// There are two kinds of ktree nodes- branch nodes and leaf nodes
// Both contain a signature matrix, plus their own signature
// (the root node signature does not matter and can be blank)
// Branch nodes then contain 'order' links to other nodes
// Leaf nodes do not.
// However, as leaf nodes may become branch nodes, we allocate
// the space anyway.
// As the space to be used is determined at runtime, we use
// parallel arrays, not structs

#ifndef INCLUDE_test_tree_HPP
#define INCLUDE_test_tree_HPP

#include <omp.h>
#include <unordered_map>
#include <unordered_set>
#include <experimental/filesystem>
#include "bloom_filter.hpp"

// #include "read.hpp"

using namespace std;
typedef pair<float, float> data_type;

double split_threshold = 5;
double stay_threshold = 1;
size_t minimiser_match_threshold = 4;
size_t partree_capacity = 100;

data_type createMeanSig(const vector<data_type> clusterSigs)
{
    data_type meanSig = make_pair(0, 0);

    for (data_type signature : clusterSigs)
    {
        meanSig.first += signature.first;
        meanSig.second += signature.second;
    }

    meanSig.first = meanSig.first / clusterSigs.size();
    meanSig.second = meanSig.second / clusterSigs.size();

    return meanSig;
}

double calcDistance(data_type a, data_type b)
{
    float x = a.first - b.first;
    float y = a.second - b.second;
    return sqrt(x * x + y * y);
}

// RMSD
double calcDistortion(const vector<data_type> clusterSigs)
{
    data_type meanSig = createMeanSig(clusterSigs);
    double sumSquareDistance = 0;

    for (data_type signature : clusterSigs)
    {
        double distance = calcDistance(meanSig, signature);
        sumSquareDistance += distance * distance;
    }
    return sqrt(sumSquareDistance / clusterSigs.size());
}

class test_tree
{
public:
    size_t root = 0;                    // # of root node
    vector<size_t> childCounts;         // n entries, number of children
    vector<int> isBranchNode;           // n entries, is this a branch node
    vector<vector<size_t>> childLinks;  // n * o entries, links to children
    vector<vector<size_t>> seqIDs;      // n * o entries, links to children
    vector<size_t> parentLinks;         // n entries, links to parents
    vector<float> priority;             // n entries, links to parents
    vector<data_type> means;            // n * signatureSize entries, node signatures
    vector<vector<data_type>> matrices; // capacity * signatureSize * n
    vector<omp_lock_t> locks;           // n locks
    size_t capacity = 0;                // Set during construction, currently can't change

    void reserve(size_t capacity)
    {
        // For safety, only call this at startup currently
        if (this->capacity != 0)
        {
            fprintf(stderr, "Reserve can only be called from 0 capacity\n");
            exit(1);
        }
        this->capacity = capacity;

        //#pragma omp parallel
        {
            //#pragma omp single
            {
                childCounts.resize(capacity);
            }
            //#pragma omp single
            {
                isBranchNode.resize(capacity);
            }
            //#pragma omp single
            {
                childLinks.resize(capacity);
            }
            //#pragma omp single
            {
                seqIDs.resize(capacity);
            }
            //#pragma omp single
            {
                priority.resize(capacity);
            }
            //#pragma omp single
            {
                parentLinks.resize(capacity);
            }
            //#pragma omp single
            {
                locks.resize(capacity);
            }
            //#pragma omp single
            {
                matrices.resize(capacity);
            }
            //#pragma omp single
            {
                means.resize(capacity);
            }
        }
    }

    test_tree(size_t capacity)
    {
        reserve(capacity);
        childCounts[root] = 0;
    }

    size_t getNewNodeIdx(vector<size_t> &insertionList)
    {
        if (insertionList.empty())
        {
            fprintf(stderr, "ERROR: ran out of insertion points\n");
            exit(1);
        }
        size_t idx = insertionList.back();
        insertionList.pop_back();

        // Initialise lock
        omp_init_lock(&locks[idx]);
        return idx;
    }

    size_t findAncestor(size_t node)
    {
        while (parentLinks[node] != root)
        {
            node = parentLinks[node];
        }
        return node;
    }

    tuple<size_t, size_t> findAncestorNlevel(size_t node)
    {
        size_t level = 0;
        while (parentLinks[node] != root)
        {
            node = parentLinks[node];
            level++;
        }
        return make_tuple(node, level);
    }

    size_t findLevel(size_t node)
    {
        size_t level = 0;
        while (node != root)
        {
            node = parentLinks[node];
            level++;
        }
        return level;
    }

    void printNodeJson(FILE *stream, size_t tnode)
    {
        fprintf(stream, "{\"node\":\"%zu\",\"priority\":\"%.2f\",\"childCount\":\"%zu\",\"content\":\"*", tnode, priority[tnode], seqIDs[tnode].size());
        for (size_t seq : seqIDs[tnode])
        {
            fprintf(stream, "%zu,", seq);
        }
        fprintf(stream, "\",\"children\":[");
    }

    void printSubTreeJson(FILE *stream, size_t tnode)
    {
        if (childCounts[tnode]>0)
        {
            printNodeJson(stream, tnode);
            printSubTreeJson(stream, childLinks[tnode][0]);

            for (size_t i = 1; i < childLinks[tnode].size(); i++)
            {
                fprintf(stream, ",");
                printSubTreeJson(stream, childLinks[tnode][i]);
            }
            fprintf(stream, "]}");
        }
        else
        {
            fprintf(stream, "{\"node\":\"%zu\",\"priority\":\"%.2f\",\"childCount\":\"%zu\",\"content\":\"*", tnode, priority[tnode], seqIDs[tnode].size());
            for (size_t seq : seqIDs[tnode])
            {
                fprintf(stream, "%zu,", seq);
            }
            fprintf(stream, "\"}");
        }
    }

    void printTreeJson(FILE *stream)
    {
        fprintf(stream, "var treeData = ");
        printSubTreeJson(stream, root);
        fprintf(stream, ";\n");
    }

    void addSigToMatrix(size_t node, data_type signature)
    {
        matrices[node].push_back(signature);
    }

    void clearNode(size_t node)
    {
        isBranchNode[node] = 0;
        childCounts[node] = 0;
        childLinks[node].clear();
        seqIDs[node].clear();
        matrices[node].clear();
        means[node] = createMeanSig(matrices[node]);
        priority[node] = 0;
    }

    // merge if number of children is greater than tree order
    // decrease split threshold
    // to allow node with larger distortion
    // the local split threshold should be gradual, but we just make it 110% of the og threshold
    // stay threshold should have nothing to do here, if not split just stay
    // at this stage, we don't do k-means
    // implement as if this is an independant tree
    void mergeChildren(size_t node, vector<size_t> &insertionList)
    {
        // double local_split = split_threshold * 1.1;

        double local_split = minimiser_match_threshold;

        vector<data_type> temp_means;
        vector<size_t> temp_childLink;

        // insert first child
        size_t child = childLinks[node][0];
        temp_childLink.push_back(child);
        temp_means.push_back(means[child]);

        for (size_t i = 1; i < childCounts[node]; i++)
        {
            child = childLinks[node][i];
            size_t n = 0;
            for (data_type mean : temp_means)
            {
                float distance = calcDistance(means[child], mean);
                // fprintf(stderr, "Distance:%.2f, %.2f, %.2f\n", distance,split_threshold, local_split);
                if (distance < local_split)
                {
                    fprintf(stderr, "Merge %zu into %zu\n", child, childLinks[node][n]);
                    break;
                }
                n++;
            }
            fprintf(stderr, "n:%zu\n", n);
            if (n == temp_means.size())
            {
                temp_childLink.push_back(child);
                temp_means.push_back(means[child]);
            }
        }

        fprintf(stderr, "Before: %zu, After:%zu\n", childCounts[node], temp_childLink.size());
    }

    // union mean of children
    inline void updateNodeMean(size_t node)
    {
        // if (childCounts[node] > 1)
        // {
        //     size_t size = signatureSize * bits_per_char;
        //     data_type meanSig(size);
        //     fill(&meanSig[0], &meanSig[0] + size, 0);

        //     for (size_t child : childLinks[node])
        //     {
        //         for (size_t i = 0; i < size; i++)
        //         {
        //             meanSig[i] |= means[child][i];
        //         }
        //     }
        //     means[node] = meanSig;
        //     priority[node] = countSingleSetBits(means[node]);
        // }
        means[node] = createMeanSig(matrices[node]);
        //? p
        priority[node] = calcDistortion(matrices[node]);
    }

    // return if found same, and the destination node
    inline tuple<bool, size_t> traverse(data_type signature, vector<size_t> &insertionList)
    {
        size_t node = root;
        double local_stay = stay_threshold;
        double offset = 0.1;

        while (childCounts[node]>0)
        {
            vector<size_t> mismatch;
            vector<size_t> matching;

            for (size_t i = 0; i < childCounts[node]; i++)
            {
                size_t child = childLinks[node][i];
                double distance = calcDistance(means[child], signature);

                fprintf(stderr, " <%zu,%.2f> ", child, distance);

                // found same, move on with the next seq
                if (distance <= (local_stay))
                {
                    // priority[child]++;
                    return make_tuple(true, child);
                }

                // count how many nodes mismatch
                if (distance >= (split_threshold))
                {
                    // mismatch++;
                    mismatch.push_back(child);
                }
                else
                {
                    matching.push_back(child);
                }
            }

            // nothing is close enough, spawn new child under parent
            if (mismatch.size() == childCounts[node])
            {
                fprintf(stderr, " -no match:%zu, ", node);
                return make_tuple(false, node);
            }
            else if (matching.size() > 1) // matching with multiple
            {
                size_t temp = createNode(signature, insertionList, node);

                parentLinks[temp] = node;

                childCounts[node] = mismatch.size() + 1;
                childLinks[node].clear();
                childLinks[node] = mismatch;
                childLinks[node].push_back(temp);

                fprintf(stderr, "\nxxx multiple: %zu,%zu\n", temp, node);

                for (size_t n : matching)
                {
                    for (data_type sig : matrices[n])
                    {
                        addSigToMatrix(temp, sig);
                    }

                    for (size_t grandchild : childLinks[n])
                    {
                        parentLinks[grandchild] = temp;
                        childLinks[temp].push_back(grandchild);
                    }

                    for (size_t id : seqIDs[n])
                    {
                        seqIDs[temp].push_back(id);
                    }

                    childCounts[temp] += childCounts[n];
                    fprintf(stderr, "leaves: %zu,%zu\n", n, childCounts[temp]);
                    clearNode(n);
                    insertionList.push_back(n);
                }
                // fprintf(stderr, "leaves: %zu\n",childCounts[temp]);
                // printSubTreeJson(stderr,root);

                //?
                priority[temp] = childCounts[temp] - 1;

                return make_tuple(true, temp);
            }
            // else if (matching.size() == 1 && childCounts[node] == 1) // matching with multiple
            // {
            //     size_t temp = createNode(signature, insertionList, node);

            //     parentLinks[temp] = node;
            //     isBranchNode[temp] = 0;

            //     childCounts[node] = 1;
            //     childLinks[node].clear();
            //     childLinks[node].push_back(temp);

            //     childCounts[temp] = 2;
            //     childLinks[node].push_back(matching[0]);
            //     childLinks[temp].push_back(temp);

            //     fprintf(stderr, "\nxxx here: %zu,%zu\n", temp, node);

            //     for (size_t n : matching)
            //     {
            //         for (data_type sig : matrices[n])
            //         {
            //             addSigToMatrix(temp, sig);
            //         }

            //         for (size_t grandchild : childLinks[n])
            //         {
            //             parentLinks[grandchild] = temp;
            //             childLinks[temp].push_back(grandchild);
            //         }

            //         for (size_t id : seqIDs[n])
            //         {
            //             seqIDs[temp].push_back(id);
            //         }

            //         childCounts[temp] += childCounts[n];
            //         fprintf(stderr, "leaves: %zu,%zu\n", n, childCounts[temp]);
            //         clearNode(n);
            //         insertionList.push_back(n);
            //     }
            //     // fprintf(stderr, "leaves: %zu\n",childCounts[temp]);
            //     // printSubTreeJson(stderr,root);

            //     //?
            //     priority[temp] = childCounts[temp] - 1;

            //     return make_tuple(true, temp);
            // }

            // updateNodeMean(node);
            node = matching[0]; // choose the first matching for now, may change to best
        }

        return std::make_tuple(false, node);
    }

    inline size_t createNode(data_type signature, vector<size_t> &insertionList, size_t parent)
    {
        // size_t parent = root;
        size_t node = getNewNodeIdx(insertionList);
        parentLinks[node] = parent;
        childLinks[parent].push_back(node);
        childCounts[parent]++;
        // priority[node] = 1;
        addSigToMatrix(node, signature);
        means[node] = signature;
        //? p
        // priority[node] = countSingleSetBits(means[node]);
        // priority[node] = 1;
        return node;
    }

    inline size_t first_insert(data_type signature, vector<size_t> &insertionList, size_t idx, size_t parent = 0)
    {
        size_t node = createNode(signature, insertionList, parent);
        seqIDs[node].push_back(idx);
        return node;
    }

    inline size_t insert(data_type signature, vector<size_t> &insertionList, size_t idx)
    {
        bool stay = false;
        size_t parent = root;
        tie(stay, parent) = traverse(signature, insertionList);
        fprintf(stderr, "inserting seq %zu at node %zu; stay: %d\n", idx, parent, stay);

        if (!stay) // add new node
        {
            size_t node = createNode(signature, insertionList, parent);
            seqIDs[node].push_back(idx);
            return node;
        }
        else // add seq without adding new node
        {
            seqIDs[parent].push_back(idx);
            addSigToMatrix(parent, signature);
            means[parent] = createMeanSig(matrices[parent]);
            // priority[parent] = countSingleSetBits(means[parent]);
            //? p
            priority[parent] = calcDistortion(matrices[parent]);
            // priority[parent]++;

            // priority[parent]++;
            // return parent;
        }

        fprintf(stderr, "\ninserting %zu at %zu\n", idx, parent);
        return parent;
    }

    inline size_t search(data_type signature, size_t idx = 0)
    {
        size_t node = root;

        size_t best_child = node;
        double best_distance = 0;

        while (childCounts[node]>0)
        {
            size_t local_best_child = node;
            double local_best_distance = 0;

            for (size_t i = 0; i < childCounts[node]; i++)
            {
                size_t child = childLinks[node][i];
                double distance = calcDistance(means[child], signature);

                fprintf(stderr, " <%zu,%.2f> ", child, distance);

                // found same, move on with the next seq
                if (distance <= (stay_threshold))
                {
                    // seqIDs[child].push_back(idx);
                    return child;
                }

                if (distance <= local_best_distance)
                {
                    local_best_distance = distance;
                    local_best_child = child;
                }
            }
            if (local_best_distance <= best_distance)
            {
                best_distance = local_best_distance;
                best_child = local_best_child;
            }

            node = local_best_child;
        }

        // seqIDs[best_child].push_back(idx);
        return best_child;
    }

    void clearSeqId(size_t lastNode)
    {
        for (size_t i = 0; i < lastNode; i++)
        {
            seqIDs[i].clear();
        }
    }

    // make sure cluster centroid is up-to-date before clearing matrices
    // clear matrices and seqID for the next insertion cycle
    // cluster means still remain
    void prepReinsert(size_t node)
    {

        updateNodeMean(node);
        matrices[node].clear();
        seqIDs[node].clear();

        if (childCounts[node]>0)
        {
            for (size_t child : childLinks[node])
            {
                prepReinsert(child);
            }
        }
    }

    size_t reinsert(data_type signature, size_t idx)
    {
        size_t node = search(signature);
        seqIDs[node].push_back(idx);
        addSigToMatrix(node, signature);
        return node;
    }

    void destroyLocks(size_t node)
    {
        omp_destroy_lock(&locks[node]);
        if (childCounts[node]>0)
        {
            for (size_t i = 0; i < childCounts[node]; i++)
            {
                destroyLocks(childLinks[node][i]);
            }
        }
    }

    inline void destroyLocks()
    {
        destroyLocks(root);
    }
};

#endif