// There are two kinds of ktree nodes- branch nodes and leaf nodes
// Both contain a signature matrix, plus their own signature
// (the root node signature does not matter and can be blank)
// Branch nodes then contain 'order' links to other nodes
// Leaf nodes do not.
// However, as leaf nodes may become branch nodes, we allocate
// the space anyway.
// As the space to be used is determined at runtime, we use
// parallel arrays, not structs

#ifndef INCLUDE_temp_tree_HPP
#define INCLUDE_temp_tree_HPP

#include <omp.h>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include "read.hpp"
#include "distance.hpp"

using namespace std;

// typedef unsigned char cell_type;
bloom_parameters parameters;
size_t partree_capacity = 100;

vector<cell_type> createMeanSig(const vector<seq_type> clusterSigs)
{
    vector<cell_type> meanSig(signatureSize * bits_per_char);
    vector<size_t> counter(signatureSize * bits_per_char);

    for (seq_type signature : clusterSigs)
    {
        seq_type signatureData = getMinimiseSet(signature);
        for (size_t i = 0; i < signatureSize; i++)
        {
            for (int n = 0; n < bits_per_char; n++)
            {
                if ((signatureData[0][i] >> n) & 1)
                {
                    counter[i * bits_per_char + n]++;
                }
            }
        }
    }

    for (int i = 0; i < counter.size(); i++)
    {
        if (counter[i] > clusterSigs.size() / 2)
        {
            meanSig[i / bits_per_char] |= (cell_type)1 << (i % bits_per_char);
        }
    }

    return meanSig;
}

class temp_tree
{
public:
    size_t root = 0;                   // # of root node
    vector<size_t> childCounts;        // n entries, number of children
    vector<int> isBranchNode;          // n entries, is this a branch node
    vector<vector<size_t>> childLinks; // n * o entries, links to children
    vector<vector<size_t>> seqIDs;     // n * o entries, links to children
    vector<size_t> parentLinks;        // n entries, links to parents
    vector<distance_type> priority;    // n entries, links to parents
    vector<vector<cell_type>> means;   // n * signatureSize entries, node signatures
    vector<vector<seq_type>> matrices; // capacity * signatureSize * n
    vector<omp_lock_t> locks;          // n locks
    size_t capacity = 0;               // Set during construction, currently can't change

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
                means.resize(capacity * signatureSize);
            }
        }
    }

    temp_tree(size_t capacity)
    {
        reserve(capacity);
        childCounts[root] = 0;
        isBranchNode[root] = 0;
    }

    void printMatrix(FILE *stream, size_t node)
    {
        fprintf(stream, ">>>printing matrix at node %zu\n", node);
        for (seq_type seq : matrices[node])
        {
            dbgPrintSignature(stream, seq);
        }
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
        if (isBranchNode[tnode])
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
    }

    void addSigToMatrix(size_t node, seq_type signature)
    {
        matrices[node].push_back(signature);

        // // debug
        // printMatrix(stderr, node);
    }

    // calculate mean form a vector of sigs
    inline vector<cell_type> createMeanSig(const vector<seq_type> clusterSigs)
    {
        vector<cell_type> meanSig(signatureSize * bits_per_char);
        vector<size_t> counter(signatureSize * bits_per_char);

        for (seq_type signature : clusterSigs)
        {
            seq_type signatureData = getMinimiseSet(signature);
            for (size_t i = 0; i < signatureSize; i++)
            {
                for (int n = 0; n < bits_per_char; n++)
                {
                    if ((signatureData[0][i] >> n) & 1)
                    {
                        counter[i * bits_per_char + n]++;
                    }
                }
            }
        }

        for (int i = 0; i < counter.size(); i++)
        {
            if (counter[i] > clusterSigs.size() / 2)
            {
                meanSig[i / bits_per_char] |= (cell_type)1 << (i % bits_per_char);
            }
        }

        return meanSig;
    }

    // return if found same, and the destination node
    inline tuple<bool, size_t> traverseHD(seq_type signature) const
    {
        size_t node = root;
        size_t a = countSetBits(signature);
        fprintf(stderr, " \n(%zu,%f )", a, split_threshold * (a + a));

        while (isBranchNode[node])
        {
            fprintf(stderr, " \n%zu: ", node);
            size_t lowestHD = numeric_limits<size_t>::max();
            size_t lowestHDchild = childLinks[node][0];
            size_t mismatch = 0;

            for (size_t i = 0; i < childCounts[node]; i++)
            {
                size_t child = childLinks[node][i];
                size_t hd = calcDistance(matrices[child][0], signature);
                size_t b = countSetBits(matrices[child][0]);

                fprintf(stderr, " <%zu,%zu,%zu> ", child, hd, b);

                // found same, move on with the next seq
                if (hd <= stay_threshold * (a + b))
                {
                    return make_tuple(true, child);
                }
                else if (hd < lowestHD)
                {
                    lowestHD = hd;
                    lowestHDchild = child;
                }

                // count how many nodes mismatch
                if (hd >= split_threshold * (a + b))
                {
                    mismatch++;
                }
            }

            // nothing is close enough, spawn new child under parent
            if (mismatch == childCounts[node])
            {
                return make_tuple(false, node);
            }

            // matching with all, grow on top
            if (mismatch == 0)
            {
                return make_tuple(false, parentLinks[node]);
            }

            node = lowestHDchild;
        }

        return std::make_tuple(false, node);
    }

    // return if found same, and the destination node
    inline tuple<bool, size_t> traverseInter(seq_type signature) const
    {
        size_t node = root;
        size_t a = countSetBits(signature);
        fprintf(stderr, " \n(%zu, %f, ,%f)", a, split_threshold * a, stay_threshold * a);

        while (isBranchNode[node])
        {
            fprintf(stderr, " \n%zu: ", node);
            size_t maxInter = 0;
            size_t maxInterchild = childLinks[node][0];
            size_t mismatch = 0;

            for (size_t i = 0; i < childCounts[node]; i++)
            {
                size_t child = childLinks[node][i];
                size_t inter = calcDistance(matrices[child][0], signature);
                size_t b = countSetBits(matrices[child][0]);

                fprintf(stderr, " <%zu,%zu,%zu> ", child, inter, b);

                size_t len = max(a, b);

                // found same, move on with the next seq
                if (inter >= stay_threshold * len)
                {
                    return make_tuple(true, child);
                }
                else if (inter < maxInter)
                {
                    maxInter = inter;
                    maxInterchild = child;
                }

                // count how many nodes mismatch
                if (inter <= split_threshold * len)
                {
                    mismatch++;
                }
            }

            // nothing is close enough, spawn new child under parent
            if (mismatch == childCounts[node])
            {
                return make_tuple(false, node);
            }

            node = maxInterchild;
        }

        return std::make_tuple(false, node);
    }

    // if there is only 1 path and the input is closer to the parent than the child, rotate them
    inline size_t rotate(size_t node, seq_type signature, vector<size_t> &insertionList)
    {
        size_t parent = parentLinks[node];
        size_t old_node = node;

        // find where to insert
        while (node != root)
        {
            size_t match = 0;
            for (int i = 0; i < childCounts[parent]; i++)
            {
                size_t temp = childLinks[parent][i];
                size_t distance = calcDistance(signature, matrices[temp][0]);
                //? change length later
                if (distance > stay_threshold * signature.size() && distance < split_threshold * signature.size())
                {
                    match++;
                }
            }

            if (match != childCounts[parent])
            {
                break;
            }

            old_node = node;
            node = parent;
            parent = parentLinks[node];
        }

        // insert new node as child of parent and the parent and its siblings are now the children of the new node
        size_t temp = getNewNodeIdx(insertionList);

        isBranchNode[temp] = 1;
        parentLinks[temp] = node;
        childCounts[temp] = 1;
        childLinks[temp].push_back(old_node);
        addSigToMatrix(temp, signature);

        // update parent of old_node and siblings
        for (size_t sibling : childLinks[node])
        {
            parentLinks[sibling] = temp;
        }

        // update parent
        replace(childLinks[node].begin(), childLinks[node].end(), old_node, temp);

        fprintf(stderr, "***Adding %zu between %zu, %zu, %zu", temp, parent, node, old_node);
        return temp;
    }

    inline bool isUnipath(size_t node)
    {
        size_t parent = parentLinks[node];
        size_t grandparent = parentLinks[parent];

        // // return the highest node in the unipath
        // // if output = input node => there is no unipath
        // // else the highest point is the parent of the output
        // while (childCounts[parent] == 1 && childCounts[grandparent] == 1)
        // {
        //     node = parent;
        //     parent = grandparent;
        //     grandparent = parentLinks[parent];
        // }
        // return node;

        return childCounts[parent] == 1 && childCounts[grandparent] == 1;
    }

    // priority = distance to ancestor
    inline size_t rotateAnc(size_t node)
    {
        if (node == root)
        {
            return 0;
        }
        size_t parent = parentLinks[node];

        if (parent == root)
        {
            return 0;
        }

        if (priority[parent] < priority[node])
        {
            return 0; // heap order, larger priority on top
        }

        fprintf(stderr, ">%zu; priority: %.1f,%.1f\n", parent, priority[parent], priority[node]);

        size_t grandparent = parentLinks[parent];
        vector<size_t> &siblings = childLinks[parent];
        fprintf(stderr, "\nrotate: %zu,%zu,%zu\n", grandparent, parent, node);

        // remove(siblings.begin(), siblings.end(), node);
        siblings.erase(remove(siblings.begin(), siblings.end(), node), siblings.end());

        vector<size_t> &uncles = childLinks[grandparent];
        uncles.erase(remove(uncles.begin(), uncles.end(), parent), uncles.end());

        childLinks[grandparent].push_back(node);
        childCounts[grandparent] += siblings.size();
        for (size_t sibling : siblings)
        {
            fprintf(stderr, "\nSibling: %zu\n", sibling);
            childLinks[grandparent].push_back(sibling);
            parentLinks[sibling] = grandparent;
        }

        parentLinks[node] = grandparent;
        childCounts[node]++;

        for (size_t grand : childLinks[node])
        {
            fprintf(stderr, "%zu,", grand);
        }
        fprintf(stderr, "\n");
        childLinks[node].push_back(parent);
        for (size_t grand : childLinks[node])
        {
            fprintf(stderr, "%zu,", grand);
        }
        fprintf(stderr, "\n");

        parentLinks[parent] = node;

        int temp = isBranchNode[parent];
        isBranchNode[parent] = isBranchNode[node];
        isBranchNode[node] = temp;

        fprintf(stderr, "finish: %zu\n", grandparent);

        for (size_t child : childLinks[grandparent])
        {
            fprintf(stderr, "%zu>", child);
            for (size_t grand : childLinks[child])
            {
                fprintf(stderr, "%zu,", grand);
            }
            fprintf(stderr, "\n");
        }

        return parent;
    }

    // return if found same, and the destination node
    inline tuple<bool, size_t> traverse(seq_type signature, vector<size_t> &insertionList)
    {
        size_t node = root;
        double local_stay = stay_threshold;
        double offset = 0.1;

        while (isBranchNode[node])
        {
            vector<size_t> mismatch;
            vector<size_t> matching;

            for (size_t i = 0; i < childCounts[node]; i++)
            {
                size_t child = childLinks[node][i];
                double sim = calcDistance(matrices[child][0], signature);

                fprintf(stderr, " <%zu,%.2f> ", child, sim);

                // found same, move on with the next seq
                if (sim >= (local_stay))
                {
                    // priority[child]++;
                    // rotateAnc(child);
                    return make_tuple(true, child);
                }

                // count how many nodes mismatch
                if (sim <= (split_threshold))
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
                // node = rotateAnc(node);
                return make_tuple(false, node);
            }
            else if (matching.size() > 1)
            {
                size_t temp = createNode(signature, insertionList, node);

                parentLinks[temp] = node;
                isBranchNode[temp] = 1;
                childCounts[temp] = matching.size();
                childLinks[temp] = matching;

                childCounts[node] = mismatch.size() + 1;
                childLinks[node].clear();
                childLinks[node] = mismatch;
                childLinks[node].push_back(temp);

                fprintf(stderr, "\nxxx multiple: %zu,%zu\n", temp, node);

                for (size_t n : matching)
                {
                    parentLinks[n] = temp;
                    fprintf(stderr, "leaves: %zu\n", n);
                }

                return make_tuple(true, temp);
            }
            node = matching[0]; // choose the first matching for now, may change to best
        }

        return std::make_tuple(false, node);
    }

    // return if found same, and the destination node
    inline tuple<bool, size_t> traverseMatching(seq_type signature, vector<size_t> &insertionList)
    {
        size_t node = root;
        size_t a = signature.size();
        fprintf(stderr, " \n(%zu, %f, ,%f)", a, split_threshold * a, stay_threshold * a);
        double offset = 0;

        while (isBranchNode[node])
        {
            fprintf(stderr, " \n%zu: ", node);
            vector<size_t> mismatch;
            vector<size_t> matching_leaves;
            vector<size_t> matching_branch;

            for (size_t i = 0; i < childCounts[node]; i++)
            {
                size_t child = childLinks[node][i];
                size_t inter = calcDistance(matrices[child][0], signature);
                size_t b = matrices[child][0].size();

                size_t len = max(a, b);
                fprintf(stderr, " <%zu,%zu,%.2f,%.2f> ", child, inter, (stay_threshold - offset) * len, split_threshold * len);

                // found same, move on with the next seq
                if (inter >= (stay_threshold - offset) * len)
                {
                    return make_tuple(true, child);
                }

                // count how many nodes mismatch
                if (inter <= (split_threshold)*len)
                {
                    // mismatch++;
                    mismatch.push_back(child);
                }
                else if (!isBranchNode[child])
                {
                    matching_leaves.push_back(child);
                    fprintf(stderr, " -m%zu, ", child);
                }
                else
                {
                    //? treat matching branch as mismatch for now
                    matching_branch.push_back(child);
                    fprintf(stderr, " -b%zu, ", child);
                }
            }

            // nothing is close enough, spawn new child under parent
            if (mismatch.size() == childCounts[node])
            {
                fprintf(stderr, " -s%zu, ", node);
                return make_tuple(false, node);
            }
            // else if (mismatch == 0)
            // {
            //     size_t temp = rotate(node, signature, insertionList);
            //     return make_tuple(true, temp);
            // }

            else if (matching_leaves.size() == 1)
            {
                node = matching_leaves[0];
            }
            else if (matching_branch.size() == 1)
            {
                node = matching_branch[0];
            }
            else
            {
                size_t temp = getNewNodeIdx(insertionList);
                addSigToMatrix(temp, signature);
                // priority[temp] = calcDistance(signature, matrices[findAncestor(node)][0]);
                size_t ancestor = findAncestor(node);
                if (ancestor == root)
                {
                    priority[temp] = signature.size();
                }
                else
                {
                    priority[temp] = calcDistance(signature, matrices[ancestor][0]);
                }

                parentLinks[temp] = node;
                isBranchNode[temp] = 1;
                childCounts[temp] = matching_leaves.size() + matching_branch.size();
                childLinks[temp] = matching_branch;
                childLinks[temp].insert(childLinks[temp].end(), matching_leaves.begin(), matching_leaves.end());

                childCounts[node] = mismatch.size() + 1;
                childLinks[node].clear();
                childLinks[node] = mismatch;
                childLinks[node].push_back(temp);

                for (size_t n : matching_leaves)
                {
                    parentLinks[n] = temp;
                }

                for (size_t n : matching_branch)
                {
                    parentLinks[n] = temp;
                }

                fprintf(stderr, "\nxxx multiple: %zu,%zu,%zu\n", temp, node, childLinks[node][childCounts[node] - 1]);

                return make_tuple(true, temp);
            }
            // offset += 0.1;
        }

        // grow height, should I go on top or bottom
        // return rotateAnc(node, signature, insertionList);
        return std::make_tuple(false, node);
        // return make_tuple(true, rotateAnc(node, signature, insertionList));
    }

    inline size_t createNode(seq_type signature, vector<size_t> &insertionList, size_t parent)
    {
        // size_t parent = root;
        size_t node = getNewNodeIdx(insertionList);
        parentLinks[node] = parent;
        childLinks[parent].push_back(node);
        isBranchNode[parent] = 1;
        childCounts[parent]++;
        // priority[node] = 1;
        addSigToMatrix(node, signature);
        means[node] = getMinimiseSet(signature)[0];
        priority[node] = countSingleSetBits(means[node]);
        return node;
    }

    inline size_t first_insert(seq_type signature, vector<size_t> &insertionList, size_t idx, size_t parent = 0)
    {
        size_t node = createNode(signature, insertionList, parent);
        seqIDs[node].push_back(idx);
        return node;
    }

    inline size_t insert(seq_type signature, vector<size_t> &insertionList, size_t idx)
    {
        bool stay = false;
        size_t parent = root;
        tie(stay, parent) = traverse(signature, insertionList);
        fprintf(stderr, "inserting seq %zu at node %zu; stay: %d\n", idx, parent, stay);

        if (!stay) // add new node
        {
            parent = createNode(signature, insertionList, parent);
            seqIDs[parent].push_back(idx);
        }
        else // add seq without adding new node
        {
            seqIDs[parent].push_back(idx);
            addSigToMatrix(parent, signature);
            means[parent] = createMeanSig(matrices[parent]);
            priority[parent] = countSingleSetBits(means[parent]);
            // priority[parent]++;
            // return parent;
        }

        fprintf(stderr, "\ninserting %zu at %zu\n", idx, parent);
        return parent;

        // size_t insertionPoint = traverse(signature);

        // // fprintf(stdout, " at %zu\n", insertionPoint);
        // omp_set_lock(&locks[insertionPoint]);

        // // if (childCounts[insertionPoint] < order)
        // {
        //     addSigToMatrix(insertionPoint, signature);

        //     childLinks[insertionPoint].push_back(idx);
        //     fprintf(stdout, " at %zu\n", insertionPoint);

        //     childCounts[insertionPoint]++;
        //     recalculateSig(insertionPoint);
        //     if (isBranchNode[root])
        //     {
        //         recalculateUp(parentLinks[insertionPoint]);
        //     }
        // }
        // omp_unset_lock(&locks[insertionPoint]);
    }

    inline void removeSingleton(vector<tuple<size_t, size_t>> &clusters, vector<size_t> &insertionList)
    {
        vector<size_t> tempChildren;
        vector<size_t> singletons;
        for (size_t i = 0; i < childCounts[root]; i++)
        {
            size_t node = childLinks[root][i];
            if (childCounts[node] > 1)
            {
                tempChildren.push_back(node);
            }
            else
            {
                singletons.push_back(node);
            }
        }
        childLinks[root].clear();
        childLinks[root] = tempChildren;
        childCounts[root] = tempChildren.size();

        for (size_t node : singletons)
        {
            size_t i = seqIDs[node][0];
            // insertionList.push_back(node);
            size_t clu = insert(matrices[node][0], insertionList, i);
            clusters[i] = make_tuple(clu, findAncestor(clu));
        }
    }

    inline size_t search(seq_type signature, size_t idx)
    {
        size_t node = root;
        double local_stay = stay_threshold;

        size_t best_child = node;
        double best_sim = 0;

        while (isBranchNode[node])
        {
            size_t local_best_child = node;
            double local_best_sim = 0;

            for (size_t i = 0; i < childCounts[node]; i++)
            {
                size_t child = childLinks[node][i];
                double sim = calcDistance(matrices[child][0], signature);

                fprintf(stderr, " <%zu,%.2f> ", child, sim);

                // found same, move on with the next seq
                if (sim >= (local_stay))
                {
                    fprintf(stderr, "\n found %zu at %zu\n>", idx, child);

                    for (size_t id : seqIDs[child])
                    {
                        if (id == idx)
                        {
                            fprintf(stderr, "***match\n");
                        }
                    }
                    return child;
                }

                if (sim >= local_best_sim)
                {
                    local_best_sim = sim;
                    local_best_child = child;
                }
            }
            if (local_best_sim >= best_sim)
            {
                best_sim = local_best_sim;
                best_child = local_best_child;
            }

            node = local_best_child;
        }

        fprintf(stderr, "\n found %zu at %zu\n>", idx, best_child);

        for (size_t id : seqIDs[best_child])
        {
            if (id == idx)
            {
                fprintf(stderr, "***match\n");
            }
        }
        return best_child;
    }

    // // return if found same, and the destination node
    // inline size_t search(seq_type signature, size_t idx) const
    // {
    //     size_t node = root;
    //     size_t a = countSetBits(signature);

    //     while (isBranchNode[node])
    //     {
    //         // fprintf(stderr, " \n%zu: ", node);
    //         size_t lowestHD = numeric_limits<size_t>::max();
    //         size_t lowestHDchild = childLinks[node][0];
    //         size_t mismatch = 0;

    //         for (size_t i = 0; i < childCounts[node]; i++)
    //         {
    //             size_t child = childLinks[node][i];
    //             size_t hd = calcDistance(matrices[child][0], signature);
    //             size_t b = countSetBits(matrices[child][0]);

    //             // found same, move on with the next seq
    //             if (hd <= stay_threshold * max(a, b))
    //             {
    //                 return child;
    //             }
    //             else if (hd < lowestHD)
    //             {
    //                 lowestHD = hd;
    //                 lowestHDchild = child;
    //             }

    //             // count how many nodes mismatch
    //             if (hd >= split_threshold * max(a, b))
    //             {
    //                 mismatch++;
    //             }

    //             // fprintf(stderr, " <%zu,%zu> ", child, hd);
    //         }

    //         //? nothing is close enough, return parent
    //         if (mismatch == childCounts[node])
    //         {
    //             fprintf(stderr, ">%zu mismatch\n", idx);
    //             return parentLinks[node];
    //         }

    //         node = lowestHDchild;
    //     }

    //     return node;
    // }

    void destroyLocks(size_t node)
    {
        omp_destroy_lock(&locks[node]);
        if (isBranchNode[node])
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