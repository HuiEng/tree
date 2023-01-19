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

// vector<cell_type> createMeanSig(const vector<seq_type> clusterSigs)
// {
//     vector<cell_type> meanSig(signatureSize * bits_per_char);
//     vector<size_t> counter(signatureSize * bits_per_char);

//     for (seq_type signature : clusterSigs)
//     {
//         seq_type signatureData = getMinimiseSet(signature);
//         for (size_t i = 0; i < signatureSize; i++)
//         {
//             for (int n = 0; n < bits_per_char; n++)
//             {
//                 if ((signatureData[0][i] >> n) & 1)
//                 {
//                     counter[i * bits_per_char + n]++;
//                 }
//             }
//         }
//     }

//     for (int i = 0; i < counter.size(); i++)
//     {
//         if (counter[i] >= clusterSigs.size() / 2)
//         {
//             meanSig[i / bits_per_char] |= (cell_type)1 << (i % bits_per_char);
//         }
//     }

//     return meanSig;
// }

// RMSD
double calcDistortion(const vector<seq_type> clusterSigs)
{
    return 0;
    // vector<cell_type> meanSig = createMeanSig(clusterSigs);
    // double sumSquaresimilarity = 0;

    // for (seq_type signature : clusterSigs)
    // {
    //     vector<cell_type> signatureData = getMinimiseSet(signature)[0];
    //     double similarity = calcJaccardGlobal(meanSig, signatureData);
    //     sumSquaresimilarity += similarity * similarity;
    // }
    // return sqrt(sumSquaresimilarity / clusterSigs.size());
}

class temp_tree
{
public:
    size_t root = 0;                   // # of root node
    vector<size_t> childCounts;        // n entries, number of children
    vector<int> isBranchNode;          // n entries, is this a branch node
    vector<int> rootNodes;             // n entries, is this root of subtree
    vector<vector<size_t>> childLinks; // n * o entries, links to children
    vector<int> isAmbiNode;            // n entries, is this a branch node
    vector<vector<size_t>> ambiLinks;  // n * o entries, links to children
    vector<vector<size_t>> seqIDs;     // n * o entries, links to children
    vector<size_t> parentLinks;        // n entries, links to parents
    vector<distance_type> priority;    // n entries, links to parents
    // vector<vector<cell_type>> means;   // n * signatureSize entries, node signatures
    vector<seq_type> means;            // n * signatureSize entries, node signatures
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
                isAmbiNode.resize(capacity);
                ambiLinks.resize(capacity);
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
                // means.resize(capacity * signatureSize);
                means.resize(capacity);
            }
        }
    }

    temp_tree(size_t capacity)
    {
        reserve(capacity);
        childCounts[root] = 0;
        isBranchNode[root] = 0;
        rootNodes.push_back(root);
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
        fprintf(stream, "{\"node\":\"%zu\",", tnode);
        fprintf(stream, "\"branch\":\"%zu\",", isBranchNode[tnode]);
        fprintf(stream, "\"ambi\":\"%zu\",", isAmbiNode[tnode]);
        fprintf(stream, "\"priority\":\"%.2f\",", priority[tnode]);
        fprintf(stream, "\"childCount\":\"%zu\",\"content\":\"*", seqIDs[tnode].size());
        // fprintf(stream, "{\"node\":\"%zu\",\"branch\":\"%zu\",\"priority\":\"%.2f\",\"childCount\":\"%zu\",\"content\":\"*", tnode, isBranchNode[tnode], calcNodeMaxsimilarity(tnode), seqIDs[tnode].size());

        for (size_t seq : seqIDs[tnode])
        {
            fprintf(stream, "%zu,", seq);
        }
        // fprintf(stream, "\",\"children\":[");
    }

    void printSubTreeJson(FILE *stream, size_t tnode)
    {
        if (childCounts[tnode] > 0)
        {
            printNodeJson(stream, tnode);
            fprintf(stream, "\",\"children\":[");

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
            // fprintf(stream, "{\"node\":\"%zu\",\"priority\":\"%.2f\",\"childCount\":\"%zu\",\"content\":\"*", tnode, priority[tnode], seqIDs[tnode].size());
            // // fprintf(stream, "{\"node\":\"%zu\",\"priority\":\"%.2f\",\"childCount\":\"%zu\",\"content\":\"*", tnode, calcNodeMaxsimilarity(tnode), seqIDs[tnode].size());

            // for (size_t seq : seqIDs[tnode])
            // {
            //     fprintf(stream, "%zu,", seq);
            // }
            printNodeJson(stream, tnode);

            fprintf(stream, "\"}");
        }
    }

    void printTreeJson(FILE *stream, size_t node = 0)
    {
        fprintf(stream, "var treeData = ");
        printSubTreeJson(stream, node);
        fprintf(stream, ";\n");
    }

    void addSigToMatrix(size_t node, seq_type signature)
    {
        matrices[node].push_back(signature);

        // // debug
        // printMatrix(stderr, node);
    }

    int getNodeIdx(size_t node)
    {
        size_t parent = parentLinks[node];
        int idx = -1;
        for (size_t i = 0; i < childCounts[parent]; i++)
        {
            if (childLinks[parent][i] == node)
            {
                idx = i;
                break;
            }
        }
        return idx;
    }

    // delete a child from its parent, need to format the child separately
    inline void deleteNode(size_t node)
    {
        fprintf(stderr, "deleting %zu\n", node);
        size_t parent = parentLinks[node];

        int idx = getNodeIdx(node);

        if (node == 71)
        {
            printTreeJson(stderr, root);
            fprintf(stderr, "ERROR  %zu  %zu!!\n", node, matrices[parent].size());
        }

        if (idx == -1)
        {
            fprintf(stderr, "ERROR deleting %zu from %zu!!\n", node, parent);
        }
        else
        {
            matrices[parent].erase(matrices[parent].begin() + idx);
            childLinks[parent].erase(childLinks[parent].begin() + idx);
            childCounts[parent]--;
            // parentLinks[node] = 0;
            // isBranchNode[node] = 0;
        }
    }

    void clearNode(size_t node)
    {
        isBranchNode[node] = 0;
        childCounts[node] = 0;
        childLinks[node].clear();
        seqIDs[node].clear();
        matrices[node].clear();
        means[node].clear(); //?
        priority[node] = 0;
        parentLinks[node] = 0;
    }

    // union mean of children
    inline void updateNodeMean(size_t node)
    {
        if (isBranchNode[node])
        {
            seq_type meanSig = means[childLinks[node][0]];

            for (size_t c = 1; c < childCounts[node]; c++)
            {
                size_t child = childLinks[node][c];

                // skip ambi
                if (isAmbiNode[child])
                {
                    continue;
                }
                meanSig = doUnion(meanSig, means[child]);
            }

            // size_t winNum = meanSig.size();

            // vector<vector<size_t>> counters;
            // for (size_t w = 0; w < winNum; w++)
            // {
            //     vector<size_t> counter(signatureSize * bits_per_char);
            //     counters.push_back(counter);
            // }

            // for (size_t c = 1; c < childCounts[node]; c++)
            // {
            //     seq_type signatureData = means[childLinks[node][c]];
            //     for (size_t w = 0; w < min(winNum, signatureData.size()); w++)
            //     {
            //         for (size_t i = 0; i < signatureSize; i++)
            //         {
            //             for (int n = 0; n < bits_per_char; n++)
            //             {
            //                 if ((signatureData[w][i] >> n) & 1)
            //                 {
            //                     counters[w][i * bits_per_char + n]++;
            //                 }
            //             }
            //         }
            //     }
            // }

            // for (size_t w = 0; w < winNum; w++)
            // {
            //     vector<size_t> counter = counters[w];
            //     for (int i = 0; i < counter.size(); i++)
            //     {
            //         // make it upperbound
            //         if (counter[i] >= (childCounts[node] + 1) / 2)
            //         {
            //             meanSig[w][i / bits_per_char] |= (cell_type)1 << (i % bits_per_char);
            //         }
            //     }
            // }
            means[node] = meanSig;
        }
    }

    // add priority to the branch child
    double calcNodeMaxsimilarity(size_t node)
    {
        seq_type meanSig = means[node];

        //?
        if (isBranchNode[node])
        {
            meanSig = means[childLinks[node][0]];
        }

        double min_similarity = 1;
        size_t i = 0;

        // fprintf(stderr, ">>%f,%f\n", meanSig.first, meanSig.second);
        for (seq_type sig : matrices[node])
        {
            double similarity = calcSimilarity(meanSig, sig);

            // if (isBranchNode[node])
            // {
            //     fprintf(stderr, "%zu,%f,%f,%f\n", childLinks[node][i], sig.first, sig.second, similarity);
            // }
            // else
            // {
            //     fprintf(stderr, "%zu,%f,%f,%f\n", seqIDs[node][i], sig.first, sig.second, similarity);
            // }
            if (similarity < min_similarity)
            {
                min_similarity = similarity;
            }
            i++;
        }

        return -min_similarity;
    }

    // union mean of children
    inline void updatePriority(size_t node)
    {
        // priority[node] = calcDistortion(matrices[node]);
        // priority[node] = calcNodeMaxsimilarity(node);

        // priority[node] = calcNodeDistortion(node);
    }

    // priority = similarity to ancestor
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

    // return if match with branch centroid
    // return false if node is not branch
    // user set thresholds for leaves
    // thresholds for branch = threshold set by user + priority or distortion of children
    // if stay is leave, just add new seq in; if stay is branch, can choose to check it again
    // else return node, the vectors will be updated as well
    inline bool checkNode(size_t node, seq_type signature, vector<size_t> &insertionList, vector<size_t> &mismatch, vector<size_t> &NN, vector<size_t> &stay)
    {
        // if (childCounts[node] > 5)
        // {
        //     forcesplitBranch(node, insertionList);
        // }

        while (splitBranch(node, insertionList))
        {
            // splitBranch(node, insertionList);
            // return checkNode(node, signature, insertionList, mismatch, NN, stay);
        }

        bool matchBranchCentroid = false;
        size_t i = 0;
        if (isBranchNode[node])
        {
            // skip first child later;
            i++;

            // first child is the centroid, deal with it separately
            size_t child = childLinks[node][0];
            fprintf(stderr, "checking %zu \n", child);
            splitBranch(child, insertionList);
            double local_stay_t = stay_threshold;

            // increase stay threshold of a node
            if (isBranchNode[child])
            {
                local_stay_t += priority[child];
            }

            double similarity = calcSimilarity(means[child], signature);

            if (similarity >= local_stay_t)
            {
                // do not need to merge NN, should have done so in previous operation
                fprintf(stderr, " stay in centroid %zu,%.2f \n", child, similarity);
                matchBranchCentroid = true;
            }
        }

        // check the other children
        for (i; i < childCounts[node]; i++)
        {
            size_t child = childLinks[node][i];
            double local_stay_t = stay_threshold;
            double local_split_t = split_threshold;

            // increase stay threshold of a node
            if (isBranchNode[child])
            {
                local_stay_t += priority[child];
                local_split_t += priority[child];
            }

            double similarity = calcSimilarity(means[child], signature);

            if (similarity >= local_stay_t)
            {
                fprintf(stderr, "<%zu,%.2f>: %f stay\n", child, similarity, local_stay_t);
                stay.push_back(child);
                // dest = child;
            }

            else if (similarity < local_split_t)
            {
                fprintf(stderr, "<%zu,%.2f>: %f miss\n", child, similarity, local_split_t);
                mismatch.push_back(child);
            }
            else
            {
                fprintf(stderr, "<%zu,%.2f>: NN\n", child, similarity);
                NN.push_back(child);
            }
        }

        return matchBranchCentroid;
    }

    inline size_t stayNode(seq_type signature, vector<size_t> &insertionList, size_t idx, size_t node)
    {
        if (isBranchNode[node])
        {
            return tt_branch(signature, insertionList, idx, node);
        }
        else
        {
            seqIDs[node].push_back(idx);
            addSigToMatrix(node, signature);
            updatePriority(node);
            return node;
        }
    }

    void moveParent(size_t child, size_t new_parent, bool d = true)
    {
        childCounts[new_parent]++;
        childLinks[new_parent].push_back(child);
        addSigToMatrix(new_parent, means[child]);
        if (d)
        {
            deleteNode(child);
        }
        parentLinks[child] = new_parent;
    }

    // haven't update childLinks[t_parent] yet
    inline size_t createParent(size_t node, vector<size_t> &insertionList)
    {
        size_t t_parent = getNewNodeIdx(insertionList);
        isBranchNode[t_parent] = 1;
        parentLinks[t_parent] = node;

        // update node
        childCounts[node]++;
        childLinks[node].push_back(t_parent);

        fprintf(stderr, "create t_parent %zu\n", t_parent);
        return t_parent;
    }

    inline size_t createNode(seq_type signature, vector<size_t> &insertionList, size_t t_parent, size_t idx)
    {
        size_t new_node = getNewNodeIdx(insertionList);
        parentLinks[new_node] = t_parent;
        addSigToMatrix(new_node, signature);
        seqIDs[new_node].push_back(idx);
        means[new_node] = signature;

        // then add new node to t_parent
        childCounts[t_parent]++;
        childLinks[t_parent].push_back(new_node);
        addSigToMatrix(t_parent, means[new_node]);

        if (isBranchNode[t_parent])
        {
            updatePriority(t_parent);
        }

        fprintf(stderr, "create new node %zu\n", new_node);
        return new_node;
    }

    // Ambi Node cannot be the first child or else update node mean will have problem
    inline size_t createAmbiNode(seq_type signature, vector<size_t> &insertionList, size_t node, size_t idx)
    {
        size_t dest = createNode(signature, insertionList, node, idx);
        fprintf(stderr, ">>Ambi\n");
        ambiLinks[node].push_back(dest);
        isAmbiNode[dest] = 1;
        return dest;
    }

    size_t splitBranch(size_t node, vector<size_t> &insertionList)
    {

        if (!isBranchNode[node])
        {
            return 0;
        }
        else if (priority[node] < split_node_threshold)
        {
            return 0;
        }
        fprintf(stderr, "splitbranch node %zu, %f\n", node, priority[node]);

        vector<size_t> temp_centroids = {childLinks[node][0]};
        vector<vector<size_t>> clusters;
        vector<size_t> temp;
        clusters.push_back(temp);
        vector<size_t> candidates_idx;
        for (size_t i = childCounts[node] - 1; i > 0; i--)
        {
            size_t child = childLinks[node][i];
            bool add = true;
            for (size_t centroid : temp_centroids)
            {
                double similarity = calcSimilarity(means[centroid], matrices[node][i]);
                // fprintf(stderr, "%zu, %zu, %f\n", child, centroid, similarity);
                if (similarity > split_node_threshold)
                {
                    add = false;
                    break;
                }
            }
            if (add)
            {
                // fprintf(stderr, "----- %zu\n", child);
                temp_centroids.push_back(child);
                vector<size_t> t;
                clusters.push_back(t);
            }
            else
            {
                candidates_idx.push_back(i);
            }
        }

        if (temp_centroids.size() == 1)
        {
            fprintf(stderr, "??cannot splitbranch node %zu, %f\n", node, priority[node]);
            return 0;
        }

        for (size_t n : candidates_idx)
        {
            size_t child = childLinks[node][n];

            size_t dest = 0;
            double max_similarity = 0;
            for (size_t i = 0; i < temp_centroids.size(); i++)
            {
                double similarity = calcSimilarity(means[temp_centroids[i]], matrices[node][n]);
                // fprintf(stderr, "%zu, %zu, %f\n", child, temp_centroids[i], similarity);
                if (similarity >= max_similarity)
                {
                    max_similarity = similarity;
                    dest = i;
                }
            }

            // fprintf(stderr, "--- %zu goes to %zu\n", child, dest);
            clusters[dest].push_back(child);
        }

        if (temp_centroids.size() == 2)
        {
            for (vector<size_t> cluster : clusters)
            {
                if (cluster.size() == 0)
                {
                    fprintf(stderr, ">>??something is wrong cannot split evenly\n");
                    return 0;
                }
            }
        }
        // prepare newBranch for valid clusters
        size_t parent = parentLinks[node];
        for (size_t i = 1; i < temp_centroids.size(); i++)
        {
            size_t centroid = temp_centroids[i];
            // fprintf(stderr, "*** %zu\n", centroid);
            if (isBranchNode[centroid])
            {
                moveParent(centroid, parent);
            }
            else if (clusters[i].size() > 0)
            {
                size_t t_parent = createParent(parent, insertionList);
                moveParent(centroid, t_parent);
                means[t_parent] = means[centroid];
                addSigToMatrix(parent, means[t_parent]);
            }
            else
            {
                // leaf & no NN
                moveParent(centroid, parent);
            }
        }

        // new centroid should have been move out at this stage
        for (size_t i = 1; i < temp_centroids.size(); i++)
        {
            size_t t_parent = parentLinks[temp_centroids[i]];
            for (size_t n = 0; n < clusters[i].size(); n++)
            {
                moveParent(clusters[i][n], t_parent);
            }
            updatePriority(t_parent);
        }

        updatePriority(node);
        if (isBranchNode[parent])
        {
            updatePriority(parent);
        }

        fprintf(stderr, ">> %f\n", priority[node]);
        return 1;
    }

    void prep()
    {
        for (size_t i = 0; i < childCounts[root]; i++)
        {
            size_t child = childLinks[root][i];
            // priority>split_threshold good enough?
            if (isBranchNode[child] && priority[child] > split_node_threshold)
            {
                size_t centroid = childLinks[child][0];
                if (isBranchNode[centroid])
                {
                    // dissolve this branch
                    childLinks[root][i] = centroid;
                    parentLinks[centroid] = root;

                    for (size_t j = 1; j < childCounts[child]; j++)
                    {
                        moveParent(childLinks[child][j], root, false);
                    }
                    clearNode(child); // doesn't do anything now, just clear memory; can insert back to insertionList
                }
                else
                {
                    for (size_t grandchild : childLinks[child])
                    {
                        if (isBranchNode[grandchild] && priority[grandchild] > split_node_threshold)
                        {
                            moveParent(grandchild, root);
                        }
                    }
                    updatePriority(child);
                }
            }
        }
    }

    size_t forcesplitBranch(size_t node, vector<size_t> &insertionList)
    {
        if (node != root)
        {
            return 0;
        }

        // printTreeJson(stderr);
        prep();
        // printTreeJson(stderr);

        fprintf(stderr, "force split node %zu, %f\n", node, priority[node]);
        vector<size_t> temp_centroids = {childLinks[node][0]};
        vector<vector<size_t>> clusters;
        vector<size_t> temp;
        clusters.push_back(temp);
        vector<size_t> t;
        clusters.push_back(t);

        double min_similarity = 1;
        size_t candidate = 0;
        seq_type mean0 = means[childLinks[node][0]];
        for (size_t i = childCounts[node] - 1; i > 0; i--)
        {
            double similarity = calcSimilarity(mean0, matrices[node][i]);
            // fprintf(stderr, "%zu, %f\n", childLinks[node][i], similarity);
            if (similarity < min_similarity)
            {
                min_similarity = similarity;
                candidate = childLinks[node][i];
            }
        }

        if (candidate == 0)
        {
            fprintf(stderr, "??something is wrong %zu, %f\n", node, priority[node]);
            return 0;
        }

        temp_centroids.push_back(candidate);
        // fprintf(stderr, ">??%zu, %zu\n", matrices[node].size(), childCounts[node]);

        // cluster the rest of the matrices, centroid will not be in clusters
        for (size_t n = 1; n < matrices[node].size(); n++)
        {
            size_t child = childLinks[node][n];
            if (child == candidate)
            {
                continue;
            }
            size_t dest = 0;
            double max_similarity = 0;
            for (size_t i = 0; i < temp_centroids.size(); i++)
            {
                double similarity = calcSimilarity(means[temp_centroids[i]], matrices[node][n]);
                // fprintf(stderr, "%zu, %f\n", temp_centroids[i], similarity);
                if (similarity > max_similarity)
                {
                    max_similarity = similarity;
                    dest = i;
                }
            }
            // fprintf(stderr, "--- %zu goes to %zu\n", child, dest);
            clusters[dest].push_back(child);
        }

        if (clusters[0].size() == 0 || clusters[1].size() == 0)
        {
            fprintf(stderr, "??something is wrong cannot split evenly\n");
            return 0;
        }

        // prepare newBranch for valid clusters

        // size_t parent = createParent(node, insertionList);
        // means[parent] = means[childLinks[node][0]];
        if (isBranchNode[node])
        {
            // promote the other clusters
            fprintf(stderr, "??promote\n");
            size_t parent = parentLinks[node];
            for (size_t i = 1; i < temp_centroids.size(); i++)
            {
                size_t centroid = temp_centroids[i];
                size_t t_parent = createParent(parent, insertionList);
                means[t_parent] = means[centroid];
                moveParent(centroid, t_parent);
                addSigToMatrix(parent, means[t_parent]);
            }
            // //?
            // for (size_t i = 1; i < temp_centroids.size(); i++)
            // {
            //     size_t centroid = temp_centroids[i];
            //     fprintf(stderr, "*** %zu\n", centroid);
            //     if (isBranchNode[centroid])
            //     {
            //         moveParent(centroid, parent);
            //     }
            //     else if (clusters[i].size() > 1)
            //     {
            //         size_t t_parent = createParent(parent, insertionList);
            //         moveParent(centroid, t_parent);
            //         means[t_parent] = means[centroid];
            //         addSigToMatrix(parent, means[t_parent]);
            //     }
            //     else
            //     {
            //         // leaf & no NN
            //         moveParent(centroid, parent);
            //     }
            // }
        }
        else
        {
            for (size_t i = 0; i < temp_centroids.size(); i++)
            {
                size_t centroid = temp_centroids[i];
                size_t t_parent = createParent(node, insertionList);
                means[t_parent] = means[centroid];
                moveParent(centroid, t_parent);
                addSigToMatrix(node, means[t_parent]);
            }
        }

        // new centroid should have been move out at this stage
        for (size_t i = 0; i < temp_centroids.size(); i++)
        {

            size_t centroid = temp_centroids[i];
            size_t t_parent = parentLinks[centroid];
            fprintf(stderr, "*** %zu\n", centroid);
            for (size_t n = 0; n < clusters[i].size(); n++)
            {
                moveParent(clusters[i][n], t_parent);
                fprintf(stderr, "--- %zu,%zu\n", clusters[i][n], t_parent);
            }
            updatePriority(t_parent);
        }

        updatePriority(node);

        // fprintf(stderr, ">?? done %zu, %zu\n", matrices[node].size(), childCounts[node]);

        // for (size_t i = 0; i < childCounts[root]; i++)
        // {
        //     if (matrices[root][i].first != means[childLinks[root][i]].first || matrices[root][i].second != means[childLinks[root][i]].second)
        //     {
        //         fprintf(stderr, "++ %zu, %zu\n", i, childLinks[root][i]);
        //         fprintf(stderr, "> %f, %f\n", matrices[root][i].first, matrices[root][i].second);
        //         fprintf(stderr, "* %f, %f\n\n", means[childLinks[root][i]].first, means[childLinks[root][i]].second);
        //     }
        // }

        return 1;
    }

    size_t forceSplitRoot(vector<size_t> &insertionList)
    {
        fprintf(stderr, ">> Force split root\n");
        // rootNodes.resize(1);
        size_t node = root;
        vector<size_t> temp_centroids = {childLinks[node][0]};
        vector<vector<size_t>> clusters;
        //?
        vector<size_t> temp;
        clusters.push_back(temp);
        vector<size_t> t;
        clusters.push_back(t);

        double min_similarity = 1;
        size_t candidate = 0;
        seq_type mean0 = means[childLinks[node][0]];
        for (size_t i = childCounts[node] - 1; i > 0; i--)
        {
            size_t child = childLinks[node][i];
            if (isAmbiNode[child])
            {
                continue;
            }
            double similarity = calcSimilarity(mean0, matrices[node][i]);
            // fprintf(stderr, "%zu, %f\n", child, similarity);
            if (similarity < min_similarity)
            {
                min_similarity = similarity;
                candidate = child;
            }
        }

        if (candidate == 0)
        {
            fprintf(stderr, "??something is wrong %zu, %f\n", node, priority[node]);
            return 0;
        }

        temp_centroids.push_back(candidate);

        fprintf(stderr, "??something is wrong %zu,%zu,%zu \n",  matrices[node].size(), childLinks[node].size(),childCounts[node]);

        for (size_t n = 0; n < matrices[node].size(); n++)
        {
            size_t child = childLinks[node][n];
            size_t dest = 0;
            double max_similarity = 0;
            for (size_t i = 0; i < temp_centroids.size(); i++)
            {
                double similarity = calcSimilarity(means[temp_centroids[i]], matrices[node][n]);
                // fprintf(stderr, "%zu, %f\n", temp_centroids[i], similarity);
                if (similarity > max_similarity)
                {
                    max_similarity = similarity;
                    dest = i;
                }
            }
            // fprintf(stderr, "--- %zu goes to %zu\n", child, dest);
            clusters[dest].push_back(child);
        }

        if (clusters[0].size() == 0 || clusters[1].size() == 0)
        {
            fprintf(stderr, "??something is wrong cannot split evenly\n");
            return 0;
        }

        for (size_t i = 0; i < temp_centroids.size(); i++)
        {
            size_t t_parent = createParent(node, insertionList);
            rootNodes.push_back(t_parent);
            for (size_t n = 0; n < clusters[i].size(); n++)
            {
                moveParent(clusters[i][n], t_parent);
                fprintf(stderr, "--- %zu,%zu\n", clusters[i][n], t_parent);
            }
            updatePriority(t_parent);
            updateNodeMean(t_parent);
            addSigToMatrix(node, means[t_parent]);
        }

        return 1;
    }

    inline size_t similarityStatus(size_t child, size_t rank, seq_type signature)
    {
        double offset = 1.2;
        // check the other children
        double local_stay_t = stay_threshold;
        // size_t rank = findLevel(child);
        for (size_t t = 0; t < rank; t++)
        {
            local_stay_t = local_stay_t + stay_threshold * offset;
        }
        double NN_t = local_stay_t + stay_threshold * offset;

        double local_split_t = split_threshold;
        double similarity = calcSimilarity(means[child], signature);

        if (similarity >= local_stay_t)
        {
            fprintf(stderr, "<%zu,%.2f>: %f stay\n", child, similarity, local_stay_t);
            return STAY_F;
        }
        else if (similarity <= local_split_t)
        {
            fprintf(stderr, "<%zu,%.2f>: %.2f miss\n", child, similarity, split_threshold);
            return SPLIT_F;
        }
        else
        {
            fprintf(stderr, "<%zu,%.2f>: NN\n", child, similarity);
            if (isBranchNode[child])
            {
                return NN_BRANCH_F;
            }
            else
            {
                return NN_LEAVE_F;
            }
        }
    }

    inline size_t checkRoot(seq_type signature, vector<size_t> &insertionList, vector<size_t> &mismatch, vector<size_t> &NN_leaves, vector<size_t> &NN_branches, vector<size_t> &stay, size_t node = 0)
    {
        size_t dest = 0;
        size_t best_rank = numeric_limits<size_t>::max();
        double offset = 1.2;
        // check the other children
        for (size_t child : childLinks[node])
        {
            // skip ambiNode
            if (isAmbiNode[child])
            {
                fprintf(stderr, "Ambi: %zu\n", child);
                continue;
            }
            // size_t rank = findLevel(child);
            size_t rank = 0;
            size_t status = similarityStatus(child, rank, signature);
            switch (status)
            {
            case STAY_F:
                stay.push_back(child);
                if (rank < best_rank)
                {
                    best_rank = rank;
                    dest = child;
                }
                break;

            case SPLIT_F:
                mismatch.push_back(child);
                break;
            case NN_LEAVE_F:
                NN_leaves.push_back(child);
                break;
            case NN_BRANCH_F:
                NN_branches.push_back(child);
                break;
            default:
                break;
            }
        }
        return dest;
    }

    // put NN leaves and branches under the same supercluster
    // create an ambiNode to store seqs that are NN to some or all the NN
    // all children in this supercluster except for the ambiNode should be "split"
    inline size_t superCluster(seq_type signature, vector<size_t> &insertionList, size_t idx, size_t node, vector<size_t> NN_leaves, vector<size_t> NN_branches)
    {
        size_t t_parent = createParent(node, insertionList);
        fprintf(stderr, "*** super\n");
        for (size_t child : NN_leaves)
        {
            moveParent(child, t_parent);
        }

        for (size_t child : NN_branches)
        {
            moveParent(child, t_parent);
        }

        addSigToMatrix(node, signature); // this should be updated later

        size_t dest = createAmbiNode(signature, insertionList, t_parent, idx);
        updateNodeMean(t_parent);
        updateNodeMean(node);
        if (ambiLinks[node].size() > 0)
        {
            // shouldn't happen at root
            fprintf(stderr, "check ambi\n");
        }
        return dest;
    }

    inline size_t tt_root(seq_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        // size_t node = root;
        vector<size_t> mismatch;
        vector<size_t> NN_leaves;
        vector<size_t> NN_branches;
        vector<size_t> stay;

        size_t dest = checkRoot(signature, insertionList, mismatch, NN_leaves, NN_branches, stay, node);

        if (dest != 0)
        {
            fprintf(stderr, "#stay in %zu from %zu\n", dest, stay.size());
            return stayNode(signature, insertionList, idx, dest);
        }

        if (mismatch.size() == childCounts[node] - ambiLinks[node].size())
        {
            fprintf(stderr, "#mismatch\n");
            return createNode(signature, insertionList, node, idx);
        }

        size_t NN_total = NN_leaves.size() + NN_branches.size();

        if (NN_branches.size() == 0)
        {
            // if (NN_leaves.size() > 1)
            {
                fprintf(stderr, "NN leaves\n");
                return superCluster(signature, insertionList, idx, node, NN_leaves, NN_branches);

                // dbgPrintSignatureIdx(stderr, means[t_parent]);
                // for (size_t c : childLinks[t_parent])
                // {
                //     dbgPrintSignatureIdx(stderr, means[c]);
                // }
            }
        }
        else if (NN_branches.size() == 1)
        {
            fprintf(stderr, "#match ONE branch > ");
            dest = NN_branches[0];
            if (NN_total == 1)
            {
                // no NN leaves
                fprintf(stderr, "?no NN leaves\n");
                return tt_branch(signature, insertionList, idx, dest);
            }
            else
            {
                fprintf(stderr, "?with NN leaves\n");
                return superCluster(signature, insertionList, idx, node, NN_leaves, NN_branches);
            }
        }
        else
        {
            fprintf(stderr, "#match multiple branches > ");
            double max_sim = 0;
            for (size_t c : NN_branches)
            {
                double similarity = calcSimilarity(means[c], signature);
                if (c > max_sim)
                {
                    max_sim = similarity;
                    dest = c;
                }
            }

            if (NN_leaves.size() == 0)
            {
                // continue with nearest, might change
                fprintf(stderr, "?without leaves\n");
                return tt_branch(signature, insertionList, idx, dest);
            }
            else
            {
                // if there is any leaf has higher similarity than all the branches,
                // ignore NN branches, supercluster on the NN leaves only
                // or else, continue with the nearest branch and ignore leaves
                // might change
                fprintf(stderr, "?with leaves\n");
                for (size_t c : NN_leaves)
                {
                    double similarity = calcSimilarity(means[c], signature);
                    if (c > max_sim)
                    {
                        max_sim = similarity;
                        dest = c;
                        NN_branches.clear();
                        break;
                    }
                }

                if (isBranchNode[dest])
                {
                    return tt_branch(signature, insertionList, idx, dest);
                }
                else
                {
                    return superCluster(signature, insertionList, idx, node, NN_leaves, NN_branches);
                }
            }
        }
        // if (NN_total == 1)
        // {
        //     fprintf(stderr, "??NN with one without stay\n");
        //     return createNode(signature, insertionList, node, idx);
        // }
        // else if (NN_total > 1)
        // {
        //     fprintf(stderr, "NN with multiple without stay\n");
        //     size_t t_parent = createParent(node, insertionList);
        //     dest = createNode(signature, insertionList, t_parent, idx);
        //     // means[t_parent] = signature;
        //     // addSigToMatrix(node, means[t_parent]);
        //     for (size_t child : NN_leaves)
        //     {
        //         moveParent(child, t_parent);
        //     }

        //     // means[t_parent] = createMeanSig(matrices[t_parent]);
        //     // updatePriority(t_parent);
        //     updateNodeMean(t_parent);
        //     return dest;
        // }

        return createNode(signature, insertionList, node, idx);
    }

    inline size_t tt_branch(seq_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        if (find(rootNodes.begin(), rootNodes.end(), node) != rootNodes.end())
        {
            fprintf(stderr, "is subtree %zu\n", node);
            return tt_root(signature, insertionList, idx, node);
        }

        size_t current_childCount = childCounts[node];
        vector<size_t> mismatch;
        vector<size_t> NN_leaves;
        vector<size_t> NN_branches;
        vector<size_t> stay;

        size_t dest = checkRoot(signature, insertionList, mismatch, NN_leaves, NN_branches, stay, node);
        fprintf(stderr, "tt_branch %zu\n", node);

        if (dest != 0)
        {
            fprintf(stderr, "#b- stay in %zu from %zu\n", dest, stay.size());
            return stayNode(signature, insertionList, idx, dest);
        }

        if (mismatch.size() == childCounts[node] - ambiLinks[node].size())
        {
            //? this might be a wrong supercluster, dissolve this branch and check from parent again
            fprintf(stderr, "//?#b- mismatch all\n");
            dest = createNode(signature, insertionList, node, idx);
            updateNodeMean(node);
            return dest;
        }
        else if (mismatch.size() > 0)
        {
            fprintf(stderr, "//?#b- mismatch, promote\n");
        }

        size_t NN_total = NN_leaves.size() + NN_branches.size();
        if (mismatch.size() == 0)
        {
            // NN with everything
            fprintf(stderr, "#??NN with all\n");
            if (ambiLinks[node].size() == 0)
            {
                return createAmbiNode(signature, insertionList, node, idx);
            }
            else
            {
                return stayNode(signature, insertionList, idx, ambiLinks[node][0]);
            }
        }
        else
        {
            // some mismatches
            fprintf(stderr, "#NN with some\n");
            if (ambiLinks[node].size() == 0)
            {
                return superCluster(signature, insertionList, idx, node, NN_leaves, NN_branches);
            }
            else
            {
                return stayNode(signature, insertionList, idx, ambiLinks[node][0]);
            }
        }
        fprintf(stderr, "//?#b- wrong\n");
        return createNode(signature, insertionList, node, idx);
    }

    inline size_t tt(seq_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        if (node == root)
        {
            return tt_root(signature, insertionList, idx);
        }
        else if (isBranchNode[node])
        {
            return tt_branch(signature, insertionList, idx, node);
        }
        else
        {
            fprintf(stderr, "Something went wrong\n");
        }
    }

    inline size_t first_insert(seq_type signature, vector<size_t> &insertionList, size_t idx, size_t parent = 0)
    {
        return createNode(signature, insertionList, root, idx);
    }

    inline size_t insert(seq_type signature, vector<size_t> &insertionList, size_t idx)
    {
        size_t node = tt(signature, insertionList, idx);

        fprintf(stderr, "inserted %zu at %zu\n\n", idx, node);

        if (childCounts[root] > 5)
        {
            printTreeJson(stderr);
            forceSplitRoot(insertionList);
            printTreeJson(stderr);
        }
        return node;
    }

    inline size_t search(seq_type signature, size_t idx = 0)
    {
        size_t node = root;
        size_t best_child = node;
        double best_similarity = 0;

        do
        {
            size_t local_best_child = node;
            double local_best_similarity = 0;

            for (size_t i = 0; i < childCounts[node]; i++)
            {
                size_t child = childLinks[node][i];
                if (isAmbiNode[child])
                {
                    continue;
                }

                double similarity = calcSimilarity(means[child], signature);
                // fprintf(stderr, " <%zu,%.2f> ", child, similarity);

                // found same, move on with the next seq
                if (similarity >= (stay_threshold))
                {
                    // while (isBranchNode[child])
                    // {
                    //     child = childLinks[child][0];
                    // }
                    return child;
                }

                if (similarity >= local_best_similarity)
                {
                    local_best_similarity = similarity;
                    local_best_child = child;
                }
            }
            if (local_best_similarity >= best_similarity)
            {
                best_similarity = local_best_similarity;
                best_child = local_best_child;
            }

            node = local_best_child;

        } while (isBranchNode[node]);

        // seqIDs[best_child].push_back(idx);
        return best_child;
    }

    // cluster means still remain
    void prepReinsert(size_t node = 0)
    {
        if (!isBranchNode[node])
        {
            seqIDs[node].clear();
        }

        for (size_t child : childLinks[node])
        {
            prepReinsert(child);
        }
    }

    size_t reinsert(seq_type signature, size_t idx)
    {
        size_t node = search(signature);
        seqIDs[node].push_back(idx);
        addSigToMatrix(node, signature);
        return node;
    }

    void getEmptyLeaves(vector<size_t> &leaves, size_t node = 0)
    {
        for (size_t child : childLinks[node])
        {
            if (isBranchNode[child])
            {
                getEmptyLeaves(leaves, child);
            }
            else if (isAmbiNode[child] || seqIDs[child].size() == 0)
            {
                leaves.push_back(child);
            }
            else
            {
                // updatePriority(child);
                updateNodeMean(node);
            }
        }
    }

    bool redundant(size_t parent)
    {
        fprintf(stderr, "redundant %zu\n", parent);
        if (childCounts[parent] == 0)
        {
            return true;
        }
        else if (childCounts[parent] == 1)
        {
            size_t child = childLinks[parent][0];
            if (isBranchNode[child])
            {
                return childCounts[child] == 0;
            }
        }
        return false;
    }

    void dropBranch(size_t node = 0)
    {
        if (isBranchNode[node])
        {
            size_t centroid = node;
            size_t parent = parentLinks[node];
            size_t idx = getNodeIdx(node);
            while (isBranchNode[centroid] && childCounts[centroid] == 1)
            {
                size_t temp = centroid;
                centroid = childLinks[centroid][0];
                clearNode(temp);
            }

            if (centroid != node)
            {
                childLinks[parent][idx] = centroid;
                parentLinks[centroid] = parent;
            }

            for (size_t child : childLinks[centroid])
            {
                dropBranch(child);
            }
        }
    }

    void trim(size_t last_idx)
    {
        std::set<size_t> branches;
        fprintf(stderr, "\nTrimming\n");

        vector<size_t> leaves;

        getEmptyLeaves(leaves);

        // do leaves
        for (size_t i : leaves)
        {
            if (seqIDs[i].size() == 0)
            {
                size_t parent = parentLinks[i];
                deleteNode(i);
                if (childCounts[parent] == 0)
                {
                    branches.insert(parent);
                }
            }
        }

        // do parents
        for (auto &b : branches)
        {
            size_t node = b;
            size_t parent = parentLinks[node];

            while (redundant(node))
            {
                size_t temp = parent;
                parent = parentLinks[parent];
                deleteNode(node);
                clearNode(node);
                node = temp;
            }
            // updatePriority(node);
            updateNodeMean(node);
        }

        // remove unitig for non-empty leaves
        for (size_t child : childLinks[root])
        {
            dropBranch(child);
        }
    }

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