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
#include <bitset>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include "read.hpp"
#include "distance.hpp"
#include <stdarg.h>
#include <stdio.h>

using namespace std;

// typedef unsigned char cell_type;
bloom_parameters parameters;
size_t partree_capacity = 100;
size_t singleton = 2;
size_t tree_order = 5;
bool print_ = false;

void printMsg(const char *format, ...)
{
    if (print_)
    {
        va_list args;
        va_start(args, format);
        vfprintf(stderr, format, args);
        va_end(args);
    }
}

seq_type createMeanSig(const vector<seq_type> clusterSigs)
{
    seq_type meanSig;
    // find the smallest windows count
    size_t winNum = clusterSigs[0].size();
    for (seq_type matrix : clusterSigs)
    {
        if (matrix.size() < winNum)
        {
            winNum = matrix.size();
        }
    }
    vector<vector<size_t>> counters;
    for (size_t w = 0; w < winNum; w++)
    {
        vector<size_t> counter(signatureSize * bits_per_char);
        counters.push_back(counter);

        vector<cell_type> temp(signatureSize * bits_per_char);
        meanSig.push_back(temp);
    }

    for (seq_type signatureData : clusterSigs)
    {
        for (size_t w = 0; w < winNum; w++)
        {
            for (size_t i = 0; i < signatureSize; i++)
            {
                for (int n = 0; n < bits_per_char; n++)
                {
                    if ((signatureData[w][i] >> n) & 1)
                    {
                        counters[w][i * bits_per_char + n]++;
                    }
                }
            }
        }
    }
    for (size_t w = 0; w < winNum; w++)
    {
        vector<size_t> counter = counters[w];
        for (int i = 0; i < counter.size(); i++)
        {
            // make it upperbound
            if (counter[i] >= (clusterSigs.size() + 1) / 2)
            {
                meanSig[w][i / bits_per_char] |= (cell_type)1 << (i % bits_per_char);
            }
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
    vector<int> isRootNode;            // n entries, is this root of subtree
    vector<int> isSuperNode;           // n entries, is this root of subtree
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
                isRootNode.reserve(capacity);
                isSuperNode.reserve(capacity);
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
        isRootNode[root] = 1;
        isSuperNode[root] = 0;
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
        fprintf(stream, "\"subtree\":\"%zu\",", isRootNode[tnode]);
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

    void calcNodeDistance(FILE *stream, size_t tnode, size_t idx, seq_type sig)
    {
        double distance = calcHD(means[tnode], sig);
        fprintf(stream, "%zu,%zu,%.4f", tnode, idx, distance);
    }

    void printNodeDistance(FILE *stream, const vector<seq_type> &seqs, vector<size_t> &clusters)
    {
        fprintf(stream, "seq_id,clu,HD\n");
        for (size_t i = 0; i < clusters.size(); i++)
        {
            size_t tnode = clusters[i];
            double distance = calcHD(means[tnode], seqs[i]);
            fprintf(stream, "%zu,%zu,%.4f\n", i, tnode, distance);
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
        for (size_t i = 0; i < childLinks[parent].size(); i++)
        {
            if (childLinks[parent][i] == node)
            {
                idx = i;
                break;
            }
        }
        return idx;
    }

    seq_type createRandomSig(vector<size_t> children)
    {
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine rng(seed);
        uniform_int_distribution<size_t> dist(0, children.size() - 1);

        seq_type randomSig; // = means[children[dist(rng)]];
        // find the smallest windows count
        size_t winNum = means[children[0]].size();
        for (size_t child : children)
        {
            seq_type matrix = means[child];
            if (matrix.size() < winNum)
            {
                winNum = matrix.size();
            }
        }

        randomSig.resize(winNum);

        vector<vector<size_t>> counters;
        for (size_t w = 0; w < winNum; w++)
        {
            size_t s = dist(rng);
            randomSig[w] = means[children[s]][w];
        }

        return randomSig;
    }

    // delete a child from its parent, need to format the child separately
    inline void deleteNode(size_t node)
    {
        // printMsg("deleting %zu\n", node);
        size_t parent = parentLinks[node];

        int idx = getNodeIdx(node);

        if (idx == -1)
        {
            printMsg("ERROR deleting %zu from %zu!!\n", node, parent);
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

    inline vector<seq_type> getNonAmbiMatrix(size_t node)
    {
        vector<seq_type> temp_matrix;
        for (size_t c = 0; c < childCounts[node]; c++)
        {
            size_t child = childLinks[node][c];
            // skip ambi
            if (!isAmbiNode[child])
            {
                temp_matrix.push_back(means[child]);
                // temp_matrix.push_back(matrices[node][c]);
            }
        }

        return temp_matrix;
    }

    // union mean of children
    inline void updateNodeMean(size_t node)
    {
        if (isBranchNode[node] || isRootNode[node])
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
            means[node] = meanSig;

            // means[node] = createMeanSig(getNonAmbiMatrix(node));
        }
        else
        {
            means[node] = createMeanSig(matrices[node]);
        }

        // size_t parent = parentLinks[node];
        // matrices[parent][getNodeIdx(node)] = means[node];
    }

    inline void updateParentMean(size_t node)
    {
        while (node != root)
        {
            int idx = getNodeIdx(node);
            updateNodeMean(node);
            updatePriority(node);
            size_t parent = parentLinks[node];

            if (idx == -1)
            {
                printMsg("ERROR updating mean %zu from %zu!!\n", node, parent);
            }
            else
            {
                matrices[parent][idx] = means[node];
            }
            node = parent;
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

        // printMsg(">>%f,%f\n", meanSig.first, meanSig.second);
        for (seq_type sig : matrices[node])
        {
            double similarity = calcSimilarity(meanSig, sig);

            // if (isBranchNode[node])
            // {
            //     printMsg("%zu,%f,%f,%f\n", childLinks[node][i], sig.first, sig.second, similarity);
            // }
            // else
            // {
            //     printMsg("%zu,%f,%f,%f\n", seqIDs[node][i], sig.first, sig.second, similarity);
            // }
            if (similarity < min_similarity)
            {
                min_similarity = similarity;
            }
            i++;
        }

        return -min_similarity;
    }

    // RMSD
    double test(size_t node)
    {
        vector<seq_type> temp_matrix = getNonAmbiMatrix(node);
        if (temp_matrix.size() <= 1)
        {
            return 0;
        }

        double sumSquaredistance = 0;

        seq_type meanSig = createMeanSig(temp_matrix);
        dbgPrintSignatureIdx(stderr, meanSig);
        auto children = childLinks[node];
        size_t i = 0;

        for (seq_type signature : temp_matrix)
        {
            double distance = calcSimilarity(meanSig, signature);
            dbgPrintSignatureIdx(stderr, signature);
            while (!isAmbiNode[children[i]])
            {
                i++;
            }

            printMsg(">>%zu,%.2f\n", children[i], distance);
            sumSquaredistance += distance * distance;
        }

        return sqrt(sumSquaredistance / temp_matrix.size());
    }

    // RMSD
    double calcDistortion(size_t node)
    {
        vector<seq_type> temp_matrix = getNonAmbiMatrix(node);
        if (temp_matrix.size() <= 1)
        {
            return 0;
        }

        double sumSquaredistance = 0;
        // vector<cell_type> meanSig = getMinimiseSet(createMeanSig(temp_matrix))[0];

        // for (seq_type signature : temp_matrix)
        // {
        //     vector<cell_type> signatureData = getMinimiseSet(signature)[0];
        //     double distance = calcJaccardGlobal_cell(meanSig, signatureData);
        //     sumSquaredistance += distance * distance;
        // }

        seq_type meanSig = createMeanSig(temp_matrix);

        for (seq_type signature : temp_matrix)
        {
            double distance = calcSimilarity(meanSig, signature);
            sumSquaredistance += distance * distance;
        }

        return sqrt(sumSquaredistance / temp_matrix.size());
    }

    // RMSD
    double calcDistortionHD(size_t node)
    {

        vector<seq_type> temp_matrix = getNonAmbiMatrix(node);
        if (temp_matrix.size() <= 1)
        {
            return 0;
        }

        seq_type meanSig = createMeanSig(temp_matrix);

        double sumSquareHD = 0;

        for (seq_type signature : temp_matrix)
        {
            double similarity = calcHD(meanSig, signature);
            sumSquareHD += similarity * similarity;
        }
        return sqrt(sumSquareHD / temp_matrix.size());
    }

    inline double calcAvgSim(size_t node)
    {
        vector<seq_type> temp_matrix = getNonAmbiMatrix(node);
        if (temp_matrix.size() <= 1)
        {
            return 0;
        }

        double sumDistance = 0;
        seq_type meanSig = createMeanSig(temp_matrix);

        for (seq_type signature : temp_matrix)
        {
            double distance = calcSimilarity(meanSig, signature);
            sumDistance += distance;
        }

        return sumDistance / temp_matrix.size();
    }

    // union mean of children
    inline void updatePriority(size_t node)
    {
        // priority[node] = calcDistortion(node);
        priority[node] = calcAvgSim(node);

        // priority[node] = calcDistortionHD(node);
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

        printMsg(">%zu; priority: %.1f,%.1f\n", parent, priority[parent], priority[node]);

        size_t grandparent = parentLinks[parent];
        vector<size_t> &siblings = childLinks[parent];
        printMsg("\nrotate: %zu,%zu,%zu\n", grandparent, parent, node);

        // remove(siblings.begin(), siblings.end(), node);
        siblings.erase(remove(siblings.begin(), siblings.end(), node), siblings.end());

        vector<size_t> &uncles = childLinks[grandparent];
        uncles.erase(remove(uncles.begin(), uncles.end(), parent), uncles.end());

        childLinks[grandparent].push_back(node);
        childCounts[grandparent] += siblings.size();
        for (size_t sibling : siblings)
        {
            printMsg("\nSibling: %zu\n", sibling);
            childLinks[grandparent].push_back(sibling);
            parentLinks[sibling] = grandparent;
        }

        parentLinks[node] = grandparent;
        childCounts[node]++;

        for (size_t grand : childLinks[node])
        {
            printMsg("%zu,", grand);
        }
        printMsg("\n");
        childLinks[node].push_back(parent);
        for (size_t grand : childLinks[node])
        {
            printMsg("%zu,", grand);
        }
        printMsg("\n");

        parentLinks[parent] = node;

        int temp = isBranchNode[parent];
        isBranchNode[parent] = isBranchNode[node];
        isBranchNode[node] = temp;

        printMsg("finish: %zu\n", grandparent);

        for (size_t child : childLinks[grandparent])
        {
            printMsg("%zu>", child);
            for (size_t grand : childLinks[child])
            {
                printMsg("%zu,", grand);
            }
            printMsg("\n");
        }

        return parent;
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
            // updatePriority(node);
            return node;
        }
    }

    // check against all ambiNode
    // create "dest" (non-ambi) and move all matrices that stay with input sig from the ambi nodes
    // no nothing stay, delete dest AND
    // create new ambi if there is at least one split in all ambi nodes
    // else store in the first ambiNode
    inline size_t stayAmbi(seq_type signature, vector<size_t> &insertionList, size_t idx, size_t node)
    {
        return stayNode(signature, insertionList, idx, ambiLinks[node][0]);

        // // size_t ambi = ambiLinks[node][0];
        // size_t newAmbi = 0;

        // vector<seq_type> staySigs;
        // vector<size_t> staySeqIDs;

        // for (size_t ambi : ambiLinks[node])
        // {
        //     printMsg("*** stayAmbi %zu\n", ambi);
        //     bool newAmbi_F = false;
        //     vector<seq_type> tempSigs;
        //     vector<size_t> tempSeqIDs;
        //     size_t i = 0;

        //     for (seq_type matrix : matrices[ambi])
        //     {
        //         size_t status = similarityStatus(matrix, signature);
        //         if (status == STAY_F)
        //         {
        //             staySigs.push_back(matrix);
        //             staySeqIDs.push_back(seqIDs[ambi][i]);
        //         }
        //         else
        //         {
        //             if (status == SPLIT_F)
        //             {
        //                 newAmbi_F = true;
        //             }
        //             tempSigs.push_back(matrix);
        //             tempSeqIDs.push_back(seqIDs[ambi][i]);
        //         }
        //         i++;
        //     }

        //     if (tempSeqIDs.size() != seqIDs[ambi].size())
        //     {
        //         matrices[ambi] = tempSigs;
        //         seqIDs[ambi] = tempSeqIDs;
        //     }

        //     if (newAmbi_F)
        //     {
        //         newAmbi++;
        //     }
        // }

        // for (size_t ambi : ambiLinks[node])
        // {
        //     if (seqIDs[ambi].size() == 0)
        //     {
        //         deleteNode(ambi);
        //     }
        // }

        // if (staySeqIDs.size() == 0)
        // {
        //     if (newAmbi == ambiLinks[node].size())
        //     {

        //         printMsg("***???\n");
        //         return createAmbiNode(signature, insertionList, node, idx);
        //     }
        //     else
        //     {
        //         return stayNode(signature, insertionList, idx, ambiLinks[node][0]);
        //     }
        // }
        // else
        // {
        //     printMsg("***---\n");
        //     size_t dest = createNode(signature, insertionList, node, idx);
        //     matrices[dest] = staySigs;
        //     seqIDs[dest] = staySeqIDs;
        //     matrices[dest].push_back(signature);
        //     seqIDs[dest].push_back(idx);
        //     means[dest] = matrices[dest][0];

        //     updateParentMean(node);
        //     return dest;
        // }
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
        isRootNode[t_parent] = 0;

        // update node
        childCounts[node]++;
        childLinks[node].push_back(t_parent);

        printMsg("create t_parent %zu\n", t_parent);
        return t_parent;
    }

    inline size_t createNode(seq_type signature, vector<size_t> &insertionList, size_t t_parent, size_t idx)
    {
        size_t new_node = getNewNodeIdx(insertionList);
        parentLinks[new_node] = t_parent;
        addSigToMatrix(new_node, signature);
        seqIDs[new_node].push_back(idx);
        means[new_node] = signature;
        isRootNode[new_node] = 0;

        // then add new node to t_parent
        childCounts[t_parent]++;
        childLinks[t_parent].push_back(new_node);
        addSigToMatrix(t_parent, means[new_node]);

        if (isBranchNode[t_parent])
        {
            updatePriority(t_parent);
        }

        printMsg("create new node %zu\n", new_node);
        return new_node;
    }

    // Ambi Node cannot be the first child or else update node mean will have problem
    inline size_t createAmbiNode(seq_type signature, vector<size_t> &insertionList, size_t node, size_t idx)
    {
        size_t dest = createNode(signature, insertionList, node, idx);

        printMsg(">>Ambi\n");
        // size_t dest = getNewNodeIdx(insertionList);
        // parentLinks[dest] = node;
        // childLinks[node].push_back(dest);
        // seqIDs[dest].push_back(idx);
        ambiLinks[node].push_back(dest);
        isAmbiNode[dest] = 1;
        return dest;
    }

    // size_t findFurthest(size_t node, seq_type mean0)
    // {
    //     double min_similarity = 1;
    //     size_t candidate = 0;
    //     // seq_type mean0 = means[childLinks[node][0]];
    //     for (size_t i = 0; i < childCounts[node]; i++)
    //     {
    //         size_t child = childLinks[node][i];
    //         if (isAmbiNode[child])
    //         {
    //             continue;
    //         }
    //         double similarity = calcSimilarity(mean0, matrices[node][i]);
    //         // double similarity = calcSimilarity(mean0, means[child]);

    //         if (similarity < min_similarity)
    //         {
    //             min_similarity = similarity;
    //             candidate = child;
    //         }
    //     }

    //     if (candidate == 0)
    //     {
    //         printMsg("??cannot find candidate %zu, %f\n", node, priority[node]);
    //     }
    //     return candidate;
    // }

    // // this code is not working with unpacked
    // size_t forceSplitRoot(vector<size_t> &insertionList, bool unpack = false, size_t node = 0)
    // {
    //     // printTreeJson(stderr);
    //     // rootNodes.resize(1);
    //     // size_t node = root;

    //     // unpack previous subtree
    //     unpack = false;
    //     if (unpack)
    //     {
    //         printMsg(">> Unpacking...\n");
    //         vector<size_t> tempChildLinks = childLinks[node];
    //         for (size_t subtree : tempChildLinks)
    //         {
    //             if (isRootNode[subtree])
    //             {
    //                 for (size_t child : childLinks[subtree])
    //                 {
    //                     moveParent(child, node, false);
    //                 }
    //                 deleteNode(subtree);
    //             }
    //         }
    //     }

    //     //? assume the first child is never ambi
    //     vector<size_t> temp_centroids;
    //     vector<vector<size_t>> clusters;
    //     //?
    //     vector<size_t> temp;
    //     clusters.push_back(temp);
    //     vector<size_t> t;
    //     clusters.push_back(t);

    //     // size_t candidate = findFurthest(node, createMeanSig(getNonAmbiMatrix(node)));

    //     size_t candidate = findFurthest(node, means[node]);
    //     temp_centroids.push_back(candidate);
    //     temp_centroids.push_back(findFurthest(node, means[candidate]));

    //     // printMsg("??something is wrong %zu,%zu,%zu \n", matrices[node].size(), childLinks[node].size(), childCounts[node]);

    //     for (size_t n = 0; n < matrices[node].size(); n++)
    //     {
    //         size_t child = childLinks[node][n];
    //         size_t dest = 0;
    //         double max_similarity = 0;
    //         for (size_t i = 0; i < temp_centroids.size(); i++)
    //         {
    //             // double similarity = calcSimilarity(means[temp_centroids[i]], matrices[node][n]);
    //             double similarity = calcSimilarity(means[temp_centroids[i]], means[child]);

    //             printMsg("%zu, %f\n", temp_centroids[i], similarity);
    //             if (similarity > max_similarity)
    //             {
    //                 max_similarity = similarity;
    //                 dest = i;
    //             }
    //         }
    //         printMsg("--- %zu goes to %zu\n", child, dest);
    //         clusters[dest].push_back(child);
    //     }

    //     if (clusters[0].size() <= 1 || clusters[1].size() <= 1)
    //     {
    //         printMsg("??something is wrong cannot split evenly\n");
    //         return 0;
    //     }

    //     printMsg(">> Force split root %zu\n", node);

    //     for (size_t i = 0; i < temp_centroids.size(); i++)
    //     {
    //         size_t t_parent = createParent(node, insertionList);
    //         isRootNode[t_parent] = 1;
    //         for (size_t n = 0; n < clusters[i].size(); n++)
    //         {
    //             moveParent(clusters[i][n], t_parent);
    //             // printMsg("--- %zu,%zu\n", clusters[i][n], t_parent);
    //         }
    //         updateNodeMean(t_parent);
    //         updatePriority(t_parent);
    //         addSigToMatrix(node, means[t_parent]);
    //     }
    //     updateParentMean(node);
    //     // printTreeJson(stderr);

    //     return 1;
    // }

    vector<size_t> findFurthest(vector<size_t> children)
    {
        double min_similarity = 1;
        size_t a, b;

        for (size_t i = 0; i < children.size(); i++)
        {
            for (size_t j = i + 1; j < children.size(); j++)
            {
                double similarity = calcSimilarity(means[children[i]], means[children[j]]);
                printMsg("## %zu, %zu, %f\n", children[i], children[j], similarity);

                if (similarity < min_similarity)
                {
                    min_similarity = similarity;
                    a = children[i];
                    b = children[j];
                }
            }
        }

        return {a, b};
    }

    size_t forceSplitRoot(vector<size_t> &insertionList, bool unpack = false, size_t node = 0)
    {
        size_t clusterCount = 2;
        // unpack previous subtree
        vector<size_t> children;
        // vector<size_t> grandchildren;
        vector<size_t> to_remove_children;

        if (unpack)
        {
            printMsg(">> Unpacking...\n");
            for (size_t subtree : childLinks[node])
            {
                if (isRootNode[subtree])
                {
                    for (size_t child : childLinks[subtree])
                    {
                        // moveParent(child, node, false);
                        // grandchildren.push_back(child);
                        children.push_back(child);
                    }
                    // deleteNode(subtree);
                    to_remove_children.push_back(subtree);
                }
                else if (!isAmbiNode[subtree])
                {
                    children.push_back(subtree);
                }
            }
            // children.insert(children.end(), grandchildren.begin(), grandchildren.end());
        }
        else
        {
            children = childLinks[node];
        }

        vector<size_t> clusters(children.size());
        vector<seq_type> temp_centroids(clusterCount);
        for (size_t i = 0; i < clusterCount; i++)
        {
            temp_centroids[i] = createRandomSig(children);
        }

        vector<size_t> clustersSize(clusterCount);

        for (size_t n = 0; n < children.size(); n++)
        {
            size_t child = children[n];
            size_t dest = 0;
            double max_similarity = 0;
            for (size_t i = 0; i < temp_centroids.size(); i++)
            {
                // double similarity = calcSimilarity(means[temp_centroids[i]], matrices[node][n]);
                double similarity = calcSimilarity(temp_centroids[i], means[child]);
                printMsg("%zu, %f\n", i, similarity);
                if (similarity > max_similarity)
                {
                    max_similarity = similarity;
                    dest = i;
                }
            }
            printMsg("--- %zu goes to %zu\n", child, dest);
            clusters[n] = dest;
            clustersSize[dest]++;
        }

        for (size_t size : clustersSize)
        {
            if (size <= 1)
            {
                printMsg("??something is wrong cannot split evenly\n");
                return 0;
            }
        }

        printTreeJson(stderr);
        printMsg(">> Force split root %zu\n", node);

        // reuse clusterSize to store the new t_parents
        for (size_t i = 0; i < temp_centroids.size(); i++)
        {
            size_t t_parent = createParent(node, insertionList);
            isRootNode[t_parent] = 1;
            clustersSize[i] = t_parent;
        }

        for (size_t n = 0; n < clusters.size(); n++)
        {
            moveParent(children[n], clustersSize[clusters[n]]);
            // printMsg("--- %zu,%zu\n", clusters[i][n], t_parent);
        }

        for (size_t child : to_remove_children)
        {
            deleteNode(child);
        }

        for (size_t t_parent : clustersSize)
        {
            updateNodeMean(t_parent);
            updatePriority(t_parent);
            addSigToMatrix(node, means[t_parent]);
        }
        updateParentMean(node);
        printTreeJson(stderr);

        return 1;
    }

    size_t addSubtree(vector<size_t> &insertionList, bool unpack = false, size_t node = 0)
    {
        size_t clusterCount = 2;
        // unpack previous subtree
        vector<size_t> children;
        vector<size_t> to_remove_children;

        if (unpack)
        {
            printMsg(">> Unpacking...\n");
            for (size_t subtree : childLinks[node])
            {
                if (isRootNode[subtree])
                {
                    for (size_t child : childLinks[subtree])
                    {
                        children.push_back(child);
                    }
                    to_remove_children.push_back(subtree);
                }
                else if (!isAmbiNode[subtree])
                {
                    children.push_back(subtree);
                }
            }
        }
        else
        {
            children = childLinks[node];
        }
        printTreeJson(stderr);

        vector<vector<size_t>> clusters(clusterCount);
        vector<seq_type> temp_centroids(clusterCount);
        temp_centroids[0] = means[node];

        // find another sibling centroid if any subtree has no child
        while (true)
        {
            temp_centroids[1] = createRandomSig(children);
            for (size_t child : children)
            {
                size_t dest = 0;
                double max_similarity = 0;
                for (size_t i = 0; i < temp_centroids.size(); i++)
                {
                    double similarity = calcSimilarity(temp_centroids[i], means[child]);
                    printMsg("%zu, %f\n", i, similarity);
                    if (similarity > max_similarity)
                    {
                        max_similarity = similarity;
                        dest = i;
                    }
                }
                printMsg("--- %zu goes to %zu\n", child, dest);
                clusters[dest].push_back(child);
            }
            if (clusters[0].size() > 1 && clusters[1].size() > 1)
            {
                break;
            }
            clusters.clear();
            clusters.resize(2);
        }

        printMsg(">> Add subtree %zu\n", node);

        size_t grandparent = parentLinks[node];
        size_t t_parent = createParent(grandparent, insertionList);
        isRootNode[t_parent] = 1;

        for (size_t child : clusters[1])
        {
            moveParent(child, t_parent);
        }

        // if unpacking
        for (size_t child : to_remove_children)
        {
            deleteNode(child);
        }

        updateNodeMean(node);
        updatePriority(node);

        updateNodeMean(t_parent);
        updatePriority(t_parent);
        addSigToMatrix(grandparent, means[t_parent]);

        updateParentMean(grandparent);
        printTreeJson(stderr);

        return 1;
    }

    inline size_t similarityStatus(seq_type sig1, seq_type sig2, double local_stay_t, double offset = 0)
    {
        double local_split_t = split_threshold * offset;
        double similarity = calcSimilarity(sig1, sig2);

        printMsg("%.2f", similarity);

        if (similarity >= local_stay_t)
        {
            printMsg(": STAY> b\n");
            return STAY_F;
        }
        else if (similarity <= local_split_t)
        {
            printMsg(": SPLIT> b\n");
            return SPLIT_F;
        }
        else
        {
            printMsg(": NN> b\n");
            return NN_F;
        }
    }

    inline size_t similarityStatus(size_t child, seq_type signature)
    {
        printMsg("<%zu, ", child);
        size_t status;
        //?
        double offset = 0;

        if (isRootNode[child])
        {
            //? root can be stay or split only
            double similarity = calcOverlap(signature, means[child]);
            printMsg("%.2f> r\n ", similarity);
            if (similarity >= stay_threshold)
            {
                return STAY_F;
            }
            else
            {
                return SPLIT_F;
            }
        }
        else if (isSuperNode[child])
        {
            // super can only be stay or split
            status = similarityStatus(means[child], signature, priority[child], 100);
        }
        else
        {
            status = similarityStatus(means[child], signature, stay_threshold);
        }

        if (status == NN_F)
        {
            if (isBranchNode[child])
            {
                return NN_BRANCH_F;
            }
            else
            {
                return NN_LEAVE_F;
            }
        }
        return status;
    }

    template <typename bset>
    inline size_t getMethod(bset statuses, bool multiple = false)
    {
        if (statuses.none())
        {
            return 0;
        }

        // bset temp(std::string("00100010"));

        if ((statuses & bset("00001010")).any())
        {
            if (statuses.test(5))
            {
                return 2;
            }
            else
            {
                return 1;
            }
        }
        else if (statuses == bset().set(5))
        {
            if (multiple)
            {
                return 4;
            }
            else
            {
                return 3;
            }
        }
        else if (statuses == bset().set(2))
        {
            if (multiple)
            {
                return 7;
            }
            else
            {
                return 5;
            }
        }
        else if (statuses == bset().set(4))
        {
            if (multiple)
            {
                return 7;
            }
            else
            {
                return 6;
            }
        }
        else if (statuses == bset("00010100"))
        {
            return 7;
        }
        else if (statuses == bset("00100100"))
        {
            return 8;
        }
        else if (statuses == bset("00110000"))
        {
            return 9;
        }
        else if (statuses == bset("00110100"))
        {
            return 10;
        }
        else
        {
            printMsg("ERROR - status\n");
            return 11;
        }
        return 0;
    }

    inline pair<size_t, size_t> checkRootP(seq_type signature, vector<size_t> &insertionList, vector<size_t> &mismatch, vector<size_t> &NN_leaves, vector<size_t> &NN_branches, vector<size_t> &stay, size_t node = 0)
    {
        std::bitset<8> statuses;
        size_t dest = 0;
        double max_similarity = 0;

        // check the other children
        for (size_t child : childLinks[node])
        {
            // skip ambiNode
            if (isAmbiNode[child])
            {
                printMsg("Ambi: %zu\n", child);
                continue;
            }
            double similarity = 0;
            size_t status = similarityStatus(child, signature);
            switch (status)
            {
            case STAY_F:
                stay.push_back(child);

                if (isSuperNode[child] || isRootNode[child])
                {
                    statuses.set(5);
                }
                else
                {
                    similarity = calcSimilarity(means[child], signature);
                    if (similarity > max_similarity)
                    {
                        max_similarity = similarity;
                        dest = child;
                    }

                    if (isBranchNode[child])
                    {
                        statuses.set(3);
                    }
                    else
                    {
                        statuses.set(1);
                    }
                }

                break;
            case SPLIT_F:
                mismatch.push_back(child);
                break;
            case NN_LEAVE_F:
                NN_leaves.push_back(child);
                statuses.set(2);
                break;
            case NN_BRANCH_F:
                NN_branches.push_back(child);
                statuses.set(4);
                break;
            default:
                break;
            }
        }

        // cerr << statuses << "\n";
        // cerr << "Status: " << getMethod(statuses, (NN_branches.size() > 1 || NN_leaves.size() > 1)) << "\n";

        size_t method = getMethod(statuses, (NN_branches.size() > 1 || NN_leaves.size() > 1));

        return make_pair(method, dest);
    }

    inline size_t checkRoot(seq_type signature, vector<size_t> &insertionList, vector<size_t> &mismatch, vector<size_t> &NN_leaves, vector<size_t> &NN_branches, vector<size_t> &stay, size_t node = 0)
    {
        std::bitset<8> statuses;
        size_t dest = 0;
        double max_similarity = 0;

        // check the other children
        for (size_t child : childLinks[node])
        {
            // skip ambiNode
            if (isAmbiNode[child])
            {
                printMsg("Ambi: %zu\n", child);
                continue;
            }
            double similarity = 0;
            size_t status = similarityStatus(child, signature);
            switch (status)
            {
            case STAY_F:
                stay.push_back(child);
                similarity = calcSimilarity(means[child], signature);
                if (similarity > max_similarity)
                {
                    max_similarity = similarity;
                    dest = child;
                }

                if (isSuperNode[child] || isRootNode[child])
                {
                    statuses.set(5);
                }
                else if (isBranchNode[child])
                {
                    statuses.set(3);
                }
                else
                {
                    statuses.set(1);
                }
                break;
            case SPLIT_F:
                mismatch.push_back(child);
                break;
            case NN_LEAVE_F:
                NN_leaves.push_back(child);
                statuses.set(2);
                break;
            case NN_BRANCH_F:
                NN_branches.push_back(child);
                statuses.set(4);
                break;
            default:
                break;
            }
        }

        // cerr << statuses << "\n";
        // cerr << "Status: " << getMethod(statuses, (NN_branches.size() > 1 || NN_leaves.size() > 1)) << "\n";

        size_t method = getMethod(statuses, (NN_branches.size() > 1 || NN_leaves.size() > 1));

        return dest;
    }

    // put NN leaves and branches under the same supercluster
    // create an ambiNode to store seqs that are NN to some or all the NN
    // all children in this supercluster except for the ambiNode should be "split"
    inline size_t superCluster(seq_type signature, vector<size_t> &insertionList, size_t idx, size_t node, vector<size_t> NN_leaves, vector<size_t> NN_branches)
    {
        size_t t_parent = createParent(node, insertionList);
        printMsg("*** super\n");
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
        updateParentMean(node);
        if (ambiLinks[node].size() > 0)
        {
            // shouldn't happen at root
            printMsg("check ambi\n");
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

        size_t method, dest;
        std::tie(method, dest) = checkRootP(signature, insertionList, mismatch, NN_leaves, NN_branches, stay, node);

        switch (method)
        {
        case 0:
            printMsg("#mismatch\n");
            return createNode(signature, insertionList, node, idx);
            break;
        case 1:
            printMsg("#stay in %zu from %zu\n", dest, stay.size());
            return stayNode(signature, insertionList, idx, dest);
            break;
        case 2:
            printMsg("#ignore super, stay in %zu from %zu\n", dest, stay.size());
            return stayNode(signature, insertionList, idx, dest);
            break;
        case 3:
            printMsg("#stay in one super %zu from %zu\n", dest, stay.size());
            return tt_root(signature, insertionList, idx, stay[0]);
            break;
        case 4:
            printMsg("#stay multiple super %zu\n", stay.size());
            break;
        case:
            printMsg("# %zu\n", );
            break;
        case:
            printMsg("# %zu\n", );
            break;
        case:
            printMsg("# %zu\n", );
            break;
        case:
            printMsg("# %zu\n", );
            break;
        case:
            printMsg("# %zu\n", );
            break;
        default:
            printMsg("ERROR!!!\n");
            break;
        }

        return 0;

        if (dest != 0)
        {
            printMsg("#stay in %zu from %zu\n", dest, stay.size());
            return stayNode(signature, insertionList, idx, dest);
        }

        if (mismatch.size() == childCounts[node] - ambiLinks[node].size())
        {
            printMsg("#mismatch\n");
            return createNode(signature, insertionList, node, idx);
        }

        size_t NN_total = NN_leaves.size() + NN_branches.size();

        if (NN_branches.size() == 0)
        {
            {
                printMsg("NN leaves\n");
                return superCluster(signature, insertionList, idx, node, NN_leaves, NN_branches);
            }
        }
        else if (NN_branches.size() == 1)
        {
            printMsg("#match ONE branch > ");
            dest = NN_branches[0];
            if (NN_total == 1)
            {
                // no NN leaves
                printMsg("?no NN leaves\n");
                return tt_branch(signature, insertionList, idx, dest);
            }
            else if (isRootNode[dest])
            {
                printMsg("?with NN leaves & is subtree\n");
                // ignore leaves for now
                return tt_root(signature, insertionList, idx, dest);
            }
            else
            {
                printMsg("?with NN leaves\n");
                return superCluster(signature, insertionList, idx, node, NN_leaves, NN_branches);
            }
        }
        else
        {
            printMsg("#match multiple branches > ");
            double max_sim = 0;
            for (size_t c : NN_branches)
            {
                double similarity = calcSimilarity(means[c], signature);
                if (similarity > max_sim)
                {
                    max_sim = similarity;
                    dest = c;
                }
            }

            if (NN_leaves.size() == 0)
            {
                // continue with nearest, might change
                printMsg("?without leaves\n");
                return tt_branch(signature, insertionList, idx, dest);
            }
            else
            {
                // if there is any leaf has higher similarity than all the branches,
                // ignore NN branches, supercluster on the NN leaves only
                // or else, continue with the nearest branch and ignore leaves
                // might change
                printMsg("?with leaves\n");
                for (size_t c : NN_leaves)
                {
                    double similarity = calcSimilarity(means[c], signature);
                    if (similarity > max_sim)
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

        return createNode(signature, insertionList, node, idx);
    }

    // promote subtree only if the subtree and its parent are not root(0)
    // AND
    // if distortion[node]<distortion[parent]
    // OR
    // if the similarity btw the subtree and its parent is smaller than the split threshold
    inline size_t promoteSubtree(size_t node)
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

        if (childCounts[parent] < tree_order)
        {
            return 0;
        }

        if (isRootNode[node])
        {
            if (priority[node] >= priority[parent])
            {
                return 0;
            }
        }
        else
        {
            double similarity = calcSimilarity(means[node], means[parent]);
            if (similarity > split_threshold)
            {
                return 0;
            }
        }

        size_t grandparent = parentLinks[parent];
        printMsg("promoting %zu from %zu to %zu\n", node, parent, grandparent);
        printTreeJson(stderr);
        moveParent(node, grandparent);
        updateNodeMean(parent);
        updatePriority(parent);

        updateNodeMean(grandparent);
        updatePriority(grandparent);
        printTreeJson(stderr);
        return promoteSubtree(node);
    }

    // check the children in stay_subtree to see if any of them is closer to sig
    // sig is stored in t_parent
    // if yes, move the child to t_parent
    // return child with higest similarity with the input sig
    inline size_t moveFromSubtree(vector<size_t> stay_subtree, size_t t_parent, seq_type signature)
    {
        size_t best_child = 0;
        double max_similarity = stay_threshold;
        // create new subtree and move overlapping children out
        //? maybe make it as branch rather than subtree if distance is close
        // t_parent = createParent(node, insertionList);
        for (size_t subtree : stay_subtree)
        {
            for (size_t child : childLinks[subtree])
            {
                double sim1 = calcSimilarity(signature, means[child]);
                double sim2 = calcSimilarity(means[subtree], means[child]);

                if (sim1 > sim2)
                {
                    printMsg("moving %zu\n", child);
                    moveParent(child, t_parent);
                }

                if (sim1 > max_similarity)
                {
                    max_similarity = sim1;
                    best_child = child;
                }

                // if (sim1 > priority[subtree])
                // {
                //     moveParent(child, t_parent);
                // }
            }
            updateNodeMean(subtree);
            updatePriority(subtree);
        }
        return best_child;
    }

    inline size_t insertStaySubtree(vector<size_t> stay_subtree, seq_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        size_t dest = 0;
        if (stay_subtree.size() == 0)
        {
            printMsg("@@ Stay 0\n");
            dest = createNode(signature, insertionList, node, idx);
            //? check if need to move children from other subtree to minimise its distortion
            return dest;
        }
        else if (stay_subtree.size() == 1)
        {
            printMsg("## Stay 1\n");
            return tt_rootSplit(signature, insertionList, idx, stay_subtree[0]);
        }
        else
        {
            printMsg("$$ Stay many %zu\n", stay_subtree.size());
            size_t t_parent = createParent(node, insertionList);
            addSigToMatrix(node, means[t_parent]);
            size_t best_child = moveFromSubtree(stay_subtree, t_parent, signature);
            if (best_child == 0)
            {
                isRootNode[t_parent] = 1;
                dest = createNode(signature, insertionList, t_parent, idx);
            }
            else
            {
                printMsg("^^ proceed with  %zu\n", best_child);
                dest = tt_rootSplit(signature, insertionList, idx, best_child);
            }

            updateParentMean(t_parent);
            return dest;
        }
    }

    //? for now just move into the first level of the subtree
    //? later need to change to find the best level
    inline void moveIntoSubtree(size_t node, vector<size_t> stay_subtree)
    {
        // moveParent(node, subtree);
        double best_similarity = 0;
        size_t dest = 0;
        for (size_t subtree : stay_subtree)
        {
            double similarity = calcSimilarity(means[node], means[subtree]);
            if (similarity > best_similarity)
            {
                best_similarity = similarity;
                dest = subtree;
            }
        }
        moveParent(node, dest);
    }

    // use this when children are all subtrees
    // find the nearest subtree
    // split subtree if distortion is too big
    inline size_t tt_rootSplit(seq_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        // if (node != root && isRootNode[node])
        // {
        //     if (childCounts[node] > findLevel(node) * tree_order)
        //     {
        //         forceSplitRoot(insertionList, false, node);
        //     }

        //     // promote subtree if distance to parent greater than split threshold
        //     promoteSubtree(node);
        // }

        if (childCounts[node] > 0)
        {
            if (!isRootNode[childLinks[node][0]])
            {
                return tt_root(signature, insertionList, idx, node);
            }
        }

        vector<size_t> mismatch_subtree;
        vector<size_t> stay_subtree;
        vector<size_t> nonSubTree;

        for (size_t subtree : childLinks[node])
        {
            if (isAmbiNode[subtree])
            {
                continue;
            }

            if (isRootNode[subtree])
            {
                double similarity = calcSimilarity(signature, means[subtree]);
                if (similarity >= priority[subtree])
                {
                    stay_subtree.push_back(subtree);
                }
                else
                {
                    mismatch_subtree.push_back(subtree);
                }
            }
            else
            {
                nonSubTree.push_back(subtree);
            }
        }

        size_t t_parent = 0;
        size_t dest = 0;
        if (nonSubTree.size() > 0)
        {
            printMsg("!! Contains %zu non subtree, stay in %zu subtree\n", nonSubTree.size(), stay_subtree.size());
            //? do this for now, later change to check against all non subtree
            // return tt_root(signature, insertionList, idx, nonSubTree[0]);
            vector<size_t> mismatch;
            vector<size_t> NN_leaves;
            vector<size_t> NN_branches;
            vector<size_t> stay;

            for (size_t child : nonSubTree)
            {
                size_t status = similarityStatus(child, signature);
                switch (status)
                {
                case STAY_F:
                    stay.push_back(child);
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

            size_t NN_total = NN_branches.size() + NN_leaves.size();
            if (mismatch.size() == nonSubTree.size())
            {
                printMsg("--- Non subtrees all mismatch\n");
                return insertStaySubtree(stay_subtree, signature, insertionList, idx, node);
            }
            else if (nonSubTree.size() == 1)
            {
                // if the only nonsubtree is mismatch, would have caught it in the prev if block
                size_t temp_node = nonSubTree[0];
                if (stay.size() == 1)
                {
                    printMsg("!!! stay 1 leave\n");
                    if (stay_subtree.size() == 0)
                    {
                        return stayNode(signature, insertionList, idx, stay[0]);
                    }
                    else
                    {
                        printMsg("      moving\n");
                        t_parent = createParent(node, insertionList);
                        addSigToMatrix(node, means[temp_node]); // this should be updated later
                        moveParent(temp_node, t_parent);

                        size_t best_child = moveFromSubtree(stay_subtree, t_parent, signature);
                        //?
                        if (best_child == 0)
                        {
                            return tt_branch(signature, insertionList, idx, temp_node);
                        }
                        else
                        {
                            return tt_rootSplit(signature, insertionList, idx, best_child);
                        }
                    }
                }
                else if (NN_branches.size() == 1)
                {
                    printMsg("@@@ NN 1 branch\n");
                    if (stay_subtree.size() == 0)
                    {
                        return tt_branch(signature, insertionList, idx, temp_node);
                    }
                    else
                    {
                        printMsg("      moving\n");
                        t_parent = temp_node;
                        size_t best_child = moveFromSubtree(stay_subtree, t_parent, signature);
                        //?
                        if (best_child == 0)
                        {
                            return tt_branch(signature, insertionList, idx, temp_node);
                        }
                        else
                        {
                            return tt_rootSplit(signature, insertionList, idx, best_child);
                        }
                    }
                }
                else
                {
                    //?
                    printMsg("$$$ NN 1 leaf\n");
                    return superCluster(signature, insertionList, idx, t_parent, NN_leaves, NN_branches);
                }
            }
            // below all means nonsubtree > 1

            else if (NN_total == 1 || NN_leaves.size() == NN_total)
            {
                if (NN_total > 1 && NN_leaves.size() == NN_total)
                {
                    printMsg("Groun NN leaves first\n");
                    t_parent = createParent(node, insertionList);
                    addSigToMatrix(node, signature); // this should be updated later
                    for (size_t child : NN_leaves)
                    {
                        moveParent(child, t_parent);
                    }
                    updateParentMean(t_parent);
                    NN_leaves.clear();
                    NN_branches.push_back(t_parent);
                }

                if (stay.size() == 0)
                {
                    if (stay_subtree.size() == 0)
                    {
                        printMsg("!!!! only 1 NN (leaf or branch), no stay in subtree or nonsubtree\n");
                        if (NN_branches.size() == 1)
                        {
                            return tt_branch(signature, insertionList, idx, NN_branches[0]);
                        }
                        else
                        {
                            return createAmbiNode(signature, insertionList, idx, NN_leaves[0]);
                        }
                    }
                    else
                    {
                        printMsg("!!!! only 1 NN (leaf or branch), no stay nonsubtree,but stay in subtrees\n");

                        double max_similarty = 0;
                        size_t best_subtree = 0;
                        if (NN_branches.size() == 1)
                        {
                            moveIntoSubtree(NN_branches[0], stay_subtree);
                        }
                        else
                        {
                            moveIntoSubtree(NN_leaves[0], stay_subtree);
                        }

                        for (size_t subtree : stay_subtree)
                        {
                            updateNodeMean(subtree);
                            updatePriority(subtree);
                        }
                        updateParentMean(node);
                        return insertStaySubtree(stay_subtree, signature, insertionList, idx, node);
                        // return tt_rootSplit(signature, insertionList, idx, t_parent);
                    }
                }
                else
                {
                    if (stay_subtree.size() == 0)
                    {
                        if (t_parent != 0)
                        {
                        }
                        else if (NN_branches.size() > 0)
                        {
                            t_parent = NN_branches[0];
                            // NN 1 branch
                        }
                        else
                        { // NN 1 leave
                            t_parent = createParent(node, insertionList);
                            moveParent(NN_leaves[0], t_parent);
                        }

                        for (size_t child : stay)
                        {
                            moveParent(child, t_parent);
                        }
                        updateParentMean(t_parent);
                        return tt_branch(signature, insertionList, idx, t_parent);
                    }
                }
            }

            else if (NN_total > 0) // there might still be mismatch and stay
            {
                if (stay.size() == 0)
                {
                    if (NN_total == 1)
                    {
                        if (stay_subtree.size() == 0)
                        {
                            printMsg("!!!! only 1 NN (leaf or branch), no stay in subtree or nonsubtree\n");
                            if (NN_branches.size() == 1)
                            {
                                return tt_branch(signature, insertionList, idx, NN_branches[0]);
                            }
                            else
                            {
                                return createAmbiNode(signature, insertionList, idx, NN_leaves[0]);
                            }
                        }
                        else
                        {
                            printMsg("!!!! only 1 NN (leaf or branch), no stay nonsubtree,but stay in subtrees\n");

                            double max_similarty = 0;
                            size_t best_subtree = 0;
                            if (NN_branches.size() == 1)
                            {
                                moveIntoSubtree(NN_branches[0], stay_subtree);
                            }
                            else
                            {
                                moveIntoSubtree(NN_leaves[0], stay_subtree);
                            }

                            for (size_t subtree : stay_subtree)
                            {
                                updateNodeMean(subtree);
                                updatePriority(subtree);
                            }
                            updateParentMean(node);
                            return insertStaySubtree(stay_subtree, signature, insertionList, idx, node);
                            // return tt_rootSplit(signature, insertionList, idx, t_parent);
                        }
                    }
                    else if (NN_leaves.size() == NN_total)
                    {
                    }
                }
                if (stay.size() > 0)
                {
                    if (NN_total == 1)
                    {
                        if (NN_branches.size() == 1)
                        {
                            t_parent == NN_branches[0];
                        }
                        else
                        {
                            size_t temp_node = NN_leaves[0];
                            t_parent = createParent(node, insertionList);
                            addSigToMatrix(node, means[temp_node]); // this should be updated later
                            moveParent(temp_node, t_parent);
                        }
                    }
                    else
                    {
                        size_t temp_node = nonSubTree[0];
                        t_parent = createParent(node, insertionList);
                        addSigToMatrix(node, means[temp_node]); // this should be updated later
                        moveParent(temp_node, t_parent);
                    }

                    printMsg("!!! stay 1 leave\n");
                    if (stay_subtree.size() == 0)
                    {
                        return stayNode(signature, insertionList, idx, stay[0]);
                    }
                    else
                    {
                        printMsg("      moving\n");
                        t_parent = createParent(node, insertionList);
                        addSigToMatrix(node, signature); // this should be updated later

                        size_t best_child = moveFromSubtree(stay_subtree, t_parent, signature);
                        printMsg("      best child %zu\n", best_child);
                        //?
                        if (best_child == 0)
                        {
                            return tt_branch(signature, insertionList, idx, t_parent);
                        }
                        else
                        {
                            return tt_rootSplit(signature, insertionList, idx, best_child);
                        }
                    }
                }
                else // no stay, but might have mismatch. Here we just ignore mismatch
                {
                }
            }
            // by now means some or all of the nonsubtrees are stay, no NN here
            else
            {
                //? is it possible to stay in multiple??

                t_parent = createParent(node, insertionList);
                addSigToMatrix(node, signature); // this should be updated later

                if (stay_subtree.size() == 0)
                {
                    for (size_t child : stay)
                    {
                        moveParent(child, t_parent);
                    }
                    return tt_branch(signature, insertionList, idx, t_parent);
                }

                isRootNode[t_parent] = 1;
                // move nonsubtree into the stay_subtree if the similarity > priority
                double dest_similarity = 0;
                for (size_t child : stay)
                {
                    double child_similarity = calcSimilarity(means[child], signature);
                    if (child_similarity > dest_similarity)
                    {
                        dest_similarity = child_similarity;
                        dest = child;
                    }

                    double max_similarity = 0;
                    size_t best_loc = 0;
                    for (size_t subtree : stay_subtree)
                    {
                        double similarity = calcSimilarity(means[child], means[subtree]);
                        if (similarity > max_similarity)
                        {
                            max_similarity = similarity;
                            best_loc = subtree;
                        }
                    }

                    if (max_similarity > priority[best_loc])
                    {
                        printMsg("      moving %zu to %zu (%.2f)\n", child, best_loc, max_similarity);
                        moveParent(child, best_loc);
                    }
                    else
                    {
                        printMsg("      not worth moving %zu to %zu (%.2f)\n", child, best_loc, max_similarity);
                    }

                    //? is it fair to compare subtree with nonsubtree?
                    if (max_similarity > dest_similarity)
                    {
                        dest_similarity = max_similarity;
                        dest = best_loc;
                    }
                }
                if (isRootNode[dest])
                {
                    return tt_rootSplit(signature, insertionList, idx, dest);
                }
                else
                {
                    return stayNode(signature, insertionList, idx, dest);
                }
            }
            printMsg("ERROR insert went wrong\n");
            return 0;
        }
        else
        {
            printMsg("--- Subtrees only\n");
            return insertStaySubtree(stay_subtree, signature, insertionList, idx, node);
        }

        // if (maxOverlap < maxSimilarity || maxOverlap < 0.5)
        // {
        //     printMsg("#Retry subtree %zu (%.2f,%.2f)\n", node, maxOverlap, maxSimilarity);
        //     return tt_root(signature, insertionList, idx, node);
        // }

        // printMsg("#stay in subtree %zu (%.2f)\n", bestSubtree, maxOverlap);
        // if (childCounts[bestSubtree] > 4 && priority[bestSubtree] < split_threshold)
        // {
        //     addSubtree(insertionList, false, bestSubtree);
        // }

        // return tt_rootSplit(signature, insertionList, idx, bestSubtree);
    }

    inline size_t forceSplitSubtree(vector<size_t> &insertionList, size_t node)
    {
        updatePriority(node);
        return 0;
        if (childCounts[node] <= tree_order && priority[node] > split_threshold)
        {
            return 0;
        }

        printMsg(">> Force split subtree %zu\n", node);
        forceSplitRoot(insertionList, false, node);
        size_t parent = parentLinks[node];
        for (size_t child : childLinks[node])
        {
            moveParent(child, parent, false);
        }
        deleteNode(node);

        return 1;
    }

    inline size_t tt_branch(seq_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        if (isRootNode[node])
        {
            printMsg("is subtree %zu\n", node);
            forceSplitSubtree(insertionList, node);
            return tt_root(signature, insertionList, idx, node);
        }

        size_t current_childCount = childCounts[node];
        vector<size_t> mismatch;
        vector<size_t> NN_leaves;
        vector<size_t> NN_branches;
        vector<size_t> stay;

        size_t dest = checkRoot(signature, insertionList, mismatch, NN_leaves, NN_branches, stay, node);
        printMsg("tt_branch %zu\n", node);

        if (dest != 0)
        {
            printMsg("#b- stay in %zu from %zu\n", dest, stay.size());
            return stayNode(signature, insertionList, idx, dest);
        }

        if (mismatch.size() == childCounts[node] - ambiLinks[node].size())
        {
            //? this might be a wrong supercluster, dissolve this branch and check from parent again
            printMsg("//?#b- mismatch all\n");
            dest = createNode(signature, insertionList, node, idx);
            updateParentMean(node);
            return dest;
        }
        else if (mismatch.size() > 0)
        {
            printMsg("//?#b- mismatch, promote\n");
        }

        size_t NN_total = NN_leaves.size() + NN_branches.size();
        if (mismatch.size() == 0)
        {
            // NN with everything
            printMsg("#??NN with all\n");
            if (ambiLinks[node].size() == 0)
            {
                return createAmbiNode(signature, insertionList, node, idx);
            }
            else
            {
                return stayAmbi(signature, insertionList, idx, node);
            }
        }
        else
        {
            // some mismatches
            printMsg("#NN with some\n");
            if (ambiLinks[node].size() == 0)
            {
                return superCluster(signature, insertionList, idx, node, NN_leaves, NN_branches);
            }
            else
            {
                return stayAmbi(signature, insertionList, idx, node);
            }
        }
        printMsg("//?#b- wrong\n");
        return createNode(signature, insertionList, node, idx);
    }

    inline size_t tt(seq_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        if (isRootNode[node])
        {
            return tt_root(signature, insertionList, idx);
        }
        else if (isBranchNode[node])
        {
            return tt_branch(signature, insertionList, idx, node);
        }
        else
        {
            printMsg("Something went wrong\n");
        }
    }

    inline size_t first_insert(seq_type signature, vector<size_t> &insertionList, size_t idx, size_t parent = 0)
    {
        return createNode(signature, insertionList, root, idx);
    }

    inline size_t insert(seq_type signature, vector<size_t> &insertionList, size_t idx)
    {
        size_t node = tt(signature, insertionList, idx);
        printMsg("inserted %zu at %zu\n\n", idx, node);
        return node;
    }

    inline size_t insertSplit(seq_type signature, vector<size_t> &insertionList, size_t idx)
    {
        size_t node = tt_rootSplit(signature, insertionList, idx);
        printMsg("inserted %zu at %zu\n\n", idx, node);
        return node;
    }

    inline size_t insertSplitRoot(seq_type signature, vector<size_t> &insertionList, size_t idx)
    {
        size_t node = insertSplit(signature, insertionList, idx);
        if (childCounts[root] > tree_order)
        {
            forceSplitRoot(insertionList, false);

            // //?
            // vector<size_t> subtrees = childLinks[root];
            // for (size_t t_parent : subtrees)
            // {
            //     forceSplitSubtree(insertionList, t_parent);
            // }
        }
        return node;
    }

    // inline size_t search(seq_type signature, size_t idx = 0, size_t node = 0)
    // {
    //     // size_t node = root;

    //     size_t local_best_child = node;
    //     double local_best_similarity = 0;
    //     double local_best_similarity_b = 0;
    //     size_t local_best_child_b = node;

    //     for (size_t i = 0; i < childCounts[node]; i++)
    //     {
    //         size_t child = childLinks[node][i];
    //         if (isAmbiNode[child])
    //         {
    //             continue;
    //         }

    //         // double similarity = calcSimilarity(means[child], signature);
    //         double similarity = calcOverlap(signature, means[child]);
    //         printMsg(" <%zu,%.2f> ", child, similarity);

    //         // found same, move on with the next seq
    //         if (similarity >= (stay_threshold))
    //         {
    //             // while (isBranchNode[child])
    //             // {
    //             //     child = childLinks[child][0];
    //             // }
    //             if (isBranchNode[child])
    //             {
    //                 if (local_best_similarity_b == 1)
    //                 {
    //                     return search(signature, idx, child);
    //                 }
    //                 else if (similarity >= local_best_similarity_b)
    //                 {
    //                     local_best_similarity_b = similarity;
    //                     local_best_child_b = child;
    //                 }
    //             }
    //             else
    //             {
    //                 return child;
    //             }
    //         }
    //         else
    //         {

    //             if (isBranchNode[child])
    //             {
    //                 if (similarity >= local_best_similarity_b)
    //                 {
    //                     local_best_similarity_b = similarity;
    //                     local_best_child_b = child;
    //                 }
    //             }
    //             else
    //             {
    //                 if (similarity >= local_best_similarity)
    //                 {
    //                     local_best_similarity = similarity;
    //                     local_best_child = child;
    //                 }
    //             }
    //         }
    //     }

    //     // // stay in branch, proceed with children
    //     // if (local_best_child_b != node)
    //     // {
    //     //     return search(signature, idx, local_best_child_b);
    //     // }

    //     if (local_best_similarity >= local_best_similarity_b)
    //     {
    //         return local_best_child;
    //     }
    //     else
    //     {
    //         return search(signature, idx, local_best_child_b);
    //     }
    // }

    // find highest overlap among children
    // break tie by checking similarity of the finals
    inline size_t search(seq_type signature, size_t idx = 0, size_t node = 0)
    {
        double best_similarity = 0;
        vector<size_t> candidates;

        for (size_t i = 0; i < childCounts[node]; i++)
        {
            size_t child = childLinks[node][i];
            if (isAmbiNode[child])
            {
                continue;
            }

            double similarity = calcSimilarity(means[child], signature);
            similarity += calcOverlap(signature, means[child]);
            printMsg(" <%zu,%.2f> ", child, similarity);
            printMsg(" (%.2f, %.2f, %.2f) ", calcSimilarity(means[child], signature), calcOverlap(signature, means[child]), priority[child]);

            if (similarity > best_similarity)
            {
                best_similarity = similarity;
                candidates.clear();
                // best_child = child;
            }

            if (similarity == best_similarity)
            {
                candidates.push_back(child);
            }
        }

        size_t best_child = candidates[0];

        if (candidates.size() == 1)
        {
            if (isBranchNode[best_child])
            {
                printMsg("\n>%zu\n", best_child);
                return search(signature, idx, best_child);
            }
            else
            {
                return best_child;
            }
        }
        else
        {

            printMsg("\n*** here\n");
            double best_dest_similarity = 0;
            for (size_t child : candidates)
            {
                printMsg("child %zu\n", child);
                size_t leaf = child;
                double similarity = best_similarity;
                if (isBranchNode[child])
                {
                    leaf = search(signature, idx, child);
                    similarity = calcSimilarity(means[child], signature);
                    similarity += calcOverlap(signature, means[leaf]);
                    printMsg("> leaf %zu\n", leaf);
                }

                if (similarity > best_dest_similarity)
                {
                    best_dest_similarity = similarity;
                    best_child = leaf;
                }
            }

            return best_child;
        }
    }

    // find highest overlap among children
    // break tie by checking similarity of the finals
    inline size_t search2(seq_type signature, size_t idx = 0, size_t node = 0)
    {
        printMsg("\n");
        double best_similarity = 0;
        size_t best_branch = 0;

        double best_overlap = 0;
        size_t best_subtree = 0;
        vector<size_t> overlapAll;

        for (size_t i = 0; i < childCounts[node]; i++)
        {
            size_t child = childLinks[node][i];
            if (isAmbiNode[child])
            {
                continue;
            }
            else if (isRootNode[child])
            {
                double overlap = calcOverlap(signature, means[child]);
                printMsg(" <%zu,o-%.2f> ", child, overlap);
                if (overlap == 1)
                {
                    overlapAll.push_back(child);
                }
                else if (overlap == best_overlap)
                {
                    double s1 = calcSimilarity(means[best_subtree], signature);
                    double s2 = calcSimilarity(means[child], signature);
                    if (s2 > s1)
                    {
                        printMsg(" <%zu,s2-%.2f,%.2f> ", child, s1, s2);
                        best_subtree = child;
                    }
                }
                else if (overlap > best_overlap)
                {
                    best_overlap = overlap;
                    best_subtree = child;
                }
            }
            else
            {
                double similarity = calcSimilarity(means[child], signature);
                printMsg(" <%zu,s-%.2f> ", child, similarity);
                if (similarity > best_similarity)
                {
                    best_similarity = similarity;
                    best_branch = child;
                }
            }
        }

        if (overlapAll.size() > 0)
        {
            best_subtree = overlapAll[0];
        }
        printMsg("@@@ %zu, %zu", best_subtree, overlapAll.size());
        // if node has no subtree
        if (best_subtree == 0)
        {
            if (isBranchNode[best_branch])
            {
                return search2(signature, idx, best_branch);
            }
            else
            {
                return best_branch;
            }
        }
        else if (overlapAll.size() <= 1)
        {
            return search2(signature, idx, best_subtree);
        }
        else
        {
            printMsg("\n---");
            double best_subtree_similarity = 0;
            for (size_t child : overlapAll)
            {
                double similarity = calcSimilarity(means[child], signature);
                printMsg(" <%zu,s-%.2f> ", child, similarity);
                if (similarity > best_subtree_similarity)
                {
                    best_subtree_similarity = similarity;
                    best_subtree = child;
                }
            }

            if (best_branch == 0)
            {
                return search2(signature, idx, best_subtree);
            }
            else
            {
                if (best_subtree_similarity > best_similarity)
                {
                    return search2(signature, idx, best_subtree);
                }
                if (isBranchNode[best_branch])
                {
                    return search2(signature, idx, best_branch);
                }
                else
                {
                    return best_branch;
                }
            }
        }
    }

    // find highest overlap among children
    // break tie by checking similarity of the finals
    inline size_t search3(seq_type signature, size_t idx = 0, size_t node = 0, double best_score = 0, double dest = 0)
    {
        printMsg("\n");
        size_t temp_dest = 0;
        double temp_best_similarity = best_score;

        for (size_t child : childLinks[node])
        {
            if (isAmbiNode[child])
            {
                continue;
            }
            else if (isRootNode[child] || isBranchNode[child])
            {
                size_t best_leaf = search3(signature, idx, child, temp_best_similarity);
                if (best_leaf != 0)
                {
                    double similarity = calcSimilarity(means[best_leaf], signature);
                    printMsg("<%zu,o-%.2f>\n", child, similarity);
                    if (best_leaf != 0 && similarity > temp_best_similarity)
                    {
                        temp_best_similarity = similarity;
                        temp_dest = child;
                    }
                }
            }
            else
            {
                double similarity = calcSimilarity(means[child], signature);
                printMsg(" <%zu,s-%.2f> ", child, similarity);
                if (similarity > temp_best_similarity)
                {
                    temp_best_similarity = similarity;
                    temp_dest = child;
                }
            }
        }

        return temp_dest;
    }

    // cluster means still remain
    void prepReinsert(size_t node = 0)
    {
        if (!isAmbiNode[node])
        {
            updatePriority(node);
            if (!isBranchNode[node])
            {
                updateNodeMean(node);
            }
        }
        seqIDs[node].clear();
        for (size_t child : childLinks[node])
        {
            prepReinsert(child);
        }
    }

    void updateTree(size_t node = 0)
    {
        if (!isAmbiNode[node])
        {
            updatePriority(node);
            if (!isBranchNode[node])
            {
                updateNodeMean(node);
            }
        }
        for (size_t child : childLinks[node])
        {
            updateTree(child);
        }
    }

    size_t reinsert(seq_type signature, size_t idx)
    {
        // size_t node = search(signature);
        size_t node = search2(signature);
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
            else if (isAmbiNode[child] | seqIDs[child].size() < singleton) // remove singleton
            {
                leaves.push_back(child);
            }
            else
            {
                updatePriority(child);
                updateParentMean(node);
            }
        }
    }

    bool redundant(size_t parent)
    {
        printMsg("redundant %zu\n", parent);
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
        printMsg("\nTrimming\n");

        vector<size_t> leaves;
        getEmptyLeaves(leaves);

        // do leaves
        for (size_t i : leaves)
        {
            size_t parent = parentLinks[i];
            deleteNode(i);
            if (childCounts[parent] == 0)
            {
                branches.insert(parent);
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
            updatePriority(node);
            updateParentMean(node);
        }

        // remove unitig for non-empty leaves
        for (size_t child : childLinks[root])
        {
            dropBranch(child);
        }
    }

    void outputHierarchy(FILE *pFile, size_t node = 0, size_t rank = 0)
    {
        for (size_t child : childLinks[node])
        {
            fprintf(pFile, "%zu,%zu,%zu\n", node, child, rank);
            if (isBranchNode[child])
            {
                outputHierarchy(pFile, child, rank + 1);
            }
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