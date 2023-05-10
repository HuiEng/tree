
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
size_t singleton = 1;
size_t minClusSize = partree_capacity;
size_t tree_order = 5;
bool print_ = false;

template <typename T>
void removeVecValue(vector<T> &vec, T value)
{
    vec.erase(std::remove(vec.begin(), vec.end(), value), vec.end());
}

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
        fprintf(stream, "\"branch\":\"%d\",", isBranchNode[tnode]);
        fprintf(stream, "\"ambi\":\"%d\",", isAmbiNode[tnode]);
        fprintf(stream, "\"subtree\":\"%d\",", isRootNode[tnode]);
        fprintf(stream, "\"super\":\"%d\",", isSuperNode[tnode]);
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

    seq_type createRandomSig(vector<size_t> children, size_t s = 0)
    {
        // unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        unsigned seed = s;
        if (seed == 0)
        {
            seed = chrono::system_clock::now().time_since_epoch().count();
        }
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
        if (isBranchNode[node] || isSuperNode[node] || isRootNode[node])
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
        updatePriority(node);
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
    double test(size_t node = 0)
    {
        for (size_t i = 0; i < childCounts[node]; i++)
        {
            for (size_t j = i + 1; j < childCounts[node]; j++)
            {
                double similarity = calcSimilarity(means[childLinks[node][i]], means[childLinks[node][j]]);
                fprintf(stderr, "%zu,%zu,%.2f\n", childLinks[node][i], childLinks[node][j], similarity);
            }
        }
        return 0;

        // printMsg(">>>\n");
        // node = 10;

        // size_t nodeb = 12;
        // for (size_t j = 0; j < childCounts[node]; j++)
        // {
        //     double similarity = calcSimilarity(means[nodeb], means[childLinks[node][j]]);
        //     printMsg("10,%zu,%.2f\n", childLinks[node][j], similarity);
        // }

        printMsg(">>>\n");
        double similarity = calcSimilarity(means[31], means[32]);
        printMsg("%zu,%zu,%.2f\n", 31, 32, similarity);

        return 0;

        vector<size_t> t{11, 8, 4, 5, 7, 38};
        for (size_t i = 0; i < t.size(); i++)
        {
            for (size_t j = i + 1; j < t.size(); j++)
            {
                double similarity = calcSimilarity(means[t[i]], means[t[j]]);
                printMsg("%zu,%zu,%.2f\n", t[i], t[j], similarity);
            }
        }
        return 0;

        node = 73;
        size_t nodeb = 61;
        for (size_t i = 0; i < childCounts[node]; i++)
        {
            for (size_t j = 0; j < childCounts[nodeb]; j++)
            {
                double similarity = calcSimilarity(means[childLinks[node][i]], means[childLinks[nodeb][j]]);
                printMsg("%zu,%zu,%.2f\n", childLinks[node][i], childLinks[nodeb][j], similarity);
            }
        }
        printMsg(">>>\n");

        nodeb = 67;
        for (size_t i = 0; i < childCounts[node]; i++)
        {
            for (size_t j = 0; j < childCounts[nodeb]; j++)
            {
                double similarity = calcSimilarity(means[childLinks[node][i]], means[childLinks[nodeb][j]]);
                printMsg("%zu,%zu,%.2f\n", childLinks[node][i], childLinks[nodeb][j], similarity);
            }
        }
        return 0;
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
        vector<seq_type> temp_matrix = matrices[node];
        if (isBranchNode[node])
        {
            temp_matrix = getNonAmbiMatrix(node);
        }
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
        ;
        // double avg_sim = sumDistance / temp_matrix.size();
        // if (isRootNode[node])
        // {
        //     double max_width = 0;
        //     for (size_t child : childLinks[node])
        //     {
        //         if (priority[child] > max_width)
        //         {
        //             max_width = priority[child];
        //         }
        //     }
        //     avg_sim += max_width / 2;
        // }
        // return avg_sim;
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
            return tt(signature, insertionList, idx, node);
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
        // printMsg("moving %zu to %zu \n", child, new_parent);

        if (d)
        {
            deleteNode(child);
        }
        parentLinks[child] = new_parent;
        // printMsg("deleted %zu \n", child);
    }

    // haven't update childLinks[t_parent] yet
    inline size_t createParent(size_t node, vector<size_t> &insertionList)
    {
        size_t t_parent = getNewNodeIdx(insertionList);
        isBranchNode[t_parent] = 1;
        parentLinks[t_parent] = node;
        isRootNode[t_parent] = 0;
        isSuperNode[t_parent] = 0;

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
        isSuperNode[new_node] = 0;

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

    // rotate subtree if it has at least 1 subtree child with priority lower than it
    // rotate with the lowest priority
    size_t rotate(size_t t_parent)
    {
        bool need_rotate = false;
        double min_priority = 1;
        size_t candidate = 0;
        size_t candidate_idx = 0;
        for (size_t i = 0; i < childCounts[t_parent]; i++)
        {
            size_t child = childLinks[t_parent][i];
            if (isRootNode[child])
            {
                if (priority[child] < priority[t_parent])
                {
                    need_rotate = true;
                    printMsg("---%zu, %.2f\n", child, priority[child]);
                    if (priority[child] < min_priority)
                    {
                        min_priority = priority[child];
                        candidate = child;
                        candidate_idx = i;
                    }
                }
            }
        }

        if (!need_rotate)
        {
            return 0;
        }

        printMsg("Rotating %zu with %zu\n", t_parent, candidate);
        printTreeJson(stderr);

        size_t grandparent = parentLinks[t_parent];
        size_t t_parent_idx = getNodeIdx(t_parent);

        // remove candidate from t_parent then make candidate the parent of t_parent
        matrices[t_parent].erase(matrices[t_parent].begin() + candidate_idx);
        childLinks[t_parent].erase(childLinks[t_parent].begin() + candidate_idx);
        childCounts[t_parent]--;
        updateNodeMean(t_parent);
        updatePriority(t_parent);

        // add t_parent to candidate, make grandparent the parent of candidate then update grandparent
        moveParent(t_parent, candidate, false);
        matrices[grandparent].erase(matrices[grandparent].begin() + t_parent_idx);
        childLinks[grandparent].erase(childLinks[grandparent].begin() + t_parent_idx);
        moveParent(candidate, grandparent, false);

        if (childCounts[t_parent] <= 1)
        {
            dropBranch(t_parent);
            printMsg("Dropping unitig %zu \n", t_parent);
        }

        updateParentMean(candidate);
        return 1;
    }

    inline pair<size_t, double> similarityStatusF(seq_type sig1, seq_type sig2, double local_stay_t, double offset = 1, bool isBranch = false)
    {
        //?
        if (local_stay_t > stay_threshold)
        {
            local_stay_t = stay_threshold;
        }
        double local_split_t = split_threshold * offset;
        double similarity = calcSimilarity(sig1, sig2);

        printMsg("%.2f (%.2f)", similarity, local_split_t);

        if (similarity >= local_stay_t)
        {
            printMsg(": STAY>\n");
            // return STAY_F;
            return make_pair(STAY_F, similarity);
        }
        else if (similarity <= local_split_t)
        {
            printMsg(": SPLIT>\n");
            // return SPLIT_F;
            return make_pair(SPLIT_F, similarity);
        }
        else
        {
            printMsg(": NN>\n");
            // return NN_F;
            if (isBranch)
            {
                return make_pair(NN_BRANCH_F, similarity);
            }
            else
            {
                return make_pair(NN_LEAVE_F, similarity);
            }
        }
    }

    inline pair<size_t, double> similarityStatus(size_t child, seq_type signature)
    {
        printMsg("<%zu, ", child);

        if (isSuperNode[child])
        {
            // super can only be stay or split
            printMsg("(super %.2f)", priority[child]);
            return similarityStatusF(means[child], signature, priority[child], 100);
        }
        else if (isBranchNode[child])
        {
            // super can only be stay or split
            printMsg("(branch %.2f)", priority[child]);
            return similarityStatusF(means[child], signature, priority[child], 1, true);
        }
        else
        {
            return similarityStatusF(means[child], signature, stay_threshold);
        }
    }

    template <typename bset>
    inline size_t getMethod1(bset statuses)
    {
        size_t pos;
        for (pos = 0; pos < statuses.size(); ++pos)
        {
            if (statuses.test(pos))
            {
                break;
            }
        }
        std::cerr << pos << '\n';
        return 0;
    }

    template <typename bset>
    inline size_t getMethod2(bset statuses)
    {
        return 0;
    }

    template <typename bset>
    inline size_t getMethod3(bset statuses)
    {
        return 0;
    }

    template <typename bset>
    inline size_t getMethod4(bset statuses)
    {
        return 0;
    }

    template <typename bset>
    inline size_t getMethod5(bset statuses)
    {
        return 0;
    }

    template <typename bset>
    inline size_t getMethod(bset statuses)
    {
        if (statuses.none())
        {
            return 0;
        }

        switch (statuses.count())
        {
        case 1:
            return getMethod1(statuses);
            break;
        case 2:
            return getMethod2(statuses);
            break;
        case 3:
            return getMethod3(statuses);
            break;
        case 4:
            return getMethod4(statuses);
            break;
        case 5:
            return getMethod5(statuses);
            break;

        default:
            break;
        }
    }

    // template <typename bset>
    // inline size_t getMethod(bset statuses, size_t staySuper, size_t totalSuper, size_t NNBranch)
    // {
    //     if (statuses.none())
    //     {
    //         return 0;
    //     }

    //     if (statuses == bset().set(1))
    //     {
    //         return 1;
    //     }
    //     else if (statuses == bset().set(2) || statuses == bset().set(1).set(2))
    //     {
    //         return 2;
    //     }
    //     else if (statuses == bset().set(3))
    //     {
    //         if (staySuper == 1)
    //         {
    //             return 3;
    //         }
    //         else if (staySuper == totalSuper)
    //         {
    //             return 4;
    //         }
    //         else
    //         {
    //             return 5;
    //         }
    //     }
    //     else if (statuses == bset().set(2).set(3))
    //     {
    //         return 6;
    //     }
    //     else if (statuses == bset().set(1).set(2).set(3))
    //     {
    //         return 7;
    //     }
    //     else if (statuses == bset().set(4))
    //     {
    //         if (NNBranch == 1)
    //         {
    //             return 8;
    //         }
    //         else
    //         {
    //             return 9;
    //         }
    //     }
    //     else
    //     {
    //         printMsg("ERROR - status\n");
    //         return 11;
    //     }
    //     return 0;
    // }

    // create ambi and insert to it if there is none
    // else, add to suitable ambi and turn the ambi to leaf if the prioroty >stay threshold
    // else create new ambi
    inline size_t insertBranch(seq_type signature, vector<size_t> &insertionList, size_t idx, size_t branch)
    {
        if (ambiLinks[branch].size() == 0)
        {
            return createAmbiNode(signature, insertionList, branch, idx);
        }
        size_t dest = 0;
        size_t ambi_split = 0;
        for (size_t ambi : ambiLinks[branch])
        {
            // preload sig into matrices to check for priority
            addSigToMatrix(ambi, signature);
            updatePriority(ambi);
            // release preloaded matrix
            matrices[ambi].pop_back();

            printMsg(">>ambi %zu, %.2f\n", ambi, priority[ambi]);

            if (priority[ambi] >= stay_threshold)
            {
                dest = stayNode(signature, insertionList, idx, ambi);
                printMsg("@@ambi %zu, %.2f\n", dest, priority[ambi]);
                isAmbiNode[dest] = 0;
                updateParentMean(dest);
                //? also remove from ambi link;
                // ambiLinks[branch].clear();
                removeVecValue(ambiLinks[branch], ambi);
                return dest;
            }
            else if (priority[ambi] < split_threshold)
            {
                ambi_split++;
            }
            else
            {
                dest = ambi;
            }
        }
        if (ambi_split == ambiLinks[branch].size())
        {
            return createAmbiNode(signature, insertionList, branch, idx);
        }
        else
        {
            return stayNode(signature, insertionList, idx, dest);
        }
    }

    // rearrange grandchildren to their best parent
    inline size_t recluster(size_t node)
    {
        printMsg("reclustering %zu\n", node);
        bool moved = false;
        for (size_t child : childLinks[node])
        {
            if (childCounts[child] - ambiLinks[child].size() <= 1)
            {
                continue;
            }
            for (size_t grandchild : childLinks[child])
            {
                if (isAmbiNode[grandchild])
                {
                    continue;
                }
                double max_simimlarity = 0;
                size_t dest = child;
                for (size_t t_parent : childLinks[node])
                {
                    double similarity = calcSimilarity(means[t_parent], means[grandchild]);
                    if (similarity > max_simimlarity)
                    {
                        max_simimlarity = similarity;
                        dest = t_parent;
                    }
                }

                if (dest != child)
                {
                    moveParent(grandchild, dest);
                    moved = true;
                    for (size_t ambi : ambiLinks[child])
                    {
                        moveParent(ambi, dest);
                    }
                }
            }
        }

        if (!moved)
        {
            return 0;
        }

        // tidy
        for (size_t child : childLinks[node])
        {
            if (childCounts[child] == 0)
            {
                deleteNode(child);
            }
            else if (childCounts[child] == 1)
            {
                size_t grandchild = childLinks[child][0];
                size_t idx = getNodeIdx(child);

                childLinks[node][idx] = grandchild;
                parentLinks[grandchild] = node;
                clearNode(child);
            }
            else
            {
                updateNodeMean(child);
            }
        }
        updateParentMean(node);
        return 1;
    }

    // dest cannot be super or root
    // in the case of stay size = 1
    // dest = stay[0]
    inline size_t checkRootP(seq_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        vector<size_t> stay_branch, stay_leaf;
        vector<size_t> NN_leaves;
        vector<size_t> NN_branches;
        vector<size_t> stay_super;

        std::bitset<8> statuses;
        size_t dest = 0;
        double max_similarity = 0;
        size_t dest_branch = 0;
        double max_branch_similarity = 0;
        size_t branchCount = 0;
        size_t dest_super = 0;
        double max_super_similarity = 0;
        size_t superCount = 0;
        size_t mismatch = 0;
        // combineChild(node, insertionList);

        // check the other children
        for (size_t child : childLinks[node])
        {
            // skip ambiNode
            if (isAmbiNode[child])
            {
                // printMsg("Ambi: %zu\n", child);
                continue;
            }
            if (isSuperNode[child])
            {
                superCount++;
            }
            else if (isBranchNode[child])
            {
                branchCount++;
            }

            double similarity = 0;
            size_t status = 0;
            std::tie(status, similarity) = similarityStatus(child, signature);
            switch (status)
            {
            case STAY_F:
                // stay.push_back(child);
                if (isSuperNode[child])
                {
                    stay_super.push_back(child);
                    if (similarity > max_super_similarity)
                    {
                        max_super_similarity = similarity;
                        dest_super = child;
                    }
                    statuses.set(5);
                }
                else if (isBranchNode[child])
                {
                    stay_branch.push_back(child);
                    if (similarity > max_branch_similarity)
                    {
                        max_branch_similarity = similarity;
                        dest_branch = child;
                    }
                    statuses.set(3);
                }
                else
                {
                    stay_leaf.push_back(child);
                    if (similarity > max_similarity)
                    {
                        max_similarity = similarity;
                        dest = child;
                    }
                    // printMsg("child %zu\n", dest);
                    statuses.set(1);
                }

                break;
            case SPLIT_F:
                // mismatch.push_back(child);
                mismatch++;
                break;
            case NN_LEAVE_F:
                NN_leaves.push_back(child);
                statuses.set(2);
                break;
            default:
                NN_branches.push_back(child);
                statuses.set(4);
                break;
            }
        }

        if (statuses.none())
        {
            return createNode(signature, insertionList, node, idx);
        }

        size_t t_parent, t_branch;
        size_t setbits = statuses.count();
        if (setbits == 1)
        {
            size_t pos;
            for (pos = 0; pos < statuses.size(); ++pos)
            {
                if (statuses.test(pos))
                {
                    break;
                }
            }
            switch (pos)
            {
            case 1:
                printMsg("S1 or SM Leaf\n");
                if (stay_leaf.size() > 1)
                {
                    createBranch(node, insertionList, stay_leaf);
                }
                return stayNode(signature, insertionList, idx, dest);
                break;
            case 2:
                printMsg("N1 or NM Leaf\n");
                t_parent = createBranch(node, insertionList, NN_leaves);
                return createAmbiNode(signature, insertionList, t_parent, idx);
                break;
            case 3:
                printMsg("S1 or SM Branch\n");
                if (stay_branch.size() > 1)
                {
                    t_parent = createSuper(node, insertionList, stay_branch);
                    recluster(t_parent);
                }
                return insertBranch(signature, insertionList, idx, dest_branch);
                break;
            case 4:
                printMsg("N1 or NM Branch %zu\n", mismatch);
                if (NN_branches.size() > 1)
                {
                    if (mismatch == 0)
                    {
                        t_parent = node;
                    }
                    else
                    {
                        t_parent = createSuper(node, insertionList, NN_branches);
                    }
                    dest = createUniBranch(t_parent, insertionList, signature, idx);
                    recluster(t_parent);
                    return dest;
                }
                else
                {
                    return insertBranch(signature, insertionList, idx, NN_branches[0]);
                }
                break;
            case 5:
                printMsg("S1 or SM Super\n");
                if (stay_super.size() > 1)
                {
                    t_parent = createSuper(node, insertionList, stay_super);
                    recluster(t_parent);
                }
                return stayNode(signature, insertionList, idx, dest_super);
                break;
            default:
                printMsg("ERROR 1 bit!!\n");
                break;
            }
        }
        else if (setbits == 2)
        {
            size_t tempbits = (statuses & bitset<8>("00101010")).count();
            // only stay
            if (tempbits == 2)
            {
                if (statuses.test(5))
                {
                    if (statuses.test(1))
                    {
                        printMsg("Stay Leaf and Super\n");
                        dest = stayNode(signature, insertionList, idx, dest);
                        relocate(stay_super, stay_leaf);
                        return dest;
                    }
                    else
                    {
                        printMsg("Stay Branch and Super\n");
                        dest = insertBranch(signature, insertionList, idx, dest_branch);
                        relocate(stay_super, stay_branch);
                        return dest;
                    }
                }
                else
                {
                    printMsg("Stay Leaf and Branch\n");
                    dest = stayNode(signature, insertionList, idx, dest);
                    relocate(stay_branch, stay_leaf);
                    return dest;
                }
            }
            else if (tempbits == 1)
            {
                printMsg("NN Branch and ");
                if (statuses.test(1))
                {
                    printMsg("Stay Leaf\n");
                    dest = stayNode(signature, insertionList, idx, dest);
                    t_parent = createSuper(node, insertionList, NN_branches);
                    createBranch(t_parent, insertionList, stay_leaf);
                    recluster(t_parent);
                    return dest;
                }
                else if (statuses.test(3))
                {
                    printMsg("Stay Branch\n");
                    dest = insertBranch(signature, insertionList, idx, dest_branch);
                    t_parent = createSuper(node, insertionList, stay_branch);
                    for (size_t b : NN_branches)
                    {
                        moveParent(b, t_parent);
                    }
                    recluster(t_parent);
                    return dest;
                }
                else if (statuses.test(5))
                {
                    printMsg("Stay Super\n");
                    relocate(stay_super, NN_leaves);
                    if (stay_super.size() > 1)
                    {
                        t_parent = createSuper(node, insertionList, stay_super);
                        // recluster(t_parent);
                        // return stayNode(signature, insertionList, idx, t_parent);
                    }

                    return stayNode(signature, insertionList, idx, dest_super);
                }
            }
            else
            {
                printMsg("NN Leaf and NN Branch\n");
                t_branch = createBranch(node, insertionList, NN_leaves);
                dest = createNode(signature, insertionList, t_branch, idx);
                t_parent = createSuper(node, insertionList, NN_branches);
                moveParent(t_branch, t_parent);
                recluster(t_parent);
                return dest;
            }
        }
        else if (setbits == 3)
        {
            if (statuses.test(5))
            {
                printMsg("Stay Super and ");
                if (statuses.test(1))
                {
                    printMsg("Stay Leaf and ");
                    dest = stayNode(signature, insertionList, idx, dest);
                    if (statuses.test(2))
                    {
                        printMsg("NN Leaf\n");
                        relocate(stay_super, stay_leaf);
                        relocate(stay_super, NN_leaves);
                    }
                    else if (statuses.test(3))
                    {
                        printMsg("Stay Branch\n");
                        relocate(stay_branch, stay_leaf);
                        relocate(stay_super, stay_leaf);
                    }
                    else
                    {
                        printMsg("NN branch\n");
                        relocate(stay_super, stay_leaf);
                        relocate(stay_super, NN_branches);
                    }
                    if (stay_super.size() > 1)
                    {
                        t_parent = createSuper(node, insertionList, stay_super);
                    }
                    else
                    {
                        t_parent = stay_super[0];
                    }
                    recluster(t_parent);
                    return dest;
                }
                else if (statuses.test(2))
                {
                    printMsg("NN Leaf and ");
                    if (statuses.test(3))
                    {
                        printMsg("Stay Branch\n");
                        dest = insertBranch(signature, insertionList, idx, dest_branch);
                        relocate(stay_branch, NN_leaves);

                        if (stay_super.size() > 1)
                        {
                            t_parent = createSuper(node, insertionList, stay_super);
                            t_branch = createSuper(t_parent, insertionList, stay_branch);
                            recluster(t_parent);
                        }
                        else
                        {
                            relocate(stay_super, stay_branch);
                        }
                        return dest;
                    }
                    else
                    {
                        //?
                        printMsg("NN branch\n");
                        t_branch = createBranch(node, insertionList, NN_leaves);
                        dest = createNode(signature, insertionList, t_branch, idx);

                        moveParent(t_branch, dest_super);
                        relocate(stay_super, NN_branches);
                        if (stay_super.size() > 1)
                        {
                            t_parent = createSuper(node, insertionList, stay_super);
                            recluster(t_parent);
                        }
                        return dest;
                    }
                }
                else
                {
                    printMsg("Stay Branch and NN branch\n");
                    dest = insertBranch(signature, insertionList, idx, dest_branch);
                    relocate(stay_super, stay_branch);
                    relocate(stay_super, NN_branches);

                    if (stay_super.size() > 1)
                    {
                        t_parent = createSuper(node, insertionList, stay_super);
                    }
                    else
                    {
                        t_parent = stay_super[0];
                    }
                    recluster(t_parent);
                    return dest;
                }
            }
            else
            {
                if (statuses.test(1))
                {
                    printMsg("Stay Leaf and ");
                    if (statuses.test(2))
                    {
                        printMsg("NN Leaf and ");
                        if (statuses.test(3))
                        {
                            printMsg("Stay Branch\n");
                            dest = stayNode(signature, insertionList, idx, dest);
                            relocate(stay_branch, stay_leaf);
                            relocate(stay_branch, NN_leaves);
                            if (stay_branch.size() > 1)
                            {
                                t_parent = createSuper(node, insertionList, NN_branches);
                                recluster(t_parent);
                            }
                            return dest;
                        }
                        else
                        {
                            printMsg("NN Branch\n");
                            dest = stayNode(signature, insertionList, idx, dest);
                            if (NN_branches.size() > 1)
                            {
                                t_parent = createSuper(node, insertionList, NN_branches);
                                t_branch = createBranch(t_parent, insertionList, NN_leaves);
                                for (size_t l : stay_leaf)
                                {
                                    moveParent(l, t_branch);
                                }
                                recluster(t_parent);
                            }
                            else
                            {
                                relocate(NN_branches, stay_leaf);
                                relocate(NN_branches, NN_leaves);
                            }
                            return dest;
                        }
                    }
                    else
                    {
                        printMsg("Stay Branch and NN branch\n");
                        dest = stayNode(signature, insertionList, idx, dest);
                        relocate(stay_branch, stay_leaf);
                        t_parent = createSuper(node, insertionList, NN_branches);
                        for (size_t b : stay_branch)
                        {
                            moveParent(b, t_parent);
                        }
                        recluster(t_parent);
                        return dest;
                    }
                }
                else
                {
                    printMsg("NN Leaf and Stay Branch and NN branch\n");
                    dest = insertBranch(signature, insertionList, idx, dest_branch);
                    relocate(stay_branch, NN_leaves);
                    t_parent = createSuper(node, insertionList, NN_branches);
                    for (size_t b : stay_branch)
                    {
                        moveParent(b, t_parent);
                    }
                    recluster(t_parent);
                    return dest;
                }
            }
        }
        else if (setbits == 4)
        {
            if (statuses.test(1))
            {
                printMsg("Stay Leaf and ");
                if (statuses.test(2))
                {
                    printMsg("NN Leaf and ");
                    if (statuses.test(3))
                    {
                        printMsg("Stay Branch and ");
                        if (statuses.test(4))
                        {
                            printMsg("NN Branch\n");
                            dest = stayNode(signature, insertionList, idx, dest);
                            t_branch = createBranch(node, insertionList, stay_leaf);
                            for (size_t l : NN_leaves)
                            {
                                moveParent(l, t_branch);
                            }
                            stay_branch.push_back(t_branch);
                            t_parent = createSuper(node, insertionList, stay_branch);
                            for (size_t b : NN_branches)
                            {
                                moveParent(b, t_parent);
                            }
                            recluster(t_parent);
                            return dest;
                        }
                        else
                        {
                            printMsg("Stay Super\n");
                            dest = stayNode(signature, insertionList, idx, dest);
                            relocate(stay_branch, stay_leaf);
                            relocate(stay_branch, NN_leaves);
                            relocate(stay_super, stay_branch);
                            if (stay_super.size() > 1)
                            {
                                t_parent = createSuper(node, insertionList, stay_super);
                            }
                            else
                            {
                                t_parent = stay_super[0];
                            }
                            recluster(t_parent);
                            return dest;
                        }
                    }
                    else
                    {
                        printMsg("NN Branch and Stay Super\n");
                        dest = stayNode(signature, insertionList, idx, dest);
                        t_branch = createBranch(node, insertionList, stay_leaf);
                        for (size_t l : NN_leaves)
                        {
                            moveParent(l, t_branch);
                        }
                        NN_branches.push_back(t_branch);
                        relocate(stay_super, NN_branches);
                        if (stay_super.size() > 1)
                        {
                            t_parent = createSuper(node, insertionList, stay_super);
                        }
                        else
                        {
                            t_parent = stay_super[0];
                        }
                        recluster(t_parent);
                        return dest;
                    }
                }
                else
                {
                    printMsg("Stay Branch and NN Branch and Stay Super\n");
                    dest = stayNode(signature, insertionList, idx, dest);
                    relocate(stay_branch, stay_leaf);
                    relocate(stay_super, stay_branch);
                    relocate(stay_super, NN_branches);
                    if (stay_super.size() > 1)
                    {
                        t_parent = createSuper(node, insertionList, stay_super);
                    }
                    else
                    {
                        t_parent = stay_super[0];
                    }
                    recluster(t_parent);
                    return dest;
                }
            }
            else
            {
                printMsg("NN Leaf and Stay Branch and NN Branch and Stay Super\n");
                dest = insertBranch(signature, insertionList, idx, dest_branch);
                relocate(stay_branch, NN_leaves);
                relocate(stay_super, stay_branch);
                relocate(stay_super, NN_branches);
                if (stay_super.size() > 1)
                {
                    t_parent = createSuper(node, insertionList, stay_super);
                }
                else
                {
                    t_parent = stay_super[0];
                }
                recluster(t_parent);
                return dest;
            }
        }

        printMsg("ERROR!!\n");
        return 0;

        // // cerr << statuses << "\n";
        // // cerr << "Status: " << getMethod(statuses, (NN_branches.size() > 1 || NN_leaves.size() > 1)) << "\n";

        // size_t method = getMethod(statuses);//, stay_super.size(), superCount);

        // // return nearest branch if not staying in any leaf
        // if (dest == 0)
        // {
        //     dest = dest_branch;
        // }
        // return make_pair(method, dest);
    }

    inline size_t ambiBranch(seq_type signature, vector<size_t> &insertionList, size_t idx, size_t node, vector<size_t> NN)
    {
        size_t t_parent = createParent(node, insertionList);
        printMsg("*** ambiBranch\n");
        for (size_t child : NN)
        {
            moveParent(child, t_parent);
        }

        addSigToMatrix(node, signature);

        return createAmbiNode(signature, insertionList, t_parent, idx);
    }

    inline size_t createUniBranch(size_t node, vector<size_t> &insertionList, seq_type signature, size_t idx)
    {
        size_t t_parent = createParent(node, insertionList);
        size_t dest = createNode(signature, insertionList, t_parent, idx);
        addSigToMatrix(node, signature);
        updateParentMean(t_parent);
        return dest;
    }

    // create a new branch under "node"
    // then move candidates to the new branch
    inline size_t createBranch(size_t node, vector<size_t> &insertionList, vector<size_t> candidates)
    {
        size_t t_parent = createParent(node, insertionList);

        for (size_t child : candidates)
        {
            moveParent(child, t_parent);
        }
        addSigToMatrix(node, means[t_parent]);
        updateParentMean(t_parent);
        return t_parent;
    }

    // create a new super under "node"
    // then move candidates to the new super
    inline size_t createSuper(size_t node, vector<size_t> &insertionList, vector<size_t> candidates)
    {
        size_t t_parent = createBranch(node, insertionList, candidates);
        isSuperNode[t_parent] = 1;
        return t_parent;
    }

    //?
    inline void relocate(vector<size_t> parents, vector<size_t> children)
    {
        printMsg("Relocating...\n");
        for (size_t child : children)
        {
            double max_simimlarity = split_threshold;
            size_t dest = 0;
            for (size_t parent : parents)
            {
                double similarity = calcSimilarity(means[child], means[parent]);
                if (similarity > max_simimlarity)
                {
                    max_simimlarity = similarity;
                    dest = parent;
                }
            }
            //? maybe find the best level to move
            if (dest != 0)
            {
                moveParent(child, dest);
            }
        }
    }

    inline size_t tt(seq_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        size_t t_parent = 0;
        size_t temp_dest = 0;
        size_t ambi_split = 0;
        double max_similarity = 0;

        size_t method, dest;

        // size_t node = root;

        return checkRootP(signature, insertionList, idx, node);

        // switch (method)
        // {
        // case 0:
        //     printMsg("#mismatch\n");
        //     return createNode(signature, insertionList, node, idx);
        //     break;
        // case 1:
        //     printMsg("#stay in leaf %zu\n", dest);
        //     return stayNode(signature, insertionList, idx, dest);
        //     break;
        // case 2:
        //     printMsg("#NN with leaf: %zu\n", NN_leaves.size());
        //     t_parent = node;
        //     if (ambiLinks[t_parent].size() == 0)
        //     {
        //         t_parent = createBranch(node, insertionList, NN_leaves);
        //         dest = createAmbiNode(signature, insertionList, t_parent, idx);
        //     }
        //     else
        //     {
        //         for (size_t i = 0; i < ambiLinks[t_parent].size(); i++)
        //         {
        //             size_t ambi = ambiLinks[t_parent][i];
        //             // dest = ambiLinks[t_parent][0];
        //             // dest = ambi;

        //             // dest = stayNode(signature, insertionList, idx, dest);
        //             // updatePriority(dest);

        //             // preload sig into matrices to check for priority
        //             addSigToMatrix(ambi, signature);
        //             double temp_priority = priority[ambi];
        //             updatePriority(ambi);
        //             // release preloaded matrix
        //             matrices[ambi].pop_back();

        //             printMsg(">>ambi %zu, %.2f\n", ambi, priority[ambi]);

        //             if (priority[ambi] > stay_threshold)
        //             {
        //                 dest = stayNode(signature, insertionList, idx, ambi);
        //                 printMsg("@@ambi %zu, %.2f\n", dest, priority[dest]);
        //                 isAmbiNode[ambi] = 0;
        //                 updateParentMean(dest);
        //                 //? also remove from ambi link;
        //                 // ambiLinks[t_parent].clear();
        //                 ambiLinks[t_parent].erase(ambiLinks[t_parent].begin() + i);
        //                 return dest;
        //             }
        //             else if (priority[ambi] < split_threshold)
        //             {

        //                 printMsg("??ambi split %zu, %.2f\n", ambi, priority[ambi]);
        //                 ambi_split++;
        //                 priority[ambi] = temp_priority;
        //                 // dest = createAmbiNode(signature, insertionList, t_parent, idx);
        //             }
        //         }
        //         if (ambi_split == ambiLinks[t_parent].size())
        //         {
        //             dest = createAmbiNode(signature, insertionList, t_parent, idx);
        //         }
        //         else
        //         {
        //             dest = stayNode(signature, insertionList, idx, ambiLinks[t_parent][0]);
        //         }
        //     }

        //     return dest;
        //     break;
        // case 3:
        //     dest = stay_super[0];
        //     printMsg("#stay with ONE branch: %zu\n", dest);
        //     return tt(signature, insertionList, idx, dest);
        //     break;
        // case 4:
        //     printMsg("#stay with ALL branch: %zu\n", dest);
        //     return tt(signature, insertionList, idx, dest);
        //     break;
        // case 5:
        //     printMsg("#stay with SOME branch: %zu\n", dest);
        //     t_parent = createSuper(node, insertionList, stay_supers);
        //     return tt(signature, insertionList, idx, dest);
        //     break;
        // case 6:
        //     printMsg("#stay with branch(s) and NN with %zu leaves: %zu\n", dest, NN_leaves.size());
        //     dest = tt(signature, insertionList, idx, dest);
        //     relocate(NN_branches, NN_leaves);
        //     return dest;
        //     break;
        // case 7:
        //     printMsg("ERROR 7!!!\n");
        //     break;
        // case 8:
        //     printMsg("ERROR 8!!!\n");
        //     break;
        // case 9:
        //     t_parent = createSuper(node, insertionList, stay_supers);
        //     break;
        // default:
        //     printMsg("ERROR!!!\n");
        //     break;
        // }

        // return 0;
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
        size_t node = 0;
        // size_t node = tt_rootSplit(signature, insertionList, idx);
        printMsg("inserted %zu at %zu\n\n", idx, node);
        return node;
    }

    inline size_t insertSplitRoot(seq_type signature, vector<size_t> &insertionList, size_t idx)
    {
        // size_t node = insertSplit(signature, insertionList, idx);
        size_t node = insert(signature, insertionList, idx);
        if (childCounts[root] > tree_order)
        {
            // forceSplitRoot(insertionList, false);

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

    inline size_t searchBest(seq_type signature, size_t node = 0)
    {
        size_t dest = 0;
        double max_similarity = 0;

        for (size_t child : childLinks[node])
        {
            double similarity = calcSimilarity(means[child], signature);
            printMsg(" (%zu, %.2f) ", child, calcSimilarity(means[child], signature));

            if (similarity > max_similarity)
            {
                max_similarity = similarity;
                dest = child;
            }
        }
        printMsg("\n>>>%zu\n", dest);

        if (isBranchNode[dest])
        {
            return searchBest(signature, dest);
        }
        else
        {
            return dest;
        }
    }

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
        // if (!isAmbiNode[node])
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
        // if (!isAmbiNode[node])
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

    // update priority of ambi nodes, for debugging only
    void prepareAmbi(size_t node = 0)
    {
        for (size_t child : childLinks[node])
        {
            if (isBranchNode[child])
            {
                prepareAmbi(child);
            }
            else if (isAmbiNode[child])
            {
                if (matrices[child].size() > 1)
                {
                    updatePriority(child);
                    if (priority[child] > stay_threshold)
                    {
                        isAmbiNode[child] = 0;
                        //? also remove from ambi link;
                    }
                }
            }
        }
    }

    size_t reinsert(seq_type signature, size_t idx)
    {
        // size_t node = search(signature);
        size_t node = searchBest(signature);
        // size_t node = search2(signature);
        seqIDs[node].push_back(idx);
        addSigToMatrix(node, signature);
        return node;
        // return parentLinks[node];
    }

    void getEmptyLeaves(vector<size_t> &leaves, size_t node = 0)
    {
        for (size_t child : childLinks[node])
        {
            if (isBranchNode[child])
            {
                // getEmptyLeaves(leaves, child);

                // unify tight branch
                if (isRootNode[child] || isSuperNode[child])
                {
                    getEmptyLeaves(leaves, child);
                }
                else if (priority[child] > stay_threshold)
                {
                    for (size_t i = 1; i < childCounts[child]; i++)
                    {
                        leaves.push_back(childLinks[child][i]);
                    }
                }
                else
                {
                    getEmptyLeaves(leaves, child);
                }
            }
            else if (isAmbiNode[child] | seqIDs[child].size() <= singleton) // remove singleton
            {
                leaves.push_back(child);
            }
            else
            {
                updatePriority(child);
                updateParentMean(node);
                // size_t size = seqIDs[child].size();
                // if (size <= minClusSize)
                // {
                //     minClusSize = size;
                // }
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

    void trim()
    {
        std::set<size_t> branches;

        vector<size_t> leaves;
        getEmptyLeaves(leaves);
        // singleton = minClusSize + 1;
        printMsg("\nTrimming, singleton size = %zu\n", singleton);
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