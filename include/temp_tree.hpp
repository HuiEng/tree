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
size_t singleton = 2;
size_t tree_order = 5;

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

    // delete a child from its parent, need to format the child separately
    inline void deleteNode(size_t node)
    {
        // fprintf(stderr, "deleting %zu\n", node);
        size_t parent = parentLinks[node];

        int idx = getNodeIdx(node);

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
            updateNodeMean(node);
            updatePriority(node);
            node = parentLinks[node];
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

    // RMSD
    double calcDistortion(size_t node)
    {
        vector<seq_type> temp_matrix = getNonAmbiMatrix(node);
        if (temp_matrix.size() <= 1)
        {
            return 0;
        }

        double sumSquaredistance = 0;
        vector<cell_type> meanSig = getMinimiseSet(createMeanSig(temp_matrix))[0];

        for (seq_type signature : temp_matrix)
        {
            vector<cell_type> signatureData = getMinimiseSet(signature)[0];
            double distance = calcJaccardGlobal(meanSig, signatureData);
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

    // union mean of children
    inline void updatePriority(size_t node)
    {
        priority[node] = calcDistortion(node);
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
        //     fprintf(stderr, "*** stayAmbi %zu\n", ambi);
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

        //         fprintf(stderr, "***???\n");
        //         return createAmbiNode(signature, insertionList, node, idx);
        //     }
        //     else
        //     {
        //         return stayNode(signature, insertionList, idx, ambiLinks[node][0]);
        //     }
        // }
        // else
        // {
        //     fprintf(stderr, "***---\n");
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
        // size_t dest = getNewNodeIdx(insertionList);
        // parentLinks[dest] = node;
        // childLinks[node].push_back(dest);
        // seqIDs[dest].push_back(idx);
        ambiLinks[node].push_back(dest);
        isAmbiNode[dest] = 1;
        return dest;
    }

    size_t findFurthest(size_t node, seq_type mean0)
    {
        double min_similarity = 1;
        size_t candidate = 0;
        // seq_type mean0 = means[childLinks[node][0]];
        for (size_t i = 0; i < childCounts[node]; i++)
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
            fprintf(stderr, "??cannot find candidate %zu, %f\n", node, priority[node]);
        }
        return candidate;
    }

    size_t forceSplitRoot(vector<size_t> &insertionList, size_t node = 0, bool unpack = false)
    {
        // printTreeJson(stderr);
        fprintf(stderr, ">> Force split root %zu\n", node);
        // rootNodes.resize(1);
        // size_t node = root;

        // unpack previous subtree
        if (unpack)
        {
            vector<size_t> tempChildLinks = childLinks[node];
            for (size_t subtree : tempChildLinks)
            {
                if (isRootNode[subtree])
                {
                    for (size_t child : childLinks[subtree])
                    {
                        moveParent(child, node, false);
                    }
                    deleteNode(subtree);
                }
            }
        }

        //? assume the first child is never ambi
        vector<size_t> temp_centroids;
        vector<vector<size_t>> clusters;
        //?
        vector<size_t> temp;
        clusters.push_back(temp);
        vector<size_t> t;
        clusters.push_back(t);

        // size_t candidate = findFurthest(node, createMeanSig(getNonAmbiMatrix(node)));

        size_t candidate = findFurthest(node, means[node]);
        temp_centroids.push_back(candidate);
        temp_centroids.push_back(findFurthest(node, means[candidate]));

        // fprintf(stderr, "??something is wrong %zu,%zu,%zu \n", matrices[node].size(), childLinks[node].size(), childCounts[node]);

        for (size_t n = 0; n < matrices[node].size(); n++)
        {
            size_t child = childLinks[node][n];
            size_t dest = 0;
            double max_similarity = 0;
            for (size_t i = 0; i < temp_centroids.size(); i++)
            {
                double similarity = calcSimilarity(means[temp_centroids[i]], matrices[node][n]);
                fprintf(stderr, "%zu, %f\n", temp_centroids[i], similarity);
                if (similarity > max_similarity)
                {
                    max_similarity = similarity;
                    dest = i;
                }
            }
            fprintf(stderr, "--- %zu goes to %zu\n", child, dest);
            clusters[dest].push_back(child);
        }

        if (clusters[0].size() <= 1 || clusters[1].size() <= 1)
        {
            fprintf(stderr, "??something is wrong cannot split evenly\n");
            return 0;
        }

        for (size_t i = 0; i < temp_centroids.size(); i++)
        {
            size_t t_parent = createParent(node, insertionList);
            isRootNode[t_parent] = 1;
            for (size_t n = 0; n < clusters[i].size(); n++)
            {
                moveParent(clusters[i][n], t_parent);
                // fprintf(stderr, "--- %zu,%zu\n", clusters[i][n], t_parent);
            }
            updateNodeMean(t_parent);
            updatePriority(t_parent);
            addSigToMatrix(node, means[t_parent]);
        }
        updateParentMean(node);
        // printTreeJson(stderr);

        return 1;
    }

    inline size_t similarityStatus(seq_type sig1, seq_type sig2, double offset = 0)
    {
        // check the other children
        double local_stay_t = stay_threshold - offset;

        // double offset = 1.2;
        // // size_t rank = findLevel(child);
        // for (size_t t = 0; t < rank; t++)
        // {
        //     local_stay_t = local_stay_t + stay_threshold * offset;
        // }
        // // double NN_t = local_stay_t + stay_threshold * offset;

        double local_split_t = split_threshold;
        double similarity = calcSimilarity(sig1, sig2);

        fprintf(stderr, "%.2f", similarity);

        if (similarity >= local_stay_t)
        {
            fprintf(stderr, ": STAY>\n");
            return STAY_F;
        }
        else if (similarity <= local_split_t)
        {
            fprintf(stderr, ": SPLIT>\n");
            return SPLIT_F;
        }
        else
        {
            fprintf(stderr, ": NN>\n");
            return NN_F;
        }
    }

    inline size_t similarityStatus(size_t child, size_t rank, seq_type signature)
    {
        fprintf(stderr, "<%zu, ", child);
        size_t status;

        if (isRootNode[child])
        {
            // root can be stay or split only
            double similarity = calcOverlap(signature, means[child]);
            fprintf(stderr, "%.2f> r\n ", similarity);
            if (similarity >= stay_threshold)
            {
                return STAY_F;
            }
            else
            {
                return SPLIT_F;
            }
        }
        else if (isBranchNode[child])
        {
            // double offset = 0.1 * rank;
            double offset = 0;
            status = similarityStatus(means[child], signature, offset);
        }
        else
        {
            status = similarityStatus(means[child], signature);
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
            size_t rank = findLevel(child);
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
        updateParentMean(node);
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
            else if (isRootNode[dest])
            {
                fprintf(stderr, "?with NN leaves & is subtree\n");
                // ignore leaves for now
                return tt_root(signature, insertionList, idx, dest);
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

        return createNode(signature, insertionList, node, idx);
    }

    inline size_t forceSplitSubtree(vector<size_t> &insertionList, size_t node)
    {
        updatePriority(node);
        if (childCounts[node] <= tree_order && priority[node] < 0.5)
        {
            return 0;
        }

        fprintf(stderr, ">> Force split subtree %zu\n", node);
        forceSplitRoot(insertionList, node, false);
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
            fprintf(stderr, "is subtree %zu\n", node);
            forceSplitSubtree(insertionList, node);
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
            updateParentMean(node);
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
                return stayAmbi(signature, insertionList, idx, node);
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
                return stayAmbi(signature, insertionList, idx, node);
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
        return node;
    }

    inline size_t insertSplitRoot(seq_type signature, vector<size_t> &insertionList, size_t idx)
    {
        size_t node = insert(signature, insertionList, idx);
        if (childCounts[root] > tree_order)
        {
            forceSplitRoot(insertionList);

            //?
            vector<size_t> subtrees = childLinks[root];
            for (size_t t_parent : subtrees)
            {
                forceSplitSubtree(insertionList, t_parent);
            }
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
    //         fprintf(stderr, " <%zu,%.2f> ", child, similarity);

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
            fprintf(stderr, " <%zu,%.2f> ", child, similarity);
            fprintf(stderr, " (%.2f, %.2f, %.2f) ", calcSimilarity(means[child], signature), calcOverlap(signature, means[child]), priority[child]);

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
                fprintf(stderr, "\n>%zu\n",best_child);
                return search(signature, idx, best_child);
            }
            else
            {
                return best_child;
            }
        }
        else
        {

            fprintf(stderr, "\n*** here\n");
            double best_dest_similarity = 0;
            for (size_t child : candidates)
            {
                fprintf(stderr, "child %zu\n", child);
                size_t leaf = child;
                double similarity = best_similarity;
                if (isBranchNode[child])
                {
                    leaf = search(signature, idx, child);
                    similarity = calcSimilarity(means[child], signature);
                    similarity += calcOverlap(signature, means[leaf]);
                    fprintf(stderr, "> leaf %zu\n", leaf);
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
            if (seqIDs[i].size() < singleton)
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