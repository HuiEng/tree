// There are two kinds of ktree nodes- branch nodes and leaf nodes
// Both contain a signature matrix, plus their own signature
// (the root node signature does not matter and can be blank)
// Branch nodes then contain 'order' links to other nodes
// Leaf nodes do not.
// However, as leaf nodes may become branch nodes, we allocate
// the space anyway.
// As the space to be used is determined at runtime, we use
// parallel arrays, not structs

#ifndef INCLUDE_self_tree_HPP
#define INCLUDE_self_tree_HPP

#include <omp.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include "bloom_filter.hpp"
#include "read.hpp"

typedef vector<vector<cell_type>> seq_type;
double split_threshold = 1;
double stay_threshold = 0;
using namespace std;

// typedef unsigned char cell_type;
bloom_parameters parameters;
size_t signatureSize; // Signature size (depends on element in BF, obtained while read binary)
size_t partree_capacity = 100;

void toBinary(FILE *stream, vector<cell_type> sig)
{
    //   fprintf(stderr, "%p: ", sig);
    for (size_t i = 0; i < signatureSize; i++)
    {
        int binary[bits_per_char];
        for (int n = 0; n < bits_per_char; n++)
            binary[bits_per_char - 1 - n] = (sig[i] >> n) & 1;

        for (int n = 0; n < bits_per_char; n++)
            fprintf(stream, "%d", binary[n]);
    }
    fprintf(stream, "\n");
}

// print each window in one line
void dbgPrintSignature(FILE *stream, vector<vector<cell_type>> seq)
{
    for (auto window : seq)
    {
        toBinary(stream, window);
    }
    fprintf(stream, "\n");
}

class self_tree
{
public:
    size_t root = 0;                   // # of root node
    vector<size_t> childCounts;        // n entries, number of children
    vector<int> isBranchNode;          // n entries, is this a branch node
    vector<vector<size_t>> childLinks; // n * o entries, links to children
    vector<vector<size_t>> seqIDs;     // n * o entries, links to children
    vector<size_t> parentLinks;        // n entries, links to parents
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
                childLinks.resize(capacity);
            }
            //#pragma omp single
            {
                seqIDs.resize(capacity);
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

    self_tree(size_t capacity)
    {
        reserve(capacity);
        childCounts[root] = 0;
        isBranchNode[root] = 0;
    }

    size_t calcHDwrapper(seq_type shorter, seq_type longer) const
    {
        size_t c = 0;
        // treat tail subseq as mismatch
        for (int w = 0; w < shorter.size(); w++)
        {
            for (size_t i = 0; i < signatureSize; i++)
            {
                c += __builtin_popcountll(shorter[w][i] ^ longer[w][i]);
            }
        }
        // treat tail subseq as mismatch
        for (int w = shorter.size(); w < longer.size(); w++)
        {
            for (size_t i = 0; i < signatureSize; i++)
            {
                c += __builtin_popcountll(longer[w][i]);
            }
        }
        return c;
    }

    size_t calcHD(seq_type a, seq_type b) const
    {
        if (a.size() < b.size())
        {
            return calcHDwrapper(a, b);
        }
        else
        {
            return calcHDwrapper(b, a);
        }
    }

    /*
        // breadth first, printing top down
        void find_leaves(size_t parent, vector<size_t> &leaves)
        {
            for (size_t child : childLinks[parent])
            {
                if (isBranchNode[child])
                {
                    find_leaves(child, leaves);
                }
                else
                {
                    leaves.push_back(child);
                }
            }
        }

        inline size_t traverse(const sig_type *signature) const
        {
            size_t node = root;

            while (isBranchNode[node])
            {
                size_t lowestHD = numeric_limits<size_t>::max();
                size_t lowestHDchild = childLinks[node][0];

                for (size_t i = 0; i < childCounts[node]; i++)
                {
                    size_t child = childLinks[node][i];
                    size_t hd = calcHD(&means[child * signatureSize], signature);
                    if (hd < lowestHD)
                    {
                        lowestHD = hd;
                        lowestHDchild = child;
                    }
                }

                node = lowestHDchild;
            }

            return node;
        }


        void addSigToMatrix(size_t node, const sig_type *sig)
        {
            matri[node].insert(matri[node].end(), sig, sig + signatureSize);
        }

    */

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

    size_t findAncestor(size_t node){
        while(parentLinks[node]!=root){
            node = parentLinks[node];
        }
        return node;
    }
    
    /*
        // lock parent while updating matrix, return parent node for unlocking later
        size_t findParent(size_t node, sig_type *meanSig)
        {
            size_t parent = parentLinks[node];
            // Lock the parent
            omp_set_lock(&locks[parent]);

            size_t idx = numeric_limits<size_t>::max();
            for (size_t i = 0; i < childCounts[parent]; i++)
            {
                if (childLinks[parent][i] == node)
                {
                    idx = i;
                    break;
                }
            }

            // get the wrong parent, try again
            size_t count = 0;
            while (idx == numeric_limits<size_t>::max())
            {
                // unset old lock
                omp_unset_lock(&locks[parent]);

                // find and lock new parent
                parent = parentLinks[node];
                omp_set_lock(&locks[parent]);

                // find idx
                for (size_t i = 0; i < childCounts[parent]; i++)
                {
                    if (childLinks[parent][i] == node)
                    {
                        idx = i;
                        break;
                    }
                }
                count++;
                // after 10 trials, skip update, return parent for unlocking
                if (count > 10)
                {
                    return parent;
                }
            }

            memcpy(&matri[parent][idx * signatureSize], meanSig, sizeof(sig_type) * signatureSize);

            return parent;
        }

        void update(size_t node)
        {
            for (int i = 0; i < childCounts[node]; i++)
            {
                size_t child = childLinks[node][i];
                memcpy(&matri[node][i * signatureSize], &means[child * signatureSize], sizeof(sig_type) * signatureSize);
            }
            if (isBranchNode[parentLinks[node]])
            {
                // fprintf(stdout,"again %zu %zu\n",node, parentLinks[node]);
                update(parentLinks[node]);
            }
        }

        template <class RNG>
        void splitNode(RNG &&rng, size_t node, const sig_type *sig, vector<size_t> &insertionList, size_t link)
        {
            size_t nodeSigs = childCounts[node] + 1;
            vector<sig_type> sigs(nodeSigs * signatureSize);
            //?
            copy(matri[node].begin(), matri[node].end(), sigs.begin());
            memcpy(&sigs[childCounts[node] * signatureSize], sig, sizeof(sig_type) * signatureSize);

            vector<sig_type> meanSigs = createRandomSigs(rng, sigs);
            vector<size_t> clusters(nodeSigs);
            vector<vector<size_t>> clusterLists;

            // dbgPrintSignature(&meanSigs[0]);
            // dbgPrintSignature(&meanSigs[signatureSize]);
            vector<sig_type> temp(2 * signatureSize);

            for (int iteration = 0; iteration < 4; iteration++)
            {
                // fprintf(stderr, "Iteration %d\n", iteration);
                reclusterSignatures(clusters, meanSigs, sigs);
                clusterLists = createClusterLists(clusters);
                memcpy(&temp[0], &meanSigs[0], sizeof(meanSigs));
                meanSigs = createClusterSigs(clusterLists, sigs, nodeSigs);
            }

            // Create the sibling node
            size_t sibling = getNewNodeIdx(insertionList);
            size_t newlyAddedIdx = childCounts[node];

            childCounts[sibling] = clusterLists[1].size();
            isBranchNode[sibling] = isBranchNode[node];
            {
                size_t siblingIdx = 0;
                for (size_t seqIdx : clusterLists[1])
                {

                    // If this is a branch node, relink the child to the new parent
                    if (isBranchNode[sibling])
                    {
                        if (seqIdx < newlyAddedIdx)
                        {
                            childLinks[sibling].push_back(childLinks[node][seqIdx]);
                        }
                        else
                        {
                            childLinks[sibling].push_back(link);
                        }
                        parentLinks[childLinks[sibling][siblingIdx]] = sibling;
                    }
                    else
                    {
                        if (seqIdx < newlyAddedIdx)
                        {
                            childLinks[sibling].push_back(childLinks[node][seqIdx]);
                        }
                        else
                        {
                            childLinks[sibling].push_back(link);
                        }
                    }
                    addSigToMatrix(sibling, &sigs[seqIdx * signatureSize]);
                    siblingIdx++;
                }
            }

            memcpy(&means[node * signatureSize], &meanSigs[0], sizeof(sig_type) * signatureSize);
            memcpy(&means[sibling * signatureSize], &meanSigs[signatureSize], sizeof(sig_type) * signatureSize);

            childCounts[node] = clusterLists[0].size();
            // Fill the current node with the other cluster of signatures
            {
                matri[node].clear();
                size_t nodeIdx = 0;
                for (size_t seqIdx : clusterLists[0])
                {

                    // If this is a branch node, relink the child to the new parent
                    if (isBranchNode[node])
                    {
                        if (seqIdx < newlyAddedIdx)
                        {
                            childLinks[node][nodeIdx] = childLinks[node][seqIdx];
                        }
                        else
                        {
                            childLinks[node][nodeIdx] = link;
                        }
                        parentLinks[childLinks[node][nodeIdx]] = node;
                    }
                    else
                    {
                        if (seqIdx < newlyAddedIdx)
                        {
                            childLinks[node][nodeIdx] = childLinks[node][seqIdx];
                        }
                        else
                        {
                            childLinks[node][nodeIdx] = link;
                        }
                    }
                    addSigToMatrix(node, &sigs[seqIdx * signatureSize]);
                    nodeIdx++;
                }
                if (isBranchNode[node])
                {
                    childLinks[node].resize(childCounts[node]);
                }
            }

            if (!isBranchNode[node])
            {
                childLinks[node].resize(childCounts[node]);
                fprintf(stdout, "splitting at %zu\n", node);

                for (size_t child : childLinks[node])
                {
                    fprintf(stdout, "%zu,", child);
                }
                fprintf(stdout, "\n");

                fprintf(stdout, "sibling %zu\n", sibling);

                for (size_t child : childLinks[sibling])
                {
                    fprintf(stdout, "%zu,", child);
                }
                fprintf(stdout, "\n");

                // childLinks[node].clear();
                // childLinks[sibling].clear();
                // for (size_t i = 0; i < clusters.size(); i++)
                // {
                //     fprintf(stderr,"%zu\n", clusters[i]);
                //     if (clusters[i]==0){
                //         childLinks[node].push_back(i);

                //     }else{
                //         childLinks[sibling].push_back(i);
                //     }
                // }
            }

            // Is this the root level?
            if (node == root)
            {
                // Create a new root node
                size_t newRoot;
                newRoot = getNewNodeIdx(insertionList);

                // Link this node and the sibling to it
                parentLinks[node] = newRoot;
                parentLinks[sibling] = newRoot;

                childCounts[newRoot] = 2;
                isBranchNode[newRoot] = 1;
                childLinks[newRoot].push_back(node);
                childLinks[newRoot].push_back(sibling);
                addSigToMatrix(newRoot, &meanSigs[0 * signatureSize]);
                addSigToMatrix(newRoot, &meanSigs[1 * signatureSize]);

                root = newRoot;

                recalculateSig(root);
            }
            else
            {
                // // First, update the reference to this node in the parent with the new mean
                // size_t parent = findParent(node, &meanSigs[0]);

                // // Connect sibling node to parent
                // parentLinks[sibling] = parent;

                // if (childCounts[parent] + 1 < order)
                // {
                //     childLinks[parent].push_back(sibling);

                //     childCounts[parent]++;
                //     addSigToMatrix(parent, &meanSigs[1 * signatureSize]);

                //     // Update signatures (may change?)
                //     recalculateUp(parent);
                // }
                // else
                // {
                //     // fprintf(stdout, "after\n");
                //     splitNode(rng, parent, &meanSigs[1 * signatureSize], insertionList, sibling);
                // }
                // update(parent);
                // // Unlock the parent
                // omp_unset_lock(&locks[parent]);
            }
            // fprintf(stdout, ">>>\n");
            // for (int i = 0; i < childCounts[root]; i++)
            // {
            //     dbgPrintSignature(&matri[root][i * signatureSize]);
            //     dbgPrintSignature(&means[childLinks[root][i] * signatureSize]);
            //     fprintf(stdout, "\n");
            // }
            // fprintf(stdout, ">>>\n");
            // printTree();
        }

    */

    void printNodeJson(FILE *stream, size_t tnode)
    {
        fprintf(stream, "{\"node\":\"%zu\",\"childCount\":\"%zu\",\"content\":\"*", tnode, seqIDs[tnode].size());
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
            fprintf(stream, "{\"node\":\"%zu\",\"childCount\":\"%zu\",\"content\":\"*", tnode, seqIDs[tnode].size());
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

    inline size_t first_insert(seq_type signature, vector<size_t> &insertionList, size_t idx)
    {
        size_t parent = root;
        size_t node = getNewNodeIdx(insertionList);
        parentLinks[node] = parent;
        childLinks[parent].push_back(node);
        isBranchNode[parent] = 1;
        childCounts[parent]++;

        seqIDs[node].push_back(idx);
        addSigToMatrix(node, signature);
        return node;
    }

    // return if found same, and the destination node
    inline tuple<bool, size_t> traverse(seq_type signature) const
    {
        size_t node = root;
        size_t a = signature.size();

        while (isBranchNode[node])
        {
            // fprintf(stderr, " \n%zu: ", node);
            size_t lowestHD = numeric_limits<size_t>::max();
            size_t lowestHDchild = childLinks[node][0];
            size_t mismatch = 0;

            for (size_t i = 0; i < childCounts[node]; i++)
            {
                size_t child = childLinks[node][i];
                size_t hd = calcHD(matrices[child][0], signature);
                size_t b = matrices[child][0].size();

                // found same, move on with the next seq
                if (hd <= stay_threshold * max(a, b))
                {
                    return make_tuple(true, child);
                }
                else if (hd < lowestHD)
                {
                    lowestHD = hd;
                    lowestHDchild = child;
                }

                // count how many nodes mismatch
                if (hd >= split_threshold * a)
                {
                    mismatch++;
                }

                // fprintf(stderr, " <%zu,%zu> ", child, hd);
            }


            // nothing is close enough, spawn new child under parent
            if (mismatch == childCounts[node])
            {
                return make_tuple(false, node);
            }

            node = lowestHDchild;
        }

        return std::make_tuple(false, node);
    }

    inline size_t insert(seq_type signature, vector<size_t> &insertionList, size_t idx)
    {
        bool stay;
        size_t parent;
        tie(stay, parent) = traverse(signature);
        // fprintf(stderr, "inserting seq %zu at node %zu; stay: %d\n", idx, parent, stay);

        if (!stay)
        {
            size_t node = getNewNodeIdx(insertionList);
            parentLinks[node] = parent;
            childLinks[parent].push_back(node);
            isBranchNode[parent] = 1;
            childCounts[parent]++;
            addSigToMatrix(node, signature);
            seqIDs[node].push_back(idx);
            return node;
        }
        else
        {
            seqIDs[parent].push_back(idx);
            return parent;
        }

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

    // return if found same, and the destination node
    inline size_t search(seq_type signature, size_t idx) const
    {
        size_t node = root;
        size_t a = signature.size();

        while (isBranchNode[node])
        {
            // fprintf(stderr, " \n%zu: ", node);
            size_t lowestHD = numeric_limits<size_t>::max();
            size_t lowestHDchild = childLinks[node][0];
            size_t mismatch = 0;

            for (size_t i = 0; i < childCounts[node]; i++)
            {
                size_t child = childLinks[node][i];
                size_t hd = calcHD(matrices[child][0], signature);
                size_t b = matrices[child][0].size();

                // found same, move on with the next seq
                if (hd <= stay_threshold * max(a, b))
                {
                    return child;
                }
                else if (hd < lowestHD)
                {
                    lowestHD = hd;
                    lowestHDchild = child;
                }

                // count how many nodes mismatch
                if (hd >= split_threshold * a)
                {
                    mismatch++;
                }

                // fprintf(stderr, " <%zu,%zu> ", child, hd);
            }


            //? nothing is close enough, return parent
            if (mismatch == childCounts[node])
            {
                fprintf(stderr, ">%zu mismatch\n", idx);
                return parentLinks[node];
            }

            node = lowestHDchild;
        }

        return node;
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