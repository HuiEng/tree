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

typedef vector<sig_type> seq_type;
using namespace std;

// typedef unsigned char cell_type;
bloom_parameters parameters;
size_t signatureSize; // Signature size (depends on element in BF, obtained while read binary)
size_t partree_capacity = 1000000;

void dbgPrintSignature(const sig_type *sig)
{
    //   fprintf(stderr, "%p: ", sig);
    for (size_t i = 0; i < signatureSize; i++)
    {
        int binary[bits_per_char];
        for (int n = 0; n < bits_per_char; n++)
            binary[bits_per_char - 1 - n] = (sig[i] >> n) & 1;

        for (int n = 0; n < bits_per_char; n++)
            fprintf(stdout, "%d", binary[n]);
    }
    fprintf(stdout, "\n");
}

vector<vector<size_t>> createClusterLists(const vector<size_t> &clusters, size_t clusterCount = 2)
{
    // constexpr size_t clusterCount = 2;
    vector<vector<size_t>> clusterLists(clusterCount);
    for (size_t i = 0; i < clusters.size(); i++)
    {
        clusterLists[clusters[i]].push_back(i);
    }
    return clusterLists;
}

//?
template <class RNG>
vector<sig_type> createRandomSigs(RNG &&rng, const vector<sig_type> &sigs, size_t clusterCount = 2)
{
    // constexpr size_t clusterCount = 2;
    vector<sig_type> clusterSigs(signatureSize * clusterCount);
    size_t signatureCount = sigs.size() / signatureSize;
    uniform_int_distribution<size_t> dist(0, signatureCount - 1);
    bool finished = false;

    unordered_set<string> uniqueSigs;
    for (size_t i = 0; i < signatureCount; i++)
    {
        size_t sig = dist(rng);
        string sigData(signatureSize * sizeof(sig_type), ' ');
        memcpy(&sigData[0], &sigs[sig * signatureSize], signatureSize * sizeof(sig_type));
        uniqueSigs.insert(sigData);
        if (uniqueSigs.size() >= clusterCount)
        {
            finished = true;
            break;
        }
    }

    size_t i = 0;
    for (const auto &sig : uniqueSigs)
    {
        memcpy(&clusterSigs[i * signatureSize], sig.data(), signatureSize * sizeof(sig_type));
        i++;
    }

    if (!finished)
    {
        if (uniqueSigs.size() != 1)
        {
            fprintf(stderr, "This should not happen\n");
            exit(1);
        }
        for (size_t i = 0; i < signatureSize; i++)
        {
            clusterSigs.push_back(clusterSigs[i]);
        }
    }
    return clusterSigs;
}

void reclusterSignatures(vector<size_t> &clusters, const vector<sig_type> &meanSigs, const vector<sig_type> &sigs, size_t clusterCount = 2)
{
    set<size_t> allClusters;
    for (size_t sig = 0; sig < clusters.size(); sig++)
    {
        const sig_type *sourceSignature = &sigs[sig * signatureSize];
        size_t minHdCluster = 0;
        size_t minHd = numeric_limits<size_t>::max();

        for (size_t cluster = 0; cluster < clusterCount; cluster++)
        {
            const sig_type *clusterSignature = &meanSigs[cluster * signatureSize];
            size_t hd = 0;
            for (size_t i = 0; i < signatureSize; i++)
            {
                hd += __builtin_popcountll(sourceSignature[i] ^ clusterSignature[i]);
            }
            if (hd < minHd)
            {
                minHd = hd;
                minHdCluster = cluster;
            }
        }
        clusters[sig] = minHdCluster;
        allClusters.insert(minHdCluster);
    }

    if (allClusters.size() == 1)
    {
        // We can't have everything in the same cluster.
        // If this did happen, just split them evenly
        for (size_t sig = 0; sig < clusters.size(); sig++)
        {
            clusters[sig] = sig % clusterCount;
        }
    }
}

vector<sig_type> createClusterSigs(const vector<vector<size_t>> &clusterLists, const vector<sig_type> &sigs, size_t count, size_t clusterCount = 2)
{
    // constexpr size_t clusterCount = 2;
    vector<sig_type> clusterSigs(signatureSize * clusterCount);
    //#pragma omp parallel
    {
        //#pragma omp for
        for (size_t cluster = 0; cluster < clusterLists.size(); cluster++)
        {
            sig_type *flattenedSignature = &clusterSigs[cluster * signatureSize];
            fill(flattenedSignature, flattenedSignature + signatureSize, 0);

            vector<size_t> counter(signatureSize * bits_per_char);

            for (size_t signature : clusterLists[cluster])
            {
                const sig_type *signatureData = &sigs[signatureSize * signature];
                for (size_t i = 0; i < signatureSize; i++)
                {
                    for (int n = 0; n < bits_per_char; n++)
                    {
                        if ((signatureData[i] >> n) & 1)
                        {
                            counter[i * bits_per_char + n]++;
                        }
                    }
                    flattenedSignature[i] |= signatureData[i];
                }
            }
        }
    }
    return clusterSigs;
}

template <class RNG>
void kMeans(RNG &&rng, const vector<sig_type> &sigs, size_t clusterCount)
{
    size_t nodeSigs = sigs.size() / signatureSize;
    vector<sig_type> meanSigs = createRandomSigs(rng, sigs, clusterCount);
    vector<size_t> clusters(nodeSigs);
    vector<vector<size_t>> clusterLists;

    for (int iteration = 0; iteration < 4; iteration++)
    {
        reclusterSignatures(clusters, meanSigs, sigs, clusterCount);
        clusterLists = createClusterLists(clusters, clusterCount);
        meanSigs = createClusterSigs(clusterLists, sigs, nodeSigs, clusterCount);
    }

    for (int i = 0; i < clusterCount; i++)
    {
        fprintf(stderr, ">");
        for (size_t seq : clusterLists[i])
        {
            fprintf(stderr, ",%zu", seq);
        }
        fprintf(stderr, "\n");
    }
}

class self_tree
{
public:
    size_t root = 0;                   // # of root node
    vector<size_t> childCounts;        // n entries, number of children
    vector<int> isBranchNode;          // n entries, is this a branch node
    vector<int> isLeafNode;            // n entries, is this a leaf node
    vector<vector<size_t>> childLinks; // n * o entries, links to children
    vector<vector<size_t>> dest;       // n * o entries, links to seq idx
    vector<size_t> parentLinks;        // n entries, links to parents
    vector<sig_type> means;            // n * signatureSize entries, node signatures
    vector<vector<sig_type>> matri;    // capacity * signatureSize * n
    vector<sig_type> matrices;         // capacity * signatureSize * n
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
                isLeafNode.resize(capacity);
            }
            //#pragma omp single
            {
                childLinks.resize(capacity);
            }
            //#pragma omp single
            {
                dest.resize(capacity);
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
                matri.resize(capacity);
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

    size_t calcHD(const sig_type *a, const sig_type *b) const
    {
        size_t c = 0;
        for (size_t i = 0; i < signatureSize; i++)
        {
            c += __builtin_popcountll(a[i] ^ b[i]);
        }
        return c;
    }

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

    //?
    inline size_t traverse(size_t parent, sig_type signature, vector<size_t> &insertionList)
    {
        for (size_t i = 0; i < childCounts[parent]; i++)
        {

            size_t child = childLinks[parent][i];
            if (matrices[child] == signature)
            {
                // std::cout << child << "found\n";
                // childCounts[child]++;
                return child;
            }
        }

        size_t node = getNewNodeIdx(insertionList);
        addSigToMatrix(node, signature);
        childLinks[parent].push_back(node);
        parentLinks[node] = parent;
        childCounts[parent]++;
        isBranchNode[parent] = 1;

        return node;
    }

    // void addSigToMatrix(size_t node, const seq_type *seq)
    // {
    //     matri[node].insert(matri[node].end(), seq, seq + signatureSize);
    // }

    void addSigToMatrix(size_t node, const sig_type *sig)
    {
        matri[node].insert(matri[node].end(), sig, sig + signatureSize);
    }

    void addSigToMatrix(size_t node, sig_type sig)
    {
        matrices[node] = sig;
    }

    void printMatrix(size_t node)
    {
        for (sig_type val : matri[node])
        {
            std::cout << val << "n";
        }
    }

    //?
    void recalculateSig(size_t node)
    {
        size_t children = childCounts[node];
        vector<sig_type> matrix = matri[node];
        sig_type *sig = &means[node * signatureSize];
        fill(sig, sig + signatureSize, static_cast<unsigned char>(0x00));

        for (size_t i = 0; i < children - 1; i++)
        {
            // dbgPrintSignature(&matrix[i * signatureSize]);
            size_t k = 0;
            for (size_t j = 0; j < signatureSize; j++)
            {
                for (int n = bits_per_char - 1; n >= 0; n--)
                {
                    k++;
                }

                sig[j] |= matrix[i * signatureSize + j];
            }
            // std::cerr << "\n----------------------\n";
        }

        size_t i = children - 1;
        size_t k = 0;
        size_t counts_t = 0;
        for (size_t j = 0; j < signatureSize; j++)
        {
            for (int n = bits_per_char - 1; n >= 0; n--)
            {
                k++;
            }

            sig[j] |= matrix[i * signatureSize + j];
        }
    }

    void recalculateUp(size_t node)
    {
        size_t limit = 10;
        // fprintf(stderr, "RecalculateUp %zu\n", node);
        while (node != root)
        {
            recalculateSig(node);
            node = parentLinks[node];
            if (omp_test_lock(&locks[node]))
            {
                omp_unset_lock(&locks[node]);
            }
            else
            {
                break;
            }

            // Put a limit on how far we go up
            // At some point it stops mattering, plus this helps avoid inf loops
            // caused by cycles getting into the tree structure
            limit--;
            if (limit == 0)
                return;
            // fprintf(stderr, "-> %zu\n", node);
        }
        // update root
        recalculateSig(node);
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

    void printNode(size_t node)
    {
        fprintf(stdout, "%zu", node);
        for (size_t i = 0; i < childCounts[node]; i++)
        {
            size_t child = childLinks[node][i];
            fprintf(stdout, "\t%zu", child);
        }
        fprintf(stdout, "\n");
    }

    void printLeave(size_t node)
    {
        fprintf(stdout, "-%zu: \n", node);
        for (size_t leaf : childLinks[node])
        {
            fprintf(stdout, "%zu,", leaf);
        }
        fprintf(stdout, "\n");
    }

    // breadth first, printing top down
    void printTree()
    {
        if (isBranchNode[root])
        {

            set<size_t> parents = {root};
            vector<size_t> leaves;
            while (parents.size() > 0)
            {
                vector<size_t> temp_parents;
                fprintf(stdout, ">\n");
                for (size_t parent : parents)
                {
                    printNode(parent);
                    for (size_t child : childLinks[parent])
                    {
                        if (isBranchNode[child])
                        {
                            temp_parents.push_back(child);
                        }
                        else
                        {
                            leaves.push_back(child);
                        }
                    }
                }
                parents = set<size_t>(temp_parents.begin(), temp_parents.end());
                if (leaves.size() > 0)
                {
                    fprintf(stdout, ">\n");
                    for (size_t leaf : leaves)
                    {
                        printLeave(leaf);
                    }
                    fprintf(stdout, "---\n");
                }
            }
        }
    }

    void printSubTree(size_t tnode)
    {
        if (!isLeafNode[tnode])
        {
            // node idx;minimiser;{child idx,...}
            std::cout << tnode << ";" << matrices[tnode] << ";";
            for (size_t child : childLinks[tnode])
            {
                std::cout << child << ",";
            }
            std::cout << "\n";

            for (size_t child : childLinks[tnode])
            {
                printSubTree(child);
            }
        }
        else
        {
            std::cout << tnode << ";" << matrices[tnode] << "\n>>>";
            for (size_t seq : dest[tnode])
            {
                std::cout << seq << ",";
            }
            std::cout << "\n";

            // a leaf node can also be a branch if there is a shorter seqs that share the same set of minimisers
            if (isBranchNode[tnode])
            {
                // node idx;minimiser;{child idx,...}
                std::cout << tnode << ";" << matrices[tnode] << ";";
                for (size_t child : childLinks[tnode])
                {
                    std::cout << child << ",";
                }
                std::cout << "\n";
                for (size_t child : childLinks[tnode])
                {
                    printSubTree(child);
                }
            }
        }
    }

    void printNodeJson(size_t tnode)
    {
        cout << "{\"node\":\"" << tnode << "\",\"content\":\"*" << matrices[tnode] << "\",\"children\":[";
        // for (size_t i = 0; i<childLinks[tnode].size() -1;i++)
        // {
        //     std::cout << childLinks[tnode][i] << ",";
        // }
        // std::cout << childLinks[tnode].back() << "]}\n";
        // printSubTreeJson
    }

    void printSubTreeJson(size_t tnode)
    {
        if (!isLeafNode[tnode])
        {
            printNodeJson(tnode);
            printSubTreeJson(childLinks[tnode][0]);

            for (size_t i = 1; i<childLinks[tnode].size();i++)
            {
                std::cout << ",";
                printSubTreeJson(childLinks[tnode][i]);
            }

            // std::cout << "]}";
        }
        else
        {
            cout << "{\"node\":\"" << tnode << "\",\"content\":\"*" << matrices[tnode] << "\",\"children\":[{\"node\":\">>>";
            for (size_t seq : dest[tnode])
            {
                std::cout << seq << ",";
            }
            std::cout << "\"}";

            // a leaf node can also be a branch if there is a shorter seqs that share the same set of minimisers
            if (isBranchNode[tnode])
            {
                std::cout << ",";
                printSubTreeJson(childLinks[tnode][0]);

                // for (size_t i = 0; i<childLinks[tnode].size() -1;i++)
                // {
                //     std::cout << ",";
                //     printSubTreeJson(childLinks[tnode][i]);
                // }

                // std::cout << "]}";

                // cout << "{\"content\":\"---\"}";
            }
            // else
            // {
            //     std::cout << "]}";
            // }
        }
        std::cout << "]}";
    }

    // inline void first_insert(vector<sig_type> signature, vector<size_t> &insertionList, size_t idx)
    // {
    //     sig_type temp = 0;
    //     size_t size = sizeof(sig_type);

    //     addSigToMatrix(root, (size_t)-1);
    //     size_t parent = root;
    //     size_t node = 0;
    //     for (int i = 0; i < signature.size(); i++)
    //     {
    //         // memcpy(&temp, &signature[i * size], size);
    //         temp = signature[i];
    //         node = getNewNodeIdx(insertionList);
    //         addSigToMatrix(node, signature[i]);
    //         // std::cout << i << "," << temp << "\n";
    //         parentLinks[node] = parent;
    //         childLinks[parent].push_back(node);
    //         isBranchNode[parent] = 1;
    //         childCounts[parent]++;
    //         parent = node;
    //     }

    //     // last node record the idx
    //     // isBranchNode[parent] = 1;
    //     isLeafNode[node] = 1;
    //     dest[node].push_back(idx);
    //     // printMatrix(root);
    //     // std::cout << "print tree\n";
    //     // size_t tnode = root;
    //     // while (isBranchNode[tnode])
    //     // {
    //     //     std::cout << tnode << "," << matrices[tnode] << "\n";
    //     //     tnode = childLinks[tnode][0];
    //     // }
    //     // std::cout << tnode << "," << matrices[tnode] << "\n";
    // }

    
    inline void first_insert(seq_type signature, vector<size_t> &insertionList, size_t idx)
    {
        size_t parent = root;
        size_t node = getNewNodeIdx(insertionList);
        parentLinks[node] = parent;
        childLinks[parent].push_back(node);
        isBranchNode[parent] = 1;
        childCounts[parent]++;
        addSigToMatrix(node, &signature[0]);
    }
    
    
    inline void next_insert(vector<sig_type> signature, vector<size_t> &insertionList, size_t idx)
    {
        size_t parent = traverse(root, signature[0], insertionList);

        // omp_set_lock(&locks[insertionPoint]);

        // size_t parent = root;

        size_t node = 0;

        for (int i = 1; i < signature.size(); i++)
        {
            // size_t node = getNewNodeIdx(insertionList);
            // addSigToMatrix(node, signature[i]);
            // parentLinks[node] = parent;
            // childLinks[parent].push_back(node);
            // isBranchNode[parent] = 1;
            // childCounts[parent]++;
            node = traverse(parent, signature[i], insertionList);
            parent = node;
        }
        // isBranchNode[parent] = 1;

        isLeafNode[node] = 1;
        dest[node].push_back(idx);
        // std::cout << ">>print tree\n";
        // printSubTree(root);
        // omp_unset_lock(&locks[insertionPoint]);
    }

    size_t search(size_t parent, vector<sig_type> signature, size_t counter, vector<size_t> &matches)
    {

        // reached the end of seq, no match is found
        if (counter == signature.size())
            return 0;

        for (size_t i = 0; i < childCounts[parent]; i++)
        {

            size_t child = childLinks[parent][i];
            if (matrices[child] == signature[counter])
            {
                // // this can be partial match
                // if (isLeafNode[child])
                // {
                //     matches.push_back(dest[child][0]);
                // }

                if (counter == signature.size() - 1 && isLeafNode[child])
                {
                    matches.insert(matches.end(), dest[child].begin(), dest[child].end());
                }
                else
                {
                    search(child, signature, counter + 1, matches);
                }
            }
        }
    }

    void search(vector<sig_type> signature)
    {
        vector<size_t> matches;
        search(root, signature, 0, matches);
        for (size_t match : matches)
        {
            std::cout << match << ",";
        }
        std::cout << "\n";
    }

    template <class RNG>
    inline void insert(RNG &&rng, const sig_type *signature, vector<size_t> &insertionList, size_t idx)
    {
        size_t insertionPoint = traverse(signature);

        // fprintf(stdout, " at %zu\n", insertionPoint);
        omp_set_lock(&locks[insertionPoint]);

        // if (childCounts[insertionPoint] < order)
        {
            addSigToMatrix(insertionPoint, signature);

            childLinks[insertionPoint].push_back(idx);
            fprintf(stdout, " at %zu\n", insertionPoint);

            childCounts[insertionPoint]++;
            recalculateSig(insertionPoint);
            if (isBranchNode[root])
            {
                recalculateUp(parentLinks[insertionPoint]);
            }
        }
        // else
        // {
        //     fprintf(stdout, " \n");
        //     splitNode(rng, insertionPoint, signature, insertionList, idx);
        //     // splitNode(rng, insertionPoint, signature, insertionList, 0);
        // }
        omp_unset_lock(&locks[insertionPoint]);
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