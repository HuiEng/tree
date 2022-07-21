// There are two kinds of ktree nodes- branch nodes and leaf nodes
// Both contain a signature matrix, plus their own signature
// (the root node signature does not matter and can be blank)
// Branch nodes then contain 'order' links to other nodes
// Leaf nodes do not.
// However, as leaf nodes may become branch nodes, we allocate
// the space anyway.
// As the space to be used is determined at runtime, we use
// parallel arrays, not structs

#ifndef INCLUDE_BF_KTree_HPP
#define INCLUDE_BF_KTree_HPP

#include <omp.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include "bloom_filter.hpp"

using namespace std;

typedef unsigned char cell_type;
bloom_parameters parameters;
size_t signatureSize; // Signature size (depends on element in BF, obtained while read binary)
size_t ktree_order = 10;
size_t ktree_capacity = 1000000;
bool intersect = false;
double threshold = 0;
size_t countThreshold = 2;

// void hashToCounter(cell_type hash, vector<size_t> &counter)
// {
//     int binary[bits_per_char];
//     for (int n = 0; n < bits_per_char; n++)
//     {
//         if ((hash >> n) & 1)
//             counter[bits_per_char - 1 - n]++;
//     }
// }

void dbgPrintSignature(const cell_type *sig)
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
vector<cell_type> createRandomSigs(RNG &&rng, const vector<cell_type> &sigs, size_t clusterCount = 2)
{
    // constexpr size_t clusterCount = 2;
    vector<cell_type> clusterSigs(signatureSize * clusterCount * 2); // front 2 for union, back 2 for inter
    size_t signatureCount = sigs.size() / signatureSize;
    uniform_int_distribution<size_t> dist(0, signatureCount - 1);
    bool finished = false;

    // memcpy(&clusterSigs[0], &sigs[0], signatureSize * sizeof(cell_type));
    // memcpy(&clusterSigs[signatureSize], &sigs[(signatureCount - 1) * signatureSize], signatureSize * sizeof(cell_type));

    unordered_set<string> uniqueSigs;
    for (size_t i = 0; i < signatureCount; i++)
    {
        size_t sig = dist(rng);
        string sigData(signatureSize * sizeof(cell_type), ' ');
        memcpy(&sigData[0], &sigs[sig * signatureSize], signatureSize * sizeof(cell_type));
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
        memcpy(&clusterSigs[i * signatureSize], sig.data(), signatureSize * sizeof(cell_type));
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

    // append a copy to the end for inter
    memcpy(&clusterSigs[(clusterCount - 1) * 2 * signatureSize], &clusterSigs[0], 2 * signatureSize * sizeof(cell_type));

    return clusterSigs;
}

void reclusterSignatures(vector<size_t> &clusters, const vector<cell_type> &meanSigs, const vector<cell_type> &sigs, size_t clusterCount = 2)
{
    set<size_t> allClusters;
    for (size_t sig = 0; sig < clusters.size(); sig++)
    {
        const cell_type *sourceSignature = &sigs[sig * signatureSize];
        size_t minHdCluster = 0;
        size_t minHd = numeric_limits<size_t>::max();

        for (size_t cluster = 0; cluster < clusterCount; cluster++)
        {
            const cell_type *clusterSignature = &meanSigs[cluster * signatureSize];
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

void reclusterSignatures(vector<size_t> &clusters, const vector<cell_type> &meanSigs, const vector<cell_type> &sigs_union, const vector<cell_type> &sigs_inter)
{
    size_t clusterCount = 2;
    set<size_t> allClusters;
    if (intersect)
    {
        for (size_t sig = 0; sig < clusters.size(); sig++)
        {
            const cell_type *sourceSignature_union = &sigs_union[sig * signatureSize];
            const cell_type *sourceSignature_inter = &sigs_inter[sig * signatureSize];
            size_t maxBitsCluster = 0;
            size_t maxBitsCluster_union = 0;
            size_t max_cmn_set_bits_union = 0;
            size_t maxBitsCluster_inter = 0;
            size_t max_cmn_set_bits_inter = 0; // threshold;
            size_t max_cmn_set_bits_inter_union = 0;

            for (size_t cluster = 0; cluster < clusterCount; cluster++)
            {
                size_t cmn_union = BFintersect(sourceSignature_union, &meanSigs[cluster * signatureSize], signatureSize);
                size_t cmn_inter = BFintersect(sourceSignature_inter, &meanSigs[(cluster + clusterCount) * signatureSize], signatureSize);

                if (cmn_union >= max_cmn_set_bits_union)
                {
                    if (cmn_inter > max_cmn_set_bits_inter)
                    {
                        max_cmn_set_bits_union = cmn_union;
                        max_cmn_set_bits_inter = cmn_inter;
                        maxBitsCluster_union = cluster;
                    }
                }
            }
            maxBitsCluster = maxBitsCluster_union;

            // if (max_cmn_set_bits_inter > signatureSize * threshold)
            // {
            //     maxBitsCluster = maxBitsCluster_inter;
            // }

            clusters[sig] = maxBitsCluster;
            allClusters.insert(maxBitsCluster);
        }
    }
    else
    {
        for (size_t sig = 0; sig < clusters.size(); sig++)
        {
            const cell_type *sourceSignature = &sigs_union[sig * signatureSize];
            size_t minHdCluster = 0;
            size_t minHd = numeric_limits<size_t>::max();

            for (size_t cluster = 0; cluster < 2; cluster++)
            {
                const cell_type *clusterSignature = &meanSigs[cluster * signatureSize];
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
    }

    if (allClusters.size() == 1)
    {
        // We can't have everything in the same cluster.
        // If this did happen, just split them evenly
        for (size_t sig = 0; sig < clusters.size(); sig++)
        {
            clusters[sig] = sig % 2;
        }
    }
}

void reclusterSignatures2(vector<size_t> &clusters, const vector<cell_type> &meanSigs, const vector<cell_type> &sigs_union, const vector<cell_type> &sigs_inter)
{
    size_t clusterCount = 2;
    set<size_t> allClusters;
    if (intersect)
    {
        for (size_t sig = 0; sig < clusters.size(); sig++)
        {
            const cell_type *sourceSignature_union = &sigs_union[sig * signatureSize];
            const cell_type *sourceSignature_inter = &sigs_inter[sig * signatureSize];
            size_t maxBitsCluster = 0;
            size_t maxBitsCluster_union = 0;
            size_t max_cmn_set_bits_union = 0;
            size_t maxBitsCluster_inter = 0;
            size_t max_cmn_set_bits_inter = 0; // threshold;
            size_t max_cmn_set_bits_inter_union = 0;

            for (size_t cluster = 0; cluster < clusterCount; cluster++)
            {
                size_t cmn_union = BFintersect(sourceSignature_union, &meanSigs[cluster * signatureSize], signatureSize);
                size_t cmn_inter = BFintersect(sourceSignature_inter, &meanSigs[(cluster + clusterCount) * signatureSize], signatureSize);

                if (cmn_union >= max_cmn_set_bits_union)
                {
                    if (cmn_inter > max_cmn_set_bits_inter)
                    {
                        max_cmn_set_bits_union = cmn_union;
                        max_cmn_set_bits_inter = cmn_inter;
                        maxBitsCluster_union = cluster;
                    }
                }
            }
            // fprintf(stdout, "cmn_union:%zu\n", max_cmn_set_bits_union);
            maxBitsCluster = maxBitsCluster_union;

            // if (max_cmn_set_bits_inter > signatureSize * threshold)
            // {
            //     maxBitsCluster = maxBitsCluster_inter;
            // }

            clusters[sig] = maxBitsCluster;
            allClusters.insert(maxBitsCluster);
        }
    }
    else
    {
        for (size_t sig = 0; sig < clusters.size(); sig++)
        {
            const cell_type *sourceSignature = &sigs_union[sig * signatureSize];
            size_t minHdCluster = 0;
            size_t minHd = numeric_limits<size_t>::max();

            for (size_t cluster = 0; cluster < 2; cluster++)
            {
                const cell_type *clusterSignature = &meanSigs[cluster * signatureSize];
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
    }

    if (allClusters.size() == 1)
    {
        // We can't have everything in the same cluster.
        // If this did happen, just split them evenly
        for (size_t sig = 0; sig < clusters.size(); sig++)
        {
            clusters[sig] = sig % 2;
        }
    }
}

vector<cell_type> createClusterSigs(const vector<vector<size_t>> &clusterLists, const vector<cell_type> &sigs, size_t count, size_t clusterCount = 2)
{
    // constexpr size_t clusterCount = 2;
    vector<cell_type> clusterSigs(signatureSize * clusterCount * 2); // front 2 for union, back 2 for inter
    //#pragma omp parallel
    {
        //#pragma omp for
        for (size_t cluster = 0; cluster < clusterLists.size(); cluster++)
        {
            cell_type *flattenedSignature_union = &clusterSigs[cluster * signatureSize];
            fill(flattenedSignature_union, flattenedSignature_union + signatureSize, 0);

            cell_type *flattenedSignature_inter = &clusterSigs[(cluster + clusterCount) * signatureSize];
            fill(flattenedSignature_inter, flattenedSignature_inter + signatureSize, 1);

            vector<size_t> counter(signatureSize * bits_per_char);

            for (size_t signature : clusterLists[cluster])
            {
                const cell_type *signatureData = &sigs[signatureSize * signature];
                for (size_t i = 0; i < signatureSize; i++)
                {
                    for (int n = 0; n < bits_per_char; n++)
                    {
                        if ((signatureData[i] >> n) & 1)
                        {
                            counter[i * bits_per_char + n]++;
                        }
                    }
                    flattenedSignature_union[i] |= signatureData[i];
                    flattenedSignature_inter[i] &= signatureData[i];
                }
            }

            // for (int i = 0; i < counter.size(); i++)
            // {
            //     if (counter[i] > ktree_order * threshold)
            //     {
            //         flattenedSignature_union[i / bits_per_char] |= (cell_type)1 << (i % bits_per_char);
            //     }
            // }
        }
    }
    return clusterSigs;
}

template <class RNG>
void kMeans(RNG &&rng, const vector<cell_type> &sigs, size_t clusterCount)
{
    size_t nodeSigs = sigs.size() / signatureSize;
    vector<cell_type> meanSigs = createRandomSigs(rng, sigs, clusterCount);
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

class BF_KTree
{
public:
    size_t root = numeric_limits<size_t>::max(); // # of root node
    vector<size_t> childCounts;                  // n entries, number of children
    vector<int> isBranchNode;                    // n entries, is this a branch node
    vector<vector<size_t>> childLinks;           // n * o entries, links to children
    vector<vector<size_t>> counts;               // n * o entries, links to children
    vector<size_t> counts_threshold;             // n * o entries, links to children
    vector<size_t> parentLinks;                  // n entries, links to parents
    vector<cell_type> means_union;               // n * signatureSize entries, node signatures
    vector<cell_type> means_inter;               // n * signatureSize entries, node signatures
    vector<vector<cell_type>> matrices_union;    // capacity * signatureSize * n
    vector<vector<cell_type>> matrices_inter;    // capacity * signatureSize * n
    vector<omp_lock_t> locks;                    // n locks
    size_t order;
    size_t capacity = 0; // Set during construction, currently can't change

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
                counts.resize(capacity);
            }
            //#pragma omp single
            {
                counts_threshold.resize(capacity);
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
                //?
                matrices_union.resize(capacity);
            }
            //#pragma omp single
            {
                matrices_inter.resize(capacity);
            }
            //#pragma omp single
            {
                means_union.resize(capacity * signatureSize);
            }
            //#pragma omp single
            {
                means_inter.resize(capacity * signatureSize);
            }
        }
    }

    BF_KTree(size_t order_, size_t capacity) : order{order_}
    {
        reserve(capacity);
    }

    size_t calcHD(const cell_type *a, const cell_type *b) const
    {
        size_t c = 0;
        for (size_t i = 0; i < signatureSize; i++)
        {
            c += __builtin_popcountll(a[i] ^ b[i]);
        }
        return c;
    }

    // void printNode(size_t node)
    // {
    //     fprintf(stdout, "%zu", node);
    //     for (size_t i = 0; i < childCounts[node]; i++)
    //     {
    //         size_t child = childLinks[node][i];
    //         fprintf(stdout, "\t%zu", child);
    //     }
    //     fprintf(stdout, "\n");
    // }

    // void printLeave(size_t node)
    // {
    //     vector<cell_type> matrix = matrices_inter[node];
    //     size_t size = matrix.size() / signatureSize;
    //     fprintf(stdout, "%zu\t-%zu\n", node, size);

    //     cell_type *a = &means_inter[node];
    //     size_t c = 0;
    //     for (size_t i = 0; i < signatureSize; i++)
    //     {
    //         c += __builtin_popcountll(a[i]);
    //     }
    //     fprintf(stderr, "%zu\t: %zu children with %zu bits", node, size, c);
    //     for (size_t i = 0; i < size; i++)
    //     {
    //         cell_type *a = &matrix[i * signatureSize];
    //         c = 0;
    //         for (size_t i = 0; i < signatureSize; i++)
    //         {
    //             c += __builtin_popcountll(a[i]);
    //         }
    //         fprintf(stderr, ", %zu", c);
    //     }
    //     fprintf(stderr, "\n");
    // }

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

    size_t traverse_all(vector<size_t> &leaves, const cell_type *signature)
    {
        size_t distance_union = 0;
        size_t distance_inter = 0;
        size_t count_threshold = 0;
        size_t final = 0;

        for (size_t leaf : leaves)
        {
            size_t cmn_union = BFintersect(&means_union[leaf * signatureSize], signature, signatureSize);
            size_t cmn_inter = BFintersect(&means_inter[leaf * signatureSize], signature, signatureSize);
            size_t count = counts_threshold[leaf];

            if (cmn_union >= cmn_union)
            {
                // if (cmn_inter >= distance_inter)
                if (count >= count_threshold)
                {
                    distance_union = cmn_union;
                    distance_inter = cmn_inter;
                    final = leaf;
                    count_threshold = count;
                    fprintf(stderr, "%zu,%zu\n", distance_union, count_threshold);
                }
            }
        }
        // fprintf(stderr,"\n\n");
        return final;
    }
    // only prune when no common set bit
    // inline size_t traverse_all(const cell_type *signature, vector<size_t> nodes, size_t final, size_t distance) const
    // {
    //     // size_t node = root;
    //     // size_t max = countSetBits(signature, signatureSize);
    //     // size_t min = BFintersect(&means_union[final * signatureSize], signature, signatureSize);
    //     vector<size_t> candidates;

    //     for (size_t node : nodes)
    //     {
    //         if (isBranchNode[node])
    //         {
    //             fprintf(stderr,"%zu\n",node);
    //             for (size_t i = 0; i < childCounts[node]; i++)
    //             {
    //                 size_t child = childLinks[node][i];
    //                 size_t cmn_union = BFintersect(&means_union[child * signatureSize], signature, signatureSize);
    //                 // fprintf(stderr,"***%zu,%zu\n",child,cmn_union);
    //                 if (cmn_union > threshold)
    //                 {
    //                     candidates.push_back(child);
    //                 }
    //             }
    //         }else{
    //             size_t temp_dist = BFintersect(&means_union[node * signatureSize], signature, signatureSize);
    //             if (temp_dist > distance){
    //                 final = node;
    //                 distance = temp_dist;
    //                 fprintf(stderr,"%zu,%zu\n",node,distance);
    //             }else{
    //                 // fprintf(stderr,"---%zu,%zu\n",node,distance);
    //             }
    //         }
    //     }
    //     if (candidates.size() > 0)
    //     {
    //         traverse_all(signature, candidates, final, distance);
    //     }

    //     fprintf(stderr,">*-final %zu,%zu\n",final,distance);
    //     return final;

    // }

    inline size_t traverse(const cell_type *signature) const
    {
        size_t node = root;
        if (intersect)
        {
            while (isBranchNode[node])
            {
                size_t maxBitschild = childLinks[node][0];
                size_t max_cmn_set_bits_union = 0;
                size_t max_cmn_set_bits_inter = 0; // threshold;
                size_t max_cmn_set_bits_inter_union = 0;
                size_t maxBitsChild_union = maxBitschild;
                size_t maxBitsChild_inter = maxBitschild;
                size_t count_t = 0;

                for (size_t i = 0; i < childCounts[node]; i++)
                {
                    size_t child = childLinks[node][i];
                    size_t cmn_union = BFintersect(&means_union[child * signatureSize], signature, signatureSize);
                    size_t cmn_inter = BFintersect(&means_inter[child * signatureSize], signature, signatureSize);
                    size_t count = counts_threshold[child];

                    if (cmn_union >= max_cmn_set_bits_union)
                    {
                        // if (count >= count_t)
                        if (cmn_inter >= max_cmn_set_bits_inter)
                        {
                            max_cmn_set_bits_union = cmn_union;
                            max_cmn_set_bits_inter = cmn_inter;
                            // count_t = count;
                            maxBitsChild_union = child;
                            // fprintf(stdout, ",%zu,", count_t);
                        }
                    }
                }
                maxBitschild = maxBitsChild_union;
                // if (max_cmn_set_bits_inter > signatureSize * threshold)
                // {
                //     maxBitschild = maxBitsChild_inter;
                // }
                node = maxBitschild;
            }
        }
        else
        {
            while (isBranchNode[node])
            {
                size_t lowestHD = numeric_limits<size_t>::max();
                size_t lowestHDchild = childLinks[node][0];

                for (size_t i = 0; i < childCounts[node]; i++)
                {
                    size_t child = childLinks[node][i];
                    size_t hd = calcHD(&means_union[child * signatureSize], signature);
                    if (hd < lowestHD)
                    {
                        lowestHD = hd;
                        lowestHDchild = child;
                    }
                }

                node = lowestHDchild;
            }
        }

        return node;
    }

    void addSigToMatrix_union(size_t node, const cell_type *sig)
    {
        matrices_union[node].insert(matrices_union[node].end(), sig, sig + signatureSize);
    }

    void addSigToMatrix_inter(size_t node, const cell_type *sig)
    {
        matrices_inter[node].insert(matrices_inter[node].end(), sig, sig + signatureSize);
    }

    //?
    void recalculateSig(size_t node)
    {
        size_t children = childCounts[node];
        vector<cell_type> matrix_union = matrices_union[node];
        cell_type *sig_union = &means_union[node * signatureSize];
        fill(sig_union, sig_union + signatureSize, static_cast<unsigned char>(0x00));

        vector<cell_type> matrix_inter = matrices_inter[node];
        cell_type *sig_inter = &means_inter[node * signatureSize];
        fill(sig_inter, sig_inter + signatureSize, 1);

        // vector<size_t> counter(signatureSize * bits_per_char);
        counts[node].clear();
        counts[node].resize(signatureSize * bits_per_char);

        for (size_t i = 0; i < children - 1; i++)
        {
            // dbgPrintSignature(&matrix_union[i * signatureSize]);
            size_t k = 0;
            for (size_t j = 0; j < signatureSize; j++)
            {
                for (int n = bits_per_char - 1; n >= 0; n--)
                {
                    if ((matrix_union[i * signatureSize + j] >> n) & 1)
                    {
                        counts[node][k]++;
                    }
                    k++;
                }

                sig_union[j] |= matrix_union[i * signatureSize + j];
                sig_inter[j] &= matrix_inter[i * signatureSize + j];
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
                if ((matrix_union[i * signatureSize + j] >> n) & 1)
                {
                    counts[node][k]++;
                }

                if (counts[node][k] > countThreshold)
                {
                    counts_t++;
                }
                k++;
            }

            sig_union[j] |= matrix_union[i * signatureSize + j];
            sig_inter[j] &= matrix_inter[i * signatureSize + j];
        }
        counts_threshold[node] = counts_t;

        // for (int i = 0; i < counter.size(); i++)
        // {
        //     if (counter[i] > order * threshold)
        //     {
        //         sig_union[i / bits_per_char] |= (cell_type)1 << (i % bits_per_char);
        //     }
        // }

        // std::cerr << "\n";
        // dbgPrintSignature(sig_union);
        // for (size_t count : counts[node])
        // {
        //     std::cerr << count << ",";
        // }
        // std::cerr << "\n----------------------\n";
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
    size_t findParent(size_t node, cell_type *meanSig_union, cell_type *meanSig_inter)
    {
        size_t parent = parentLinks[node];
        // Lock the parent
        omp_set_lock(&locks[parent]);

        size_t idx = numeric_limits<size_t>::max();
        for (size_t i = 0; i < childCounts[parent]; i++)
        {
            // if (childLinks[parent * order + i] == node) {
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
                // if (childLinks[parent * order + i] == node) {
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

        memcpy(&matrices_union[parent][idx * signatureSize], meanSig_union, sizeof(cell_type) * signatureSize);
        memcpy(&matrices_inter[parent][idx * signatureSize], meanSig_inter, sizeof(cell_type) * signatureSize);

        return parent;
    }

    void update(size_t node)
    {
        for (int i = 0; i < childCounts[node]; i++)
        {
            size_t child = childLinks[node][i];
            memcpy(&matrices_union[node][i * signatureSize], &means_union[child * signatureSize], sizeof(cell_type) * signatureSize);
            memcpy(&matrices_inter[node][i * signatureSize], &means_inter[child * signatureSize], sizeof(cell_type) * signatureSize);
        }
        if (isBranchNode[parentLinks[node]])
        {
            // fprintf(stdout,"again %zu %zu\n",node, parentLinks[node]);
            update(parentLinks[node]);
        }
    }

    template <class RNG>
    void splitNode(RNG &&rng, size_t node, const cell_type *sig, vector<size_t> &insertionList, size_t link)
    {
        size_t nodeSigs = childCounts[node] + 1;
        vector<cell_type> sigs_union(nodeSigs * signatureSize);
        //?
        copy(matrices_union[node].begin(), matrices_union[node].end(), sigs_union.begin());
        memcpy(&sigs_union[childCounts[node] * signatureSize], sig, sizeof(cell_type) * signatureSize);

        vector<cell_type> sigs_inter(nodeSigs * signatureSize);
        //?
        copy(matrices_inter[node].begin(), matrices_inter[node].end(), sigs_inter.begin());
        memcpy(&sigs_inter[childCounts[node] * signatureSize], sig, sizeof(cell_type) * signatureSize);

        vector<cell_type> meanSigs = createRandomSigs(rng, sigs_union);
        vector<size_t> clusters(nodeSigs);
        vector<vector<size_t>> clusterLists;

        // dbgPrintSignature(&meanSigs[0]);
        // dbgPrintSignature(&meanSigs[signatureSize]);
        vector<cell_type> temp(2 * signatureSize);

        for (int iteration = 0; iteration < 4; iteration++)
        {
            // fprintf(stderr, "Iteration %d\n", iteration);
            reclusterSignatures(clusters, meanSigs, sigs_union, sigs_inter);
            clusterLists = createClusterLists(clusters);
            memcpy(&temp[0], &meanSigs[0], sizeof(meanSigs));
            meanSigs = createClusterSigs(clusterLists, sigs_union, nodeSigs);
        }

        reclusterSignatures2(clusters, temp, sigs_union, sigs_inter);

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
                addSigToMatrix_union(sibling, &sigs_union[seqIdx * signatureSize]);
                addSigToMatrix_inter(sibling, &sigs_inter[seqIdx * signatureSize]);
                siblingIdx++;
            }
        }

        memcpy(&means_union[node * signatureSize], &meanSigs[0], sizeof(cell_type) * signatureSize);
        memcpy(&means_union[sibling * signatureSize], &meanSigs[signatureSize], sizeof(cell_type) * signatureSize);

        memcpy(&means_inter[node * signatureSize], &meanSigs[2 * signatureSize], sizeof(cell_type) * signatureSize);
        memcpy(&means_inter[sibling * signatureSize], &meanSigs[3 * signatureSize], sizeof(cell_type) * signatureSize);

        childCounts[node] = clusterLists[0].size();
        // Fill the current node with the other cluster of signatures
        {
            matrices_union[node].clear();
            matrices_inter[node].clear();
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
                addSigToMatrix_union(node, &sigs_union[seqIdx * signatureSize]);
                addSigToMatrix_inter(node, &sigs_inter[seqIdx * signatureSize]);
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
            addSigToMatrix_union(newRoot, &meanSigs[0 * signatureSize]);
            addSigToMatrix_union(newRoot, &meanSigs[1 * signatureSize]);

            addSigToMatrix_inter(newRoot, &meanSigs[2 * signatureSize]);
            addSigToMatrix_inter(newRoot, &meanSigs[3 * signatureSize]);

            root = newRoot;

            recalculateSig(root);
        }
        else
        {
            // First, update the reference to this node in the parent with the new mean
            size_t parent = findParent(node, &meanSigs[0], &meanSigs[2 * signatureSize]);

            // Connect sibling node to parent
            parentLinks[sibling] = parent;

            if (childCounts[parent] + 1 < order)
            {
                childLinks[parent].push_back(sibling);

                childCounts[parent]++;
                addSigToMatrix_union(parent, &meanSigs[1 * signatureSize]);
                addSigToMatrix_inter(parent, &meanSigs[3 * signatureSize]);

                // Update signatures (may change?)
                recalculateUp(parent);
            }
            else
            {
                // fprintf(stdout, "after\n");
                splitNode(rng, parent, &meanSigs[1 * signatureSize], insertionList, sibling);
            }
            update(parent);
            // Unlock the parent
            omp_unset_lock(&locks[parent]);
        }
        // fprintf(stdout, ">>>\n");
        // for (int i = 0; i < childCounts[root]; i++)
        // {
        //     dbgPrintSignature(&matrices_union[root][i * signatureSize]);
        //     dbgPrintSignature(&means_union[childLinks[root][i] * signatureSize]);
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
        // vector<cell_type> matrix = matrices_inter[node];
        // size_t size = matrix.size() / signatureSize;
        // fprintf(stdout, "%zu\t-%zu\n", node, size);

        // cell_type *a = &means_inter[node];
        // size_t c = 0;
        // for (size_t i = 0; i < signatureSize; i++)
        // {
        //     c += __builtin_popcountll(a[i]);
        // }
        // fprintf(stderr, "%zu\t: %zu children with %zu bits", node, size, c);
        // for (size_t i = 0; i < size; i++)
        // {
        //     cell_type *a = &matrix[i * signatureSize];
        //     c = 0;
        //     for (size_t i = 0; i < signatureSize; i++)
        //     {
        //         c += __builtin_popcountll(a[i]);
        //     }
        //     fprintf(stderr, ", %zu", c);
        // }
        // fprintf(stderr, "\n");
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

    //?
    void get_counts()
    {
        vector<size_t> leaves;
        find_leaves(root, leaves);
        fprintf(stderr, ">>>*\n");

        for (size_t leaf : leaves)
        {
            size_t count = 0;
            for (size_t cnt : counts[leaf])
            {
                if (cnt > 2)
                {
                    count++;
                }
            }
            counts_threshold[leaf] = count;
            fprintf(stderr, "%zu,%zu\n", leaf, count);
        }
        fprintf(stderr, ">>>\n");
    }

    template <class RNG>
    inline void insert(RNG &&rng, const cell_type *signature, vector<size_t> &insertionList, size_t idx)
    {
        // Warning: ALWAYS INSERT THE FIRST NODE SINGLE-THREADED
        // We don't have any protection from this because it would slow everything down to do so
        if (root == numeric_limits<size_t>::max())
        {
            root = getNewNodeIdx(insertionList);
            childCounts[root] = 0;
            isBranchNode[root] = 0;
        }

        // vector<size_t> leaves;
        // find_leaves(root, leaves);

        // for (size_t leaf : leaves)
        // {
        //     size_t count = 0;
        //     for (size_t cnt : counts[leaf])
        //     {
        //         if (cnt > 2)
        //         {
        //             count++;
        //         }
        //     }
        //     counts_threshold[leaf] = count;
        // }
        size_t insertionPoint = traverse(signature);
        // size_t insertionPoint = traverse_all(leaves, signature);
        // size_t insertionPoint = traverse(signature, insertionList);

        // fprintf(stdout, " at %zu\n", insertionPoint);
        omp_set_lock(&locks[insertionPoint]);
        if (childCounts[insertionPoint] < order)
        {
            addSigToMatrix_union(insertionPoint, signature);
            addSigToMatrix_inter(insertionPoint, signature);

            childLinks[insertionPoint].push_back(idx);
            fprintf(stdout, " at %zu\n", insertionPoint);

            childCounts[insertionPoint]++;
            recalculateSig(insertionPoint);
            if (isBranchNode[root])
            {
                recalculateUp(parentLinks[insertionPoint]);
            }
        }
        else
        {
            get_counts();
            fprintf(stdout, " \n");
            splitNode(rng, insertionPoint, signature, insertionList, idx);
            // splitNode(rng, insertionPoint, signature, insertionList, 0);
        }
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