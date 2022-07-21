// There are two kinds of ktree nodes- branch nodes and leaf nodes
// Both contain a signature matrix, plus their own signature
// (the root node signature does not matter and can be blank)
// Branch nodes then contain 'order' links to other nodes
// Leaf nodes do not.
// However, as leaf nodes may become branch nodes, we allocate
// the space anyway.
// As the space to be used is determined at runtime, we use
// parallel arrays, not structs

#ifndef INCLUDE_part_KTree_MULTI2_HPP
#define INCLUDE_part_KTree_MULTI2_HPP

#include <omp.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include "bloom_filter.hpp"
#include "read.hpp"

using namespace std;

// typedef unsigned char cell_type;
bloom_parameters parameters;
size_t signatureSize; // Signature size (depends on element in BF, obtained while read binary)
size_t maxWindow = 0;
size_t partree_order = 10;
size_t HD_threshold = 0;

void toBinary(FILE *stream, cell_type letter)
{
    int binary[bits_per_char];
    for (int n = 0; n < bits_per_char; n++)
        binary[bits_per_char - 1 - n] = (letter >> n) & 1;

    for (int n = 0; n < bits_per_char; n++)
        fprintf(stream, "%d", binary[n]);
}

void toBinaryIdx(FILE *stream, cell_type letter, int offset)
{
    int binary[bits_per_char];
    for (int n = 0; n < bits_per_char; n++)
        binary[bits_per_char - 1 - n] = (letter >> n) & 1;

    for (int n = 0; n < bits_per_char; n++)
    {
        if (binary[n] > 0)
        {
            fprintf(stream, "%d,", n + offset);
        }
    }
}

void printBF(FILE *stream, vector<cell_type> bf, bool idxOnly = true)
{
    if (idxOnly)
    {
        for (int i = 0; i < bf.size(); i++)
        {
            toBinaryIdx(stream, bf[i], i * bits_per_char);
        }
    }
    else
    {
        for (auto val : bf)
        {
            toBinary(stream, val);
        }
        fprintf(stream, "-\n");
    }
}

// only OR bits with count greater than threshold
vector<vector<cell_type>> calcUnion(vector<vector<vector<cell_type>>> sigs, size_t threshold = 0)
{
    // default if no threshold given, half the sigs size
    if (threshold == 0)
    {
        threshold = sigs.size() / 2;
    }

    vector<vector<cell_type>> result;
    for (int w = 0; w < maxWindow; w++)
    {
        vector<cell_type> a(signatureSize);
        for (size_t i = 0; i < signatureSize; i++)
        {
            vector<size_t> c(bits_per_char);
            fill(c.begin(), c.end(), 0);
            for (size_t j = 0; j < sigs.size(); j++)
            {
                // fprintf(stderr, "--%zu\n", sigs[j][w][i]);
                for (int n = 0; n < bits_per_char; n++)
                {
                    if ((sigs[j][w][i] >> n) & 1 == 1)
                    {
                        c[n]++;
                        // fprintf(stderr, "--%d\n", n);
                    }
                }
            }
            for (int n = 0; n < bits_per_char; n++)
            {
                if (c[n] > threshold)
                {
                    a[i] |= (1 << n);
                    // fprintf(stderr, "--**%zu\n", c[n]);
                }
            }
        }
        result.push_back(a);
    }
    return result;
}

size_t calcInter(vector<cell_type> a, vector<cell_type> b)
{
    size_t c = 0;
    for (size_t i = 0; i < a.size(); i++)
    {
        c += __builtin_popcountll(a[i] & b[i]);
    }
    return c;
}

size_t calcUnion(vector<cell_type> a, vector<cell_type> b)
{
    size_t c = 0;
    for (size_t i = 0; i < a.size(); i++)
    {
        c += __builtin_popcountll(a[i] | b[i]);
    }
    return c;
}

size_t calcHD(vector<cell_type> a, vector<cell_type> b)
{
    size_t c = 0;
    for (size_t i = 0; i < a.size(); i++)
    {
        c += __builtin_popcountll(a[i] ^ b[i]);
    }
    return c;
}

size_t countSetBits(vector<cell_type> a)
{
    size_t c = 0;
    for (size_t i = 0; i < a.size(); i++)
    {
        c += __builtin_popcountll(a[i]);
    }
    return c;
}

vector<vector<cell_type>> do_union(vector<vector<cell_type>> a, vector<vector<cell_type>> b)
{
    if (a.size() == 0)
    {
        return b;
    }
    size_t size = min(a.size(), b.size());
    if (a.size() > b.size())
    {
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < a[i].size(); j++)
            {
                a[i][j] |= b[i][j];
            }
        }
        return a;
    }
    else
    {
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < b[i].size(); j++)
            {
                b[i][j] |= a[i][j];
            }
        }
        return b;
    }
}


vector<vector<cell_type>> do_inter(vector<vector<cell_type>> a, vector<vector<cell_type>> b)
{
    if (a.size() == 0)
    {
        return b;
    }
    size_t size = min(a.size(), b.size());
    if (a.size() > b.size())
    {
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < a[i].size(); j++)
            {
                a[i][j] &= b[i][j];
            }
        }
        return a;
    }
    else
    {
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < b[i].size(); j++)
            {
                b[i][j] &= a[i][j];
            }
        }
        return b;
    }
}

size_t compareSeqsHD(vector<vector<cell_type>> a, vector<vector<cell_type>> b)
{
    size_t distance = 0;
    for (int i = 0; i < min(a.size(), b.size()); i++)
    {
        distance += calcHD(a[i], b[i]);
    }
    return distance;
}

size_t compareSeqsWindows(vector<vector<cell_type>> mean, vector<vector<cell_type>> x)
{
    size_t matchingWindows = 0;

    for (int w = 0; w < min(mean.size(), x.size()); w++) // compare window by window
    {
        size_t setBits = countSetBits(x[w]);
        size_t threshold = min(HD_threshold, setBits);
        size_t distance = calcInter(mean[w], x[w]);
        if (distance > threshold)
        {
            matchingWindows++;
        }
    }
    return matchingWindows;
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
vector<vector<vector<cell_type>>> createRandomSigs(RNG &&rng, vector<vector<vector<cell_type>>> sigs, size_t clusterCount = 2)
{
    // constexpr size_t clusterCount = 2;
    vector<vector<vector<cell_type>>> clusterSigs;
    size_t signatureCount = sigs.size();
    uniform_int_distribution<size_t> dist(0, signatureCount - 1);
    bool finished = false;

    for (int n = 0; n < clusterCount; n++)
    {
        vector<vector<cell_type>> uniqueSigs;
        for (size_t i = 0; i < maxWindow; i++)
        {
            size_t sig = dist(rng);
            size_t idx = i;
            if (sigs[sig].size() < idx)
            {
                idx = 0;
            }
            uniqueSigs.push_back(sigs[sig][idx]);
        }

        clusterSigs.push_back(uniqueSigs);
    }

    return clusterSigs;
}

// vector<vector<vector<cell_type>>> createRandomSigs(vector<vector<vector<cell_type>>> sigs, size_t clusterCount = 2)
// {
//     vector<vector<vector<cell_type>>> clusterSigs;
//     for (int n = 0; n < clusterCount; n++)
//     {
//         clusterSigs.push_back(sigs[n * 26]);
//     }

//     return clusterSigs;
// }

void reclusterSignatures(vector<size_t> &clusters, vector<vector<vector<cell_type>>> meanSigs, vector<vector<vector<cell_type>>> sigs, size_t clusterCount = 2)
{
    set<size_t> allClusters;
    for (size_t sig = 0; sig < clusters.size(); sig++)
    {
        vector<vector<cell_type>> sourceSignature = sigs[sig];
        size_t maxInterCluster = 0;
        double maxInter_overall = -1;

        for (size_t cluster = 0; cluster < clusterCount; cluster++)
        {
            // check how many windows above the threshold with cluster sig,
            // pick the one with most window match
            //? pointer &
            vector<vector<cell_type>> clusterSignature = meanSigs[cluster];
            double windows_above_threshold = 0;
            size_t minSize = min(clusterSignature.size(), sourceSignature.size());
            for (size_t i = 0; i < minSize; i++)
            {
                size_t threshold = min(HD_threshold, countSetBits(clusterSignature[i]));
                size_t dist = calcInter(sourceSignature[i], clusterSignature[i]);

                if (dist > threshold)
                {
                    windows_above_threshold++;
                }
            }
            windows_above_threshold = windows_above_threshold / minSize;

            if (windows_above_threshold > maxInter_overall)
            {
                maxInter_overall = windows_above_threshold;
                maxInterCluster = cluster;
            }
        }
        clusters[sig] = maxInterCluster;
        allClusters.insert(maxInterCluster);
    }

    // for (size_t sig = 0; sig < clusters.size(); sig++)
    // {
    //     vector<vector<cell_type>> sourceSignature = sigs[sig];
    //     size_t minHDCluster = 0;
    //     size_t minHD = numeric_limits<size_t>::max();

    //     for (size_t cluster = 0; cluster < clusterCount; cluster++)
    //     {
    //         // check how many windows above the threshold with cluster sig,
    //         // pick the one with most window match
    //         //? pointer &
    //         vector<vector<cell_type>> clusterSignature = meanSigs[cluster];
    //         size_t HD = compareSeqsHD(clusterSignature, sourceSignature);

    //         if (HD < minHD)
    //         {
    //             minHD = HD;
    //             minHDCluster = cluster;
    //         }
    //     }
    //     clusters[sig] = minHDCluster;
    //     allClusters.insert(minHDCluster);
    // }

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

void reclusterSignaturesHD(vector<size_t> &clusters, vector<vector<vector<cell_type>>> meanSigs, vector<vector<vector<cell_type>>> sigs, size_t clusterCount = 2)
{
    set<size_t> allClusters;

    for (size_t sig = 0; sig < clusters.size(); sig++)
    {
        vector<vector<cell_type>> sourceSignature = sigs[sig];
        size_t minHDCluster = 0;
        size_t minHD = numeric_limits<size_t>::max();

        for (size_t cluster = 0; cluster < clusterCount; cluster++)
        {
            // check how many windows above the threshold with cluster sig,
            // pick the one with most window match
            //? pointer &
            vector<vector<cell_type>> clusterSignature = meanSigs[cluster];
            size_t HD = compareSeqsHD(clusterSignature, sourceSignature);

            if (HD < minHD)
            {
                minHD = HD;
                minHDCluster = cluster;
            }
        }
        clusters[sig] = minHDCluster;
        allClusters.insert(minHDCluster);
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

vector<vector<vector<cell_type>>> createClusterSigs(const vector<vector<size_t>> &clusterLists, vector<vector<vector<cell_type>>> sigs)
{
    vector<vector<vector<cell_type>>> clusterSigs;

    // //#pragma omp parallel
    // {
    //     size_t n = 0;
    //     //#pragma omp for
    //     for (size_t cluster = 0; cluster < clusterLists.size(); cluster++)
    //     {
    //         vector<vector<cell_type>> flattenedSignature_union; //? windows, bf
    //         if (clusterLists[cluster].size() == 0)
    //         {
    //             for (int i = 0; i < maxWindow; i++)
    //             {
    //                 vector<cell_type> bf(signatureSize);
    //                 flattenedSignature_union.push_back(bf);
    //             }
    //         }

    //         else
    //         {
    //             flattenedSignature_union = sigs[clusterLists[cluster][0]];
    //         }

    //         for (size_t signature : clusterLists[cluster])
    //         {
    //             n++;
    //             vector<vector<cell_type>> signatureData = sigs[signature];
    //             for (size_t i = 1; i < maxWindow; i++)
    //             {
    //                 vector<cell_type> window = signatureData[i];

    //                 for (size_t j = 0; j < signatureSize; j++)
    //                 {
    //                     flattenedSignature_union[i][j] |= window[j];
    //                 }
    //             }
    //             clusterSigs.push_back(flattenedSignature_union);
    //         }

    //         // // debug
    //         // for (size_t i = 0; i < maxWindow; i++)
    //         // {

    //         //     printBF(stderr, flattenedSignature_union[i]);
    //         //     fprintf(stderr, "\n");
    //         //     // for (size_t signature : clusterLists[cluster])
    //         //     // {
    //         //     //     vector<vector<cell_type>> signatureData = sigs[signature];
    //         //     //     printBF(stderr, signatureData[i]);
    //         //     //     fprintf(stderr, "\n");
    //         //     // }
    //         // }
    //     }
    // }

    // // // for (auto clu : clusterSigs)
    // // // {
    // // //     for (auto bf : clu)
    // // //     {
    // // //         printBF(stderr, bf);
    // // //     }
    // // //     fprintf(stderr, "\n");
    // // // }

    for (size_t cluster = 0; cluster < clusterLists.size(); cluster++)
    {
        vector<vector<vector<cell_type>>> cluSeqs;
        for (size_t signature : clusterLists[cluster])
        {
            cluSeqs.push_back(sigs[signature]);
        }
        clusterSigs.push_back(calcUnion(cluSeqs));
    }
    return clusterSigs;
}

vector<double> calcStdev(const vector<vector<size_t>> &clusterLists, vector<vector<vector<cell_type>>> meanSigs, vector<vector<vector<cell_type>>> sigs, size_t clusterCount)
{
    vector<double> stdevs(clusterCount);
#pragma omp parallel
    {
#pragma omp for
        for (size_t cluster = 0; cluster < clusterLists.size(); cluster++)
        {
            vector<vector<cell_type>> clusterSignature = meanSigs[cluster];
            size_t squareTotal = 0;

            // calc std
            for (size_t signature : clusterLists[cluster])
            {
                vector<vector<cell_type>> sourceSignature = sigs[signature];
                size_t hd = compareSeqsHD(sourceSignature, clusterSignature);
                squareTotal += hd * hd;
            }
            double stdev = sqrt(squareTotal / (double)(clusterLists[cluster].size()));
            stdevs[cluster] = stdev;
        }
    }
    return stdevs;
}

double indexDB(const vector<vector<size_t>> &clusterLists, vector<vector<vector<cell_type>>> meanSigs, vector<vector<vector<cell_type>>> sigs, size_t clusterCount)
{
    vector<double> stdevs = calcStdev(clusterLists, meanSigs, sigs, clusterCount);
    double DB;

    // calc cluster similarity
    for (int i = 0; i < clusterCount; i++)
    {
        double maxSim = 0;
        for (int j = 0; j < clusterCount; j++)
        {
            if (i != j)
            {
                double similarity = (stdevs[i] + stdevs[j]) / compareSeqsHD(meanSigs[i], meanSigs[j]);
                if (similarity > maxSim)
                {
                    maxSim = similarity;
                }
            }
        }
        DB += maxSim;
    }
    DB = DB / clusterCount;
    return DB;
}

template <class RNG>
void kMeans(RNG &&rng, vector<vector<vector<cell_type>>> sigs, size_t clusterCount)
{
    size_t nodeSigs = sigs.size();
    vector<vector<vector<cell_type>>> meanSigs = createRandomSigs(rng, sigs, clusterCount);
    // vector<vector<vector<cell_type>>> meanSigs = createRandomSigs(sigs, clusterCount);
    vector<size_t> clusters(nodeSigs);
    vector<vector<size_t>> clusterLists;

    for (int iteration = 0; iteration < 4; iteration++)
    {
        // reclusterSignatures(clusters, meanSigs, sigs, clusterCount);
        reclusterSignaturesHD(clusters, meanSigs, sigs, clusterCount);
        clusterLists = createClusterLists(clusters, clusterCount);
        meanSigs = createClusterSigs(clusterLists, sigs);
    }

    auto DB = indexDB(clusterLists, meanSigs, sigs, clusterCount);
    fprintf(stderr, "DB index: %f\n", DB);

    // for (int i = 0; i < clusterCount; i++)
    // {
    //     fprintf(stderr, ">");
    //     for (size_t seq : clusterLists[i])
    //     {
    //         fprintf(stderr, ",%zu", seq);
    //     }
    //     fprintf(stderr, "\n");
    // }

    for (int i = 0; i < nodeSigs; i++)
    {
        fprintf(stdout, "%d,%zu\n", i, clusters[i]);
    }
}

class part_KTree_multi2
{
public:
    size_t root = 0;                                    // # of root node
    vector<size_t> childCounts;                         // n entries, number of children
    vector<int> isBranchNode;                           // n entries, is this a branch node
    vector<int> isLeafNode;                             // n entries, is this a leaf node
    vector<vector<size_t>> childLinks;                  // n * o entries, links to children
    vector<vector<size_t>> dest;                        // n * o entries, links to seq idx
    vector<size_t> parentLinks;                         // n entries, links to parents
    vector<vector<sig_type>> matri;                     // capacity * signatureSize * n
    vector<vector<vector<vector<cell_type>>>> matrices; // list of BF in each node
    vector<vector<vector<cell_type>>> overallSigs;      // BF of windows under a branch, store in leaf
    vector<size_t> leafNodes;                           // n entries
    vector<vector<vector<cell_type>>> meanSigs;         // mean of each node
    vector<omp_lock_t> locks;                           // n locks
    size_t capacity = 0;                                // Set during construction, currently can't change

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
                meanSigs.resize(capacity);
            }
            //#pragma omp single
            {
                overallSigs.resize(capacity);
            }
        }
    }

    // each node is a window, represented by a BF (vector<cell_type>) that stores a number of minimisers
    void addSigToMatrix(size_t node, vector<vector<cell_type>> sig)
    {
        // matrices[node] = sig;

        // for (int i = 0; i < sig.size(); i++)
        // {
        //     matrices[node][i] |= sig[i];
        // }
        matrices[node].push_back(sig);
    }

    part_KTree_multi2(size_t capacity, size_t threshold)
    {
        reserve(capacity);
        childCounts[root] = 0;
        isBranchNode[root] = 0;
        // vector<cell_type> temp;
        // addSigToMatrix(root, temp);
        // // matrices[root] = temp;
        HD_threshold = threshold;
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

    //?
    void printMatrix(size_t node, bool idxOnly = false)
    {
        if (idxOnly)
        {
            for (auto sig : meanSigs[node])
            {
                printBF(stdout, sig, idxOnly);
                fprintf(stdout, "</br>");
            }
        }
        else
        {
            for (auto sig : meanSigs[node])
            {
                for (cell_type val : sig)
                {
                    if (val > (cell_type)0)
                    {
                        fprintf(stdout, "%d,", val);
                    }
                }
                fprintf(stdout, "</br>");
            }
        }
    }

    void printNodeJson(size_t tnode)
    {
        fprintf(stdout, "{\"node\":\"%zu\",\"content\":\"*", tnode);
        printMatrix(tnode, true);
        fprintf(stdout, "\",\"children\":[");
    }

    void printSubTreeJson(size_t tnode)
    {
        // fprintf(stderr, "%zu, %zu, %zu, %zu\n", tnode, isLeafNode[tnode], childCounts[tnode], childLinks[tnode].size());
        if (!isLeafNode[tnode])
        {
            printNodeJson(tnode);
            if (childLinks[tnode].size() > 0)
            {
                printSubTreeJson(childLinks[tnode][0]);

                for (size_t i = 1; i < childCounts[tnode]; i++)
                {
                    fprintf(stdout, ",");
                    printSubTreeJson(childLinks[tnode][i]);
                }
            }
        }
        else
        {
            fprintf(stdout, "{\"node\":\"%zu\",\"content\":\"*", tnode);
            printMatrix(tnode, true);
            fprintf(stdout, "\",\"children\":[{\"node\":\">>>");
            for (size_t seq : dest[tnode])
            {
                fprintf(stdout, "%zu,", seq);
            }
            fprintf(stdout, "\"}");

            // a leaf node can also be a branch if there is a shorter seqs that share the same set of minimisers
            if (isBranchNode[tnode])
            {
                for (size_t child : childLinks[tnode])
                {
                    fprintf(stdout, ",");
                    printSubTreeJson(child);
                }
            }
        }
        fprintf(stdout, "]}");
    }

    //?
    void updateMeanSig(size_t node)
    {
        // vector<vector<cell_type>> mean = matrices[node][0];

        // size_t BFsize = mean[0].size();
        // for (auto sig : matrices[node])
        // {
        //     for (size_t w = 0; w < min(sig.size(), mean.size()); w++)
        //     {
        //         // printBF(stderr, sig);
        //         // fprintf(stderr, "-\n");
        //         for (int i = 0; i < BFsize; i++)
        //         {
        //             mean[w][i] |= sig[w][i];
        //         }
        //     }
        // }
        // meanSigs[node] = mean;

        size_t threshold = matrices[node].size() / 2;

        size_t minW = numeric_limits<size_t>::max();
        size_t shortest = 0;
        for (int i = 0; i < matrices[node].size(); i++)
        {
            size_t temp = matrices[node][i].size();
            if (temp < minW)
            {
                minW = temp;
                shortest = i;
            }
        }
        vector<vector<cell_type>> mean = matrices[node][shortest];
        size_t size = mean[0].size();

        //? do something to window size
        for (int w = 0; w < mean.size(); w++)
        {
            for (size_t i = 0; i < size; i++)
            {
                vector<size_t> c(bits_per_char);
                fill(c.begin(), c.end(), 0);
                for (size_t j = 0; j < matrices[node].size(); j++)
                {
                    // fprintf(stderr, "--%zu\n", sigs[j][w][i]);
                    for (int n = 0; n < bits_per_char; n++)
                    {
                        if ((matrices[node][j][w][i] >> n) & 1 == 1)
                        {
                            c[n]++;
                            // fprintf(stderr, "--%d\n", n);
                        }
                    }
                }
                for (int n = 0; n < bits_per_char; n++)
                {
                    if (c[n] > threshold)
                    {
                        mean[w][i] |= (1 << n);
                        // fprintf(stderr, "--**%zu\n", c[n]);
                    }
                }
            }
        }

        meanSigs[node] = mean;
    }

    inline size_t first_insert(vector<vector<cell_type>> signature, vector<size_t> &insertionList, size_t idx, bool store = true)
    {
        size_t node = getNewNodeIdx(insertionList);
        size_t parent = root;

        while (isBranchNode[parent])
        {
            size_t maxWindows = signature.size() / 2;
            size_t bestLeaf = 0;

            // for (size_t child : childLinks[parent])
            for (int i = 0; i < childCounts[parent]; i++)
            {
                size_t child = childLinks[parent][i];
                size_t matching = compareSeqsWindows(meanSigs[child], signature);
                if (matching > maxWindows)
                {
                    maxWindows = matching;
                    bestLeaf = child;
                }
            }
            if (bestLeaf != 0)
            {
                parent = bestLeaf;
            }
            else
            {
                fprintf(stderr, "not close\n");
                break;
            }
        }

        parentLinks[node] = parent;
        meanSigs[parent] = do_union(meanSigs[parent], signature);
        childLinks[parent].push_back(node);
        isBranchNode[parent] = 1;
        childCounts[parent]++;
        isLeafNode[node] = 1;
        overallSigs[node] = signature;
        addSigToMatrix(parent, signature);
        leafNodes.push_back(node);

        if (store)
        {
            dest[node].push_back(idx);
            fprintf(stderr, "%zu\n", idx);
        }

        return node;
    }

    inline size_t grow(size_t parent, vector<vector<cell_type>> signature, vector<size_t> &insertionList, size_t idx, bool store = true)
    {
        size_t node = getNewNodeIdx(insertionList);
        // grow new root;
        if (parent == root)
        {
            fprintf(stderr,"grow root\n");
            size_t newRoot = getNewNodeIdx(insertionList);
            
            

            meanSigs[newRoot] = meanSigs[root]; // will be union later
            childCounts[newRoot]=1; // will be added later
            childLinks[newRoot].push_back(root);
            addSigToMatrix(newRoot, meanSigs[root]);
            addSigToMatrix(newRoot, signature);

            root = newRoot;
            parent = root;
        }

        parentLinks[node] = parent;
        meanSigs[parent] = do_inter(meanSigs[parent], signature);
        childLinks[parent].push_back(node);
        isBranchNode[parent] = 1;
        childCounts[parent]++;
        isLeafNode[node] = 1;
        overallSigs[node] = signature;
        addSigToMatrix(parent, signature);
        leafNodes.push_back(node);

        if (store)
        {
            dest[node].push_back(idx);
            fprintf(stderr, "%zu\n", idx);
        }

        return node;
    }

    int traverseRoot(vector<vector<cell_type>> signature)
    {
        size_t node = root;
        size_t maxWindows = signature.size() / 2;
        size_t bestLeaf = 0;

        while (isBranchNode[node])
        {
            size_t oldBest = bestLeaf;
            maxWindows = signature.size() / 2;
            for (size_t leaf : childLinks[node])
            {
                size_t matching = compareSeqsWindows(meanSigs[leaf], signature);
                if (matching > maxWindows)
                {
                    maxWindows = matching;
                    bestLeaf = leaf;
                }
            }
            // nothing better in children, escape
            if (oldBest == bestLeaf)
            {
                fprintf(stderr, "---%zu\n", bestLeaf);
                return -bestLeaf;
            }
            node = bestLeaf;
        }
        return bestLeaf;
    }

    size_t traverse(vector<vector<cell_type>> signature)
    {
        size_t maxWindows = signature.size() / 2;
        size_t bestLeaf = 0;

        for (size_t leaf : leafNodes)
        {
            size_t matching = compareSeqsWindows(overallSigs[leaf], signature);
            if (matching > maxWindows)
            {
                maxWindows = matching;
                bestLeaf = leaf;
            }
        }
        return bestLeaf;
    }

    size_t traverseHD(vector<vector<cell_type>> signature)
    {
        size_t minHD = numeric_limits<size_t>::max();
        size_t bestLeaf = 0;

        for (size_t leaf : leafNodes)
        {
            size_t HD = compareSeqsHD(overallSigs[leaf], signature);
            if (HD < minHD)
            {
                minHD = HD;
                bestLeaf = leaf;
            }
        }
        return bestLeaf;
    }

    size_t traverseHDRoot(vector<vector<cell_type>> signature)
    {
        size_t node = root;

        while (isBranchNode[node])
        {
            size_t minHD = numeric_limits<size_t>::max();
            size_t bestLeaf = 0;
            for (size_t leaf : childLinks[node])
            {
                size_t HD = compareSeqsHD(meanSigs[leaf], signature);
                if (HD < minHD)
                {
                    minHD = HD;
                    bestLeaf = leaf;
                }
            }
            node = bestLeaf;
        }
        return node;
    }

    // lock parent while updating matrix, return parent node for unlocking later
    size_t updateParentMatrix(size_t node, vector<vector<cell_type>> meanSig)
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

        //?
        matrices[parent][idx] = meanSig;

        return parent;
    }

    void recalculateUp(size_t node)
    {
        size_t limit = 10;
        // fprintf(stderr, "RecalculateUp %zu\n", node);
        while (node != root)
        {
            updateMeanSig(node);

            // update parent matrix with new mean
            size_t parent = updateParentMatrix(node, meanSigs[node]);
            omp_unset_lock(&locks[parent]);

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
    }

    inline size_t insert(vector<vector<cell_type>> signature, vector<size_t> &insertionList, size_t idx, bool store = true)
    {

        // size_t node = traverse(signature);
        int node = traverseRoot(signature);

        if (node == 0)
        {
            // fprintf(stderr, "here\n");
            // node = first_insert(signature, insertionList, idx, !store);

            node = grow(node,signature, insertionList, idx, !store);
        }
        else if (node < 0)
        {
            fprintf(stderr, "here\n");
            node = grow(-node, signature, insertionList, idx, !store);
        }

        if (store)
        {
            dest[node].push_back(idx);
        }
        return node;
    }

    template <class RNG>
    void splitNode(RNG &&rng, size_t node, vector<vector<cell_type>> signature, vector<size_t> &insertionList, size_t link)
    {
        size_t nodeSigs = childCounts[node] + 1;
        vector<vector<vector<cell_type>>> sigs = matrices[node];
        sigs.push_back(signature);

        size_t clusterCount = 2;
        vector<vector<vector<cell_type>>> kmeanSigs = createRandomSigs(rng, sigs, clusterCount);
        vector<size_t> clusters(nodeSigs);
        vector<vector<size_t>> clusterLists;
        for (int iteration = 0; iteration < 4; iteration++)
        {
            // fprintf(stderr, "Iteration %d\n", iteration);
            reclusterSignatures(clusters, kmeanSigs, sigs);
            clusterLists = createClusterLists(clusters);
            kmeanSigs = createClusterSigs(clusterLists, sigs);
        }

        // for (auto cluster : clusterLists)
        // {
        //     fprintf(stderr, "\n>\n");
        //     for (auto val : cluster)
        //     {
        //         fprintf(stderr, "%zu,", val);
        //     }
        // }

        // Create the sibling node
        size_t sibling = getNewNodeIdx(insertionList);

        size_t newlyAddedIdx = childCounts[node];

        childCounts[sibling] = clusterLists[1].size();
        isBranchNode[sibling] = isBranchNode[node];
        isLeafNode[sibling] = isLeafNode[node];
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
                addSigToMatrix(sibling, sigs[seqIdx]);
                siblingIdx++;
            }
        }

        // memcpy(&meanSigs[node], &kmeanSigs[0], sizeof(kmeanSigs[0]));
        // memcpy(&meanSigs[sibling], &kmeanSigs[1], sizeof(kmeanSigs[0]));

        //?
        meanSigs[node] = kmeanSigs[0];
        meanSigs[sibling] = kmeanSigs[1];

        //     printBF(stderr, w);
        // for (auto w : meanSigs[0])
        // {
        //     printBF(stderr, w);
        //     fprintf(stderr, ">");
        // }
        // fprintf(stderr, "\n");

        // for (auto w : meanSigs[sibling])
        // {
        //     printBF(stderr, w);
        //     fprintf(stderr, ">");
        // }
        // fprintf(stderr, "\n");

        childCounts[node] = clusterLists[0].size();
        // Fill the current node with the other cluster of signatures
        {
            matrices[node].clear();
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
                addSigToMatrix(node, sigs[seqIdx]);
                nodeIdx++;
            }
            if (isBranchNode[node])
            {
                childLinks[node].resize(childCounts[node]);
            }
        }

        // Is this the root level?
        if (node == root)
        {
            // fprintf(stderr, "Node being split is root node\n");

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
            addSigToMatrix(newRoot, kmeanSigs[0]);
            addSigToMatrix(newRoot, kmeanSigs[1]);

            updateMeanSig(newRoot);

            root = newRoot;

            // recalculateSig(root);
        }
        else
        {
            // First, update the reference to this node in the parent with the new mean
            size_t parent = updateParentMatrix(node, kmeanSigs[0]);

            // fprintf(stderr, "Split again, r: %zu, n: %zu, s:%zu\n", parentLinks[node], node, sibling);

            // Connect sibling node to parent
            parentLinks[sibling] = parent;

            // Now add a link in the parent node to the sibling node
            if (childCounts[parent] + 1 < partree_order)
            {

                childLinks[parent].push_back(sibling);

                childCounts[parent]++;
                addSigToMatrix(parent, kmeanSigs[1]);

                // Update signatures (may change?)
                recalculateUp(parent);
            }
            else
            {
                // fprintf(stderr, ">Split again, r: %zu, n: %zu, s:%zu\n", parentLinks[node], node, sibling);
                splitNode(rng, parent, kmeanSigs[1], insertionList, sibling);
            }
            // Unlock the parent
            omp_unset_lock(&locks[parent]);
        }

        // fprintf(stderr, "Split finished, r: %zu, n: %zu, s:%zu\n", parentLinks[node], node, sibling);
    }

    template <class RNG>
    void insert_split(RNG &&rng, vector<vector<cell_type>> signature, vector<size_t> &insertionList)
    {
        size_t insertionPoint = traverseHD(signature);
        omp_set_lock(&locks[insertionPoint]);

        // fprintf(stderr, " at %zu\n", insertionPoint);

        if (childCounts[insertionPoint] > partree_order)
        {
            splitNode(rng, insertionPoint, signature, insertionList, 0);
        }
        else
        {
            addSigToMatrix(insertionPoint, signature);
            updateMeanSig(insertionPoint);
            childCounts[insertionPoint]++;
        }
        omp_unset_lock(&locks[insertionPoint]);
    }

    inline size_t insert_again(vector<vector<cell_type>> signature, size_t idx)
    {

        // vector<size_t> matches;
        // vector<size_t> depths;
        // size_t ans = traverse_searchHD(root, signature, 0, matches, depths);

        // auto max = max_element(depths.begin(), depths.end());
        // auto i = distance(depths.begin(), max);

        // size_t node = matches[i];
        // size_t parent = node;
        // size_t counter = depths[i];
        // // fprintf(stderr, "---%zu,%zu,%zu\n", parent, counter, signature.size() - 1);

        // if (counter == signature.size())
        // {
        //     fprintf(stderr, "---%zu,%zu,%zu\n", idx, counter, signature.size());
        // }

        size_t node = traverseHDRoot(signature);
        // size_t node = traverseHD(signature);
        // fprintf(stderr, "---%zu,%zu\n", idx, node);

        dest[node].push_back(idx);

        return node;
    }

    // for every leaf, output node idx;cluster size; seq idx
    void outputClusterbyLeaf(FILE *pFile, size_t node)
    {
        if (!isLeafNode[node])
        {
            for (size_t child : childLinks[node])
            {
                outputClusterbyLeaf(pFile, child);
            }
        }
        else
        {
            fprintf(pFile, "%zu;%zu;", node, dest[node].size());
            for (size_t seq : dest[node])
            {
                fprintf(pFile, "%zu,", seq);
            }
            fprintf(pFile, "\n");

            if (isBranchNode[node])
            {
                for (size_t child : childLinks[node])
                {
                    outputClusterbyLeaf(pFile, child);
                }
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