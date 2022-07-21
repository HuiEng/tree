// There are two kinds of ktree nodes- branch nodes and leaf nodes
// Both contain a signature matrix, plus their own signature
// (the root node signature does not matter and can be blank)
// Branch nodes then contain 'order' links to other nodes
// Leaf nodes do not.
// However, as leaf nodes may become branch nodes, we allocate
// the space anyway.
// As the space to be used is determined at runtime, we use
// parallel arrays, not structs

#ifndef INCLUDE_part_KTree_MULTI_HPP
#define INCLUDE_part_KTree_MULTI_HPP

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

vector<cell_type> do_union(vector<cell_type> a, vector<cell_type> b)
{
    vector<cell_type> u(a.size());
    // printBF(stderr, a);
    // fprintf(stderr, "\n");
    // printBF(stderr, b);
    // fprintf(stderr, "\n");

    for (int i = 0; i < a.size(); i++)
    {
        u[i] = a[i] | b[i];
    }
    // printBF(stderr, u);
    // fprintf(stderr, "\n");
    return u;
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

vector<vector<vector<cell_type>>> createRandomSigs(vector<vector<vector<cell_type>>> sigs, size_t clusterCount = 2)
{
    vector<vector<vector<cell_type>>> clusterSigs;
    for (int n = 0; n < clusterCount; n++)
    {
        clusterSigs.push_back(sigs[n * 26]);
    }

    return clusterSigs;
}

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
                //// fprintf(stderr,"%zu,%zu\n",i,threshold);
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

class part_KTree_multi
{
public:
    size_t root = 0;                               // # of root node
    vector<size_t> childCounts;                    // n entries, number of children
    vector<int> isBranchNode;                      // n entries, is this a branch node
    vector<int> isLeafNode;                        // n entries, is this a leaf node
    vector<vector<size_t>> childLinks;             // n * o entries, links to children
    vector<vector<size_t>> dest;                   // n * o entries, links to seq idx
    vector<size_t> parentLinks;                    // n entries, links to parents
    vector<vector<sig_type>> matri;                // capacity * signatureSize * n
    vector<vector<vector<cell_type>>> matrices;    // list of BF in each node
    vector<vector<vector<cell_type>>> overallSigs; // BF of windows under a branch, store in leaf
    vector<size_t> leafNodes;                      // n entries
    vector<vector<cell_type>> meanSigs;            // mean of each node
    vector<omp_lock_t> locks;                      // n locks
    size_t capacity = 0;                           // Set during construction, currently can't change

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
    void addSigToMatrix(size_t node, vector<cell_type> sig)
    {
        // matrices[node] = sig;

        // for (int i = 0; i < sig.size(); i++)
        // {
        //     matrices[node][i] |= sig[i];
        // }
        matrices[node].push_back(sig);
    }

    part_KTree_multi(size_t capacity, size_t threshold)
    {
        reserve(capacity);
        childCounts[root] = 0;
        isBranchNode[root] = 0;
        vector<cell_type> temp;
        addSigToMatrix(root, temp);
        // matrices[root] = temp;
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

    void printMatrix(size_t node, bool idxOnly = false)
    {
        if (idxOnly)
        {
            printBF(stdout, meanSigs[node], idxOnly);
        }
        else
        {
            for (cell_type val : meanSigs[node])
            {
                if (val > (cell_type)0)
                {
                    fprintf(stdout, "%d,", val);
                }
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
        if (!isLeafNode[tnode])
        {
            printNodeJson(tnode);
            printSubTreeJson(childLinks[tnode][0]);

            for (size_t i = 1; i < childCounts[tnode]; i++)
            {
                fprintf(stdout, ",");
                printSubTreeJson(childLinks[tnode][i]);
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

    void Trim(vector<size_t> &insertionList)
    {
        for (size_t i = 0; i < insertionList[0]; i++)
        {
            if (isLeafNode[i] && dest[i].size() == 0)
            {
                size_t node = i;
                size_t parent = parentLinks[node];

                while (parent != root)
                {
                    insertionList.push_back(node); // reuse redundant node
                    auto p = find(childLinks[parent].begin(), childLinks[parent].end(), node);
                    childLinks[parent].erase(p);
                    childCounts[parent]--;
                    node = parent;
                    parent = parentLinks[parent];
                }

                // insertionList.push_back(node); // reuse redundant node
                // auto p = find(childLinks[parent].begin(), childLinks[parent].end(), node);
                // childLinks[parent].erase(p);
                // childCounts[parent]--;
            }
        }
    }

    void updateMeanSig(size_t node)
    {
        vector<cell_type> mean = matrices[node][0];
        for (auto sig : matrices[node])
        {
            // printBF(stderr, sig);
            // fprintf(stderr, "-\n");
            for (int i = 0; i < mean.size(); i++)
            {
                mean[i] |= sig[i];
            }
        }

        // size_t threshold = matrices[node].size() / 2;
        // size_t size = matrices[node][0].size();
        // vector<cell_type> mean(size);
        // for (size_t i = 0; i < size; i++)
        // {
        //     vector<size_t> c(bits_per_char);
        //     fill(c.begin(), c.end(), 0);
        //     for (size_t j = 0; j < matrices[node].size(); j++)
        //     {
        //         // fprintf(stderr, "--%zu\n", sigs[j][w][i]);
        //         for (int n = 0; n < bits_per_char; n++)
        //         {
        //             if ((matrices[node][j][i] >> n) & 1 == 1)
        //             {
        //                 c[n]++;
        //                 // fprintf(stderr, "--%d\n", n);
        //             }
        //         }
        //     }
        //     for (int n = 0; n < bits_per_char; n++)
        //     {
        //         if (c[n] > threshold)
        //         {
        //             mean[i] |= (1 << n);
        //             // fprintf(stderr, "--**%zu\n", c[n]);
        //         }
        //     }
        // }

        meanSigs[node] = mean;
        // fprintf(stderr, ">>>\n");
        // printBF(stderr, meanSigs[node]);
        // fprintf(stderr, "\n");
    }

    size_t traverseHD(size_t parent, vector<vector<cell_type>> signature, size_t counter, vector<size_t> &candidates, vector<size_t> &max_depths)
    {
        vector<size_t> scores;
        vector<size_t> next_match(childCounts[parent]);
        vector<size_t> maxInter_children;
        size_t minHD = numeric_limits<size_t>::max();

        // printBF(stderr, signature[counter]);
        // fprintf(stderr, "\n");

        // find max score
        for (size_t i = 0; i < childCounts[parent]; i++)
        {
            size_t child = childLinks[parent][i];
            size_t distance = calcHD(meanSigs[child], signature[counter]);

            if (distance < minHD)
            {
                minHD = distance;
            }

            // // bias for next window match
            // if (nodeBits < HD_threshold * 3 && isBranchNode[child] && counter < signature.size() - 1)
            // {
            //     for (size_t grandchild : childLinks[child])
            //     {
            //         if (calcHD(meanSigs[grandchild], signature[counter + 1]) < minHD)
            //         {
            //             next_match[i] = 1;
            //             // fprintf(stderr, "next match %zu\n", grandchild);
            //             break;
            //         }
            //     }
            // }
            scores.push_back(distance);

            // printBF(stderr, meanSigs[child]);
            // fprintf(stderr, "--%zu\n", distance);
        }

        // fprintf(stderr, "distance,%zu\n", maxInter);

        if (minHD < HD_threshold)
        {
            for (size_t i = 0; i < childCounts[parent]; i++)
            {
                if (scores[i] == minHD)
                {
                    size_t child = childLinks[parent][i];

                    if (counter == signature.size() - 1)
                    {
                        {
                            candidates.push_back(child);
                            max_depths.push_back(counter + 1);
                        }
                    }

                    else
                    {
                        traverseHD(child, signature, counter + 1, candidates, max_depths);
                    }
                }
                else if (counter > 0)
                {
                    // reached mismatch, only match up to parent;
                    candidates.push_back(parent);
                    max_depths.push_back(counter - 1);
                }
            }
        }

        // no match
        candidates.push_back(parent);
        max_depths.push_back(counter);

        // return 0;
    }

    // find at which window mismatch starts
    size_t traverse_searchHD(size_t parent, vector<vector<cell_type>> signature, size_t counter, vector<size_t> &candidates, vector<size_t> &max_depths)
    {
        size_t minHD = numeric_limits<size_t>::max();
        vector<size_t> scores;
        vector<size_t> next_match(childCounts[parent]);
        vector<size_t> minHD_children;

        // printBF(stderr, signature[counter]);
        // fprintf(stderr, "\n");

        // find max score
        for (size_t i = 0; i < childCounts[parent]; i++)
        {
            size_t child = childLinks[parent][i];
            size_t distance = calcHD(meanSigs[child], signature[counter]);

            if (distance < minHD)
            {
                minHD = distance;
            }
            scores.push_back(distance);

            // printBF(stderr, meanSigs[child]);
            // fprintf(stderr, "--%zu\n", distance);
        }

        // fprintf(stderr, "distance,%zu\n", maxInter);

        for (size_t i = 0; i < childCounts[parent]; i++)
        {
            if (scores[i] == minHD)
            {
                size_t child = childLinks[parent][i];

                if (counter == signature.size() - 1)
                {
                    candidates.push_back(child);
                    max_depths.push_back(counter + 1);
                }

                else
                {
                    traverse(child, signature, counter + 1, candidates, max_depths);
                }
            }
            else if (counter > 0)
            {
                // reached mismatch, only match up to parent;
                candidates.push_back(parent);
                max_depths.push_back(counter - 1);
            }
        }

        // no match
        candidates.push_back(parent);
        max_depths.push_back(counter);
    }

    // find at which window mismatch starts
    size_t traverse_searchHD(size_t parent, vector<vector<cell_type>> signature)
    {
        size_t minHD_child = 0;
        size_t counter = 0;

        while (isBranchNode[parent] && counter < signature.size())
        {
            size_t minHD = numeric_limits<size_t>::max();
            // find max score
            for (size_t i = 0; i < childCounts[parent]; i++)
            {
                size_t child = childLinks[parent][i];
                size_t distance = calcHD(meanSigs[child], signature[counter]);

                if (distance < minHD)
                {
                    minHD = distance;
                    minHD_child = child;
                }
            }
            parent = minHD_child;
            counter++;
        }
        return parent;
    }

    // find at which window mismatch starts
    size_t traverse2(size_t parent, vector<vector<cell_type>> signature, size_t counter, vector<size_t> &candidates, vector<size_t> &max_depths, vector<size_t> &HDs, size_t totalHD)
    {
        size_t setBits = countSetBits(signature[counter]);
        size_t maxInter = min(HD_threshold, setBits);
        vector<size_t> scores;

        // find max score
        for (size_t i = 0; i < childCounts[parent]; i++)
        {
            size_t child = childLinks[parent][i];
            size_t distance = calcInter(meanSigs[child], signature[counter]);

            if (distance > maxInter)
            {
                maxInter = distance;
            }

            scores.push_back(distance);
        }
        // fprintf(stderr, "distance,%zu\n", maxInter);

        totalHD += maxInter;

        for (size_t i = 0; i < childCounts[parent]; i++)
        {
            if (scores[i] == maxInter)
            {
                size_t child = childLinks[parent][i];

                size_t BFdensity = countSetBits(meanSigs[child]);

                if (counter == signature.size() - 1)
                {
                    candidates.push_back(child);
                    max_depths.push_back(counter + 1);
                    HDs.push_back(totalHD);
                }
                else if (BFdensity > setBits * 2)
                {
                    candidates.push_back(parent);
                    max_depths.push_back(counter);
                    HDs.push_back(totalHD);
                    fprintf(stderr, "here\n");
                }

                else
                {
                    traverse2(child, signature, counter + 1, candidates, max_depths, HDs, totalHD);
                }
            }
            else if (counter > 0)
            {
                // reached mismatch, only match up to parent;
                candidates.push_back(parent);
                max_depths.push_back(counter - 1);
                HDs.push_back(totalHD);
            }
        }

        // no match
        candidates.push_back(parent);
        max_depths.push_back(counter);
        HDs.push_back(totalHD);
    }

    // find at which window mismatch starts
    size_t traverse(size_t parent, vector<vector<cell_type>> signature, size_t counter, vector<size_t> &candidates, vector<size_t> &max_depths)
    {
        size_t setBits = countSetBits(signature[counter]);
        size_t maxInter = min(HD_threshold, setBits);
        // size_t maxInter = max(HD_threshold, setBits);
        // if (counter < setBits)
        //     maxInter = max(HD_threshold, setBits - counter / 2 - 1);

        // if (counter == signature.size() - 1)
        // {
        //     fprintf(stderr, "%zu\n", maxInter);
        // }
        // if (parent == root)
        // {
        //     maxInter = setBits - 1;
        // }
        vector<size_t> scores;
        vector<size_t> next_match(childCounts[parent]);
        vector<size_t> maxInter_children;

        // printBF(stderr, signature[counter]);
        // fprintf(stderr, "\n");

        // find max score
        for (size_t i = 0; i < childCounts[parent]; i++)
        {
            size_t child = childLinks[parent][i];
            size_t nodeBits = countSetBits(meanSigs[child]);
            size_t distance = calcInter(meanSigs[child], signature[counter]);
            size_t dist_union = calcUnion(meanSigs[child], signature[counter]);

            if (nodeBits >= HD_threshold * 3)
            {
                distance = 0;
            }

            if (distance >= maxInter)
            {
                maxInter = distance;
            }

            // bias for next window match
            if (nodeBits < HD_threshold * 3 && isBranchNode[child] && counter < signature.size() - 1)
            {
                for (size_t grandchild : childLinks[child])
                {
                    if (calcInter(meanSigs[grandchild], signature[counter + 1]) >= HD_threshold)
                    {
                        next_match[i] = 1;
                        // fprintf(stderr, "next match %zu\n", grandchild);
                        break;
                    }
                }
            }
            scores.push_back(distance);

            // printBF(stderr, meanSigs[child]);
            // fprintf(stderr, "--%zu\n", distance);
        }

        // fprintf(stderr, "distance,%zu\n", maxInter);

        for (size_t i = 0; i < childCounts[parent]; i++)
        {
            if (scores[i] == maxInter || next_match[i] == 1)
            {
                size_t child = childLinks[parent][i];

                if (counter == signature.size() - 1)
                {
                    // if (isLeafNode[child])
                    {
                        candidates.push_back(child);
                        max_depths.push_back(counter + 1);
                    }
                    // // perfect match
                    // if (isLeafNode[child])
                    // {
                    //     return child;
                    // }
                    // // perfect match to a shorter seq
                    // else
                    // {
                    //     //?
                    //     // size_t node = 0;
                    //     // for (size_t i = counter; i < signature.size(); i++)
                    //     // {
                    //     //     node = getNewNodeIdx(insertionList);
                    //     //     addSigToMatrix(node, signature[i]);

                    //     //     parentLinks[node] = parent;
                    //     //     childLinks[parent].push_back(node);
                    //     //     isBranchNode[parent] = 1;
                    //     //     childCounts[parent]++;
                    //     //     parent = node;
                    //     // }
                    //     // return node;
                    //     return child;
                    // }
                }

                else
                {
                    traverse(child, signature, counter + 1, candidates, max_depths);
                }
            }
            else if (counter > 0)
            {
                // reached mismatch, only match up to parent;
                candidates.push_back(parent);
                max_depths.push_back(counter - 1);
            }
        }

        // no match
        candidates.push_back(parent);
        max_depths.push_back(counter);

        // return 0;
    }

    // find at which window mismatch starts
    size_t traverse_search(size_t parent, vector<vector<cell_type>> signature, size_t counter, vector<size_t> &candidates, vector<size_t> &max_depths)
    {
        size_t setBits = countSetBits(signature[counter]);
        size_t maxInter = min(HD_threshold, setBits);
        vector<size_t> scores;
        vector<size_t> next_match(childCounts[parent]);
        vector<size_t> maxInter_children;

        // printBF(stderr, signature[counter]);
        // fprintf(stderr, "\n");

        // find max score
        for (size_t i = 0; i < childCounts[parent]; i++)
        {
            size_t child = childLinks[parent][i];
            size_t nodeBits = countSetBits(meanSigs[child]);
            size_t distance = calcInter(meanSigs[child], signature[counter]);

            if (distance >= maxInter)
            {
                maxInter = distance;
            }

            // bias for next window match
            if (nodeBits < HD_threshold * 3 && isBranchNode[child] && counter < signature.size() - 1)
            {
                for (size_t grandchild : childLinks[child])
                {
                    if (calcInter(meanSigs[grandchild], signature[counter + 1]) >= HD_threshold)
                    {
                        next_match[i] = 1;
                        // fprintf(stderr, "next match %zu\n", grandchild);
                        break;
                    }
                }
            }

            scores.push_back(distance);

            // printBF(stderr, meanSigs[child]);
            // fprintf(stderr, "--%zu\n", distance);
        }

        // fprintf(stderr, "distance,%zu\n", maxInter);

        for (size_t i = 0; i < childCounts[parent]; i++)
        {
            if (scores[i] == maxInter || next_match[i] == 1)
            {
                size_t child = childLinks[parent][i];

                if (counter == signature.size() - 1)
                {
                    candidates.push_back(child);
                    max_depths.push_back(counter + 1);
                }

                else
                {
                    traverse(child, signature, counter + 1, candidates, max_depths);
                }
            }
            else if (counter > 0)
            {
                // reached mismatch, only match up to parent;
                candidates.push_back(parent);
                max_depths.push_back(counter - 1);
            }
        }

        // no match
        candidates.push_back(parent);
        max_depths.push_back(counter);
    }

    inline size_t first_insert(vector<vector<cell_type>> signature, vector<size_t> &insertionList, size_t idx, bool store = false)
    {
        size_t node = 0;
        size_t parent = root;
        for (int i = 0; i < signature.size(); i++)
        {
            node = getNewNodeIdx(insertionList);
            addSigToMatrix(node, signature[i]);
            updateMeanSig(node);

            parentLinks[node] = parent;
            childLinks[parent].push_back(node);
            isBranchNode[parent] = 1;
            childCounts[parent]++;
            parent = node;
        }
        isLeafNode[node] = 1;

        if (store)
        {
            dest[node].push_back(idx);
        }
        // dest[node].push_back(idx);

        //?
        overallSigs[node] = signature;
        leafNodes.push_back(node);

        return node;
    }

    inline size_t updateBranch(size_t leaf, vector<vector<cell_type>> signature, vector<size_t> &insertionList, size_t idx, bool store = false)
    {
        size_t node = leaf;
        size_t parent = parentLinks[node];
        size_t w = min(overallSigs[node].size(), signature.size()) - 1;
        while (node != root)
        {

            size_t setBits = countSetBits(signature[w]);
            size_t threshold = min(HD_threshold, setBits);
            size_t distance = calcInter(meanSigs[node], signature[w]);
            if (distance > threshold)
            {
                // matchingWindows++;
                addSigToMatrix(node, signature[w]);
                updateMeanSig(node);
            }
            else
            {
                size_t temp = getNewNodeIdx(insertionList);
                addSigToMatrix(temp, signature[w]);
                updateMeanSig(temp);

                parent = parentLinks[node];

                parentLinks[temp] = parent;
                childLinks[parent].push_back(temp);
                isBranchNode[parent] = 1;
                childCounts[parent]++;
            }
            w--;
        }
        isLeafNode[node] = 1;

        if (store)
        {
            dest[node].push_back(idx);
        }
        // dest[node].push_back(idx);

        //?
        overallSigs[node] = signature;
        leafNodes.push_back(node);

        return node;
    }

    size_t traverse3(vector<vector<cell_type>> signature)
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

    size_t traverse4(vector<vector<cell_type>> signature)
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


    inline size_t insert(vector<vector<cell_type>> signature, vector<size_t> &insertionList, size_t idx, bool store = true)
    {

        size_t node = traverse3(signature);

        if (node == 0)
        {
            // fprintf(stderr, "here\n");
            // node = first_insert(signature, insertionList, idx, store);

            size_t sibling = traverse4(signature);
        }
        else
        {
            // size_t matching = compareSeqsWindows(overallSigs[node], signature);
            // if (matching == signature.size())
            // {
            //     fprintf(stderr, "here %zu\n", idx);
            // }

            // size_t BFsize = overallSigs[node][0].size();

            // for (int w = 0; w < min(overallSigs[node].size(),signature.size()); w++)
            // {
            //     for (int i = 0; i < BFsize; i++)
            //     {
            //         overallSigs[node][w][i] |= signature[w][i];
            //     }
            // }
        }
        // fprintf(stderr, "at %zu\n", node);

        // vector<size_t> matches;
        // vector<size_t> depths;
        // vector<size_t> HDs;
        // size_t ans = traverse2(root, signature, 0, matches, depths, HDs, 0);

        // // auto max = max_element(depths.begin(), depths.end());
        // auto max = max_element(HDs.begin(), HDs.end());
        // auto i = distance(HDs.begin(), max);

        // size_t node = matches[i];
        // size_t parent = node;
        // size_t counter = depths[i];
        // // for (int i = 0; i < matches.size(); i++)
        // // {
        // //     fprintf(stderr, "%zu,%zu,%zu\n", matches[i], depths[i], HDs[i]);
        // // }
        // // fprintf(stderr, "---%zu,%zu,%zu\n", parent, counter, signature.size() - 1);

        // // merge sigs for chosen path
        // size_t temp = node;
        // for (int i = 0; i < counter; i++)
        // {
        //     addSigToMatrix(temp, signature[counter - i - 1]);
        //     updateMeanSig(temp);
        //     temp = parentLinks[temp];
        // }

        // // add new nodes (if any)
        // for (int i = counter; i < signature.size(); i++)
        // {
        //     node = getNewNodeIdx(insertionList);
        //     // fprintf(stderr, "---%zu\n", node);
        //     addSigToMatrix(node, signature[i]);
        //     // matrices[node] = signature[i];
        //     meanSigs[node] = signature[i];

        //     parentLinks[node] = parent;
        //     childLinks[parent].push_back(node);
        //     isBranchNode[parent] = 1;
        //     childCounts[parent]++;
        //     parent = node;
        // }
        // isLeafNode[node] = 1;

        if (node == 0)
        {
            fprintf(stderr, "not found %zu\n", idx);
        }

        if (store)
        {
            dest[node].push_back(idx);
        }
        return node;
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

        size_t node = traverse3(signature);
        // fprintf(stderr, "---%zu,%zu\n", idx, node);

        dest[node].push_back(idx);

        return node;
    }

    size_t search(size_t parent, vector<vector<cell_type>> signature, size_t counter, vector<size_t> &matches)
    {

        // reached the end of seq, no match is found
        if (counter == signature.size())
        {
            fprintf(stderr, "not found\n");
            return 0;
        }

        size_t threshold = min(HD_threshold, countSetBits(signature[counter]));

        for (size_t i = 0; i < childCounts[parent]; i++)
        {

            size_t child = childLinks[parent][i];
            if (calcInter(meanSigs[child], signature[counter]) >= threshold)
            {
                // // this can be partial match
                // if (isLeafNode[child])
                // {
                //     matches.push_back(dest[child][0]);
                // }

                if (counter == signature.size() - 1)
                {
                    // matches.insert(matches.end(), dest[child].begin(), dest[child].end());
                    matches.push_back(child);
                }
                else
                {
                    search(child, signature, counter + 1, matches);
                }
            }
        }

        // all children not match
        return 0;
    }

    void search(vector<vector<cell_type>> signature)
    {
        vector<size_t> matches;
        vector<size_t> depths;
        size_t ans = traverse(root, signature, 0, matches, depths);

        // auto max = max_element(depths.begin(), depths.end());
        // auto i = distance(depths.begin(), max);
        // fprintf(stderr, "%zu\n", matches[i]);

        auto max = *max_element(depths.begin(), depths.end());
        for (int i = 0; i < matches.size(); i++)
        {
            if (depths[i] == max)
            {
                fprintf(stderr, "%zu,", matches[i]);
            }
        }
        fprintf(stderr, "\n");
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