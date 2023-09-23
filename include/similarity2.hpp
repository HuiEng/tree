#ifndef INCLUDE_SIMILARITY2_HPP
#define INCLUDE_SIMILARITY2_HPP

#include "bloom_filter.hpp"
#include "read.hpp"
#include "distance.hpp"
#include <functional>
using namespace std;

vector<vector<size_t>> compressClusterList(const string file)
{
    string line;
    ifstream clusterStream(file);

    size_t id = 0;
    size_t clus = 0;
    vector<size_t> clusters;

    while (getline(clusterStream, line))
    {
        sscanf(line.c_str(), "%zu,%zu", &id, &clus);
        clusters.push_back(clus);
    }
    clusterStream.close();

    unordered_map<size_t, size_t> remap;
    vector<vector<size_t>> clusterGroup;
    size_t i = 0;
    for (size_t &clus : clusters)
    {
        if (remap.count(clus))
        {
            clus = remap[clus];
            clusterGroup[clus].push_back(i);
        }
        else
        {
            size_t newClus = remap.size();
            remap[clus] = newClus;
            clus = newClus;
            clusterGroup.push_back({i});
        }
        i++;
    }
    fprintf(stderr, "Output %zu clusters\n", remap.size());
    return clusterGroup;
}

template <typename sigType>
vector<vector<sigType>> getSigsGroup(vector<vector<size_t>> clusterGroup,
                                     vector<sigType> &seqs, size_t mul)
{
    vector<vector<sigType>> sigsGroup;
    for (vector<size_t> cluster : clusterGroup)
    {

        vector<sigType> temp;
        for (size_t id : cluster)
        {
            if (mul == 1)
            {
                temp.push_back(seqs[id]);
            }
            else
            {
                sigType *signature = &seqs[id * mul];
                temp.insert(temp.end(), signature, signature + signatureSize);
            }
        }
        sigsGroup.push_back(temp);
    }

    return sigsGroup;
}

void calcAllSimilarityWindow1(FILE *pFile,
                              vector<seq_type> seqs_a,
                              vector<size_t> ids_a)
{
    for (size_t i = 0; i < ids_a.size(); i++)
    {
        for (size_t j = i + 1; j < ids_a.size(); j++)
        {
            double similarity = calcMatchingWindows(seqs_a[i], seqs_a[j]);
            fprintf(pFile, "%zu,%zu,%.2f\n", ids_a[i], ids_a[j], similarity * 100);
        }
    }
}

// Jaccard similarity
void calcAllSimilarityKmers1(FILE *pFile,
                             vector<cell_type> seqs_a,
                             vector<size_t> ids_a)
{
    for (size_t i = 0; i < ids_a.size(); i++)
    {
        for (size_t j = i + 1; j < ids_a.size(); j++)
        {
            double similarity = calcSimilarity(&seqs_a[i * signatureSize], &seqs_a[j * signatureSize], signatureSize);

            fprintf(pFile, "%zu,%zu,%.2f\n", ids_a[i], ids_a[j], similarity * 100);
        }
    }
}

template <typename funcType, typename sigType>
void sepCluster_intra(vector<vector<size_t>> clusterGroup,
                      funcType simFunc,
                      vector<vector<sigType>> sigsGroup, size_t mul = 1)
{

    string rawname = "-all_sim.txt";
    if (mul == 1)
    {
        char buffer[50];
        sprintf(buffer, "-t%zu-window_sim.txt", minimiser_match_threshold);
        rawname = buffer;
    }
    FILE *pFile = fopen(("intra" + rawname).c_str(), "w");
    fprintf(pFile, "i,j,similarity\n");

    for (size_t i = 0; i < clusterGroup.size(); i++)
    {
        simFunc(pFile, sigsGroup[i], clusterGroup[i]);
    }

    fprintf(stderr, "done intra_sim\n");
}

void calcAllSimilarityWindow2(FILE *pFile,
                              vector<seq_type> seqs_a,
                              vector<size_t> ids_a,
                              vector<seq_type> seqs_b,
                              vector<size_t> ids_b)
{
    for (size_t i = 0; i < ids_a.size(); i++)
    {
        for (size_t j = 0; j < ids_b.size(); j++)
        {
            double similarity = calcMatchingWindows(seqs_a[i], seqs_b[j]);
            fprintf(pFile, "%zu,%zu,%.2f\n", ids_a[i], ids_b[j], similarity * 100);
        }
    }
}

// Jaccard similarity
void calcAllSimilarityKmers2(FILE *pFile,
                             vector<cell_type> seqs_a,
                             vector<size_t> ids_a,
                             vector<cell_type> seqs_b,
                             vector<size_t> ids_b)
{
    for (size_t i = 0; i < ids_a.size(); i++)
    {
        for (size_t j = 0; j < ids_b.size(); j++)
        {
            double similarity = calcSimilarity(&seqs_a[i * signatureSize], &seqs_b[j * signatureSize], signatureSize);
            fprintf(pFile, "%zu,%zu,%.2f\n", ids_a[i], ids_b[j], similarity * 100);
        }
    }
}

template <typename funcType, typename sigType>
void sepCluster_inter(vector<vector<size_t>> clusterGroup,
                      funcType simFunc,
                      vector<vector<sigType>> sigsGroup, size_t mul = 1)
{
    string rawname = "-all_sim.txt";
    if (mul == 1)
    {
        char buffer[50];
        sprintf(buffer, "-t%zu-window_sim.txt", minimiser_match_threshold);
        rawname = buffer;
    }
    FILE *pFile = fopen(("inter" + rawname).c_str(), "w");
    fprintf(pFile, "i,j,similarity\n");

    for (size_t i = 0; i < clusterGroup.size(); i++)
    {
        vector<sigType> clusterSigs_a = sigsGroup[i];
        vector<size_t> idx_a = clusterGroup[i];
        for (size_t j = i + 1; j < clusterGroup.size(); j++)
        {
            simFunc(pFile, clusterSigs_a, idx_a,
                    sigsGroup[j], clusterGroup[j]);
        }
    }
    fprintf(stderr, "done inter_sim\n");
}

template <typename sigType>
void sepCluster(vector<vector<size_t>> clusterGroup,
                vector<sigType> &seqs)
{
    vector<vector<sigType>> sigsGroup = getSigsGroup(clusterGroup, seqs, 1);
    sepCluster_intra(clusterGroup, &calcAllSimilarityWindow1, sigsGroup);
    sepCluster_inter(clusterGroup, &calcAllSimilarityWindow2, sigsGroup);
}

template <typename sigType>
void sepCluster(vector<vector<size_t>> clusterGroup,
                vector<sigType> &seqs, size_t mul)
{
    vector<vector<sigType>> sigsGroup = getSigsGroup(clusterGroup, seqs, mul);
    sepCluster_intra(clusterGroup, &calcAllSimilarityKmers1, sigsGroup, mul);
    sepCluster_inter(clusterGroup, &calcAllSimilarityKmers2, sigsGroup, mul);
}

#endif