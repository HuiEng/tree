#ifndef INCLUDE_SIMILARITY_HPP
#define INCLUDE_SIMILARITY_HPP

#include "bloom_filter.hpp"
#include "read.hpp"
#include "distance.hpp"
#include <functional>
using namespace std;

// static size_t signatureSize;         // Signature size (depends on element in BF, obtained while read binary)
bool skip = false;
size_t batch = 300;
double min_print_sim = 0.5;

// void toBinary(cell_type letter)
// {
//     int binary[8];
//     for (int n = 0; n < 8; n++)
//         binary[7 - n] = (letter >> n) & 1;

//     for (int n = 0; n < 8; n++)
//         fprintf(stdout, "%d", binary[n]);
// }

// void printBF(const cell_type *bf)
// {
//     for (std::size_t i = 0; i < signatureSize; ++i)
//     {
//         toBinary(bf[i * sizeof(cell_type)]);
//     }
//     fprintf(stdout, "\n");
// }

void calcAllSimilarityLocal(FILE *pFile, vector<seq_type> seqs, size_t offset = 0)
{
    size_t seqCount = seqs.size();
    if (skip)
    {
        for (size_t i = 0; i < seqCount; i++)
        {
            // // debug
            // fprintf(stdout, "\n%zu,%zu", i, seqs[i].size());
            // for (int w = 0; w < min(seqs[i].size(),seqs[37].size()); w++)
            // {
            //     size_t c = 0;
            //     for (size_t n = 0; n < signatureSize; n++)
            //     {
            //         c += __builtin_popcountll(seqs[i][w][n] & seqs[37][w][n]);
            //     }
            //     fprintf(stdout, ",%zu", c);
            // }
            fprintf(pFile, "%zu,%zu,%.2f\n", i + offset, i + offset, 100.0);
            for (size_t j = i + 1; j < seqCount; j++)
            {
                // double dist = calcDistance(seqs[i], seqs[j]) * 100.0 / max(countBits(seqs[i]), countBits(seqs[j]));
                double similarity = calcJaccardLocal(seqs[i], seqs[j]);
                if (similarity > min_print_sim)
                {
                    fprintf(pFile, "%zu,%zu,%.2f\n", i + offset, j + offset, similarity * 100);
                }
            }
        }
    }
    else
    {
        for (size_t i = 0; i < seqCount; i++)
        {

            for (size_t j = i; j < seqCount; j++)
            {
                // double dist = calcDistance(seqs[i], seqs[j]) * 100.0 / max(countBits(seqs[i]), countBits(seqs[j]));
                double dist = calcJaccardLocal(seqs[i], seqs[j]) * 100.0;
                fprintf(pFile, "%zu,%zu,%.2f\n", i + offset, j + offset, dist);
            }
        }
    }
}

void calcAllSimilarityWindow(FILE *pFile, vector<seq_type> seqs, size_t offset = 0)
{
    size_t seqCount = seqs.size();
    if (skip)
    {

        for (size_t i = 0; i < seqCount; i++)
        {
            fprintf(pFile, "%zu,%zu,%.2f\n", i + offset, i + offset, 100.0);
            for (size_t j = i + 1; j < seqCount; j++)
            {
                double similarity = calcMatchingWindows(seqs[i], seqs[j]);
                if (similarity > min_print_sim)
                {
                    fprintf(pFile, "%zu,%zu,%.2f\n", i + offset, j + offset, similarity * 100);
                }
            }
        }
    }
    else
    {
        for (size_t i = 0; i < seqCount; i++)
        {
            for (size_t j = 0; j < seqCount; j++)
            {
                double dist = calcMatchingWindows(seqs[i], seqs[j]) * 100.0;
                fprintf(pFile, "%zu,%zu,%.2f\n", i + offset, j + offset, dist);
            }
        }
    }
}

void calcAllSimilarityGlobal(FILE *pFile, vector<seq_type> seqs, size_t offset = 0)
{
    size_t seqCount = seqs.size();
    // fprintf(pFile, "i,j,similarity\n");
    if (skip)
    {

        for (size_t i = 0; i < seqCount; i++)
        {
            fprintf(pFile, "%zu,%zu,%.2f\n", i + offset, i + offset, 100.0);
            for (size_t j = i + 1; j < seqCount; j++)
            {
                double similarity = calcJaccardGlobal(seqs[i], seqs[j]);
                if (similarity > min_print_sim)
                {
                    fprintf(pFile, "%zu,%zu,%.2f\n", i + offset, j + offset, similarity * 100);
                }
            }
        }
    }
    else
    {
        for (size_t i = 0; i < seqCount; i++)
        {
            for (size_t j = 0; j < seqCount; j++)
            {
                double dist = calcJaccardGlobal(seqs[i], seqs[j]) * 100.0;
                fprintf(pFile, "%zu,%zu,%.2f\n", i + offset, j + offset, dist);
            }
        }
    }
}

// Jaccard similarity
void calcAllSimilarityKmers(FILE *pFile, vector<cell_type> seqs, size_t offset = 0)
{
    size_t seqCount = seqs.size() / signatureSize;
    if (skip)
    {
        for (size_t i = 0; i < seqCount; i++)
        {
            // size_t temp = countSetBits(&seqs[i * signatureSize], signatureSize);
            fprintf(pFile, "%zu,%zu,%.2f\n", i + offset, i + offset, 100.0);

            for (size_t j = i + 1; j < seqCount; j++)
            {
                // size_t bits = max(temp, countSetBits(&seqs[j * signatureSize], signatureSize));
                // fprintf(stdout, "%zu,%zu,%zu\n", i, j, bits);
                double sim = calcSimilarity(&seqs[i * signatureSize], &seqs[j * signatureSize], signatureSize);
                if (sim > min_print_sim)
                {
                    fprintf(pFile, "%zu,%zu,%.2f\n", i + offset, j + offset, sim * 100);
                }
            }
        }
    }
    else
    {
        for (size_t i = 0; i < seqCount; i++)
        {
            size_t temp = countSetBits(&seqs[i * signatureSize], signatureSize);
            for (size_t j = 0; j < seqCount; j++)
            {
                double sim = calcSimilarity(&seqs[i * signatureSize], &seqs[j * signatureSize], signatureSize);
                fprintf(pFile, "%zu,%zu,%.2f\n", i + offset, j + offset, sim * 100);
            }
        }
    }
}

void calcAllSetBits(const vector<cell_type> &sigs)
{
    size_t seqCount = sigs.size() / signatureSize;
    for (size_t i = 0; i < seqCount; i++)
    {
        // fprintf(stderr,"%zu,%zu,%zu\n",kmerLength,windowLength,countSetBits(&sigs[i*signatureSize]));
        fprintf(stderr, "BF density for the first seq: %zu\n", countSetBits(&sigs[i * signatureSize], signatureSize));
        break;
    }
}

vector<vector<vector<cell_type>>> readPartitionBF(const string file_path, size_t size, size_t &idx)
{
    ifstream rf(file_path, ios::out | ios::binary);
    if (!rf.is_open())
    {
        fprintf(stderr, "Invalid File. Please try again\n");
        exit(0);
    }

    unsigned long long int length;
    if (rf)
        rf.read(reinterpret_cast<char *>(&length), sizeof(unsigned long long int));

    vector<vector<vector<cell_type>>> seqs;
    vector<vector<cell_type>> tseq; // list of BFs for a seq

    //? 1 window
    // cout << "length: " << length << "\n";
    vector<cell_type> bf(length);
    cell_type temp = 0;
    size_t i = 0;
    size_t count = 0;
    size_t c = 0;
    cell_type t;

    // while (rf && i < idx)
    // {
    //     // fprintf(stderr, "%zu\n", i);
    //     rf.read((char *)&t, sizeof(cell_type));
    //     i++;
    // }

    // i = 0;

    rf.seekg(idx + sizeof(unsigned long long int), rf.beg);
    rf.read((char *)&bf[i], sizeof(cell_type));
    rf.seekg(idx + sizeof(unsigned long long int), rf.beg);

    while (rf)
    {
        rf.read((char *)&bf[i], sizeof(cell_type));
        i++;
        c++;
        if (i == length)
        {
            if (isEmpty(bf))
            {
                seqs.push_back(tseq);
                tseq.clear();
                count++;
                if (count == size)
                {
                    break;
                }
            }
            else
            {
                tseq.push_back(bf);
            }
            fill(bf.begin(), bf.end(), 0);
            i = 0;
        }
    }
    rf.close();
    idx += c;
    return seqs;
}

unsigned long long int readSignatures(const string file, vector<cell_type> &sigs, size_t size, size_t &idx)
{
    ifstream rf(file, ios::out | ios::binary);
    if (!rf.is_open())
    {
        fprintf(stderr, "Invalid File. Please try again\n");
        exit(0);
    }

    // get length of file:
    rf.seekg(0, rf.end);
    unsigned long long int len = rf.tellg();
    len = len - sizeof(unsigned long long int); // / sizeof(uint64_t);
    rf.seekg(0, rf.beg);

    unsigned long long int length;
    if (rf)
        rf.read(reinterpret_cast<char *>(&length), sizeof(unsigned long long int));

    size = size * length;
    sigs.resize(size);
    size_t i = 0;
    cell_type t;

    // while (rf && i < idx)
    // {
    //     // fprintf(stderr, "%zu\n", i);
    //     rf.read((char *)&t, sizeof(cell_type));
    //     i++;
    // }
    // i = 0;

    rf.seekg(idx + sizeof(unsigned long long int), rf.beg);
    rf.read((char *)&sigs[i], sizeof(cell_type));
    // if (idx >= len)
    // {
    //     sigs.clear();
    //     rf.close();
    //     return length;
    // }
    rf.seekg(idx + sizeof(unsigned long long int), rf.beg);

    while (rf)
    {
        rf.read((char *)&sigs[i], sizeof(cell_type));
        i++;

        if (i == size)
        {
            break;
        }
    }

    sigs.resize(i);
    idx += i;

    rf.close();

    return length;
}

void calcAllSimilarityBatch(FILE *pFile, vector<seq_type> seqsA, vector<seq_type> seqsB, size_t offsetA, size_t offsetB)
{
    for (size_t i = 0; i < seqsA.size(); i++)
    {
        for (size_t j = 0; j < seqsB.size(); j++)
        {
            double dist = calcJaccardLocal(seqsA[i], seqsB[j]) * 100.0;
            fprintf(pFile, "%zu,%zu,%.2f\n", i + offsetA, j + offsetB, dist);
        }
    }
}

void calcAllSimilarityLocalBatch(FILE *pFile, vector<seq_type> seqsA, vector<seq_type> seqsB, size_t offsetA, size_t offsetB)
{
    for (size_t i = 0; i < seqsA.size(); i++)
    {
        for (size_t j = 0; j < seqsB.size(); j++)
        {
            double dist = calcMatchingWindows(seqsA[i], seqsB[j]) * 100.0;
            fprintf(pFile, "%zu,%zu,%.2f\n", i + offsetA, j + offsetB, dist);
        }
    }
}

void calcAllSimilarityGlobalBatch(FILE *pFile, vector<seq_type> seqsA, vector<seq_type> seqsB, size_t offsetA, size_t offsetB)
{
    for (size_t i = 0; i < seqsA.size(); i++)
    {
        for (size_t j = 0; j < seqsB.size(); j++)
        {
            double dist = calcJaccardGlobal(seqsA[i], seqsB[j]) * 100.0;
            fprintf(pFile, "%zu,%zu,%.2f\n", i + offsetA, j + offsetB, dist);
        }
    }
}

void calcAllSimilarityKmersBatch(FILE *pFile, vector<cell_type> seqsA, vector<cell_type> seqsB, size_t offsetA, size_t offsetB)
{
    size_t seqCountA = seqsA.size() / signatureSize;
    size_t seqCountB = seqsB.size() / signatureSize;
    for (size_t i = 0; i < seqCountA; i++)
    {
        for (size_t j = 0; j < seqCountB; j++)
        {
            double dist = calcSimilarity(&seqsA[i * signatureSize], &seqsB[j * signatureSize], signatureSize) * 100.0;
            fprintf(pFile, "%zu,%zu,%.2f\n", i + offsetA, j + offsetB, dist);
        }
    }
}

template <typename seq_batch_type, typename t1, typename t2>
void batchSimFunction(FILE *pFile, vector<seq_batch_type> seqs_batch, t1 simFunc, t2 simBatchFunc)
{
    // size_t temp = seqs_batch.size() - 1;
    // size_t seqCount = temp * batch + seqs_batch[temp].size() / signatureSize;
    // temp++;
    // fprintf(stderr, "Loaded %zu seqs...\n", seqCount);
    skip = true;
    size_t seqCount = seqs_batch.size();

    // FILE *pFile = fopen("test-all_sim.txt", "w");
    fprintf(pFile, "i,j,similarity\n");
    for (size_t i = 0; i < seqCount; i++)
    {
        seq_batch_type seqsA = seqs_batch[i];
        size_t offset = i * batch;
        simFunc(pFile, seqsA, offset);
        for (size_t j = i + 1; j < seqCount; j++)
        {
            seq_batch_type seqsB = seqs_batch[j];
            simBatchFunc(pFile, seqsA, seqsB, offset, j * batch);
        }
    }
}

int batching(string bfIn, FILE *pFile,
             function<void(FILE *, vector<seq_type>, size_t)> simFunc,
             function<void(FILE *, vector<seq_type>, vector<seq_type>, size_t, size_t)> simBatchFunc)
{
    size_t idx = 0;
    vector<seq_type> seqs = readPartitionBF(bfIn, batch, idx);
    // vector<seq_type> seqs = readPartitionBF(bfIn);
    // dummy code, assume there is at least 10 input seqs
    for (int i = 0; i < 10; i++)
    {
        if (seqs[i].size() > 0)
        {
            signatureSize = seqs[i][0].size();
            break;
        }
    }
    if (signatureSize == 0)
    {
        fprintf(stderr, "Something is wrong with the input data, please generate signature with diff params\n");
        return 0;
    }

    fprintf(pFile, "i,j,similarity\n");
    // calcAllSimilarityGlobal(pFile, seqs);
    simFunc(pFile, seqs, 0);

    size_t start = seqs.size();
    vector<seq_type> temp;
    vector<seq_type> seqsA;
    size_t startA = 0;
    size_t idxA = 0;
    bool changed = false;

    do
    {
        temp = readPartitionBF(bfIn, batch, idx);
        simFunc(pFile, temp, start);
        simBatchFunc(pFile, seqs, temp, startA, start);

        if (!changed)
        {
            seqsA = temp;
            idxA = idx;
            changed = true;
        }
        start += temp.size();
    } while (temp.size() > 0);

    while (seqsA.size() > 1)
    {
        startA += seqs.size();
        start = startA + seqsA.size();
        seqs = seqsA;
        idx = idxA;

        changed = false;
        do
        {
            temp = readPartitionBF(bfIn, batch, idx);
            simBatchFunc(pFile, seqs, temp, startA, start);
            if (!changed)
            {
                seqsA = temp;
                idxA = idx;
                changed = true;
            }
            start += temp.size();
        } while (temp.size() > 0);
    }
    fprintf(stderr, "Loaded %zu seqs..\n", start);
    return 1;
}

#endif