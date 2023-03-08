// #include "bf_ktree.hpp"
#include <fstream>
#include <random>
#include "similarity_main_cmdline.hpp"
#include "bloom_filter.hpp"
#include "read.hpp"
#include "self_tree.hpp"

using namespace std;

static similarity_main_cmdline args; // Command line switches and arguments
// static size_t signatureSize;         // Signature size (depends on element in BF, obtained while read binary)
bool skip = false;
size_t batch = 300;

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

void calcAllSimilarity(FILE *pFile, vector<seq_type> seqs)
{
    size_t seqCount = seqs.size();
    fprintf(pFile, "i,j,similarity\n");
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

            for (size_t j = i + 1; j < seqCount; j++)
            {
                // double dist = calcDistance(seqs[i], seqs[j]) * 100.0 / max(countBits(seqs[i]), countBits(seqs[j]));
                double dist = calcJaccardLocal(seqs[i], seqs[j]) * 100.0;
                fprintf(pFile, "%zu,%zu,%.2f\n", i, j, dist);
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
                fprintf(pFile, "%zu,%zu,%.2f\n", i, j, dist);
            }
        }
    }
}

void calcAllSimilarityLocal(FILE *pFile, vector<seq_type> seqs)
{
    size_t seqCount = seqs.size();
    fprintf(pFile, "i,j,similarity\n");
    if (skip)
    {

        for (size_t i = 0; i < seqCount; i++)
        {
            for (size_t j = i + 1; j < seqCount; j++)
            {
                double dist = calcMatchingWindows(seqs[i], seqs[j]) * 100.0;
                fprintf(pFile, "%zu,%zu,%.2f\n", i, j, dist);
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
                fprintf(pFile, "%zu,%zu,%.2f\n", i, j, dist);
            }
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

void calcAllSimilarityGlobal(FILE *pFile, vector<seq_type> seqs, size_t offset = 0)
{
    size_t seqCount = seqs.size();
    // fprintf(pFile, "i,j,similarity\n");
    if (skip)
    {

        for (size_t i = 0; i < seqCount; i++)
        {
            for (size_t j = i + 1; j < seqCount; j++)
            {
                double dist = calcJaccardGlobal(seqs[i], seqs[j]) * 100.0;
                fprintf(pFile, "%zu,%zu,%.2f\n", i + offset, j + offset, dist);
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
void calcAllSimilarityKmers(FILE *pFile, vector<cell_type> seqs)
{
    size_t seqCount = seqs.size() / signatureSize;
    fprintf(pFile, "i,j,similarity\n");
    if (skip)
    {
        for (size_t i = 0; i < seqCount; i++)
        {
            size_t temp = countSetBits(&seqs[i * signatureSize], signatureSize);
            for (size_t j = i; j < seqCount; j++)
            {
                // size_t bits = max(temp, countSetBits(&seqs[j * signatureSize], signatureSize));
                // fprintf(stdout, "%zu,%zu,%zu\n", i, j, bits);
                double sim = calcSimilarity(&seqs[i * signatureSize], &seqs[j * signatureSize], signatureSize);
                fprintf(pFile, "%zu,%zu,%.2f\n", i, j, sim * 100);
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
                fprintf(pFile, "%zu,%zu,%.2f\n", i, j, sim * 100);
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

// vector<vector<vector<cell_type>>> readPartitionBF(const string file_path, size_t size, size_t &idx)
// {
//     ifstream rf(file_path, ios::out | ios::binary);
//     if (!rf.is_open())
//     {
//         fprintf(stderr, "Invalid File. Please try again\n");
//         exit(0);
//     }

//     unsigned long long int length;
//     if (rf)
//         rf.read(reinterpret_cast<char *>(&length), sizeof(unsigned long long int));

//     vector<vector<vector<cell_type>>> seqs;
//     vector<vector<cell_type>> tseq; // list of BFs for a seq

//     //? 1 window
//     // cout << "length: " << length << "\n";
//     vector<cell_type> bf(length);
//     cell_type temp = 0;
//     size_t i = 0;
//     size_t count = 0;
//     size_t c = sizeof(unsigned long long int);

//     while (rf)
//     {
//         rf.read((char *)&bf[i], sizeof(cell_type));
//         i++;
//         if (i == length)
//         {
//             if (isEmpty(bf))
//             {
//                 seqs.push_back(tseq);
//                 tseq.clear();
//                 count++;
//                 if (count == size)
//                 {
//                     break;
//                 }
//             }
//             else
//             {
//                 tseq.push_back(bf);
//             }
//             fill(bf.begin(), bf.end(), 0);
//             c += i;
//             i = 0;
//         }
//     }
//     c += i;
//     rf.close();
//     idx = c;
//     fprintf(stderr, "%zu,%zu\n", c, length);
//     return seqs;
// }

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

    while (rf && i < idx)
    {
        // fprintf(stderr, "%zu\n", i);
        rf.read((char *)&t, sizeof(cell_type));
        i++;
    }

    i = 0;

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

unsigned long long int readSignatures(const string file, vector<cell_type> &sigs, size_t size)
{
    ifstream rf(file, ios::out | ios::binary);
    if (!rf.is_open())
    {
        fprintf(stderr, "Invalid File. Please try again\n");
        exit(0);
    }

    unsigned long long int length;
    if (rf)
        rf.read(reinterpret_cast<char *>(&length), sizeof(unsigned long long int));

    size_t i = 0;
    size = size * length;
    sigs.resize(size);

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
    rf.close();

    return length;
}

int similarity_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    std::ios::sync_with_stdio(false); // No sync with stdio -> faster

    string bfIn = "data.bin";
    size_t bf_element_cnt = 1000;

    if (args.bf_input_given)
        bfIn = args.bf_input_arg;

    if (args.threshold_given)
        minimiser_match_threshold = args.threshold_arg;

    if (args.skip_arg)
    {
        skip = true;
    }

    size_t firstindex = bfIn.find_last_of("/") + 1;
    size_t lastindex = bfIn.find_last_of(".");
    string rawname = bfIn.substr(firstindex, lastindex - firstindex);

    if (args.output_given)
    {
        rawname = rawname + "-" + args.output_arg;
    }
    FILE *pFile;

    if (args.all_kmer_arg)
    {
        pFile = fopen((rawname + "-all_sim.txt").c_str(), "w");
        vector<cell_type> seqs;
        // signatureSize = readSignatures(bfIn, seqs, batch);
        signatureSize = readSignatures(bfIn, seqs);

        fprintf(stderr, "Loaded %zu seqs...\n", seqs.size() / signatureSize);
        calcAllSimilarityKmers(pFile, seqs);
    }
    else
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

        // fprintf(stderr,"done\n" );
        fprintf(stderr, "Loaded %zu seqs...\n", seqs.size());

        // calcAllSetBits(sigs);
        if (args.local_arg)
        {
            char buffer[50];
            sprintf(buffer, "-t%zu-local_sim.txt", minimiser_match_threshold);
            pFile = fopen((rawname + buffer).c_str(), "w");
            calcAllSimilarityLocal(pFile, seqs);
        }
        else if (args.global_arg)
        {
            pFile = fopen((rawname + "-global_sim.txt").c_str(), "w");
            fprintf(pFile, "i,j,similarity\n");
            calcAllSimilarityGlobal(pFile, seqs);
            size_t start = seqs.size();
            vector<seq_type> temp;
            vector<seq_type> seqsA;
            size_t startA, idxA;
            bool changed = false;

            do
            {
                temp = readPartitionBF(bfIn, batch, idx);
                calcAllSimilarityGlobal(pFile, temp, start);
                calcAllSimilarityGlobalBatch(pFile, seqs, temp, start - seqs.size(), start);
                start += temp.size();

                if (!changed)
                {
                    seqsA = temp;
                    startA = start;
                    idxA = idx;
                    changed = true;
                }
            } while (temp.size() > 0);

            while (seqsA.size() > 0)
            {
                seqs = seqsA;
                start = startA;
                idx = idxA;

                changed = false;
                do
                {
                    temp = readPartitionBF(bfIn, batch, idx);
                    calcAllSimilarityGlobalBatch(pFile, seqs, temp, start - seqs.size(), start);
                    start += temp.size();

                    if (!changed)
                    {
                        seqsA = temp;
                        startA = start;
                        idxA = idx;
                        changed = true;
                    }
                } while (temp.size() > 0);
            }
            fprintf(stderr, "Loaded %zu seqs..\n", start);
        }
        else
        {
            pFile = fopen((rawname + "_sim.txt").c_str(), "w");
            calcAllSimilarity(pFile, seqs);
        }
    }

    // for (int i = 0; i < sigs.size(); i += signatureSize)
    // {
    //     printBF(&sigs[i]);
    // }

    // for (auto sig:sigs){
    //     fprintf(stdout,"%zu\n",sig);
    // }

    return 0;
}
