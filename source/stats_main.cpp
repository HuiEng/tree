// #include "bf_ktree.hpp"
#include <fstream>
#include <random>
#include "stats_main_cmdline.hpp"
#include "stats.hpp"

#include <regex>
#include <vector>
#include <numeric>
#include <iostream>
#include "similarity.hpp"

using namespace std;

static stats_main_cmdline args; // Command line switches and arguments

template <typename seq_batch_type, typename t1, typename t2>
void batchFunction(vector<seq_batch_type> seqs_batch, t1 simFunc, t2 simBatchFunc)
{
    size_t temp = seqs_batch.size() - 1;
    size_t seqCount = temp * batch_size + seqs_batch[temp].size() / signatureSize;
    temp++;
    fprintf(stderr, "Loaded %zu seqs...\n", seqCount);

    FILE *pFile = fopen("test-all_sim.txt", "w");
    fprintf(pFile, "i,j,similarity\n");
    for (size_t i = 0; i < temp; i++)
    {
        seq_batch_type seqsA = seqs_batch[i];
        size_t offset = i * batch_size;
        simFunc(pFile, seqsA, offset);
        for (size_t j = i + 1; j < temp; j++)
        {
            seq_batch_type seqsB = seqs_batch[j];
            simBatchFunc(pFile, seqsA, seqsB, offset, j * batch_size);
        }
    }
}

size_t estimateSeqCount(ifstream &rf)
{
    // get length of file:
    rf.seekg(0, rf.end);
    unsigned long long int len = rf.tellg();
    len = len - sizeof(unsigned long long int); // / sizeof(uint64_t);
    rf.seekg(0, rf.beg);

    unsigned long long int length;
    if (rf)
        rf.read(reinterpret_cast<char *>(&length), sizeof(unsigned long long int));

    size_t winNum = 0;
    vector<cell_type> bf(length);
    size_t i = 0;
    while (rf)
    {
        rf.read((char *)&bf[i], sizeof(cell_type));
        i++;
        if (i == length)
        {
            if (isEmpty(bf))
            {
                rf.seekg(0, rf.beg);
                return (len * 1.0 / length) / winNum;
            }
            else
            {
                winNum++;
            }
            fill(bf.begin(), bf.end(), 0);
            i = 0;
        }
    }
}

vector<vector<vector<cell_type>>> readPartitionBFSample(const string file_path, size_t &signatureSize, size_t sampleSize)
{
    ifstream rf(file_path, ios::out | ios::binary);
    if (!rf.is_open())
    {
        fprintf(stderr, "Invalid File. Please try again\n");
        exit(0);
    }

    size_t seqCount = estimateSeqCount(rf);

    unsigned long long int length;
    if (rf)
        rf.read(reinterpret_cast<char *>(&length), sizeof(unsigned long long int));
    fprintf(stderr, "File contains approximately %zu seqs...\n", seqCount);
    signatureSize = length;

    //? 1 window
    // cout << "length: " << length << "\n";
    vector<cell_type> bf(length);
    cell_type temp = 0;
    size_t i = 0;

    vector<vector<vector<cell_type>>> seqs;
    vector<vector<cell_type>> tseq; // list of BFs for a seq

    if (sampleSize > seqCount)
    {
        while (rf)
        {
            rf.read((char *)&bf[i], sizeof(cell_type));
            i++;
            if (i == length)
            {
                if (isEmpty(bf))
                {
                    seqs.push_back(tseq);
                    tseq.clear();
                }
                else
                {
                    tseq.push_back(bf);
                }
                fill(bf.begin(), bf.end(), 0);
                i = 0;
            }
        }
    }
    else
    {
        // actual sample size might be smaller due to the estimated seqCount
        // just ignore for now
        vector<size_t> indices = getIndices(seqCount, sampleSize);
        // for (size_t i : indices)
        // {
        //     fprintf(stderr, "%zu,", i);
        // }
        // fprintf(stderr, "\n");

        size_t idx = indices.back();
        indices.pop_back();
        size_t c = 0;
        size_t last_seen = 0;
        size_t n = sizeof(unsigned long long int);

        while (rf)
        {
            rf.read((char *)&bf[i], sizeof(cell_type));
            i++;
            n++;
            if (i == length)
            {
                if (isEmpty(bf))
                {
                    if (c == idx)
                    {
                        fprintf(stderr, "%zu,", idx);
                        seqs.push_back(tseq);
                        idx = indices.back();
                        indices.pop_back();
                        last_seen = n;

                        if (seqs.size() == sampleSize)
                        {
                            break;
                        }
                    }
                    c++;
                    tseq.clear();
                }
                else
                {
                    tseq.push_back(bf);
                }
                fill(bf.begin(), bf.end(), 0);
                i = 0;
            }
        }

        if (seqs.size() < sampleSize)
        {
            fill(bf.begin(), bf.end(), 0);
            tseq.clear();
            i = 0;

            fprintf(stderr, "\nloaded %zu\n", seqs.size());
            rf.clear();
            rf.seekg(last_seen, rf.beg);
            while (rf)
            {
                rf.read((char *)&bf[i], sizeof(cell_type));
                i++;
                if (i == length)
                {
                    if (isEmpty(bf))
                    {
                        seqs.push_back(tseq);
                        tseq.clear();

                        if (seqs.size() == sampleSize)
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
        }
    }
    rf.close();

    return seqs;
}

int stats_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    ios::sync_with_stdio(false); // No sync with stdio -> faster

    string bfIn = args.bf_input_arg;
    size_t bf_element_cnt = 1000;

    size_t firstindex = bfIn.find_last_of("/") + 1;
    size_t lastindex = bfIn.find_last_of(".");
    string rawname = bfIn.substr(firstindex, lastindex - firstindex);

    if (args.output_given)
    {
        rawname = rawname + "-" + args.output_arg;
    }

    if (args.max_given)
    {
        max_seqCount = args.max_arg;
    }

    cmatch matches;
    if (regex_search(args.bf_input_arg, matches, regex("-b([0-9]+)")))
    {
        stringstream sstream(matches[1]);
        sstream >> bf_element_cnt;
    }

    vector<seq_type> seqs = readPartitionBFSample(bfIn, signatureSize, 10);
    fprintf(stderr, "Loaded %zu samples...\n", seqs.size());
    FILE *pFile = fopen("test-all_sim.txt", "w");
    fprintf(pFile, "i,j,similarity\n");
    calcAllSimilarityGlobal(pFile, seqs);

    return 0;

    if (args.all_kmer_arg)
    {
        vector<cell_type> seqs;
        signatureSize = readSignaturesSample(bfIn, seqs, max_seqCount);
        calcAllStatsKmers(seqs);

        // signatureSize = readSignatures(bfIn, seqs);

        // seqs.clear();

        // if (seqs.size() == 0)
        // {
        //     vector<vector<cell_type>> seqs_batch = readSignaturesBatch(bfIn, batch_size, signatureSize);
        //     // batchSimFunction(seqs_batch, &calcAllSimilarityKmers, &calcAllSimilarityKmersBatch);
        // }
        // else
        // {
        //     fprintf(stderr, "Loaded %zu seqs...\n", seqs.size() / signatureSize);
        //     calcAllStatsKmers(seqs);
        // }
    }
    else
    {
        vector<seq_type> seqs = readPartitionBF(bfIn, signatureSize);
        if (signatureSize == 0)
        {
            fprintf(stderr, "Something is wrong with the input data, please generate signature with diff params\n");
            return 0;
        }

        fprintf(stderr, "done\n");
        // fprintf(stderr, "Loaded %zu seqs...\n", seqs.size());

        if (seqs.size() == 0)
        {
            fprintf(stderr, "here\n");
            vector<vector<seq_type>> seqs_batch = readPartitionBFBatch(bfIn, batch_size, signatureSize);
            size_t temp = seqs_batch.size() - 1;
            size_t seqCount = temp * batch_size + seqs_batch[temp].size();
            fprintf(stderr, " Batch Loaded %zu seqs...\n", seqCount);

            if (args.local_arg)
            {
                if (args.threshold_given)
                {

                    minimiser_match_threshold = args.threshold_arg;
                }
                fprintf(stderr, "Matching at least %zu minimisers per window...\n", minimiser_match_threshold);
                calcAllStatsBatch(seqs_batch, seqCount, &calcMatchingWindows);
            }
            else if (args.global_arg)
            {
                calcAllStatsBatch(seqs_batch, seqCount, &calcJaccardGlobal);
            }
            else
            {
                calcAllStatsBatch(seqs_batch, seqCount, &calcJaccardLocal);
            }
        }
        else
        {

            if (args.local_arg)
            {
                if (args.threshold_given)
                {

                    minimiser_match_threshold = args.threshold_arg;
                }
                fprintf(stderr, "Matching at least %zu minimisers per window...\n", minimiser_match_threshold);
                calcAllStatsLocal(seqs);
            }
            else if (args.global_arg)
            {
                calcAllStatsGlobal(seqs);
            }
            else
            {
                calcAllStats(seqs);
            }
        }
    }

    return 0;
}
