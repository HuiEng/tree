// #include "bf_ktree.hpp"
#include <fstream>
#include <random>
#include "stats_main_cmdline.hpp"
#include "stats.hpp"

#include <regex>
#include <vector>
#include <numeric>
#include <iostream>
// #include "similarity.hpp"

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

    if (args.runs_given)
    {
        runs = args.runs_arg;
        max_seqCount = max_seqCount * runs;
    }

    cmatch matches;
    if (regex_search(args.bf_input_arg, matches, regex("-b([0-9]+)")))
    {
        stringstream sstream(matches[1]);
        sstream >> bf_element_cnt;
    }

    if (args.all_kmer_arg)
    {
        vector<cell_type> seqs;
        signatureSize = readSignaturesSample(bfIn, seqs, max_seqCount);
        fprintf(stderr, "Loaded %zu seqs...\n", seqs.size() / signatureSize);
        max_seqCount = args.max_arg;
        calcAllStatsKmersBatch(seqs);
    }
    else
    {
        vector<seq_type> seqs = readPartitionBFSample(bfIn, signatureSize, max_seqCount);
        fprintf(stderr, "Loaded %zu seqs...\n", seqs.size());
        max_seqCount = args.max_arg;
        if (args.local_arg)
        {
            if (args.threshold_given)
            {

                minimiser_match_threshold = args.threshold_arg;
            }
            fprintf(stderr, "Matching at least %zu minimisers per window...\n", minimiser_match_threshold);
            calcAllStatsBatch(seqs, &calcMatchingWindows);
        }
        else if (args.global_arg)
        {
            calcAllStatsBatch(seqs, &calcJaccardGlobal);
        }
        else
        {
            calcAllStatsBatch(seqs, &calcJaccardLocal);
        }
    }

    return 0;

    // if (args.all_kmer_arg)
    // {
    //     vector<cell_type> seqs;
    //     signatureSize = readSignaturesSample(bfIn, seqs, max_seqCount);
    //     calcAllStatsKmers(seqs);

    //     // signatureSize = readSignatures(bfIn, seqs);

    //     // seqs.clear();

    //     // if (seqs.size() == 0)
    //     // {
    //     //     vector<vector<cell_type>> seqs_batch = readSignaturesBatch(bfIn, batch_size, signatureSize);
    //     //     // batchSimFunction(seqs_batch, &calcAllSimilarityKmers, &calcAllSimilarityKmersBatch);
    //     // }
    //     // else
    //     // {
    //     //     fprintf(stderr, "Loaded %zu seqs...\n", seqs.size() / signatureSize);
    //     //     calcAllStatsKmers(seqs);
    //     // }
    // }
    // else
    // {

    //     vector<seq_type> seqs = readPartitionBFSample(bfIn, signatureSize, max_seqCount);

    //     if (args.local_arg)
    //     {
    //         if (args.threshold_given)
    //         {

    //             minimiser_match_threshold = args.threshold_arg;
    //         }
    //         fprintf(stderr, "Matching at least %zu minimisers per window...\n", minimiser_match_threshold);
    //         calcAllStatsLocal(seqs);
    //     }
    //     else if (args.global_arg)
    //     {
    //         calcAllStatsGlobal(seqs);
    //     }
    //     else
    //     {
    //         calcAllStats(seqs);
    //     }

    //     // vector<seq_type> seqs = readPartitionBF(bfIn, signatureSize);
    //     // if (signatureSize == 0)
    //     // {
    //     //     fprintf(stderr, "Something is wrong with the input data, please generate signature with diff params\n");
    //     //     return 0;
    //     // }

    //     // // fprintf(stderr, "Loaded %zu seqs...\n", seqs.size());

    //     // if (seqs.size() == 0)
    //     // {
    //     //     fprintf(stderr, "here\n");
    //     //     vector<vector<seq_type>> seqs_batch = readPartitionBFBatch(bfIn, batch_size, signatureSize);
    //     //     size_t temp = seqs_batch.size() - 1;
    //     //     size_t seqCount = temp * batch_size + seqs_batch[temp].size();
    //     //     fprintf(stderr, " Batch Loaded %zu seqs...\n", seqCount);

    //     //     if (args.local_arg)
    //     //     {
    //     //         if (args.threshold_given)
    //     //         {

    //     //             minimiser_match_threshold = args.threshold_arg;
    //     //         }
    //     //         fprintf(stderr, "Matching at least %zu minimisers per window...\n", minimiser_match_threshold);
    //     //         calcAllStatsBatch(seqs_batch, seqCount, &calcMatchingWindows);
    //     //     }
    //     //     else if (args.global_arg)
    //     //     {
    //     //         calcAllStatsBatch(seqs_batch, seqCount, &calcJaccardGlobal);
    //     //     }
    //     //     else
    //     //     {
    //     //         calcAllStatsBatch(seqs_batch, seqCount, &calcJaccardLocal);
    //     //     }
    //     // }
    //     // else
    //     // {

    //     //     if (args.local_arg)
    //     //     {
    //     //         if (args.threshold_given)
    //     //         {

    //     //             minimiser_match_threshold = args.threshold_arg;
    //     //         }
    //     //         fprintf(stderr, "Matching at least %zu minimisers per window...\n", minimiser_match_threshold);
    //     //         calcAllStatsLocal(seqs);
    //     //     }
    //     //     else if (args.global_arg)
    //     //     {
    //     //         calcAllStatsGlobal(seqs);
    //     //     }
    //     //     else
    //     //     {
    //     //         calcAllStats(seqs);
    //     //     }
    //     // }
    // }

    return 0;
}
