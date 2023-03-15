// #include "bf_ktree.hpp"
#include <fstream>
#include <random>
#include "stats_main_cmdline.hpp"
#include "stats.hpp"

#include <regex>
#include <vector>
#include <numeric>
#include <iostream>

using namespace std;

static stats_main_cmdline args; // Command line switches and arguments
size_t batch_size = 5;
// Jaccard similarity
void calcAllSimilarityKmers2(FILE *pFile, vector<cell_type> seqs, size_t offset = 0)
{
    size_t seqCount = seqs.size() / signatureSize;

    for (size_t i = 0; i < seqCount; i++)
    {
        // size_t temp = countSetBits(&seqs[i * signatureSize], signatureSize);
        fprintf(pFile, "%zu,%zu,%.2f\n", i + offset, i + offset, 100.0);

        for (size_t j = i + 1; j < seqCount; j++)
        {
            // size_t bits = max(temp, countSetBits(&seqs[j * signatureSize], signatureSize));
            // fprintf(stdout, "%zu,%zu,%zu\n", i, j, bits);
            double sim = calcSimilarity(&seqs[i * signatureSize], &seqs[j * signatureSize], signatureSize);
            fprintf(pFile, "%zu,%zu,%.2f\n", i + offset, j + offset, sim * 100);
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

    if (args.all_kmer_arg)
    {
        vector<cell_type> seqs;
        signatureSize = readSignatures(bfIn, seqs);

        seqs.clear();

        if (seqs.size() == 0)
        {
            fprintf(stderr, "here\n");
            vector<vector<cell_type>> seqs_batch = readSignaturesBatch(bfIn, batch_size);

            size_t temp = seqs_batch.size() - 1;
            size_t seqCount = temp * batch_size + seqs_batch[temp].size() / signatureSize;
            temp++;
            fprintf(stderr, "Loaded %zu seqs...\n", seqCount);

            FILE *pFile = fopen("test-all_sim.txt", "w");
            fprintf(pFile, "i,j,similarity\n");
            for (size_t i = 0; i < temp; i++)
            {
                vector<cell_type> seqsA = seqs_batch[i];
                size_t offset = i * batch_size;
                calcAllSimilarityKmers2(pFile, seqsA, offset);
                for (size_t j = i + 1; j < temp; j++)
                {
                    vector<cell_type> seqsB = seqs_batch[j];
                    fprintf(stderr, "%zu,%zu\n", offset, j * batch_size);
                    calcAllSimilarityKmersBatch(pFile, seqsA, seqsB, offset, j * batch_size);
                }
            }
        }
        else
        {
            fprintf(stderr, "Loaded %zu seqs...\n", seqs.size() / signatureSize);

            calcAllStatsKmers(seqs);
        }
    }
    else
    {
        vector<seq_type> seqs = readPartitionBF(bfIn);
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

    return 0;
}
