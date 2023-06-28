// #include "bf_ktree.hpp"
#include <fstream>
#include <random>
#include "similarity_main_cmdline.hpp"
#include "similarity.hpp"
static similarity_main_cmdline args; // Command line switches and arguments

void allKmersBatch(FILE *pFile, const string inputFile)
{
    // pFile = fopen((rawname + "-all_sim.txt").c_str(), "w");
    size_t idx = 0;
    vector<cell_type> seqs;
    signatureSize = readSignatures(inputFile, seqs, batch, idx);

    fprintf(pFile, "i,j,similarity\n");
    calcAllSimilarityKmers(pFile, seqs);

    size_t start = idx;
    vector<cell_type> temp;
    vector<cell_type> seqsA;
    size_t startA = start;
    size_t idxA = 0;
    bool changed = false;

    do
    {
        readSignatures(inputFile, temp, batch, idx);
        size_t offset = (idx - temp.size()) / signatureSize;
        calcAllSimilarityKmers(pFile, temp, offset);
        calcAllSimilarityKmersBatch(pFile, seqs, temp, idxA / signatureSize, offset);

        if (!changed)
        {
            seqsA = temp;
            startA = start;
            changed = true;
        }
        start = idx;

    } while (temp.size() > 0);

    size_t seqCount = idx;
    while (seqsA.size() > 1)
    {
        seqs = seqsA;
        idx = startA + seqsA.size();
        start = idx;
        idxA = startA;

        changed = false;

        do
        {
            readSignatures(inputFile, temp, batch, idx);
            size_t offset = (idx - temp.size()) / signatureSize;
            calcAllSimilarityKmersBatch(pFile, seqs, temp, idxA / signatureSize, offset);

            if (!changed)
            {
                seqsA = temp;
                startA = start;
                // idxA = idx;
                changed = true;
            }
            start = idx;

        } while (temp.size() > 0);
    }
    seqCount = seqCount / signatureSize;
    fprintf(stderr, "Loaded %zu seqs..\n", seqCount);
}

int similarity_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    std::ios::sync_with_stdio(false); // No sync with stdio -> faster

    string inputFile = args.bf_input_arg;


    if (args.threshold_given)
        minimiser_match_threshold = args.threshold_arg;

    if (args.skip_arg)
    {
        skip = true;
    }

    if (args.min_sim_given)
    {
        min_print_sim = args.min_sim_arg;
        skip = true;
    }

    if (skip)
    {
        fprintf(stderr, "Only printing entries > %.1f %%\n", min_print_sim * 100);
    }

    size_t firstindex = inputFile.find_last_of("/") + 1;
    size_t lastindex = inputFile.find_last_of(".");
    string rawname = inputFile.substr(firstindex, lastindex - firstindex);

    if (args.output_given)
    {
        rawname = rawname + "-" + args.output_arg;
    }
    FILE *pFile;

    if (args.batch_given)
    {
        skip = true;
        batch = args.batch_arg;
        if (args.all_kmer_arg)
        {
            pFile = fopen((rawname + "-all_sim.txt").c_str(), "w");
            allKmersBatch(pFile, inputFile);
            // vector<vector<cell_type>> seqs_batch = readSignaturesBatch(inputFile, batch, signatureSize);
            // size_t temp = seqs_batch.size() - 1;
            // size_t seqCount = temp * batch + seqs_batch[temp].size() / signatureSize;
            // fprintf(stderr, "Loaded %zu seqs...\n", seqCount);
            // batchSimFunction(pFile, seqs_batch, &calcAllSimilarityKmers, &calcAllSimilarityKmersBatch);
        }
        else
        {
            // vector<vector<vector<vector<cell_type>>>> seqs_batch = readPartitionBFBatch(inputFile, batch, signatureSize);
            // size_t temp = seqs_batch.size() - 1;
            // size_t seqCount = temp * batch + seqs_batch[temp].size();
            // fprintf(stderr, "Loaded %zu seqs...\n", seqCount);

            if (args.local_arg)
            {
                char buffer[50];
                sprintf(buffer, "-t%zu-local_sim.txt", minimiser_match_threshold);
                pFile = fopen((rawname + buffer).c_str(), "w");
                batching(inputFile, pFile, &calcAllSimilarityWindow, &calcAllSimilarityLocalBatch);
                // batchSimFunction(pFile, seqs_batch, &calcAllSimilarityLocal, &calcAllSimilarityLocalBatch);
            }
            else if (args.global_arg)
            {
                pFile = fopen((rawname + "-global_sim.txt").c_str(), "w");
                batching(inputFile, pFile, &calcAllSimilarityGlobal, &calcAllSimilarityGlobalBatch);
                // batchSimFunction(pFile, seqs_batch, &calcAllSimilarityGlobal, &calcAllSimilarityGlobalBatch);
            }
            else
            {
                pFile = fopen((rawname + "_sim.txt").c_str(), "w");
                batching(inputFile, pFile, &calcAllSimilarityLocal, &calcAllSimilarityBatch);
                // batchSimFunction(pFile, seqs_batch, &calcAllSimilarity, &calcAllSimilarityBatch);
            }
        }
    }
    else
    {
        if (args.all_kmer_arg)
        {
            pFile = fopen((rawname + "-all_sim.txt").c_str(), "w");
            vector<cell_type> seqs;
            if (args.multiple_arg)
            {
                signatureSize = readList(inputFile, seqs);
            }
            else
            {
                signatureSize = readSignatures(inputFile, seqs);
            }

            fprintf(stderr, "Loaded %zu seqs...\n", seqs.size() / signatureSize);
            fprintf(pFile, "i,j,similarity\n");

            calcAllSimilarityKmers(pFile, seqs);
        }
        else
        {
            vector<seq_type> seqs;
            if (args.multiple_arg)
            {
                seqs = readListPartitionBF(inputFile, signatureSize);
            }
            else
            {
                seqs = readPartitionBF(inputFile, signatureSize);
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
                sprintf(buffer, "-t%zu-window_sim.txt", minimiser_match_threshold);
                pFile = fopen((rawname + buffer).c_str(), "w");
                fprintf(pFile, "i,j,similarity\n");

                calcAllSimilarityWindow(pFile, seqs);
            }
            else if (args.global_arg)
            {
                pFile = fopen((rawname + "-global_sim.txt").c_str(), "w");
                fprintf(pFile, "i,j,similarity\n");

                calcAllSimilarityGlobal(pFile, seqs);
            }
            else
            {
                pFile = fopen((rawname + "-local_sim.txt").c_str(), "w");
                fprintf(pFile, "i,j,similarity\n");

                calcAllSimilarityLocal(pFile, seqs);
            }
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
