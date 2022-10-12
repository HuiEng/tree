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

void toBinary(cell_type letter)
{
    int binary[8];
    for (int n = 0; n < 8; n++)
        binary[7 - n] = (letter >> n) & 1;

    for (int n = 0; n < 8; n++)
        fprintf(stdout, "%d", binary[n]);
}

void printBF(const cell_type *bf)
{
    for (std::size_t i = 0; i < signatureSize; ++i)
    {
        toBinary(bf[i * sizeof(cell_type)]);
    }
    fprintf(stdout, "\n");
}

void calcAllSimilarity(FILE *pFile, vector<vector<vector<cell_type>>> seqs)
{
    size_t seqCount = seqs.size();
    fprintf(pFile, "i,j,similarity\n");
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

        for (size_t j = i; j < seqCount; j++)
        {
            // double dist = calcDistance(seqs[i], seqs[j]) * 100.0 / max(countBits(seqs[i]), countBits(seqs[j]));
            double dist = calcJaccardLocal(seqs[i], seqs[j]) * 100.0;
            fprintf(pFile, "%zu,%zu,%.2f\n", i, j, dist);
        }
    }
}

void calcAllSimilarityGlobal(FILE *pFile, vector<vector<vector<cell_type>>> seqs)
{
    size_t seqCount = seqs.size();
    fprintf(pFile, "i,j,similarity\n");
    for (size_t i = 0; i < seqCount; i++)
    {
        for (size_t j = i; j < seqCount; j++)
        {
            double dist = calcJaccardGlobal(seqs[i], seqs[j]) * 100.0;
            fprintf(pFile, "%zu,%zu,%.2f\n", i, j, dist);
        }
    }
}

// Jaccard similarity
void calcAllSimilarityKmers(FILE *pFile, vector<cell_type> seqs)
{
    size_t seqCount = seqs.size() / signatureSize;
    fprintf(pFile, "i,j,similarity\n");
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

    size_t firstindex = bfIn.find_last_of("/") + 1;
    size_t lastindex = bfIn.find_last_of(".");
    string rawname = bfIn.substr(firstindex, lastindex - firstindex);

    if (args.output_given)
    {
        rawname = rawname + "-" + args.output_arg;
    }

    if (args.all_kmer_arg)
    {
        FILE *pFile = fopen((rawname + "-all_sim.txt").c_str(), "w");
        vector<cell_type> seqs;
        signatureSize = readSignatures(bfIn, seqs);

        fprintf(stderr, "Loaded %zu seqs...\n", seqs.size() / signatureSize);
        calcAllSimilarityKmers(pFile, seqs);
    }
    else
    {
        vector<vector<vector<cell_type>>> seqs = readPartitionBF(bfIn);
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
        if (args.global_arg)
        {
            FILE *pFile = fopen((rawname + "-global_sim.txt").c_str(), "w");
            calcAllSimilarityGlobal(pFile, seqs);
        }
        else
        {
            FILE *pFile = fopen((rawname + "_sim.txt").c_str(), "w");
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
