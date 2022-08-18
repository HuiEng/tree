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

void calcAllSimilarity(vector<vector<vector<cell_type>>> seqs)
{
    size_t seqCount = seqs.size();
    fprintf(stdout, "i,j,similarity\n");
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

        for (size_t j = 0; j < seqCount; j++)
        {
            double dist = calcDistance(seqs[i], seqs[j]) * 100.0 / max(countBits(seqs[i]),countBits(seqs[j]));
            fprintf(stdout, "%zu,%zu,%.2f\n", i, j, dist);
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

    vector<vector<vector<cell_type>>> seqs = readPartitionBF(bfIn);
    signatureSize = seqs[0][0].size();
    size_t seqCount = seqs.size();
    fprintf(stderr, "Loaded seqs...\n");

    // calcAllSetBits(sigs);
    calcAllSimilarity(seqs);

    // for (int i = 0; i < sigs.size(); i += signatureSize)
    // {
    //     printBF(&sigs[i]);
    // }

    // for (auto sig:sigs){
    //     fprintf(stdout,"%zu\n",sig);
    // }

    return 0;
}
