// #include "bf_ktree.hpp"
#include <fstream>
#include <random>
#include "similarity_main_cmdline.hpp"
#include "bloom_filter.hpp"
#include "read.hpp"

using namespace std;

static similarity_main_cmdline args; // Command line switches and arguments
static size_t signatureSize;         // Signature size (depends on element in BF, obtained while read binary)

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

void calcAllSimilarity(const vector<cell_type> &sigs)
{
    size_t seqCount = sigs.size() / signatureSize;
    for (size_t i = 0; i < seqCount; i++)
    {
        for (size_t j = i + 1; j < seqCount; j++)
        {
            fprintf(stdout, "%zu,%zu,%f\n", i, j, calcSimilarity(&sigs[i * signatureSize], &sigs[j * signatureSize], signatureSize));
            // break;
        }
        // break;
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

    string bfIn = "out.bin";
    size_t bf_element_cnt = 1000;

    if (args.bf_input_given)
        bfIn = args.bf_input_arg;

    vector<cell_type> sigs;
    if (args.multiple_arg)
    {
        auto bf_len = readSignaturesMultiple(bfIn, sigs);
        signatureSize = bf_len;
    }
    else
    {
        auto bf_len = readSignatures(bfIn, sigs);
        signatureSize = bf_len;
    }
    fprintf(stderr, "Loaded signatures...\n");

    calcAllSetBits(sigs);
    calcAllSimilarity(sigs);

    for (int i = 0; i < sigs.size(); i += signatureSize)
    {
        printBF(&sigs[i]);
    }

    // for (auto sig:sigs){
    //     fprintf(stdout,"%zu\n",sig);
    // }

    return 0;
}
