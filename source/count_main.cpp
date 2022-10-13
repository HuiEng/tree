// #include "bf_ktree.hpp"
#include <fstream>
#include <random>
#include "count_main_cmdline.hpp"
#include "bloom_filter.hpp"
#include "read.hpp"
#include "self_tree.hpp"

using namespace std;

static count_main_cmdline args; // Command line switches and arguments
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

// count set bits per seq
void countWinnow(FILE *pFile, vector<cell_type> seqs)
{
    size_t seqCount = seqs.size() / signatureSize;
    fprintf(pFile, "i,count\n");
    for (size_t i = 0; i < seqCount; i++)
    {
        size_t bits = countSetBits(&seqs[i * signatureSize], signatureSize);
        fprintf(pFile, "%zu,%zu\n", i, bits);
    }
}

void countChunk(FILE *pFile, vector<vector<vector<cell_type>>> seqs)
{
    size_t seqCount = seqs.size();
    for (size_t i = 0; i < seqCount; i++)
    {
        vector<cell_type> temp = getMinimiseSet(seqs[i])[0];
        size_t bits = countSetBits(&temp[0], signatureSize);
        fprintf(pFile, "%zu,%zu\n", i, bits);
    }
}

int count_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    std::ios::sync_with_stdio(false); // No sync with stdio -> faster

    string bfIn = args.bf_input_arg;

    size_t firstindex = bfIn.find_last_of("/") + 1;
    size_t lastindex = bfIn.find_last_of(".");
    string rawname = bfIn.substr(firstindex, lastindex - firstindex);

    if (args.output_given)
    {
        rawname = args.output_arg;
    }

    if (!args.chunk_arg)
    {
        vector<cell_type> seqs;
        signatureSize = readSignatures(bfIn, seqs);

        fprintf(stderr, "Loaded %zu seqs...\n", seqs.size() / signatureSize);
        FILE *pFile = fopen((rawname + "-all_cnt.txt").c_str(), "w");
        countWinnow(pFile, seqs);
    }
    else
    {
        vector<vector<vector<cell_type>>> seqs = readPartitionBF(bfIn);
        // dummy code, assume there is at least 10 valid input seqs
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
        FILE *pFile = fopen((rawname + "-all_cnt.txt").c_str(), "w");
        countChunk(pFile, seqs);
    }

    return 0;
}
