// #include "bf_ktree.hpp"
#include <fstream>
#include <random>
#include "stats_main_cmdline.hpp"
#include "bloom_filter.hpp"
#include "read.hpp"
#include "self_tree.hpp"

using namespace std;

static stats_main_cmdline args; // Command line switches and arguments
// static size_t signatureSize;         // Signature size (depends on element in BF, obtained while read binary)
size_t max_seqCount = 100;

void calcAllstatsGlobal(vector<vector<vector<cell_type>>> seqs)
{
    size_t seqCount = seqs.size();
    for (size_t i = 0; i < seqCount; i++)
    {
        for (size_t j = 0; j < seqCount; j++)
        {
            double dist = calcJaccardGlobal(seqs[i], seqs[j]) * 100.0;
            fprintf(stderr, "%zu,%zu,%.2f\n", i, j, dist);
        }
    }
}

// Jaccard stats
void calcAllstatsKmers(vector<cell_type> seqs)
{
    size_t seqCount = seqs.size() / signatureSize;
    vector<size_t> indices(min(seqCount, max_seqCount));
    iota(indices.begin(), indices.end(), 0);

    for (size_t i : indices)
    {
        size_t temp = countSetBits(&seqs[i * signatureSize], signatureSize);
        for (size_t j : indices)
        {
            double sim = calcSimilarity(&seqs[i * signatureSize], &seqs[j * signatureSize], signatureSize);
            fprintf(stderr, "%zu,%zu,%.2f\n", i, j, sim * 100);
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

int stats_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    std::ios::sync_with_stdio(false); // No sync with stdio -> faster

    string bfIn = args.bf_input_arg;
    // size_t bf_element_cnt = 1000;

    size_t firstindex = bfIn.find_last_of("/") + 1;
    size_t lastindex = bfIn.find_last_of(".");
    string rawname = bfIn.substr(firstindex, lastindex - firstindex);

    if (args.output_given)
    {
        rawname = rawname + "-" + args.output_arg;
    }

    if (args.all_kmer_arg)
    {
        vector<cell_type> seqs;
        signatureSize = readSignatures(bfIn, seqs);

        fprintf(stderr, "Loaded %zu seqs...\n", seqs.size() / signatureSize);
        calcAllstatsKmers(seqs);
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
    }

    return 0;
}
