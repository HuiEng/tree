#include <random>
#include "tree_main_cmdline.hpp"
#include "sec_tree.hpp"
#include "cluster.hpp"
// #include "stats.hpp"

typedef sec_tree tree_type;

using namespace std;
static tree_main_cmdline args; // Command line switches and arguments



int tree_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    ios::sync_with_stdio(false); // No sync with stdio -> faster

    

    string inputFile = args.input_arg;
    string outName = setArgs(args);
    

    vector<size_t> clusters;
    vector<seq_type> seqs;
    if (args.multiple_arg)
    {
        seqs = readListPartitionBF(inputFile, signatureSize);
    }
    else
    {
        seqs = readPartitionBF(inputFile, signatureSize);
    }
    size_t seqCount = seqs.size();

    if (signatureSize == 0)
    {
        fprintf(stderr, "Something is wrong with the input data, please generate signature with diff params\n");
        return 0;
    }
    if (cap == 0)
    {
        cap = seqCount;
    }

    fprintf(stderr, "Loaded %zu seqs, signatureSize %zu...\n", seqs.size(), signatureSize);
    clusters = clusterSignatures<tree_type, seq_type>(seqs, seqCount);

    fprintf(stderr, "writing output...\n");
    FILE *pFile = fopen(outName.c_str(), "w");
    outputClusters(pFile, clusters);

    return 0;
}
