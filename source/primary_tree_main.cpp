// int primary_tree_main(int argc, char *argv[]) { return 0; }
#include <random>
#include "tree_main_cmdline.hpp"
#include "primary_tree.hpp"
#include "cluster.hpp"
// #include "stats.hpp"

typedef primary_tree primary_tree_type;
using namespace std;
static tree_main_cmdline args; // Command line switches and arguments

int primary_tree_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    ios::sync_with_stdio(false); // No sync with stdio -> faster

    string inputFile = args.input_arg;
    args.single_arg = true;
    string outName = setArgs(args);


    vector<size_t> clusters;
    vector<cell_type> seqs;

    if (args.multiple_arg)
    {
        // clusters = clusterSignaturesList(inputFile);
        signatureSize = readList(inputFile, seqs);
    }
    else
    {
        signatureSize = readSignatures(inputFile, seqs);

        // signatureWidth = signatureSize * sizeof(cell_type);
        // size_t seqCount = seqs.size() / signatureSize;
        // if (cap == 0)
        // {
        //     cap = seqCount;
        // }

        // fprintf(stderr, "Loaded %zu seqs...signatureSize %zu\n", seqCount, signatureSize);
        // clusters = clusterSignaturesPrim(seqs);
    }

    signatureWidth = signatureSize * sizeof(cell_type);
    size_t seqCount = seqs.size() / signatureSize;
    if (cap == 0)
    {
        cap = seqCount;
    }

    fprintf(stderr, "Loaded %zu seqs...signatureSize %zu\n", seqCount, signatureSize);
    clusters = clusterSignatures<primary_tree_type, cell_type>(seqs, seqCount, signatureSize);

    fprintf(stderr, "writing output...\n");
    FILE *pFile = fopen(outName.c_str(), "w");
    outputClusters(pFile, clusters);

    return 0;
}
