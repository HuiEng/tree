#include <random>
#include "tree_main_cmdline.hpp"
#include "sec_tree.hpp"
// #include "stats.hpp"

typedef sec_tree tree_type;

using namespace std;
static tree_main_cmdline args; // Command line switches and arguments



int tree_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    ios::sync_with_stdio(false); // No sync with stdio -> faster

    

    string inputFile = args.input_arg;
    setArgs(args);
    

    if (args.split_threshold_given)
    {
        split_threshold = args.split_threshold_arg;
    }
    else if (args.multiple_arg)
    {
        split_threshold = getSplitThresholdList(inputFile);
    }
    else
    {
        split_threshold = getSplitThreshold(inputFile);
    }

    split_node_threshold = split_threshold / 2;

    fprintf(stderr, "split threshold: %.2f\n", split_threshold);
    fprintf(stderr, "stay threshold: %.2f\n", stay_threshold);
    fprintf(stderr, "split_node_threshold threshold: %.2f\n", split_node_threshold);

    // signatureSize = seqs[0][0].size();
    default_random_engine rng;
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
    setArgs(args);

    fprintf(stderr, "Loaded %zu seqs, signatureSize %zu...\n", seqs.size(), signatureSize);
    clusters = clusterSignatures<tree_type, seq_type>(seqs, seqCount);

    size_t firstindex = inputFile.find_last_of("/") + 1;
    size_t lastindex = inputFile.find_last_of(".");
    string rawname = inputFile.substr(firstindex, lastindex - firstindex);

    auto fileName = rawname + "-s" + to_string((int)(stay_threshold * 100)) + "-l" + to_string((int)(split_threshold * 100));
    if (args.tag_given)
    {
        fileName = fileName + "-" + args.tag_arg;
    }
    FILE *pFile = fopen((fileName + ".txt").c_str(), "w");
    outputClusters(pFile, clusters);

    return 0;
}
