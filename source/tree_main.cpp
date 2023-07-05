#include <random>
#include "tree_main_cmdline.hpp"
#include "sec_tree.hpp"
// #include "stats.hpp"

typedef sec_tree tree_type;

using namespace std;
static tree_main_cmdline args; // Command line switches and arguments



int tree_main(int argc, char *argv[])
{
    split_threshold = 0.5;
    stay_threshold = 0.8;
    minimiser_match_threshold = 4;
    partree_capacity = 10000;

    args.parse(argc, argv);
    ios::sync_with_stdio(false); // No sync with stdio -> faster

    //?
    if (!args.input_given)
    {
        cout << "No input and/or query given! Exiting...\n";
        return 0;
    }

    string inputFile = args.input_arg;
    if (args.random_arg)
    {
        random_ = true;
        seed = args.seed_arg;
    }

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

    if (args.stay_threshold_given)
    {
        stay_threshold = args.stay_threshold_arg;
    }
    split_node_threshold = split_threshold / 2;

    fprintf(stderr, "split threshold: %.2f\n", split_threshold);
    fprintf(stderr, "stay threshold: %.2f\n", stay_threshold);
    fprintf(stderr, "split_node_threshold threshold: %.2f\n", split_node_threshold);

    if (args.minimiser_match_given)
    {
        minimiser_match_threshold = args.minimiser_match_arg;
        fprintf(stderr, "minimiser_match threshold: %zu\n", minimiser_match_threshold);
    }

    cap = args.sizeCap_arg;
    if (args.iteration_given)
    {
        iteration_given = true;
        iteration = args.iteration_arg;
    }

    if (args.capacity_given)
    {
        partree_capacity = args.capacity_arg;
        fprintf(stderr, "partree_capacity: %zu\n", partree_capacity);
    }

    debug_ = args.debug_arg;
    print_ = args.print_arg;
    force_split_ = args.force_split_arg;
    if (args.tree_order_given)
    {
        tree_order = args.tree_order_arg;
        fprintf(stderr, "tree_order: %zu\n", tree_order);
    }

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
    setTreeMeta(args);

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
