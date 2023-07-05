// int primary_tree_main(int argc, char *argv[]) { return 0; }
#include <random>
#include "primary_tree_main_cmdline.hpp"
#include "primary_tree.hpp"
// #include "stats.hpp"

typedef primary_tree primary_tree_type;
using namespace std;
static primary_tree_main_cmdline args; // Command line switches and arguments

int primary_tree_main(int argc, char *argv[])
{
    split_threshold = 0.5;
    stay_threshold = 0.8;
    minimiser_match_threshold = 4;
    // partree_capacity = 10000;

    args.parse(argc, argv);
    ios::sync_with_stdio(false); // No sync with stdio -> faster

    //?
    if (!args.input_given)
    {
        cout << "No input given! Exiting...\n";
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
        split_threshold = getSplitThresholdListSingle(inputFile);
    }
    else
    {
        split_threshold = getSplitThresholdSingle(inputFile);
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
    setTreeMeta(args);

    fprintf(stderr, "Loaded %zu seqs...signatureSize %zu\n", seqCount, signatureSize);
    clusters = clusterSignatures<primary_tree_type, cell_type>(seqs, seqCount, signatureSize);

    fprintf(stderr, "writing output...\n");

    size_t firstindex = inputFile.find_last_of("/") + 1;
    size_t lastindex = inputFile.find_last_of(".");
    string rawname = inputFile.substr(firstindex, lastindex - firstindex);

    auto fileName = rawname + "-s" + to_string((int)(stay_threshold * 100)) + "-l" + to_string((int)(split_threshold * 100));
    if (args.tag_given)
    {
        fileName = args.tag_arg + fileName;
    }
    FILE *pFile = fopen((fileName + ".txt").c_str(), "w");
    outputClusters(pFile, clusters);

    return 0;
}
