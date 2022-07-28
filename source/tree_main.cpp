#include <random>
#include "tree_main_cmdline.hpp"
#include "self_tree.hpp"
typedef self_tree tree_type;

using namespace std;

static tree_main_cmdline args; // Command line switches and arguments

void outputClusters(FILE *pFile, const vector<tuple<size_t, size_t>> &clusters)
{
    fprintf(pFile, "seqID,cluster,ancestor\n");
    for (size_t sig = 0; sig < clusters.size(); sig++)
    {
        fprintf(pFile, "%llu,%llu,%llu\n",
                static_cast<unsigned long long>(sig),
                static_cast<unsigned long long>(get<0>(clusters[sig])),
                static_cast<unsigned long long>(get<1>(clusters[sig])));
    }
}

void compressClusterList(vector<size_t> &clusters)
{
    unordered_map<size_t, size_t> remap;
    for (size_t &clus : clusters)
    {
        if (remap.count(clus))
        {
            clus = remap[clus];
        }
        else
        {
            size_t newClus = remap.size();
            remap[clus] = newClus;
            clus = newClus;
        }
    }
    fprintf(stderr, "Output %zu clusters\n", remap.size());
}

vector<tuple<size_t, size_t>> clusterSignatures(const vector<vector<vector<cell_type>>> &seqs)
{
    size_t seqCount = seqs.size();
    vector<tuple<size_t, size_t>> clusters(seqCount);
    tree_type tree(partree_capacity);

    size_t firstNodes = 1;
    if (firstNodes > seqCount)
        firstNodes = seqCount;

    vector<size_t> insertionList; // potential nodes idx except root; root is always 0

    default_random_engine rng;

    // node 0 reserved for root, node 1 reserved for leaves idx
    for (size_t i = firstNodes; i < partree_capacity; i++)
    {
        insertionList.push_back(partree_capacity - i);
    }

    // Insert first 1 nodes single-threaded
    for (size_t i = 0; i < firstNodes; i++)
    {
        size_t clus = tree.first_insert(seqs[i], insertionList, i);
        clusters[i] = make_tuple(clus, clus);
    }

    for (size_t i = firstNodes; i < seqCount; i++)
    {
        // fprintf(stdout, "inserting %zu", i);
        size_t clus = tree.insert(seqs[i], insertionList, i);
        // clusters[i] = tree.findAncestor(clus);
        clusters[i] = make_tuple(clus, tree.findAncestor(clus));
    }

    for (size_t i = 0; i < seqCount; i++)
    {
        // fprintf(stdout, "inserting %zu", i);
        size_t clus = tree.search(seqs[i], i);
        // clusters[i] = tree.findAncestor(clus);
        clusters[i] = make_tuple(clus, tree.findAncestor(clus));
    }

    // Recursively destroy all locks
    tree.destroyLocks();

    tree.printTreeJson(stdout);

    return clusters;
}

int tree_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    std::ios::sync_with_stdio(false); // No sync with stdio -> faster

    //?
    if (!args.input_given)
    {
        cout << "No input and/or query given! Exiting...\n";
        return 0;
    }

    if (args.split_threshold_given)
    {
        split_threshold = args.split_threshold_arg;
        fprintf(stderr, "split threshold: %f\n", split_threshold);
    }

    if (args.stay_threshold_given)
    {
        stay_threshold = args.stay_threshold_arg;
        fprintf(stderr, "stay threshold: %f\n", stay_threshold);
    }

    string inputFile = args.input_arg;

    vector<vector<vector<cell_type>>> seqs = readPartitionBF(inputFile);
    fprintf(stderr, "Loaded signatures...\n");

    signatureSize = seqs[0][0].size();
    fprintf(stderr, "Building Signature...\n");
    default_random_engine rng;
    vector<tuple<size_t, size_t>> clusters = clusterSignatures(seqs);

    fprintf(stderr, "writing output...\n");

    size_t firstindex = inputFile.find_last_of("/") + 1;
    size_t lastindex = inputFile.find_last_of(".");
    string rawname = inputFile.substr(firstindex, lastindex - firstindex);

    auto fileName = rawname + "-s" + to_string((int)stay_threshold) + "-l" + to_string((int)split_threshold);
    FILE *pFile = fopen((fileName + ".txt").c_str(), "w");
    outputClusters(pFile, clusters);

    return 0;
}
