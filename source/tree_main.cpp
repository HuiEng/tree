#include <random>
#include "tree_main_cmdline.hpp"
#include "self_tree.hpp"
typedef self_tree tree_type;

using namespace std;

static tree_main_cmdline args; // Command line switches and arguments

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

tree_type clusterSignatures(const vector<vector<vector<cell_type>>> &seqs)
{
    size_t seqCount = seqs.size();
    vector<size_t> clusters(seqCount);
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

    firstNodes = 2;

    // Insert first 1 nodes single-threaded
    for (size_t i = 0; i < firstNodes; i++)
    {
        tree.first_insert(seqs[i], insertionList, i);
        
    }


    for (size_t i = firstNodes; i < seqCount; i++)
    {
        // fprintf(stdout, "inserting %zu", i);
        tree.insert(seqs[i], insertionList, i, i%2+1);
    }

    tree.insert(seqs[0], insertionList, seqCount, 7);

    // cout<<"var treeData=";
    // tree.printSubTreeJson(tree.root);

    // Recursively destroy all locks
    tree.destroyLocks();

    return tree;
}

void outputClusters(FILE *pFile, const vector<size_t> &clusters)
{
    for (size_t sig = 0; sig < clusters.size(); sig++)
    {
        fprintf(pFile, "%llu,%llu\n", static_cast<unsigned long long>(sig), static_cast<unsigned long long>(clusters[sig]));
    }
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

    string inputFile = args.input_arg;

    vector<vector<vector<cell_type>>> seqs = readPartitionBF(inputFile);
    fprintf(stderr, "Loaded signatures...\n");

    signatureSize = seqs[0][0].size();
    fprintf(stderr, "Building Signature...\n");
    default_random_engine rng;
    tree_type tree = clusterSignatures(seqs);

    tree.printTreeJson(stdout);

    return 0;
}
