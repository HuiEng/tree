#include <random>
#include "partree_main_cmdline.hpp"
#include "part_ktree.hpp"

using namespace std;

static partree_main_cmdline args; // Command line switches and arguments

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

part_KTree clusterSignatures(const vector<vector<sig_type>> &seqs)
{
    size_t seqCount = seqs.size() - 1;
    vector<size_t> clusters(seqCount);
    part_KTree tree(partree_capacity);

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
        tree.first_insert(seqs[i], insertionList, i);
    }

    for (size_t i = firstNodes; i < seqCount; i++)
    {
        // fprintf(stdout, "inserting %zu", i);
        tree.next_insert(seqs[i], insertionList, i);
    }

    cout<<"var treeData=";
    tree.printSubTreeJson(tree.root);

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

void search(part_KTree tree, const vector<vector<sig_type>> &queries)
{
    // vector<sig_type> query;
    // query.push_back(0);
    // query.push_back(18);
    // query.push_back(35);
    // query.push_back(0);

    for (size_t i = 0; i < queries.size() - 1; i++)
    {
        cout << "Found " << i << " in: ";
        tree.search(queries[i]);
    }
}

int partree_main(int argc, char *argv[])
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

    vector<vector<sig_type>> seqs = readPartition(inputFile);
    // std::cout << "seqCount: " << seqCount << "\n";
    // for (auto seq : seqs)
    // {
    //     for (auto h : seq)
    //     {
    //         std::cout << h << ",";
    //     }
    //     std::cout << "\n";
    // }

    fprintf(stderr, "Loaded %zu signatures...\n", seqs.size());

    fprintf(stderr, "Building Signature...\n");
    default_random_engine rng;
    part_KTree tree = clusterSignatures(seqs);

    if (args.query_given)
    {
        fprintf(stderr, "Searching queries...\n");
        vector<vector<sig_type>> queries = readPartition(args.query_arg);
        search(tree, queries);
    }

    return 0;
}
