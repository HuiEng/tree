#include <random>
#include "ktree_main_cmdline.hpp"
#include "bf_ktree.hpp"
#include "read.hpp"

using namespace std;

static ktree_main_cmdline args; // Command line switches and arguments

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

vector<size_t> clusterSignatures(const vector<cell_type> &sigs)
{
    size_t sigCount = sigs.size() / signatureSize;
    vector<size_t> clusters(sigCount);
    BF_KTree tree(ktree_order, ktree_capacity);

    size_t firstNodes = 1;
    if (firstNodes > sigCount)
        firstNodes = sigCount;

    vector<size_t> insertionList;
    for (size_t i = 0; i <= firstNodes; i++)
    {
        insertionList.push_back(firstNodes - i);
    }

    default_random_engine rng;
    // Insert first 1 nodes single-threaded
    for (size_t i = 0; i < firstNodes; i++)
    {
        tree.insert(rng, &sigs[i * signatureSize], insertionList, i);
    }

    // What's the next free insertion point?
    size_t nextFree = insertionList.back();

    //#pragma omp parallel
    {
        default_random_engine rng;
        vector<size_t> insertionList;

        //#pragma omp for
        for (size_t i = nextFree; i < ktree_capacity; i++)
        {
            insertionList.push_back(ktree_capacity - i);
        }

        //#pragma omp for
        for (size_t i = firstNodes; i < sigCount; i++)
        {
            fprintf(stdout,"inserting %zu",i);
            tree.insert(rng, &sigs[i * signatureSize], insertionList, i);
        }
    }
    // tree.printTree();

    // We've created the tree. Now reinsert everything
    //#pragma omp parallel for

    // tree.get_counts();

    // vector<size_t> leaves;
    // tree.find_leaves(tree.root, leaves);
    for (size_t i = 0; i < sigCount; i++)
    {
        // size_t clus = tree.traverse_all(leaves, &sigs[i * signatureSize]);
        size_t clus = tree.traverse(&sigs[i * signatureSize]);
        clusters[i] = clus;
    }

    // We want to compress the cluster list down
    compressClusterList(clusters);
    tree.printTree();

    // vector<size_t> node;
    // node.push_back(tree.root);
    // size_t final = tree.traverse_all(&sigs[0],node,0,0);
    // fprintf(stderr,">>%zu,%zu\n",tree.root,final);

    // fprintf(stderr,">>%zu\n",tree.traverse_all(&sigs[0]));


    // Recursively destroy all locks
    tree.destroyLocks();

    // tree.do_kmeans(rng);


    return clusters;
}


void outputClusters(FILE *pFile, const vector<size_t> &clusters)
{
    for (size_t sig = 0; sig < clusters.size(); sig++)
    {
        fprintf(pFile, "%llu,%llu\n", static_cast<unsigned long long>(sig), static_cast<unsigned long long>(clusters[sig]));
    }
}


int ktree_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    std::ios::sync_with_stdio(false); // No sync with stdio -> faster

    if (!args.specific_arg && !args.input_given){
        cout<<"No input given! Exiting...\n";
        return 0;
    }

    if (args.order_given)
        ktree_order = args.order_arg;

    if (args.capacity_given)
        ktree_capacity = args.capacity_arg;

    if (args.threshold_given)
        threshold = args.threshold_arg;

    intersect = args.intersect_arg;

    string inputFile = args.input_arg;

    vector<cell_type> sigs;
    if (args.multiple_arg)
    {
        auto bf_len = readSignaturesMultiple(inputFile, sigs);
        signatureSize = bf_len;
    }
    else if (args.specific_given)
    {
        
        inputFile = args.specific_arg;
        auto bf_len = readSignaturesSpecific(inputFile, sigs);
        signatureSize = bf_len;
    }
    else
    {
        auto bf_len = readSignatures(inputFile, sigs);
        signatureSize = bf_len;
    }
    fprintf(stderr, "Loaded signatures...\n");

    fprintf(stderr, "Clustering signatures...\n");
    default_random_engine rng;
    kMeans(rng, sigs, 5);
    // auto clusters = clusterSignatures(sigs);
    // fprintf(stderr, "writing output...\n");

    // size_t firstindex = inputFile.find_last_of("/") + 1;
    // size_t lastindex = inputFile.find_last_of(".");
    // string rawname = inputFile.substr(firstindex, lastindex - firstindex);

    // auto fileName = rawname + "-o" + to_string(ktree_order);
    // FILE *pFile = fopen((fileName + ".txt").c_str(), "w");
    // outputClusters(pFile, clusters);

    return 0;
}
