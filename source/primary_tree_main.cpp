<<<<<<< HEAD
// int primary_tree_main(int argc, char *argv[]) { return 0; }
=======
>>>>>>> 7d7dfbb5ca924824ed232a598c0f72d870ac4d19
#include <random>
#include "primary_tree_main_cmdline.hpp"
#include "primary_tree.hpp"
#include "stats.hpp"

typedef primary_tree primary_tree_type;
typedef vector<tuple<size_t, tuple<size_t, size_t>>> output_type;
bool random_ = false;
size_t cap = 0;
size_t iteration = 0;
bool iteration_given = false;
using namespace std;
bool debug_ = false;
bool force_split_ = false;
unsigned seed;

static primary_tree_main_cmdline args; // Command line switches and arguments

// void dbgPrintMatrix(const cell_type *matrix)
// {
//     size_t ktree_csig_height = (ktree_order + (bits_per_char -1)) / bits_per_char;
//     for (size_t i = 0; i < signatureSize * bits_per_char; i++)
//     {
//         fprintf(stderr, "%03zu:", i);
//         for (size_t j = 0; j < ktree_csig_height * bits_per_char; j++)
//         {
//             auto val = matrix[i * ktree_csig_height + (j / bits_per_char)];
//             if (val & (1ull << (j % bits_per_char)))
//             {
//                 fprintf(stderr, "1");
//             }
//             else
//             {
//                 fprintf(stderr, "0");
//             }
//         }
//         fprintf(stderr, "\n");
//         if (i >= 5)
//         {
//             fprintf(stderr, "...............\n");
//             break;
//         }
//     }
// }

void outputClusters(FILE *pFile, const vector<size_t> &clusters)
{
    for (size_t sig = 0; sig < clusters.size(); sig++)
    {
        fprintf(pFile, "%llu,%llu\n", static_cast<unsigned long long>(sig), static_cast<unsigned long long>(clusters[sig]));
    }
}

void outputClusters(FILE *pFile, const output_type &clusters)
{
    fprintf(pFile, "seqID,cluster,ancestor,level\n");
    for (size_t sig = 0; sig < clusters.size(); sig++)
    {
        fprintf(pFile, "%llu,%llu,%llu,%llu\n",
                static_cast<unsigned long long>(sig),
                static_cast<unsigned long long>(get<0>(clusters[sig])),
                static_cast<unsigned long long>(get<0>(get<1>(clusters[sig]))),
                static_cast<unsigned long long>(get<1>(get<1>(clusters[sig]))));
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

<<<<<<< HEAD
// vector<size_t> clusterSignatures2(const vector<cell_type> &seqs)
// {
//     size_t seqCount = seqs.size() / signatureSize;
//     vector<size_t> clusters(cap);
//     primary_tree_type primary_tree(partree_capacity);
//     // primary_tree.means.resize(partree_capacity * signatureSize);

//     size_t firstNodes = 1;
//     if (firstNodes > seqCount)
//         firstNodes = seqCount;

//     vector<size_t> insertionList; // potential nodes idx except root; root is always 0

//     default_random_engine rng;

//     // node 0 reserved for root, node 1 reserved for leaves idx
//     for (size_t i = firstNodes; i < partree_capacity; i++)
//     {
//         insertionList.push_back(partree_capacity - i);
//     }
//     vector<size_t> foo;
//     for (int i = 0; i < seqCount; i++)
//     {
//         foo.push_back(i);
//     }

//     if (random_)
//     {
//         // unsigned seed = chrono::system_clock::now().time_since_epoch().count();
//         shuffle(foo.begin(), foo.end(), default_random_engine(seed));
//     }
//     foo.resize(cap);

//     clusters[foo[0]] = primary_tree.first_insert(&seqs[foo[0] * signatureSize], insertionList, foo[0]);
//     if (force_split_)
//     {
//         for (size_t i = 1; i < cap; i++)
//         {
//             printMsg("inserting %zu\n", foo[i]);
//             size_t clus = primary_tree.insertSplitRoot(&seqs[foo[i] * signatureSize], insertionList, foo[i]) + 1;
//         }
//     }
//     else
//     {
//         for (size_t i = 1; i < cap; i++)
//         {
//             printMsg("inserting %zu\n", foo[i]);
//             size_t clus = primary_tree.insert(&seqs[foo[i] * signatureSize], insertionList, foo[i]);
//             clusters[foo[i]] = clus;
//         }
//     }

//     // // for debugging
//     // if (iteration_given)
//     // {
//     //     // prep to remove and reinsert ambi
//     //     fprintf(stderr, "\n\n\nBefore\n");
//     //     primary_tree.printTreeJson(stderr);

//     //     singleton = 0;
//     //     // primary_tree.trim();
//     //     primary_tree.removeAmbi();
//     //     singleton = 2;
//     //     printMsg("\n\nReinserting ambi (all)\n");
//     //     primary_tree.prepReinsert();
//     //     for (size_t i = 0; i < cap; i++)
//     //     {
//     //         size_t clus = primary_tree.reinsert(&seqs[foo[i] * signatureSize], foo[i]);
//     //         printMsg("\n Reinsert %zu at %zu\n", foo[i], clus);
//     //         // clusters[foo[i]] = primary_tree.findAncestor(clus);
//     //         clusters[foo[i]] = clus;
//     //     }

//     //     primary_tree.printTreeJson(stderr);
//     // }

//     // for (size_t run = 0; run < iteration; run++)
//     // {
//     //     fprintf(stderr, "Iteration %zu (singleton = %zu)\n", run, singleton);

//     //     // primary_tree.trim();

//     //     primary_tree.removeAmbi();
//     //     primary_tree.prepReinsert();
//     //     // singleton++;

//     //     primary_tree.printTreeJson(stderr);
//     //     for (size_t i = 0; i < cap; i++)
//     //     {
//     //         size_t clus = primary_tree.reinsert(&seqs[foo[i] * signatureSize], foo[i]);

//     //         printMsg("\n found %zu at %zu\n", foo[i], clus);
//     //         // clusters[foo[i]] = primary_tree.findAncestor(clus);
//     //         clusters[foo[i]] = clus;
//     //     }

//     //     if (debug_)
//     //     {
//     //         auto fileName = "nodeDistance-r" + to_string((size_t)(run)) + ".txt";
//     //         FILE *nFile = fopen(fileName.c_str(), "w");
//     //         primary_tree.printNodeDistance(nFile, seqs, clusters);

//     //         fileName = "clusters-r" + to_string((size_t)(run)) + ".txt";
//     //         FILE *cFile = fopen(fileName.c_str(), "w");
//     //         outputClusters(cFile, clusters);
//     //     }
//     // }
//     // primary_tree.updateTree();

//     // FILE *pFile = fopen("nodeDistance.txt", "w");
//     // primary_tree.printNodeDistance(pFile, seqs, clusters);

//     // Recursively destroy all locks
//     primary_tree.destroyLocks();

//     primary_tree.printTreeJson(stdout);
//     // FILE *hFile = fopen("hierarchy.txt", "w");
//     // fprintf(hFile, "parent,child,rank\n");
//     // primary_tree.outputHierarchy(hFile);

//     return clusters;
// }
=======
vector<size_t> clusterSignatures2(const vector<cell_type> &seqs)
{
    size_t seqCount = seqs.size() / signatureSize;
    vector<size_t> clusters(cap);
    primary_tree_type primary_tree(partree_capacity);
    // primary_tree.means.resize(partree_capacity * signatureSize);

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
    vector<size_t> foo;
    for (int i = 0; i < seqCount; i++)
    {
        foo.push_back(i);
    }

    if (random_)
    {
        // unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        shuffle(foo.begin(), foo.end(), default_random_engine(seed));
    }
    foo.resize(cap);

    clusters[foo[0]] = primary_tree.first_insert(&seqs[foo[0] * signatureSize], insertionList, foo[0]);
    if (force_split_)
    {
        for (size_t i = 1; i < cap; i++)
        {
            printMsg("inserting %zu\n", foo[i]);
            size_t clus = primary_tree.insertSplitRoot(&seqs[foo[i] * signatureSize], insertionList, foo[i]) + 1;
        }
    }
    else
    {
        for (size_t i = 1; i < cap; i++)
        {
            printMsg("inserting %zu\n", foo[i]);
            size_t clus = primary_tree.insert(&seqs[foo[i] * signatureSize], insertionList, foo[i]);
            clusters[foo[i]] = clus;
        }
    }

    // // for debugging
    // if (iteration_given)
    // {
    //     // prep to remove and reinsert ambi
    //     fprintf(stderr, "\n\n\nBefore\n");
    //     primary_tree.printTreeJson(stderr);

    //     singleton = 0;
    //     // primary_tree.trim();
    //     primary_tree.removeAmbi();
    //     singleton = 2;
    //     printMsg("\n\nReinserting ambi (all)\n");
    //     primary_tree.prepReinsert();
    //     for (size_t i = 0; i < cap; i++)
    //     {
    //         size_t clus = primary_tree.reinsert(&seqs[foo[i] * signatureSize], foo[i]);
    //         printMsg("\n Reinsert %zu at %zu\n", foo[i], clus);
    //         // clusters[foo[i]] = primary_tree.findAncestor(clus);
    //         clusters[foo[i]] = clus;
    //     }

    //     primary_tree.printTreeJson(stderr);
    // }

    // for (size_t run = 0; run < iteration; run++)
    // {
    //     fprintf(stderr, "Iteration %zu (singleton = %zu)\n", run, singleton);

    //     // primary_tree.trim();

    //     primary_tree.removeAmbi();
    //     primary_tree.prepReinsert();
    //     // singleton++;

    //     primary_tree.printTreeJson(stderr);
    //     for (size_t i = 0; i < cap; i++)
    //     {
    //         size_t clus = primary_tree.reinsert(&seqs[foo[i] * signatureSize], foo[i]);

    //         printMsg("\n found %zu at %zu\n", foo[i], clus);
    //         // clusters[foo[i]] = primary_tree.findAncestor(clus);
    //         clusters[foo[i]] = clus;
    //     }

    //     if (debug_)
    //     {
    //         auto fileName = "nodeDistance-r" + to_string((size_t)(run)) + ".txt";
    //         FILE *nFile = fopen(fileName.c_str(), "w");
    //         primary_tree.printNodeDistance(nFile, seqs, clusters);

    //         fileName = "clusters-r" + to_string((size_t)(run)) + ".txt";
    //         FILE *cFile = fopen(fileName.c_str(), "w");
    //         outputClusters(cFile, clusters);
    //     }
    // }
    // primary_tree.updateTree();

    // FILE *pFile = fopen("nodeDistance.txt", "w");
    // primary_tree.printNodeDistance(pFile, seqs, clusters);

    // Recursively destroy all locks
    primary_tree.destroyLocks();

    primary_tree.printTreeJson(stdout);
    // FILE *hFile = fopen("hierarchy.txt", "w");
    // fprintf(hFile, "parent,child,rank\n");
    // primary_tree.outputHierarchy(hFile);

    return clusters;
}
>>>>>>> 7d7dfbb5ca924824ed232a598c0f72d870ac4d19

int primary_tree_main(int argc, char *argv[])
{
    split_threshold = 0.5;
    stay_threshold = 0.8;
    minimiser_match_threshold = 4;
<<<<<<< HEAD
    // partree_capacity = 10000;
=======
    partree_capacity = 10000;
>>>>>>> 7d7dfbb5ca924824ed232a598c0f72d870ac4d19

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
    else // if (args.single_arg)
    {
        split_threshold = getSplitThresholdSingle(inputFile);
    }
    // else
    // {
    //     split_threshold = getSplitThreshold(inputFile);
    // }

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

    cap = args.capacity_arg;
    if (args.iteration_given)
    {
        iteration_given = true;
        iteration = args.iteration_arg;
    }

    debug_ = args.debug_arg;
    print_ = args.print_arg;
    force_split_ = args.force_split_arg;

    // signatureSize = seqs[0][0].size();
    default_random_engine rng;
    vector<size_t> clusters;
    vector<cell_type> seqs;
    signatureSize = readSignatures(inputFile, seqs);
    signatureWidth = signatureSize * sizeof(cell_type);
    if (cap == 0)
    {
        cap = seqs.size();
    }

    fprintf(stderr, "Loaded %zu seqs...signatureSize %zu\n", seqs.size() / signatureSize, signatureSize);
<<<<<<< HEAD
    // clusters = clusterSignatures2(seqs);
=======
    clusters = clusterSignatures2(seqs);
>>>>>>> 7d7dfbb5ca924824ed232a598c0f72d870ac4d19

    fprintf(stderr, "writing output...\n");

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
