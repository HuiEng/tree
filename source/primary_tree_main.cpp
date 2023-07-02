// int primary_tree_main(int argc, char *argv[]) { return 0; }
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

string readTreeLine(string s, string folder, primary_tree_type &primary_tree)
{
    size_t parent = 0;
    size_t child = 0;
    string delimiter = ">";

    size_t pos = s.find(delimiter);
    string node = s.substr(0, pos);
    sscanf(node.c_str(), "%zu", &parent);
    s.erase(0, pos + 1);
    // cout << "Node: " << node << endl;

    delimiter = ",";
    string childStr = "";
    while ((pos = s.find(delimiter)) != string::npos)
    {
        childStr = s.substr(0, pos);
        sscanf(childStr.c_str(), "%zu", &child);
        // cout << "Child: " << childStr << endl;
        s.erase(0, pos + 1);

        vector<cell_type> signature;
        readSignatures((folder + childStr + ".bin"), signature);
        primary_tree.readNode(parent, child, &signature[0]);
    }
    primary_tree.updatePriority(parent);

    // return parent;
    return node;
}

// output signature size, if file too big => return empty "sigs", use readSignaturesBatch
size_t readTree(const string folder, primary_tree_type &primary_tree)
{
    string line;

    // Read from the text file
    ifstream listStream((folder + "tree.txt").c_str());

    // read signatureSize
    getline(listStream, line);

    // // use first line to set tree params
    // primary_tree_type primary_tree(partree_capacity);
    // getline(listStream, line);
    // sscanf(line.c_str(), "%zu", &signatureSize);
    // signatureWidth = signatureSize * sizeof(cell_type);
    // primary_tree.means.resize(partree_capacity * signatureSize);

    // read last nodeIdx
    getline(listStream, line);
    size_t offset = 0;
    sscanf(line.c_str(), "%zu", &offset);

    while (getline(listStream, line))
    {
        readTreeLine(line, folder, primary_tree);
    }

    // Close the file
    listStream.close();
    // primary_tree.updateTree();
    // primary_tree.printTreeJson(stdout);
    // return primary_tree;

    return offset;
}

vector<size_t> clusterSignaturesPrim(const vector<cell_type> &seqs)
{
    size_t seqCount = seqs.size() / signatureSize;
    vector<size_t> clusters(cap);

    size_t firstNodes = 1;
    if (firstNodes > seqCount)
        firstNodes = seqCount;

    vector<size_t> insertionList; // potential nodes idx except root; root is always 0
    size_t offset = 0;
    default_random_engine rng;

    primary_tree_type primary_tree(partree_capacity);
    primary_tree.means.resize(partree_capacity * signatureSize);

    if (args.topology_in_given)
    {
        offset = readTree(args.topology_in_arg, primary_tree);
    }

    // node 0 reserved for root, node 1 reserved for leaves idx
    for (size_t i = firstNodes; i < partree_capacity - offset; i++)
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
            size_t clus = primary_tree.insertSplitRoot(&seqs[foo[i] * signatureSize], insertionList, foo[i]);
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

    primary_tree.printTreeJson(stderr);

    // for debugging
    if (iteration_given)
    {
        // prep to remove and reinsert ambi
        fprintf(stderr, "\n\n\nBefore\n");
        primary_tree.printTreeJson(stderr);

        singleton = 0;
        // primary_tree.trim();
        primary_tree.removeAmbi();
        singleton = 2;
        printMsg("\n\nReinserting ambi (all)\n");
        primary_tree.prepReinsert();
        for (size_t i = 0; i < cap; i++)
        {
            size_t clus = primary_tree.reinsert(&seqs[foo[i] * signatureSize], foo[i]);
            printMsg("\n Reinsert %zu at %zu\n", foo[i], clus);
            // clusters[foo[i]] = primary_tree.findAncestor(clus);
            clusters[foo[i]] = clus;
        }

        primary_tree.printTreeJson(stderr);
    }

    for (size_t run = 0; run < iteration; run++)
    {
        fprintf(stderr, "Iteration %zu (singleton = %zu)\n", run, singleton);

        // primary_tree.trim();

        primary_tree.removeAmbi();
        primary_tree.prepReinsert();
        // singleton++;

        primary_tree.printTreeJson(stderr);
        for (size_t i = 0; i < cap; i++)
        {
            size_t clus = primary_tree.reinsert(&seqs[foo[i] * signatureSize], foo[i]);

            printMsg("\n found %zu at %zu\n", foo[i], clus);
            // clusters[foo[i]] = primary_tree.findAncestor(clus);
            clusters[foo[i]] = clus;
        }

        if (debug_)
        {
            auto fileName = "nodeDistance-r" + to_string((size_t)(run)) + ".txt";
            FILE *nFile = fopen(fileName.c_str(), "w");
            primary_tree.printNodeDistance(nFile, seqs, clusters);

            fileName = "clusters-r" + to_string((size_t)(run)) + ".txt";
            FILE *cFile = fopen(fileName.c_str(), "w");
            outputClusters(cFile, clusters);
        }
    }
    primary_tree.updateTree();

    FILE *pFile = fopen("nodeDistance.txt", "w");
    primary_tree.printNodeDistance(pFile, seqs, clusters);

    // Recursively destroy all locks
    primary_tree.destroyLocks();

    primary_tree.printTreeJson(stdout);
    FILE *hFile = fopen("hierarchy.txt", "w");
    fprintf(hFile, "parent,child,rank\n");
    primary_tree.outputHierarchy(hFile);

    if (args.topology_out_given)
    {
        string outFile = args.topology_out_arg;
        FILE *tFile = fopen((outFile + "tree.txt").c_str(), "w");
        primary_tree.printTree(tFile, insertionList, args.topology_out_arg);
    }

    return clusters;
}

// assume 1 file 1 seq
vector<size_t> clusterSignaturesList(const string listFile)
{

    string file;

    // Read from the text file
    ifstream listStream(listFile);
    size_t seqCount = getListLength(listStream);
    if (cap == 0)
    {
        cap = seqCount;
    }

    // read first file and set tree params;
    getline(listStream, file);
    vector<cell_type> seq;
    signatureSize = readSignatures(file, seq);
    signatureWidth = signatureSize * sizeof(cell_type);

    vector<size_t> clusters(seqCount);
    primary_tree_type primary_tree(partree_capacity);
    primary_tree.means.resize(partree_capacity * signatureSize);

    size_t firstNodes = 1;
    if (firstNodes > seqCount)
        firstNodes = seqCount;

    vector<size_t> insertionList; // potential nodes idx except root; root is always 0

    // node 0 reserved for root, node 1 reserved for leaves idx
    for (size_t i = firstNodes; i < partree_capacity; i++)
    {
        insertionList.push_back(partree_capacity - i);
    }

    clusters[0] = primary_tree.first_insert(&seq[0], insertionList, 0);

    // Use a while loop together with the getline() function to read the file line by line
    for (size_t i = 1; i < cap; i++)
    {
        getline(listStream, file);
        readSignatures(file, seq);
        printMsg("inserting %zu\n", i);
        size_t clus = 0;
        if (force_split_)
        {
            clus = primary_tree.insertSplitRoot(&seq[0], insertionList, i);
        }
        else
        {
            clus = primary_tree.insert(&seq[0], insertionList, i);
        }

        // clusters.push_back(clus);
        clusters[i] = clus;
    }

    fprintf(stderr, "Loaded %zu seqs...signatureSize %zu\n", seqCount, signatureSize);
    cap = max(cap, seqCount);

    // for debugging
    if (iteration_given)
    {
        // prep to remove and reinsert ambi
        fprintf(stderr, "\n\n\nBefore\n");
        primary_tree.printTreeJson(stderr);

        singleton = 0;
        // primary_tree.trim();
        primary_tree.removeAmbi();
        singleton = 2;
        printMsg("\n\nReinserting ambi (all)\n");
        primary_tree.prepReinsert();

        ifstream listStream_(listFile);
        for (size_t i = 0; i < cap; i++)
        {
            getline(listStream_, file);
            readSignatures(file, seq);
            size_t clus = primary_tree.reinsert(&seq[0], i);
            printMsg("\n Reinsert %zu at %zu\n", i, clus);
            clusters[i] = clus;
        }

        primary_tree.printTreeJson(stderr);
    }

    for (size_t run = 0; run < iteration; run++)
    {
        fprintf(stderr, "Iteration %zu (singleton = %zu)\n", run, singleton);
        // primary_tree.trim();

        primary_tree.removeAmbi();
        primary_tree.prepReinsert();
        // singleton++;

        primary_tree.printTreeJson(stderr);
        ifstream listStream_(listFile);
        for (size_t i = 0; i < cap; i++)
        {
            getline(listStream_, file);
            readSignatures(file, seq);
            size_t clus = primary_tree.reinsert(&seq[0], i);

            printMsg("\n found %zu at %zu\n", i, clus);
            clusters[i] = clus;
        }
    }
    primary_tree.updateTree();
    primary_tree.printTreeJson(stdout);
    if (args.topology_out_given)
    {
        string outFile = args.topology_out_arg;
        FILE *tFile = fopen((outFile + "tree.txt").c_str(), "w");
        primary_tree.printTree(tFile, insertionList, args.topology_out_arg);
    }
    return clusters;
}

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

    fprintf(stderr, "Loaded %zu seqs...signatureSize %zu\n", seqCount, signatureSize);
    clusters = clusterSignaturesPrim(seqs);
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
