#ifndef INCLUDE_CLUSTER_HPP
#define INCLUDE_CLUSTER_HPP
#include <omp.h>
using namespace std;

bool random_ = false;
size_t cap = 0;
size_t iteration = 0;
bool iteration_given = false;
bool debug_ = false;
bool force_split_ = false;
bool skip_ = false;
unsigned seed;

struct tree_meta_st
{
    bool readTree_ = false;
    bool writeTree_ = false;
    const char *inputTreePath;
    const char *outputTreePath;
};

tree_meta_st tree_meta;

template <typename cmdline_type>
string setArgs(cmdline_type args)
{
    split_threshold = 0.5;
    stay_threshold = 0.8;
    minimiser_match_threshold = 4;

    //?
    if (!args.input_given)
    {
        cout << "No input and/or query given! Exiting...\n";
        exit(1);
    }

    string inputFile = args.input_arg;

    if (args.random_arg)
    {
        random_ = true;
        seed = args.seed_arg;
    }

    if (args.iteration_given)
    {
        iteration_given = true;
        iteration = args.iteration_arg;
    }

    if (args.stay_threshold_given)
    {
        stay_threshold = args.stay_threshold_arg;
    }

    if (args.split_threshold_given)
    {
        split_threshold = args.split_threshold_arg;
    }
    else if (args.skip_arg)
    {
        skip_ = true;
        fprintf(stderr, "Skip inital tree construction, thresholds have no use\n");
        // if (iteration == 0)
        // {
        //     iteration = 1;
        //     fprintf(stderr, "No iteration given, default = 1\n");
        // }
        // fprintf(stderr, "Iteration = %zu\n", iteration);
    }
    else if (args.single_arg)
    {
        if (args.multiple_arg)
        {
            split_threshold = getSplitThresholdListSingle(inputFile);
        }
        else
        {
            split_threshold = getSplitThresholdSingle(inputFile);
        }
    }
    else
    {
        split_threshold = getSplitThreshold(inputFile, args.multiple_arg);
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

    tree_meta.readTree_ = args.topology_in_given;
    tree_meta.writeTree_ = args.topology_out_given;
    tree_meta.inputTreePath = args.topology_in_arg;
    tree_meta.outputTreePath = args.topology_out_arg;

    size_t firstindex = inputFile.find_last_of("/") + 1;
    size_t lastindex = inputFile.find_last_of(".");
    string rawname = inputFile.substr(firstindex, lastindex - firstindex);

    auto fileName = rawname + "-s" + to_string((int)(stay_threshold * 100)) + "-l" + to_string((int)(split_threshold * 100));
    if (args.tag_given)
    {
        fileName = args.tag_arg + fileName;
    }

    return (fileName + ".txt");
}

template <typename tree_type>
string readTreeLine(string s, string folder, tree_type &tree)
{
    int parent_temp = 0;
    size_t parent = 0;
    double priority = 0;
    size_t child = 0;
    string delimiter = ">";

    size_t pos = s.find(delimiter);
    string node = s.substr(0, pos);
    // sscanf(node.c_str(), "%zu", &parent);
    sscanf(node.c_str(), "%d(%lf)", &parent_temp, &priority);
    parent = parent_temp;
    if (parent_temp < 0)
    {
        parent = parent_temp * -1;
        tree.isRootNode[parent] = 1;
    }
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

        // vector<cell_type> signature;
        // readSignatures((folder + childStr + ".bin"), signature);
        // tree.readNode(parent, child, &signature[0]);
        // auto signature = tree.readInput((folder + childStr + ".bin").c_str());
        // tree.readNode(parent, child, signature, priority);
        tree.readNode(parent, child, (folder + childStr + ".bin").c_str(), priority);
    }
    // tree.updatePriority(parent);

    // return parent;
    return node;
}

// output signature size, if file too big => return empty "sigs", use readSignaturesBatch
template <typename tree_type>
size_t readTree(const string folder, tree_type &tree)
{
    string line;

    // Read from the text file
    ifstream listStream((folder + "tree.txt").c_str());

    // read signatureSize
    getline(listStream, line);

    // // use first line to set tree params
    // tree_type tree(partree_capacity);
    // getline(listStream, line);
    // sscanf(line.c_str(), "%zu", &signatureSize);
    // signatureWidth = signatureSize * sizeof(cell_type);
    // tree.means.resize(partree_capacity * signatureSize);

    // read last nodeIdx
    getline(listStream, line);
    size_t offset = 0;
    sscanf(line.c_str(), "%zu", &offset);

    while (getline(listStream, line))
    {
        readTreeLine(line, folder, tree);
    }

    // Close the file
    listStream.close();
    // tree.updateTree();
    // tree.printTreeJson(stdout);
    // return tree;

    return offset;
}

void outputClusters(FILE *pFile, const vector<size_t> &clusters)
{
    for (size_t sig = 0; sig < clusters.size(); sig++)
    {
        fprintf(pFile, "%llu,%llu\n", static_cast<unsigned long long>(sig), static_cast<unsigned long long>(clusters[sig]));
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

size_t getSingleton(vector<size_t> &clusters)
{
    unordered_map<size_t, size_t> remap;
    for (size_t clus : clusters)
    {
        if (remap.count(clus))
        {
            remap[clus]++;
        }
        else
        {
            remap[clus] = 1;
        }
    }

    // for (auto &p : remap)
    //     cerr << " " << p.first << " => " << p.second << '\n';

    auto p = min_element(remap.begin(), remap.end(),
                         [](const auto &l, const auto &r)
                         { return l.second < r.second; });

    // cerr << "min " << p->second << '\n';
    // fprintf(stderr, "Output %zu clusters\n", remap.size());

    return p->second;
}

seq_type getSeq(const vector<seq_type> &seqs, size_t i)
{
    return seqs[i];
}

const cell_type *getSeq(const vector<cell_type> &seqs, size_t i)
{
    return &seqs[i];
}

template <typename tree_type, typename signature_type>
vector<size_t> clusterSignatures(const vector<signature_type> &seqs, size_t seqCount, size_t mul = 1)
{
    // size_t seqCount = seqs.size() / signatureSize;
    vector<size_t> clusters(cap);

    size_t firstNodes = 1;
    if (firstNodes > seqCount)
        firstNodes = seqCount;

    size_t offset = 0;
    default_random_engine rng;

    tree_type tree(partree_capacity);
    tree.means.resize(partree_capacity * mul);

    if (tree_meta.readTree_)
    {
        offset = readTree(tree_meta.inputTreePath, tree);
        tree.printTreeJson(stderr);
    }
    offset += firstNodes;
    vector<size_t> insertionList(partree_capacity - offset); // potential nodes idx except root; root is always 0
    vector<size_t> foo(seqCount);

// node 0 reserved for root, node 1 reserved for leaves idx
#pragma omp parallel
    {
#pragma omp for
        for (size_t i = 0; i < partree_capacity - offset; i++)
        {
            // insertionList.push_back(partree_capacity - i);
            insertionList[i] = partree_capacity - i - 1;
        }

#pragma omp for
        for (int i = 0; i < seqCount; i++)
        {
            foo[i] = i;
        }
    }

    if (random_)
    {
        // unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        shuffle(foo.begin(), foo.end(), default_random_engine(seed));
    }
    foo.resize(cap);

    if (skip_)
    {
#pragma omp parallel for
        for (size_t i = 0; i < cap; i++)
        {
            size_t clus = tree.reinsert(getSeq(seqs, i * mul), foo[i]);
            printMsg("\n Reinsert %zu at %zu\n", foo[i], clus);
            // clusters[foo[i]] = tree.findAncestor(clus);
            clusters[foo[i]] = clus;
        }
        // tree.removeAmbi();
        // tree.printTreeJson(stderr);
    }
    else
    {
        try
        {

            if (force_split_)
            {
                for (size_t i = 0; i < cap; i++)
                {
                    printMsg("inserting %zu\n", foo[i]);
                    size_t clus = tree.insertSplitRoot(getSeq(seqs, i * mul), insertionList, foo[i]);
                    clusters[foo[i]] = clus;
                }
            }
            else
            {
                for (size_t i = 0; i < cap; i++)
                {
                    printMsg("inserting %zu\n", foo[i]);
                    size_t clus = tree.insert(getSeq(seqs, i * mul), insertionList, foo[i]);
                    clusters[foo[i]] = clus;
                }
            }
        }
        catch (const std::exception &e) // caught by reference to base
        {
            std::cout << " a standard exception was caught, with message '"
                      << e.what() << "'\n";
            tree.printTreeJson(stdout);
            
            // Recursively destroy all locks
            tree.destroyLocks();
            return clusters;
        }
        // for debugging
        if (iteration_given)
        {
            // prep to remove and reinsert ambi
            // fprintf(stderr, "\n\n\nBefore\n");
            // tree.printTreeJson(stderr);

            singleton = 0;
            // tree.trim();
            tree.removeAmbi();
            singleton = 1;
            printMsg("\n\nReinserting ambi (all)\n");
            tree.prepReinsert();

            if (tree_meta.writeTree_)
            {
                // string outFile = tree_meta.outputTreePath;
                // FILE *tFile = fopen((outFile + "tree.txt").c_str(), "w");
                tree.printTree(tree_meta.outputTreePath, insertionList);
            }
            // #pragma omp parallel for
            for (size_t i = 0; i < cap; i++)
            {
                size_t clus = tree.reinsert(getSeq(seqs, i * mul), foo[i]);
                printMsg("\n Reinsert %zu at %zu\n", foo[i], clus);
                // clusters[foo[i]] = tree.findAncestor(clus);
                clusters[foo[i]] = clus;
            }

            // tree.printTreeJson(stderr);
        }
        // else if (tree_meta.readTree_)
        // {
        //     singleton = 0;
        //     tree.removeAmbi();
        // }
    }

    singleton = 1;

    for (size_t run = 0; run < iteration; run++)
    {
        if (debug_)
        {

            // singleton = getSingleton(clusters) + 1;
            singleton++;
        }
        fprintf(stderr, "Iteration %zu (singleton = %zu)\n", run, singleton);

        // tree.trim();

        tree.removeAmbi();
        tree.prepReinsert();
// singleton++;

// tree.printTreeJson(stderr);
#pragma omp parallel for
        for (size_t i = 0; i < cap; i++)
        {
            size_t clus = tree.reinsert(getSeq(seqs, i * mul), foo[i]);

            printMsg("\n found %zu at %zu\n", foo[i], clus);
            // clusters[foo[i]] = tree.findAncestor(clus);
            clusters[foo[i]] = clus;
        }

        // if (debug_)
        // {
        //     auto fileName = "nodeDistance-r" + to_string((size_t)(run)) + ".txt";
        //     FILE *nFile = fopen(fileName.c_str(), "w");
        //     tree.printNodeDistance(nFile, seqs, clusters);

        //     fileName = "clusters-r" + to_string((size_t)(run)) + ".txt";
        //     FILE *cFile = fopen(fileName.c_str(), "w");
        //     outputClusters(cFile, clusters);
        // }
    }
    tree.updateTree();

    // FILE *pFile = fopen("nodeDistance.txt", "w");
    // tree.printNodeDistance(pFile, seqs, clusters);

    // tree.printTreeJson(stdout);
    // FILE *hFile = fopen("hierarchy.txt", "w");
    // fprintf(hFile, "parent,child,rank\n");
    // tree.outputHierarchy(hFile);

    // if (tree_meta.writeTree_)
    // {
    //     string outFile = tree_meta.outputTreePath;
    //     FILE *tFile = fopen((outFile + "tree.txt").c_str(), "w");
    //     tree.printTree(tFile, insertionList, tree_meta.outputTreePath);
    // }

    if (tree_meta.writeTree_)
    {
        tree.printTreeJson(stdout);
    }

    // Recursively destroy all locks
    tree.destroyLocks();
    return clusters;
}

#endif