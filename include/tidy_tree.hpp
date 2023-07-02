
// There are two kinds of ktree nodes- branch nodes and leaf nodes
// Both contain a signature matrix, plus their own signature
// (the root node signature does not matter and can be blank)
// Branch nodes then contain 'order' links to other nodes
// Leaf nodes do not.
// However, as leaf nodes may become branch nodes, we allocate
// the space anyway.
// As the space to be used is determined at runtime, we use
// parallel arrays, not structs

#ifndef INCLUDE_TIDY_TREE_HPP
#define INCLUDE_TIDY_TREE_HPP
#include <omp.h>
#include <bitset>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include "read.hpp"
#include "distance.hpp"
#include <stdarg.h>
#include <stdio.h>

using namespace std;
bloom_parameters parameters;
size_t partree_capacity = 10000;
size_t singleton = 1;
size_t minClusSize = partree_capacity;
size_t tree_order = 5;
bool print_ = false;

struct tt_data
{
    std::bitset<8> statuses;
    size_t dest = 0;
    double max_similarity = 0;
    size_t dest_branch = 0;
    double max_branch_similarity = 0;
    size_t dest_super = 0;
    double max_super_similarity = 0;
    size_t dest_root = 0;
    double max_root_similarity = 0;
    size_t mismatch = 0;

    vector<size_t> stay_branch;
    vector<size_t> stay_leaf;
    vector<size_t> nn_leaf;
    vector<size_t> nn_branch;
    vector<size_t> stay_super;
    vector<size_t> stay_root;
};

template <typename T>
void removeVecValue(vector<T> &vec, T value)
{
    vec.erase(std::remove(vec.begin(), vec.end(), value), vec.end());
}

template <typename T>
void removeVecIdx(vector<T> &vec, size_t idx)
{
    vec.erase(vec.begin() + idx);
}

template <typename T>
void insertVecRange(vector<T> &to, vector<T> &from)
{
    to.insert(to.end(), from.begin(), from.end());
}

void printMsg(const char *format, ...)
{
    if (print_)
    {
        va_list args;
        va_start(args, format);
        vfprintf(stderr, format, args);
        va_end(args);
    }
}

template <typename T>
T returnEmpy()
{
    T empty;
    return empty;
}

// Derived class
template <typename s_type, typename const_s_type, typename sVec_type>
class tidy_tree
{
public:
    size_t root = 0;                   // # of root node
    vector<size_t> childCounts;        // n entries, number of children
    vector<int> isBranchNode;          // n entries, is this a branch node
    vector<int> isRootNode;            // n entries, is this root of subtree
    vector<int> isSuperNode;           // n entries, is this root of subtree
    vector<vector<size_t>> childLinks; // n * o entries, links to children
    vector<int> isAmbiNode;            // n entries, is this a branch node
    vector<vector<size_t>> ambiLinks;  // n * o entries, links to children
    vector<vector<size_t>> seqIDs;     // n * o entries, links to children
    vector<size_t> parentLinks;        // n entries, links to parents
    vector<distance_type> priority;    // n entries, links to parents
    // vector<vector<cell_type>> means;   // n * signatureSize entries, node signatures
    sVec_type means;            // n * signatureSize entries, node signatures
    vector<sVec_type> matrices; // capacity * signatureSize * n
    vector<omp_lock_t> locks;   // n locks
    size_t capacity = 0;        // Set during construction, currently can't change

    virtual s_type getMeanSig(size_t node) { return returnEmpy<s_type>(); }

    virtual void updateMeanSig(size_t node, const_s_type signature) {}

    virtual void addSigToMatrix(size_t node, const_s_type signature) {}

    virtual double calcSimilaritySigToNode(size_t node, sVec_type signatures, size_t i) { return 0; }

    virtual double calcSimilarityWrap(const_s_type a, const_s_type b, size_t signatureSize_ = 0) { return 0; }

    virtual sVec_type createRandomSigs(size_t node, size_t clusterCount, size_t s = 0) { return returnEmpy<sVec_type>(); }

    virtual inline sVec_type getNonAmbiMatrix(size_t node) { return returnEmpy<sVec_type>(); }

    virtual sVec_type unionNodeMean(size_t node) { return returnEmpy<sVec_type>(); }

    // union mean of children
    virtual inline void updateNodeMean(size_t node) {}

    virtual inline void updateMatrixIdx(size_t parent, size_t idx, size_t node) {}

    virtual inline double calcAvgSim(size_t node) { return 0; }

    virtual void clearMean(size_t node) {}

    virtual double preloadPriority(size_t node, const_s_type signature) { return 0; }

    void reserve(size_t capacity)
    {
        // For safety, only call this at startup currently
        if (this->capacity != 0)
        {
            fprintf(stderr, "Reserve can only be called from 0 capacity\n");
            exit(1);
        }
        this->capacity = capacity;

        //#pragma omp parallel
        {
            //#pragma omp single
            {
                childCounts.resize(capacity);
            }
            //#pragma omp single
            {
                isBranchNode.resize(capacity);
            }
            //#pragma omp single
            {
                isAmbiNode.resize(capacity);
                ambiLinks.resize(capacity);
                isRootNode.reserve(capacity);
                isSuperNode.reserve(capacity);
            }
            //#pragma omp single
            {
                childLinks.resize(capacity);
            }
            //#pragma omp single
            {
                seqIDs.resize(capacity);
            }
            //#pragma omp single
            {
                priority.resize(capacity);
            }
            //#pragma omp single
            {
                parentLinks.resize(capacity);
            }
            //#pragma omp single
            {
                locks.resize(capacity);
            }
            //#pragma omp single
            {
                matrices.resize(capacity);
            }
            //#pragma omp single
            {
                // means.resize(capacity * signatureSize);
                means.resize(capacity);
            }
        }
    }

    tidy_tree(size_t capacity)
    {
        reserve(capacity);
        childCounts[root] = 0;
        isBranchNode[root] = 0;
        isRootNode[root] = 1;
        isSuperNode[root] = 0;
    }

    void printMatrix(FILE *stream, size_t node)
    {
        fprintf(stream, ">>>printing matrix at node %zu\n", node);
        for (s_type seq : matrices[node])
        {
            dbgPrintSignature(stream, seq);
        }
    }

    size_t getNewNodeIdx(vector<size_t> &insertionList)
    {
        if (insertionList.empty())
        {
            fprintf(stderr, "ERROR: ran out of insertion points\n");
            exit(1);
        }
        size_t idx = insertionList.back();
        insertionList.pop_back();

        // Initialise lock
        omp_init_lock(&locks[idx]);
        return idx;
    }

    size_t findAncestor(size_t node)
    {
        while (parentLinks[node] != root)
        {
            node = parentLinks[node];
        }
        return node;
    }

    tuple<size_t, size_t> findAncestorNlevel(size_t node)
    {
        size_t level = 0;
        while (parentLinks[node] != root)
        {
            node = parentLinks[node];
            level++;
        }
        return make_tuple(node, level);
    }

    size_t findLevel(size_t node)
    {
        size_t level = 0;
        while (node != root)
        {
            node = parentLinks[node];
            level++;
        }
        return level;
    }

    void printNodeJson(FILE *stream, size_t tnode)
    {
        fprintf(stream, "{\"node\":\"%zu\",", tnode);
        fprintf(stream, "\"branch\":\"%d\",", isBranchNode[tnode]);
        fprintf(stream, "\"ambi\":\"%d\",", isAmbiNode[tnode]);
        fprintf(stream, "\"subtree\":\"%d\",", isRootNode[tnode]);
        fprintf(stream, "\"super\":\"%d\",", isSuperNode[tnode]);
        fprintf(stream, "\"priority\":\"%.2f\",", priority[tnode]);
        fprintf(stream, "\"childCount\":\"%zu\",\"content\":\"*", seqIDs[tnode].size());
        // fprintf(stream, "{\"node\":\"%zu\",\"branch\":\"%zu\",\"priority\":\"%.2f\",\"childCount\":\"%zu\",\"content\":\"*", tnode, isBranchNode[tnode], calcNodeMaxsimilarity(tnode), seqIDs[tnode].size());

        for (size_t seq : seqIDs[tnode])
        {
            fprintf(stream, "%zu,", seq);
        }
        // fprintf(stream, "\",\"children\":[");
    }

    void printSubTreeJson(FILE *stream, size_t tnode)
    {
        if (childCounts[tnode] > 0)
        {
            printNodeJson(stream, tnode);
            fprintf(stream, "\",\"children\":[");

            printSubTreeJson(stream, childLinks[tnode][0]);

            for (size_t i = 1; i < childLinks[tnode].size(); i++)
            {
                fprintf(stream, ",");
                printSubTreeJson(stream, childLinks[tnode][i]);
            }
            fprintf(stream, "]}");
        }
        else
        {
            // fprintf(stream, "{\"node\":\"%zu\",\"priority\":\"%.2f\",\"childCount\":\"%zu\",\"content\":\"*", tnode, priority[tnode], seqIDs[tnode].size());
            // // fprintf(stream, "{\"node\":\"%zu\",\"priority\":\"%.2f\",\"childCount\":\"%zu\",\"content\":\"*", tnode, calcNodeMaxsimilarity(tnode), seqIDs[tnode].size());

            // for (size_t seq : seqIDs[tnode])
            // {
            //     fprintf(stream, "%zu,", seq);
            // }
            printNodeJson(stream, tnode);

            fprintf(stream, "\"}");
        }
    }

    void calcNodeDistance(FILE *stream, size_t tnode, size_t idx, const_s_type sig)
    {
        double distance = calcHD(getMeanSig(tnode), sig);
        fprintf(stream, "%zu,%zu,%.4f", tnode, idx, distance);
    }

    void printNodeDistance(FILE *stream, sVec_type seqs, vector<size_t> &clusters) {}

    void printTreeJson(FILE *stream, size_t node = 0)
    {
        fprintf(stream, "var treeData = ");
        printSubTreeJson(stream, node);
        fprintf(stream, ";\n");
    }

    virtual void printSignature(ostream &wf, size_t node) {}

    void printNodeMean(string outfolder, size_t node)
    {
        ofstream wf(outfolder + to_string(node) + ".bin", ios::out | ios::binary);
        writeInt(wf, signatureSize);
        // wf.write(reinterpret_cast<const char *>(&i), sizeof(cell_type));
        printSignature(wf, node);
        wf.close();
    }

    void printTree(FILE *stream, string outfolder, size_t node)
    {
        fprintf(stream, "%zu>", node);
        for (size_t child : childLinks[node])
        {
            fprintf(stream, "%zu,", child);
            printNodeMean(outfolder, child);
        }
        fprintf(stream, "\n");
        for (size_t child : childLinks[node])
        {
            if (childCounts[child] > 0)
            {
                printTree(stream, outfolder, child);
            }
        }
    }

    void printTree(FILE *stream, vector<size_t> &insertionList, string outfolder = "")
    {
        fprintf(stream, "%zu\n%zu\n", signatureSize, getNewNodeIdx(insertionList));
        printTree(stream, outfolder, 0);
    }

    void readNode(size_t parent, size_t child, const_s_type signature)
    {
        isBranchNode[parent] = 1;
        parentLinks[child] = parent;
        childLinks[parent].push_back(child);
        childCounts[parent]++;
        updateMeanSig(child, signature);
        addSigToMatrix(parent, signature);
    }

    int getNodeIdx(size_t node)
    {
        size_t parent = parentLinks[node];
        int idx = -1;
        for (size_t i = 0; i < childLinks[parent].size(); i++)
        {
            if (childLinks[parent][i] == node)
            {
                idx = i;
                break;
            }
        }
        return idx;
    }

    s_type createRandomSig(vector<size_t> children, size_t s = 0)
    {
        // unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        unsigned seed = s;
        if (seed == 0)
        {
            seed = chrono::system_clock::now().time_since_epoch().count();
        }
        default_random_engine rng(seed);
        uniform_int_distribution<size_t> dist(0, children.size() - 1);

        const_s_type randomSig; // = getMeanSig(children[dist(rng)])
        // find the smallest windows count
        size_t winNum = getMeanSig(children[0]).size();
        for (size_t child : children)
        {
            s_type matrix = getMeanSig(child);
            if (matrix.size() < winNum)
            {
                winNum = matrix.size();
            }
        }

        randomSig.resize(winNum);

        vector<vector<size_t>> counters;
        for (size_t w = 0; w < winNum; w++)
        {
            size_t s = dist(rng);
            randomSig[w] = getMeanSig(children[s])[w];
        }

        return randomSig;
    }

    // delete a child from its parent, need to format the child separately
    inline void deleteNode(size_t node)
    {
        // printMsg("deleting %zu\n", node);
        size_t parent = parentLinks[node];

        int idx = getNodeIdx(node);

        if (idx == -1)
        {
            printMsg("ERROR deleting %zu from %zu!!\n", node, parent);
        }
        else
        {

            removeVecIdx(matrices[parent], idx);
            removeVecIdx(childLinks[parent], idx);
            childCounts[parent]--;
            // parentLinks[node] = 0;
            // isBranchNode[node] = 0;
        }
    }

    void clearNode(size_t node)
    {
        isBranchNode[node] = 0;
        childCounts[node] = 0;
        childLinks[node].clear();
        ambiLinks[node].clear();
        seqIDs[node].clear();
        matrices[node].clear();
        clearMean(node);
        priority[node] = 0;
        parentLinks[node] = 0;
    }

    void deleteUnitig(size_t node)
    {
        size_t parent = parentLinks[node];
        size_t child = childLinks[node][0];
        size_t idx = getNodeIdx(node);

        childLinks[parent][idx] = child;
        parentLinks[child] = parent;
        // deleteNode(node);
        clearNode(node);
    }

    inline bool checkdeleteUniSuper(size_t node)
    {
        if (!isRootNode[node] && childCounts[node] == 1)
        {
            deleteUnitig(node);
            return true;
        }
        return false;
    }

    inline void updateParentMean(size_t node)
    {
        while (node != root)
        {
            int idx = getNodeIdx(node);
            updateNodeMean(node);
            // updatePriority(node);
            size_t parent = parentLinks[node];

            if (idx == -1)
            {
                printMsg("ERROR updating mean %zu from %zu!!\n", node, parent);
            }
            else
            {
                updateMatrixIdx(parent, idx, node);
            }
            node = parent;
        }
    }

    // union mean of children
    inline void updatePriority(size_t node)
    {
        priority[node] = calcAvgSim(node);
    }

    inline size_t stayNode(const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t node)
    {
        if (isBranchNode[node])
        {
            return tt_root(signature, insertionList, idx, node);
        }
        else
        {
            seqIDs[node].push_back(idx);
            addSigToMatrix(node, signature);
            // updatePriority(node);
            return node;
        }
    }

    // check against all ambiNode
    // create "dest" (non-ambi) and move all matrices that stay with input sig from the ambi nodes
    // no nothing stay, delete dest AND
    // create new ambi if there is at least one split in all ambi nodes
    // else store in the first ambiNode
    inline size_t stayAmbi(const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t node)
    {
        return stayNode(signature, insertionList, idx, ambiLinks[node][0]);
    }

    void moveParent(size_t child, size_t new_parent, bool d = true)
    {
        childCounts[new_parent]++;
        childLinks[new_parent].push_back(child);
        addSigToMatrix(new_parent, getMeanSig(child));
        // printMsg("moving %zu to %zu \n", child, new_parent);

        if (d)
        {
            deleteNode(child);
        }
        parentLinks[child] = new_parent;
        // printMsg("deleted %zu \n", child);
    }

    // haven't update childLinks[t_parent] yet
    inline size_t createParent(size_t node, vector<size_t> &insertionList)
    {
        size_t t_parent = getNewNodeIdx(insertionList);
        isBranchNode[t_parent] = 1;
        parentLinks[t_parent] = node;
        isRootNode[t_parent] = 0;
        isSuperNode[t_parent] = 0;

        // update node
        childCounts[node]++;
        childLinks[node].push_back(t_parent);

        printMsg("create t_parent %zu\n", t_parent);
        return t_parent;
    }

    inline size_t createNode(const_s_type signature, vector<size_t> &insertionList, size_t t_parent, size_t idx)
    {
        size_t new_node = getNewNodeIdx(insertionList);
        parentLinks[new_node] = t_parent;
        addSigToMatrix(new_node, signature);
        seqIDs[new_node].push_back(idx);
        updateMeanSig(new_node, signature);
        isRootNode[new_node] = 0;
        isSuperNode[new_node] = 0;

        // then add new node to t_parent
        childCounts[t_parent]++;
        childLinks[t_parent].push_back(new_node);
        addSigToMatrix(t_parent, getMeanSig(new_node));

        if (isBranchNode[t_parent])
        {
            printMsg("--------------b %zu, %.2f\n", new_node, t_parent);
            updatePriority(t_parent);
        }

        printMsg("create new node %zu\n", new_node);
        return new_node;
    }

    // Ambi Node cannot be the first child or else update node mean will have problem
    inline size_t createAmbiNode(const_s_type signature, vector<size_t> &insertionList, size_t node, size_t idx)
    {
        size_t dest = createNode(signature, insertionList, node, idx);

        printMsg(">>Ambi\n");
        ambiLinks[node].push_back(dest);
        isAmbiNode[dest] = 1;
        updatePriority(node);
        return dest;
    }

    inline pair<size_t, double> similarityStatusF(const_s_type sig1, const_s_type sig2, double local_stay_t, double offset = 1, bool isBranch = false)
    {
        //?
        if (local_stay_t > stay_threshold)
        {
            local_stay_t = stay_threshold;
        }

        double local_split_t = split_threshold * offset;
        double similarity = calcSimilarityWrap(sig1, sig2);

        printMsg("%.2f (%.2f)", similarity, local_split_t);

        if (similarity >= local_stay_t)
        {
            printMsg(": STAY>\n");
            // return STAY_F;
            return make_pair(STAY_F, similarity);
        }
        else if (similarity <= local_split_t)
        {
            printMsg(": SPLIT>\n");
            // return SPLIT_F;
            return make_pair(SPLIT_F, similarity);
        }
        else
        {
            printMsg(": NN>\n");
            // return NN_F;
            if (isBranch)
            {
                return make_pair(NN_BRANCH_F, similarity);
            }
            else
            {
                return make_pair(NN_LEAVE_F, similarity);
            }
        }
    }

    inline pair<size_t, double> similarityStatus(size_t child, const_s_type signature)
    {
        printMsg("<%zu, ", child);

        if (isSuperNode[child])
        {
            // super can only be stay or split
            printMsg("(super %.2f)", priority[child]);
            return similarityStatusF(getMeanSig(child), signature, priority[child], 100);
        }
        else if (isBranchNode[child])
        {
            // super can only be stay or split

            if (priority[child] == 0)
            {
                // branch only 1 child
                printMsg("(branch)");
                return similarityStatusF(getMeanSig(child), signature, stay_threshold, 1, true);
            }
            else
            {
                printMsg("(branch %.2f)", priority[child]);
                return similarityStatusF(getMeanSig(child), signature, priority[child], 1, true);
            }
        }
        else
        {
            return similarityStatusF(getMeanSig(child), signature, stay_threshold);
        }
    }

    inline void mergeBranch(size_t branch)
    {
        size_t parent = parentLinks[branch];
        vector<size_t> children = childLinks[branch];
        clearNode(branch);
        parentLinks[branch] = parent;

        printMsg("merging branch %zu\n", branch);
        for (size_t child : children)
        {
            if (!isAmbiNode[child])
            {
                insertVecRange(matrices[branch], matrices[child]);
                insertVecRange(seqIDs[branch], seqIDs[child]);
            }
            clearNode(child);
        }

        updateParentMean(branch);
    }

    // // use it when all children are singletons
    // void dissolveBranch(size_t branch)
    // {
    //     matrices[branch].clear();
    //     for (size_t child:childLinks[branch]){
    //         insertVecRange(matrices[branch], )
    //     }
    // }

    void mergeAmbi(size_t branch)
    {
        vector<size_t> candidates;
        for (size_t ambi : ambiLinks[branch])
        {
            if (isSingleton(ambi))
            {
                candidates.push_back(ambi);
            }
        }
        if (candidates.size() > 1)
        {
            size_t dest = candidates.back();
            candidates.pop_back();
            printMsg("Merging ambis: ");
            for (size_t s_ambi : candidates)
            {
                printMsg("%zu,", s_ambi);
                insertVecRange(matrices[dest], matrices[s_ambi]);
                insertVecRange(seqIDs[dest], seqIDs[s_ambi]);
                removeVecValue(ambiLinks[branch], s_ambi);
                deleteNode(s_ambi);
                clearNode(s_ambi);
            }
            printMsg("\n To: %zu\n", dest);
            updateNodeMean(dest);
        }
    }

    inline size_t insertAmbi(const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t branch)
    {

        if (ambiLinks[branch].size() == 0)
        {
            return createAmbiNode(signature, insertionList, branch, idx);
        }

        printMsg("InsertAmbi %zu\n", branch);

        size_t dest = 0;
        size_t ambi_split = 0;
        for (size_t ambi : ambiLinks[branch])
        {
            // // preload sig into matrices to check for priority
            // addSigToMatrix(ambi, signature);
            // updatePriority(ambi);
            // // release preloaded matrix
            // matrices[ambi].pop_back();

            double preload_priority = preloadPriority(ambi, signature);

            printMsg(">>ambi %zu, %.2f\n", ambi, preload_priority);

            if (preload_priority >= stay_threshold)
            {
                dest = stayNode(signature, insertionList, idx, ambi);
                isAmbiNode[dest] = 0;
                updateParentMean(dest);
                printMsg("@@ambi %zu, %.2f\n", dest, priority[ambi]);

                //? also remove from ambi link;
                // ambiLinks[branch].clear();
                removeVecValue(ambiLinks[branch], ambi);
                return dest;
            }
            else if (preload_priority < split_threshold)
            {
                ambi_split++;
            }
            else
            {
                dest = ambi;
            }
        }

        if (ambi_split == ambiLinks[branch].size())
        {
            return createAmbiNode(signature, insertionList, branch, idx);
        }
        else
        {
            return stayNode(signature, insertionList, idx, dest);
        }
    }

    // create ambi and insert to it if there is none
    // else, add to suitable ambi and turn the ambi to leaf if the prioroty >stay threshold
    // else create new ambi
    inline size_t insertBranch(const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t branch)
    {
        double max_similarity = stay_threshold;
        size_t dest = 0;
        for (size_t child : childLinks[branch])
        {
            double similarity = calcSimilarityWrap(getMeanSig(child), signature);
            printMsg("b %zu, %.2f\n", child, similarity);
            if (similarity >= max_similarity)
            {
                max_similarity = similarity;
                dest = child;
            }
        }

        if (dest != 0)
        {
            dest = stayNode(signature, insertionList, idx, dest);
            if (isAmbiNode[dest])
            {
                isAmbiNode[dest] = 0;
                updateParentMean(dest);
                removeVecValue(ambiLinks[branch], dest);
            }
            return dest;
        }

        return insertAmbi(signature, insertionList, idx, branch);
    }

    // find best leaf from the NN branches, if dest = 0, create new node
    inline size_t doNNBranch(const_s_type signature, size_t node, tt_data dt, vector<size_t> &insertionList, size_t idx)
    {
        // xxx, should be ok
        size_t t_parent = node;
        vector<size_t> target = dt.nn_branch;
        if (checkCreateSuper(node, true, dt.mismatch, target.size()))
        {
            t_parent = createSuper(node, insertionList, target);
            updateParentMean(t_parent);
        }

        printMsg("doNNBranch %zu\n", t_parent);

        double max_similarity = stay_threshold;
        size_t dest = 0;
        double max_ambi_similarity = split_threshold;
        size_t dest_ambi = 0;
        for (size_t branch : target)
        {
            for (size_t child : childLinks[branch])
            {
                double similarity = calcSimilarityWrap(getMeanSig(child), signature);
                if (isAmbiNode[child])
                {
                    printMsg("b_ambi %zu, %.2f\n", child, similarity);
                    if (similarity >= max_ambi_similarity)
                    {
                        max_ambi_similarity = similarity;
                        dest_ambi = child;
                    }
                }
                else
                {
                    printMsg("b %zu, %.2f\n", child, similarity);
                    if (similarity >= max_similarity)
                    {
                        max_similarity = similarity;
                        dest = child;
                    }
                }
            }
        }

        if (dest + dest_ambi == 0)
        {
            dest = createNode(signature, insertionList, t_parent, idx);
        }
        else if (dest != 0)
        {
            dest = stayNode(signature, insertionList, idx, dest);
        }
        else
        {
            dest = stayNode(signature, insertionList, idx, dest_ambi);
            updateNodeMean(dest_ambi);
            if (priority[dest_ambi] >= stay_threshold)
            {
                printMsg("@@ambi %zu, %.2f\n", dest_ambi, priority[dest_ambi]);
                isAmbiNode[dest_ambi] = 0;
                // updateParentMean(dest_ambi);// will do below
                removeVecValue(ambiLinks[parentLinks[dest_ambi]], dest_ambi);
            }
        }

        updateParentMean(dest);
        return dest;
    }

    // rearrange grandchildren to their best parent
    inline size_t recluster(size_t node)
    {
        printMsg("reclustering %zu\n", node);
        bool moved = false;
        for (size_t child : childLinks[node])
        {
            if (childCounts[child] - ambiLinks[child].size() <= 1)
            {
                continue;
            }
            // printMsg(">>>%zu\n",child);
            for (size_t grandchild : childLinks[child])
            {
                if (isAmbiNode[grandchild])
                {
                    continue;
                }
                double max_simimlarity = 0;
                size_t dest = child;
                for (size_t t_parent : childLinks[node])
                {
                    double similarity = calcSimilarityWrap(getMeanSig(t_parent), getMeanSig(grandchild));
                    // printMsg("%zu, %zu, %.2f\n", grandchild, t_parent, similarity);
                    if (similarity > max_simimlarity)
                    {
                        max_simimlarity = similarity;
                        dest = t_parent;
                    }
                }

                if (dest != child)
                {
                    printMsg("success %zu, %zu\n", grandchild, dest);
                    moveParent(grandchild, dest);
                    moved = true;
                    for (size_t ambi : ambiLinks[child])
                    {
                        moveParent(ambi, dest);
                    }
                }
            }
        }

        if (!moved)
        {
            for (size_t child : childLinks[node])
            {
                if (childCounts[child] == 1)
                {
                    deleteUnitig(child);
                }
            }
            updateNodeMean(node);
            return 0;
        }

        // tidy
        for (size_t child : childLinks[node])
        {
            if (childCounts[child] == 0)
            {
                deleteNode(child);
            }
            else if (childCounts[child] == 1)
            {
                deleteUnitig(child);
            }
            else
            {
                updateNodeMean(child);
            }
        }
        updateParentMean(node);
        // printMsg("success\n");
        return 1;
    }

    inline size_t createUniBranch(size_t node, vector<size_t> &insertionList, const_s_type signature, size_t idx)
    {
        size_t t_parent = createParent(node, insertionList);
        size_t dest = createNode(signature, insertionList, t_parent, idx);
        addSigToMatrix(node, signature);
        updateParentMean(t_parent);
        return dest;
    }

    // create a new branch under "node"
    // then move candidates to the new branch
    inline size_t createBranch(size_t node, vector<size_t> &insertionList, vector<size_t> candidates)
    {
        size_t t_parent = createParent(node, insertionList);

        for (size_t child : candidates)
        {
            moveParent(child, t_parent);
        }

        addSigToMatrix(node, getMeanSig(t_parent));
        updateParentMean(t_parent);
        return t_parent;
    }

    // create a new super under "node"
    // then move candidates to the new super
    inline size_t createSuper(size_t node, vector<size_t> &insertionList, vector<size_t> candidates)
    {
        size_t t_parent = createBranch(node, insertionList, candidates);
        isSuperNode[t_parent] = 1;
        return t_parent;
    }

    //? return if all children moved
    inline bool relocate(vector<size_t> parents, vector<size_t> children)
    {
        size_t moved = 0;
        printMsg("Relocating...\n");
        for (size_t child : children)
        {
            double max_simimlarity = split_threshold;
            size_t dest = 0;
            for (size_t parent : parents)
            {
                double similarity = calcSimilarityWrap(getMeanSig(child), getMeanSig(parent));
                printMsg("reloc %zu, %.2f\n", parent, similarity);
                if (similarity > max_simimlarity)
                {
                    max_simimlarity = similarity;
                    dest = parent;
                }
            }
            //? maybe find the best level to move
            if (dest != 0 && dest != parentLinks[child])
            {
                moveParent(child, dest);
                moved++;
            }
        }

        if (moved > 0)
        {
            for (size_t node : parents)
            {
                updateNodeMean(node);
            }
        }
        return moved == children.size();
    }

    // check if first child of branch is a singleton
    inline bool isSingleton(size_t child)
    {
        return matrices[child].size() == 1;
    }

    // separate roots and non roots children of node
    // first list contains roots
    // second list contains non-roots
    inline vector<vector<size_t>> separateRootChildren(size_t node)
    {
        vector<vector<size_t>> children(2);
        for (size_t child : childLinks[node])
        {
            if (isRootNode[child])
            {
                children[0].push_back(child);
            }
            else
            {
                children[1].push_back(child);
            }
        }
        return children;
    }

    inline tt_data getStatusRoot(const_s_type signature, vector<size_t> roots)
    {
        tt_data dt;
        for (size_t child : roots)
        {
            ///???
            double similarity = calcSimilarityWrap(getMeanSig(child), signature);
            printMsg("(root %zu, %.2f, %.2f)", child, priority[child], similarity);
            if (similarity > priority[child])
            {
                if (similarity > dt.max_root_similarity)
                {
                    dt.max_root_similarity = similarity;
                    dt.dest_root = child;
                }
                dt.statuses.set(6);
            }
        }
        return dt;
    }

    inline tt_data getStatus(const_s_type signature, size_t node)
    {
        tt_data dt;

        for (size_t child : childLinks[node])
        {
            // skip ambiNode
            if (isAmbiNode[child])
            {
                // printMsg("Ambi: %zu\n", child);
                continue;
            }

            double similarity = 0;
            size_t status = 0;

            if (isRootNode[child])
            { ///???
                similarity = calcSimilarityWrap(getMeanSig(child), signature);
                printMsg("(root %zu, %.2f, %.2f)\n", child, priority[child], similarity);
                if (similarity > priority[child])
                {
                    if (similarity > dt.max_root_similarity)
                    {
                        dt.max_root_similarity = similarity;
                        dt.dest_root = child;
                    }
                    dt.statuses.set(6);
                }
                continue;
            }
            std::tie(status, similarity) = similarityStatus(child, signature);
            switch (status)
            {
            case STAY_F:
                // stay.push_back(child);
                if (isSuperNode[child])
                {
                    if (priority[child] <= split_threshold)
                    {
                        printMsg("split super\n");
                        dissolveSuper(child);
                        dt.statuses.set(6);
                        return getStatus(signature, node);
                    }
                    dt.stay_super.push_back(child);
                    if (similarity > dt.max_super_similarity)
                    {
                        dt.max_super_similarity = similarity;
                        dt.dest_super = child;
                    }
                    dt.statuses.set(5);
                }
                else if (isBranchNode[child])
                {
                    dt.stay_branch.push_back(child);
                    if (similarity > dt.max_branch_similarity)
                    {
                        dt.max_branch_similarity = similarity;
                        dt.dest_branch = child;
                    }
                    dt.statuses.set(3);
                }
                else
                {
                    dt.stay_leaf.push_back(child);
                    if (similarity > dt.max_similarity)
                    {
                        dt.max_similarity = similarity;
                        dt.dest = child;
                    }
                    // printMsg("child %zu\n", dt.dest);
                    dt.statuses.set(1);
                }

                break;
            case SPLIT_F:
                // mismatch.push_back(child);
                dt.mismatch++;
                break;
            case NN_LEAVE_F:
                dt.nn_leaf.push_back(child);
                dt.statuses.set(2);
                break;
            default:
                dt.nn_branch.push_back(child);
                dt.statuses.set(4);
                break;
            }
        }
        return dt;
    }

    // return true if it will not create a unitigrelocateRemaining(node,
    inline bool checkCreateSuper(size_t node, bool moved_all, size_t mismatch, size_t target_size)
    {
        bool doRoot = false;
        bool doSuper = false;

        if (isRootNode[node])
        {
            if (target_size > 1)
            {
                doRoot = true;
            }
            else if (!moved_all)
            {
                doRoot = true;
            }
        }
        else
        {
            if (!moved_all && mismatch != 0)
            {
                doSuper = true;
            }
            else if (target_size > 1 && mismatch != 0)
            {
                doSuper = true;
            }
        }

        return doRoot || doSuper;
    }

    inline bool relocateRemaining(size_t node, bool moved_all_current, vector<size_t> candidates, vector<size_t> targets)
    {
        bool moved_all_candidates = moved_all_current;
        if (!moved_all_candidates)
        {
            // try to relocate the remaining stay leaves to super
            vector<size_t> remaining_candidates;
            for (size_t child : candidates)
            {
                if (parentLinks[child] == node)
                {
                    remaining_candidates.push_back(child);
                }
            }

            moved_all_candidates = relocate(targets, remaining_candidates);
        }
        return moved_all_candidates;
    }

    inline void movedRemaining(size_t node, size_t t_parent, bool moved_all, vector<size_t> candidates)
    {
        if (!moved_all)
        {
            for (size_t child : candidates)
            {
                if (parentLinks[child] == node)
                {
                    moveParent(child, t_parent);
                }
            }
        }
    }

    // relocate candidates to target, then createSuper accordingly
    // level = 3 = super [default]
    // level = 2 = nn branch
    // level = 1 = stay branch
    inline size_t doTarget_candidates(tt_data dt, size_t node, vector<size_t> candidates, vector<size_t> &insertionList, size_t level = 3)
    {
        vector<size_t> target;
        switch (level)
        {
        case 1:
            target = dt.stay_branch;
            break;
        case 2:
            target = dt.nn_branch;
            break;
        case 3:
            target = dt.stay_super;
            break;
        default:
            fprintf(stderr, "ERROR - doTarget_StayNN\n");
            return node;
            break;
        }
        bool moved_all_nn = relocate(target, candidates);

        if (!checkdeleteUniSuper(node))
        {
            if (checkCreateSuper(node, moved_all_nn, dt.mismatch, target.size()))
            { // there is some nn branches can't be moved to stay super, create new super
                size_t t_parent = createSuper(node, insertionList, target);
                movedRemaining(node, t_parent, moved_all_nn, candidates);
                updateParentMean(t_parent);
                return t_parent;
            }
        }
        return node;
    }

    // relocate stay_candidates and nn_candidates to target, then createSuper accordingly
    // level = 3 = super [default]
    // level = 2 = nn branch
    // level = 1 = stay branch
    inline size_t doTarget_StayNN(tt_data dt, size_t node, vector<size_t> stay_candidates, vector<size_t> nn_candidates, vector<size_t> &insertionList, size_t level = 3)
    {
        vector<size_t> target;
        switch (level)
        {
        case 1:
            target = dt.stay_branch;
            break;
        case 2:
            target = dt.nn_branch;
            break;
        case 3:
            target = dt.stay_super;
            break;
        default:
            fprintf(stderr, "ERROR - doTarget_StayNN\n");
            return node;
            break;
        }

        bool moved_all_stay = relocate(target, stay_candidates);
        bool moved_all_nn = relocate(target, nn_candidates);
        bool moved_all = moved_all_stay && moved_all_nn;

        if (checkdeleteUniSuper(node))
        {
            return node;
        }

        if (checkCreateSuper(node, moved_all, dt.mismatch, target.size()))
        { // there is some stay and/or nn candidates can't be moved to stay X, create new Super
            size_t t_parent = createSuper(node, insertionList, target);
            movedRemaining(node, t_parent, moved_all_stay, stay_candidates);
            movedRemaining(node, t_parent, moved_all_nn, nn_candidates);
            updateParentMean(t_parent);
            return t_parent;
        }

        return node;
    }

    // move low to high
    // then remaining low to target (stay_super)
    // then move high to target
    // step 2&3 reversible, same effect
    inline size_t doStaySuperLowHigh(tt_data dt, size_t node, vector<size_t> low_candidates, vector<size_t> high_candidates, vector<size_t> &insertionList)
    {
        vector<size_t> target = dt.stay_super;
        bool moved_all_low = relocate(high_candidates, low_candidates);
        moved_all_low = relocateRemaining(node, moved_all_low, low_candidates, target);
        bool moved_all_high = relocate(target, high_candidates);
        bool moved_all = moved_all_low && moved_all_high;

        if (!checkdeleteUniSuper(node))
        {
            if (checkCreateSuper(node, moved_all, dt.mismatch, target.size()))
            { // there is some nn highes can't be moved to stay super, create new super
                size_t t_parent = createSuper(node, insertionList, target);
                movedRemaining(node, t_parent, moved_all_low, low_candidates);
                movedRemaining(node, t_parent, moved_all_high, high_candidates);
                updateParentMean(t_parent);
                return t_parent;
            }
        }
        return node;
    }

    // relocate stay (level 1, default) or nn leaf (level 2) to stay branch, or else stay super, or else nn branch
    // then stay branch to super
    // then nn branch to super
    inline size_t doLeaf_others(tt_data dt, size_t node, vector<size_t> &insertionList, size_t level = 1)
    {
        vector<size_t> candidates;
        if (level == 1)
        {
            candidates = dt.stay_leaf;
        }
        else
        {
            candidates = dt.nn_leaf;
        }

        bool moved_all_leaf = relocate(dt.stay_branch, candidates);
        moved_all_leaf = relocateRemaining(node, moved_all_leaf, candidates, dt.stay_super);
        moved_all_leaf = relocateRemaining(node, moved_all_leaf, candidates, dt.nn_branch);
        bool moved_all_stay_branch = relocate(dt.stay_super, dt.stay_branch);
        bool moved_all_nn_branch = relocate(dt.stay_super, dt.nn_branch);
        bool moved_all = moved_all_leaf && moved_all_stay_branch && moved_all_nn_branch;

        if (!checkdeleteUniSuper(node))
        {
            if (checkCreateSuper(node, moved_all, dt.mismatch, candidates.size()))
            { // there is some nn highes can't be moved to stay super, create new super
                size_t t_parent = createSuper(node, insertionList, candidates);
                movedRemaining(node, t_parent, moved_all_leaf, candidates);
                movedRemaining(node, t_parent, moved_all_stay_branch, dt.stay_branch);
                movedRemaining(node, t_parent, moved_all_nn_branch, dt.nn_branch);

                updateParentMean(t_parent);
                return t_parent;
            }
        }
        return node;
    }

    // relocate candidates to high, then remaining to low
    inline size_t doCandidates_HighLow(tt_data dt, size_t node, vector<size_t> candidates, vector<size_t> high_parents, vector<size_t> low_parents, vector<size_t> &insertionList)
    {
        bool moved_all = relocate(high_parents, candidates);
        moved_all = relocateRemaining(node, moved_all, candidates, low_parents);

        if (!checkdeleteUniSuper(node))
        {
            if (isRootNode[node] || dt.mismatch != 0)
            {
                // create super for high and low
                size_t t_parent = createSuper(node, insertionList, high_parents);
                for (size_t child : low_parents)
                {
                    moveParent(child, t_parent);
                }

                // there is some candidates can't be moved to high or low, move to super
                movedRemaining(node, t_parent, moved_all, candidates);
                updateParentMean(t_parent);
                return t_parent;
            }
        }
        return node;
    }

    inline size_t doBit1(tt_data dt, const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t node)
    {
        size_t t_parent, t_branch;
        size_t pos;
        for (pos = 0; pos < dt.statuses.size(); ++pos)
        {
            if (dt.statuses.test(pos))
            {
                break;
            }
        }

        switch (pos)
        {
        case 1:
            printMsg(">>Stat>> Stay Leaf\n");
            if (dt.stay_leaf.size() > 1)
            {
                createBranch(node, insertionList, dt.stay_leaf);
            }
            return stayNode(signature, insertionList, idx, dt.dest);
            break;
        case 2:
            printMsg(">>Stat>> NN Leaf\n");
            t_branch = createBranch(node, insertionList, dt.nn_leaf);
            return insertBranch(signature, insertionList, idx, t_branch);
            // return createAmbiNode(signature, insertionList, t_branch, idx);
            break;
        case 3:
            printMsg(">>Stat>> Stay Branch\n");
            if (dt.stay_branch.size() > 1)
            {
                if (isRootNode[node] || dt.mismatch > 0)
                {
                    t_parent = createSuper(node, insertionList, dt.stay_branch);
                    recluster(t_parent);
                }
            }
            return insertBranch(signature, insertionList, idx, dt.dest_branch);
            break;
        case 4:
            printMsg(">>Stat>> NN Branch %zu\n", dt.mismatch);

            if (dt.nn_branch.size() > 1)
            {

                return doNNBranch(signature, node, dt, insertionList, idx);
            }
            else
            {
                return insertBranch(signature, insertionList, idx, dt.nn_branch[0]);
            }
            break;
        case 5:
            printMsg(">>Stat>> Stay Super\n");
            if (dt.stay_super.size() > 1)
            {
                if (isRootNode[node] || dt.mismatch != 0)
                {
                    t_parent = createSuper(node, insertionList, dt.stay_super);
                    recluster(t_parent);
                }
            }
            // return stayNode(signature, insertionList, idx, dt.dest_super);
            return tt(signature, insertionList, idx, dt.dest_super);
            break;
        default:
            fprintf(stderr, "ERROR 1 bit %zu!!", idx);
            return 0;
            break;
        }
    }

    inline size_t doBit2(tt_data dt, const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t node)
    {
        size_t t_parent, t_branch;
        size_t tempbits = (dt.statuses & bitset<8>("00101010")).count();
        // only stay
        if (tempbits == 2)
        {
            if (dt.statuses.test(5))
            {
                if (dt.statuses.test(1))
                {
                    printMsg(">>Stat>> Stay Leaf and Stay Super\n");
                    dt.dest = stayNode(signature, insertionList, idx, dt.dest);
                    doTarget_candidates(dt, node, dt.stay_leaf, insertionList);
                }
                else
                {
                    printMsg(">>Stat>> Stay Branch and Stay Super\n");
                    dt.dest = insertBranch(signature, insertionList, idx, dt.dest_branch);
                    doTarget_candidates(dt, node, dt.stay_branch, insertionList);
                }
            }
            else
            {
                printMsg(">>Stat>> Stay Leaf and Stay Branch\n");
                dt.dest = stayNode(signature, insertionList, idx, dt.dest);
                doTarget_candidates(dt, node, dt.stay_leaf, insertionList, 1);
            }
            return dt.dest;
        }
        else if (tempbits == 1)
        {
            if (dt.statuses.test(4))
            {
                if (dt.statuses.test(1))
                {
                    printMsg(">>Stat>> Stay Leaf and NN Branch\n");
                    dt.dest = stayNode(signature, insertionList, idx, dt.dest);
                    doTarget_candidates(dt, node, dt.stay_leaf, insertionList, 2);
                    return dt.dest;
                }
                else if (dt.statuses.test(3))
                {
                    printMsg(">>Stat>> Stay Branch and NN Branch\n");
                    dt.dest = insertBranch(signature, insertionList, idx, dt.dest_branch);
                    if (isRootNode[node] || dt.mismatch != 0)
                    {
                        t_parent = createSuper(node, insertionList, dt.stay_branch);
                        for (size_t b : dt.nn_branch)
                        {
                            moveParent(b, t_parent);
                        }
                        // recluster(t_parent);
                    }
                    return dt.dest;
                }
                else if (dt.statuses.test(5))
                { //?
                    printMsg(">>Stat>> NN Branch and Stay Super\n");
                    doTarget_candidates(dt, node, dt.nn_branch, insertionList);
                    return tt(signature, insertionList, idx, dt.dest_super);
                }
            }
            else
            {
                // printMsg("NN Leaf and ");
                if (dt.statuses.test(1))
                {
                    printMsg(">>Stat>> Stay Leaf and NN Leaf\n");
                    dt.dest = stayNode(signature, insertionList, idx, dt.dest);
                    t_branch = createBranch(node, insertionList, dt.stay_leaf);
                    for (size_t l : dt.nn_leaf)
                    {
                        moveParent(l, t_branch);
                    }
                    return dt.dest;
                }
                else if (dt.statuses.test(3))
                {
                    printMsg(">>Stat>> NN Leaf and Stay Branch\n");
                    dt.dest = insertBranch(signature, insertionList, idx, dt.dest_branch);
                    doTarget_candidates(dt, node, dt.nn_leaf, insertionList, 1);
                    return dt.dest;
                }
                else if (dt.statuses.test(5))
                {
                    printMsg(">>Stat>> NN Leaf and Stay Super\n");
                    doTarget_candidates(dt, node, dt.nn_leaf, insertionList);
                    return tt(signature, insertionList, idx, dt.dest_super);
                }
            }
        }
        else
        {
            printMsg(">>Stat>> NN Leaf and NN Branch\n");

            t_parent = doTarget_candidates(dt, node, dt.nn_leaf, insertionList, 2);

            for (size_t child : dt.nn_leaf)
            {
                size_t parent = parentLinks[child];
                if (parent != node)
                {
                    isAmbiNode[child] = 1;
                    ambiLinks[parent].push_back(child);
                }
            }

            for (size_t branch : dt.nn_branch)
            {
                mergeAmbi(branch);
            }
            mergeAmbi(node);
            updateParentMean(node);
            // this works better for non force split
            return createNode(signature, insertionList, t_parent, idx);

            // // this works better for force split
            // return tt(signature, insertionList, idx, node);
        }

        fprintf(stderr, "ERROR 2 bit %zu!!", idx);
        return 0;
    }

    inline size_t doBit3(tt_data dt, const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t node)
    {
        size_t t_parent, t_branch;
        if (dt.statuses.test(5))
        {
            if (dt.statuses.test(1))
            {
                dt.dest = stayNode(signature, insertionList, idx, dt.dest);
                if (dt.statuses.test(2))
                {
                    printMsg(">>Stat>> Stay Leaf and NN Leaf and Stay Super\n");
                    doTarget_StayNN(dt, node, dt.stay_leaf, dt.nn_leaf, insertionList);
                }
                else if (dt.statuses.test(3))
                {
                    printMsg(">>Stat>> Stay Leaf and Stay Branch and Stay Super\n");
                    doStaySuperLowHigh(dt, node, dt.stay_leaf, dt.stay_branch, insertionList);
                }
                else
                {
                    printMsg(">>Stat>> Stay Leaf and NN branch and Stay Super\n");
                    doTarget_StayNN(dt, node, dt.stay_leaf, dt.nn_branch, insertionList);
                }

                return dt.dest;
            }
            else if (dt.statuses.test(2))
            {
                // printMsg("Stay Super and NN Leaf and ");
                if (dt.statuses.test(3))
                {
                    printMsg(">>Stat>> NN Leaf and Stay Branch and Stay Super\n");
                    dt.dest = insertBranch(signature, insertionList, idx, dt.dest_branch);
                    doStaySuperLowHigh(dt, node, dt.nn_leaf, dt.stay_branch, insertionList);
                    return dt.dest;
                }
                else
                {
                    printMsg(">>Stat>> NN Leaf and NN branch and Stay Super\n");
                    doStaySuperLowHigh(dt, node, dt.nn_leaf, dt.nn_branch, insertionList);
                    return tt(signature, insertionList, idx, dt.dest_super);
                }
            }
            else
            {
                printMsg(">>Stat>> Stay Branch and NN branch and Stay Super\n");
                dt.dest = insertBranch(signature, insertionList, idx, dt.dest_branch);
                doTarget_StayNN(dt, node, dt.stay_branch, dt.nn_branch, insertionList);
                return dt.dest;
            }
        }
        else
        {
            if (dt.statuses.test(1))
            {
                if (dt.statuses.test(2))
                {
                    if (dt.statuses.test(3))
                    {
                        printMsg(">>Stat>> Stay Leaf and NN Leaf and Stay Branch\n");
                        dt.dest = stayNode(signature, insertionList, idx, dt.dest);
                        doTarget_StayNN(dt, node, dt.stay_leaf, dt.nn_leaf, insertionList, 1);
                    }
                    else
                    {
                        printMsg(">>Stat>> Stay Leaf and NN Leaf and NN Branch\n");
                        dt.dest = stayNode(signature, insertionList, idx, dt.dest);
                        doTarget_StayNN(dt, node, dt.stay_leaf, dt.nn_leaf, insertionList, 2);
                    }
                }
                else
                {
                    printMsg(">>Stat>> Stay Leaf and Stay Branch and NN branch\n");
                    dt.dest = stayNode(signature, insertionList, idx, dt.dest);
                    doCandidates_HighLow(dt, node, dt.stay_leaf, dt.stay_branch, dt.nn_branch, insertionList);
                }
            }
            else
            {
                printMsg(">>Stat>> NN Leaf and Stay Branch and NN branch\n");
                dt.dest = insertBranch(signature, insertionList, idx, dt.dest_branch);
                doCandidates_HighLow(dt, node, dt.nn_leaf, dt.stay_branch, dt.nn_branch, insertionList);
            }

            return dt.dest;
        }

        fprintf(stderr, "ERROR 3 bit %zu!!", idx);
        return 0;
    }

    inline size_t doBit4(tt_data dt, const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t node)
    {
        size_t t_parent, t_branch;

        if (dt.statuses.test(1))
        {
            dt.dest = stayNode(signature, insertionList, idx, dt.dest);
            if (dt.statuses.test(2))
            {
                insertVecRange(dt.stay_leaf, dt.nn_leaf);
                if (dt.statuses.test(3))
                {
                    if (dt.statuses.test(4))
                    {
                        printMsg(">>Stat>> Stay Leaf and NN Leaf and Stay Branch and NN Branch\n");
                        doCandidates_HighLow(dt, node, dt.stay_leaf, dt.stay_branch, dt.nn_branch, insertionList);
                    }
                    else
                    {
                        printMsg(">>Stat>> Stay Leaf and NN Leaf and Stay Branch and Stay Super\n");
                        doStaySuperLowHigh(dt, node, dt.stay_leaf, dt.stay_branch, insertionList);
                    }
                }
                else
                {
                    printMsg(">>Stat>> Stay Leaf and NN Leaf and NN Branch and Stay Super\n");
                    doTarget_StayNN(dt, node, dt.stay_leaf, dt.nn_branch, insertionList);
                }
            }
            else
            {
                printMsg(">>Stat>> Stay Leaf and Stay Branch and NN Branch and Stay Super\n");
                dt.dest = stayNode(signature, insertionList, idx, dt.dest);
                doLeaf_others(dt, node, insertionList);
            }
            return dt.dest;
        }
        else
        {
            printMsg(">>Stat>> NN Leaf and Stay Branch and NN Branch and Stay Super\n");
            dt.dest = insertBranch(signature, insertionList, idx, dt.dest_branch);
            doLeaf_others(dt, node, insertionList, 2);
            return dt.dest;
        }

        fprintf(stderr, "ERROR 4 bit %zu!!", idx);
        return 0;
    }

    // dest cannot be super or root
    // in the case of stay size = 1
    // dest = stay[0]
    inline size_t growtree_without_root(tt_data dt, const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        if (dt.statuses.none())
        {
            return createNode(signature, insertionList, node, idx);
        }

        size_t t_parent, t_branch;
        size_t setbits = dt.statuses.count();
        if (setbits == 1)
        {
            return doBit1(dt, signature, insertionList, idx, node);
        }
        else if (setbits == 2)
        {
            return doBit2(dt, signature, insertionList, idx, node);
        }
        else if (setbits == 3)
        {
            return doBit3(dt, signature, insertionList, idx, node);
        }
        else if (setbits == 4)
        {
            return doBit4(dt, signature, insertionList, idx, node);
        }
        else if (setbits == 5)
        {
            printMsg(">>Stat>> Stay Leaf and NN Leaf and Stay Branch and NN Branch and Stay Super\n");
            t_branch = createBranch(node, insertionList, dt.stay_leaf);
            for (size_t l : dt.nn_leaf)
            {
                moveParent(l, t_branch);
            }
            printMsg("combined leaves to branch %zu\n", t_branch);
            dt.stay_branch.push_back(t_branch);
            dt.dest_branch = t_branch;
            dt.statuses.reset(1);
            dt.statuses.reset(2);
            return doBit3(dt, signature, insertionList, idx, node);
        }

        printMsg(">>Stat>> ERROR!!\n");
        return 0;
    }

    inline size_t tt(const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        tt_data dt = getStatus(signature, node);
        return growtree_without_root(dt, signature, insertionList, idx, node);
    }

    inline void dissolveSingleton(size_t singleton, size_t dest)
    {
        printMsg(">>> Dissolve singleton %zu to %zu\n", singleton, dest);
        seqIDs[dest].push_back(seqIDs[singleton][0]);
        addSigToMatrix(dest, getMeanSig(singleton));
        deleteNode(singleton);
        clearNode(singleton);
    }

    inline void dissolveSuper(size_t parent)
    {
        printMsg(">>> Dissolving super %zu\n", parent);
        size_t grandparent = parentLinks[parent];
        vector<size_t> candidates = childLinks[parent];
        for (size_t child : candidates)
        {
            moveParent(child, grandparent);
        }
        deleteNode(parent);
        clearNode(parent);
        updateParentMean(grandparent);
    }

    inline size_t tt_root(const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        if (node != root && priority[node] < split_threshold)
        {
            printMsg("add subtree %zu, %.2f\n", node, priority[node]);
            size_t parent = parentLinks[node];
            if (addSubtree(node, insertionList) == 0)
            {
                printMsg("parent %zu\n", parent);
                size_t last_child = childLinks[node][childCounts[node] - 1];
                printMsg("last_child %zu\n", last_child);
                if (relocate(childLinks[parent], vector<size_t>{last_child}))
                {
                    printMsg(">> Relocated %zu\n", last_child);
                    return tt_root(signature, insertionList, idx, parent);
                }
            }
            else
            {
                return tt_root(signature, insertionList, idx, parent);
            }
        }

        if (!isRootNode[childLinks[node][0]])
        {
            return tt(signature, insertionList, idx, node);
        }

        double max_similarity = 0;
        size_t best_root = 0;
        for (size_t child : childLinks[node])
        {
            double similarity = calcSimilarityWrap(getMeanSig(child), signature);

            if (similarity + split_node_threshold >= priority[child])
            {
                printMsg("***(root %zu, %.2f, %.2f)\n", child, priority[child], similarity);
            }
            else
            {
                printMsg("(root %zu, %.2f, %.2f)\n", child, priority[child], similarity);
            }
            if (similarity > max_similarity)
            {
                max_similarity = similarity;
                best_root = child;
            }
        }
        size_t dest = 0;
        printMsg("traverse root %zu\n", best_root);
        dest = tt_root(signature, insertionList, idx, best_root);

        if (childCounts[best_root] > tree_order)
        {

            if (forceSplitRoot(insertionList, best_root) == 1)
            {
                dissolveSuper(best_root);
            }
        }
        return dest;
    }

    inline size_t first_insert(const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t parent = 0)
    {
        return createNode(signature, insertionList, root, idx);
    }

    inline size_t insert(const_s_type signature, vector<size_t> &insertionList, size_t idx)
    {
        size_t node = tt(signature, insertionList, idx);
        printMsg("inserted %zu at %zu\n\n", idx, node);
        return node;
    }

    size_t addSubtree(size_t node, vector<size_t> &insertionList)
    {
        size_t clusterCount = 2;
        vector<size_t> children = childLinks[node];
        vector<vector<size_t>> clusters(clusterCount);

        sVec_type temp_centroids = createRandomSigs(node, clusterCount);

        for (size_t child : children)
        {
            size_t dest = 0;
            double max_similarity = 0;
            for (size_t i = 0; i < clusterCount; i++)
            {
                // double similarity = calcSimilarityWrap(temp_centroids[i], getMeanSig(child));
                double similarity = calcSimilaritySigToNode(child, temp_centroids, i);
                printMsg("%zu, %f\n", i, similarity);
                if (similarity > max_similarity)
                {
                    max_similarity = similarity;
                    dest = i;
                }
            }
            printMsg("--- %zu goes to %zu\n", child, dest);
            clusters[dest].push_back(child);
        }

        for (vector<size_t> cluster : clusters)
        {
            if (cluster.size() <= 1)
            {
                // printTreeJson(stderr);
                printMsg("??something is wrong cannot add subtree evenly\n");
                return 0;
            }
        }

        if (print_)
        {
            printTreeJson(stderr);
        }
        printMsg(">> addSubtree %zu\n", node);

        size_t t_parent = createSuper(parentLinks[node], insertionList, clusters[1]);
        isSuperNode[t_parent] = 0;
        isRootNode[t_parent] = 1;

        updateParentMean(node);
        if (print_)
        {
            printTreeJson(stderr);
        }

        return 1;
    }

    size_t forceSplitRoot(vector<size_t> &insertionList, size_t node = 0, size_t clusterCount = 2)
    {
        // size_t clusterCount = 2;
        vector<size_t> children = childLinks[node];
        vector<vector<size_t>> clusters(clusterCount);

        size_t child = children[0];
        children = childLinks[node];
        sVec_type temp_centroids = createRandomSigs(node, clusterCount);

        for (size_t child : children)
        {
            size_t dest = 0;
            double max_similarity = 0;
            for (size_t i = 0; i < clusterCount; i++)
            {
                // double similarity = calcSimilarityWrap(temp_centroids[i], getMeanSig(child));
                double similarity = calcSimilaritySigToNode(child, temp_centroids, i);
                printMsg("%zu, %f\n", i, similarity);
                if (similarity > max_similarity)
                {
                    max_similarity = similarity;
                    dest = i;
                }
            }
            printMsg("--- %zu goes to %zu\n", child, dest);
            clusters[dest].push_back(child);
        }

        for (vector<size_t> cluster : clusters)
        {
            if (cluster.size() <= 1)
            {
                printMsg("??something is wrong cannot split evenly\n");
                return 0;
            }
        }

        if (print_)
        {
            printTreeJson(stderr);
        }
        printMsg(">> forceSplitRoot %zu\n", node);

        // reuse clusterSize to store the new t_parents
        for (vector<size_t> cluster : clusters)
        {
            size_t t_parent = createParent(node, insertionList);
            isRootNode[t_parent] = 1;
            for (size_t child : cluster)
            {
                moveParent(child, t_parent);
            }
            updateNodeMean(t_parent);
            addSigToMatrix(node, getMeanSig(t_parent));
        }
        updateParentMean(node);
        if (print_)
        {
            printTreeJson(stderr);
        }

        if (node != 0)
        {
            recluster(node);
            if (print_)
            {
                printTreeJson(stderr);
            }
        }

        return 1;
    }

    inline size_t insertSplitRoot(const_s_type signature, vector<size_t> &insertionList, size_t idx)
    {
        size_t node = tt_root(signature, insertionList, idx);
        printMsg("inserted %zu at %zu,%zu,%zu\n\n", idx, node, childCounts[root], tree_order);
        if (childCounts[root] > tree_order)
        {
            forceSplitRoot(insertionList);
        }
        return node;
    }

    inline size_t searchBest(const_s_type signature, size_t node = 0)
    {
        size_t dest = 0;
        double max_similarity = 0;

        for (size_t child : childLinks[node])
        {
            double similarity = calcSimilarityWrap(getMeanSig(child), signature);
            printMsg(" (%zu, %.2f) ", child, similarity);

            if (similarity > max_similarity)
            {
                max_similarity = similarity;
                dest = child;
            }
        }
        printMsg("\n>>>%zu\n", dest);

        if (isBranchNode[dest])
        {
            return searchBest(signature, dest);
        }
        else
        {
            return dest;
        }
    }

    inline size_t deepFirstSearch(const_s_type signature, size_t node = 0)
    {
        if (isBranchNode[node])
        {
            return searchBest(signature, node);
        }
        size_t dest = 0;
        double max_similarity = 0;
        printMsg("\n------------- %zu\n", node);

        for (size_t subtree : childLinks[node])
        {
            size_t candidate = searchBestSubtree(signature, subtree);
            double similarity = calcSimilarityWrap(getMeanSig(candidate), signature);
            printMsg(" (%zu, %.2f) ", candidate, similarity);

            if (similarity > max_similarity)
            {
                max_similarity = similarity;
                dest = candidate;
            }
        }

        return dest;
    }

    inline size_t searchBestSubtree(const_s_type signature, size_t node = 0)
    {
        vector<vector<size_t>> children = separateRootChildren(node);
        if (children[0].size() == 0)
        {
            return search(signature, node);
        }

        size_t dest = 0;
        double max_similarity = 0;

        // do roots
        for (size_t subtree : children[0])
        {
            size_t candidate = searchBestSubtree(signature, subtree);
            double similarity = calcSimilarityWrap(getMeanSig(candidate), signature);

            if (similarity > max_similarity)
            {
                max_similarity = similarity;
                dest = candidate;
            }
        }
        // do non-roots
        if (children[1].size() > 0)
        {
            size_t best_non_root = selectiveSearch(signature, children[1]);
            double similarity = calcSimilarityWrap(getMeanSig(best_non_root), signature);

            if (similarity > max_similarity)
            {
                max_similarity = similarity;
                dest = best_non_root;
            }
        }
        return dest;
    }

    // find highest overlap among children
    // break tie by checking similarity of the finals
    inline size_t search(const_s_type signature, size_t node = 0)
    {
        double best_similarity = 0;
        vector<size_t> candidates;

        for (size_t child : childLinks[node])
        {
            double similarity = calcSimilarityWrap(getMeanSig(child), signature);
            printMsg(" <%zu,%.2f> ", child, similarity);

            if (similarity > best_similarity)
            {
                best_similarity = similarity;
                candidates.clear();
            }

            if (similarity == best_similarity)
            {
                candidates.push_back(child);
            }
        }

        size_t best_child = candidates[0];

        if (candidates.size() == 1)
        {
            if (isBranchNode[best_child])
            {
                printMsg("\n>%zu\n", best_child);
                return search(signature, best_child);
            }
            else
            {
                return best_child;
            }
        }
        else
        {

            printMsg("\n*** here\n");
            double best_dest_similarity = 0;
            for (size_t child : candidates)
            {
                printMsg("child %zu\n", child);
                size_t leaf = child;
                double similarity = best_similarity;
                if (isBranchNode[child])
                {
                    leaf = search(signature, child);
                    similarity = calcSimilarityWrap(getMeanSig(child), signature);
                    // similarity += calcOverlap(signature, getMeanSig(leaf));
                    printMsg("> leaf %zu\n", leaf);
                }

                if (similarity > best_dest_similarity)
                {
                    best_dest_similarity = similarity;
                    best_child = leaf;
                }
            }

            return best_child;
        }
    }

    // find highest overlap among children
    // break tie by checking similarity of the finals
    inline size_t selectiveSearch(const_s_type signature, vector<size_t> children)
    {
        double best_similarity = 0;
        vector<size_t> candidates;

        for (size_t child : children)
        {
            double similarity = calcSimilarityWrap(getMeanSig(child), signature);
            if (similarity > best_similarity)
            {
                best_similarity = similarity;
                candidates.clear();
            }

            if (similarity == best_similarity)
            {
                candidates.push_back(child);
            }
        }

        size_t best_child = candidates[0];

        if (candidates.size() == 1)
        {
            if (isBranchNode[best_child])
            {
                printMsg("\n>%zu\n", best_child);
                return search(signature, best_child);
            }
            else
            {
                return best_child;
            }
        }
        else
        {

            printMsg("\n*** here\n");
            double best_dest_similarity = 0;
            for (size_t child : candidates)
            {
                printMsg("child %zu\n", child);
                size_t leaf = child;
                double similarity = best_similarity;
                if (isBranchNode[child])
                {
                    leaf = search(signature, child);
                    similarity = calcSimilarityWrap(getMeanSig(child), signature);
                    // similarity += calcOverlap(signature, getMeanSig(leaf));
                    printMsg("> leaf %zu\n", leaf);
                }

                if (similarity > best_dest_similarity)
                {
                    best_dest_similarity = similarity;
                    best_child = leaf;
                }
            }

            return best_child;
        }
    }

    // cluster means still remain
    void prepReinsert(size_t node = 0)
    {
        // if (!isAmbiNode[node])
        {
            updatePriority(node);
            if (!isBranchNode[node])
            {
                updateNodeMean(node);
            }
        }
        seqIDs[node].clear();
        for (size_t child : childLinks[node])
        {
            prepReinsert(child);
        }
    }

    void updateTree(size_t node = 0)
    {
        // if (!isAmbiNode[node])
        {
            updatePriority(node);
            if (!isBranchNode[node])
            {
                updateNodeMean(node);
            }
        }
        for (size_t child : childLinks[node])
        {
            updateTree(child);
        }
    }

    size_t reinsert(const_s_type signature, size_t idx)
    {
        // size_t node = search(signature);
        size_t node = searchBestSubtree(signature);
        // size_t node = search2(signature);
        seqIDs[node].push_back(idx);
        addSigToMatrix(node, signature);
        return node;
    }

    void removeAmbi(size_t node = 0)
    {
        size_t cleared = 0;

        vector<size_t> children = childLinks[node];
        if (isRootNode[node] || isSuperNode[node])
        {
            for (size_t child : children)
            {
                removeAmbi(child);
            }
        }
        else if (isBranchNode[node])
        {
            if (priority[node] >= stay_threshold)
            {
                mergeBranch(node);
            }
            else
            {
                for (size_t child : children)
                {
                    removeAmbi(child);
                }
            }
        }
        else if (isAmbiNode[node] || seqIDs[node].size() <= singleton)
        {
            deleteNode(node);
            cleared = node;
        }
        else
        {
            updateParentMean(node);
        }
        size_t parent = parentLinks[node];
        while (childCounts[parent] <= 1 && parent != root)
        {
            deleteUnitig(parent);
            parent = parentLinks[parent];
        }

        if (cleared != 0)
        {
            clearNode(cleared);
        }
    }

    void outputHierarchy(FILE *pFile, size_t node = 0, size_t rank = 0)
    {
        for (size_t child : childLinks[node])
        {
            fprintf(pFile, "%zu,%zu,%zu\n", node, child, rank);
            if (isBranchNode[child])
            {
                outputHierarchy(pFile, child, rank + 1);
            }
        }
    }

    void destroyLocks(size_t node)
    {
        omp_destroy_lock(&locks[node]);
        if (isBranchNode[node])
        {
            for (size_t i = 0; i < childCounts[node]; i++)
            {
                destroyLocks(childLinks[node][i]);
            }
        }
    }

    inline void destroyLocks()
    {
        destroyLocks(root);
    }
};

#endif