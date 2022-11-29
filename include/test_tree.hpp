// There are two kinds of ktree nodes- branch nodes and leaf nodes
// Both contain a signature matrix, plus their own signature
// (the root node signature does not matter and can be blank)
// Branch nodes then contain 'order' links to other nodes
// Leaf nodes do not.
// However, as leaf nodes may become branch nodes, we allocate
// the space anyway.
// As the space to be used is determined at runtime, we use
// parallel arrays, not structs

#ifndef INCLUDE_test_tree_HPP
#define INCLUDE_test_tree_HPP

#include <omp.h>
#include <unordered_map>
#include <unordered_set>
#include <experimental/filesystem>
#include "bloom_filter.hpp"

// #include "read.hpp"

using namespace std;
typedef pair<float, float> data_type;

double split_threshold = 5;
double stay_threshold = 1;
size_t minimiser_match_threshold = 4;
size_t partree_capacity = 10000;

data_type createMeanSig(const vector<data_type> clusterSigs)
{
    data_type meanSig = make_pair(0, 0);

    for (data_type signature : clusterSigs)
    {
        meanSig.first += signature.first;
        meanSig.second += signature.second;
    }

    meanSig.first = meanSig.first / clusterSigs.size();
    meanSig.second = meanSig.second / clusterSigs.size();

    return meanSig;
}

double calcDistance(data_type a, data_type b)
{
    float x = a.first - b.first;
    float y = a.second - b.second;
    return sqrt(x * x + y * y);
}

// RMSD
double calcDistortion(const vector<data_type> clusterSigs)
{
    data_type meanSig = createMeanSig(clusterSigs);
    double sumSquareDistance = 0;

    for (data_type signature : clusterSigs)
    {
        double distance = calcDistance(meanSig, signature);
        sumSquareDistance += distance * distance;
    }
    return sqrt(sumSquareDistance / clusterSigs.size());
}

class test_tree
{
public:
    size_t root = 0;                    // # of root node
    vector<size_t> childCounts;         // n entries, number of children
    vector<size_t> branchCounts;        // n entries, number of children which are a branch
    vector<int> isBranchNode;           // n entries, is this a branch node
    vector<vector<size_t>> childLinks;  // n * o entries, links to children
    vector<vector<size_t>> seqIDs;      // n * o entries, links to children
    vector<size_t> parentLinks;         // n entries, links to parents
    vector<float> priority;             // n entries, links to parents
    vector<data_type> means;            // n * signatureSize entries, node signatures
    vector<vector<data_type>> matrices; // capacity * signatureSize * n
    vector<omp_lock_t> locks;           // n locks
    size_t capacity = 0;                // Set during construction, currently can't change

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
                branchCounts.resize(capacity);
            }
            //#pragma omp single
            {
                isBranchNode.resize(capacity);
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
                means.resize(capacity);
            }
        }
    }

    test_tree(size_t capacity)
    {
        reserve(capacity);
        childCounts[root] = 0;
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
        fprintf(stream, "{\"node\":\"%zu\",\"branch\":\"%zu\",\"priority\":\"%.2f\",\"childCount\":\"%zu\",\"content\":\"*", tnode, isBranchNode[tnode], priority[tnode], seqIDs[tnode].size());
        for (size_t seq : seqIDs[tnode])
        {
            fprintf(stream, "%zu,", seq);
        }
        fprintf(stream, "\",\"children\":[");
    }

    void printSubTreeJson(FILE *stream, size_t tnode)
    {
        if (childCounts[tnode] > 0)
        {
            printNodeJson(stream, tnode);
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
            fprintf(stream, "{\"node\":\"%zu\",\"priority\":\"%.2f\",\"childCount\":\"%zu\",\"content\":\"*", tnode, priority[tnode], seqIDs[tnode].size());
            for (size_t seq : seqIDs[tnode])
            {
                fprintf(stream, "%zu,", seq);
            }
            fprintf(stream, "\"}");
        }
    }

    void printTreeJson(FILE *stream)
    {
        fprintf(stream, "var treeData = ");
        printSubTreeJson(stream, root);
        fprintf(stream, ";\n");
    }

    void addSigToMatrix(size_t node, data_type signature)
    {
        matrices[node].push_back(signature);
    }

    void clearNode(size_t node)
    {
        isBranchNode[node] = 0;
        childCounts[node] = 0;
        childLinks[node].clear();
        seqIDs[node].clear();
        matrices[node].clear();
        means[node] = createMeanSig(matrices[node]);
        priority[node] = 0;
    }

    // merge if number of children is greater than tree order
    // decrease split threshold
    // to allow node with larger distortion
    // the local split threshold should be gradual, but we just make it 110% of the og threshold
    // stay threshold should have nothing to do here, if not split just stay
    // at this stage, we don't do k-means
    // implement as if this is an independant tree
    void mergeChildren(size_t node, vector<size_t> &insertionList)
    {
        // double local_split = split_threshold * 1.1;

        double local_split = minimiser_match_threshold;

        vector<data_type> temp_means;
        vector<size_t> temp_childLink;

        // insert first child
        size_t child = childLinks[node][0];
        temp_childLink.push_back(child);
        temp_means.push_back(means[child]);

        for (size_t i = 1; i < childCounts[node]; i++)
        {
            child = childLinks[node][i];
            size_t n = 0;
            for (data_type mean : temp_means)
            {
                float distance = calcDistance(means[child], mean);
                // fprintf(stderr, "Distance:%.2f, %.2f, %.2f\n", distance,split_threshold, local_split);
                if (distance < local_split)
                {
                    fprintf(stderr, "Merge %zu into %zu\n", child, childLinks[node][n]);
                    break;
                }
                n++;
            }
            fprintf(stderr, "n:%zu\n", n);
            if (n == temp_means.size())
            {
                temp_childLink.push_back(child);
                temp_means.push_back(means[child]);
            }
        }

        fprintf(stderr, "Before: %zu, After:%zu\n", childCounts[node], temp_childLink.size());
    }

    // union mean of children
    inline void updateNodeMean(size_t node)
    {
        means[node] = createMeanSig(matrices[node]);
        //? p
        priority[node] = calcDistortion(matrices[node]);

        size_t parent = parentLinks[node];

        if (isBranchNode[parent])
        {
            int idx = -1;
            for (size_t i = 0; i < childCounts[parent]; i++)
            {
                if (childLinks[parent][i] == node)
                {
                    idx = i;
                    break;
                }
            }

            if (idx == -1)
            {
                fprintf(stderr, "ERROR!!\n");
            }
            else
            {
                matrices[parent][idx] == means[node];
                updateNodeMean(parent);
            }
        }
    }

    // return if found same, and the destination node
    inline tuple<bool, size_t> traverse(data_type signature, vector<size_t> &insertionList)
    {
        size_t node = root;
        double local_stay = stay_threshold;
        double offset = 0.1;

        while (childCounts[node] > 0)
        {
            vector<size_t> mismatch;
            vector<size_t> matching;

            for (size_t i = 0; i < childCounts[node]; i++)
            {
                size_t child = childLinks[node][i];
                double distance = calcDistance(means[child], signature);

                fprintf(stderr, " <%zu,%.2f> ", child, distance);

                // found same, move on with the next seq
                if (distance <= (local_stay))
                {
                    // priority[child]++;
                    return make_tuple(true, child);
                }

                // count how many nodes mismatch
                if (distance >= (split_threshold))
                {
                    // mismatch++;
                    mismatch.push_back(child);
                }
                else
                {
                    matching.push_back(child);
                }
            }

            // nothing is close enough, spawn new child under parent
            if (mismatch.size() == childCounts[node])
            {
                fprintf(stderr, " -no match:%zu, ", node);
                return make_tuple(false, node);
            }
            else if (matching.size() > 1) // matching with multiple
            {
                size_t temp = createNode(signature, insertionList, node);

                parentLinks[temp] = node;

                childCounts[node] = mismatch.size() + 1;
                childLinks[node].clear();
                childLinks[node] = mismatch;
                childLinks[node].push_back(temp);

                childCounts[temp] = matching.size();

                fprintf(stderr, "\nxxx multiple: %zu,%zu\n", temp, node);

                for (size_t n : matching)
                {
                    for (data_type sig : matrices[n])
                    {
                        addSigToMatrix(temp, means[n]);
                    }

                    parentLinks[n] = temp;
                    childLinks[temp].push_back(n);
                    fprintf(stderr, "leaves: %zu\n", n);
                }
                // fprintf(stderr, "leaves: %zu\n",childCounts[temp]);
                // printSubTreeJson(stderr,root);

                //?
                priority[temp] = childCounts[temp] - 1;

                return make_tuple(true, temp);
            }

            else if (matching.size() == 1 && childCounts[node] == 1) // matching with multiple
            {
                size_t match = matching[0];

                if (isBranchNode[match])
                {
                    node = match;
                    fprintf(stderr, "matching with branch %zu\n", node);
                }
                else
                {
                    childCounts[node] = 0;
                    childLinks[node].clear();
                    size_t tempUnion = createNode(signature, insertionList, node);
                    parentLinks[tempUnion] = node;
                    isBranchNode[tempUnion] = 1;
                    branchCounts[node]++;
                    childLinks[node][0] = tempUnion;

                    size_t newNode = createNode(signature, insertionList, tempUnion);
                    childCounts[tempUnion]++;
                    childLinks[tempUnion].push_back(match);

                    fprintf(stderr, "\nxxx here: %zu>%zu,%zu\n", tempUnion, newNode, match);

                    addSigToMatrix(tempUnion, matrices[match][0]);
                    updateNodeMean(tempUnion);
                    return make_tuple(true, newNode);
                }
            }
            // updateNodeMean(node);
            node = matching[0]; // choose the first matching for now, may change to best
        }

        return std::make_tuple(false, node);
    }

    inline bool hasBranch(size_t node)
    {
        for (size_t child : childLinks[node])
        {
            if (isBranchNode[child])
            {
                return true;
            }
        }
        return false;
    }

    inline size_t traverse(data_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        size_t current_childCount = childCounts[node];

        vector<size_t> mismatch;
        // near neighbour
        vector<size_t> NN_branches;
        vector<size_t> NN_leaves;

        for (size_t i = 0; i < childCounts[node]; i++)
        {
            size_t child = childLinks[node][i];
            double distance = calcDistance(means[child], signature);

            fprintf(stderr, " <%zu,%.2f> ", child, distance);

            if (distance <= stay_threshold)
            {
                // stay in branch, check its children
                if (isBranchNode[child])
                {
                    fprintf(stderr, "stay in branch %zu\n", child);
                    return traverse(signature, insertionList, idx, child);
                }
                else
                {
                    seqIDs[child].push_back(idx);
                    addSigToMatrix(child, signature);
                    updateNodeMean(child);
                    return child;
                }
            }

            if (distance > split_threshold)
            {
                mismatch.push_back(child);
            }
            else if (isBranchNode[child])
            {
                NN_branches.push_back(child);
            }
            else
            {
                NN_leaves.push_back(child);
            }
        }

        if (mismatch.size() == current_childCount) // All diff => split
        {
            // // maybe proceed with children if there is only 1 branch in this subtree and its a NN
            // if (NN_branches.size() == 1 && current_childCount == 1)
            // {
            //     node = NN_branches[0];
            // }

            size_t newNode = createNode(signature, insertionList, node);
            seqIDs[newNode].push_back(idx);
            return newNode;
        }
        // all children of this node are leaves, no branch at all
        else if (!hasBranch(node))
        {
            // match with all nodes (leaves)
            if (mismatch.size() == 0)
            {
                if (isBranchNode[node])
                {
                    size_t newChild = createNode(signature, insertionList, node);
                    seqIDs[newChild].push_back(idx);
                    return newChild;
                }
                else
                {
                    size_t tempParent = getNewNodeIdx(insertionList);
                    fprintf(stderr, "1 created tempParent %zu\n", tempParent);
                    isBranchNode[tempParent] = 1;
                    branchCounts[node]++;
                    parentLinks[tempParent] = node;
                    childCounts[tempParent] = childCounts[node];
                    childLinks[tempParent] = childLinks[node];

                    childCounts[node] = 1;
                    childLinks[node].clear();
                    childLinks[node].push_back(tempParent);

                    for (size_t n : childLinks[tempParent])
                    {
                        addSigToMatrix(tempParent, means[n]);
                        parentLinks[n] = tempParent;
                    }

                    size_t newChild = createNode(signature, insertionList, tempParent);
                    seqIDs[newChild].push_back(idx);
                    return newChild;
                }
            }
            // match with some leaves only
            else
            {

                size_t tempParent = getNewNodeIdx(insertionList);
                fprintf(stderr, "2 created tempParent %zu\n", tempParent);
                isBranchNode[tempParent] = 1;
                branchCounts[node]++;
                parentLinks[tempParent] = node;
                childCounts[tempParent] = NN_leaves.size();
                childLinks[tempParent] = NN_leaves;
                for (size_t n : NN_leaves)
                {
                    addSigToMatrix(tempParent, means[n]);
                    parentLinks[n] = tempParent;
                }

                childLinks[node].clear();
                childLinks[node].push_back(tempParent);
                matrices[node].clear();
                addSigToMatrix(node, means[tempParent]);

                size_t newChild = createNode(signature, insertionList, tempParent);
                seqIDs[newChild].push_back(idx);

                fprintf(stderr, "here %zu, %zu, %zu\n", node, isBranchNode[node], NN_leaves.size());

                if (isBranchNode[node])
                {
                    fprintf(stderr, "a\n");
                    childCounts[node] = 2;

                    // merge mismatch leaves
                    // put NN in one branch, mismatch in another branch
                    if (mismatch.size() > 1)
                    {
                        fprintf(stderr, "a\n");
                        size_t mismatchBranch = getNewNodeIdx(insertionList);
                        fprintf(stderr, "created mismatchBranch %zu\n", mismatchBranch);
                        isBranchNode[mismatchBranch] = 1;
                        branchCounts[node]++;
                        parentLinks[mismatchBranch] = node;
                        childCounts[mismatchBranch] = mismatch.size();
                        childLinks[mismatchBranch] = mismatch;
                        for (size_t n : mismatch)
                        {
                            addSigToMatrix(mismatchBranch, means[n]);
                            parentLinks[n] = mismatchBranch;
                        }
                        childLinks[node].push_back(mismatchBranch);
                        updateNodeMean(mismatchBranch);
                    }
                    else
                    {
                        fprintf(stderr, "b\n");
                        // should just be 1 mismatch
                        for (size_t n : mismatch)
                        {
                            addSigToMatrix(node, means[n]);
                            childLinks[node].push_back(n);
                        }
                    }
                    updateNodeMean(node);
                }
                else
                {
                    fprintf(stderr, "b\n");
                    childCounts[node] = mismatch.size() + 1;
                    for (size_t n : mismatch)
                    {
                        addSigToMatrix(node, means[n]);
                        childLinks[node].push_back(n);
                    }
                }

                return newChild;
            }
        }
        else // NN with branch
        {
            // NN with leaves only => new node or new seq then merge all into new branch
            if (NN_branches.size() == 0)
            {
                // create new child
                size_t newChild = getNewNodeIdx(insertionList);
                fprintf(stderr, "created node %zu\n", newChild);
                seqIDs[newChild].push_back(idx);
                addSigToMatrix(newChild, signature);
                updateNodeMean(newChild);

                // NN to all leaves
                if (isBranchNode[node] && mismatch.size() + NN_leaves.size() == childCounts[node])
                {
                    fprintf(stderr, "match with all leaves %zu\n", node);
                    parentLinks[newChild] = node;
                    childCounts[node]++;
                    childLinks[node].push_back(newChild);
                    addSigToMatrix(node, means[newChild]);
                    updateNodeMean(node);
                }
                // some leaves not NN
                else
                {
                    // create new branch and let NN leaves be its children
                    size_t tempParent = getNewNodeIdx(insertionList);
                    fprintf(stderr, "3 created tempParent %zu\n", tempParent);
                    isBranchNode[tempParent] = 1;
                    branchCounts[node]++;
                    childCounts[tempParent] = NN_leaves.size() + 1;
                    for (size_t n : NN_leaves)
                    {
                        addSigToMatrix(tempParent, means[n]);
                        childLinks[tempParent].push_back(n);

                        // update NN children's parent to temp
                        parentLinks[n] = tempParent;
                    }

                    // add newChild to the back
                    childLinks[tempParent].push_back(newChild);
                    addSigToMatrix(tempParent, means[newChild]);
                    parentLinks[newChild] = tempParent;

                    updateNodeMean(tempParent);

                    // format node (previous parent), replace NN children with temp
                    vector<data_type> temp_matrice;
                    childCounts[node] = mismatch.size() + 1;
                    childLinks[node].clear();
                    for (size_t n : mismatch)
                    {
                        childLinks[node].push_back(n);
                        temp_matrice.push_back(means[n]);
                    }
                    childLinks[node].push_back(tempParent);
                    temp_matrice.push_back(means[tempParent]);
                    parentLinks[tempParent] = node;

                    //? check if vector replace is good
                    matrices[node] = temp_matrice;
                    updateNodeMean(node);
                }
                return newChild;
            }
            else if (NN_leaves.size() == 0)
            {
                // NN to 1 branch only, check its children
                if (NN_branches.size() == 1)
                {
                    fprintf(stderr, ">>>\n");
                    return traverse(signature, insertionList, idx, NN_branches[0]);
                }
                else
                {
                    // NN with multiple branches, no leaf
                    // create sibling
                    size_t newChild = createNode(signature, insertionList, node);
                    seqIDs[newChild].push_back(idx);
                    return newChild;
                }
            }
            else // mix case
            {

                // NN with everything, including branch
                if (mismatch.size() == 0)
                {
                    if (isBranchNode[node])
                    {
                        size_t newChild = createNode(signature, insertionList, node);
                        seqIDs[newChild].push_back(idx);
                        return newChild;
                    }
                    else
                    {
                        size_t tempParent = getNewNodeIdx(insertionList);
                        fprintf(stderr, "0 created tempParent %zu\n", tempParent);
                        isBranchNode[tempParent] = 1;
                        branchCounts[node]++;
                        parentLinks[tempParent] = node;
                        childCounts[tempParent] = childCounts[node];
                        childLinks[tempParent] = childLinks[node];

                        childCounts[node] = 1;
                        childLinks[node].clear();
                        childLinks[node].push_back(tempParent);

                        for (size_t n : childLinks[tempParent])
                        {
                            addSigToMatrix(tempParent, means[n]);
                            parentLinks[n] = tempParent;
                        }

                        size_t newChild = createNode(signature, insertionList, tempParent);
                        seqIDs[newChild].push_back(idx);
                        return newChild;
                    }
                }

                else if (NN_branches.size() == 1 && NN_leaves.size() == 1)
                {
                    fprintf(stderr, ">>>\n");
                    return traverse(signature, insertionList, idx, NN_branches[0]);
                }
                else
                {
                    // create new child
                    size_t newChild = getNewNodeIdx(insertionList);
                    fprintf(stderr, "created node %zu\n", newChild);
                    seqIDs[newChild].push_back(idx);
                    addSigToMatrix(newChild, signature);
                    updateNodeMean(newChild);

                    // create new branch and let NN leaves be its children
                    size_t tempParent = getNewNodeIdx(insertionList);
                    fprintf(stderr, "4 created tempParent %zu\n", tempParent);
                    isBranchNode[tempParent] = 1;
                    childCounts[tempParent] = NN_leaves.size() + NN_branches.size() + 1;
                    for (size_t n : NN_leaves)
                    {
                        addSigToMatrix(tempParent, means[n]);
                        childLinks[tempParent].push_back(n);
                        parentLinks[n] = tempParent;
                    }
                    for (size_t n : NN_branches)
                    {
                        addSigToMatrix(tempParent, means[n]);
                        childLinks[tempParent].push_back(n);
                        parentLinks[n] = tempParent;
                    }

                    // add newChild to the back
                    childLinks[tempParent].push_back(newChild);
                    addSigToMatrix(tempParent, means[newChild]);
                    parentLinks[newChild] = tempParent;

                    updateNodeMean(tempParent);

                    // format node (previous parent), replace NN children with temp
                    vector<data_type> temp_matrice;
                    childCounts[node] = mismatch.size() + 1;
                    childLinks[node].clear();
                    for (size_t n : mismatch)
                    {
                        childLinks[node].push_back(n);
                        temp_matrice.push_back(means[n]);
                    }
                    childLinks[node].push_back(tempParent);
                    temp_matrice.push_back(means[tempParent]);

                    //? check if vector replace is good
                    matrices[node] = temp_matrice;
                    updateNodeMean(node);

                    // // create new branch and let NN leaves be its children
                    // size_t tempParent = getNewNodeIdx(insertionList);
                    // isBranchNode[tempParent] = 1;
                    // childCounts[tempParent] = NN_leaves.size();
                    // for (size_t n : NN_leaves)
                    // {
                    //     addSigToMatrix(tempParent, means[n]);
                    //     childLinks[tempParent].push_back(n);

                    //     // update NN children's parent to temp
                    //     parentLinks[n] = tempParent;
                    // }

                    // updateNodeMean(tempParent);

                    // // format node (previous parent), replace NN children with temp
                    // vector<data_type> temp_matrice;
                    // childCounts[node] = NN_branches.size() + 2;
                    // childLinks[node].clear();
                    // for (size_t n : NN_branches)
                    // {
                    //     childLinks[node].push_back(n);
                    //     temp_matrice.push_back(means[n]);
                    // }
                    // childLinks[node].push_back(tempParent);
                    // temp_matrice.push_back(means[tempParent]);

                    // // add newChild to the back
                    // childLinks[node].push_back(newChild);
                    // addSigToMatrix(node, means[newChild]);
                    // parentLinks[newChild] = node;

                    // //? check if vector replace is good
                    // matrices[node] = temp_matrice;
                    // updateNodeMean(node);
                    fprintf(stderr, "Mix > %zu,%zu,%zu\n", node, tempParent, newChild);
                    return newChild;
                }
            }
        }
    }

    inline size_t createNode(data_type signature, vector<size_t> &insertionList, size_t parent)
    {
        // size_t parent = root;
        size_t node = getNewNodeIdx(insertionList);

        parentLinks[node] = parent;
        childLinks[parent].push_back(node);
        childCounts[parent]++;
        addSigToMatrix(node, signature);
        updateNodeMean(node);
        if (isBranchNode[parent])
        {
            addSigToMatrix(parent, signature);
            updateNodeMean(parent);
        }

        fprintf(stderr, "\ncreated node %zu\n", node);

        //? p
        // priority[node] = countSingleSetBits(means[node]);
        // priority[node] = 1;

        return node;
    }

    inline size_t first_insert(data_type signature, vector<size_t> &insertionList, size_t idx, size_t parent = 0)
    {
        size_t node = createNode(signature, insertionList, parent);
        seqIDs[node].push_back(idx);
        return node;
    }

    inline size_t insert(data_type signature, vector<size_t> &insertionList, size_t idx)
    {
        // bool stay = false;
        // size_t parent = root;
        // tie(stay, parent) = traverse(signature, insertionList);
        // fprintf(stderr, "inserting seq %zu at node %zu; stay: %d\n", idx, parent, stay);

        // if (!stay) // add new node
        // {
        //     size_t node = createNode(signature, insertionList, parent);
        //     seqIDs[node].push_back(idx);
        //     return node;
        // }
        // else // add seq without adding new node
        // {
        //     seqIDs[parent].push_back(idx);
        //     addSigToMatrix(parent, signature);
        //     means[parent] = createMeanSig(matrices[parent]);
        //     // priority[parent] = countSingleSetBits(means[parent]);
        //     //? p
        //     priority[parent] = calcDistortion(matrices[parent]);
        //     // priority[parent]++;

        //     // priority[parent]++;
        //     // return parent;
        // }

        size_t parent = traverse(signature, insertionList, idx);

        fprintf(stderr, "\ninserting %zu at %zu\n", idx, parent);
        return parent;
    }

    inline size_t search(data_type signature, size_t idx = 0)
    {
        size_t node = root;

        size_t best_child = node;
        double best_distance = numeric_limits<double>::max();

        while (childCounts[node] > 0)
        {
            size_t local_best_child = node;
            // double local_best_distance = 0;
            double local_best_distance = numeric_limits<double>::max();

            for (size_t i = 0; i < childCounts[node]; i++)
            {
                size_t child = childLinks[node][i];
                double distance = calcDistance(means[child], signature);

                fprintf(stderr, " <%zu,%.2f> ", child, distance);

                // found same, move on with the next seq
                if (distance <= (stay_threshold))
                {
                    // seqIDs[child].push_back(idx);
                    return child;
                }

                if (distance <= local_best_distance)
                {
                    local_best_distance = distance;
                    local_best_child = child;
                }
            }
            if (local_best_distance <= best_distance)
            {
                best_distance = local_best_distance;
                best_child = local_best_child;
            }

            node = local_best_child;
        }

        // seqIDs[best_child].push_back(idx);
        return best_child;
    }

    inline size_t searchBranch(data_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        return 0;
    }

    void clearSeqId(size_t lastNode)
    {
        for (size_t i = 0; i < lastNode; i++)
        {
            seqIDs[i].clear();
        }
    }

    // make sure cluster centroid is up-to-date before clearing matrices
    // clear matrices and seqID for the next insertion cycle
    // cluster means still remain
    void prepReinsert(size_t node)
    {

        updateNodeMean(node);
        matrices[node].clear();
        seqIDs[node].clear();

        if (childCounts[node] > 0)
        {
            for (size_t child : childLinks[node])
            {
                prepReinsert(child);
            }
        }
    }

    size_t reinsert(data_type signature, size_t idx)
    {
        size_t node = search(signature);
        seqIDs[node].push_back(idx);
        addSigToMatrix(node, signature);
        return node;
    }

    void destroyLocks(size_t node)
    {
        omp_destroy_lock(&locks[node]);
        if (childCounts[node] > 0)
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