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

double split_threshold;
double stay_threshold;
size_t minimiser_match_threshold;
size_t partree_capacity;
double split_node_threshold;

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
    // return sqrt(sumSquareDistance);
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

    void printTreeJson(FILE *stream, size_t node = 0)
    {
        fprintf(stream, "var treeData = ");
        printSubTreeJson(stream, node);
        fprintf(stream, ";\n");
    }

    void printSubTreeMatrices(FILE *stream, const vector<data_type> &seqs, size_t node = 0)
    {
        if (isBranchNode[node])
        {
            for (data_type sig : matrices[node])
            {
                fprintf(stream, "%zu,%zu,%zu,%.2f,%.2f,%.2f\n", node, parentLinks[node], isBranchNode[node], priority[node], sig.first, sig.second);
            }

            for (size_t child : childLinks[node])
            {
                printSubTreeMatrices(stream, seqs, child);
            }
        }
        else if (node == root)
        {

            for (size_t child : childLinks[node])
            {
                printSubTreeMatrices(stream, seqs, child);
            }
        }
        else
        {
            for (size_t n : seqIDs[node])
            {
                data_type sig = seqs[n];
                fprintf(stream, "%zu,%zu,%zu,%.2f,%f,%f\n", node, parentLinks[node], isBranchNode[node], priority[node], sig.first, sig.second);
            }
        }
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

        if (isBranchNode[node])
        {
            updateParentMean(node);
        }
        // fprintf(stderr, "done\n");
    }

    // union mean of children
    inline size_t updateParentMean(size_t node)
    {
        fprintf(stderr, "updating %zu\n", node);
        size_t parent = parentLinks[node];
        // do nothing if parent is root
        if (parent == root)
        {
            return 0;
        }

        if (parent == node)
        {
            fprintf(stderr, "something is wrong\n");
            return 0;
        }
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
            fprintf(stderr, "ERROR updating %zu!!\n", parent);
        }
        else
        {
            fprintf(stderr, "here %zu %zu!!\n", node, parent);
            // updateNodeMean(node); //wrong
            matrices[parent][idx] == means[node];
            updateParentMean(parent);
        }
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

    // inline size_t forceSplit(size_t node, vector<size_t> &insertionList, double threshold)
    // {
    //     size_t childCount = childCounts[node];

    //     vector<vector<size_t>> clusters;
    //     // initialise cluster
    //     vector<size_t> temp = {childLinks[node][0]};
    //     clusters.push_back(temp);

    //     for (size_t i = 1; i < childCount; i++)
    //     {
    //         size_t child = childLinks[node][i];
    //         bool found = false;
    //         // for (int c = 0; c < clusters.size(); c++)
    //         for (vector<size_t> &cluster : clusters)
    //         {
    //             // vector<size_t> &cluster = clusters[c];
    //             for (size_t n : cluster)
    //             {
    //                 if (calcDistance(means[n], means[child]) < threshold)
    //                 {
    //                     found = true;
    //                     break;
    //                 }
    //             }
    //             if (found)
    //             {
    //                 cluster.push_back(child);
    //                 break;
    //             }
    //         }
    //         if (!found)
    //         {
    //             vector<size_t> t = {child};
    //             clusters.push_back(t);
    //         }
    //     }

    //     fprintf(stderr, "split %zu, %.2f\n", node, priority[node]);

    //     childLinks[node].clear();
    //     childCounts[node] = clusters.size();
    //     matrices[node].clear();
    //     for (vector<size_t> &cluster : clusters)
    //     {
    //         for (size_t child : cluster)
    //         {
    //             fprintf(stderr, "%zu,", child);
    //         }
    //         fprintf(stderr, "\n>\n");
    //         if (cluster.size() == 1)
    //         {
    //             // no need extra node for singleton, add it back to the node
    //             childLinks[node].push_back(cluster[0]);
    //             addSigToMatrix(node, means[cluster[0]]);

    //             // // may further split
    //             // if (isBranchNode[cluster[0]])
    //             // {
    //             //     splitNode(cluster[0], insertionList)
    //             // }
    //         }
    //         else
    //         {
    //             size_t newNode = getNewNodeIdx(insertionList);
    //             parentLinks[newNode] = node;

    //             childCounts[newNode] = cluster.size();
    //             childLinks[newNode] = cluster;
    //             for (size_t child : cluster)
    //             {
    //                 addSigToMatrix(newNode, means[child]);
    //             }
    //             updateNodeMean(newNode);
    //             isBranchNode[newNode] = 1;
    //             childLinks[node].push_back(newNode);
    //             addSigToMatrix(node, means[newNode]);
    //         }
    //     }
    //     updateNodeMean(node);
    //     fprintf(stderr, "nowF %zu, %.2f\n", node, priority[node]);

    //     return 1;
    // }

    // inline size_t splitNode(size_t node, vector<size_t> &insertionList)
    // {
    //     if (hasBranch(node))
    //     {
    //         return 0;
    //     }
    //     size_t childCount = childCounts[node];

    //     if (priority[node] < split_node_threshold | childCount <= 1)
    //     {
    //         return 0;
    //     }

    //     vector<vector<size_t>> clusters;
    //     // initialise cluster
    //     vector<size_t> temp = {childLinks[node][0]};
    //     clusters.push_back(temp);

    //     for (size_t i = 1; i < childCount; i++)
    //     {
    //         size_t child = childLinks[node][i];
    //         bool found = false;
    //         // for (int c = 0; c < clusters.size(); c++)
    //         for (vector<size_t> &cluster : clusters)
    //         {
    //             // vector<size_t> &cluster = clusters[c];
    //             for (size_t n : cluster)
    //             {
    //                 if (calcDistance(means[n], means[child]) < split_node_threshold)
    //                 {
    //                     found = true;
    //                     break;
    //                 }
    //             }
    //             if (found)
    //             {
    //                 cluster.push_back(child);
    //                 break;
    //             }
    //         }
    //         if (!found)
    //         {
    //             vector<size_t> t = {child};
    //             clusters.push_back(t);
    //         }
    //     }

    //     // no need to split
    //     if (clusters.size() == 1)
    //     {
    //         return 0;
    //     }

    //     fprintf(stderr, "split %zu, %.2f\n", node, priority[node]);

    //     childLinks[node].clear();
    //     childCounts[node] = clusters.size();
    //     matrices[node].clear();
    //     for (vector<size_t> &cluster : clusters)
    //     {
    //         for (size_t child : cluster)
    //         {
    //             fprintf(stderr, "%zu,", child);
    //         }
    //         fprintf(stderr, "\n>\n");
    //         if (cluster.size() == 1)
    //         {
    //             // no need extra node for singleton, add it back to the node
    //             childLinks[node].push_back(cluster[0]);
    //             addSigToMatrix(node, means[cluster[0]]);

    //             // // may further split
    //             // if (isBranchNode[cluster[0]])
    //             // {
    //             //     splitNode(cluster[0], insertionList)
    //             // }
    //         }
    //         else
    //         {
    //             size_t newNode = getNewNodeIdx(insertionList);
    //             isBranchNode[newNode] = 1;
    //             parentLinks[newNode] = node;

    //             childCounts[newNode] = cluster.size();
    //             childLinks[newNode] = cluster;
    //             for (size_t child : cluster)
    //             {
    //                 addSigToMatrix(newNode, means[child]);
    //             }
    //             updateNodeMean(newNode);
    //             childLinks[node].push_back(newNode);
    //             addSigToMatrix(node, means[newNode]);
    //         }
    //     }
    //     updateNodeMean(node);
    //     fprintf(stderr, "now %zu %.2f\n", node, priority[node]);
    //     return 1;
    // }

    // haven't update childLinks[t_parent] yet
    inline size_t createTempParent(size_t node, vector<size_t> &insertionList, bool copy = false)
    {
        size_t t_parent = getNewNodeIdx(insertionList);
        isBranchNode[t_parent] = 1;
        parentLinks[t_parent] = node;

        // copy everything from node, will be updated later
        if (copy)
        {
            childCounts[t_parent] = childCounts[node];
            childLinks[t_parent] = childLinks[node];
            matrices[t_parent] = matrices[node];
            means[t_parent] = means[node];
            priority[t_parent] = priority[node];

            for (size_t child : childLinks[t_parent])
            {
                parentLinks[child] = t_parent;
            }
        }

        // update node
        childCounts[node] = 1;
        childLinks[node].clear();
        childLinks[node].push_back(t_parent);
        matrices[node].clear();
        // node mean will be updated later once t_parent is ready
        addSigToMatrix(node, means[t_parent]);

        fprintf(stderr, "create t_parent %zu\n", t_parent);
        return t_parent;
    }

    inline size_t createNode(data_type signature, vector<size_t> &insertionList, size_t t_parent, size_t idx)
    {
        size_t new_node = getNewNodeIdx(insertionList);
        parentLinks[new_node] = t_parent;
        addSigToMatrix(new_node, signature);
        seqIDs[new_node].push_back(idx);
        updateNodeMean(new_node);

        // then add new node to t_parent
        childCounts[t_parent]++;
        childLinks[t_parent].push_back(new_node);
        addSigToMatrix(t_parent, means[new_node]);

        //? maybe call it outside
        updateNodeMean(t_parent);

        fprintf(stderr, "create new node %zu\n", new_node);
        return new_node;
    }

    // return stay child
    // if stay is leave, just add new seq in; if stay is branch, can choose to check it again
    // else return node, the vectors will be updated as well
    inline size_t checkNode(size_t node, data_type signature, vector<size_t> &mismatch, vector<size_t> &NN_branches, vector<size_t> &NN_leaves)
    {

        for (size_t i = 0; i < childCounts[node]; i++)
        {
            size_t child = childLinks[node][i];
            double distance = calcDistance(means[child], signature);

            fprintf(stderr, " <%zu,%.2f> ", child, distance);

            if (distance <= stay_threshold)
            {
                return child;
            }

            if (distance > split_threshold)
            {
                mismatch.push_back(child);
            }
            else if (isBranchNode[child])
            {
                // NN_branches.push_back(child);
                // split earlier if it's branch, change offset later
                if (distance > split_node_threshold)
                {
                    fprintf(stderr, " -*%zu*- ", child);
                    mismatch.push_back(child);
                }
                else
                {
                    fprintf(stderr, " -%zu- ", child);
                    NN_branches.push_back(child);
                }
            }
            else
            {
                NN_leaves.push_back(child);
            }
        }
        
        return node;
    }

    // unpack child branches and recluster again
    // create new t_parent if necessary
    size_t reclusterNode(size_t node, vector<size_t> &insertionList)
    {
        if (!hasBranch(node))
        {
            return 0;
        }
        fprintf(stderr, "reclustering %zu\n", node);
    }

    inline size_t traverse(data_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        size_t current_childCount = childCounts[node];

        vector<size_t> mismatch;
        // near neighbour
        vector<size_t> NN_branches;
        vector<size_t> NN_leaves;

        size_t child = checkNode(node, signature, mismatch, NN_branches, NN_leaves);

        if (child != node)
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

        // mismatch with everything, spawn new sibling
        if (mismatch.size() == current_childCount)
        {
            fprintf(stderr, ">>mismatch\n");
            return createNode(signature, insertionList, node, idx);
        }

        // // continue with the only NN branch, merge any NN leaves
        // if (NN_branches.size() == 1)
        // {
        //     fprintf(stderr, ">b\n");
        //     // continue with NN branch, create t_parent if there is any NN_leaves
        //     size_t NN_branch = NN_branches[0];

        //     // create t_parent to merge NN branch with NN leaves
        //     // maybe check distortion before doing this
        //     // eg if (NN_leaves.size() > 0 && priority[NN_branch] < split_node_threshold)
        //     if (NN_leaves.size() > 0)
        //     {
        //         size_t t_parent = node;

        //         // create supercluster if node is not already one
        //         if (!isBranchNode[node])
        //         {
        //             t_parent = createTempParent(node, insertionList);

        //             // add mismatches back to node
        //             childCounts[node] += mismatch.size();
        //             for (size_t m : mismatch)
        //             {
        //                 childLinks[node].push_back(m);
        //                 addSigToMatrix(node, means[m]);
        //             }
        //         }

        //         // add matches to t_parent
        //         childCounts[t_parent] += NN_leaves.size() + 1;
        //         childLinks[t_parent].push_back(NN_branch);
        //         addSigToMatrix(t_parent, means[NN_branch]);

        //         for (size_t l : NN_leaves)
        //         {
        //             childLinks[t_parent].push_back(l);
        //             addSigToMatrix(t_parent, means[l]);
        //         }

        //         updateNodeMean(t_parent);

        //         // //? may not need this here
        //         // updateNodeMean(node);
        //     }
        //     return traverse(signature, insertionList, idx, NN_branch);
        // }

        // matches all, add level => create new branch
        if (mismatch.size() == 0)
        {
            fprintf(stderr, ">>match all");
            if (current_childCount == 1 && hasBranch(node))
            {
                fprintf(stderr, ", traverse %zu\n", childLinks[node][0]);
                return traverse(signature, insertionList, idx, childLinks[node][0]);
            }
            else
            {
                size_t t_parent = node;
fprintf(stderr, ">>1\n");
                // if node alr a branch, use node
                // this is to avoid unipath
                if (!isBranchNode[node])
                {
                    fprintf(stderr, ">>2\n");
                    t_parent = createTempParent(node, insertionList, true);
                }
fprintf(stderr, ">>3\n");
                return createNode(signature, insertionList, t_parent, idx);
            }

            // // this node does not have branch child
            // if (NN_branches.size() == 0)
            // {
            //     fprintf(stderr, "a\n");
            //     size_t t_parent = node;

            //     // if node alr a branch, use node
            //     // this is to avoid unipath
            //     if (!isBranchNode[node])
            //     {
            //         t_parent = createTempParent(node, insertionList, true);
            //     }

            //     return createNode(signature, insertionList, t_parent, idx);
            // }
            // // this node does not have leaves, at least 2 branches and they are all NN
            // else if (NN_leaves.size() == 0)
            // {
            //     fprintf(stderr, "b\n");
            //     size_t t_parent = node;

            //     // create supercluster if node is not already one
            //     if (!isBranchNode[node])
            //     {
            //         t_parent = createTempParent(node, insertionList, true);
            //     }

            //     double temp_priority = priority[t_parent];
            //     reclusterNode(t_parent, insertionList);
            //     if (temp_priority == priority[t_parent])
            //     {
            //         fprintf(stderr, "b1\n");
            //         return createNode(signature, insertionList, t_parent, idx);
            //     }
            //     else
            //     {
            //         fprintf(stderr, "b2\n");
            //         return traverse(signature, insertionList, idx, t_parent);
            //     }
            // }
            // // mix case
            // else
            // {
            //     size_t t_parent = node;

            //     // create supercluster if node is not already one
            //     if (!isBranchNode[node])
            //     {
            //         t_parent = createTempParent(node, insertionList, true);
            //     }

            //     double temp_priority = priority[t_parent];
            //     reclusterNode(t_parent, insertionList);
            //     if (temp_priority == priority[t_parent])
            //     {
            //         fprintf(stderr, "c1\n");
            //         return createNode(signature, insertionList, t_parent, idx);
            //     }
            //     else
            //     {
            //         fprintf(stderr, "c2\n");
            //         return traverse(signature, insertionList, idx, t_parent);
            //     }
            // }
        }

        // add NN leaves to t_parent, leave mismatches in node
        // if (NN_branches.size() == 0)
        // {
        //     fprintf(stderr, ">>a\n");
        //     size_t t_parent = node;

        //     // create supercluster if node is not already one
        //     if (!isBranchNode[node])
        //     {
        //         t_parent = createTempParent(node, insertionList);

        //         // add mismatches back to node
        //         childCounts[node] += mismatch.size();
        //         for (size_t m : mismatch)
        //         {
        //             childLinks[node].push_back(m);
        //             addSigToMatrix(node, means[m]);
        //         }
        //     }

        //     // add matches to t_parent
        //     childCounts[t_parent] += NN_leaves.size();
        //     for (size_t l : NN_leaves)
        //     {
        //         childLinks[t_parent].push_back(l);
        //         addSigToMatrix(t_parent, means[l]);
        //     }

        //     // prep new node then add to t_parent
        //     size_t dest = createNode(signature, insertionList, t_parent, idx);

        //     // //? maybe don't need this here
        //     // reclusterNode(node, insertionList);

        //     return dest;
        // }
        // else //?
        // {
        //     fprintf(stderr, ">>c\n");
        //     // create t_parent to merge NN branch with NN leaves
        //     size_t t_parent = node;

        //     // create supercluster if node is not already one
        //     if (!isBranchNode[node])
        //     {
        //         t_parent = createTempParent(node, insertionList);

        //         // add mismatches back to node
        //         childCounts[node] += mismatch.size();
        //         for (size_t m : mismatch)
        //         {
        //             childLinks[node].push_back(m);
        //             addSigToMatrix(node, means[m]);
        //         }
        //     }

        //     // add matches to t_parent
        //     childCounts[t_parent] += NN_branches.size() + NN_leaves.size();
        //     for (size_t b : NN_branches)
        //     {
        //         childLinks[t_parent].push_back(b);
        //         addSigToMatrix(t_parent, means[b]);
        //     }
        //     for (size_t l : NN_leaves)
        //     {
        //         childLinks[t_parent].push_back(l);
        //         addSigToMatrix(t_parent, means[l]);
        //     }

        //     // prep new node then add to t_parent
        //     size_t dest = createNode(signature, insertionList, t_parent, idx);

        //     reclusterNode(t_parent, insertionList);

        //     return dest;
        // }

        // NN with some children
        size_t t_parent = node;

        if (isBranchNode[node])
        {
            // promote mismatches
            // makes mismatch the sibling of node
            size_t parent = parentLinks[node];
            fprintf(stderr, "promote %zu to %zu\n", node, parent);
            // add mismatches back to node
            childCounts[parent] += mismatch.size();
            for (size_t m : mismatch)
            {
                childLinks[parent].push_back(m);
                addSigToMatrix(parent, means[m]);
                fprintf(stderr, ">%zu\n", m);
            }

            // clear node, will be updated later
            childCounts[node] = 0;
            childLinks[node].clear();
            matrices[node].clear();

            //?
            updateNodeMean(parent);
        }
        else
        {
            fprintf(stderr, "merging some NN to %zu\n", t_parent);
            // create new branch to merge NN, mismatches stay in node
            t_parent = createTempParent(node, insertionList);
            // add mismatches back to node
            childCounts[node] += mismatch.size();
            for (size_t m : mismatch)
            {
                childLinks[node].push_back(m);
                addSigToMatrix(node, means[m]);
            }
        }

        // add matches to t_parent
        childCounts[t_parent] += NN_branches.size() + NN_leaves.size();
        for (size_t b : NN_branches)
        {
            childLinks[t_parent].push_back(b);
            addSigToMatrix(t_parent, means[b]);
        }
        for (size_t l : NN_leaves)
        {
            childLinks[t_parent].push_back(l);
            addSigToMatrix(t_parent, means[l]);
        }
        // // prep new node then add to t_parent
        // size_t dest = createNode(signature, insertionList, t_parent, idx);

        // reclusterNode(t_parent, insertionList);

        double temp_priority = priority[t_parent];
        reclusterNode(t_parent, insertionList);
        if (temp_priority == priority[t_parent])
        {
            fprintf(stderr, "c1\n");
            return createNode(signature, insertionList, t_parent, idx);
        }
        else
        {
            fprintf(stderr, "c2\n");
            return traverse(signature, insertionList, idx, t_parent);
        }
    }

    inline size_t insert(data_type signature, vector<size_t> &insertionList, size_t idx)
    {
        size_t node = traverse(signature, insertionList, idx);

        fprintf(stderr, "\ninserting %zu at %zu\n", idx, node);
        return node;
    }

    inline size_t superCluster(size_t cluster)
    {
        while (isBranchNode[parentLinks[cluster]])
        {
            cluster = parentLinks[cluster];
        }
        return cluster;
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