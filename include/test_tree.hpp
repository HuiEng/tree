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
    if (clusterSigs.size() <= 1)
    {
        return 0;
    }
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
        // fprintf(stream, "{\"node\":\"%zu\",\"branch\":\"%zu\",\"priority\":\"%.2f\",\"childCount\":\"%zu\",\"content\":\"*", tnode, isBranchNode[tnode], calcNodeMaxDistance(tnode), seqIDs[tnode].size());

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
            // fprintf(stream, "{\"node\":\"%zu\",\"priority\":\"%.2f\",\"childCount\":\"%zu\",\"content\":\"*", tnode, calcNodeMaxDistance(tnode), seqIDs[tnode].size());

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

    // add priority to the branch child
    double calcNodeDistortion(size_t node)
    {
        data_type meanSig = createMeanSig(matrices[node]);
        double sumSquareDistance = 0;

        for (size_t child : childLinks[node])
        {
            double distance = calcDistance(meanSig, means[child]);
            if (isBranchNode[child])
            {
                distance += priority[child];
            }
            sumSquareDistance += distance * distance;
        }
        return sqrt(sumSquareDistance / childCounts[node]);
    }

    // add priority to the branch child
    double calcNodeMaxDistance(size_t node)
    {
        data_type meanSig = means[node];
        double max_distance = 0;

        for (data_type sig : matrices[node])
        {
            double distance = calcDistance(meanSig, sig);
            if (distance > max_distance)
            {
                max_distance = distance;
            }
        }
        return max_distance;
    }

    // union mean of children
    inline void updatePriority(size_t node)
    {
        // priority[node] = calcDistortion(matrices[node]);
        priority[node] = calcNodeMaxDistance(node);

        // priority[node] = calcNodeDistortion(node);
    }

    // union mean of children
    inline void updateNodeMean(size_t node)
    {
        // fprintf(stderr, "updating node %zu\n", node);
        // means[node] = createMeanSig(matrices[node]);
        // //? p
        // updatePriority(node);

        if (isBranchNode[node])
        {
            fprintf(stderr, "updating node %zu\n", node);
            means[node] = createMeanSig(matrices[node]);
        }

        updatePriority(node);
        updateParentMean(node);
    }

    int getNodeIdx(size_t node)
    {
        size_t parent = parentLinks[node];
        int idx = -1;
        for (size_t i = 0; i < childCounts[parent]; i++)
        {
            if (childLinks[parent][i] == node)
            {
                idx = i;
                break;
            }
        }
        return idx;
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
        int idx = getNodeIdx(node);

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

    // delete a child from its parent, need to format the child separately
    inline void deleteNode(size_t node)
    {
        fprintf(stderr, "deleting %zu\n", node);
        size_t parent = parentLinks[node];

        int idx = getNodeIdx(node);

        if (idx == -1)
        {
            fprintf(stderr, "ERROR deleting %zu!!\n", parent);
        }
        else
        {
            matrices[parent].erase(matrices[parent].begin() + idx);
            childLinks[parent].erase(childLinks[parent].begin() + idx);
            childCounts[parent]--;
            // parentLinks[node] = 0;
            // isBranchNode[node] = 0;
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

    // haven't update childLinks[t_parent] yet
    inline size_t createParent(size_t node, vector<size_t> &insertionList)
    {
        size_t t_parent = getNewNodeIdx(insertionList);
        isBranchNode[t_parent] = 1;
        parentLinks[t_parent] = node;

        // update node
        childCounts[node]++;
        childLinks[node].push_back(t_parent);

        fprintf(stderr, "create t_parent %zu\n", t_parent);
        return t_parent;
    }

    inline size_t createNode(data_type signature, vector<size_t> &insertionList, size_t t_parent, size_t idx)
    {
        size_t new_node = getNewNodeIdx(insertionList);
        parentLinks[new_node] = t_parent;
        addSigToMatrix(new_node, signature);
        seqIDs[new_node].push_back(idx);
        // updateNodeMean(new_node);
        means[new_node] = signature;

        // then add new node to t_parent
        childCounts[t_parent]++;
        childLinks[t_parent].push_back(new_node);
        addSigToMatrix(t_parent, means[new_node]);

        // //? maybe call it outside
        // updateNodeMean(t_parent);

        if (isBranchNode[t_parent])
        {
            updatePriority(t_parent);
            // updateNodeMean(t_parent);
        }

        fprintf(stderr, "create new node %zu\n", new_node);
        return new_node;
    }

    // return stay child
    // if stay is leave, just add new seq in; if stay is branch, can choose to check it again
    // else return node, the vectors will be updated as well
    inline size_t checkNode(size_t node, data_type signature, vector<size_t> &mismatch, vector<size_t> &NN_branches, vector<size_t> &NN_leaves)
    {

        for (size_t child : childLinks[node])
        {
            double local_stay_t = stay_threshold;

            // increase stay threshold of a node
            if (isBranchNode[child])
            {
                local_stay_t += priority[child];
            }
            double local_split_t = local_stay_t * 4;

            double distance = calcDistance(means[child], signature);

            fprintf(stderr, " <%zu,%.2f> ", child, distance);

            if (distance <= local_stay_t)
            {
                return child;
            }

            if (distance > local_split_t)
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

        // get children of branches
        // candidates for reclustering
        vector<size_t> temp_nodes;
        for (size_t child : childLinks[node])
        {
            if (isBranchNode[child])
            {
                for (size_t n : childLinks[child])
                {
                    temp_nodes.push_back(n);
                }
                childCounts[child] = 0;
                childLinks[child].clear();
                matrices[child].clear();
            }
        }

        size_t changed = 0;

        // recluster
        vector temp_centroids = childLinks[node];
        for (size_t n : temp_nodes)
        {
            double best_dist = numeric_limits<double>::max();
            size_t best_centroid = 0;
            size_t best_centroid_idx = 0; // for updating node
            data_type signature = means[n];

            for (size_t i = 0; i < childCounts[node]; i++)
            {
                size_t child = childLinks[node][i];
                double distance = calcDistance(means[child], signature);

                fprintf(stderr, " <%zu,%.2f> ", child, distance);

                // found same, move on with the next seq
                if (distance <= best_dist)
                {
                    best_dist = distance;
                    best_centroid = child;
                    best_centroid_idx = i;
                }
            }

            fprintf(stderr, "\n");
            if (best_centroid != parentLinks[n])
            {
                fprintf(stderr, "child %zu from %zu to %zu\n", n, parentLinks[n], best_centroid);
                changed = 1;

                // create branch node if it is a leave
                if (!isBranchNode[best_centroid])
                {
                    size_t t_parent = getNewNodeIdx(insertionList);
                    childLinks[node][best_centroid_idx] = t_parent;

                    isBranchNode[t_parent] = 1;
                    parentLinks[t_parent] = node;
                    childCounts[t_parent] = 1;
                    childLinks[t_parent].push_back(best_centroid);
                    means[t_parent] = means[best_centroid];
                    addSigToMatrix(t_parent, means[best_centroid]);
                    parentLinks[best_centroid] = t_parent;
                    best_centroid = t_parent;
                }
                parentLinks[n] = best_centroid;
            }

            childCounts[best_centroid]++;
            childLinks[best_centroid].push_back(n);
            addSigToMatrix(best_centroid, signature);
        }

        for (size_t p : childLinks[node])
        {
            fprintf(stderr, "\n>>>%zu\n", p);
            for (size_t child : childLinks[p])
            {
                fprintf(stderr, "%zu\n", child);
            }
        }

        if (changed)
        {
            matrices[node].clear();
            for (size_t child : childLinks[node])
            {
                // updateNodeMean(child);
                means[child] = createMeanSig(matrices[child]);
                priority[child] = calcDistortion(matrices[child]);
                addSigToMatrix(node, means[child]);
            }
            updateNodeMean(node);
        }

        return changed;
    }

    // delete 1 branch, move all its children to its parent
    void splitNode(size_t node)
    {
        size_t parent = parentLinks[node];
        childCounts[parent] += childCounts[node];

        for (size_t child : childLinks[node])
        {
            childLinks[parent].push_back(child);
            addSigToMatrix(parent, means[child]);
            parentLinks[child] = parent;
        }

        deleteNode(node);

        if (isBranchNode[parent])
        {
            // updatePriority(parent);
            updateNodeMean(parent);
        }

        clearNode(node);
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
                // updateNodeMean(child);
                updatePriority(child);

                return child;
            }
        }

        // mismatch with everything, spawn new sibling
        if (mismatch.size() == current_childCount)
        {
            fprintf(stderr, ">>mismatch\n");
            size_t dest = createNode(signature, insertionList, node, idx);

            if (isBranchNode[node])
            {
                fprintf(stderr, ">>break cluster\n");
                splitNode(node);
            }
            return dest;
        }

        // proceed if only 1 branch matching
        if (NN_branches.size() == 1 && NN_leaves.size() == 0)
        {
            fprintf(stderr, "match 1 branch\n");
            return traverse(signature, insertionList, idx, NN_branches[0]);
        }

        size_t t_parent = node;
        if (isBranchNode[node])
        {
            // promote mismatches
            // makes mismatch the sibling of node
            size_t parent = parentLinks[node];

            if (mismatch.size() > 0)
            {
                fprintf(stderr, "promote mismatch from %zu to %zu\n", node, parent);
            }

            // add mismatches back to node
            childCounts[parent] += mismatch.size();
            for (size_t m : mismatch)
            {
                childLinks[parent].push_back(m);
                addSigToMatrix(parent, means[m]);
                parentLinks[m] = parent;
                fprintf(stderr, ">%zu\n", m);
            }

            // clear node, will be updated later
            childCounts[node] = 0;
            childLinks[node].clear();
            matrices[node].clear();

            //?
            // updateNodeMean(parent);
        }
        else
        {
            fprintf(stderr, "merging some NN to %zu\n", node);
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

        // matches with leaves only
        if (NN_branches.size() == 0)
        {
            childCounts[t_parent] += NN_leaves.size();
            childLinks[t_parent] = NN_leaves;
            for (size_t l : NN_leaves)
            {
                addSigToMatrix(t_parent, means[l]);
                parentLinks[l] = t_parent;
            }

            fprintf(stderr, "match leaves only %zu\n", t_parent);
            size_t dest = createNode(signature, insertionList, t_parent, idx);
            updateNodeMean(t_parent);
            updatePriority(node);
            return dest;
        }

        // add matches to t_parent
        childCounts[t_parent] += NN_branches.size() + NN_leaves.size();

        double best_dist = numeric_limits<double>::max();
        size_t best_node = 0;
        bool cont = false;

        for (size_t b : NN_branches)
        {
            childLinks[t_parent].push_back(b);
            addSigToMatrix(t_parent, means[b]);
            parentLinks[b] = t_parent;

            // get distance of branch child
            double distance = calcDistance(means[b], signature);

            if (distance < best_dist)
            {
                fprintf(stderr, "\n>b<%zu,%.2f>\n", b, distance);
                best_dist = distance;
                best_node = b;
                cont = false;
            }

            // get distance of grandchildren
            for (size_t grandC : childLinks[b])
            {
                double distance = calcDistance(means[grandC], signature);

                fprintf(stderr, " <%zu,%.2f> ", grandC, distance);

                if (distance < best_dist)
                {
                    best_dist = distance;
                    best_node = grandC;
                    cont = true;
                }

                if (distance <= stay_threshold)
                {
                    break;
                }
            }
        }

        for (size_t l : NN_leaves)
        {
            childLinks[t_parent].push_back(l);
            addSigToMatrix(t_parent, means[l]);
            parentLinks[l] = t_parent;
        }

        // updatePriority(t_parent);
        updateNodeMean(t_parent);
        updatePriority(node);

        //? promote mismatch grandchildren
        for (size_t b : NN_branches)
        {
            vector<size_t> mismatch_grandchildren;
            vector<size_t> match_grandchildren;
            // get distance of grandchildren
            for (size_t grandC : childLinks[b])
            {
                double distance = calcDistance(means[node], means[grandC]);
                fprintf(stderr, "\n>grand <%zu,%.2f>\n", grandC, distance);
                if (distance > split_threshold)
                {
                    mismatch_grandchildren.push_back(grandC);
                }
                else
                {
                    match_grandchildren.push_back(grandC);
                }
            }

            // if grandchild doesn not match with grandparent, promote
            if (mismatch_grandchildren.size() == childCounts[b])
            {
                fprintf(stderr, "promote mismatch branch and grandchildren %zu\n", b);
                deleteNode(b);
                size_t grandparent = parentLinks[node];
                parentLinks[b] = grandparent;
                childCounts[grandparent]++;
                childLinks[grandparent].push_back(b);
                addSigToMatrix(grandparent, means[b]);
            }
            else if (mismatch_grandchildren.size() > 0)
            {
                fprintf(stderr, "promoting some child from branch %zu\n", b);
                size_t mismatch_parent = getNewNodeIdx(insertionList);
                isBranchNode[mismatch_parent] = 1;
                childCounts[node]++;
                childLinks[node].push_back(mismatch_parent);
                parentLinks[mismatch_parent] = node;
                addSigToMatrix(node, means[mismatch_parent]);

                // // do this to avoid unipath, ignore for now
                // if (mismatch_grandchildren.size() > 1 ||
                //     (mismatch_grandchildren.size() == 1 && !isBranchNode[mismatch_grandchildren[0]]))
                // {
                // }

                // update mismatch_parent
                childCounts[mismatch_parent] += mismatch_grandchildren.size();
                childLinks[mismatch_parent] = mismatch_grandchildren;
                for (size_t grandC : mismatch_grandchildren)
                {
                    addSigToMatrix(mismatch_parent, means[grandC]);
                    parentLinks[grandC] = mismatch_parent;
                }
                updateNodeMean(mismatch_parent);

                // update current NN_branch
                childCounts[b] = match_grandchildren.size();
                childLinks[b] = match_grandchildren;
                matrices[b].clear();
                for (size_t grandC : match_grandchildren)
                {
                    addSigToMatrix(b, means[grandC]);
                }
                updateNodeMean(b);
            }
        }

        if (childCounts[node] == 0)
        {
            fprintf(stderr, "drop all children of node, and itself %zu\n", node);
            deleteNode(node);
            return traverse(signature, insertionList, idx, parentLinks[node]);
        }

        // this only possible for grandchildren
        if (best_dist <= stay_threshold)
        {
            fprintf(stderr, "stay in grandchild %zu\n", best_node);
            if (isBranchNode[best_node])
            {
                return traverse(signature, insertionList, idx, best_node);
            }
            else
            {
                seqIDs[best_node].push_back(idx);
                addSigToMatrix(best_node, signature);
                updatePriority(best_node);
                return best_node;
            }
        }
        else
        {
            // check if it it closer or NN leaves
            for (size_t l : NN_leaves)
            {
                double distance = calcDistance(means[l], signature);
                if (distance < best_dist)
                {
                    fprintf(stderr, "\n>l<%zu,%.2f>\n", l, distance);
                    best_dist = distance;
                    best_node = l;
                    cont = false;
                }
            }

            if (cont)
            {
                if (isBranchNode[best_node])
                {
                    fprintf(stderr, "cont. traverse child %zu\n", best_node);
                    return traverse(signature, insertionList, idx, best_node);
                }
                else
                {
                    // add new level
                    // copy content of the best node to a new node
                    size_t new_node = getNewNodeIdx(insertionList);
                    parentLinks[new_node] = best_node;
                    matrices[new_node] = matrices[best_node];
                    means[new_node] = means[best_node];
                    priority[new_node] = priority[best_node];
                    seqIDs[new_node] = seqIDs[best_node];

                    // make best node a branch
                    isBranchNode[best_node] = 1;
                    seqIDs[best_node].clear();
                    matrices[best_node].clear();
                    addSigToMatrix(best_node, means[best_node]);
                    childCounts[best_node] = 1;
                    childLinks[best_node].push_back(new_node);

                    // then add new node to t_parent
                    fprintf(stderr, "create sibling to leave %zu\n", new_node);
                    return createNode(signature, insertionList, best_node, idx);
                }
            }
            else
            {
                fprintf(stderr, "create sibling to %zu\n", best_node);
                return createNode(signature, insertionList, t_parent, idx);
            }
        }

        // // // prep new node then add to t_parent
        // // size_t dest = createNode(signature, insertionList, t_parent, idx);

        // // reclusterNode(t_parent, insertionList);

        // if (NN_branches.size() == 1)
        // {
        //     fprintf(stderr, "c0\n");
        //     // updateNodeMean(t_parent);
        //     return traverse(signature, insertionList, idx, NN_branches[0]);
        // }

        // double temp_priority = priority[t_parent];
        // if (!reclusterNode(t_parent, insertionList))
        // {
        //     fprintf(stderr, "c1\n");
        //     return createNode(signature, insertionList, t_parent, idx);
        // }
        // else
        // {
        //     fprintf(stderr, "c2\n");
        //     return traverse(signature, insertionList, idx, t_parent);
        // }
    }

    inline size_t tt(data_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        size_t current_childCount = childCounts[node];

        vector<size_t> mismatch;
        // near neighbour
        vector<size_t> NN_branches;
        vector<size_t> NN_leaves;

        size_t child = checkNode(node, signature, mismatch, NN_branches, NN_leaves);

        // stay
        if (child != node)
        {
            if (isBranchNode[child])
            {
                fprintf(stderr, "continue %zu\n", child);
                return tt(signature, insertionList, idx, child);
            }
            else
            {
                fprintf(stderr, "same %zu\n", child);
                seqIDs[child].push_back(idx);
                addSigToMatrix(child, signature);
                updatePriority(child);
                return child;
            }
        }

        // match none and split
        if (mismatch.size() == current_childCount)
        {
            return createNode(signature, insertionList, node, idx);
        }

        // match all
        if (mismatch.size() == 0)
        {
            size_t t_parent = node;
            if (!isBranchNode[node])
            {
                size_t t_parent = getNewNodeIdx(insertionList);
                isBranchNode[t_parent] = 1;
                parentLinks[t_parent] = node;
                size_t new_node = createNode(signature, insertionList, t_parent, idx);

                fprintf(stderr, "merging all to new parent %zu\n", t_parent);
                childCounts[t_parent] += childCounts[node];
                for (size_t n : childLinks[node])
                {
                    childLinks[t_parent].push_back(n);
                    addSigToMatrix(t_parent, means[n]);
                    parentLinks[n] = t_parent;
                }
                updatePriority(t_parent);
                means[t_parent] = signature;

                // format node
                childCounts[node] = 1;
                childLinks[node].clear();
                childLinks[node].push_back(t_parent);
                matrices[node].clear();
                matrices[node].push_back(means[t_parent]);
                updatePriority(node);
                return new_node;
            }
            fprintf(stderr, "match all, add to %zu\n", t_parent);
            return createNode(signature, insertionList, t_parent, idx);
        }

        // NN with leaves only
        if (NN_branches.size() == 0 && NN_leaves.size() > 1)
        {
            size_t t_parent = createParent(node, insertionList);
            size_t new_node = createNode(signature, insertionList, t_parent, idx);

            // means[t_parent] should be means[new_node]
            addSigToMatrix(node, means[t_parent]);

            fprintf(stderr, "merging NN leaves to %zu\n", t_parent);
            isBranchNode[t_parent] = 1;
            for (size_t leaf : NN_leaves)
            {
                deleteNode(leaf);
                childCounts[t_parent]++;
                childLinks[t_parent].push_back(leaf);
                addSigToMatrix(t_parent, means[leaf]);
                parentLinks[leaf] = t_parent;
            }
            updatePriority(t_parent);
            return t_parent;
        }

        // recluster
        // insert seq as new child of node
        // unpack grandchildren (chldren of NN branches)
        // see if any grandchildren nearer to the new seq
        if (NN_branches.size() > 0)
        {
        }

        return createNode(signature, insertionList, node, idx);
    }

    inline size_t insert(data_type signature, vector<size_t> &insertionList, size_t idx)
    {
        // size_t node = traverse(signature, insertionList, idx);
        size_t node = tt(signature, insertionList, idx);

        // //?
        // updateNodeMean(node);

        fprintf(stderr, "inserted %zu at %zu\n\n", idx, node);
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

    inline size_t findNearest(data_type signature, size_t node, size_t best_centroid, double best_dist = numeric_limits<double>::max())
    {
        size_t best_child = 0;
        for (size_t child : childLinks[node])
        {
            double distance = calcDistance(means[child], signature);
            fprintf(stderr, "%zu,%f\n", child, distance);
            if (distance < stay_threshold)
            {
                fprintf(stderr, "found stay %zu, isBranch %zu\n", child, isBranchNode[child]);
                if (isBranchNode[child])
                {
                    fprintf(stderr, "S cont. find %zu\n", best_child);
                    return findNearest(signature, child, child, distance);
                }
                return child;
            }

            if (distance < best_dist)
            {
                best_child = child;
                best_dist = distance;
                fprintf(stderr, ">>>%zu,%f\n", best_child, best_dist);
            }
        }

        if (best_child != 0)
        {
            if (isBranchNode[best_child])
            {
                fprintf(stderr, "cont. find %zu\n", best_child);
                return findNearest(signature, best_child, best_centroid, best_dist);
            }
            else
            {
                fprintf(stderr, "Best Child %zu\n", best_child);
                return best_child;
            }
        }

        fprintf(stderr, "Best centroid %zu\n", best_centroid);
        return best_centroid;
    }

    void nearest(FILE *stream, data_type signature, size_t node = 0)
    {

        fprintf(stream, "%zu,%.2f\n", node, calcDistance(means[node], signature));
        if (isBranchNode[node])
        {
            for (size_t child : childLinks[node])
            {
                nearest(stream, signature, child);
            }
        }
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