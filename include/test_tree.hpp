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

    void destroyTree()
    {
        childCounts.clear();
        isBranchNode.clear();
        childLinks.clear();
        seqIDs.clear();
        parentLinks.clear();
        priority.clear();
        means.clear();
        matrices.clear();
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
        size_t i = 0;

        // fprintf(stderr, ">>%f,%f\n", meanSig.first, meanSig.second);
        for (data_type sig : matrices[node])
        {
            double distance = calcDistance(meanSig, sig);

            if (isBranchNode[node])
            {
                // fprintf(stderr, "%zu,%f,%f,%f\n", childLinks[node][i], sig.first, sig.second, distance);
            }
            else
            {
                // fprintf(stderr, "%zu,%f,%f,%f\n", seqIDs[node][i], sig.first, sig.second, distance);
            }
            if (distance > max_distance)
            {
                max_distance = distance;
            }
            i++;
        }

        // double max_priority = 0;
        // for (size_t child : childLinks[node])
        // {
        //     if (priority[child] > max_priority)
        //     {
        //         max_priority = priority[child];
        //     }
        // }
        // return max_distance + max_priority;
        // fprintf(stderr, "done %zu\n\n\n", node);
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
            fprintf(stderr, "ERROR deleting %zu from %zu!!\n", node, parent);
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
        means[new_node] = signature;

        // then add new node to t_parent
        childCounts[t_parent]++;
        childLinks[t_parent].push_back(new_node);
        addSigToMatrix(t_parent, means[new_node]);

        if (isBranchNode[t_parent])
        {
            updatePriority(t_parent);
        }

        fprintf(stderr, "create new node %zu\n", new_node);
        return new_node;
    }

    // return if match with branch centroid
    // return false if node is not branch
    // user set thresholds for leaves
    // thresholds for branch = threshold set by user + priority or distortion of children
    // if stay is leave, just add new seq in; if stay is branch, can choose to check it again
    // else return node, the vectors will be updated as well
    inline bool checkNode(size_t node, data_type signature, vector<size_t> &insertionList, vector<size_t> &mismatch, vector<size_t> &NN, vector<size_t> &stay)
    {
        // if (childCounts[node] > 5)
        // {
        //     forcesplitBranch(node, insertionList);
        // }

        while (splitBranch(node, insertionList))
        {
            // splitBranch(node, insertionList);
            // return checkNode(node, signature, insertionList, mismatch, NN, stay);
        }

        bool matchBranchCentroid = false;
        size_t i = 0;
        if (isBranchNode[node])
        {
            // skip first child later;
            i++;

            // first child is the centroid, deal with it separately
            size_t child = childLinks[node][0];
            fprintf(stderr, "checking %zu \n", child);
            splitBranch(child, insertionList);
            double local_stay_t = stay_threshold;

            // increase stay threshold of a node
            if (isBranchNode[child])
            {
                local_stay_t += priority[child];
            }

            double distance = calcDistance(means[child], signature);

            if (distance <= local_stay_t)
            {
                // do not need to merge NN, should have done so in previous operation
                fprintf(stderr, " stay in centroid %zu,%.2f \n", child, distance);
                matchBranchCentroid = true;
            }
        }

        // check the other children
        for (i; i < childCounts[node]; i++)
        {
            size_t child = childLinks[node][i];
            // for (size_t child : childLinks[node])
            // {
            double local_stay_t = stay_threshold;
            double local_split_t = split_threshold;

            // increase stay threshold of a node
            if (isBranchNode[child])
            {
                local_stay_t += priority[child];
                local_split_t += priority[child];
            }
            // else if (priority[child] > local_stay_t)
            // {
            //     local_stay_t = priority[child]; //?
            // }

            double distance = calcDistance(means[child], signature);

            if (distance <= local_stay_t)
            {
                fprintf(stderr, "<%zu,%.2f>: %f stay\n", child, distance, local_stay_t);
                stay.push_back(child);
                // dest = child;
            }

            else if (distance > local_split_t)
            {
                fprintf(stderr, "<%zu,%.2f>: %f miss\n", child, distance, local_split_t);
                mismatch.push_back(child);
            }
            else
            {
                fprintf(stderr, "<%zu,%.2f>: NN\n", child, distance);
                NN.push_back(child);
            }
        }

        return matchBranchCentroid;
    }

    // delete branch and move its chldren up
    void deleteBranch(size_t node)
    {
        size_t parent = parentLinks[node];

        if (isBranchNode[parent] && node == childLinks[parent][0])
        {
            //?
            fprintf(stderr, "??? need to delete parent branch\n");
        }
        for (size_t child : childLinks[node])
        {
            addSigToMatrix(parent, means[child]);
            childLinks[parent].push_back(child);
            childCounts[parent]++;
            deleteNode(child);
            parentLinks[child] = parent;
        }
        deleteNode(node);
        if (isBranchNode[parent])
        {
            fprintf(stderr, "update parent\n");
            updatePriority(parent);
            fprintf(stderr, "done\n");
        }

        fprintf(stderr, "parent %zu has %zu children\n", parent, childCounts[parent]);
        for (size_t child : childLinks[parent])
        {
            fprintf(stderr, ">%zu\n", child);
        }
    }

    inline size_t stayNode(data_type signature, vector<size_t> &insertionList, size_t idx, size_t node)
    {
        if (isBranchNode[node])
        {
            return tt(signature, insertionList, idx, node);
        }
        else
        {
            seqIDs[node].push_back(idx);
            addSigToMatrix(node, signature);
            updatePriority(node);
            return node;
        }
    }

    void moveParent(size_t child, size_t new_parent, bool d = true)
    {
        childCounts[new_parent]++;
        childLinks[new_parent].push_back(child);
        addSigToMatrix(new_parent, means[child]);
        if (d)
        {
            deleteNode(child);
        }
        parentLinks[child] = new_parent;
    }

    inline size_t tt_node(data_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        size_t current_childCount = childCounts[node];
        vector<size_t> mismatch;
        vector<size_t> NN;
        vector<size_t> stay;

        bool matchCentroid = checkNode(node, signature, insertionList, mismatch, NN, stay);
        size_t dest = node;

        fprintf(stderr, "tt_node %zu\n", node);
        if (mismatch.size() == current_childCount)
        {
            fprintf(stderr, "#mismatch\n");
            return createNode(signature, insertionList, node, idx);
        }

        if (stay.size() == 1)
        {
            fprintf(stderr, "#stay in one child\n", stay[0]);
            return stayNode(signature, insertionList, idx, stay[0]);
        }
        else if (stay.size() > 1)
        {
            // // what if there is branch
            // fprintf(stderr, "stay in multiple children %zu\n", stay.size());
            // size_t t_parent = createParent(node, insertionList);
            // dest = createNode(signature, insertionList, t_parent, idx);
            // means[t_parent] = signature;
            // addSigToMatrix(node, means[t_parent]);
            // for (size_t child : stay)
            // {
            //     moveParent(child, t_parent);
            // }
            // updatePriority(t_parent);
            // return dest;

            double min_distance = numeric_limits<double>::max();
            size_t best_child = 0;
            for (size_t child : stay)
            {
                double distance = calcDistance(means[child], signature);
                if (distance < min_distance)
                {
                    min_distance = distance;
                    best_child = child;
                }
            }
            if (isBranchNode[best_child])
            {
                fprintf(stderr, "best child is branch %zu\n", best_child);
            }
            return stayNode(signature, insertionList, idx, best_child);
        }

        if (NN.size() == 1)
        {
            fprintf(stderr, "??NN with one without stay\n");
            return createNode(signature, insertionList, node, idx);
        }
        else if (NN.size() > 1)
        {
            fprintf(stderr, "NN with multiple without stay\n");
            size_t t_parent = createParent(node, insertionList);
            dest = createNode(signature, insertionList, t_parent, idx);
            means[t_parent] = signature;
            addSigToMatrix(node, means[t_parent]);
            for (size_t child : NN)
            {
                moveParent(child, t_parent);
            }
            updatePriority(t_parent);
            return dest;
        }

        fprintf(stderr, "??wrong\n");
        return createNode(signature, insertionList, node, idx);
    }

    size_t splitBranch(size_t node, vector<size_t> &insertionList)
    {

        if (!isBranchNode[node])
        {
            return 0;
        }
        else if (priority[node] < split_node_threshold)
        {
            return 0;
        }
        fprintf(stderr, "splitbranch node %zu, %f\n", node, priority[node]);

        vector<size_t> temp_centroids = {childLinks[node][0]};
        vector<vector<size_t>> clusters;
        vector<size_t> temp;
        clusters.push_back(temp);
        vector<size_t> candidates_idx;
        for (size_t i = childCounts[node] - 1; i > 0; i--)
        {
            size_t child = childLinks[node][i];
            bool add = true;
            for (size_t centroid : temp_centroids)
            {
                double distance = calcDistance(means[centroid], matrices[node][i]);
                // fprintf(stderr, "%zu, %zu, %f\n", child, centroid, distance);
                if (distance < split_node_threshold)
                {
                    add = false;
                    break;
                }
            }
            if (add)
            {
                // fprintf(stderr, "----- %zu\n", child);
                temp_centroids.push_back(child);
                vector<size_t> t;
                clusters.push_back(t);
            }
            else
            {
                candidates_idx.push_back(i);
            }
        }

        if (temp_centroids.size() == 1)
        {
            fprintf(stderr, "??cannot splitbranch node %zu, %f\n", node, priority[node]);
            return 0;
        }

        for (size_t n : candidates_idx)
        {
            size_t child = childLinks[node][n];

            size_t dest = 0;
            double min_distance = priority[node];
            for (size_t i = 0; i < temp_centroids.size(); i++)
            {
                double distance = calcDistance(means[temp_centroids[i]], matrices[node][n]);
                // fprintf(stderr, "%zu, %zu, %f\n", child, temp_centroids[i], distance);
                if (distance <= min_distance)
                {
                    min_distance = distance;
                    dest = i;
                }
            }

            // fprintf(stderr, "--- %zu goes to %zu\n", child, dest);
            clusters[dest].push_back(child);
        }

        if (temp_centroids.size() == 2)
        {
            for (vector<size_t> cluster : clusters)
            {
                if (cluster.size() == 0)
                {
                    fprintf(stderr, ">>??something is wrong cannot split evenly\n");
                    return 0;
                }
            }
        }
        // prepare newBranch for valid clusters
        size_t parent = parentLinks[node];
        for (size_t i = 1; i < temp_centroids.size(); i++)
        {
            size_t centroid = temp_centroids[i];
            // fprintf(stderr, "*** %zu\n", centroid);
            if (isBranchNode[centroid])
            {
                moveParent(centroid, parent);
            }
            else if (clusters[i].size() > 0)
            {
                size_t t_parent = createParent(parent, insertionList);
                moveParent(centroid, t_parent);
                means[t_parent] = means[centroid];
                addSigToMatrix(parent, means[t_parent]);
            }
            else
            {
                // leaf & no NN
                moveParent(centroid, parent);
            }
        }

        // new centroid should have been move out at this stage
        for (size_t i = 1; i < temp_centroids.size(); i++)
        {
            size_t t_parent = parentLinks[temp_centroids[i]];
            for (size_t n = 0; n < clusters[i].size(); n++)
            {
                moveParent(clusters[i][n], t_parent);
            }
            updatePriority(t_parent);
        }

        updatePriority(node);
        if (isBranchNode[parent])
        {
            updatePriority(parent);
        }

        fprintf(stderr, ">> %f\n", priority[node]);
        return 1;
    }

    void prep()
    {
        for (size_t i = 0; i < childCounts[root]; i++)
        {
            size_t child = childLinks[root][i];
            // priority>split_threshold good enough?
            if (isBranchNode[child] && priority[child] > split_node_threshold)
            {
                size_t centroid = childLinks[child][0];
                if (isBranchNode[centroid])
                {
                    // dissolve this branch
                    childLinks[root][i] = centroid;
                    parentLinks[centroid] = root;

                    for (size_t j = 1; j < childCounts[child]; j++)
                    {
                        moveParent(childLinks[child][j], root, false);
                    }
                    clearNode(child); // doesn't do anything now, just clear memory; can insert back to insertionList
                }
                else
                {
                    for (size_t grandchild : childLinks[child])
                    {
                        if (isBranchNode[grandchild] && priority[grandchild] > split_node_threshold)
                        {
                            moveParent(grandchild, root);
                        }
                    }
                    updatePriority(child);
                }
            }
        }
    }

    size_t forcesplitBranch(size_t node, vector<size_t> &insertionList)
    {
        if (node != root)
        {
            return 0;
        }

        // printTreeJson(stderr);
        prep();
        // printTreeJson(stderr);

        fprintf(stderr, "force split node %zu, %f\n", node, priority[node]);
        vector<size_t> temp_centroids = {childLinks[node][0]};
        vector<vector<size_t>> clusters;
        vector<size_t> temp;
        clusters.push_back(temp);
        vector<size_t> t;
        clusters.push_back(t);

        double max_distance = 0;
        size_t candidate = 0;
        data_type mean0 = means[childLinks[node][0]];
        for (size_t i = childCounts[node] - 1; i > 0; i--)
        {
            double distance = calcDistance(mean0, matrices[node][i]);
            // fprintf(stderr, "%zu, %f\n", childLinks[node][i], distance);
            if (distance > max_distance)
            {
                max_distance = distance;
                candidate = childLinks[node][i];
            }
        }

        if (candidate == 0)
        {
            fprintf(stderr, "??something is wrong %zu, %f\n", node, priority[node]);
            return 0;
        }

        temp_centroids.push_back(candidate);
        // fprintf(stderr, ">??%zu, %zu\n", matrices[node].size(), childCounts[node]);

        // cluster the rest of the matrices, centroid will not be in clusters
        for (size_t n = 1; n < matrices[node].size(); n++)
        {
            size_t child = childLinks[node][n];
            if (child == candidate)
            {
                continue;
            }
            size_t dest = 0;
            double min_distance = numeric_limits<double>::max();
            for (size_t i = 0; i < temp_centroids.size(); i++)
            {
                double distance = calcDistance(means[temp_centroids[i]], matrices[node][n]);
                // fprintf(stderr, "%zu, %f\n", temp_centroids[i], distance);
                if (distance < min_distance)
                {
                    min_distance = distance;
                    dest = i;
                }
            }
            // fprintf(stderr, "--- %zu goes to %zu\n", child, dest);
            clusters[dest].push_back(child);
        }

        if (clusters[0].size() == 0 || clusters[1].size() == 0)
        {
            fprintf(stderr, "??something is wrong cannot split evenly\n");
            return 0;
        }

        // prepare newBranch for valid clusters

        // size_t parent = createParent(node, insertionList);
        // means[parent] = means[childLinks[node][0]];
        if (isBranchNode[node])
        {
            // promote the other clusters
            fprintf(stderr, "??promote\n");
            size_t parent = parentLinks[node];
            for (size_t i = 1; i < temp_centroids.size(); i++)
            {
                size_t centroid = temp_centroids[i];
                size_t t_parent = createParent(parent, insertionList);
                means[t_parent] = means[centroid];
                moveParent(centroid, t_parent);
                addSigToMatrix(parent, means[t_parent]);
            }
            // //?
            // for (size_t i = 1; i < temp_centroids.size(); i++)
            // {
            //     size_t centroid = temp_centroids[i];
            //     fprintf(stderr, "*** %zu\n", centroid);
            //     if (isBranchNode[centroid])
            //     {
            //         moveParent(centroid, parent);
            //     }
            //     else if (clusters[i].size() > 1)
            //     {
            //         size_t t_parent = createParent(parent, insertionList);
            //         moveParent(centroid, t_parent);
            //         means[t_parent] = means[centroid];
            //         addSigToMatrix(parent, means[t_parent]);
            //     }
            //     else
            //     {
            //         // leaf & no NN
            //         moveParent(centroid, parent);
            //     }
            // }
        }
        else
        {
            for (size_t i = 0; i < temp_centroids.size(); i++)
            {
                size_t centroid = temp_centroids[i];
                size_t t_parent = createParent(node, insertionList);
                means[t_parent] = means[centroid];
                moveParent(centroid, t_parent);
                addSigToMatrix(node, means[t_parent]);
            }
        }

        // new centroid should have been move out at this stage
        for (size_t i = 0; i < temp_centroids.size(); i++)
        {

            size_t centroid = temp_centroids[i];
            size_t t_parent = parentLinks[centroid];
            fprintf(stderr, "*** %zu\n", centroid);
            for (size_t n = 0; n < clusters[i].size(); n++)
            {
                moveParent(clusters[i][n], t_parent);
                fprintf(stderr, "--- %zu,%zu\n", clusters[i][n], t_parent);
            }
            updatePriority(t_parent);
        }

        updatePriority(node);

        // fprintf(stderr, ">?? done %zu, %zu\n", matrices[node].size(), childCounts[node]);

        // for (size_t i = 0; i < childCounts[root]; i++)
        // {
        //     if (matrices[root][i].first != means[childLinks[root][i]].first || matrices[root][i].second != means[childLinks[root][i]].second)
        //     {
        //         fprintf(stderr, "++ %zu, %zu\n", i, childLinks[root][i]);
        //         fprintf(stderr, "> %f, %f\n", matrices[root][i].first, matrices[root][i].second);
        //         fprintf(stderr, "* %f, %f\n\n", means[childLinks[root][i]].first, means[childLinks[root][i]].second);
        //     }
        // }

        return 1;
    }

    inline size_t tt_branch(data_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        size_t current_childCount = childCounts[node] - 1;
        vector<size_t> mismatch;
        vector<size_t> NN;
        vector<size_t> stay;

        bool matchCentroid = checkNode(node, signature, insertionList, mismatch, NN, stay);
        size_t dest = node;

        fprintf(stderr, "tt_branch %zu\n", node);
        if (matchCentroid)
        {
            fprintf(stderr, "#stay in centroid\n");
            // dest = stayNode(signature, insertionList, idx, childLinks[node][0]);

            if (mismatch.size() > 0)
            {
                fprintf(stderr, "has mismatch\n");
                size_t parent = parentLinks[node];
                for (size_t child : mismatch)
                {
                    fprintf(stderr, "moving %zu to %zu\n", child, parent);
                    moveParent(child, parent);
                }
                if (isBranchNode[parent])
                {
                    updatePriority(parent);
                }
            }

            return stayNode(signature, insertionList, idx, childLinks[node][0]);

            // if (st
        }
        else
        {

            if (current_childCount == 1)
            {
                if (mismatch.size() == 1)
                {
                    fprintf(stderr, "??wrong branch mismatch\n");
                    return createNode(signature, insertionList, parentLinks[node], idx);
                }
                else if (stay.size() == 1)
                {
                    fprintf(stderr, "breaking branch\n");
                    size_t parent = parentLinks[node];

                    for (size_t child : childLinks[node])
                    {
                        moveParent(child, parent, false);
                    }

                    deleteNode(node);
                    if (isBranchNode[parent])
                    {
                        updatePriority(parent);
                    }

                    return stayNode(signature, insertionList, idx, stay[0]);
                }
                else if (!isBranchNode[NN[0]])
                {
                    // only NN to centroid and the other child
                    return createNode(signature, insertionList, node, idx);
                    // fprintf(stderr, "??move centroid\n");
                    // // promote centroid to grandparent, make the input seq as the centroid of current branch instead
                    // size_t t_parent = parentLinks[node];
                    // moveParent(childLinks[node][0], t_parent, false);

                    // size_t new_node = getNewNodeIdx(insertionList);
                    // parentLinks[new_node] = node;
                    // seqIDs[new_node].push_back(idx);
                    // means[new_node] = signature;

                    // // correct node (branch)
                    // matrices[node][0] = signature;
                    // means[node] = signature;
                    // childLinks[node][0] = new_node;

                    // updatePriority(node);

                    // if (isBranchNode[t_parent])
                    // {
                    //     updatePriority(t_parent);
                    // }

                    // return new_node;

                    // // return tt_node(signature, insertionList, idx, NN[0]);
                    // fprintf(stderr, "??move centroid\n");
                    // size_t new_node = createNode(means[dest], insertionList, node, 0); // idx will be changed

                    // // copy centroid to new node
                    // matrices[new_node] = matrices[dest];
                    // seqIDs[new_node] = seqIDs[dest];
                    // priority[new_node] = priority[dest];

                    // // change mean
                    // means[dest] = signature;
                    // matrices[node][0] = signature;
                    // means[node] = signature;
                    // matrices[dest].clear();
                    // matrices[dest].push_back(signature);
                    // seqIDs[dest].clear();
                    // seqIDs[dest].push_back(idx);
                    // priority[dest] = 0;

                    // updatePriority(node);
                    // return dest;
                }
                else
                {
                    fprintf(stderr, "continue branch %zu\n", NN[0]);
                    return tt_branch(signature, insertionList, idx, NN[0]);
                }
            }

            if (mismatch.size() > 0)
            {
                // promote
                fprintf(stderr, "has mismatch\n");
                size_t parent = parentLinks[node];
                for (size_t child : mismatch)
                {
                    moveParent(child, parent);
                }
                if (isBranchNode[parent])
                {
                    updatePriority(parent);
                }
            }

            if (stay.size() == 1)
            {
                fprintf(stderr, "#b-stay in one child\n", stay[0]);
                return stayNode(signature, insertionList, idx, stay[0]);
            }
            else if (stay.size() > 1)
            {
                fprintf(stderr, "#b-stay in multiple children %zu\n", stay.size());

                double min_distance = numeric_limits<double>::max();
                size_t best_child = 0;
                for (size_t child : stay)
                {
                    double distance = calcDistance(means[child], signature);
                    if (distance < min_distance)
                    {
                        min_distance = distance;
                        best_child = child;
                    }
                }
                if (isBranchNode[best_child])
                {
                    fprintf(stderr, "best child is branch %zu\n", best_child);
                }
                return stayNode(signature, insertionList, idx, best_child);
                // // return createNode(signature, insertionList, node, idx);

                // // what if there is branch
                // size_t t_parent = createParent(node, insertionList);
                // dest = createNode(signature, insertionList, t_parent, idx);
                // means[t_parent] = signature;
                // addSigToMatrix(node, means[t_parent]);
                // for (size_t child : stay)
                // {
                //     moveParent(child, t_parent);
                // }
                // updatePriority(t_parent);

                // return dest;
            }
            else
            {
                fprintf(stderr, "??b replace mean of branch\n");
                return createNode(signature, insertionList, node, idx);
            }
        }

        fprintf(stderr, "??b deal with tt_branch later\n");
        return createNode(signature, insertionList, node, idx);
    }

    inline size_t tt(data_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        if (isBranchNode[node])
        {
            return tt_branch(signature, insertionList, idx, node);
        }
        else
        {
            return tt_node(signature, insertionList, idx, node);
        }
    }

    inline size_t first_insert(data_type signature, vector<size_t> &insertionList, size_t idx)
    {
        return createNode(signature, insertionList, root, idx);
    }

    inline size_t insert(data_type signature, vector<size_t> &insertionList, size_t idx)
    {
        // size_t node = traverse(signature, insertionList, idx);
        size_t node = tt(signature, insertionList, idx);

        fprintf(stderr, "inserted %zu at %zu\n\n", idx, node);
        for (size_t i = 0; i < childCounts[root]; i++)
        {

            if (matrices[root][i].first != means[childLinks[root][i]].first || matrices[root][i].second != means[childLinks[root][i]].second)
            {
                fprintf(stderr, "+ %zu, %zu\n", i, childLinks[root][i]);
                fprintf(stderr, "> %f, %f\n", matrices[root][i].first, matrices[root][i].second);
                fprintf(stderr, "* %f, %f\n\n", means[childLinks[root][i]].first, means[childLinks[root][i]].second);
            }
        }
        // printTreeJson(stderr);

        size_t parent = parentLinks[node];
        if (childCounts[parent] > 5)
        {
            forcesplitBranch(parent, insertionList);
        }
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