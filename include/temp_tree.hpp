
// There are two kinds of ktree nodes- branch nodes and leaf nodes
// Both contain a signature matrix, plus their own signature
// (the root node signature does not matter and can be blank)
// Branch nodes then contain 'order' links to other nodes
// Leaf nodes do not.
// However, as leaf nodes may become branch nodes, we allocate
// the space anyway.
// As the space to be used is determined at runtime, we use
// parallel arrays, not structs

#ifndef INCLUDE_temp_tree_HPP
#define INCLUDE_temp_tree_HPP
#include "tidy_tree.hpp"

typedef seq_type s_type;
typedef seq_type const_s_type;
typedef vector<seq_type> sVec_type;

s_type createMeanSig(sVec_type clusterSigs)
{
    s_type meanSig;
    // find the smallest windows count
    size_t winNum = clusterSigs[0].size();
    for (s_type matrix : clusterSigs)
    {
        if (matrix.size() < winNum)
        {
            winNum = matrix.size();
        }
    }
    vector<vector<size_t>> counters;
    for (size_t w = 0; w < winNum; w++)
    {
        vector<size_t> counter(signatureSize * bits_per_char);
        counters.push_back(counter);

        vector<cell_type> temp(signatureSize * bits_per_char);
        meanSig.push_back(temp);
    }

    for (s_type signatureData : clusterSigs)
    {
        for (size_t w = 0; w < winNum; w++)
        {
            for (size_t i = 0; i < signatureSize; i++)
            {
                for (int n = 0; n < bits_per_char; n++)
                {
                    if ((signatureData[w][i] >> n) & 1)
                    {
                        counters[w][i * bits_per_char + n]++;
                    }
                }
            }
        }
    }
    for (size_t w = 0; w < winNum; w++)
    {
        vector<size_t> counter = counters[w];
        for (int i = 0; i < counter.size(); i++)
        {
            // make it upperbound
            if (counter[i] >= (clusterSigs.size() + 1) / 2)
            {
                meanSig[w][i / bits_per_char] |= (cell_type)1 << (i % bits_per_char);
            }
        }
    }
    return meanSig;
}

// Derived class
class temp_tree : public tidy_tree<s_type, const_s_type, sVec_type>
{
public:
    using tidy_tree::tidy_tree;

    s_type getMeanSig(size_t node) { return means[node]; }

    void updateMeanSig(size_t node, const_s_type signature)
    {
        means[node] = signature;
    }

    void addSigToMatrix(size_t node, const_s_type signature)
    {
        matrices[node].push_back(signature);
    }

    double calcSimilaritySigToNode(size_t node, sVec_type signatures, size_t i)
    {
        return calcSimilarityWrap(means[node], signatures[i]);
    }

    double calcSimilarityWrap(s_type a, s_type b, size_t signatureSize_ = 0)
    {
        return calcSimilarity(a, b);
    }

    sVec_type createRandomSigs(size_t node, size_t clusterCount, size_t s)
    {
        unsigned seed = s;
        if (seed == 0)
        {
            seed = chrono::system_clock::now().time_since_epoch().count();
        }
        default_random_engine rng(seed);

        sVec_type clusterSigs(clusterCount);
        
        return clusterSigs;
    }

    inline sVec_type getNonAmbiMatrix(size_t node)
    {
        sVec_type temp_matrix;
        for (size_t c = 0; c < childCounts[node]; c++)
        {
            size_t child = childLinks[node][c];
            // skip ambi
            if (!isAmbiNode[child])
            {
                temp_matrix.push_back(getMeanSig(child));
                // temp_matrix.push_back(matrices[node][c]);
            }
        }

        return temp_matrix;
    }

    ///???
    // union mean of children
    inline void updateNodeMean(size_t node)
    {
        if (isRootNode[node])
        {
            const_s_type meanSig = getMeanSig(childLinks[node][0]);

            for (size_t child : childLinks[node])
            {
                meanSig = doUnion(meanSig, getMeanSig(child));
            }
            updateMeanSig(node, meanSig);
        }
        else if (isBranchNode[node] || isSuperNode[node])
        {
            updateMeanSig(node, createMeanSig(getNonAmbiMatrix(node)));
        }
        else
        {
            updateMeanSig(node, createMeanSig(matrices[node]));
        }
        updatePriority(node);
    }

    inline void updateMatrixIdx(size_t parent, size_t idx, size_t node)
    {
        matrices[parent][idx] = getMeanSig(node);
    }

    //?
    inline double calcAvgSim(size_t node)
    {
        sVec_type temp_matrix = matrices[node];
        if (isBranchNode[node])
        {
            temp_matrix = getNonAmbiMatrix(node);
        }
        if (temp_matrix.size() <= 1)
        {
            return 0;
        }

        double sumDistance = 0;
        s_type meanSig = createMeanSig(temp_matrix);

        for (s_type signature : temp_matrix)
        {
            double distance = calcSimilarityWrap(meanSig, signature);
            sumDistance += distance;
        }

        return sumDistance / temp_matrix.size();
    }

    void clearMean(size_t node)
    {
        means[node].clear();
    }
};

#endif