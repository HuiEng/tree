
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

     double calcSimilarityWrap(s_type a, s_type b, size_t signatureSize_ = 0)
    {
        return calcSimilarity(a, b);
    }

    double calcSimilaritySigToNode(size_t node, sVec_type signatures, size_t i)
    {
        return calcSimilarityWrap(means[node], signatures[i]);
    }

    sVec_type createRandomSigs(size_t node, size_t clusterCount, size_t s)
    {
        // unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        unsigned seed = s;
        if (seed == 0)
        {
            seed = chrono::system_clock::now().time_since_epoch().count();
        }

        // find the smallest windows count
        vector<size_t> children = childLinks[node];
        size_t winNum = means[children[0]].size();
        for (size_t child : children)
        {
            s_type matrix = means[child];
            if (matrix.size() < winNum)
            {
                winNum = matrix.size();
            }
        }

        s_type randomSig(winNum); // = means[children[dist(rng)]];randomSig.resize(winNum);

        sVec_type clusterSigs(clusterCount);
        for (size_t i = 0; i < clusterCount; i++)
        {
            default_random_engine rng(seed + (i + 1) * 100);
            uniform_int_distribution<size_t> dist(0, children.size() - 1);
            for (size_t w = 0; w < winNum; w++)
            {
                size_t s = dist(rng);
                randomSig[w] = means[children[s]][w];
            }
            clusterSigs[i] = randomSig;
        }

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
    inline double calcAvgSim(sVec_type temp_matrix)
    {
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

    inline double calcAvgSim(size_t node)
    {
        sVec_type temp_matrix = matrices[node];
        if (isBranchNode[node])
        {
            temp_matrix = getNonAmbiMatrix(node);
        }
        
        return calcAvgSim(temp_matrix);
    }

    void clearMean(size_t node)
    {
        means[node].clear();
    }

    void printNodeDistance(FILE *stream, sVec_type seqs, vector<size_t> &clusters)
    {
        fprintf(stream, "seq_id,clu,HD\n");
        for (size_t i = 0; i < clusters.size(); i++)
        {
            size_t tnode = clusters[i];
            double distance = calcHD(getMeanSig(tnode), seqs[i]);
            fprintf(stream, "%zu,%zu,%.4f\n", i, tnode, distance);
        }
    }

    double preloadPriority(size_t ambi, const_s_type signature)
    {
        // preload sig into matrices to check for priority
        sVec_type temp_matrix = matrices[ambi];
        temp_matrix.push_back(signature);

        return calcAvgSim(temp_matrix);
    }
};

#endif