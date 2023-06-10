
// There are two kinds of ktree nodes- branch nodes and leaf nodes
// Both contain a signature matrix, plus their own signature
// (the root node signature does not matter and can be blank)
// Branch nodes then contain 'order' links to other nodes
// Leaf nodes do not.
// However, as leaf nodes may become branch nodes, we allocate
// the space anyway.
// As the space to be used is determined at runtime, we use
// parallel arrays, not structs

#ifndef INCLUDE_primary_tree_HPP
#define INCLUDE_primary_tree_HPP
#include "tidy_tree.hpp"
// #include "temp_tree.hpp"

using namespace std;
size_t signatureWidth = signatureSize * sizeof(cell_type);

// typedef const cell_type *s_type;
typedef cell_type *s_type;
typedef const cell_type *const_s_type;
typedef vector<cell_type> sVec_type;

sVec_type createMeanSig(const vector<cell_type> &clusterSigs)
{
    sVec_type meanSig(signatureSize);
    size_t seqCount = clusterSigs.size() / signatureSize;

    if (seqCount == 1)
    {
    }
    vector<int> unflattenedSignature(signatureWidth);
    
    for (size_t i = 0; i < seqCount; i++)
    {
        const cell_type *signatureData = &clusterSigs[i * signatureSize];
        for (size_t i = 0; i < signatureWidth; i++)
        {
            cell_type signatureMask = (cell_type)1 << (i % bits_per_char);
            printMsg("%zu,");
            if (signatureMask & signatureData[i / bits_per_char])
            {
                unflattenedSignature[i] += 1;
            }
            else
            {
                unflattenedSignature[i] -= 1;
            }
        }
    }
    printMsg("\n>>>\n");
    cell_type *flattenedSignature = &meanSig[0];
    for (size_t i = 0; i < signatureWidth; i++)
    {
        if (unflattenedSignature[i] > 0)
        {
            printMsg("%zu,");
            flattenedSignature[i / bits_per_char] |= (cell_type)1 << (i % bits_per_char);
        }
    }
    printMsg("\n>>>\n");
    return meanSig;
}

// Derived class
class primary_tree : public tidy_tree<s_type, const_s_type, sVec_type>
{
public:
    using tidy_tree::tidy_tree;

    s_type getMeanSig(size_t node) { return &means[node * signatureSize]; }

    void updateMeanSig(size_t node, const_s_type signature)
    {
        memcpy(&means[node * signatureSize], signature, signatureWidth);
    }

    ///?
    void addSigToMatrix(size_t node, const_s_type signature)
    {
        matrices[node].insert(matrices[node].end(), signature, signature + signatureSize);
    }

    double calcSimilaritySigToNode(size_t node, sVec_type signatures, size_t i)
    {
        return calcSimilarityWrap(&means[node * signatureSize], &signatures[i * signatureSize]);
    }

    double calcSimilarityWrap(const_s_type a, const_s_type b, size_t signatureSize_ = 0)
    {
        return calcSimilarity(a, b, signatureSize);
    }

    sVec_type createRandomSigs(size_t node, size_t clusterCount, size_t s)
    {
        unsigned seed = s;
        if (seed == 0)
        {
            seed = chrono::system_clock::now().time_since_epoch().count();
        }
        default_random_engine rng(seed);

        sVec_type clusterSigs(signatureSize * clusterCount);
        sVec_type sigs = matrices[node];
        size_t signatureCount = childCounts[node]; // sigs.size() / signatureSize;
        uniform_int_distribution<size_t> dist(0, signatureCount - 1);
        bool finished = false;

        unordered_set<string> uniqueSigs;
        for (size_t i = 0; i < signatureCount; i++)
        {
            size_t sig = dist(rng);
            string sigData(signatureWidth, ' ');
            memcpy(&sigData[0], &sigs[sig * signatureSize], signatureWidth);
            uniqueSigs.insert(sigData);
            if (uniqueSigs.size() >= clusterCount)
            {
                finished = true;
                break;
            }
        }

        size_t i = 0;
        for (const auto &sig : uniqueSigs)
        {
            memcpy(&clusterSigs[i * signatureSize], sig.data(), signatureWidth);
            i++;
        }

        if (!finished)
        {
            if (uniqueSigs.size() != 1)
            {
                fprintf(stderr, "This should not happen\n");
                exit(1);
            }
            for (size_t i = 0; i < signatureSize; i++)
            {
                clusterSigs.push_back(clusterSigs[i]);
            }
        }

        return clusterSigs;
    }

    inline sVec_type getNonAmbiMatrix(size_t node)
    {
        sVec_type temp_matrix;
        for (size_t child : childLinks[node])
        {
            if (!isAmbiNode[child])
            {
                s_type sig = getMeanSig(child);
                temp_matrix.insert(temp_matrix.end(), sig, sig + signatureSize);
            }
        }

        return temp_matrix;
    }

    //?
    sVec_type unionNodeMean(size_t node)
    {
        sVec_type meanSig(signatureSize, 0);
        for (size_t child : childLinks[node])
        {
            for (size_t i = 0; i < signatureSize; i++)
            {
                meanSig[i] |= means[child * signatureSize + i];
            }
        }
        return meanSig;
    }

    // union mean of children
    inline void updateNodeMean(size_t node)
    {
        if (isRootNode[node])
        {
            updateMeanSig(node, &unionNodeMean(node)[0]);
        }
        else if (isBranchNode[node] || isSuperNode[node])
        {
            updateMeanSig(node, &createMeanSig(getNonAmbiMatrix(node))[0]);
        }
        else
        {
            updateMeanSig(node, &createMeanSig(matrices[node])[0]);
        }

        updatePriority(node);
    }

    inline void updateMatrixIdx(size_t parent, size_t idx, size_t node)
    {
        memcpy(&matrices[parent][idx * signatureSize], getMeanSig(node), signatureWidth);
    }

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
        s_type meanSig = &createMeanSig(temp_matrix)[0];

        for (size_t i = 0; i < temp_matrix.size(); i += signatureSize)
        {
            double distance = calcSimilarityWrap(meanSig, &temp_matrix[i]);
            sumDistance += distance;
        }

        return sumDistance / (temp_matrix.size() / signatureSize);
    }

    void clearMean(size_t node)
    {
        s_type sig = getMeanSig(node);
        fill(sig, sig + signatureSize, 0);
    }

    void printNodeDistance(FILE *stream, sVec_type seqs, vector<size_t> &clusters)
    {
        fprintf(stream, "seq_id,clu,HD\n");
        for (size_t i = 0; i < clusters.size(); i++)
        {
            size_t tnode = clusters[i];
            double distance = calcHD(getMeanSig(tnode), &seqs[i * signatureSize]);
            fprintf(stream, "%zu,%zu,%.4f\n", i, tnode, distance);
        }
    }
};

#endif