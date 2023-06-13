
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

void printMatrtix(sVec_type temp_matrix)
{
    printMsg("printing matrices\n");
    for (size_t i = 0; i < temp_matrix.size(); i += signatureSize)
    {
        toBinaryIdx(stderr, &temp_matrix[i]);
    }
    printMsg("done print\n");
}

sVec_type createMeanSig(const vector<cell_type> &clusterSigs)
{
    size_t seqCount = clusterSigs.size() / signatureSize;
    // if (seqCount==1){
    //     return clusterSigs;
    // }

    sVec_type meanSig(signatureSize);
    vector<int> unflattenedSignature(signatureSize * bits_per_char);
    size_t c = 0;
    for (size_t i = 0; i < seqCount; i++)
    {
        const cell_type *signatureData = &clusterSigs[i * signatureSize];
        for (size_t i = 0; i < unflattenedSignature.size(); i++)
        {
            cell_type signatureMask = (cell_type)1 << (i % bits_per_char);

            if (signatureMask & signatureData[i / bits_per_char])
            {
                unflattenedSignature[i] += 1;
                c++;
            }
            // else
            // {
            //     unflattenedSignature[i] -= 1;
            // }
        }
    }
    size_t d = 0;

    cell_type *flattenedSignature = &meanSig[0];
    for (size_t i = 0; i < unflattenedSignature.size(); i++)
    {
        if (unflattenedSignature[i] >= (seqCount + 1) / 2)
        {
            d++;
            flattenedSignature[i / bits_per_char] |= (cell_type)1 << (i % bits_per_char);
        }
    }

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
        sVec_type meanSig;
        if (isRootNode[node])
        {
            meanSig = unionNodeMean(node);
            // meanSig = createMeanSig(matrices[node]);
        }
        else if (isBranchNode[node] || isSuperNode[node])
        {
            meanSig = createMeanSig(getNonAmbiMatrix(node));
        }
        else
        {
            meanSig = createMeanSig(matrices[node]);
        }

        updateMeanSig(node, &meanSig[0]);
        updatePriority(node);
    }

    inline void updateMatrixIdx(size_t parent, size_t idx, size_t node)
    {
        memcpy(&matrices[parent][idx * signatureSize], getMeanSig(node), signatureWidth);
    }

    inline double calcAvgSim(sVec_type temp_matrix)
    {
        size_t seqCount = temp_matrix.size() / signatureSize;
        if (seqCount <= 1)
        {
            return 0;
        }

        double sumDistance = 0;
        sVec_type meanVec = createMeanSig(temp_matrix);
        s_type meanSig = &meanVec[0];

        for (size_t i = 0; i < temp_matrix.size(); i += signatureSize)
        {
            double distance = calcSimilarityWrap(meanSig, &temp_matrix[i]);
            sumDistance += distance;
        }
        return sumDistance / seqCount;
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

    double preloadPriority(size_t ambi, const_s_type signature)
    {
        // preload sig into matrices to check for priority
        sVec_type temp_matrix = matrices[ambi];
        temp_matrix.insert(temp_matrix.end(), signature, signature + signatureSize);

        return calcAvgSim(temp_matrix);
    }

};

#endif