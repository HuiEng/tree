#ifndef INCLUDE_DISTANCE_HPP
#define INCLUDE_DISTANCE_HPP

#include <vector>
#include "bloom_filter.hpp"

typedef vector<vector<cell_type>> seq_type;
// typedef size_t distance_type;
typedef double distance_type;

double split_threshold = 1;
double stay_threshold = 0;
double split_node_threshold;

size_t minimiser_match_threshold = 4;
using namespace std;

enum
{
    STAY_F = 0,
    SPLIT_F = 1,
    NN_LEAVE_F = 2,
    NN_BRANCH_F = 3,
    NN_F = 4
};

size_t signatureSize = 0; // Signature size (depends on element in BF, obtained while read binary)

void toBinaryIdx(FILE *stream, vector<cell_type> sig)
{
    for (int i = 0; i < signatureSize; i++)
    {
        toBinaryIdxCell(stream, sig[i], i * bits_per_char);
    }
    fprintf(stream, "\n");
}

void toBinaryIdx(FILE *stream, const cell_type *sig)
{
    for (int i = 0; i < signatureSize; i++)
    {
        toBinaryIdxCell(stream, *(sig + i), i * bits_per_char);
    }
    fprintf(stream, "\n");
}

void toBinary(FILE *stream, cell_type *sig)
{
    //   fprintf(stderr, "%p: ", sig);
    for (size_t i = 0; i < signatureSize; i++)
    {
        int binary[bits_per_char];
        for (int n = 0; n < bits_per_char; n++)
            binary[bits_per_char - 1 - n] = (*(sig + i) >> n) & 1;

        for (int n = 0; n < bits_per_char; n++)
            fprintf(stream, "%d", binary[n]);
    }
    fprintf(stream, "\n");
}

// print each window in one line
void dbgPrintSignature(FILE *stream, seq_type seq)
{
    for (auto window : seq)
    {
        toBinary(stream, &window[0]);
    }
    fprintf(stream, "\n");
}

void dbgPrintSignatureIdx(FILE *stream, seq_type seq)
{
    for (auto window : seq)
    {
        toBinaryIdx(stream, window);
    }
    fprintf(stream, "\n");
}

size_t countSingleSetBits(vector<cell_type> seq)
{
    size_t c = 0;
    for (size_t i = 0; i < signatureSize; i++)
    {
        c += __builtin_popcountll(seq[i]);
    }
    return c;
}

size_t countSetBits(seq_type seq)
{
    size_t c = 0;
    // treat tail subseq as mismatch
    for (int w = 0; w < seq.size(); w++)
    {
        for (size_t i = 0; i < signatureSize; i++)
        {
            c += __builtin_popcountll(seq[w][i]);
        }
    }
    return c;
}

size_t calcHDWrap(seq_type shorter, seq_type longer)
{
    size_t c = 0;
    // treat tail subseq as mismatch
    for (int w = 0; w < shorter.size(); w++)
    {
        for (size_t i = 0; i < signatureSize; i++)
        {
            c += __builtin_popcountll(shorter[w][i] ^ longer[w][i]);
        }
    }
    // treat tail subseq as mismatch
    for (int w = shorter.size(); w < longer.size(); w++)
    {
        for (size_t i = 0; i < signatureSize; i++)
        {
            c += __builtin_popcountll(longer[w][i]);
        }
    }
    return c;
}

size_t calcHD(const cell_type *a, const cell_type *b)
{
    return calcHDBF(a, b, signatureSize);
}

size_t calcHD(seq_type a, seq_type b)
{
    if (a.size() < b.size())
    {
        return calcHDWrap(a, b);
    }
    else
    {
        return calcHDWrap(b, a);
    }
}

size_t calcSingleInter(vector<cell_type> a, vector<cell_type> b)
{
    size_t c = 0;
    for (size_t i = 0; i < signatureSize; i++)
    {
        c += __builtin_popcountll(a[i] & b[i]);
    }
    return c;
}

size_t calcSingleUnion(vector<cell_type> a, vector<cell_type> b)
{
    size_t c = 0;
    for (size_t i = 0; i < signatureSize; i++)
    {
        c += __builtin_popcountll(a[i] | b[i]);
    }
    return c;
}

size_t calcInter(seq_type a, seq_type b)
{
    size_t c = 0;
    // treat tail subseq as mismatch
    for (int w = 0; w < min(a.size(), b.size()); w++)
    {
        c += calcSingleInter(a[w], b[w]);
    }
    return c;
}

size_t calcUnion(seq_type a, seq_type b)
{
    size_t c = 0;
    seq_type shorter;
    seq_type longer;
    if (a.size() > b.size())
    {
        longer = a;
        shorter = b;
    }
    else
    {
        longer = b;
        shorter = a;
    }
    // treat tail subseq as mismatch
    for (int w = 0; w < shorter.size(); w++)
    {
        c += calcSingleUnion(longer[w], shorter[w]);
    }

    for (int w = shorter.size(); w < longer.size(); w++)
    {
        for (size_t i = 0; i < signatureSize; i++)
        {
            c += __builtin_popcountll(longer[w][i]);
        }
    }
    return c;
}

// matching window pair must have at least t minimisers in common
// return frequency of matching window
double calcMatchingWindows(seq_type a, seq_type b)
{
    double match = 0;
    // treat tail subseq as mismatch
    for (int w = 0; w < min(a.size(), b.size()); w++)
    {
        size_t c = calcSingleInter(a[w], b[w]);
        if (c >= minimiser_match_threshold)
        {
            match++;
        }
        else if (c > 0)
        {
            if (c == calcSingleUnion(a[w], b[w]))
            {
                match++;
            }
        }
    }
    return match / max(a.size(), b.size());
}

//?
double calcMatchingMinimisers(seq_type a, seq_type b)
{
    double size = min(countSetBits(a), countSetBits(b));
    return calcInter(a, b) / size;
}

double calcJaccard(seq_type a, seq_type b)
{
    // return calcInter(a, b) * 1.0 / calcUnion(a, b);

    // calc jaccard per window then average, assume extra windows has 0 jaccard similarity
    double jaccard = 0;
    // treat tail subseq as mismatch
    for (int w = 0; w < min(a.size(), b.size()); w++)
    {
        size_t inter_ = calcSingleInter(a[w], b[w]);
        size_t union_ = calcSingleUnion(a[w], b[w]);
        // fprintf(stderr, "%zu,%zu\n", inter_, union_);
        jaccard += inter_ * 1.0 / union_;
    }

    return jaccard / max(a.size(), b.size());
}

double calcJaccardLocalWrap(seq_type shorter, seq_type longer)
{
    // slide the shorter seq to find starting point
    size_t start = 0;
    for (start = 0; start < shorter.size(); start++)
    {
        if (calcSingleInter(shorter[start], longer[start]) >= stay_threshold)
        {
            break;
        }
    }
    double jaccard = 0;
    // treat tail subseq as mismatch
    for (int w = start; w < shorter.size(); w++)
    {
        size_t inter_ = calcSingleInter(shorter[w], longer[w]);
        size_t union_ = calcSingleUnion(shorter[w], longer[w]);
        jaccard += inter_ * 1.0 / union_;
    }

    return jaccard / (shorter.size() - start);
}

// if 1 seq is significantly shorter than the other one, cut start or tail of the longer seq
double calcJaccardLocal(seq_type a, seq_type b)
{
    if (a.size() < b.size())
    {
        return calcJaccardLocalWrap(a, b);
    }
    else
    {
        return calcJaccardLocalWrap(b, a);
    }
}

seq_type getMinimiseSet(seq_type a)
{
    seq_type output;
    if (a.size() > 0)
    {
        output.push_back(a[0]);
        for (int w = 1; w < a.size(); w++)
        {
            for (size_t i = 0; i < signatureSize; i++)
            {
                output[0][i] |= a[w][i];
            }
        }
    }
    return output;
}

double calcJaccardGlobal(seq_type a, seq_type b)
{
    seq_type x = getMinimiseSet(a);
    seq_type y = getMinimiseSet(b);
    return calcInter(x, y) * 1.0 / calcUnion(x, y);
}

double calcJaccardGlobal_cell(vector<cell_type> a, vector<cell_type> b)
{
    return calcSingleInter(a, b) * 1.0 / calcSingleUnion(a, b);
}

size_t calcPartitionBitsGlobal(seq_type a, seq_type b)
{
    seq_type x = getMinimiseSet(a);
    seq_type y = getMinimiseSet(b);
    return calcInter(x, y);
}

distance_type calcSimilarity(seq_type a, seq_type b)
{
    // if (a.size() < b.size())
    // {
    //     return calcHD(a, b);
    // }
    // else
    // {
    //     return calcHD(b, a);
    // }

    // return calcInter(a, b);
    // return calcJaccard(a,b);
    return calcJaccardGlobal(a, b);

    // return calcMatchingWindows(a, b);

    // normalised to density of b
    // return calcInter(a, b) * 1.0 / countSetBits(b);
}

// find the number of overlapping bits divide by the number of bits inthe refer_sig
double calcOverlap(seq_type query_sig, seq_type refer_sig)
{
    seq_type x = getMinimiseSet(query_sig);
    seq_type y = getMinimiseSet(refer_sig);
    return calcInter(x, y) * 1.0 / countSetBits(y);
}

// find the number of overlapping bits divide by the number of bits inthe refer_sig
double calcOverlapSingle(vector<cell_type> query_sig, vector<cell_type> refer_sig)
{
    return calcSingleInter(query_sig, refer_sig) * 1.0 / countSingleSetBits(refer_sig);
}

vector<cell_type> doSingleUnion(vector<cell_type> a, vector<cell_type> b)
{
    vector<cell_type> c = a;
    for (size_t i = 0; i < signatureSize; i++)
    {
        c[i] = a[i] | b[i];
    }
    return c;
}

seq_type doUnion(seq_type a, seq_type b)
{
    seq_type c;
    seq_type shorter;
    seq_type longer;
    if (a.size() > b.size())
    {
        longer = a;
        shorter = b;
    }
    else
    {
        longer = b;
        shorter = a;
    }
    // treat tail subseq as mismatch
    for (int w = 0; w < shorter.size(); w++)
    {
        c.push_back(doSingleUnion(longer[w], shorter[w]));
    }

    for (int w = shorter.size(); w < longer.size(); w++)
    {
        c.push_back(longer[w]);
    }
    return c;
}

#endif