#ifndef INCLUDE_STATS_HPP
#define INCLUDE_STATS_HPP

#include <vector>
#include "bloom_filter.hpp"
#include "read.hpp"
#include "distance.hpp"

using namespace std;
size_t max_seqCount = 100;
size_t batch_size = 5;

class stats
{
public:
    double mean, stdev, q1, q3;
    bool isEmpty = true;

    stats() {}

    void update(double m, double s)
    {
        mean = m;
        stdev = s;
        q1 = mean - 0.675 * stdev;
        q3 = mean + 0.675 * stdev;
        isEmpty = false;
    }

    stats(double mean, double stdev)
    {
        update(mean, stdev);
    }

    void printStats()
    {
        size_t col = 4;
        string s(1 + 11 * col, '-');
        fprintf(stderr, "|%-10s|%-10s|%-10s|%-10s|\n", "mean", "stdev", "q1", "q3");
        fprintf(stderr, "%s\n", s.c_str());
        fprintf(stderr, "|%-10.2f|%-10.2f|%-10.2f|%-10.2f|\n", mean, stdev, q1, q3);
    }
};

vector<size_t> getIndices(size_t seqCount)
{
    vector<size_t> indices(seqCount);
    iota(indices.begin(), indices.end(), 0);
    if (seqCount > max_seqCount)
    {
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        shuffle(indices.begin(), indices.end(), default_random_engine(seed));
        indices.resize(max_seqCount);
        fprintf(stderr, "Random sampling %zu seqs...\n", max_seqCount);
    }

    return indices;
}

stats summarise(vector<double> const &v)
{
    stats st;
    if (v.empty())
    {
        return st;
    }

    auto const count = static_cast<double>(v.size());

    double sum = accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / count;

    vector<double> diff(count);
    transform(v.begin(), v.end(), diff.begin(),
              bind2nd(minus<double>(), mean));
    double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev = sqrt(sq_sum / count);
    st.update(mean, stdev);

    st.printStats();
    return st;
}

template <typename funcType>
void calcAllStatsBatch(vector<vector<seq_type>> seqs, size_t seqCount, funcType simFunc)
{
    vector<size_t> indices = getIndices(seqCount);
    size_t sample_size = indices.size();
    vector<double> data;

    for (size_t i = 0; i < sample_size; i++)
    {
        size_t idx = indices[i];
        seq_type seqA = seqs[floor(idx / batch_size)][idx % batch_size];
        for (size_t j = i + 1; j < sample_size; j++)
        {
            idx = indices[j];
            double sim = simFunc(seqA, seqs[floor(idx / batch_size)][idx % batch_size]) * 100.0;
            data.push_back(sim);
        }
    }
    summarise(data);
}

void calcAllStatsLocal(vector<seq_type> seqs)
{
    vector<size_t> indices = getIndices(seqs.size());
    size_t sample_size = indices.size();
    vector<double> data;

    for (size_t i = 0; i < sample_size; i++)
    {
        for (size_t j = i + 1; j < sample_size; j++)
        {
            double sim = calcMatchingWindows(seqs[indices[i]], seqs[indices[j]]) * 100.0;
            data.push_back(sim);
        }
    }
    summarise(data);
}

void calcAllStatsGlobal(vector<seq_type> seqs)
{
    vector<size_t> indices = getIndices(seqs.size());
    size_t sample_size = indices.size();
    vector<double> data;

    for (size_t i = 0; i < sample_size; i++)
    {
        for (size_t j = i + 1; j < sample_size; j++)
        {
            double sim = calcJaccardGlobal(seqs[indices[i]], seqs[indices[j]]) * 100;
            // fprintf(stdout, "%zu,%zu,%.2f\n", i, j, sim);
            data.push_back(sim);
        }
    }
    summarise(data);
}

void calcAllStats(vector<seq_type> seqs)
{
    vector<size_t> indices = getIndices(seqs.size());
    size_t sample_size = indices.size();
    vector<double> data;

    for (size_t i = 0; i < sample_size; i++)
    {
        for (size_t j = i + 1; j < sample_size; j++)
        {
            double sim = calcJaccardLocal(seqs[indices[i]], seqs[indices[j]]) * 100.0;
            data.push_back(sim);
        }
    }
    summarise(data);
}

// Jaccard stats
void calcAllStatsKmers(vector<cell_type> seqs)
{
    vector<size_t> indices = getIndices(seqs.size() / signatureSize);
    size_t sample_size = indices.size();
    vector<double> data;

    for (size_t i = 0; i < sample_size; i++)
    {
        size_t temp = countSetBits(&seqs[indices[i] * signatureSize], signatureSize);
        for (size_t j = i + 1; j < sample_size; j++)
        {
            double sim = calcSimilarity(&seqs[indices[i] * signatureSize], &seqs[indices[j] * signatureSize], signatureSize) * 100;
            data.push_back(sim);
        }
    }
    summarise(data);
}

// Jaccard stats
void calcAllStatsKmersBatch(vector<cell_type> seqs, size_t seqCount)
{
    vector<size_t> indices = getIndices(seqCount);
    size_t sample_size = indices.size();
    vector<double> data;

    for (size_t i = 0; i < sample_size; i++)
    {
        size_t temp = countSetBits(&seqs[indices[i] * signatureSize], signatureSize);
        for (size_t j = i + 1; j < sample_size; j++)
        {
            double sim = calcSimilarity(&seqs[indices[i] * signatureSize], &seqs[indices[j] * signatureSize], signatureSize) * 100;
            data.push_back(sim);
        }
    }
    summarise(data);
}

// void calcAllSetBits(const vector<cell_type> &sigs)
// {
//     size_t seqCount = sigs.size() / signatureSize;
//     for (size_t i = 0; i < seqCount; i++)
//     {
//         // fprintf(stderr,"%zu,%zu,%zu\n",kmerLength,windowLength,countSetBits(&sigs[i*signatureSize]));
//         fprintf(stderr, "BF density for the first seq: %zu\n", countSetBits(&sigs[i * signatureSize], signatureSize));
//         break;
//     }
// }

#endif