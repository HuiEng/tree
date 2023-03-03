// #include "bf_ktree.hpp"
#include <fstream>
#include <random>
#include "stats_main_cmdline.hpp"
#include "bloom_filter.hpp"
#include "read.hpp"
#include "self_tree.hpp"

#include <regex>
#include <vector>
#include <numeric>
#include <iostream>

using namespace std;

static stats_main_cmdline args; // Command line switches and arguments
// static size_t signatureSize;         // Signature size (depends on element in BF, obtained while read binary)
size_t max_seqCount = 100;

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

void calcAllStatsLocal(vector<seq_type> seqs)
{
    size_t seqCount = seqs.size() / signatureSize;
    vector<size_t> indices = getIndices(seqCount);
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

void calcAllSetBits(const vector<cell_type> &sigs)
{
    size_t seqCount = sigs.size() / signatureSize;
    for (size_t i = 0; i < seqCount; i++)
    {
        // fprintf(stderr,"%zu,%zu,%zu\n",kmerLength,windowLength,countSetBits(&sigs[i*signatureSize]));
        fprintf(stderr, "BF density for the first seq: %zu\n", countSetBits(&sigs[i * signatureSize], signatureSize));
        break;
    }
}

int stats_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    ios::sync_with_stdio(false); // No sync with stdio -> faster

    string bfIn = args.bf_input_arg;
    size_t bf_element_cnt = 1000;

    size_t firstindex = bfIn.find_last_of("/") + 1;
    size_t lastindex = bfIn.find_last_of(".");
    string rawname = bfIn.substr(firstindex, lastindex - firstindex);

    if (args.output_given)
    {
        rawname = rawname + "-" + args.output_arg;
    }

    if (args.max_given)
    {
        max_seqCount = args.max_arg;
    }

    cmatch matches;
    if (regex_search(args.bf_input_arg, matches, regex("-b([0-9]+)")))
    {
        stringstream sstream(matches[1]);
        sstream >> bf_element_cnt;
    }

    if (args.all_kmer_arg)
    {
        vector<cell_type> seqs;
        signatureSize = readSignatures(bfIn, seqs);

        fprintf(stderr, "Loaded %zu seqs...\n", seqs.size() / signatureSize);
        calcAllStatsKmers(seqs);
    }
    else
    {
        vector<seq_type> seqs = readPartitionBF(bfIn);
        // dummy code, assume there is at least 10 input seqs
        for (int i = 0; i < 10; i++)
        {
            if (seqs[i].size() > 0)
            {
                signatureSize = seqs[i][0].size();
                break;
            }
        }
        if (signatureSize == 0)
        {
            fprintf(stderr, "Something is wrong with the input data, please generate signature with diff params\n");
            return 0;
        }

        // fprintf(stderr,"done\n" );
        fprintf(stderr, "Loaded %zu seqs...\n", seqs.size());

        if (args.local_arg)
        {
            if (args.threshold_given)
            {

                minimiser_match_threshold = args.threshold_arg;
            }
            fprintf(stderr, "Matching at least %zu minimisers per window...\n", minimiser_match_threshold);
            calcAllStatsLocal(seqs);
        }
        else if (args.global_arg)
        {
            calcAllStatsGlobal(seqs);
        }
        else
        {
            calcAllStats(seqs);
        }
    }

    return 0;
}
