#include <seqan3/search/views/partition_multi_hash.hpp>
#include <seqan3/search/views/partition_multi.hpp>
#include "minimiser.hpp"
#include "bloom_filter.hpp"
// #include "part_ktree.hpp"
#include <filesystem>
#include "build_partition_main_cmdline.hpp"

using namespace std;
namespace fs = std::filesystem;

static build_partition_main_cmdline args; // Command line switches and arguments
static uint8_t kmerLength = 9;            // Kmer length
static uint32_t windowLength = 50;        // window length
static uint32_t step_size = windowLength; // window length
size_t minimiser_size = 3;
size_t bf_element_cnt = 1000;

void doWork(bloom_parameters parameters, string inputFile)
{

    size_t temp = windowLength - kmerLength + 1;
    auto partition_view = seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmerLength}}) | seqan3::views::partition_multi(temp, kmerLength, minimiser_size, step_size);

    seqan3::sequence_file_input<dna4_traits> file_in{inputFile};

    bloom_filter bf(parameters);
    bloom_filter global_bf(parameters);
    size_t seqCount = 0;
    size_t indi_sum = 0;

    // Retrieve the sequences and ids.
    for (auto &[seq, id, qual] : file_in)
    {
        for (auto &&hashes : seq | partition_view)
        {
            for (size_t hash : hashes)
            {
                bf.insert(hash);
                global_bf.insert(hash);
            }
        }
        indi_sum += bf.setBits();
        seqCount++;
        bf.clear();
    }

    size_t total = global_bf.setBits();
    fprintf(stdout, "%zu,%zu,%zu,%.2f\n", seqCount, indi_sum, total, indi_sum * 1.0 / total);
}

int etr_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    std::ios::sync_with_stdio(false); // No sync with stdio -> faster

    if (args.kmer_given)
        kmerLength = args.kmer_arg;
    if (args.window_given)
        windowLength = args.window_arg;

    if (kmerLength > windowLength)
    {
        fprintf(stderr, "Error: kmer length must be smaller or equal to window length\n");
        return 1;
    }

    if (args.size_arg > windowLength - kmerLength + 1)
    {
        fprintf(stderr, "Error: number of minimisers per window must be smaller than %d (w - k + 1)\n", windowLength - kmerLength + 1);
        return 1;
    }

    if (args.element_given)
        bf_element_cnt = args.element_arg;

    bloom_parameters parameters;
    // How many elements roughly do we expect to insert?
    parameters.projected_element_count = bf_element_cnt;

    // Maximum tolerable false positive probability? (0,1)
    parameters.false_positive_probability = 0.001; // 1 in 10000

    // Simple randomizer (optional)
    parameters.random_seed = 0xA5A5A5A5;
    parameters.maximum_number_of_hashes = 1;

    if (!parameters)
    {
        std::cout << "Error - Invalid set of bloom filter parameters!" << std::endl;
        return 1;
    }
    parameters.compute_optimal_parameters();

    if (args.size_given)
    {
        minimiser_size = args.size_arg;
    }
    fprintf(stderr, "kmerLength,windowLength,minimiser_size,seqCount,indi_sum,total,etr\n");
    fprintf(stdout, "%u,%u,%zu,", kmerLength, windowLength, minimiser_size);

    string inputFile = args.input_arg;

    if (args.element_given)
    {
        bf_element_cnt = args.element_arg;
    }

    if (args.step_given)
    {
        step_size = args.step_arg;
    }
    else
    {
        step_size = windowLength;
    }

    if (args.folder_arg)
    {
        string delimiter = "/*";
        string folder = inputFile.substr(0, inputFile.find(delimiter));
        string ext = inputFile.substr(inputFile.find(delimiter) + delimiter.size(), inputFile.size() - 1);
        fprintf(stderr, "Reading folder %s\n", folder.c_str());

        bloom_filter bf(parameters);

        for (const auto &entry : fs::directory_iterator(folder))
        {
            if (entry.path().extension() == ext)
            {
                doWork(parameters, entry.path());
            }
        }
        return 0;
    }

    else
    {
        bloom_filter bf(parameters);
        doWork(parameters, inputFile);
    }
    return 0;
}
