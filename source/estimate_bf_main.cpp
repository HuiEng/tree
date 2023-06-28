#include <seqan3/search/views/partition_multi_hash.hpp>
#include <seqan3/search/views/partition_multi.hpp>
#include "minimiser.hpp"
#include <filesystem>
#include "build_partition_main_cmdline.hpp"

using namespace std;
namespace fs = std::filesystem;

static build_partition_main_cmdline args; // Command line switches and arguments
static uint8_t kmerLength = 9;            // Kmer length
static uint32_t windowLength = 50;        // window length
static uint32_t step_size = windowLength; // window length
size_t minimiser_size = 3;
bool multipleOut = false;

template <typename view>
void getMinimisers(view minimiser_view, string filename, set<size_t> &minimisers)
{
    seqan3::sequence_file_input<dna4_traits> file_in{filename};

    for (auto &[seq, id, qual] : file_in)
    {
        for (auto &&hashes : seq | minimiser_view)
        {
            for (size_t hash : hashes)
            {
                minimisers.insert(hash);
            }
        }
    }
}

// only consider the first window
template <typename view>
void getPartitionMinimisers(view minimiser_view, string filename, set<size_t> &minimisers)
{
    seqan3::sequence_file_input<dna4_traits> file_in{filename};
    for (auto &[seq, id, qual] : file_in)
    {
        for (auto &&hashes : seq | minimiser_view)
        {
            for (size_t hash : hashes)
            {
                minimisers.insert(hash);
            }
            break;
        }
    }
}

void estimateBFSize(string inputFile, set<size_t> &minimisers)
{
    size_t temp = windowLength - kmerLength + 1;
    auto partition_view = seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmerLength}}) | seqan3::views::partition_multi(temp, kmerLength, minimiser_size, step_size);

    if (args.toSingle_arg)
    {
        getMinimisers(partition_view, inputFile, minimisers);
    }
    else
    {
        getPartitionMinimisers(partition_view, inputFile, minimisers);
    }
}

int estimate_bf_main(int argc, char *argv[])
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

    if (args.size_given)
    {
        minimiser_size = args.size_arg;
    }
    fprintf(stderr, "Partition - Generating %zu minimisers per window...\n", minimiser_size);
    fprintf(stderr, "kmerLength= %u, windowLength = %u\n", kmerLength, windowLength);

    if (args.step_given)
    {
        step_size = args.step_arg;
    }
    else
    {
        step_size = windowLength;
    }

    set<size_t> minimisers;
    string inputFile = args.input_arg;
    if (args.folder_arg)
    {
        string delimiter = "/*";
        string folder = inputFile.substr(0, inputFile.find(delimiter));
        string ext = inputFile.substr(inputFile.find(delimiter) + delimiter.size(), inputFile.size() - 1);
        fprintf(stderr, "Reading folder %s\n", folder.c_str());

        for (const auto &entry : fs::directory_iterator(folder))
        {
            if (entry.path().extension() == ext)
            {
                estimateBFSize(entry.path(), minimisers);
            }
        }
    }
    else
    {
        estimateBFSize(inputFile, minimisers);
    }

    fprintf(stderr, "Estimate BF size = %zu \n", minimisers.size());
    return 0;
}
