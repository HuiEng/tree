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
bool debug = false;
bool compressReads = false;
bool multipleOut = false;
string outputFile = "";

template <typename view>
void getMinimisers(view minimiser_view, bloom_parameters parameters, string filename, FILE *wf)
{
    seqan3::sequence_file_input<dna4_traits> file_in{filename};
    size_t max = 8;
    size_t end = -1;
    size_t i = 0;

    bloom_filter bf(parameters);

    if (compressReads)
    {
        // Retrieve the sequences and ids.
        for (auto &[seq, id, qual] : file_in)
        {
            for (auto &&hashes : seq | minimiser_view)
            {
                for (size_t hash : hashes)
                {
                    bf.insert(hash);
                }
            }
        }
        fprintf(wf, "%zu,%zu\n", i, bf.setBits());
        i++;
        bf.clear();
    }
    else
    {
        // Retrieve the sequences and ids.
        for (auto &[seq, id, qual] : file_in)
        {
            for (auto &&hashes : seq | minimiser_view)
            {
                for (size_t hash : hashes)
                {
                    bf.insert(hash);
                }
            }
            fprintf(wf, "%zu,%zu\n", i, bf.setBits());
            i++;
            bf.clear();
        }
    }
}

template <typename view>
void getPartitionMinimisers(view minimiser_view, bloom_parameters parameters, string filename, FILE *wf)
{
    seqan3::sequence_file_input<dna4_traits> file_in{filename};
    size_t max = 8;
    size_t end = -1;

    bloom_filter bf(parameters);
    size_t i = 0;

    if (compressReads)
    {
        auto first_seq = (*file_in.begin()).sequence();
        auto windows = first_seq | minimiser_view;
        vector<bloom_filter> wBFL;
        for (auto &&hashes : windows)
        {
            bloom_filter temp(parameters);
            wBFL.push_back(temp);
        }

        for (auto &[seq, id, qual] : file_in)
        {
            size_t i = 0;
            // fprintf(stdout, ">\n");
            for (auto &&hashes : seq | minimiser_view)
            {
                for (size_t hash : hashes)
                {
                    wBFL[i].insert(hash);
                }
                i++;
            }
        }

        // print the list
        size_t s_bit = 0;
        for (bloom_filter wbf : wBFL)
        {
            s_bit += wbf.setBits();
        }
        fprintf(wf, "%zu,%zu\n", i, s_bit);
        i++;
    }

    else
    {

        // Retrieve the sequences and ids.
        for (auto &[seq, id, qual] : file_in)
        {
            size_t s_bit = 0;
            for (auto &&hashes : seq | minimiser_view)
            {
                for (size_t hash : hashes)
                {
                    bf.insert(hash);
                }
                s_bit += bf.setBits();
                bf.clear();
            }
            fprintf(wf, "%zu,%zu\n", i, s_bit);
            i++;
        }
    }
    // wf.close();
}

void doWork(FILE *wf, bloom_parameters parameters, string inputFile)
{

    if (args.canonical_arg)
    {
        auto partition_view = seqan3::views::partition_multi_hash(seqan3::shape{seqan3::ungapped{kmerLength}},
                                                                  seqan3::window_size{windowLength},
                                                                  seqan3::minimiser_size{minimiser_size},
                                                                  seqan3::step_size{step_size},
                                                                  seqan3::seed{0});
        getPartitionMinimisers(partition_view, parameters, inputFile, wf);
    }
    else
    {
        // to get minimisers with w=8,k=4
        // input param for the minimiser view is calculated by: window size - k-mer size + 1, here: 8 - 4 + 1 = 5)
        size_t temp = windowLength - kmerLength + 1;
        auto partition_view = seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmerLength}}) | seqan3::views::partition_multi(temp, kmerLength, minimiser_size, step_size);

        if (args.toSingle_arg)
        {
            getMinimisers(partition_view, parameters, inputFile, wf);
        }
        else
        {
            getPartitionMinimisers(partition_view, parameters, inputFile, wf);
        }
    }

    // wf.close();
}

int etr_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    std::ios::sync_with_stdio(false); // No sync with stdio -> faster

    //

    debug = args.debug;

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

    compressReads = args.compress_arg;
    if (compressReads)
    {
        std::cout << "Compressing Reads" << std::endl;
    }
    multipleOut = args.multiple_arg;
    if (multipleOut)
    {
        std::cout << "multipleOut" << std::endl;
    }

    if (args.size_given)
    {
        minimiser_size = args.size_arg;
    }
    fprintf(stderr, "Partition - Generating %zu minimisers per window...\n", minimiser_size);
    fprintf(stderr, "kmerLength= %u, windowLength = %u\n", kmerLength, windowLength);

    string inputFile = args.input_arg;
    size_t firstindex = inputFile.find_last_of("/") + 1;
    size_t lastindex = inputFile.find_last_of(".");
    outputFile = inputFile.substr(firstindex, lastindex - firstindex);
    string buffer = "";
    char bufferArr[50];
    if (args.element_given)
    {
        bf_element_cnt = args.element_arg;
    }

    if (args.step_given)
    {
        step_size = args.step_arg;
        sprintf(bufferArr, "-k%u-w%u-s%zu-b%zu--step%u", kmerLength, windowLength, minimiser_size, bf_element_cnt, step_size);
    }
    else
    {
        step_size = windowLength;
        sprintf(bufferArr, "-k%u-w%u-s%zu-b%zu", kmerLength, windowLength, minimiser_size, bf_element_cnt);
    }
    buffer = buffer + bufferArr;
    if (args.toSingle_arg)
    {
        buffer = buffer + "-single";
    }

    if (args.folder_arg)
    {
        if (args.output_given)
        {
            outputFile = args.output_arg;
        }
        else
        {
            fprintf(stderr, "Please provide output file path\n");
            return 0;
        }

        string delimiter = "/*";
        string folder = inputFile.substr(0, inputFile.find(delimiter));
        string ext = inputFile.substr(inputFile.find(delimiter) + delimiter.size(), inputFile.size() - 1);
        fprintf(stderr, "Reading folder %s\n", folder.c_str());

        outputFile = outputFile + buffer;

        bloom_filter bf(parameters);

        if (multipleOut)
        {
            size_t i = 0;
            for (const auto &entry : fs::directory_iterator(folder))
            {
                if (entry.path().extension() == ext)
                {
                    outputFile = outputFile + "_" + to_string(i) + ".txt";
                    FILE *wf = fopen(outputFile.c_str(), "w");
                    doWork(wf, parameters, entry.path());
                    i++;
                }
            }
        }
        else
        {
            outputFile = outputFile + ".txt";
            FILE *wf = fopen(outputFile.c_str(), "w");
            for (const auto &entry : fs::directory_iterator(folder))
            {
                if (entry.path().extension() == ext)
                {
                    // cout << entry.path().stem().string() << '\n';
                    doWork(wf, parameters, entry.path());
                }
            }
            return 0;
        }
    }

    else
    {
        if (args.output_given)
        {
            outputFile = args.output_arg;
        }
        outputFile = outputFile + buffer + ".txt";
        FILE *wf = fopen(outputFile.c_str(), "w");
        bloom_filter bf(parameters);
        doWork(wf, parameters, inputFile);
    }
    return 0;
}
