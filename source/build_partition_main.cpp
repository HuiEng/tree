#include <seqan3/search/views/partition_multi_hash.hpp>
#include <seqan3/search/views/partition_multi.hpp>
#include "minimiser.hpp"
#include "bloom_filter.hpp"
// #include "part_ktree.hpp"

#include "build_partition_main_cmdline.hpp"

using namespace std;

static build_partition_main_cmdline args; // Command line switches and arguments
static uint8_t kmerLength = 4;            // Kmer length
static uint32_t windowLength = 8;         // window length
static uint32_t step_size = windowLength; // window length
size_t minimiser_size = 3;
size_t bf_element_cnt = 1000;
bool debug = false;

void writeInt(std::ostream &os, unsigned long long int i)
{
    os.write(reinterpret_cast<const char *>(&i), sizeof(i));
}

template <typename view>
void getPartitionMinimisers(view minimiser_view, bloom_parameters parameters, string filename, string outname)
{
    seqan3::sequence_file_input<dna4_traits> file_in{filename};
    // ofstream outfile(outname);
    size_t max = 8;
    ofstream wf(outname, ios::out | ios::binary);
    size_t end = -1;

    bloom_filter bf(parameters);
    writeInt(wf, bf.table_size());

    if (debug)
    {
        size_t i = 0;

        // Retrieve the sequences and ids.
        for (auto &[seq, id, qual] : file_in)
        {
            // fprintf(stdout, ">\n");
            cout << ">0|" << i << "|0\n";
            for (auto &&hashes : seq | minimiser_view)
            {
                for (size_t hash : hashes)
                {
                    bf.insert(hash);
                    cout << hashToMer_str(kmerLength, hash) << ";";
                }
                cout << ",";
                bf.print(wf);
                // bf.printBFIdx(stderr);

                // ///// debug
                // bf.print();

                bf.clear();
            }
            cout << "\n";

            // fprintf(stderr,">s\n");

            // end of seq flag, print empty bf
            bf.print(wf);
            i++;
        }
    }
    else
    {

        // Retrieve the sequences and ids.
        for (auto &[seq, id, qual] : file_in)
        {
            // fprintf(stdout, ">\n");
            for (auto &&hashes : seq | minimiser_view)
            {
                for (size_t hash : hashes)
                {
                    bf.insert(hash);
                }
                bf.print(wf);
                // bf.printBFIdx(stderr);

                // ///// debug
                // bf.print();

                bf.clear();
            }

            // fprintf(stderr,">s\n");

            // end of seq flag, print empty bf
            bf.print(wf);
        }
    }
    wf.close();
}

template <typename view>
void getPartitionMinimisersSet(view minimiser_view, bloom_parameters parameters, string filename, string outname)
{
    seqan3::sequence_file_input<dna4_traits> file_in{filename};
    // ofstream outfile(outname);
    size_t max = 8;
    ofstream wf(outname, ios::out | ios::binary);
    size_t end = -1;

    bloom_filter bf(parameters);
    writeInt(wf, bf.table_size());
    // Retrieve the sequences and ids.
    for (auto &[seq, id, qual] : file_in)
    {
        // fprintf(stdout, ">\n");
        for (auto &&hashes : seq | minimiser_view)
        {
            for (size_t hash : hashes)
            {
                bf.insert(hash);
            }
        }

        // fprintf(stderr,">s\n");

        // end of seq flag, print empty bf
        bf.print(wf);
        bf.clear();
    }

    wf.close();
}

int build_partition_main(int argc, char *argv[])
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

    if (args.size_given)
    {
        minimiser_size = args.size_arg;
    }
    fprintf(stderr, "Partition - Generating %zu minimisers per window...\n", minimiser_size);
    fprintf(stderr, "kmerLength = %zu, windowLength = %zu\n", kmerLength, windowLength);

    string inputFile = args.input_arg;
    size_t firstindex = inputFile.find_last_of("/") + 1;
    size_t lastindex = inputFile.find_last_of(".");
    string rawname = inputFile.substr(firstindex, lastindex - firstindex);
    char buffer[50];

    if (args.step_given)
    {
        step_size = args.step_arg;
        sprintf(buffer, "-k%zu-w%zu-s%zu--step%zu.bin", kmerLength, windowLength, minimiser_size, step_size);
    }
    else
    {
        step_size = windowLength;
        sprintf(buffer, "-k%zu-w%zu-s%zu.bin", kmerLength, windowLength, minimiser_size);
    }

    string outfile = rawname + buffer;

    if (args.output_given)
        outfile = args.output_arg;

    if (args.canonical_arg)
    {
        auto partition_view = seqan3::views::partition_multi_hash(seqan3::shape{seqan3::ungapped{kmerLength}},
                                                                  seqan3::window_size{windowLength},
                                                                  seqan3::minimiser_size{minimiser_size},
                                                                  seqan3::step_size{step_size},
                                                                  seqan3::seed{0});
        getPartitionMinimisers(partition_view, parameters, args.input_arg, outfile);
    }
    else
    {
        // to get minimisers with w=8,k=4
        // input param for the minimiser view is calculated by: window size - k-mer size + 1, here: 8 - 4 + 1 = 5)
        size_t temp = windowLength - kmerLength + 1;
        auto partition_view = seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmerLength}}) | seqan3::views::partition_multi(temp, kmerLength, minimiser_size, step_size);

        if (args.set_arg){
            getPartitionMinimisersSet(partition_view, parameters, args.input_arg, outfile);
        }else{
            getPartitionMinimisers(partition_view, parameters, args.input_arg, outfile);
        }
    }

    return 0;
}
