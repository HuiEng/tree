#include "bloom_filter.hpp"
#include "minimiser.hpp"
#include <random>

#include "build_main_cmdline.hpp"
// #include "bf/all.hpp"

// using namespace bf;

using namespace std;

static build_main_cmdline args;   // Command line switches and arguments
static uint8_t kmerLength = 4;    // Kmer length
static uint32_t windowLength = 8; // Kmer length
size_t bf_element_cnt = 1000;
uint64_t seed = 0; //0x8F3F73B5CF1C9ADE; // The default seed from minimiser_hash

void writeInt(std::ostream &os, unsigned long long int i)
{
    os.write(reinterpret_cast<const char *>(&i), sizeof(i));
}

template <typename view>
void generateSig(view minimiser_view, bloom_parameters parameters, string filename, string outname, bool multiple)
{
    seqan3::sequence_file_input<dna4_traits> file_in{filename};

    bloom_filter bf(parameters);
    // basic_bloom_filter b(0.8, 100);

    if (multiple)
    {
        // Retrieve the sequences and ids.
        for (auto &[seq, id, qual] : file_in)
        {
            std::vector<std::string> textVector;
            string textString = id;
            std::istringstream iss(textString);
            std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), std::back_inserter(textVector));

            ofstream wf(outname + textVector[0] + ".bin", ios::out | ios::binary);
            writeInt(wf, bf.table_size());
            // for each seq, Instantiate Bloom Filter
            for (auto &&hash : seq | minimiser_view)
            {
                bf.insert(hash);
            }

            bf.print(wf);
            bf.clear();
            wf.close();
        }
    }
    else
    {
        ofstream wf(outname, ios::out | ios::binary);
        writeInt(wf, bf.table_size());

        // Retrieve the sequences and ids.
        for (auto &[seq, id, qual] : file_in)
        {
            // for each seq, Instantiate Bloom Filter
            for (auto &&hash : seq | minimiser_view)
            {
                bf.insert(hash);
            }

            bf.print(wf);
            bf.clear();
        }
        wf.close();
    }
}

void generateSig(bloom_parameters parameters, string filename, string outname, bool multiple)
{
    // to get minimisers with w=8,k=4
    // input param for the minimiser view is calculated by: window size - k-mer size + 1, here: 8 - 4 + 1 = 5)
    size_t temp = windowLength - kmerLength + 1;
    // auto minimiser_view = seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmerLength}}) | seqan3::views::minimiser(temp);

    // random ordering

    // Use XOR on all minimiser values
    auto minimiser_view = seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmerLength}}) | std::views::transform([seed](uint64_t i)
                                                                                                                        { return i ^ seed; }) |
                          seqan3::views::minimiser(temp);

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(1, numeric_limits<size_t>::max()/2); // distribution in range [1, 6]
    seed = dist(rng);

    auto minimiser_view_rand = seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmerLength}}) | std::views::transform([seed](uint64_t i)
                                                                                                                             { return i ^ seed; }) |
                               seqan3::views::minimiser(temp);

    seqan3::sequence_file_input<dna4_traits> file_in{filename};

    bloom_filter bf(parameters);
    // basic_bloom_filter b(0.8, 100);

    if (multiple)
    {
        // Retrieve the sequences and ids.
        for (auto &[seq, id, qual] : file_in)
        {
            std::vector<std::string> textVector;
            string textString = id;
            std::istringstream iss(textString);
            std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), std::back_inserter(textVector));

            ofstream wf(outname + textVector[0] + ".bin", ios::out | ios::binary);
            writeInt(wf, bf.table_size());
            // for each seq, Instantiate Bloom Filter
            for (auto &&hash : seq | minimiser_view)
            {
                bf.insert(hash);
            }

            bf.print(wf);
            bf.clear();
            wf.close();
        }
    }
    else
    {
        // auto &&hashes = seq | kmer_view | seqan3::views::minimiser(temp);

        // seed = (seeds[0]);
        // auto &&hashes_1 = seq | kmer_view
        // | std::views::transform([seed](uint64_t i) { return i ^ seed; })  | seqan3::views::minimiser(temp);

        // auto start = hashes.begin();
        // auto end = hashes.end();

        // while (start != end)
        // {
        //     // std::cout<< *start <<",";
        //     // bloom_filter bf(parameters);
        //     bf.insert(*start);
        //     ++start;
        // }

        // for (auto ptr = hashes.begin(); ptr < hashes.end(); ptr++){
        //     std::cout<< ptr <<"\n";
        // }

        ofstream wf(outname, ios::out | ios::binary);
        writeInt(wf, bf.table_size());
        std::cout<< "bf size: " << bf.table_size() * bits_per_char <<"\n";

        // Retrieve the sequences and ids.
        for (auto &[seq, id, qual] : file_in)
        {

            
            // for each seq, Instantiate Bloom Filter
            for (auto &&hash : seq | minimiser_view)
            {
                bf.insert(hash);
            }

            // for (auto &&hash : seq | minimiser_view_rand)
            // {
            //     bf.insert(hash);
            // }

            bf.print(wf);
            bf.clear();
        }
        wf.close();
    }
}

int build_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    std::ios::sync_with_stdio(false); // No sync with stdio -> faster

    string bfOut = "out.bin";

    if (args.bf_output_given)
        bfOut = args.bf_output_arg;
    if (args.kmer_given)
        kmerLength = args.kmer_arg;
    if (args.window_given)
        windowLength = args.window_arg;
    if (args.element_given)
        bf_element_cnt = args.element_arg;
    if (args.seed_given)
        seed = args.seed_arg;

    if (kmerLength > windowLength)
    {
        fprintf(stderr, "Error: kmer length must be smaller or equal to window length\n");
        return 1;
    }

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

    

    std::cout << "done setting bf param" << std::endl;

    if (kmerLength == windowLength)
    {
        fprintf(stderr, "Generating all k-mers...\n");
        if (args.canonical_arg)
        {
            auto minimiser_view = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{kmerLength}},
                                                                seqan3::window_size{windowLength},
                                                                seqan3::seed{0});
            generateSig(minimiser_view, parameters, args.input_arg, bfOut, args.multiple_arg);
        }
        else
        {
            auto minimiser_view = seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmerLength}});
            generateSig(minimiser_view, parameters, args.input_arg, bfOut, args.multiple_arg);
        }
    }
    else
    {
        fprintf(stderr, "Generating minimisers...\n");

        if (args.canonical_arg)
        {
            auto minimiser_view = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{kmerLength}},
                                                                seqan3::window_size{windowLength},
                                                                seqan3::seed{0});
            generateSig(minimiser_view, parameters, args.input_arg, bfOut, args.multiple_arg);
        }
        else
        {

            generateSig(parameters, args.input_arg, bfOut, args.multiple_arg);
        }
    }

    return 0;
}
