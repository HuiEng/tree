#include "bloom_filter.hpp"
#include "minimiser.hpp"
#include <random>
#include <filesystem>
#include "build_main_cmdline.hpp"
// #include "bf/all.hpp"

// using namespace bf;

using namespace std;
namespace fs = std::filesystem;

static build_main_cmdline args;   // Command line switches and arguments
static uint8_t kmerLength = 4;    // Kmer length
static uint32_t windowLength = 8; // Kmer length
size_t bf_element_cnt = 1000;
uint64_t seed = 0; // 0x8F3F73B5CF1C9ADE; // The default seed from minimiser_hash
bool debug = false;

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
    auto minimiser_view = seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmerLength}}) | seqan3::views::minimiser(temp);

    // random ordering

    // // Use XOR on all minimiser values
    // auto minimiser_view = seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmerLength}}) | std::views::transform([seed](uint64_t i)
    //                                                                                                                     { return i ^ seed; }) |
    //                       seqan3::views::minimiser(temp);

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(1, numeric_limits<size_t>::max() / 2); // distribution in range [1, 6]
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
    else if (debug)
    {

        ofstream wf(outname, ios::out | ios::binary);
        writeInt(wf, bf.table_size());
        std::cout << "bf size: " << bf.table_size() * bits_per_char << "\n";
        size_t i = 0;

        // Retrieve the sequences and ids.
        for (auto &[seq, id, qual] : file_in)
        {
            cout << ">0|" << i << "|0\n";
            // for each seq, Instantiate Bloom Filter
            for (auto &&hash : seq | minimiser_view)
            {
                bf.insert(hash);
                cout << hashToMer_str(kmerLength, hash) << ",";
            }
            cout << "\n";
            i++;

            bf.print(wf);
            bf.clear();
        }
        wf.close();
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
        std::cout << "bf size: " << bf.table_size() * bits_per_char << "\n";

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

template <typename view>
void generateSigFolder(view minimiser_view, bloom_parameters parameters, string filename, ofstream &wf)
{
    seqan3::sequence_file_input<dna4_traits> file_in{filename};

    bloom_filter bf(parameters);
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
}

void generateSigFolder(bloom_parameters parameters, string filename, ofstream &wf)
{
    // to get minimisers with w=8,k=4
    // input param for the minimiser view is calculated by: window size - k-mer size + 1, here: 8 - 4 + 1 = 5)
    size_t temp = windowLength - kmerLength + 1;
    auto minimiser_view = seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmerLength}}) | seqan3::views::minimiser(temp);
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(1, numeric_limits<size_t>::max() / 2); // distribution in range [1, 6]
    seed = dist(rng);

    auto minimiser_view_rand = seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmerLength}}) | std::views::transform([seed](uint64_t i)
                                                                                                                             { return i ^ seed; }) |
                               seqan3::views::minimiser(temp);

    seqan3::sequence_file_input<dna4_traits> file_in{filename};

    bloom_filter bf(parameters);

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
}

int build_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    std::ios::sync_with_stdio(false); // No sync with stdio -> faster

    if (args.kmer_given)
        kmerLength = args.kmer_arg;
    if (args.window_given)
        windowLength = args.window_arg;
    if (args.seed_given)
        seed = args.seed_arg;

    if (kmerLength > windowLength)
    {
        fprintf(stderr, "Error: kmer length must be smaller or equal to window length\n");
        return 1;
    }
    debug = args.debug;

    string inputFile = args.input_arg;
    size_t firstindex = inputFile.find_last_of("/") + 1;
    size_t lastindex = inputFile.find_last_of(".");
    string rawname = inputFile.substr(firstindex, lastindex - firstindex);
    char buffer[50];
    if (args.element_given)
    {
        bf_element_cnt = args.element_arg;
        if (kmerLength == windowLength)
        {
            sprintf(buffer, "-k%zu-b%zu.bin", kmerLength, bf_element_cnt);
        }
        else
        {
            sprintf(buffer, "-k%zu-w%zu-b%zu.bin", kmerLength, windowLength, bf_element_cnt);
        }
    }
    else
    {
        if (kmerLength == windowLength)
        {
            sprintf(buffer, "-k%zu.bin", kmerLength);
        }
        else
        {
            sprintf(buffer, "-k%zu-w%zu.bin", kmerLength, windowLength);
        }
    }
    string bfOut = rawname + buffer;

    if (args.bf_output_given)
        bfOut = args.bf_output_arg;

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

    if (args.folder_arg)
    {
        string inputFile = args.input_arg;
        string delimiter = "/*";
        string folder = inputFile.substr(0, inputFile.find(delimiter));
        string ext = inputFile.substr(inputFile.find(delimiter) + delimiter.size(), inputFile.size() - 1);
        fprintf(stderr, "Reading folder %s\n", folder.c_str());

        string outfile = folder + buffer;

        ofstream wf(outfile, ios::out | ios::binary);

        bloom_filter bf(parameters);
        writeInt(wf, bf.table_size());
        std::cout << "bf size: " << bf.table_size() * bits_per_char << "\n";

        if (kmerLength == windowLength)
        {
            fprintf(stderr, "Generating all k-mers...\n");
            if (args.canonical_arg)
            {
                auto minimiser_view = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{kmerLength}},
                                                                    seqan3::window_size{windowLength},
                                                                    seqan3::seed{0});

                for (const auto &entry : fs::directory_iterator(folder))
                {
                    if (entry.path().extension() == ext)
                    {
                        cout << entry.path().stem().string() << '\n';
                        generateSigFolder(minimiser_view, parameters, entry.path(), wf);
                    }
                }
            }
            else
            {
                auto minimiser_view = seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmerLength}});
                for (const auto &entry : fs::directory_iterator(folder))
                {
                    if (entry.path().extension() == ext)
                    {
                        cout << entry.path().stem().string() << '\n';
                        generateSigFolder(minimiser_view, parameters, entry.path(), wf);
                    }
                }
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
                for (const auto &entry : fs::directory_iterator(folder))
                {
                    if (entry.path().extension() == ext)
                    {
                        cout << entry.path().stem().string() << '\n';
                        generateSigFolder(minimiser_view, parameters, entry.path(), wf);
                    }
                }
            }
            else
            {
                for (const auto &entry : fs::directory_iterator(folder))
                {
                    if (entry.path().extension() == ext)
                    {
                        cout << entry.path().stem().string() << '\n';
                        generateSigFolder(parameters, entry.path(), wf);
                    }
                }
            }
        }

        wf.close();
        return 0;
    }

    if (kmerLength == windowLength)
    {
        fprintf(stderr, "Generating all k-mers...\n");
        if (args.canonical_arg)
        {
            auto minimiser_view = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{kmerLength}},
                                                                seqan3::window_size{windowLength},
                                                                seqan3::seed{0});
            generateSig(minimiser_view, parameters, inputFile, bfOut, args.multiple_arg);
        }
        else
        {
            auto minimiser_view = seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmerLength}});
            generateSig(minimiser_view, parameters, inputFile, bfOut, args.multiple_arg);
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
            generateSig(minimiser_view, parameters, inputFile, bfOut, args.multiple_arg);
        }
        else
        {

            generateSig(parameters, inputFile, bfOut, args.multiple_arg);
        }
    }

    return 0;
}
