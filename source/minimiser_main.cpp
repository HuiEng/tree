#include "minimiser.hpp"
// #include "part_ktree.hpp"

#include "minimiser_main_cmdline.hpp"

using namespace std;

static minimiser_main_cmdline args; // Command line switches and arguments
static uint8_t kmerLength = 4;      // Kmer length
static uint32_t windowLength = 8;   // window length

template <typename view>
void getMinimisers(view minimiser_view, string filename, string outname)
{
    mer_table minimiser_table{};

    seqan3::sequence_file_input<dna4_traits> file_in{filename};
    ofstream outfile(outname);

    // Retrieve the sequences and ids.
    for (auto &[seq, id, qual] : file_in)
    {
        for (auto &&hash : seq | minimiser_view)
        {
            minimiser_table[hash] = minimiser_table[hash] + 1;
        }
    }
    printToFile(kmerLength, outfile, minimiser_table);
}

int minimiser_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    std::ios::sync_with_stdio(false); // No sync with stdio -> faster

    string outfile = "minimiser.csv";

    if (args.output_given)
        outfile = args.output_arg;
    if (args.kmer_given)
        kmerLength = args.kmer_arg;
    if (args.window_given)
        windowLength = args.window_arg;

    if (kmerLength > windowLength)
    {
        fprintf(stderr, "Error: kmer length must be smaller or equal to window length\n");
        return 1;
    }

    if (kmerLength == windowLength)
    {
        fprintf(stderr, "Generating all k-mers...\n");
        if (args.canonical_arg)
        {
            auto minimiser_view = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{kmerLength}},
                                                                seqan3::window_size{windowLength},
                                                                seqan3::seed{0});
            getMinimisers(minimiser_view, args.input_arg, outfile);
        }
        else
        {
            auto minimiser_view = seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmerLength}});
            getMinimisers(minimiser_view, args.input_arg, outfile);
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
            getMinimisers(minimiser_view, args.input_arg, outfile);
        }
        else
        {
            // to get minimisers with w=8,k=4
            // input param for the minimiser view is calculated by: window size - k-mer size + 1, here: 8 - 4 + 1 = 5)
            size_t temp = windowLength - kmerLength + 1;
            auto minimiser_view = seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmerLength}}) | seqan3::views::minimiser(temp);

            getMinimisers(minimiser_view, args.input_arg, outfile);
        }
    }
    return 0;
}
