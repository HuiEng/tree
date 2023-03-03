// #include "bf_ktree.hpp"
#include <fstream>
#include <random>
#include "stats_main_cmdline.hpp"
#include "stats.hpp"

#include <regex>
#include <vector>
#include <numeric>
#include <iostream>

using namespace std;

static stats_main_cmdline args; // Command line switches and arguments

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
