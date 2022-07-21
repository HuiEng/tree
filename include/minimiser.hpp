/*
 *********************************************************************
 *                                                                   *
 *                           Minimisers                              *
 *                                                                   *
 *********************************************************************
*/

#ifndef INCLUDE_MINIMISER_HPP
#define INCLUDE_MINIMISER_HPP

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/search/views/minimiser.hpp>
// #include <seqan3/search/views/partition_hash.hpp>
// #include <seqan3/search/views/partition.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/utility/views/pairwise_combine.hpp>

#include <math.h>
#include <unordered_map>
#include <fstream>
#include <robin_hood.h>


using namespace seqan3::literals;
using namespace std;
typedef robin_hood::unordered_map<size_t, uint8_t> mer_table;


struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4; // instead of dna5
};
 

string hashToMer_str(int k, size_t hash)
{
    vector<char> mer;
    for (int i = 0; i < k; i++)
    {
        double temp = pow(4, (k - i - 1));
        int rank = (int)floor(hash / temp);
        hash -= temp * rank;

        seqan3::dna4 dna;
        dna.assign_rank(rank);
        mer.push_back(dna.to_char());
    }
    return string(mer.begin(), mer.end());
}

void printToFile(int k, ofstream &outfile, mer_table minimiser_table)
{
    for (auto &&hash : minimiser_table)
    {
        auto mer = hashToMer_str(k, hash.first);
        outfile<< ">\n" << mer << "\n";
        // outfile << mer << ",";
        // outfile << (int)hash.second << "\n";
        // seqan3::debug_stream << hash.second << '\n';
    }
}
void edit_distance()
{
    std::vector vec{std::pair{"AGTGCTACG"_dna4, "ACGTGCGACTAG"_dna4},
                    std::pair{"AGTAGACTACG"_dna4, "ACGTACGACACG"_dna4},
                    std::pair{"AGTTACGAC"_dna4, "AGTAGCGATCG"_dna4}};

    // Compute the alignment of a single pair.
    auto edit_config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme;

    // for (auto const &res : seqan3::align_pairwise(vec[0], edit_config))
    //     seqan3::debug_stream << "The score: " << res.score() << "\n";

    int i = 0;
    int j = i + 1;
    // Compute the alignment over a range of pairs.
    for (auto const &res : seqan3::align_pairwise(vec, edit_config))
    {
        if (j == 3)
        {
            i++;
            j = i + 1;
        }
        seqan3::debug_stream << i << ", " << j << ", " << res.score() << "\n";

        j++;
    }
}

#endif