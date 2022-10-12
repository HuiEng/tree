// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Mitra Darvish <mitra.darvish AT fu-berlin.de>
 * \brief Provides seqan3::views::partition_multi_hash.
 */

#pragma once

#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/partition_multi.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

namespace seqan3
{
    //!\brief strong_type for the minimiser_size.
    //!\ingroup search_views
    struct minimiser_size : seqan3::detail::strong_type<size_t, minimiser_size>
    {
        using seqan3::detail::strong_type<size_t, minimiser_size>::strong_type;
    };

    //!\brief strong_type for the step_size.
    //!\ingroup search_views
    struct step_size : seqan3::detail::strong_type<size_t, step_size>
    {
        using seqan3::detail::strong_type<size_t, step_size>::strong_type;
    };
} // namespace seqan3

namespace seqan3::detail
{
    //!\brief seqan3::views::partition_multi_hash's range adaptor object type (non-closure).
    //!\ingroup search_views
    struct partition_multi_hash_fn
    {
        /*!\brief Store the shape and the window size and return a range adaptor closure object.
         * \param[in] shape       The seqan3::shape to use for hashing.
         * \param[in] window_size The windows size to use.
         * \param[in] minimiser_size The number of minimiser per window to use.
         * \param[in] step_size   The sliding step
         * \throws std::invalid_argument if the size of the shape is greater than the `window_size`.
         * \returns               A range of converted elements.
         */
        constexpr auto operator()(shape const &shape, window_size const window_size, minimiser_size const minimiser_size, step_size const step_size) const
        {
            return seqan3::detail::adaptor_from_functor{*this, shape, window_size, minimiser_size, step_size};
        }

        /*!\brief Store the shape, the window size and the seed and return a range adaptor closure object.
         * \param[in] shape       The seqan3::shape to use for hashing.
         * \param[in] window_size The size of the window.
         * \param[in] minimiser_size The number of minimiser per window to use.
         * \param[in] step_size   The sliding step
         * \param[in] seed        The seed to use.
         * \throws std::invalid_argument if the size of the shape is greater than the `window_size`.
         * \returns               A range of converted elements.
         */
        constexpr auto operator()(shape const &shape, window_size const window_size, minimiser_size const minimiser_size, step_size const step_size, seed const seed) const
        {
            return seqan3::detail::adaptor_from_functor{*this, shape, window_size, minimiser_size, step_size, seed};
        }

        /*!\brief Call the view's constructor with the underlying view, a seqan3::shape and a window size as argument.
         * \param[in] urange      The input range to process. Must model std::ranges::viewable_range and the reference type
         *                        of the range must model seqan3::semialphabet.
         * \param[in] shape       The seqan3::shape to use for hashing.
         * \param[in] window_size The size of the window.
         * \param[in] minimiser_size The number of minimiser per window to use.
         * \param[in] step_size   The sliding step
         * \param[in] seed        The seed to use.
         * \throws std::invalid_argument if the size of the shape is greater than the `window_size`.
         * \returns               A range of converted elements.
         */
        template <std::ranges::range urng_t>
        constexpr auto operator()(urng_t &&urange,
                                  shape const &shape,
                                  window_size const window_size,
                                  minimiser_size const minimiser_size, 
                                  step_size const step_size,
                                  seed const seed = seqan3::seed{0x8F3F73B5CF1C9ADE}) const
        {
            static_assert(std::ranges::viewable_range<urng_t>,
                          "The range parameter to views::partition_multi_hash cannot be a temporary of a non-view range.");
            static_assert(std::ranges::forward_range<urng_t>,
                          "The range parameter to views::partition_multi_hash must model std::ranges::forward_range.");
            static_assert(semialphabet<std::ranges::range_reference_t<urng_t>>,
                          "The range parameter to views::partition_multi_hash must be over elements of seqan3::semialphabet.");

            if (shape.size() > window_size.get())
                throw std::invalid_argument{"The size of the shape cannot be greater than the window size."};

            auto forward_strand = std::forward<urng_t>(urange) | seqan3::views::kmer_hash(shape) | std::views::transform([seed](uint64_t i)
                                                                                                                         { return i ^ seed.get(); });

            auto reverse_strand = std::forward<urng_t>(urange) | seqan3::views::complement | std::views::reverse | seqan3::views::kmer_hash(shape) | std::views::transform([seed](uint64_t i)
                                                                                                                                                                           { return i ^ seed.get(); }) |
                                  std::views::reverse;

            return seqan3::detail::partition_multi_view(forward_strand, reverse_strand, window_size.get() - shape.size() + 1, shape.size(), minimiser_size.get(), step_size.get());
        }
    };

} // namespace seqan3::detail

namespace seqan3::views
{

    /*!\name Alphabet related views
     * \{
     */

    /*!\brief                    Computes minimisers for a range with a given shape, window size and seed.
     * \tparam urng_t            The type of the range being processed. See below for requirements. [template parameter is
     *                           omitted in pipe notation]
     * \param[in] urange         The range being processed. [parameter is omitted in pipe notation]
     * \param[in] shape          The seqan3::shape that determines how to compute the hash value.
     * \param[in] window_size    The window size to use.
     * \param[in] minimiser_size    The number of minimiser per window to use.
     * \param[in] step_size   The sliding step
     * \param[in] seed           The seed used to skew the hash values. Default: 0x8F3F73B5CF1C9ADE.
     * \returns                  A range of `size_t` where each value is the minimiser of the resp. window.
     *                           See below for the properties of the returned range.
     * \ingroup views
     *
     * \details
     *
     * A sequence can be presented by a small number of k-mers (minimisers). For a given shape and window size all k-mers
     * are determined in the forward strand and the backward strand and only the lexicographically smallest k-mer is
     * returned for one window. This process is repeated over every possible partition_multi of a sequence.
     * For example, in the sequence "TAAAGTGCTAAA" for an ungapped shape of length 3 and a window size of 5 the first
     * window will be "TAAAG", the second window will be "TGCTA".
     *
     * ### Non-lexicographical Minimisers by skewing the hash value with a seed
     * The user can change the seed to any other value he or she thinks is useful. A seed of 0 is returning the
     * lexicographical order.
     *
     * \sa seqan3::views::partition_multi_view
     *
     * \attention
     * Be aware of the requirements of the seqan3::views::kmer_hash view.
     *
     * \experimentalapi
     *
     * ### View properties
     *
     * | Concepts and traits              | `urng_t` (underlying range type)   | `rrng_t` (returned range type)   |
     * |----------------------------------|:----------------------------------:|:--------------------------------:|
     * | std::ranges::input_range         | *required*                         | *preserved*                      |
     * | std::ranges::forward_range       | *required*                         | *preserved*                      |
     * | std::ranges::bidirectional_range |                                    | *lost*                           |
     * | std::ranges::random_access_range |                                    | *lost*                           |
     * | std::ranges::contiguous_range    |                                    | *lost*                           |
     * |                                  |                                    |                                  |
     * | std::ranges::viewable_range      | *required*                         | *guaranteed*                     |
     * | std::ranges::view                |                                    | *guaranteed*                     |
     * | std::ranges::sized_range         |                                    | *lost*                           |
     * | std::ranges::common_range        |                                    | *lost*                           |
     * | std::ranges::output_range        |                                    | *lost*                           |
     * | seqan3::const_iterable_range     |                                    | *preserved*                      |
     * |                                  |                                    |                                  |
     * | std::ranges::range_reference_t   | seqan3::semialphabet               | std::size_t                      |
     *
     * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
     *
     * ### Example
     *
     * \include test/snippet/search/views/partition_multi_hash.cpp
     *
     * \hideinitializer
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    inline constexpr auto partition_multi_hash = detail::partition_multi_hash_fn{};

    //!\}

} // namespace seqan3::views
