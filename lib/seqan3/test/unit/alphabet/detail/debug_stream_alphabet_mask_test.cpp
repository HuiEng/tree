// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/mask/mask.hpp>

TEST(debug_stream_test, mask)
{
    {
        std::ostringstream o;
        seqan3::debug_stream_type my_stream{o};
        my_stream << seqan3::mask::masked;
        EXPECT_EQ(o.str(), "MASKED");
    }

    {
        std::ostringstream o;
        seqan3::debug_stream_type my_stream{o};
        my_stream << seqan3::mask::unmasked;
        EXPECT_EQ(o.str(), "UNMASKED");
    }
}
