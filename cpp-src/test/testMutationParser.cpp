/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2017 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#include "gtest/gtest.h"
#include "MutationParser.cpp"

namespace BF = boost::filesystem;
using namespace mutation;
using namespace mutation_parser;
using namespace mutation_parser::detail;

std::string FILEPATH = __FILE__;

BF::path getTestFileDir() {
    BF::path PARENTPATH = BF::path(FILEPATH).parent_path();
    BF::path FILEDIR = PARENTPATH / "files";
    return FILEDIR;
}

/*std::string getTestFilePath() {
    return (getTestFileDir() / "2_bowtie_converted.bam").string();
}*/



/*TEST(MutationParsing, LoadBam) {
    std::string file_in = getTestFilePath();
    //std::cout << "About to test loading BAM file " << file_in << std::endl;
    //EXPECT_NO_THROW(mutation_parser::detail::loadBam(file_in));
    try {
        mutation_parser::detail::loadBam(file_in);
        //FAIL();
    }
    catch (const std::runtime_error &err) {
        //std::cout << "runtime_error: " << err.what() << std::endl; // doesn't output
        //ASSERT_STREQ
    }
}*/

// MD tag parsing
TEST(MdOp, ToString) {
    MdOp m{MD_MATCH, 1, ""};
    EXPECT_NO_THROW(toString(m));
}

TEST(MdOpVector, ToString) {
    std::vector<MdOp> m{{MD_MATCH,    1, ""},
                        {MD_DELETION, 3, "ATG"}};
    EXPECT_NO_THROW(toString(m));
}

TEST(ParseMDtag, Match) {
    std::string tag = "137";
    std::vector<MdOp> output = parseMDtag(tag);
    std::vector<MdOp> expected{{MD_MATCH, 137, ""}};
    // Note: if we used Boost's testing lib, could use BOOST_CHECK_EQUAL_COLLECTIONS in these cases
    // instead of converting to strings or testing each element
    EXPECT_EQ(toString(expected), toString(output));
}

TEST(ParseMDtag, Deletion) {
    std::string tag = "^ATGCATGC";
    std::vector<MdOp> output = parseMDtag(tag);
    std::vector<MdOp> expected{{MD_DELETION, 8, "ATGCATGC"}};
    EXPECT_EQ(toString(expected), toString(output));
}

TEST(ParseMDtag, Mismatch) {
    std::string tag = "C";
    std::vector<MdOp> output = parseMDtag(tag);
    std::vector<MdOp> expected{{MD_MISMATCH, 1, "C"}};
    EXPECT_EQ(toString(expected), toString(output));
}

TEST(ParseMDtag, Complex) {
    std::string tag = "85G16G8^A0T2A0A20";
    std::vector<MdOp> output = parseMDtag(tag);
    std::vector<MdOp> expected{{MD_MATCH,    85, ""},
                               {MD_MISMATCH, 1,  "G"},
                               {MD_MATCH,    16, ""},
                               {MD_MISMATCH, 1,  "G"},
                               {MD_MATCH,    8,  ""},
                               {MD_DELETION, 1,  "A"},
                               {MD_MISMATCH, 1,  "T"},
                               {MD_MATCH,    2,  ""},
                               {MD_MISMATCH, 1,  "A"},
                               {MD_MISMATCH, 1,  "A"},
                               {MD_MATCH,    20, ""}};
    EXPECT_EQ(toString(expected), toString(output));
}

// CIGAR + MD parsing to locate mutations and reconstruct local alignment
// target sequence and aligned read over region
TEST(LocateMutations, OnlyMatch) {
    /*
     *
     *      inserts:
     *  insert locs:
     * aligned+gaps: ATGCATGCATGCATGC
     *       target: ATGCATGCATGCATGC
     */

    int left_target_pos = 0;
    std::string query_bases = "ATGCATGCATGCATGC";
    std::string query_qual =  "ABCDEFGHIJKLMNOP";
    std::vector<BAM::CigarOp> cigar_data{{'M', 16}};
    std::vector<MdOp> md = {{MD_MATCH, 16, ""}};


    std::vector<Mutation> exp_mutations = {};
    std::string exp_local_target_seq = "ATGCATGCATGCATGC";
    std::string exp_local_target_qual = "ABCDEFGHIJKLMNOP";
    std::string exp_aligned_query_seq = "ATGCATGCATGCATGC";
    std::string exp_aligned_query_qual = "ABCDEFGHIJKLMNOP";


    std::vector<Mutation> mutations;
    std::string local_target_seq;
    std::string local_target_qual;
    std::string aligned_query_seq;
    std::string aligned_query_qual;
    boost::tie(mutations,
               local_target_seq,
               local_target_qual,
               aligned_query_seq,
               aligned_query_qual) = locateMutations(left_target_pos,
                                                     query_bases,
                                                     query_qual,
                                                     cigar_data,
                                                     md,
                                                     true);

    EXPECT_EQ(toString(exp_mutations), toString(mutations));
    EXPECT_EQ(exp_local_target_seq, local_target_seq);
    EXPECT_EQ(exp_local_target_qual, local_target_qual);
    EXPECT_EQ(exp_aligned_query_seq, aligned_query_seq);
    EXPECT_EQ(exp_aligned_query_qual, aligned_query_qual);
}


TEST(LocateMutations, CigarMatchWithMdMismatch) {
/*
 *
 *      inserts:
 *  insert locs:
 * aligned+gaps: ATGCATGCGTGCATGC
 *       target: ATGCATGCATGCATGC
 */

    int left_target_pos = 0;
    std::string query_bases = "ATGCATGCGTGCATGC";
    std::string query_qual =  "ABCDEFGHIJKLMNOP";
    std::vector<BAM::CigarOp> cigar_data{{'M', 16}};
    std::vector<MdOp> md = {{MD_MATCH,    8, ""},
                            {MD_MISMATCH, 1, "A"},
                            {MD_MATCH,    7, ""}};


    std::vector<Mutation> exp_mutations = {{7, 9, "G", "I"}};
    std::string exp_local_target_seq =  "ATGCATGCATGCATGC";
    std::string exp_local_target_qual = "ABCDEFGHIJKLMNOP";
    std::string exp_aligned_query_seq = "ATGCATGCGTGCATGC";
    std::string exp_aligned_query_qual = "ABCDEFGHIJKLMNOP";


    std::vector<Mutation> mutations;
    std::string local_target_seq;
    std::string local_target_qual;
    std::string aligned_query_seq;
    std::string aligned_query_qual;
    boost::tie(mutations,
               local_target_seq,
               local_target_qual,
               aligned_query_seq,
               aligned_query_qual) = locateMutations(left_target_pos,
                                                     query_bases,
                                                     query_qual,
                                                     cigar_data,
                                                     md,
                                                     true);

    EXPECT_EQ(toString(exp_mutations), toString(mutations));
    EXPECT_EQ(exp_local_target_seq, local_target_seq);
    EXPECT_EQ(exp_local_target_qual, local_target_qual);
    EXPECT_EQ(exp_aligned_query_seq, aligned_query_seq);
    EXPECT_EQ(exp_aligned_query_qual, aligned_query_qual);
}


TEST(LocateMutations, InsertAfterGap) {
/*
 *        quals:           123
 *      inserts:           AAA
 *  insert locs:          /
 *        quals: ab cdefghijklmno
 * aligned+gaps: AT-CATGCATGCATGC
 *       target: ATGCATGCATGCATGC
 */

    int left_target_pos = 0;
    std::string query_bases = "ATCATGCAAAATGCATGC";
    std::string query_qual =  "abcdefgh123ijklmno";
    std::vector<BAM::CigarOp> cigar_data{{'M', 2},
                                         {'D', 1},
                                         {'M', 6},
                                         {'I', 3},
                                         {'M', 7}};
    std::vector<MdOp> md = {{MD_MATCH,    2,  ""},
                            {MD_DELETION, 1,  "G"},
                            {MD_MATCH,    13, ""}};


    std::vector<Mutation> exp_mutations = {{1, 3, "", ""},
                                           {8, 9, "AAA", "123"}};
    std::string exp_local_target_seq =  "ATGCATGCATGCATGC";
    std::string exp_local_target_qual = "ab!cdefghijklmno";
    std::string exp_aligned_query_seq  = "AT-CATGCATGCATGC";
    std::string exp_aligned_query_qual = "ab!cdefghijklmno";

    std::vector<Mutation> mutations;
    std::string local_target_seq;
    std::string local_target_qual;
    std::string aligned_query_seq;
    std::string aligned_query_qual;
    boost::tie(mutations,
               local_target_seq,
               local_target_qual,
               aligned_query_seq,
               aligned_query_qual) = locateMutations(left_target_pos,
                                                     query_bases,
                                                     query_qual,
                                                     cigar_data,
                                                     md,
                                                     true);

    EXPECT_EQ(toString(exp_mutations), toString(mutations));
    EXPECT_EQ(exp_local_target_seq, local_target_seq);
    EXPECT_EQ(exp_local_target_qual, local_target_qual);
    EXPECT_EQ(exp_aligned_query_seq, aligned_query_seq);
    EXPECT_EQ(exp_aligned_query_qual, aligned_query_qual);

}


TEST(LocateMutations, GapAfterInsert) {
/*
 *        quals:     123
 *      inserts:     AAA
 *  insert locs:    /
 *        quals: abcdefghijk lmno
 * aligned+gaps: ATGCATGCATG-ATGC
 *       target: ATGCATGCATGCATGC
 */

    int left_target_pos = 0;
    std::string query_bases = "ATGAAACATGCATGATGC";
    std::string query_qual =  "abc123defghijklmno";
    std::vector<BAM::CigarOp> cigar_data{{'M', 3},
                                         {'I', 3},
                                         {'M', 8},
                                         {'D', 1},
                                         {'M', 4}};
    std::vector<MdOp> md = {{MD_MATCH,    11, ""},
                            {MD_DELETION, 1,  "C"},
                            {MD_MATCH,    4,  ""}};


    std::vector<Mutation> exp_mutations = {{2,  3,  "AAA", "123"},
                                           {10, 12, "", ""}};
    std::string exp_local_target_seq =   "ATGCATGCATGCATGC";
    std::string exp_local_target_qual =  "abcdefghijk!lmno";
    std::string exp_aligned_query_seq =  "ATGCATGCATG-ATGC";
    std::string exp_aligned_query_qual = "abcdefghijk!lmno";

    std::vector<Mutation> mutations;
    std::string local_target_seq;
    std::string local_target_qual;
    std::string aligned_query_seq;
    std::string aligned_query_qual;
    boost::tie(mutations,
               local_target_seq,
               local_target_qual,
               aligned_query_seq,
               aligned_query_qual) = locateMutations(left_target_pos,
                                                     query_bases,
                                                     query_qual,
                                                     cigar_data,
                                                     md,
                                                     true);

    EXPECT_EQ(toString(exp_mutations), toString(mutations));
    EXPECT_EQ(exp_local_target_seq, local_target_seq);
    EXPECT_EQ(exp_local_target_qual, local_target_qual);
    EXPECT_EQ(exp_aligned_query_seq, aligned_query_seq);
    EXPECT_EQ(exp_aligned_query_qual, aligned_query_qual);
}


TEST(LocateMutations, SoftClippedWithMismatchMissingFirstNuc) {
/*
 *
 *      inserts:
 *  insert locs:
 *        quals:  abcdefghijklmno
 * aligned+gaps:  TGCATGCGTGCATGC
 *       target: ATGCATGCATGCATGC
 */

	int left_target_pos = 1;
	std::string query_bases = "GGGGGTGCATGCGTGCATGCGGGGG";
	std::string query_qual =  "HHHHHabcdefghijklmnoHHHHH";
	std::vector<BAM::CigarOp> cigar_data{{'S', 5},
	                                     {'M', 15},
	                                     {'S', 5}};
	std::vector<MdOp> md = {{MD_MATCH,    7, ""},
	                        {MD_MISMATCH, 1, "A"},
	                        {MD_MATCH,    7, ""}};


	std::vector<Mutation> exp_mutations = {{7, 9, "G", "h"}};
	std::string exp_local_target_seq =   "TGCATGCATGCATGC";
	std::string exp_local_target_qual =  "abcdefghijklmno";
	std::string exp_aligned_query_seq =  "TGCATGCGTGCATGC";
	std::string exp_aligned_query_qual = "abcdefghijklmno";


    std::vector<Mutation> mutations;
    std::string local_target_seq;
    std::string local_target_qual;
    std::string aligned_query_seq;
    std::string aligned_query_qual;
    boost::tie(mutations,
               local_target_seq,
               local_target_qual,
               aligned_query_seq,
               aligned_query_qual) = locateMutations(left_target_pos,
                                                     query_bases,
                                                     query_qual,
                                                     cigar_data,
                                                     md,
                                                     true);

    EXPECT_EQ(toString(exp_mutations), toString(mutations));
    EXPECT_EQ(exp_local_target_seq, local_target_seq);
    EXPECT_EQ(exp_local_target_qual, local_target_qual);
    EXPECT_EQ(exp_aligned_query_seq, aligned_query_seq);
    EXPECT_EQ(exp_aligned_query_qual, aligned_query_qual);
}


TEST(LocateMutations, InsertNearLeftEnd) {
/*
 *        quals:   12
 *      inserts:   GC
 *  insert locs:  /
 *        quals: abcdefghijklmnop
 * aligned+gaps: ATGCATGCATGCATGC
 *       target: ATGCATGCATGCATGC
 */

	int left_target_pos = 0;
	std::string query_bases = "AGCTGCATGCATGCATGC";
	std::string query_qual =  "a12bcdefghijklmnop";
	std::vector<BAM::CigarOp> cigar_data{{'M', 1},
	                                     {'I', 2},
	                                     {'M', 15}};
	std::vector<MdOp> md = {{MD_MATCH, 16, ""}};


	std::vector<Mutation> exp_mutations = {{0, 1, "GC", "12"}};
	std::string exp_local_target_seq  =  "ATGCATGCATGCATGC";
	std::string exp_local_target_qual =  "abcdefghijklmnop";
	std::string exp_aligned_query_seq =  "ATGCATGCATGCATGC";
	std::string exp_aligned_query_qual = "abcdefghijklmnop";

    std::vector<Mutation> mutations;
    std::string local_target_seq;
    std::string local_target_qual;
    std::string aligned_query_seq;
    std::string aligned_query_qual;
    boost::tie(mutations,
               local_target_seq,
               local_target_qual,
               aligned_query_seq,
               aligned_query_qual) = locateMutations(left_target_pos,
                                                     query_bases,
                                                     query_qual,
                                                     cigar_data,
                                                     md,
                                                     true);

    EXPECT_EQ(toString(exp_mutations), toString(mutations));
    EXPECT_EQ(exp_local_target_seq, local_target_seq);
    EXPECT_EQ(exp_local_target_qual, local_target_qual);
    EXPECT_EQ(exp_aligned_query_seq, aligned_query_seq);
    EXPECT_EQ(exp_aligned_query_qual, aligned_query_qual);
}


TEST(LocateMutations, InsertNearRightEnd) {
/*
 *        quals:                 12
 *      inserts:                 GC
 *  insert locs:                /
 *        quals: abcdefghijklmnop
 * aligned+gaps: ATGCATGCATGCATGC
 *       target: ATGCATGCATGCATGC
 */

	int left_target_pos = 0;
	std::string query_bases = "ATGCATGCATGCATGGCC";
	std::string query_qual =  "abcdefghijklmno12p";
	std::vector<BAM::CigarOp> cigar_data{{'M', 15},
	                                     {'I', 2},
	                                     {'M', 1}};
	std::vector<MdOp> md = {{MD_MATCH, 16, ""}};


	std::vector<Mutation> exp_mutations = {{14, 15, "GC", "12"}};
	std::string exp_local_target_seq =   "ATGCATGCATGCATGC";
	std::string exp_local_target_qual =  "abcdefghijklmnop";
	std::string exp_aligned_query_seq =  "ATGCATGCATGCATGC";
	std::string exp_aligned_query_qual = "abcdefghijklmnop";

    std::vector<Mutation> mutations;
    std::string local_target_seq;
    std::string local_target_qual;
    std::string aligned_query_seq;
    std::string aligned_query_qual;
    boost::tie(mutations,
               local_target_seq,
               local_target_qual,
               aligned_query_seq,
               aligned_query_qual) = locateMutations(left_target_pos,
                                                     query_bases,
                                                     query_qual,
                                                     cigar_data,
                                                     md,
                                                     true);

    EXPECT_EQ(toString(exp_mutations), toString(mutations));
    EXPECT_EQ(exp_local_target_seq, local_target_seq);
    EXPECT_EQ(exp_local_target_qual, local_target_qual);
    EXPECT_EQ(exp_aligned_query_seq, aligned_query_seq);
    EXPECT_EQ(exp_aligned_query_qual, aligned_query_qual);
}


TEST(LocateMutations, Complex) {
// Complex example with weird intervening inserts, soft-clipping, anything and everything
/*
 *        quals:       12
 *      inserts:       GG
 *  insert locs:      /
 *        quals:   abcd   efghi
 * aligned+gaps:   GCCC---CGCAT
 *       target: ATGCATGCATGCATGC
 */

	int left_target_pos = 2;
	std::string query_bases = "AAGCCGGCCGCATAA";
	std::string query_qual =  "HHabc12defghiHH";

	std::vector<BAM::CigarOp> cigar_data{{'S', 2},
	                                     {'M', 3},
	                                     {'I', 2},
	                                     {'M', 1},
	                                     {'D', 3},
	                                     {'M', 5},
	                                     {'S', 2}};
	std::vector<MdOp> md = {{MD_MATCH,    2, ""},
	                        {MD_MISMATCH, 1, "A"},
	                        {MD_MISMATCH, 1, "T"},
	                        {MD_DELETION, 3, "GCA"},
	                        {MD_MISMATCH, 1, "T"},
	                        {MD_MATCH,    4, ""}};

	std::vector<Mutation> exp_mutations = {{3, 5,  "C", "c"},
	                                       {4, 5,  "GG", "12"},
	                                       {4, 6,  "C", "d"},
	                                       {5, 9,  "", ""},
	                                       {8, 10, "C", "e"}};
	std::string exp_local_target_seq =   "GCATGCATGCAT";
	std::string exp_local_target_qual =  "abcd!!!efghi";
	std::string exp_aligned_query_seq =  "GCCC---CGCAT";
	std::string exp_aligned_query_qual = "abcd!!!efghi";

    std::vector<Mutation> mutations;
    std::string local_target_seq;
    std::string local_target_qual;
    std::string aligned_query_seq;
    std::string aligned_query_qual;
    boost::tie(mutations,
               local_target_seq,
               local_target_qual,
               aligned_query_seq,
               aligned_query_qual) = locateMutations(left_target_pos,
                                                     query_bases,
                                                     query_qual,
                                                     cigar_data,
                                                     md,
                                                     true);

    EXPECT_EQ(toString(exp_mutations), toString(mutations));
    EXPECT_EQ(exp_local_target_seq, local_target_seq);
    EXPECT_EQ(exp_local_target_qual, local_target_qual);
    EXPECT_EQ(exp_aligned_query_seq, aligned_query_seq);
    EXPECT_EQ(exp_aligned_query_qual, aligned_query_qual);
}


// Single Mutation toString()
TEST(Mutation, ToString) {
	Mutation m{5, 9, "CCT", "HHH"};
	EXPECT_NO_THROW(m.toString());
}

// Vector of Mutation toString()
TEST(MutationVector, ToString) {
	std::vector<Mutation> m{{1, 3,   "", ""},
	                        {2, 4,   "T", "H"},
	                        {8, 100, "", ""}};
	EXPECT_NO_THROW(toString(m));
}


// Ambiguous gap with two possible placements, right aligned
TEST(AdjustAmbiguousMutations, AmbigGapRightAligned) {
	int pos = 0;
	std::string local_target_seq =  "ATGGAT";
	std::string local_target_qual = "abc!de";
	std::string aligned_seq =  "ATG-AT";
	std::string aligned_qual = "abc!de";

	std::vector<Mutation> mutations{{2, 4, "", ""}};
	std::vector<Mutation> adjusted = adjustAmbiguousMutations(pos,
	                                                          local_target_seq,
	                                                          local_target_qual,
	                                                          aligned_seq,
	                                                          aligned_qual,
	                                                          mutations);
	std::vector<Mutation> expected{{1, 4, "G", "c"}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// Ambiguous gap with two possible placements, left aligned
TEST(AdjustAmbiguousMutations, AmbigGapLeftAligned) {
	int pos = 0;
	std::string local_target_seq = "ATGGAT";
	std::string local_target_qual= "ab!cde";
	std::string aligned_seq = "AT-GAT";
	std::string aligned_qual= "ab!cde";
	std::vector<Mutation> mutations{{1, 3, "", ""}};
	std::vector<Mutation> adjusted = adjustAmbiguousMutations(pos,
	                                                          local_target_seq,
	                                                          local_target_qual,
	                                                          aligned_seq,
	                                                          aligned_qual,
	                                                          mutations);
	std::vector<Mutation> expected{{1, 4, "G", "c"}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// Large ambiguous gap with 3 possible placements, right aligned
TEST(AdjustAmbiguousMutations, LargeAmbigGapRightAligned) {
	int pos = 0;
	std::string local_target_seq = "ATGGGGAT";
	std::string local_target_qual= "abcd!!ef";
	std::string aligned_seq =      "ATGG--AT";
	std::string aligned_qual =     "abcd!!ef";
	std::vector<Mutation> mutations{{3, 6, "", ""}};
	std::vector<Mutation> adjusted = adjustAmbiguousMutations(pos,
	                                                          local_target_seq,
	                                                          local_target_qual,
	                                                          aligned_seq,
	                                                          aligned_qual,
	                                                          mutations);
	std::vector<Mutation> expected{{1, 6, "GG", "cd"}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// Large ambiguous gap with 3 possible placements, left aligned
TEST(AdjustAmbiguousMutations, LargeAmbigGapLeftAligned) {
	int pos = 0;
	std::string local_target_seq = "ATGGGGAT";
	std::string local_target_qual= "ab!!cdef";
	std::string aligned_seq =      "AT--GGAT";
	std::string aligned_qual=      "ab!!cdef";
	std::vector<Mutation> mutations{{1, 4, "", ""}};
	std::vector<Mutation> adjusted = adjustAmbiguousMutations(pos,
	                                                          local_target_seq,
	                                                          local_target_qual,
	                                                          aligned_seq,
	                                                          aligned_qual,
	                                                          mutations);
	std::vector<Mutation> expected{{1, 6, "GG", "cd"}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// Test ambiguous gap adjacent to mismatch handling
TEST(AdjustAmbiguousMutations, AmbigGapAdjacentMismatch) {
	int pos = 0;
	std::string local_target_seq = "ATGGAT";
	std::string local_target_qual= "abc!de";
	std::string aligned_seq =      "ATG-CT";
	std::string aligned_qual =     "abc!de";
	// should be 3 possible placements, but merging would obscure this info
	// ATG-CT
	// ATGC-T
	// AT-GCT
	//
	// Treat as two separate mutations, handle merging in MutationCounter
	std::vector<Mutation> mutations{{2, 4, "", ""},
	                                {3, 5, "C", "d"}};
	std::vector<Mutation> adjusted = adjustAmbiguousMutations(pos,
	                                                          local_target_seq,
	                                                          local_target_qual,
	                                                          aligned_seq,
	                                                          aligned_qual,
	                                                          mutations);
	std::vector<Mutation> expected{{1, 4, "G", "c"},
	                               {3, 5, "C", "d"}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// ambiguous insert with two possible placements, right aligned
TEST(AdjustAmbiguousMutations, AmbigInsertRightAligned) {
	int pos = 0;
	std::string local_target_seq = "ATGAT";
	std::string local_target_qual= "abcde";
	std::string aligned_seq =      "ATGAT";
	std::string aligned_qual =     "abcde";
	std::vector<Mutation> mutations{{2, 3, "G", "1"}};
	std::vector<Mutation> adjusted = adjustAmbiguousMutations(pos,
	                                                          local_target_seq,
	                                                          local_target_qual,
	                                                          aligned_seq,
	                                                          aligned_qual,
	                                                          mutations);
	std::vector<Mutation> expected{{1, 3, "GG", "c1"}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// ambiguous insert with two possible placements, left aligned
TEST(AdjustAmbiguousMutations, AmbigInsertLeftAligned) {
	int pos = 0;
	std::string local_target_seq = "ATGAT";
	std::string local_target_qual= "abcde";
	std::string aligned_seq =      "ATGAT";
	std::string aligned_qual=      "abcde";
	std::vector<Mutation> mutations{{1, 2, "G", "1"}};
	std::vector<Mutation> adjusted = adjustAmbiguousMutations(pos,
	                                                          local_target_seq,
	                                                          local_target_qual,
	                                                          aligned_seq,
	                                                          aligned_qual,
	                                                          mutations);
	std::vector<Mutation> expected{{1, 3, "GG", "1c"}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// ambiguous gap near a different unambiguous gap on the right
// (not actually a likely output, since most aligners would group into a single gap of length 2)
TEST(AdjustAmbiguousMutations, AmbigGapWithUnambigGapOnRight) {
	int pos = 0;
	std::string local_target_seq = "ATGGATC";
    std::string local_target_qual= "ab!c!de";
	std::string aligned_seq =      "AT-G-TC";
	std::string aligned_qual =     "ab!c!de";
	std::vector<Mutation> mutations{{1, 3, "", ""},
	                                {3, 5, "", ""}};
	std::vector<Mutation> adjusted = adjustAmbiguousMutations(pos,
	                                                          local_target_seq,
	                                                          local_target_qual,
	                                                          aligned_seq,
	                                                          aligned_qual,
	                                                          mutations);
	std::vector<Mutation> expected{{1, 4, "G", "c"},
	                               {3, 5, "", ""}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// ambiguous gap near a different unambiguous gap on the left
TEST(AdjustAmbiguousMutations, AmbigGapWithUnambigGapOnLeft) {
	int pos = 0;
	std::string local_target_seq = "ATAGGTC";
	std::string local_target_qual= "ab!c!de";
	std::string aligned_seq =      "AT-G-TC";
	std::string aligned_qual =     "ab!c!de";
	std::vector<Mutation> mutations{{1, 3, "", ""},
	                                {3, 5, "", ""}};
	std::vector<Mutation> adjusted = adjustAmbiguousMutations(pos,
	                                                          local_target_seq,
	                                                          local_target_qual,
	                                                          aligned_seq,
	                                                          aligned_qual,
	                                                          mutations);
	std::vector<Mutation> expected{{1, 3, "", ""},
	                               {2, 5, "G", "c"}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// ambiguous gap near an unambiguous insertion on right
TEST(AdjustAmbiguousMutations, AmbigGapWithUnambigInsertOnRight) {
	int pos = 0;
	std::string local_target_seq = "ATGGATC";
	std::string local_target_qual= "ab!cdef";
	std::string aligned_seq =      "AT-GATC";
	std::string aligned_qual =     "ab!cdef";
	std::vector<Mutation> mutations{{1, 3, "", ""},
	                                {3, 4, "C", "1"}};
	std::vector<Mutation> adjusted = adjustAmbiguousMutations(pos,
	                                                          local_target_seq,
	                                                          local_target_qual,
	                                                          aligned_seq,
	                                                          aligned_qual,
	                                                          mutations);
	std::vector<Mutation> expected{{1, 4, "G", "c"},
	                               {3, 4, "C", "1"}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// ambiguous gap near an unambiguous insertion on right
TEST(AdjustAmbiguousMutations, AmbigGapWithUnambigInsertOnLeft) {
	int pos = 0;
	std::string local_target_seq = "ATAGGTC";
	std::string local_target_qual= "abcd!ef";
	std::string aligned_seq =      "ATAG-TC";
	std::string aligned_qual =     "abcd!ef";
	std::vector<Mutation> mutations{{1, 3, "C", "1"},
	                                {3, 5, "", ""}};
	std::vector<Mutation> adjusted = adjustAmbiguousMutations(pos,
	                                                          local_target_seq,
	                                                          local_target_qual,
	                                                          aligned_seq,
	                                                          aligned_qual,
	                                                          mutations);
	std::vector<Mutation> expected{{1, 3, "C", "1"},
	                               {2, 5, "G", "d"}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// What happens with a poorly aligned region with two ambiguous gaps nearby?
TEST(AdjustAmbiguousMutations, AmbigGapWithAmbigGapOnRight) {
	int pos = 0;
	std::string local_target_seq = "ATGGGTC";
	std::string local_target_qual= "ab!c!de";
	std::string aligned_seq =      "AT-G-TC";
	std::string aligned_qual=      "ab!c!de";
	
	std::vector<Mutation> mutations{{1, 3, "", ""},
	                                {3, 5, "", ""}};
	std::vector<Mutation> adjusted = adjustAmbiguousMutations(pos,
	                                                          local_target_seq,
	                                                          local_target_qual,
	                                                          aligned_seq,
	                                                          aligned_qual,
	                                                          mutations);

	// OLD BEHAVIOR: same nuc appended to both adjacent gaps (bad if doing sequence reconstruction)
	//std::vector<Mutation> expected{{1, 4, "G"},
	//                               {2, 5, "G"}};
	// NEW BEHAVIOR: chain adjacent ambiguous mutations into single mutations
	std::vector<Mutation> expected{{1, 5, "G", "c"}};

	EXPECT_EQ(toString(expected), toString(adjusted));
}


// ambiguous gap near sequence end on right
TEST(AdjustAmbiguousMutations, AmbigGapNearRightEnd) {
	int pos = 0;
	std::string local_target_seq = "ATGAA";
	std::string local_target_qual= "abc!d";
	std::string aligned_seq = 	   "ATG-A";
	std::string aligned_qual=      "abc!d";
	std::vector<Mutation> mutations{{2, 4, "", ""}};
	std::vector<Mutation> adjusted = adjustAmbiguousMutations(pos,
	                                                          local_target_seq,
	                                                          local_target_qual,
	                                                          aligned_seq,
	                                                          aligned_qual,
	                                                          mutations);
	std::vector<Mutation> expected{{2, 5, "A", "d"}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// ambiguous gap near sequence end on right
TEST(AdjustAmbiguousMutations, AmbigGapNearLeftEnd) {
	int pos = 0;
	std::string local_target_seq = "TTGCA";
	std::string local_target_qual= "a!bcd";
	std::string aligned_seq = 	   "T-GCA";
	std::string aligned_qual =     "a!bcd";
	std::vector<Mutation> mutations{{0, 2, "", ""}};
	std::vector<Mutation> adjusted = adjustAmbiguousMutations(pos,
	                                                          local_target_seq,
	                                                          local_target_qual,
	                                                          aligned_seq,
	                                                          aligned_qual,
	                                                          mutations);
	std::vector<Mutation> expected{{-1, 2, "T", "a"}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// ambiguous gap in read not aligned to first nuc in target
TEST(AdjustAmbiguousMutations, LargeAmbigGapLeftAlignedMissingFirstNuc) {
	int pos = 1;
	std::string local_target_seq = "ATGGGGAT";
	std::string local_target_qual= "ab!!cdef";
	std::string aligned_seq =      "AT--GGAT";
	std::string aligned_qual =     "ab!!cdef";
	std::vector<Mutation> mutations{{2, 5, "", ""}};
	std::vector<Mutation> adjusted = adjustAmbiguousMutations(pos,
	                                                          local_target_seq,
	                                                          local_target_qual,
	                                                          aligned_seq,
	                                                          aligned_qual,
	                                                          mutations);
	std::vector<Mutation> expected{{2, 7, "GG", "cd"}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


TEST(AdjustAmbiguousMutations, LargeAmbigGapLeftAlignedMissingFirstTwoNucs) {
	int pos = 2;
	std::string local_target_seq = "ATGGGGAT";
    std::string local_target_qual= "ab!!cdef";
	std::string aligned_seq =      "AT--GGAT";
	std::string aligned_qual=      "ab!!cdef";
	std::vector<Mutation> mutations{{3, 6, "", ""}};
	std::vector<Mutation> adjusted = adjustAmbiguousMutations(pos,
	                                                          local_target_seq,
	                                                          local_target_qual,
	                                                          aligned_seq,
	                                                          aligned_qual,
	                                                          mutations);
	std::vector<Mutation> expected{{3, 8, "GG", "cd"}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


/*
 * Single-nuc dels:
 *  "A-", "T-", "G-", "C-",
 * Single-nuc inserts:
 *  "-A", "-T", "-G", "-C",
 * Single-nuc mismatches:
 *  "AT", "AG", "AC",
 *  "TA", "TG", "TC",
 *  "GA", "GT", "GC",
 *  "CA", "CT", "CG",
 * Others:
 *  "multinuc_deletion",
 *  "multinuc_insertion",
 *  "multinuc_mismatch",
 *
 *  "complex_deletion",
 *  "complex_insertion"
 */
std::string local_target_seq = "ATGCATGC";
std::string local_target_qual= "abcdefgh";
std::vector<Mutation> m{{3, 5, "", ""},
                        {4, 6, "", ""},
                        {5, 7, "", ""},
                        {6, 8, "", ""},

                        {1, 2, "A", "1"},
                        {1, 2, "T", "1"},
                        {1, 2, "G", "1"},
                        {1, 2, "C", "1"},

                        {3, 5, "T", "1"},
                        {3, 5, "G", "1"},
                        {3, 5, "C", "1"},
                        {4, 6, "A", "1"},
                        {4, 6, "G", "1"},
                        {4, 6, "C", "1"},
                        {5, 7, "A", "1"},
                        {5, 7, "T", "1"},
                        {5, 7, "C", "1"},
                        {6, 8, "A", "1"},
                        {6, 8, "T", "1"},
                        {6, 8, "G", "1"},

                        {3, 6, "", ""},
                        {1, 2, "AA", "12"},
                        {1, 4, "TG", "12"}
};
std::vector<std::string> expected{"A-", "T-", "G-", "C-",
                                  "-A", "-T", "-G", "-C",
                                  "AT", "AG", "AC",
                                  "TA", "TG", "TC",
                                  "GA", "GT", "GC",
                                  "CA", "CT", "CG",
                                  "multinuc_deletion",
                                  "multinuc_insertion",
                                  "multinuc_mismatch",
};

TEST(MutationClassify, All) {
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_EQ(expected[i], m[i].classify(local_target_seq, 0));
    }
}

TEST(MutationClassify, NonzeroTargetPos) {
    std::string local_target_seq = "TGCATGC";
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_EQ(expected[i], m[i].classify(local_target_seq, 1));
    }
}


TEST(ParseSamMutations, Simple) {
    std::vector<std::string> fields = {
            "M00236:2:000000000-A21YG:1:1106:15774:10066",
            "16",
            "TPP_riboswitch",
            "1",
            "42",
            "8S20M2I79M1D37M5S",
            "*",
            "0",
            "0",
            "ATCAGAACGGCCTTCGGGCCAAGGACTCAAGGACTCCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTCGGGCG",
            "*:CECC>>EC828'AEEC*::*C???*CEGE?C>GDGGGGGGGGGGCGEEGGGGGGGGGEB>GGGGGGGGGGGGGGGHGGGGGGGGGEGGGGGGHHHHHHHHHHHHHHDHHHHHHHHIIHHHIIIHHHHGGGGGGDDEDDDEDBAAAAAAA",
            "AS:i:239",
            "XN:i:0",
            "XM:i:3",
            "XO:i:2",
            "XG:i:3",
            "NM:i:6",
            "MD:Z:22G0G1G73^C37",
            "YT:Z:UU"
    };

	std::string read_id;
    int left_target_pos;
    int right_target_pos;
    std::string local_target_seq;
    std::string local_target_qual;
    std::vector<Mutation> adjusted_mutations;
    EXPECT_NO_THROW(
            boost::tie(read_id,
					   left_target_pos,
                       right_target_pos,
                       local_target_seq,
                       local_target_qual,
                       adjusted_mutations) = mutation_parser::parseSamMutations(fields);
    );
/*    std::cout << "Parsed SAM mutations: ";
    std::cout << mutation::serializeReadInfo(left_target_pos,
                                             right_target_pos,
                                             local_target,
                                             adjusted_mutations);*/

}


// parseMutations() - grab read from bam file, check output matches expected

// TODO: test reconstructed target sequences from a number of reads mapped to the same location are identical? cover more cases than just the handful here

// bunch of exception handling conditions?

// Test serialization/deserialization

// TODO: do some random stack stops in this module and find rate limiting steps
