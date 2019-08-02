/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2018 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#include "gtest/gtest.h"
#include "MutationParser.cpp"
#include "util.h"

namespace BF = boost::filesystem;
using namespace mutation;
using namespace mutation_parser;
using namespace mutation_parser::detail;
using namespace util;

std::string FILEPATH = __FILE__;
std::string BASEPATH = "";

BF::path getTestFileDir() {
    BF::path filedir;
    if (BASEPATH == "") {
        // data file location if explicit path to main shapemapper directory
        // not provided by test runner (assumes software location has not
        // changed since compilation)
        filedir = BF::path(FILEPATH).parent_path() / "files";
    } else {
        filedir = BF::path(BASEPATH) / "internals" / "cpp-src" / "test" / "files";
        //std::cout << "Using BASEPATH " << BASEPATH << " passed by argument" << std::endl;
    }
    BF::create_directory(filedir / "tmp");
    return filedir;
}


// NOTE: for now, can't use debug_outname in tests, since ofstream mutation::debug_out
// is shared across everything importing the mutation namespace

void processMutations_wrapper(const std::string &line) {
	Read r = parseTestRead(line);
	std::vector<Read> reads;
	reads.push_back(r);

	Read processed_read = processMutations(reads,
										   1, //direction
										   false, // right_align_ambig_dels
										   false, // right_align_ambig_ins
										   6, // max_internal_match
										   0, // min_qual
										   1, // exclude_3prime
										   "", // mutation_type (empty = count all mutations)
										   false, // variant_mode
										   false, // trim_primers
										   PrimerPair(), // primer_pair
										   false); // print debug info
}

boost::tuple<std::string, std::string, std::string>
processMutations_wrapper_min_qual(const std::string &line,
								  const int min_qual) {
	Read r = parseTestRead(line);
	std::vector<Read> reads;
	reads.push_back(r);

	Read processed_read = processMutations(reads,
										   1, //direction
										   false, // right_align_ambig_dels
										   false, // right_align_ambig_ins
										   0, // max_internal_match
										   min_qual, // min_qual
										   1, // exclude_3prime
										   "", // mutation_type (empty = count all mutations)
										   false, // variant_mode
										   false, // trim_primers
										   PrimerPair(), // primer_pair
										   false); // print debug info

	std::string cm = toString(processed_read.mutations);
	std::string led = toString(processed_read.depth);
	std::string lec = toString(processed_read.count);
	return boost::make_tuple(cm, led, lec);
}

void processMutations_wrapper_exclude_3prime(const std::string &line,
											 const int exclude_3prime) {
	Read r = parseTestRead(line);
	std::vector<Read> reads;
	reads.push_back(r);

	Read processed_read = processMutations(reads,
										   1, //direction
										   false, // right_align_ambig_dels
										   false, // right_align_ambig_ins
										   6, // max_internal_match
										   0, // min_qual
										   exclude_3prime, // exclude_3prime
										   "", // mutation_type (empty = count all mutations)
										   false, // variant_mode
										   false, // trim_primers
										   PrimerPair(), // primer_pair
										   false); // print debug info
}

//----------------------------------------------------------------------------------------------


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

	int left = 0;
	std::string query_bases = "ATGCATGCATGCATGC";
	std::string query_qual =  "ABCDEFGHIJKLMNOP";
	std::vector<CigarOp> cigar_data{{'M', 16}};
	std::vector<MdOp> md = {{MD_MATCH, 16, ""}};


	std::vector<Mutation> exp_mutations = {};
	std::string exp_seq = "ATGCATGCATGCATGC";
	std::string exp_qual = "ABCDEFGHIJKLMNOP";
	std::string exp_aligned_query_seq = "ATGCATGCATGCATGC";
	std::string exp_aligned_query_qual = "ABCDEFGHIJKLMNOP";


	std::vector<Mutation> mutations;
	std::string seq;
	std::string qual;
	std::string aligned_query_seq;
	std::string aligned_query_qual;
	boost::tie(mutations,
			   seq,
			   qual,
			   aligned_query_seq,
			   aligned_query_qual) = locateMutations(left,
													 query_bases,
													 query_qual,
													 cigar_data,
													 md,
													 true);

	EXPECT_EQ(toString(exp_mutations), toString(mutations));
	EXPECT_EQ(exp_seq, seq);
	EXPECT_EQ(exp_qual, qual);
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

	int left = 0;
	std::string query_bases = "ATGCATGCGTGCATGC";
	std::string query_qual =  "ABCDEFGHIJKLMNOP";
	std::vector<CigarOp> cigar_data{{'M', 16}};
	std::vector<MdOp> md = {{MD_MATCH,    8, ""},
							{MD_MISMATCH, 1, "A"},
							{MD_MATCH,    7, ""}};


	std::vector<Mutation> exp_mutations = {{7, 9, "G", "I", ""}};
	std::string exp_seq =  "ATGCATGCATGCATGC";
	std::string exp_qual = "ABCDEFGHIJKLMNOP";
	std::string exp_aligned_query_seq = "ATGCATGCGTGCATGC";
	std::string exp_aligned_query_qual = "ABCDEFGHIJKLMNOP";


	std::vector<Mutation> mutations;
	std::string seq;
	std::string qual;
	std::string aligned_query_seq;
	std::string aligned_query_qual;
	boost::tie(mutations,
			   seq,
			   qual,
			   aligned_query_seq,
			   aligned_query_qual) = locateMutations(left,
													 query_bases,
													 query_qual,
													 cigar_data,
													 md,
													 true);

	EXPECT_EQ(toString(exp_mutations), toString(mutations));
	EXPECT_EQ(exp_seq, seq);
	EXPECT_EQ(exp_qual, qual);
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

	int left = 0;
	std::string query_bases = "ATCATGCAAAATGCATGC";
	std::string query_qual =  "abcdefgh123ijklmno";
	std::vector<CigarOp> cigar_data{{'M', 2},
									{'D', 1},
									{'M', 6},
									{'I', 3},
									{'M', 7}};
	std::vector<MdOp> md = {{MD_MATCH,    2,  ""},
							{MD_DELETION, 1,  "G"},
							{MD_MATCH,    13, ""}};


	std::vector<Mutation> exp_mutations = {{1, 3, "", "", ""},
										   {8, 9, "AAA", "123", ""}};
	std::string exp_seq =  "ATGCATGCATGCATGC";
	std::string exp_qual = "ab!cdefghijklmno";
	std::string exp_aligned_query_seq  = "AT-CATGCATGCATGC";
	std::string exp_aligned_query_qual = "ab!cdefghijklmno";

	std::vector<Mutation> mutations;
	std::string seq;
	std::string qual;
	std::string aligned_query_seq;
	std::string aligned_query_qual;
	boost::tie(mutations,
			   seq,
			   qual,
			   aligned_query_seq,
			   aligned_query_qual) = locateMutations(left,
													 query_bases,
													 query_qual,
													 cigar_data,
													 md,
													 true);

	EXPECT_EQ(toString(exp_mutations), toString(mutations));
	EXPECT_EQ(exp_seq, seq);
	EXPECT_EQ(exp_qual, qual);
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

	int left = 0;
	std::string query_bases = "ATGAAACATGCATGATGC";
	std::string query_qual =  "abc123defghijklmno";
	std::vector<CigarOp> cigar_data{{'M', 3},
									{'I', 3},
									{'M', 8},
									{'D', 1},
									{'M', 4}};
	std::vector<MdOp> md = {{MD_MATCH,    11, ""},
							{MD_DELETION, 1,  "C"},
							{MD_MATCH,    4,  ""}};


	std::vector<Mutation> exp_mutations = {{2,  3,  "AAA", "123", ""},
										   {10, 12, "", "", ""}};
	std::string exp_seq =   "ATGCATGCATGCATGC";
	std::string exp_qual =  "abcdefghijk!lmno";
	std::string exp_aligned_query_seq =  "ATGCATGCATG-ATGC";
	std::string exp_aligned_query_qual = "abcdefghijk!lmno";

	std::vector<Mutation> mutations;
	std::string seq;
	std::string qual;
	std::string aligned_query_seq;
	std::string aligned_query_qual;
	boost::tie(mutations,
			   seq,
			   qual,
			   aligned_query_seq,
			   aligned_query_qual) = locateMutations(left,
													 query_bases,
													 query_qual,
													 cigar_data,
													 md,
													 true);

	EXPECT_EQ(toString(exp_mutations), toString(mutations));
	EXPECT_EQ(exp_seq, seq);
	EXPECT_EQ(exp_qual, qual);
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

	int left = 1;
	std::string query_bases = "GGGGGTGCATGCGTGCATGCGGGGG";
	std::string query_qual =  "HHHHHabcdefghijklmnoHHHHH";
	std::vector<CigarOp> cigar_data{{'S', 5},
									{'M', 15},
									{'S', 5}};
	std::vector<MdOp> md = {{MD_MATCH,    7, ""},
							{MD_MISMATCH, 1, "A"},
							{MD_MATCH,    7, ""}};


	std::vector<Mutation> exp_mutations = {{7, 9, "G", "h", ""}};
	std::string exp_seq =   "TGCATGCATGCATGC";
	std::string exp_qual =  "abcdefghijklmno";
	std::string exp_aligned_query_seq =  "TGCATGCGTGCATGC";
	std::string exp_aligned_query_qual = "abcdefghijklmno";


	std::vector<Mutation> mutations;
	std::string seq;
	std::string qual;
	std::string aligned_query_seq;
	std::string aligned_query_qual;
	boost::tie(mutations,
			   seq,
			   qual,
			   aligned_query_seq,
			   aligned_query_qual) = locateMutations(left,
													 query_bases,
													 query_qual,
													 cigar_data,
													 md,
													 true);

	EXPECT_EQ(toString(exp_mutations), toString(mutations));
	EXPECT_EQ(exp_seq, seq);
	EXPECT_EQ(exp_qual, qual);
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

	int left = 0;
	std::string query_bases = "AGCTGCATGCATGCATGC";
	std::string query_qual =  "a12bcdefghijklmnop";
	std::vector<CigarOp> cigar_data{{'M', 1},
									{'I', 2},
									{'M', 15}};
	std::vector<MdOp> md = {{MD_MATCH, 16, ""}};


	std::vector<Mutation> exp_mutations = {{0, 1, "GC", "12", ""}};
	std::string exp_seq  =  "ATGCATGCATGCATGC";
	std::string exp_qual =  "abcdefghijklmnop";
	std::string exp_aligned_query_seq =  "ATGCATGCATGCATGC";
	std::string exp_aligned_query_qual = "abcdefghijklmnop";

	std::vector<Mutation> mutations;
	std::string seq;
	std::string qual;
	std::string aligned_query_seq;
	std::string aligned_query_qual;
	boost::tie(mutations,
			   seq,
			   qual,
			   aligned_query_seq,
			   aligned_query_qual) = locateMutations(left,
													 query_bases,
													 query_qual,
													 cigar_data,
													 md,
													 true);

	EXPECT_EQ(toString(exp_mutations), toString(mutations));
	EXPECT_EQ(exp_seq, seq);
	EXPECT_EQ(exp_qual, qual);
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

	int left = 0;
	std::string query_bases = "ATGCATGCATGCATGGCC";
	std::string query_qual =  "abcdefghijklmno12p";
	std::vector<CigarOp> cigar_data{{'M', 15},
									{'I', 2},
									{'M', 1}};
	std::vector<MdOp> md = {{MD_MATCH, 16, ""}};


	std::vector<Mutation> exp_mutations = {{14, 15, "GC", "12", ""}};
	std::string exp_seq =   "ATGCATGCATGCATGC";
	std::string exp_qual =  "abcdefghijklmnop";
	std::string exp_aligned_query_seq =  "ATGCATGCATGCATGC";
	std::string exp_aligned_query_qual = "abcdefghijklmnop";

	std::vector<Mutation> mutations;
	std::string seq;
	std::string qual;
	std::string aligned_query_seq;
	std::string aligned_query_qual;
	boost::tie(mutations,
			   seq,
			   qual,
			   aligned_query_seq,
			   aligned_query_qual) = locateMutations(left,
													 query_bases,
													 query_qual,
													 cigar_data,
													 md,
													 true);

	EXPECT_EQ(toString(exp_mutations), toString(mutations));
	EXPECT_EQ(exp_seq, seq);
	EXPECT_EQ(exp_qual, qual);
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

	int left = 2;
	std::string query_bases = "AAGCCGGCCGCATAA";
	std::string query_qual =  "HHabc12defghiHH";

	std::vector<CigarOp> cigar_data{{'S', 2},
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

	std::vector<Mutation> exp_mutations = {{3, 5,  "C", "c", ""},
										   {4, 5,  "GG", "12", ""},
										   {4, 6,  "C", "d", ""},
										   {5, 9,  "", "", ""},
										   {8, 10, "C", "e", ""}};
	std::string exp_seq =   "GCATGCATGCAT";
	std::string exp_qual =  "abcd!!!efghi";
	std::string exp_aligned_query_seq =  "GCCC---CGCAT";
	std::string exp_aligned_query_qual = "abcd!!!efghi";

	std::vector<Mutation> mutations;
	std::string seq;
	std::string qual;
	std::string aligned_query_seq;
	std::string aligned_query_qual;
	boost::tie(mutations,
			   seq,
			   qual,
			   aligned_query_seq,
			   aligned_query_qual) = locateMutations(left,
													 query_bases,
													 query_qual,
													 cigar_data,
													 md,
													 true);

	EXPECT_EQ(toString(exp_mutations), toString(mutations));
	EXPECT_EQ(exp_seq, seq);
	EXPECT_EQ(exp_qual, qual);
	EXPECT_EQ(exp_aligned_query_seq, aligned_query_seq);
	EXPECT_EQ(exp_aligned_query_qual, aligned_query_qual);
}


// Single Mutation toString()
TEST(Mutation, ToString) {
	Mutation m{5, 9, "CCT", "HHH", ""};
	EXPECT_NO_THROW(m.toString());
}

// Vector of Mutation toString()
TEST(MutationVector, ToString) {
	std::vector<Mutation> m{{1, 3,   "", "", ""},
							{2, 4,   "T", "H", ""},
							{8, 100, "", "", ""}};
	EXPECT_NO_THROW(toString(m));
}


// Ambiguous gap with two possible placements, right aligned
TEST(IdentifyAmbiguousMutations, AmbigGapRightAligned) {
	int pos = 0;
	std::string seq =  "ATGGAT";
	std::string qual = "abc!de";
	std::string aligned_seq =  "ATG-AT";
	std::string aligned_qual = "abc!de";

	std::vector<Mutation> mutations{{2, 4, "", ""}};
	std::vector<Mutation> adjusted = identifyAmbiguousMutations(pos,
																seq,
																qual,
																aligned_seq,
																aligned_qual,
																mutations);
	std::vector<Mutation> expected{{1, 4, "G", "c", ""}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// Ambiguous gap with two possible placements, left aligned
TEST(IdentifyAmbiguousMutations, AmbigGapLeftAligned) {
	int pos = 0;
	std::string seq = "ATGGAT";
	std::string qual= "ab!cde";
	std::string aligned_seq = "AT-GAT";
	std::string aligned_qual= "ab!cde";
	std::vector<Mutation> mutations{{1, 3, "", "", ""}};
	std::vector<Mutation> adjusted = identifyAmbiguousMutations(pos,
																seq,
																qual,
																aligned_seq,
																aligned_qual,
																mutations);
	std::vector<Mutation> expected{{1, 4, "G", "c", ""}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// Large ambiguous gap with 3 possible placements, right aligned
TEST(IdentifyAmbiguousMutations, LargeAmbigGapRightAligned) {
	int pos = 0;
	std::string seq = "ATGGGGAT";
	std::string qual= "abcd!!ef";
	std::string aligned_seq =      "ATGG--AT";
	std::string aligned_qual =     "abcd!!ef";
	std::vector<Mutation> mutations{{3, 6, "", "", ""}};
	std::vector<Mutation> adjusted = identifyAmbiguousMutations(pos,
																seq,
																qual,
																aligned_seq,
																aligned_qual,
																mutations);
	std::vector<Mutation> expected{{1, 6, "GG", "cd", ""}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// Large ambiguous gap with 3 possible placements, left aligned
TEST(IdentifyAmbiguousMutations, LargeAmbigGapLeftAligned) {
	int pos = 0;
	std::string seq = "ATGGGGAT";
	std::string qual= "ab!!cdef";
	std::string aligned_seq =      "AT--GGAT";
	std::string aligned_qual=      "ab!!cdef";
	std::vector<Mutation> mutations{{1, 4, "", "", ""}};
	std::vector<Mutation> adjusted = identifyAmbiguousMutations(pos,
																seq,
																qual,
																aligned_seq,
																aligned_qual,
																mutations);
	std::vector<Mutation> expected{{1, 6, "GG", "cd"}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// Test ambiguous gap adjacent to mismatch handling
TEST(IdentifyAmbiguousMutations, AmbigGapAdjacentMismatch) {
	int pos = 0;
	std::string seq = "ATGGAT";
	std::string qual= "abc!de";
	std::string aligned_seq =      "ATG-CT";
	std::string aligned_qual =     "abc!de";
	// should be 3 possible placements, but merging would obscure this info
	// ATG-CT
	// ATGC-T
	// AT-GCT
	//
	// Treat as two separate mutations, handle merging in MutationCounter
	std::vector<Mutation> mutations{{2, 4, "", "", ""},
									{3, 5, "C", "d", ""}};
	std::vector<Mutation> adjusted = identifyAmbiguousMutations(pos,
																seq,
																qual,
																aligned_seq,
																aligned_qual,
																mutations);
	std::vector<Mutation> expected{{1, 4, "G", "c", ""},
								   {3, 5, "C", "d", ""}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// ambiguous insert with two possible placements, right aligned
TEST(IdentifyAmbiguousMutations, AmbigInsertRightAligned) {
	int pos = 0;
	std::string seq = "ATGAT";
	std::string qual= "abcde";
	std::string aligned_seq =      "ATGAT";
	std::string aligned_qual =     "abcde";
	std::vector<Mutation> mutations{{2, 3, "G", "1", ""}};
	std::vector<Mutation> adjusted = identifyAmbiguousMutations(pos,
																seq,
																qual,
																aligned_seq,
																aligned_qual,
																mutations);
	std::vector<Mutation> expected{{1, 3, "GG", "c1", ""}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// ambiguous insert with two possible placements, left aligned
TEST(IdentifyAmbiguousMutations, AmbigInsertLeftAligned) {
	int pos = 0;
	std::string seq = "ATGAT";
	std::string qual= "abcde";
	std::string aligned_seq =      "ATGAT";
	std::string aligned_qual=      "abcde";
	std::vector<Mutation> mutations{{1, 2, "G", "1", ""}};
	std::vector<Mutation> adjusted = identifyAmbiguousMutations(pos,
																seq,
																qual,
																aligned_seq,
																aligned_qual,
																mutations);
	std::vector<Mutation> expected{{1, 3, "GG", "1c", ""}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// ambiguous gap near a different unambiguous gap on the right
// (not actually a likely output, since most aligners would group into a single gap of length 2)
TEST(IdentifyAmbiguousMutations, AmbigGapWithUnambigGapOnRight) {
	int pos = 0;
	std::string seq = "ATGGATC";
	std::string qual= "ab!c!de";
	std::string aligned_seq =      "AT-G-TC";
	std::string aligned_qual =     "ab!c!de";
	std::vector<Mutation> mutations{{1, 3, "", "", ""},
									{3, 5, "", "", ""}};
	std::vector<Mutation> adjusted = identifyAmbiguousMutations(pos,
																seq,
																qual,
																aligned_seq,
																aligned_qual,
																mutations);
	std::vector<Mutation> expected{{1, 4, "G", "c", ""},
								   {3, 5, "", "", ""}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// ambiguous gap near a different unambiguous gap on the left
TEST(IdentifyAmbiguousMutations, AmbigGapWithUnambigGapOnLeft) {
	int pos = 0;
	std::string seq = "ATAGGTC";
	std::string qual= "ab!c!de";
	std::string aligned_seq =      "AT-G-TC";
	std::string aligned_qual =     "ab!c!de";
	std::vector<Mutation> mutations{{1, 3, "", "", ""},
									{3, 5, "", "", ""}};
	std::vector<Mutation> adjusted = identifyAmbiguousMutations(pos,
																seq,
																qual,
																aligned_seq,
																aligned_qual,
																mutations);
	std::vector<Mutation> expected{{1, 3, "", "", ""},
								   {2, 5, "G", "c", ""}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// ambiguous gap near an unambiguous insertion on right
TEST(IdentifyAmbiguousMutations, AmbigGapWithUnambigInsertOnRight) {
	int pos = 0;
	std::string seq = "ATGGATC";
	std::string qual= "ab!cdef";
	std::string aligned_seq =      "AT-GATC";
	std::string aligned_qual =     "ab!cdef";
	std::vector<Mutation> mutations{{1, 3, "", "", ""},
									{3, 4, "C", "1", ""}};
	std::vector<Mutation> adjusted = identifyAmbiguousMutations(pos,
																seq,
																qual,
																aligned_seq,
																aligned_qual,
																mutations);
	std::vector<Mutation> expected{{1, 4, "G", "c", ""},
								   {3, 4, "C", "1", ""}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// ambiguous gap near an unambiguous insertion on right
TEST(IdentifyAmbiguousMutations, AmbigGapWithUnambigInsertOnLeft) {
	int pos = 0;
	std::string seq = "ATAGGTC";
	std::string qual= "abcd!ef";
	std::string aligned_seq =      "ATAG-TC";
	std::string aligned_qual =     "abcd!ef";
	std::vector<Mutation> mutations{{1, 3, "C", "1", ""},
									{3, 5, "", "", ""}};
	std::vector<Mutation> adjusted = identifyAmbiguousMutations(pos,
																seq,
																qual,
																aligned_seq,
																aligned_qual,
																mutations);
	std::vector<Mutation> expected{{1, 3, "C", "1", ""},
								   {2, 5, "G", "d", ""}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// What happens with a poorly aligned region with two ambiguous gaps nearby?
TEST(IdentifyAmbiguousMutations, AmbigGapWithAmbigGapOnRight) {
	int pos = 0;
	std::string seq = "ATGGGTC";
	std::string qual= "ab!c!de";
	std::string aligned_seq =      "AT-G-TC";
	std::string aligned_qual=      "ab!c!de";

	std::vector<Mutation> mutations{{1, 3, "", "", ""},
									{3, 5, "", "", ""}};
	std::vector<Mutation> adjusted = identifyAmbiguousMutations(pos,
																seq,
																qual,
																aligned_seq,
																aligned_qual,
																mutations);

	// OLD BEHAVIOR: same nuc appended to both adjacent gaps (bad if doing sequence reconstruction)
	//std::vector<Mutation> expected{{1, 4, "G"},
	//                               {2, 5, "G"}};
	// NEW BEHAVIOR: chain adjacent ambiguous mutations into single mutations
	std::vector<Mutation> expected{{1, 5, "G", "c", ""}};

	EXPECT_EQ(toString(expected), toString(adjusted));
}


// ambiguous gap near sequence end on right
TEST(IdentifyAmbiguousMutations, AmbigGapNearRightEnd) {
	int pos = 0;
	std::string seq = "ATGAA";
	std::string qual= "abc!d";
	std::string aligned_seq = 	   "ATG-A";
	std::string aligned_qual=      "abc!d";
	std::vector<Mutation> mutations{{2, 4, "", "", ""}};
	std::vector<Mutation> adjusted = identifyAmbiguousMutations(pos,
																seq,
																qual,
																aligned_seq,
																aligned_qual,
																mutations);
	std::vector<Mutation> expected{{2, 5, "A", "d", ""}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// ambiguous gap near sequence end on right
TEST(IdentifyAmbiguousMutations, AmbigGapNearLeftEnd) {
	int pos = 0;
	std::string seq = "TTGCA";
	std::string qual= "a!bcd";
	std::string aligned_seq = 	   "T-GCA";
	std::string aligned_qual =     "a!bcd";
	std::vector<Mutation> mutations{{0, 2, "", "", ""}};
	std::vector<Mutation> adjusted = identifyAmbiguousMutations(pos,
																seq,
																qual,
																aligned_seq,
																aligned_qual,
																mutations);
	std::vector<Mutation> expected{{-1, 2, "T", "a", ""}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


// ambiguous gap in read not aligned to first nuc in target
TEST(IdentifyAmbiguousMutations, LargeAmbigGapLeftAlignedMissingFirstNuc) {
	int pos = 1;
	std::string seq = "ATGGGGAT";
	std::string qual= "ab!!cdef";
	std::string aligned_seq =      "AT--GGAT";
	std::string aligned_qual =     "ab!!cdef";
	std::vector<Mutation> mutations{{2, 5, "", "", ""}};
	std::vector<Mutation> adjusted = identifyAmbiguousMutations(pos,
																seq,
																qual,
																aligned_seq,
																aligned_qual,
																mutations);
	std::vector<Mutation> expected{{2, 7, "GG", "cd", ""}};
	EXPECT_EQ(toString(expected), toString(adjusted));
}


TEST(IdentifyAmbiguousMutations, LargeAmbigGapLeftAlignedMissingFirstTwoNucs) {
	int pos = 2;
	std::string seq = "ATGGGGAT";
	std::string qual= "ab!!cdef";
	std::string aligned_seq =      "AT--GGAT";
	std::string aligned_qual=      "ab!!cdef";
	std::vector<Mutation> mutations{{3, 6, "", "", ""}};
	std::vector<Mutation> adjusted = identifyAmbiguousMutations(pos,
																seq,
																qual,
																aligned_seq,
																aligned_qual,
																mutations);
	std::vector<Mutation> expected{{3, 8, "GG", "cd", ""}};
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
std::string seq = "ATGCATGC";
std::string qual= "abcdefgh";
std::vector<Mutation> m{{3, 5, "", "", ""},
						{4, 6, "", "", ""},
						{5, 7, "", "", ""},
						{6, 8, "", "", ""},

						{1, 2, "A", "1", ""},
						{1, 2, "T", "1", ""},
						{1, 2, "G", "1", ""},
						{1, 2, "C", "1", ""},

						{3, 5, "T", "1", ""},
						{3, 5, "G", "1", ""},
						{3, 5, "C", "1", ""},
						{4, 6, "A", "1", ""},
						{4, 6, "G", "1", ""},
						{4, 6, "C", "1", ""},
						{5, 7, "A", "1", ""},
						{5, 7, "T", "1", ""},
						{5, 7, "C", "1", ""},
						{6, 8, "A", "1", ""},
						{6, 8, "T", "1", ""},
						{6, 8, "G", "1", ""},

						{3, 6, "", "", ""},
						{1, 2, "AA", "12", ""},
						{1, 4, "TG", "12", ""}
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
		EXPECT_EQ(expected[i], m[i].classify(seq, 0));
	}
}

TEST(MutationClassify, NonzeroTargetPos) {
	std::string seq = "TGCATGC";
	for (int i = 0; i < expected.size(); ++i) {
		EXPECT_EQ(expected[i], m[i].classify(seq, 1));
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

	EXPECT_NO_THROW(
			Read r = mutation_parser::parseSamFields(fields, 30, true);
	);
/*    std::cout << "Parsed SAM mutations: ";
    std::cout << mutation::serializeReadInfo(left,
                                             right,
                                             local_target,
                                             adjusted_mutations);*/

}


TEST(FlagsToBits, UnpairedMappedReverseStrand) {
	std::string field = "16";
	std::bitset<12> flags;
	EXPECT_NO_THROW(
			flags = mutation_parser::detail::flagsToBits(field);
	);
	//std::cout << "Unpaired mapped reverse strand bit flags: ";
	//std::cout << flags << std::endl;
	//std::cout << "reverse strand flag 5: " << flags[4] << std::endl;
}

// FIXME: generate some unmerged paired mapped SAM reads with ribosome data
/*
TEST(ParseFlags, PairedR1Mapped) {
    std::string field = "16";
    std::bitset<11> flags;
    EXPECT_NO_THROW(
            flags = mutation_parser::detail::parseFlags(field);
    );
    std::cout << "Paired R1 mapped bit flags: ";
    std::cout << flags;
}

TEST(ParseFlags, PairedR2Mapped) {
    std::string field = "16";
    std::bitset<11> flags;
    EXPECT_NO_THROW(
            flags = mutation_parser::detail::parseFlags(field);
    );
    std::cout << "Paired R1 mapped bit flags: ";
    std::cout << flags;
}*/

// FIXME: quickly check paired end parsing -
// write some example R1 R2 pairs mapped concordantly / discordantly / unmapped to dummy SAM file,
// then parse with
/* parseSAM(const std::string &filename,
const std::string &outname,
bool paired = false,
unsigned int min_mapq = DEFAULT_MIN_MAPQ,
bool warn_on_no_mapped = false)*/

TEST(FlagsToBits, MateUnmapped) {
	std::bitset<12> exp_flags("000001011001"); // right-most is zeroeth
	std::bitset<12> flags = flagsToBits("89");
	//std::cout << "\nPaired mate unmapped:\n";
	//std::cout << "expected flag[0]: " << bool(exp_flags[0]) << std::endl;
	//std::cout << "  actual flag[0]: " << bool(flags[0]) << std::endl;
	//std::cout << "expected flags: " << exp_flags.to_string() << std::endl;
	//std::cout << "  actual flags: " << flags.to_string() << "\n" << std::endl;
	EXPECT_EQ(flags[0], true);
	EXPECT_EQ(exp_flags.to_string(), flags.to_string());
}


////////////////////////////////////////////////////////////////////////////////////////////

TEST(Debug, VectorRangeCheckCrash) {

	std::string sample_reads = R"(M01228:25:000000000-A1CW0:1:1101:11549:3441	89	TPP	16	44	54M1D11M1D55M5S	=	16	0	GACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGATAATGCCAGTTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTCTCACA	:CGGGGE?CGGGECC>EGEEGGECEEEGGC=AEEC@>GEED==,EEGFFF;GFGFDFHFDFDBHDFHHHHFHHGHHFHEBHCHCEFCHHHHIIIIIIIIHHHHGGGGGGDDDDDDEDBBB?????	AS:i:224	XN:i:0	XM:i:1	XO:i:2	XG:i:2	NM:i:3	MD:Z:54^G11^C0G54	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:10778:3437	89	TPP	8	44	109M1I21M5S	=	8	0	GGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCACAATCGGGCTTCGGTCCGGTTCTGTGG	GE?:EEGGGGD>>GGGGECC?C:A;GGGGGGGHGGGGGGGGGE8@BGGGGGGGEEEDGGGGGGHGGGGGGGGHHHHHHHHHHHHHHHHHHHHHEHHHHHHHHIHGHHFF@HHHHGGGGGGDDDDDDEDBBBAAAAA	AS:i:254	XN:i:0	XM:i:0	XO:i:1	XG:i:1	NM:i:1	MD:Z:130	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:13477:3450	89	TPP	1	42	13S25M3D5M1D103M5S	=	1	0	TTCCGATCATCGGGGCCTTCGGGCCAAGGACTCGGGGTTTTCTCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTCATCCC	CCC)GCEC>GAAGGE8'8<>GGEGEEC:*8>DGGGC?EEEG>GEGGGGGGGGGHGGGGEGEEC;GGGGGGGEBGGGGGGGGGGGGGGHHHHHHHHHHHHHHHHHHHHHHHEHHHHHHIIIIIIIIHHHHGGGGGGDDDDEEEDBB?A?AAA	AS:i:248	XN:i:0	XM:i:1	XO:i:2	XG:i:4	NM:i:5	MD:Z:25^GCC0C4^G103	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:18624:3452	89	TPP	1	42	9S96M1I15M1D25M5S	=	1	0	GATCGAAACGGCCTTCGGGCCAAGGTCTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAGATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCTAATCCGGTTCGCCGGTCCAAATCGGGCTTCGGTCCGGTTCATGGC	*CECE?).DEC>D>A<GGEEEEEEGGGGD;ECCEC8EGGGGGGGGGGGGGGCGGGGECEEGGEEGGGGEEGGGGEGGDDEGGGGGGGHFHHHHHHHHHHHHHHHHHDHHHHHHHHHHHHIHIIIHEBHHGGGGGGDDDDDDDDBBBAA???	AS:i:248	XN:i:0	XM:i:3	XO:i:2	XG:i:2	NM:i:5	MD:Z:16A31A47G14^A25	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:19396:3449	89	TPP	1	42	11S110M2D25M5S	=	1	0	CCGATCAGTGCGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGCAAGGAAGTTCTCAATCCGGTTCGCCGTCCAAATCGGGCTTCGGTCCGGTTCTAATA	8GECEEDDAGGCED?GDGGGECGGGC:GAGGGGGCCGGGEEGGGGGGGEGGGGGGGGEGCC@:GEECDDDEEGEDEGHGGGGGGGGGGEECHEFFHHHHHHHHHHHHHECHHHEHHHEHHHIIIIHHHHGGGGGGEEDDEEDDBBB?AAAA	AS:i:251	XN:i:0	XM:i:3	XO:i:1	XG:i:2	NM:i:5	MD:Z:83T1G10G13^GA25	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:16885:3462	89	TPP	1	44	9S137M5S	=	1	0	GATCGACATGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATAACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTCTGTTC	88EEA:?C??F>>D;2CCA?EAEEEAD>DDEEAEFEFFFD?DC?EFFFFFFFFFFFEFE>BEEECEEEFFEFFFFFEFFFFFFFFEFHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHIIIIIIHIHHHHFFEFFFDDDDDDDDBBB?????	AS:i:270	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:59C77	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:18360:3463	89	TPP	4	44	12S134M5S	=	4	0	CGATCACAGCGGCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTCCACGA	GEC??0'DGAEDGAA>EGEGGEE?C8DGGDC8GEGGGEEGGGGGGGGGGGGGEGGGGGEBEGGECGGGCDGEDEGGGGEDEEGGGGGGHHHHHHHHHHHHHHHHHHHHHHHHHHHHHIIIIIIIIHHHHGGGGGGDDDDEDDDB?A?????	AS:i:268	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:134	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:14481:3463	89	TPP	1	42	18S93M7D5M2D30M5S	=	1	0	TGCTCTTCCGATCGTAAGGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCGGTTCCAGATCCAAATCGGGCTTCGGTCCGGTTCAGCGT	*CCEECCCGGGGGGGGGGCEA?DAGEEGEGGGGGCDBGGGGGGGEAGGEBGGGGGGGGGGFGGGGGGFDHHHHHHFFFCHHHHHHHHHHIIIIIHHHHHEIIIIIIIIHFCHIIHIIIIIIIIIIHHHHGGGGGGDEEDDDEEB?AAAAAA	AS:i:232	XN:i:0	XM:i:1	XO:i:2	XG:i:9	NM:i:10	MD:Z:93^CTCGATC5^CG2G27	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:15537:3468	89	TPP	93	22	68S45M5S	=	93	0	ATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCAGGTGGGCCTTCGGGCCAAGGACTCTCTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTCCCTTC	EEEEGGGGEEGGEGEECGGEEECCACGGEEACCEECCEEGECGEEE@D8GGGGGGHHHHHHFFHFEHHHHHHHHHHHHHEHHFHHIIIIHIIHHHHGGGGGGDDDDDDDDB???????	AS:i:90	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:45	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:9129:3475	89	TPP	83	44	55M5S	=	83	0	GTAAGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTCAGCAG	HHHHHHHHHHFFEEADHHHHHHHCHHIIIIIIIIHHHHFFFFFDDDDDDDEDBBB?????	AS:i:106	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:3G51	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:18668:3479	89	TPP	1	44	7S55M2I82M5S	=	1	0	TCAGAAAGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCAAGGATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTCTTGTC	1C?CC?C882>><<8:C:GEGEC8GD88GCEEGCCGGGGGGGGEGEGGGGGGGGGGEC?GEEGGGGGEDEEDGGGGGGGGEGGGGGGGHHEHHHHHHHHHHHHHHHHHHHHHHHHHHHIIIIIIIHBHHGGGGGGDEDDDDEEBBBAAAAA	AS:i:263	XN:i:0	XM:i:1	XO:i:1	XG:i:2	NM:i:3	MD:Z:56T80	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:13672:3463	89	TPP	1	44	9S137M5S	=	1	0	GATCAAGAAGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTCGTCCT	?GGEC?EEGGE2>GD<?:::GEGEGGADDGGEGGGGECG;2DC?GGGEGGGGGGGGGGGE@GGEGGGEGGGGGGGGGGGGEEGGEGBDFEDGGHHHHHHHHHHHHHHHHDHHHHHHHIIIIIIIIHHHHGGFGGGDDDDDDEEAABAAA??	AS:i:274	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:137	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:19408:3485	89	TPP	1	44	9S137M5S	=	1	0	GATCATACTGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTCGAACT	:EC?C?:??EE8<>GGGGGGGGGGGGGGGGGGGGGGCGEADGGGGGGCGGECEGGGEGEEBGGCCGGEGGGGGGGGGGGGGDEGGGEGGHEFFFDHHHHHHHHHHHHHHHHHHHHHHIIIIIHIHHHHHGGGGGGDEEDDEEDBAA?????	AS:i:274	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:137	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:14136:3486	89	TPP	1	44	9S137M5S	=	1	0	CGACGGGTAGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGATTTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTCGATTC	DEC???:0:ECAGGGGGEGGGECEGEDG?GGGGECGGEGGGGGGGGGGGGGGGGGGEEAE@GEGEGGEECGGGGGGGGGGGGGGGGGGGHHHHHHHFC,HHHHHHHHHHHHHHHHHHIIIIIIIIHHHHGGGGGGDDDDDDDEBBBAAAAA	AS:i:267	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:89A0G46	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:12786:3488	89	TPP	1	44	9S137M5S	=	1	0	GATCATCACGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGTGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTCCTACT	0CEEGG>>DE8<><'>EECGCC?::CDGG>CE?CCCCGGGEECEEEEGGCGECGGGGCC:)GGGGEAEEDGGGGGGEGEGGGGGGGGHHCHHHHHHHHHHHHFHEHHHHEEFEHHHHIIHGHIHIHHHHGGGGGGDDDDDDDDBBBA????	AS:i:270	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:34C102	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:17814:3490	89	TPP	117	22	9S21M5S	=	117	0	GGGCCAAGGAATCGGGCTTCGGTCCGGTTCCTTAT	AAA/CA>>>HHHHFFEEEFEEEEEDEEBBB?????	AS:i:42	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:21	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:9656:3491	89	TPP	10	44	20M1I108M5S	=	10	0	GCCAAGGACTCGGGGTGCCCTTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTCTTTCC	GGGGGGEEGGGDGEECC??GEGGGGECGGEEEGGGGGGGGGGGGGGGGGGGGGGHGGHGGGGGGGGGGGEGHHHHHHHHHHHHGHHHHHGHHHHHHHHHHIIIIIIIIHHHHGGGGGGDEDDEEDEBBBAA???	AS:i:250	XN:i:0	XM:i:0	XO:i:1	XG:i:1	NM:i:1	MD:Z:128	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:21374:3492	89	TPP	1	44	9S137M5S	=	1	0	GATCGATACGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGTATAATGCCAGCGTAAGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTCTCATG	CEEEGCC??GEADGGGGC:*EECGGGGGGDGGGGGGEGEDDGECGGGGGGGGGGGGGGGGEGGGGEGGGGGGGGEEEGFGGGGGFGGHFHHHHHHHHHHHHHHHHHHHHHHHHHHHHIIIIIIIIHHHHGGGGGGDEDDDDDEBBBAAAAA	AS:i:266	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:70G14G51	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:21525:3465	89	TPP	1	44	11S35M2D100M5S	=	1	0	CCGATCGCCGGGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGAAGGCTGAGAAATACCCGCATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTCTCACG	.<8.?DGDDDEEE?8?GGGGGGGGGGEGGGGGGEEEGGCGGGEEGEGGGGGGGGGEEGC:)GGGGGGGGGGGGGGGGGGGGGHGGGGHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHIIIIIIIIHHHHGGGGGGEEEEEEEEBBAA??A?	AS:i:259	XN:i:0	XM:i:1	XO:i:1	XG:i:2	NM:i:3	MD:Z:35^GT19T80	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:21136:3498	89	TPP	1	24	41S79M32D26M5S	=	1	0	GATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCGGGACGGCCTTCGGGCCAAGGACTCGGGGTGCCCTCTTCTGTGAAGGCCGAGAGATACCCGTATCACCTGATCTGGATAATGCCCTACACATCGGGCTTCGGTCCGGTTCATCCC	0EGGEECEGGGGGGGGGEEGEGEGGECEEGGGG?GGGGGGGGGGGEBEEGGGGGGGEDGGGGHHHHHHFFHFHHHHHHHHHHIIHHIIIIIHHHHHHHIHHHIIIHIIHGHIIIHIIIIIIHFFAHHHHGGGGGGEEEDDDDDBBAAAAAA	AS:i:135	XN:i:0	XM:i:9	XO:i:1	XG:i:32	NM:i:41	MD:Z:30T0C1G0C8T4A30^AGCGTAGGGAAGTTCTCGATCCGGTTCGCCGG0A1C2A20	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:20482:3501	89	TPP	1	44	9S137M5S	=	1	0	GATCAGGAAGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAGATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTCGACGT	*CEEE:CGEGE8D><8EEC???CGE?2GD;C?8GGGEEE>DDECC:GEEGGGGGGGCGEBBGGEECGGGGECGGGGGGGGEGGGGEBGGGDHHHFFHHHHHHHHHHHHHHHHHHHHHHHIIHIIIHHHHGGGGGGEDDDDDDD?BB??<??	AS:i:270	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:48A88	YT:Z:UP
M01228:25:000000000-A1CW0:1:1101:21880:3502	89	TPP	1	44	9S137M5S	=	1	0	GATCACAATGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTCTGTGAAGGCTGAGAAAAACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTCCACTC	:ECCCEECEC?2DD<DGGGEGGGGGGDDGGGGCEC8?EEGEEGGGGGGGGGCGGGGGGGGGGGHGGGGGGGGGGGGGGGGGGGGGGGEHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHIIIIIIIHHHHGGGGGGDDDDDDDDBB??????	AS:i:262	XN:i:0	XM:i:3	XO:i:0	XG:i:0	NM:i:3	MD:Z:33G0C15T86	YT:Z:UP
)";
	std::string test_dir = getTestFileDir().string();
	std::ofstream dummy(test_dir+"/tmp/tmp_vrcc.sam");
	dummy << sample_reads;
	dummy.close();
	EXPECT_NO_THROW(
		parseSAM(test_dir+"/tmp/tmp_vrcc.sam",
				 test_dir+"/tmp/tmp_vrcc.mut",
				 "", // debug_outname
				 "", // primers_filename
				 800, // max_paired_fragment_length
				 10, // min_mapq
				 false, // right_align_ambig_dels,
				 false, //right_align_ambig_ins,
				 7, // max_internal_match,
				 30, //min_qual,
				 9, // exclude_3prime,
				 "", //mutation_type,
				 false, // variant_mode,
				 false, // trim_primers
				 false, // require_forward_primer_mapped
				 false, // require_reverse_primer_mapped
				 10, // max_primer_offset
				 false, // input_is_unpaired
				 false, // debug,
				 true //warn_on_no_mapped
		);
	);
}


TEST(Debug, R2Crash) {
	std::string sample_reads = R"(M01228:25:000000000-A1CW0:1:1101:12330:3450	83	TPP	101	44	16M1I21M5S	=	1	-141	CGGTTCGCCGGATCCAGAATCGGGCTTCGGTCCGGTTCCCGCC	CC=ECC>7EC5--EA5->*CEFF>CC+BBBBB@BB@@@?????	AS:i:68	XN:i:0	XM:i:0	XO:i:1	XG:i:1	NM:i:1	MD:Z:37	YS:i:204	YT:Z:CP
M01228:25:000000000-A1CW0:1:1101:12330:3450	163	TPP	1	44	4S58M1D49M	=	101	141	AGCTGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGCGAAGGCTGAGAAATACCCGTACACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGC	????ABABDDDDDDDDGGFGGFHHIHHHCHHHHIIHIIGHHDHHHCHIIIHIIIIIIFHHHHHEHH=CDFFHHHHFFHHHHHHHHGEG5DDBEDB=.D=DCAA*;8BECEE	AS:i:204	XN:i:0	XM:i:1	XO:i:1	XG:i:1	NM:i:2	MD:Z:36T21^T49	YS:i:68	YT:Z:CP
)";
	std::string test_dir = getTestFileDir().string();
	std::ofstream dummy(test_dir+"/tmp/tmp_r2c.sam");
	dummy << sample_reads;
	dummy.close();
	EXPECT_NO_THROW(
		parseSAM(test_dir+"/tmp/tmp_r2c.sam",
				 test_dir+"/tmp/tmp_r2c.mut",
				 "", // debug_outname
				 "", // primers_filename
				 800, // max_paired_fragment_length
				 10, // min_mapq
				 false, // right_align_ambig_dels,
				 false, //right_align_ambig_ins,
				 7, // max_internal_match,
				 30, //min_qual,
				 9, // exclude_3prime,
				 "", //mutation_type,
				 false, // variant_mode,
				 false, // trim_primers
				 false, // require_forward_primer_mapped
				 false, // require_reverse_primer_mapped
				 10, // max_primer_offset
				 false, // input_is_unpaired
				 false, // debug,
				 true //warn_on_no_mapped
		);
	);
}


TEST(Debug, OutOfStepPairSegfault) {
	// this failed previously because some read pairs map to different RNAs

	std::string sample_reads = R"(M01228:25:000000000-A1CW0:1:1105:14571:4577	161	16S	1262	44	32M	23S	1426	0	CCTCGCGAGAGCAAGCGGACCTCATAAAGTGC	?????@@<BB+ABBBBEBACCFHF;CC8C898	AS:i:64	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:32	YS:i:294	YT:Z:DP
M01228:25:000000000-A1CW0:1:1105:9657:4651	83	16S	1114	2	51M1D3M79S	=	1114	-55	CTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGGATGTTTTCAGAGATGAGAATGTGCCCGGGAACTCAAAGGAGACTGCCAGAAGTCTCGCAACGACGCAACCCTTATCCTTTG	EGGEGEGGGGGGEGGGEGGGGGGGGGGGGGHHHHHHHHHHHHHHHHHHHHHHHHHIHGHFHHGGFDHHIIIHEHHHHHIIIIHHGIIIIIIIIIHHHHIIIIIHECAHHIHGGGGGGDDDDDDDDBBBAAA??	AS:i:102	XS:i:67	XN:i:0	XM:i:0	XO:i:1	XG:i:1	NM:i:1	MD:Z:51^T3	YS:i:102	YT:Z:CP
M01228:25:000000000-A1CW0:1:1105:9657:4651	163	16S	1114	2	51M1D3M79S	=	1114	-55	CTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGGATGTTTTCAGAGATGAGAATGTGCCCGGGAACTCAAAGGAGACTGCCAGAAGTCTCGCAACGACGCAACCCTTATCCTTTG	?AA??BBBDDDDDDDDFFFGGGHHHHHEHHHHHHIGDHHHGFFHFHHIHHIHFFDGGHHBFBCDFDHECDGGF?FHHHHHHHHHHHFFHCDFAFFFGGGDEGGEDE.>>=BGEBEGGBGGGGGGGGGG?EGGG	AS:i:102	XS:i:71	XN:i:0	XM:i:0	XO:i:1	XG:i:1	NM:i:1	MD:Z:51^T3	YS:i:102	YT:Z:CP
M01228:25:000000000-A1CW0:1:1105:16437:4693	83	16S	203	44	146M5S	=	18	-331	GGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGTAGAT	GC?>GGGGGGGCEGEGGEGGGGGGGEEGGGGGGGGGGGGGGFGGGGGGGHHHHHHHHHHHHHHHHGFCFHEHHHHHHHHHHIIHGIIIIIIIHHHGFGFFIIIHHHHHHHHHIIIIIIIIIHIHHFHHIGGGGGGDDDDDDDDBBBAAAA?	AS:i:292	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:146	YS:i:283	YT:Z:CP
M01228:25:000000000-A1CW0:1:1105:16437:4693	163	16S	18	44	63M2D84M	=	203	331	CATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGATTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGT	?????B?BBDDDDDDDFFEEEDHHHHHHHHHHHHHHHHHFHHHHHHHHHHHHDFHDFFHHHHFGGFFHHHHHHHDGHFHFE@CDFHHEEEEEEE1@;C,3?BBCBEEEEEEEEEEEECE?CEE:CADDD8?AECAAC:AEA?A*8A8	AS:i:283	XN:i:0	XM:i:1	XO:i:1	XG:i:2	NM:i:3	MD:Z:63^AG0C83	YS:i:292	YT:Z:CP
)";
	std::string test_dir = getTestFileDir().string();
	std::ofstream dummy(test_dir+"/tmp/tmp_oosps.sam");
	dummy << sample_reads;
	dummy.close();
	EXPECT_NO_THROW(
		parseSAM(test_dir+"/tmp/tmp_oosps.sam",
				 test_dir+"/tmp/tmp_oosps.mut",
				 "", // debug_outname
				 "", // primers_filename
				 1500, // max_paired_fragment_length
				 10, // min_mapq
				 false, // right_align_ambig_dels,
				 false, //right_align_ambig_ins,
				 7, // max_internal_match,
				 30, //min_qual,
				 9, // exclude_3prime,
				 "", //mutation_type,
				 false, // variant_mode,
				 false, // trim_primers
				 false, // require_forward_primer_mapped
				 false, // require_reverse_primer_mapped
				 10, // max_primer_offset
				 false, // input_is_unpaired
				 false, // debug,
				 true //warn_on_no_mapped
		);
	);
}


TEST(Debug, RibosomeSegfault) {
	std::string sample_reads = R"(M01228:25:000000000-A1CW0:1:1101:10980:5999	163	16S	75	3	39M1I65M	=	183	255	GGAAACAGCTTGCTGTTTCGCTGACGAGTGGCGGACGGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGC	?????BBBDDD<BBDDFFFEEEHHHHFEHHHHHH+CCCHH*>5C5CFHDGDBDFHEHFFGFFHFF??CF,4CEHF8=DB,4,=BDDBBDDEE8;?B,=C==CE;<	MD:Z:4G0A8T0C2T85
M01228:25:000000000-A1CW0:1:1101:10980:5999	83	16S	183	3	147M4S	=	75	-255	CGTCGCAAGACCAAAGAGGGGGACCCTCGGGCCTCTTGCCATCGGATTTGCCCAGATGGGATTAGCTTGTTGGTGGGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACAAGAG	<GGGGEGECCEGGGGGGGGGGGGGGGGEGGDEGEDEEBDEGEEHHFHHHHHHHHHFDHHHHHGGDBHHE?HED>IIHHDHFCHGFIHHGEHHHC7IHHHHFHFHHFHHFHHFCIIIIHIIHFHF>HHHHGGGGGGDDDDDDDDBBB?????	MD:Z:25T21G19A2A19T56
)";
	std::string test_dir = getTestFileDir().string();
	std::ofstream dummy(test_dir+"/tmp/tmp_rs.sam");
	dummy << sample_reads;
	dummy.close();
	EXPECT_NO_THROW(
		parseSAM(test_dir+"/tmp/tmp_rs.sam",
				 test_dir+"/tmp/tmp_rs.mut",
				 "", // debug_outname
				 "", // primers_filename
				 800, // max_paired_fragment_length
				 10, // min_mapq
				 false, // right_align_ambig_dels,
				 false, //right_align_ambig_ins,
				 5, // max_internal_match,
				 30, //min_qual,
				 10, // exclude_3prime,
				 "", //mutation_type,
				 false, // variant_mode,
				 false, // trim_primers
				 false, // require_forward_primer_mapped
				 false, // require_reverse_primer_mapped
				 10, // max_primer_offset
				 false, // input_is_unpaired
				 false, // debug,
				 true //warn_on_no_mapped
		);
	);
}


TEST(Debug, Segfault2) {
	std::string sample_reads = R"(M01228:25:000000000-A1CW0:1:1101:17512:6032	99	16S	1323	255	36M1D41M	=	1469	167	GACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACC	?????BBBB?BBBBBB>ACC>CF>CEHFFHHH@EEFDGCGHHHHHHEAACFFBCF;EFHHHDCCCCFHDFHBACCEH	MD:Z:36^C41
	M01228:25:000000000-A1CW0:1:1101:17512:6032	147	16S	1469	255	21M10S	=	1323	-167	CTTTGTGATTCATGACTGGGGGTGAAGCGAC	=EBE8.EEECA8.>E@@@@@--@=====<<,	MD:Z:21
)";
	std::string test_dir = getTestFileDir().string();
	std::ofstream dummy(test_dir+"/tmp/tmp_s2.sam");
	dummy << sample_reads;
	dummy.close();
	EXPECT_NO_THROW(
		parseSAM(test_dir+"/tmp/tmp_s2.sam",
				 test_dir+"/tmp/tmp_s2.mut",
				 "", // debug_outname
				 "", // primers_filename
				 800, // max_paired_fragment_length
				 10, // min_mapq
				 false, // right_align_ambig_dels,
				 false, //right_align_ambig_ins,
				 5, // max_internal_match,
				 30, //min_qual,
				 10, // exclude_3prime,
				 "", //mutation_type,
				 false, // variant_mode,
				 false, // trim_primers
				 false, // require_forward_primer_mapped
				 false, // require_reverse_primer_mapped
				 10, // max_primer_offset
				 false, // input_is_unpaired
				 false, // debug,
				 true //warn_on_no_mapped
		);
	);
};


TEST(Debug, MemoryCorruption) {
	std::string sample_reads = R"(M01228:25:000000000-A1CW0:1:2101:14804:1374	163	23S	698	1	1S129M1I12M4S	=	697	141	NCAGGTTGAAGGTTGGGTAACACTAACTGGAGGACCGAACCGACTAATTTTGAAAAATTAGCGGATGACTTGTGGCTGGGGGGGAAAGGCCAATCAAACCGGGAGATAGCTGGTTCTCCCCGAAAGCTATTTTAGGTAGCGCCCATG	!5<???BBDDDDDDDDDCFFFFCFHFFFHHFHHHHHHHHHHHHHHHGHHHHHHHHHFFHHHFHH@C@GGHHG,CDF,CEEHH'44??CEEEEAEEEEEAEDDD22AEECEACA:?:AEEAA;D?DEACEAEEEEEECEEEDD?DEEE	MD:Z:47G33T59
M01228:25:000000000-A1CW0:1:2101:14804:1374	83	23S	697	1	129M1I13M4S	=	698	-141	GCAGGTTGAAGGTTGGGTAACACTAACTGGAGGACCGAACCGACTAATTTTGAAAAATTAGCGGATGACTTGTGGCTGGGGGTGAAAGGCCAATCAAACCGGGAGATAGCTGGTTCTCCCCGAAAGCTATTTTAGGTAGCGCCCATN	EEGGEEGGGGGEEGEGEEGGGGGGGGGGGGGGGGGGGC@GGGGGGGGGGGGGGGGGGGHHHHHHHHHHHHHHFHHHHHHIIIIHIIHIIIIIIIHHHHEHIIHGHHIIHIIIIIHHHHHHHIIIIGGGGGGDDDDDDDDBB???<5!	MD:Z:48G93
)";
	std::string test_dir = getTestFileDir().string();
	std::ofstream dummy(test_dir+"/tmp/tmp_mc.sam");
	dummy << sample_reads;
	dummy.close();
	EXPECT_NO_THROW(
		parseSAM(test_dir+"/tmp/tmp_mc.sam",
				 test_dir+"/tmp/tmp_mc.mut",
				 "", // debug_outname
				 "", // primers_filename
				 800, // max_paired_fragment_length
				 10, // min_mapq
				 false, // right_align_ambig_dels,
				 false, //right_align_ambig_ins,
				 5, // max_internal_match,
				 30, //min_qual,
				 10, // exclude_3prime,
				 "", //mutation_type,
				 false, // variant_mode,
				 false, // trim_primers
				 false, // require_forward_primer_mapped
				 false, // require_reverse_primer_mapped
				 10, // max_primer_offset
				 false, // input_is_unpaired
				 false, // debug,
				 true //warn_on_no_mapped
		);
	);
}







// test mergeMatePairs() on several examples from paired_mapping/difficult_read_ids.txt
TEST(OverlapResolution, Examples) {
	std::string sample_reads = R"(M01228:25:000000000-A1CW0:1:1101:17371:10707	83	16S	18	14	54S95M2S	=	16	203	GTGAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATACCCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAATGGTAACAGGAAACAGCTTGCTGTTTCGCTGACGAGTGGCGGACGGCA	:*CEGGEGGGGECEEGGEGE??GGGEAE?3;<GECA>>>>=@D;1??@??6,+HHHHHHHHHHHHFFHHHHHHHHHCHFIIHGEFHHIHGHHIIIHHHIIIIIIIIIIIIIIIIIIIIIIIIIIIHHHHGGFGGGDDDEEDDDAAAAAA?A	AS:i:163	XS:i:100	XN:i:0	XM:i:6	XO:i:0	XG:i:0	NM:i:6	MD:Z:49C11G0A8T0C2T19	YS:i:167	YT:Z:CP
M01228:25:000000000-A1CW0:1:1101:17371:10707	163	16S	16	14	96M	=	18	-203	ATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAATGGTAACAGGAAACAGCTTGCTGTTTCGCTGACGAGTGGCGGACG	AAAAABBBDEEEEDDDGGGGGGHHHHHHHHHHIIIIHIIIIIIIIHIIHIIIIHGHFDGIGIHHHIIIIHFHHGHHHHDHHHHHHHHHHHGGGB@E	AS:i:167	XN:i:0	XM:i:6	XO:i:0	XG:i:0	NM:i:6	MD:Z:51C11G0A8T0C2T18	YS:i:163	YT:Z:CP
M01228:25:000000000-A1CW0:1:1101:6728:11695	99	16S	1106	44	138M5D7M6S	=	1107	158	GCGCAACCCTTATCTTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGCCCCTTACGACCAGGGCTACACACGTGCTACAATGGCCAACGAGCGCAAC	?????@@BDDDDDDDDFFFFFFHHIIIHHHHHCEHEHHC>>EHBGCFGHFD+AFCFHFGFGFGHIFHHIIBDEHHHH4DHHDFEEFFEFEFEFEFFFEFEEEFEFFEEEFFEEEEEEEFAEFFFEDEFCCEFF?CEEEFFC?ECA?D>D>D	AS:i:268	XN:i:0	XM:i:3	XO:i:1	XG:i:5	NM:i:8	MD:Z:14C86G36^GCATA3A3	YS:i:266	YT:Z:CP
M01228:25:000000000-A1CW0:1:1101:6728:11695	147	16S	1107	44	137M5D7M1D7M	=	1106	-158	CGCAACCCTTATCTTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGCCCCTTACGACCAGGGCTACACACGTGCTACAATGGCCAACGAGCGCAACC	.*>AAAA:E?CCCAEC88:FED>>;8@E6EEEEBEFFFFFFFFFFFHHHHHHHHFHHHHHHHHHHDBFGCHHHIHHHIIHHIHHFEHHHIHIHHHHHFHFHIIHFE7CCAIIHHHGFHFHFHHEC>FFHFFFFFFDDDDDDDDBBB?????	AS:i:266	XN:i:0	XM:i:5	XO:i:2	XG:i:6	NM:i:11	MD:Z:13C86G36^GCATA3A3^A0A2G3	YS:i:268	YT:Z:CP
M01228:25:000000000-A1CW0:1:1101:11063:14857	99	16S	1226	42	96M3D13M42S	=	1305	231	CACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCTTCATGGAGTTGAGTTGCAGACTCCAAGCGGACCTCATAAAGTGCGTCGTAGTC	??AAAAB?DDDDDDDDGGGGGGHIIIIIIIIIIIHHHHHIIHHHHHEHIIHHHHHHHIHHHHHHGCGHHDEEEHHEHHEEDEGGGGGEGGGGGGGGGGGCEGGGGGGGGGGGGGGGGGGGGEGGGGGGGGDGGGGGGGGGGGGGA??EECE	AS:i:198	XN:i:0	XM:i:3	XO:i:1	XG:i:3	NM:i:6	MD:Z:96^CGA2C4A3C1	YS:i:292	YT:Z:CP
M01228:25:000000000-A1CW0:1:1101:11063:14857	147	16S	1305	42	145M1D6M	=	1226	-231	GATTGGAGTCTGCAACTCAACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTCGGGA	CC:GEECGGEGGGEECGGGGECEGGGEGGDEGGGGGDGGGGGGGGGGGGGGGGEGGGGGGGGGGGGGGGFGFFGHHHHHHHHHFHHHFFEEHHHHHHHHHGHHHHGHIIIIIIIHHHIIIIIHIIIIIIGGGGGGDDDDDDDDAAAAA???	AS:i:292	XN:i:0	XM:i:1	XO:i:1	XG:i:1	NM:i:2	MD:Z:18G126^T6	YS:i:198	YT:Z:CP
M01228:25:000000000-A1CW0:1:1102:20323:13872	99	16S	42	44	4S147M	=	182	300	CGATGCAGGCCTAACACATGCAAGTCGAATGGTAACAGGAAACAGCTTGCTGTTTCGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGC	A?AAAABBEEDDDDDDGGGGGGIIHIHHIIIHGHHIIIIIIHFHHIIIIHIHGHIIHHHHHIHFHEFHHHHHDCHH3<CD,CFFFDF,?FFFFGGGGGGGGDDEG-=CGBE:8CECCECECEGGGCEGEGCEGEG:?CECDG?EC??:?28	AS:i:268	XN:i:0XM:i:6	XO:i:0	XG:i:0	NM:i:6	MD:Z:25C11G0A8T0C2T95	YS:i:260	YT:Z:CP
M01228:25:000000000-A1CW0:1:1102:20323:13872	147	16S	182	44	30M1I50M11D65M5S	=	42	-300	ACGTCGCAAGACCAAAGAGGGGGACCTTCGATGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGATTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGCTGAG	GGGA>ECGGECGGGGGG@GGGECCGGGGEGEEDFFB@DGFDEHE@FHHHEFDFFGGFGFHHHGHHHFFIHHFFHHGFHHFGEEHC>FCHHHHIIHGFHFFHHHHGFFGFHHHHFDHFFHFFEFAHHIIIGGGGGGBDDEDDDDBBBAAAAA	AS:i:260	XN:i:0	XM:i:2	XO:i:2	XG:i:12	NM:i:14	MD:Z:30G47G1^AACGGCTCACC65	YS:i:268	YT:Z:CP
M01228:25:000000000-A1CW0:1:1103:13613:5633	99	16S	1023	44	151M	=	1187	375	TGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACT	AAA?ABBBEEEEDDDDGGFGGGIIIHIGHIIIHIIIHHFHIIIIIHHIHHFDFHFHHHHHHDFGHHHHHIGHEGGFHHHHHHHHHHHHHGGGGGGGGEGGGEGGGGGEGEEEGGGGEEDGGGGGEECEGEGCECCC?EGGEECEC?CEGC:	AS:i:302	XN:i:0XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:151	YS:i:237	YT:Z:CP
M01228:25:000000000-A1CW0:1:1103:13613:5633	147	16S	1187	44	48M60D103M	=	1023	-375	GATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACAC	GCECCEEEEGEGGECGGEEGGGEECC@EGGGGGGGGGGEEGEGGGEEDEEGGGBGHHHHHHHHHHHHHHHHHHHHHIIHHIIHHHHHHHFIIIIHHHIIIIHHFHHIIIIHHHHGFHHHHHHHIIIHHHCGGGGGDDDDDDDDBBB?????	AS:i:237	XN:i:0	XM:i:0	XO:i:1	XG:i:60	NM:i:60	MD:Z:48^TACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCG103	YS:i:302	YT:Z:CP
M01228:25:000000000-A1CW0:1:1103:3438:14114	83	16S	438	44	150M	=	277	-311	TTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACACAGGCG	FED>8DEFFEEEECEFFFEAEEFEEEEEA:1AC?AACECC?EEFEEEE@>EEEFFEEEEFFFEEDE@?FFFHDHCEHGF<FGGHHCCDHHHHHHIHFFEHIIHIHGHHHHEHHFECEDHIHHHHFEFHDC9EFFDDDDDDDDB<??????	AS:i:296	XN:i:0XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:143G6	YS:i:266	YT:Z:CP
M01228:25:000000000-A1CW0:1:1103:3438:14114	163	16S	277	44	31M3I65M1D52M	=	438	311	CGACGATCCCTAGCTGGTCTGAGAGGGTAACGGTAGGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGG	?????@<?BBB?BDDDFFFFFFHHHHEEHHHHHHHHHHHFFFHHHFHHHHHHHFDFDDCCEH=AEFGHHFGFHHHDEH:CDHFF7@CFHHEC?DFFFFFFFF=?4>BBEEECEEEEEEB?CE?=BCBD;>D>8?ACA::AAEEEACE:AED	AS:i:266	XN:i:0	XM:i:4	XO:i:2	XG:i:4	NM:i:8	MD:Z:26A1G2C0A63^A52	YS:i:296	YT:Z:CP
M01228:25:000000000-A1CW0:1:1103:8088:26355	99	16S	830	44	84M6S	=	830	-104	GAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCACCTTCC	?????BB?BDDDDDDBFFFDCCHHHFCFEHEH@FGFFH7>CHHHHC=EFH@CCEEFHHD+CDDDDHHEEEEEEE=DDEEEBB?CC;?;;=	AS:i:168	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:84	YS:i:188	YT:Z:CP
M01228:25:000000000-A1CW0:1:1103:8088:26355	147	16S	830	44	104M	=	830	-104	GAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGGGAAGACTCAAATGAATTGACGGGGCCCAG	*0EAEEEEEEEB?B:,E?BB;;DEED@:?D4,@+C:C=<DFC+5)ED<*7HDDDD@FAFHHECCEHHC@HFFA?8A888HHFBFBFFFDDA5++BD@???????	AS:i:188	XN:i:0	XM:i:5	XO:i:0	XG:i:0	NM:i:5	MD:Z:74T0T2A20G2C1	YS:i:168	YT:Z:CP
M01228:25:000000000-A1CW0:1:1104:6129:5751	97	16S	1411	28	86M2D15M1D25M	=	1415	129	CCATTTGAGTGGTTTTCAAAAGAATTATTTATCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCAACAAGGTAACCGTAGGGAACCTGCGGTTGGATCACCTCC	???A?BBBDBBD7AB<CGGGGGIHECFEEFFAEHHIIIIIIIIHHHHHHHHHHHHIIHIIIIHCFHHIHHHHIIHIIIH5CFGHGFHHHHHHHHHHHHHHHHHGGGGGGGGGGGGG@CCEGGGGGG	AS:i:207	XN:i:0	XM:i:8	XO:i:2	XG:i:3	NM:i:11	MD:Z:4G0G6G2G8G2G0G2G54^GT15^G25	YS:i:53	YT:Z:DP
M01228:25:000000000-A1CW0:1:1104:6129:5751	145	16S	1415	28	80M2D17M1D20M5S	=	1411	-129	GGGAGGGGGGGGAAAAAGAAGGAGGGAGAGGAAAAGGAGGGAGGGAGAGGAAAAAGGGGGGAGGAAGGAAGGGGGGGAAGGAAAAAAGGGAAAAGGAGGGAAAAGGAGGTTGGATAAAATAC	EEA2DEEEFEED=FFEEDBEEFDFFD?@.FD@C@7DD=HED?DF@,@,FD@CCE<HHHHFGEEC-EFA5CCHHHHHGGHGFCHHGHHHHGGGFFHHHFHHFA8CFFBDDDDDDDBBB?55,<	AS:i:53	XN:i:0	XM:i:43	XO:i:2	XG:i:3	NM:i:46	MD:Z:5T3T0T1C8T3T2C0T0T2C0C0T0T0C7C1C0T0T1C0C1C0T0T0T1T2T0T0C1T2C0T4T4^TC1T2C4T2C0C1T1^G5C0C0T1C8C1	YS:i:207	YT:Z:DP
M01228:25:000000000-A1CW0:1:1104:11684:10093	97	16S	521	22	24M1I4M5D6M	=	521	-158	GCAGCCGCGGTAATACGGAGGGTGTCGAGATCGGA	AAAAAAAAEDDDDDEDFFFGGGCFHHHIHHIIIHH	AS:i:48	XN:i:0	XM:i:1	XO:i:2	XG:i:6	NM:i:7	MD:Z:25A2^CGTTA6	YS:i:48	YT:Z:DP
M01228:25:000000000-A1CW0:1:1104:11684:10093	145	16S	521	22	67S24M3S	=	521	158	TTTTTTTTTAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTGCAGCCGCGGTAATACGGAGGGTGTCG	GGGGGGGGEEGFFF@8CCCEFDHECHHGFHFFCFDCCFFFIHHFHEFAEHHHHHIIHHIHHIIIIIIIHHHHGGGGGGDDDDDDDDA?@??A??	AS:i:48	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:24	YS:i:48	YT:Z:DP
M01228:25:000000000-A1CW0:1:1104:2330:11289	153	16S	549	44	142M	=	549	0	CGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGGTAAGCG	AA?EAA8ACEC:*?CCA??BA??EEEB>*BBBEEDEBEFFFFHF=HFFF?FFBF?HHHHFC7HHHEHHHGFHHHDFDHGFHHHG@HECGGFEBCA-G;F=GF?FFCC,EE7FBHHHHEBFFF;FFFDBDBDDBD?B??????	AS:i:272	XN:i:0	XM:i:3	XO:i:0XG:i:0	NM:i:3	MD:Z:135T0G0T4	YT:Z:UP
M01228:25:000000000-A1CW0:1:1104:4137:13019	83	16S	880	44	48M2S	=	785	-144	CGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGCA	>EE>@EEF=@EE>>CCC>D@8.>EA>CEC8A/CC@-@@AA@+<<======	AS:i:96	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:48	YS:i:82	YT:Z:CP
M01228:25:000000000-A1CW0:1:1104:4137:13019	163	16S	785	44	1S41M	=	880	144	TGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGA	9=<=9>-<<@@@<@@-A8.AACEE-6+A>+7+AE9C>EC.CE	AS:i:82	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:41	YS:i:96	YT:Z:CP
M01228:25:000000000-A1CW0:1:1104:10078:23947	97	16S	783	44	142M5D7M1S	=	760	-177	CAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACCCTGCACG	???A?BBBDDDDDDDDGGGGGGIIHHHHIHHEHHIHHIIHHIHHHHHIFHHHHHDFHHHHHIHHIIHHHEEHHHHHHEHDFDFFHHGGDGGD=DGGGGBEGGGGGG;;8CEGGGGGGGGGGGGGGGE?EEGGEEEGCEGGGGGGECCECG	AS:i:284	XN:i:0	XM:i:1	XO:i:1	XG:i:5	NM:i:6	MD:Z:142^GGGGG2C4	YS:i:121	YT:Z:DP
M01228:25:000000000-A1CW0:1:1104:10078:23947	145	16S	760	44	62M4S	=	783	177	GGTGCGAAAGCATGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGCACG	>CDC>-??C-A-HDFCC=FFGCAHHGHHHFHHGFEHHF;HHFFFGGGFGGDDDDDDDDBB??????	AS:i:121	XN:i:0	XM:i:1	XO:i:0	XG:i:0NM:i:1	MD:Z:11G50	YS:i:284	YT:Z:DP
M01228:25:000000000-A1CW0:1:1104:11671:25136	83	16S	989	44	18M5I33M	=	829	-212	TCTTGACATCCACGGAAGTTTTTTTTTTTTTTTTGAGAATGTGCCTTCAGGAACCG	?,BC,5??+CCAEBFFD@CCCHHHHEHCCHHHHHFFFFBFDDDDDBDDBB??<???	AS:i:64	XN:i:0	XM:i:7	XO:i:1	XG:i:5	NM:i:12	MD:Z:22C0A0G0A0G0A15G7	YS:i:292	YT:Z:CP
M01228:25:000000000-A1CW0:1:1104:11671:25136	163	16S	829	44	1S146M	=	989	212	CGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGA	?????@?BBDDDDDDDFFFFFFEHEEHHIHHHHHHIHH=DCHCEDFF5CEFCACEHDFFHH;C3C=DEHBB5@@BF;?EEFFEEFEFFCEFEFFECBEEDDD>;DDD>48:A8AD?;EACCEEE*:A:8AEE?EA?A?EEEF:.'42	AS:i:292	XN:i:0	XM:i:0	XO:i:0XG:i:0	NM:i:0	MD:Z:146	YS:i:64	YT:Z:CP
M01228:25:000000000-A1CW0:1:1104:19760:26407	99	16S	199	41	2S49M24D19M	=	265	144	AAAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGGCGACGATCCCTAGC	?<????@@DDDDBDDDFDFCFFHHFHHHFHHHHHHHFC=EF@CGHFFFEFFEC5AF?C>CCEDDDHFGFG	AS:i:107	XN:i:0	XM:i:0XO:i:1	XG:i:24	NM:i:24	MD:Z:49^CTAGTAGGTGGGGTAACGGCTCAC19	YS:i:113	YT:Z:CP
M01228:25:000000000-A1CW0:1:1104:19760:26407	147	16S	265	41	3S64M3D6M	=	199	-144	ATGGGATTAGCTAGGCGACGATCCCTCGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACGCCCAG	HDHFDGFA5AA>E>5>5CEFFE>>>++FCAC=CCFFCE?/FGDCA,E;FHFFFCFFFBB?DDDDB@@??????	AS:i:113	XN:i:0XM:i:5	XO:i:1	XG:i:3	NM:i:8	MD:Z:2C1C1C16A40^ACG1T4	YS:i:107	YT:Z:CP
M01228:25:000000000-A1CW0:1:1104:18521:26491	99	16S	1257	14	67M4I13M4D7M21S	=	1319	202	AGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCAGACTCCAATCCGGACTACGACGCACTTTATGAGGTCTCGCGAGAGCAAGCGGACCTCATAAAGT	?????@@??55@=@=BGGGFA@C>>CFBCEGHHGCGCCCC>5?E?FFGHFCCEHIFGHFHHDDEFFEEEC=EGGGGFDE=BEDDEGEDGEEBB;ACA;?>48>EEGCEGEGE	AS:i:100	XS:i:58	XN:i:0	XM:i:14	XO:i:2	XG:i:8	NM:i:22	MD:Z:47G2T1G0G1G2T1C1A0C0T5C0C4A3^GGAA4T2	YS:i:266	YT:Z:CP
M01228:25:000000000-A1CW0:1:1104:18521:26491	147	16S	1319	14	95M2D43M	=	1257	-202	ACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCTGGGCCTTGTACACACCGCCCGTCACACCAGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGG	E<82::?1:CGC?GGGGGGGECB>>EGGGEGE<CEDDGGDGEGGEGEGD;GGFBDFFDFHHHFDC+EHHHHFGGGCC5<*HD>HHHIHECA9IHIHFEDEDHFHF=HHHHFCHHIIGGFFFCBBBB<DD@?<??????	AS:i:266	XN:i:0	XM:i:1	XO:i:1	XG:i:2NM:i:3	MD:Z:65C29^TG43	YS:i:100	YT:Z:CP
)";

	std::string test_dir = getTestFileDir().string();
	std::ofstream dummy(test_dir+"/tmp/tmp_e.sam");
	dummy << sample_reads;
	dummy.close();
	EXPECT_NO_THROW(
		parseSAM(test_dir+"/tmp/tmp_e.sam",
				 test_dir+"/tmp/tmp_e.mut",
				 "", // debug_outname
				 "", // primers_filename
				 800, // max_paired_fragment_length
				 10, // min_mapq
				 false, // right_align_ambig_dels,
				 false, //right_align_ambig_ins,
				 5, // max_internal_match,
				 30, //min_qual,
				 7, // exclude_3prime,
				 "", //mutation_type,
				 false, // variant_mode,
				 false, // trim_primers
				 false, // require_forward_primer_mapped
				 false, // require_reverse_primer_mapped
				 10, // max_primer_offset
				 false, // input_is_unpaired
				 false, // debug,
				 true //warn_on_no_mapped
		);
	);
}



// ugh. checking opposite strand basecalls from mutation group left->right in the absence
// of mutations in one strand causes problems if a mutation in one strand falls next to
// the mapping end on the other strand
// - should be working now
TEST(OverlapResolution, DroppedNearEnd) {
	std::string sample_reads = R"(M01228:25:000000000-A1CW0:1:1101:17776:2377:R1	83	TPP	29	42	81M4I3M3I25M5S	=	1	-146	CTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTTCCCTACACGACCCTCCCCATCGGGCTTCGGTCCGGTTCTTTTC	GEDD>C0:ECC?ACA>C4=CCA>8@3;,GGGGGGED;B;DEDD;DDBF:E@;EC?CBHF?BFFEFC<+C=C<7;HFFFFDDECECC55CAEHHFC9DHHGGGGGGEEEEEEDDBAB?????	AS:i:181	XN:i:0	XM:i:5	XO:i:2	XG:i:7	NM:i:12	MD:Z:77C0G2G5A0A20	YS:i:210	YT:Z:CP
M01228:25:000000000-A1CW0:1:1101:17776:2377:R2	163	TPP	1	42	4S109M14S	=	29	146	AAAAGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTTCCCTACACGACCCTCCC	????????DDDDDDDDEFFFFFFIHHHHHHHIIIFHHFHHHHHIIIIFHHIIIHHIHIHHEDFHGFFFHHIIIHHHHHHHHHHHHHHHHHHHHFDFFFFFDEDDDEDEEEFAEAAB=A*2:@8A*AE	AS:i:210	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:105C0G2	YS:i:181	YT:Z:CP
)";
	std::string r1_line = R"([read]	PAIRED_R1	28	136	-	INCLUDED	-999	CTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	GEDD>C0:ECC?ACA>C4=CCA>8@3;,GGGGGGED;B;DEDD;DDBF:E@;EC?CBHF?BFFEFC<+C=C<7;HFFFFDDC55HHFC9DHHGGGGGGEEEEEEDDBAB				104 106 "T" "F" "" 105 107 "C" "F" "" 108 109 "TACA" "ECEC" "" 108 110 "C" "C" "" 111 112 "CCC" "CAE" "" 114 116 "C" "C" "" 115 117 "C" "9" ""
)";
	std::string r2_line = R"([read]	PAIRED_R2	0	108	+	INCLUDED	-999	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCC	????DDDDDDDDEFFFFFFIHHHHHHHIIIFHHFHHHHHIIIIFHHIIIHHIHIHHEDFHGFFFHHIIIHHHHHHHHHHHHHHHHHHHHFDFFFFFDEDDDEDEEEFAE				104 106 "T" "E" "" 105 107 "C" "F" ""
)";
	std::string exp_merged_line = R"([read]	PAIRED	0	136	+	INCLUDED	-999	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	????DDDDDDDDEFFFFFFIHHHHHHHIIIFHHFHHHHHIIIIFHHIIIHHIHIHHGGGHGGFFHHIIIHHHHHHHHHHHHHHHHHHHHFFFFFFFDEDDDEHFFFFDEC55HHFC9DHHGGGGGGEEEEEEDDBAB	11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111	11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111		104 106 "T" "F" "" 105 107 "C" "F" "" 108 109 "TACA" "ECEC" "" 108 110 "C" "C" "" 111 112 "CCC" "CAE" "" 114 116 "C" "C" "" 115 117 "C" "9" ""
)";
	Read r1 = Read(r1_line);
	Read r2 = Read(r2_line);
	std::vector<Read> reads;
	reads.push_back(r1);
	reads.push_back(r2);
	Read merged = mutation::mergeMatePairs(reads);
	//std::cout << "\n" << merged << "\n";

    // FIXME: Read constructor doesn't 
	Read exp_merged = Read(exp_merged_line);
	EXPECT_EQ(exp_merged_line, merged.toString());
}

// - should be working now
TEST(OverlapResolution, DroppedInsert) {
	std::string sample_reads = R"(M01228:25:000000000-A1CW0:1:1101:13750:2195:R1	83	TPP	1	44	3S63M1D39M1I34M5S	=	1	-146	CACGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCACCCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTTTGCCGGATCCAAATCGGGCTTCGGTCCGGTTCTACAT	DDDEA22D?>8A:1??A?C8>D?DEEEFFEC?E?EEEAEFFFFFEFEEFFFEEEBEEEEEFFEEED?EHHHHHFHHHFDHHFCHFHHIIHHHFHIHEEHHHHHEFHHHEHHIIIIHHFIHHHHFFFFFFDDDDDDDDB?/?????	AS:i:244	XN:i:0	XM:i:4	XO:i:2	XG:i:2	NM:i:6	MD:Z:63^T0G0A0T38C31	YS:i:184	YT:Z:CP
M01228:25:000000000-A1CW0:1:1101:13750:2195:R2	163	TPP	1	44	4S63M1D38M	=	1	146	GCACGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCACCCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCG	?,<?9@<@AB9A?@@BEFFFF>CEFCCH>CEHHFF@F@D>EHEHHFD?EFF-EGCG-ACD5CEFF-CE?<<?DFBDFHFHFHHFFH+=A:@BFDDD;B,@EEEEE	AS:i:184	XN:i:0	XM:i:3	XO:i:1	XG:i:1	NM:i:4	MD:Z:63^T0G0A0T35	YS:i:244	YT:Z:CP
)";
	std::string r1_line = R"([read]	PAIRED_R1	0	136	-	INCLUDED	-999	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	EA22D?>8A:1??A?C8>D?DEEEFFEC?E?EEEAEFFFFFEFEEFFFEEEBEEEEEFFEEED!?EHHHHHFHHHFDHHFCHFHHIIHHHFHIHEEHHHHHEFHHEHHIIIIHHFIHHHHFFFFFFDDDDDDDDB?/				62 64 "" "" "" 63 65 "A" "?" "" 64 66 "C" "E" "" 65 67 "C" "H" "" 102 105 "TTT" "HHH" "" 104 106 "T" "E" ""
)";
	std::string r2_line = R"([read]	PAIRED_R2	0	101	+	INCLUDED	-999	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCG	9@<@AB9A?@@BEFFFF>CEFCCH>CEHHFF@F@D>EHEHHFD?EFF-EGCG-ACD5CEFF-C!E?<<?DFBDFHFHFHHFFH+=A:@BFDDD;B,@EEEEE				62 64 "" "" "" 63 65 "A" "E" "" 64 66 "C" "?" "" 65 67 "C" "<" ""
)";
	std::string exp_merged_line = R"([read]	PAIRED	0	136	+	INCLUDED	-999	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	EA<@DB>AA@@BEFFFF>DEFEEHFFEHHFFEFEDEFHFHHFFEEFFFEGEGEEEEEFFFFED!EEHHHHHFHHHFHHHHFHHHHIIHHHFHIHEEHHHHHEFHHEHHIIIIHHFIHHHHFFFFFFDDDDDDDDB?/	11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111	11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111		62 64 "" "" "" 63 65 "A" "?" "" 64 66 "C" "E" "" 65 67 "C" "H" "" 102 105 "TTT" "HHH" "" 104 106 "T" "E" ""
)";
	Read r1 = Read(r1_line);
	Read r2 = Read(r2_line);
	std::vector<Read> reads;
	reads.push_back(r1);
	reads.push_back(r2);
	Read merged = mutation::mergeMatePairs(reads);
	//std::cout << "\n" << merged << "\n";

	Read exp_merged = Read(exp_merged_line);
	EXPECT_EQ(exp_merged_line, merged.toString());
}

// working now
TEST(OverlapResolution, DroppedMismatch) {
	std::string sample_reads = R"(M01228:25:000000000-A1CW0:1:1101:14363:1997:R1	83	TPP	84	42	30M1D23M5S	=	1	-146	TAAGGAAGTTCTCGATCCGGTTCGCCGGATCAAATCGGGCTTCGGTCCGGTTCGTGCG	=D?4,CD?C?@D@AC)CCC7EEEC=@DCA;FFHHEHFFCDDDDDBDBB@D??@?????	AS:i:96	XN:i:0	XM:i:1XO:i:1	XG:i:1	NM:i:2	MD:Z:2G27^C23	YS:i:62	YT:Z:CP
M01228:25:000000000-A1CW0:1:1101:14363:1997:R2	163	TPP	1	42	4S31M	=	84	146	CCCGGGCCTTCGGGCCAAGGACTCGGGGTGCCCTT	,5<??@@@BBBAB==@66;9CFHFFHHBEE:EFHC	AS:i:62	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:31	YS:i:96	YT:Z:CP
)";
	std::string r1_line = R"([read]	PAIRED_R1	83	136	-	INCLUDED	-999	TAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	=D?4,CD?C?@D@AC)CCC7EEEC=@DCA;!FFHHEHFFCDDDDDBDBB@D??@				84 86 "A" "?" "" 112 115 "C" "F" ""
)";
	std::string r2_line = R"([read]	PAIRED_R2	0	30	+	INCLUDED	-999	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTT	?@@@BBBAB==@66;9CFHFFHHBEE:EFHC				
)";
	std::string exp_merged_line = R"([read]	PAIRED	0	136	+	INCLUDED	-999	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTT____________________________________________________TAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	?@@@BBBAB==@66;9CFHFFHHBEE:EFHC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~=D?4,CD?C?@D@AC)CCC7EEEC=@DCA;!FFHHEHFFCDDDDDBDBB@D??@	11111111111111111111111111111110000000000000000000000000000000000000000000000000000111111111111111111111111111111111111111111111111111111	11111111111111111111111111111110000000000000000000000000000000000000000000000000000111111111111111111111111111111111111111111111111111111		84 86 "A" "?" "" 112 115 "C" "F" ""
)";
	Read r1 = Read(r1_line);
	Read r2 = Read(r2_line);
	std::vector<Read> reads;
	reads.push_back(r1);
	reads.push_back(r2);
	Read merged = mutation::mergeMatePairs(reads);
	//std::cout << "\n" << merged << "\n";

	Read exp_merged = Read(exp_merged_line);
	EXPECT_EQ(exp_merged_line, merged.toString());
}

// ambiguously aligned complex long deletion is aligned differently between R1 and R2,
// and on top of that is not identified as ambiguous, because in both cases it's aligned
// as two mutations separated by one ref nuc
// - The deletion component is then mistakenly dropped
// - working now
TEST(OverlapResolution, UnidentifiedAmbiguousMut) {
	std::string sample_reads = R"(M01228:25:000000000-A1CW0:1:1101:14903:2861:R1	83	TPP	1	36	24S35M1I59M13D2M3D25M5S	=	1	186	GACGTGTGCTCTTCCGATCTTCAAGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTTCTGTGAAGGCTGAGAAATACCCGTATCACCTGATCTAGATAATGCCAGCGTAGGGAAGTTCCCTACACATCGGGCTTCGGTCCGGTTCACTGC	?CDGGEC?*C8EECGGEEC*GGGGGEGDE:6GGEEGGGGEGEEGGGGGGFEFFDDDFBHHHHHHHHHHHHHHHHGECHHIIIIIHGHFDIIHHHIIIIIIIHHHIHHIIIIIIHHFFHIIHHHFHEHHHGGFDGGDEEEDDEEBBB?AAA?	AS:i:193	XN:i:0	XM:i:4	XO:i:3	XG:i:17	NM:i:21	MD:Z:33G35G24^TCGATCCGGTTCG2^GGA1C2A20	YS:i:193	YT:Z:CP
M01228:25:000000000-A1CW0:1:1101:14903:2861:R2	163	TPP	1	36	4S35M1I61M16D25M25S	=	1	-186	TCAAGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTTCTGTGAAGGCTGAGAAATACCCGTATCACCTGATCTAGATAATGCCAGCGTAGGGAAGTTCCCTACACATCGGGCTTCGGTCCGGTTCACTGCAGATCGGAAGAGCGTCGTGT	AAAAABAADEEDDDDDGGGGGGIIIHHHHHHIIIIIIIIIIIIIIIIIIHIHFHHIIIHIHHHIIHIIIIHIIIIIHIFHIHHHHGFHHHDHHHDEFHHGGGGGGGHGGD=8>EGGG?CEG6<BEEGEGGEGCEGCC8?)8??C8>?8'00	AS:i:193	XN:i:0	XM:i:5	XO:i:2	XG:i:17	NM:i:22	MD:Z:33G35G24T1^GATCCGGTTCGCCGGA1C2A20	YS:i:193	YT:Z:CP
)";
	std::string r1_line = R"([read]	PAIRED_R1	0	136	-	INCLUDED	-999	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	GEGDE:6GGEEGGGGEGEEGGGGGGFEFFDDDFBHHHHHHHHHHHHHHHGECHHIIIIIHGHFDIIHHHIIIIIIIHHHIHHIIIIIIHHFFHI!!!!!!!!!!!!!IH!!!HHFHEHHHGGFDGGDEEEDDEEBBB				32 34 "T" "B" "" 34 35 "T" "H" "" 68 70 "A" "I" "" 93 107 "" "" "" 108 112 "" "" "" 112 114 "A" "H" "" 115 117 "C" "E" ""
)";
	std::string r2_line = R"([read]	PAIRED_R2	0	136	+	INCLUDED	-999	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	ABAADEEDDDDDGGGGGGIIIHHHHHHIIIIIIIIIIIIIIIIIHIHFHHIIIHIHHHIIHIIIIHIIIIIHIFHIHHHHGFHHHDHHHDEFHHGG!!!!!!!!!!!!!!!!GGGGGHGGD=8>EGGG?CEG6<BEE				32 34 "T" "I" "" 34 35 "T" "I" "" 68 70 "A" "I" "" 93 95 "C" "G" "" 95 112 "" "" "" 112 114 "A" "G" "" 115 117 "C" "G" ""
)";
	std::string exp_merged_line = R"([read]	PAIRED	0	136	+	INCLUDED	-999	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	GEGDEEEGGEEGGGGGGGIIIHHHHHHIIIIIIIIIIIIIIIIIHIHHHHIIIHIIIIIIHIIIIIIIIIIIIIIIHHHIHHIIIIIIHHFFHIGG!!!!!!!!!!!IH!!!HHGHGHHHGGFDGGGGEEEGEEBEE	11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111	11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111		32 34 "T" "I" "" 34 35 "T" "I" "" 68 70 "A" "I" "" 93 107 "" "" "" 108 112 "" "" "" 112 114 "A" "H" "" 115 117 "C" "G" ""
)";
	Read r1 = Read(r1_line);
	Read r2 = Read(r2_line);
	std::vector<Read> reads;
	reads.push_back(r1);
	reads.push_back(r2);
	Read merged = mutation::mergeMatePairs(reads);
	//std::cout << "\n" << merged << "\n";

	Read exp_merged = Read(exp_merged_line);
	EXPECT_EQ(exp_merged_line, merged.toString());
}

// 20549:2603
// especially problematic. Long deletion gets dropped, despite no real coverage of deleted region.
// - ack. fixed bug that would usually favor R2, so need to update a bunch of tests
// - now need to implement a quality comparison when a mutation group is present on one strand
//   but not on the other (here a low-quality basecall in R2 is making it through)
// - should be working now
TEST(OverlapResolution, ConflictingEnds) {
	std::string sample_reads = R"(M01228:25:000000000-A1CW0:1:1101:20549:2603:R1	83	TPP	7	42	78M1D9M18D25M5S	=	1	-146	CGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTGGGAAGTTCCCTACATCGGGCTTCGGTCCGGTTCCTAAG	ECA>AEACGECEEGGGGGEE@FDGGBHD+HHHFHHFBHHHHHHHCEEHFIHFHHHIHIIIIIIIIIHFIHFHFHHFIIIIIIIIHHIHHFFEHHHGGGGGGDDDDBDEDBBBAAAAA	AS:i:183	XN:i:0	XM:i:3	XO:i:2	XG:i:19	NM:i:22	MD:Z:78^A9^TCGATCCGGTTCGCCGGA0T1C1A20YS:i:177	YT:Z:CP
M01228:25:000000000-A1CW0:1:1101:20549:2603:R2	163	TPP	1	42	4S84M1D13M3S	=	7	146	ACGGGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTNCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTGGGAAGTTCCCTACAT	?????BB@DDEDDDDDGGGGGGIIIHGHFHHHIIIII#55CDFFHHFHHIHIIIIIIIIHHEHIIFHHFHHHHHHHHHHFFHHHFHHHEA5AFGDDDDEEGGGG	AS:i:177	XN:i:0	XM:i:3	XO:i:1	XG:i:1	NM:i:4	MD:Z:33G50^A9T1G1	YS:i:183	YT:Z:CP
)";
	std::string r1_line = R"([read]	PAIRED_R1	6	136	-	INCLUDED	-999	CGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	ECA>AEACGECEEGGGGGEE@FDGGBHD+HHHFHHFBHHHHHHHCEEHFIHFHHHIHIIIIIIIIIHFIHFHFHHFII!IIIIIIHHI!!!!!!!!!!!!!!!!!!HHFFEHHHGGGGGGDDDDBDEDBBB				83 85 "" "" "" 93 114 "CC" "HH" "" 113 115 "T" "F" "" 115 117 "C" "E" ""
)";
	std::string r2_line = R"([read]	PAIRED_R2	0	97	+	INCLUDED	-999	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGA	?BB@DDEDDDDDGGGGGGIIIHGHFHHHIIIII#55CDFFHHFHHIHIIIIIIIIHHEHIIFHHFHHHHHHHHHHFFHHHFHHH!EA5AFGDDDDEEG				32 34 "N" "#" "" 83 85 "" "" "" 93 95 "C" "D" "" 95 97 "T" "E" ""
)";
	std::string exp_merged_line = R"([read]	PAIRED	0	136	+	INCLUDED	-999	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	?BB@DDEDDDDEGGGGGGIIIHGHFHHHIIIIID5HHHFHHHFHHIHIIIIIIIIIHFHIIIHIIIIIIIIIHHIHFHHHHHII!IIIIIIHHIDEEG!!!!!!!!!!!!!!HHFFEHHHGGGGGGDDDDBDEDBBB	11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111	11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111		83 85 "" "" "" 93 114 "CC" "HH" "" 113 115 "T" "F" "" 115 117 "C" "E" ""
)";
	Read r1 = Read(r1_line);
	Read r2 = Read(r2_line);
	std::vector<Read> reads;
	reads.push_back(r1);
	reads.push_back(r2);
	Read merged = mutation::mergeMatePairs(reads);
	//std::cout << "\n" << merged << "\n";

	Read exp_merged = Read(exp_merged_line);
	exp_merged.toString();
	merged.toString();
	EXPECT_EQ(exp_merged_line, merged.toString());
}

// 18159:2928 seems to pick wrong strand for second mutation
// working now
TEST(OverlapResolution, LowerQualityStrandSelected) {
	std::string sample_reads = R"(M01228:25:000000000-A1CW0:1:1101:18159:2928:R1	83	TPP	1	41	32S19M23D95M5S	=	1	176	GGAGTTCAGACGTGTGCTCTTCCGATCGCACTGGCCTTCGGGCCAAGGACTCTGAGAAATACCCGTATCACCTGACCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTCCGCTA	CGCEGGEECC?GGGGCC9GGECD?8EGGGEGGGGDA?AGGGGGEGGGGGGGGGGGGGGGGGGGGGGGGGEGGGGGGGGGGGGGGGGHGGGGGHHHHHHHHHHHHHHHHHHHHHHHHHIIIIIIIIHHHHGGGGGGDEEDEEEEAAAAAAAA	AS:i:196	XN:i:0	XM:i:1	XO:i:1	XG:i:23	NM:i:24	MD:Z:19^CGGGGTGCCCTTCTGCGTGAAGG24T70	YS:i:196	YT:Z:CP
M01228:25:000000000-A1CW0:1:1101:18159:2928:R2	163	TPP	1	41	4S19M23D95M7S	=	1	-176	CACTGGCCTTCGGGCCAAGGACTCTGAGAAATACCCGTATCACCTGACCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTCCGCTAAG	<<?<<B?BBBB?BBBBCF>FCC>EBC>ACFB9ACGHDC>FFFFHDDEDFFBFCGHGGDDGGFHE@EHHHCEEFHBEH:EFHECEHHD)CE:=)=5BBDEBBD@@6:AB)?BEEBE?;CE)??;A:	AS:i:196	XN:i:0	XM:i:1	XO:i:1	XG:i:23	NM:i:24	MD:Z:19^CGGGGTGCCCTTCTGCGTGAAGG24T70	YS:i:196	YT:Z:CP
)";
	std::string r1_line = R"([read]	PAIRED_R1	0	136	-	INCLUDED	-999	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	GGDA?AGGGGGEGGGGGGG!!!!!!!!!!!!!!!!!!!!!!!GGGGGGGGGGGGGGGGGGEGGGGGGGGGGGGGGGGHGGGGGHHHHHHHHHHHHHHHHHHHHHHHHHIIIIIIIIHHHHGGGGGGDEEDEEEEAAA	11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111	11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111		18 43 "C" "G" "" 65 67 "C" "G" ""
)";
	std::string r2_line = R"([read]	PAIRED_R2	0	136	+	INCLUDED	-999	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	<B?BBBB?BBBBCF>FCC>!!!!!!!!!!!!!!!!!!!!!!!EBC>ACFB9ACGHDC>FFFFHDDEDFFBFCGHGGDDGGFHE@EHHHCEEFHBEH:EFHECEHHD)CE:=)=5BBDEBBD@@6:AB)?BEEBE?;C	11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111	11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111		18 43 "C" "E" "" 65 67 "C" "D" ""
)";
	std::string exp_merged_line = R"([read]	PAIRED	0	136	+	INCLUDED	-999	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	GGDBBBGGGGGEGGGGGGG!!!!!!!!!!!!!!!!!!!!!!!GGGGGGGGGGGGHGGGGGFGHGGGGGGGGGGHGGGHGGGHGHHHHHHHHHHHHHHHHHHHHHHHHHIIIIIIIIHHHHGGGGGGDEEDEEEEAAC	11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111	11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111		18 43 "C" "G" "" 65 67 "C" "G" ""
)";
	Read r1 = Read(r1_line);
	Read r2 = Read(r2_line);
	std::vector<Read> reads;
	reads.push_back(r1);
	reads.push_back(r2);
	Read merged = mutation::mergeMatePairs(reads);
	//std::cout << "\n" << merged << "\n";

	Read exp_merged = Read(exp_merged_line);
	EXPECT_EQ(toString(exp_merged.mutations), toString(merged.mutations));
	EXPECT_EQ(exp_merged_line, merged.toString());
}

// 19057:2147 wrong strand for second mutation
// - should be working now
TEST(OverlapResolution, LowerQualityStrandSelectedB) {
	std::string sample_reads = R"(M01228:25:000000000-A1CW0:1:1101:19057:2147:R1	83	TPP	1	44	10S35M1D101M5S	=	1	152	CGATCTAGATGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCTGAAGGCTGAGAAATACCCGCATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTCCTACC	CFFFFEEEC?8E>?>8'EEEEA?CEECDDD?AAA8EAA:?EAEEEEAACEEEEEBEFFEEEBDEFFFEFFFFFFFEFEFFFFFFFEFEHHHFHHHHHHHHHHHEHHHHHHHHHHHHHIIIIIIIIHHHHFFFFFFDDDDDDDDBB??????	AS:i:262	XN:i:0	XM:i:1	XO:i:1	XG:i:1	NM:i:2	MD:Z:35^G20T80	YS:i:240	YT:Z:CP
M01228:25:000000000-A1CW0:1:1101:19057:2147:R2	163	TPP	1	44	4S35M1D90M	=	1	-152	AGATGGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCTGAAGGCTGAGAAATACCCGCATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTT	?AA??BBBEEEDDDDDGGGGGGIIIGHBCCHHHHHFHHHIIHIIIIIIDFHFHHHEGHHHHHEHHHHIIIIIHIHFHHHHHHHHFHDDEHFHFGFGFGGGDBDEGB@EDE@EEEGGGEGGE:C??822C	AS:i:240	XN:i:0	XM:i:1	XO:i:1	XG:i:1	NM:i:2	MD:Z:35^G20T69	YS:i:262	YT:Z:CP
)";
	std::string r1_line = R"([read]	PAIRED_R1	0	136	-	INCLUDED	-999	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	8E>?>8'EEEEA?CEECDDD?AAA8EAA:?EAEEE!EAACEEEEEBEFFEEEBDEFFFEFFFFFFFEFEFFFFFFFEFEHHHFHHHHHHHHHHHEHHHHHHHHHHHHHIIIIIIIIHHHHFFFFFFDDDDDDDDBB?	11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111	11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111		34 36 "" "" "" 55 57 "C" "F" ""
)";
	std::string r2_line = R"([read]	PAIRED_R2	0	125	+	INCLUDED	-999	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTT	?BBBEEEDDDDDGGGGGGIIIGHBCCHHHHHFHHH!IIHIIIIIIDFHFHHHEGHHHHHEHHHHIIIIIHIHFHHHHHHHHFHDDEHFHFGFGFGGGDBDEGB@EDE@EEEGGGEGGE:C??822C	111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111	111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111		34 36 "" "" "" 55 57 "C" "H" ""
)";
	std::string exp_merged_line = R"([read]	PAIRED	0	136	+	INCLUDED	-999	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	?EBBEEEEEEEDGGGGGGIIIGHBCEHHHHHFHHH!IIHIIIIIIDFHFHHHEGHHHHHFHHHHIIIIIHIHFHHHHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHHHIIIIIIIIHHHHFFFFFFDDDDDDDDBB?	11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111	11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111		34 36 "" "" "" 55 57 "C" "H" ""
)";
	Read r1 = Read(r1_line);
	Read r2 = Read(r2_line);
	std::vector<Read> reads;
	reads.push_back(r1);
	reads.push_back(r2);
	Read merged = mutation::mergeMatePairs(reads);
	//std::cout << "\n" << merged << "\n";

	Read exp_merged = Read(exp_merged_line);
	EXPECT_EQ(toString(exp_merged.mutations), toString(merged.mutations));
	EXPECT_EQ(exp_merged_line, merged.toString());
}

// note: no more exclude_3prime parameter, and needs an initial depth
// - initial depths calced in stripEnd() and stripPrimers()
// - depths are now one longer on the end compared to previous implementation

TEST(filterQscoresCountDepths, simple){
	/*
     *
     *        q-scores: HHHHHHHHHHHHHHHH
     *    aligned read: AATTGGCCATGCCGTA
     *          target: AATTGGCCATGCCGTA
     *
     *
     * "H" = phred score 39
     * "!" = phred score 0
     * "#" = phred score 2
     */

	std::string seq =  "AATTGGCCATGCCGTA";
	std::string qual = "HHHHHHHHHHHHHHHH";
	std::vector<Mutation> mutations{};
	int left = 0;
	int min_qual = 30;

	std::string exp_depth = "1111111111111111";
	std::string exp_count = "0000000000000000";
	std::vector<Mutation> exp_included_mutations = {};
	std::vector<Mutation> exp_excluded_mutations = {};

	std::vector<bool> depth;
	std::vector<bool> count;
	std::vector<bool> initial_depth = stringToBoolVect("1111111111111111");
	std::vector<Mutation> included_mutations;
	std::vector<Mutation> excluded_mutations;
	boost::tie(depth,
			   count,
			   included_mutations,
			   excluded_mutations) = filterQscoresCountDepths(mutations,
															  seq,
															  qual,
															  initial_depth,
															  left,
															  min_qual,
															  "",
															  false); // variant_mode off

	EXPECT_EQ(exp_depth, toString(depth));
	EXPECT_EQ(exp_count, toString(count));
}


TEST(filterQscoresCountDepths, variant_mode_on){
	/*
     *
     *        q-scores: HHHHHHHHHHHHHHHH
     *    aligned read: AATTGGCC--GCCGTA
     *          target: AATTGGCCATGCCGTA
     *
     *
     * "H" = phred score 39
     * "!" = phred score 0
     * "#" = phred score 2
     */

	std::string seq =  "AATTGGCCATGCCGTA";
	std::string qual = "HHHHHHHH!!HHHHHH";
	std::vector<Mutation> mutations{{7, 10, "", "", ""}};
	int left = 0;
	int min_qual = 30;

	std::string exp_depth = "1111111111111111";
	std::string exp_count = "0000000001000000";
	std::vector<Mutation> exp_included_mutations = {};
	std::vector<Mutation> exp_excluded_mutations = {};

	std::vector<bool> depth;
	std::vector<bool> count;
	std::vector<bool> initial_depth = stringToBoolVect("1111111111111111");
	std::vector<Mutation> included_mutations;
	std::vector<Mutation> excluded_mutations;
	boost::tie(depth,
			   count,
			   included_mutations,
			   excluded_mutations) = filterQscoresCountDepths(mutations,
															  seq,
															  qual,
															  initial_depth,
															  left,
															  min_qual,
															  "",
															  true); // variant_mode on
	EXPECT_EQ(exp_depth, toString(depth));
	EXPECT_EQ(exp_count, toString(count));
}

TEST(filterQscoresCountDepths, variant_mode_off){
	/*
     *
     *        q-scores: HHHHHHHHHHHHHHHH
     *    aligned read: AATTGGCC--GCCGTA
     *          target: AATTGGCCATGCCGTA
     *
     *
     * "H" = phred score 39
     * "!" = phred score 0
     * "#" = phred score 2
     */

	std::string seq =  "AATTGGCCATGCCGTA";
	std::string qual = "HHHHHHHH!!HHHHHH";
	std::vector<Mutation> mutations{{7, 10, "", "", ""}};
	int left = 0;
	int exclude_3prime = 1;
	int min_qual = 30;

	std::string exp_depth = "1111111101111111";
	std::string exp_count = "0000000001000000";
	std::vector<Mutation> exp_included_mutations = {};
	std::vector<Mutation> exp_excluded_mutations = {};

	std::vector<bool> depth;
	std::vector<bool> count;
	std::vector<bool> initial_depth = stringToBoolVect("1111111111111111");
	std::vector<Mutation> included_mutations;
	std::vector<Mutation> excluded_mutations;
	boost::tie(depth,
			   count,
			   included_mutations,
			   excluded_mutations) = filterQscoresCountDepths(mutations,
															  seq,
															  qual,
															  initial_depth,
															  left,
															  min_qual,
															  "",
															  false); // variant_mode off
	EXPECT_EQ(exp_depth, toString(depth));
	EXPECT_EQ(exp_count, toString(count));
}

TEST(filterQscoresCountDepths, variant_mode_mismatch){
	/*
     *
     *        q-scores: HHHHHHHHHHHHHHHH
     *    aligned read: AATTGGCG-TGCCGTA
     *          target: AATTGGCCATGCCGTA
     *
     *
     * "H" = phred score 39
     * "!" = phred score 0
     * "#" = phred score 2
     */

	std::string seq =  "AATTGGCCATGCCGTA";
	std::string qual = "HHHHHHHH!HHHHHHH";
	std::vector<Mutation> mutations{{6, 8, "G", "H", ""},
									{7, 9, "", "", ""}};
	int left = 0;
	int min_qual = 30;

	std::string exp_depth = "1111111111111111";
	std::string exp_count = "0000000110000000";
	std::vector<Mutation> exp_included_mutations = {};
	std::vector<Mutation> exp_excluded_mutations = {};

	std::vector<bool> depth;
	std::vector<bool> count;
	std::vector<bool> initial_depth = stringToBoolVect("1111111111111111");
	std::vector<Mutation> included_mutations;
	std::vector<Mutation> excluded_mutations;
	boost::tie(depth,
			   count,
			   included_mutations,
			   excluded_mutations) = filterQscoresCountDepths(mutations,
															  seq,
															  qual,
															  initial_depth,
															  left,
															  min_qual,
															  "",
															  true); // variant_mode on
	EXPECT_EQ(exp_depth, toString(depth));
	EXPECT_EQ(exp_count, toString(count));
}
/*
boost::tuple<std::vector<Mutation>,
		std::vector<bool>>
		stripEnd(const Read &r,
				 const std::vector<bool> &depth_in,
				 const int exclude_length,
				 const int which_end,
				 const bool debug)
*/

TEST(stripEnd, largeEndStrip){
	/*
     *
     *        q-scores: HHHHHHHHHHHHHHHH
     *    aligned read: AATTGGCCATGCCGTA
     *          target: AATTGGCCATGCCGTA
     *
     *
     * "H" = phred score 39
     * "!" = phred score 0
     * "#" = phred score 2
     */

	std::string seq =  "AATTGGCCATGCCGTA";
	std::string qual = "HHHHHHHHHHHHHHHH";
	std::vector<bool> initial_depth = stringToBoolVect("1111111111111111");
	std::vector<Mutation> mutations{};
	Read r(0, 15, seq);
	r.setQual(qual).setMutations(mutations);
	int exclude_length = 30;

	std::vector<Mutation> included_mutations;
	std::vector<bool> depth;

	std::string exp_depth = "0000000000000000";
	std::vector<Mutation> exp_included_mutations = {};

	boost::tie(included_mutations,
			   depth) = stripEnd(r,
								 initial_depth,
								 exclude_length,
								 RIGHT,
								 false);

	EXPECT_EQ(exp_depth, toString(depth));
	EXPECT_EQ(exp_included_mutations, included_mutations);
}

TEST(stripEnd, endStripWithMutations){
	/*
     *
     *        q-scores: HHHHHHHHHHHHHHHH
     *    aligned read: AATTGGCCGTGCCGTA
     *          target: AATTGGCCATGCCGTA
     *
     *
     * "H" = phred score 39
     * "!" = phred score 0
     * "#" = phred score 2
     */

	std::string seq =  "AATTGGCCATGCCGTA";
	std::string qual = "HHHHHHHHHHHHHHHH";
	std::vector<bool> initial_depth = stringToBoolVect("1111111111111111");
	std::vector<Mutation> mutations{{7, 9, "A", "G", ""}};
	Read r(0, 15, seq);
	r.setQual(qual).setMutations(mutations);
	int exclude_length = 17;

	std::vector<Mutation> included_mutations;
	std::vector<bool> depth;

	std::string exp_depth = "0000000000000000";
	std::vector<Mutation> exp_included_mutations = {};

	boost::tie(included_mutations,
			   depth) = stripEnd(r,
								 initial_depth,
								 exclude_length,
								 RIGHT,
								 false);

	EXPECT_EQ(exp_depth, toString(depth));
	EXPECT_EQ(exp_included_mutations, included_mutations);
}

TEST(shiftAmbigIndels, gapWithMMShiftLeft){
/*
 * replacement seq:    CC
 *                     ^^^^^^
 *          target: TGCCGCGCGTGTA
 */
	std::string seq =  "TGCCGCGCGTGTA";
	std::string qual = "ABCDEFGHIJKLM";
	std::vector<Mutation> mutations{{2, 9, "CC", "#!"}};
	int left = 0;
	bool right_align_ambig_dels = false;
	bool right_align_ambig_ins = false;

	std::vector<Mutation> shifted_mutations;
	shifted_mutations = shiftAmbigIndels(mutations,
										 seq,
										 qual,
										 left,
										 right_align_ambig_dels,
										 right_align_ambig_ins);
	std::vector<Mutation> exp_mutations{{2, 7, "", "", "_ambig"},
										{7, 9, "C", "!", "_ambig"}};
	EXPECT_EQ(toString(exp_mutations), toString(shifted_mutations));
}

TEST(shiftAmbigIndels, gapWithMMShiftRight){
/*
 * replacement seq:    CC
 *                     ^^^^^^
 *          target: TGCCGCGCGTGTA
 */
	std::string seq =  "TGCCGCGCGTGTA";
	std::string qual = "ABCDEFGHIJKLM";
	std::vector<Mutation> mutations{{2, 9, "CC", "#!"}};
	int left = 0;
	bool right_align_ambig_dels = true;
	bool right_align_ambig_ins = true;

	std::vector<Mutation> shifted_mutations;
	shifted_mutations = shiftAmbigIndels(mutations,
										 seq,
										 qual,
										 left,
										 right_align_ambig_dels,
										 right_align_ambig_ins);
	std::vector<Mutation> exp_mutations{{3, 5, "C", "!", "_ambig"},
										{4, 9, "", "", "_ambig"}};
	EXPECT_EQ(toString(exp_mutations), toString(shifted_mutations));
}

TEST(shiftAmbigIndels, insertWithMMShiftLeft){
/*
 * replacement seq:    CGCGCG
 *                     ^^
 *          target: TGCCCTGTA
 */
	std::string seq =  "TGCCCTGTA";
	std::string qual = "ABCDEFGHI";
	std::vector<Mutation> mutations{{2, 5, "CGCGCG", "123456"}};
	int left = 0;
	bool right_align_ambig_dels = false;
	bool right_align_ambig_ins = false;

	std::vector<Mutation> shifted_mutations;
	shifted_mutations = shiftAmbigIndels(mutations,
										 seq,
										 qual,
										 left,
										 right_align_ambig_dels,
										 right_align_ambig_ins);
	std::vector<Mutation> exp_mutations{{2, 3, "CGCG", "1234", "_ambig"},
										{3, 5, "G", "6", "_ambig"}};
	EXPECT_EQ(toString(exp_mutations), toString(shifted_mutations));
}

TEST(shiftAmbigIndels, insertWithMMShiftRight){
/*
 * replacement seq:    CGCGCG
 *                     ^^
 *          target: TGCCCTGTA
 */
	std::string seq = "TGCCCTGTA";
	std::string qual = "ABCDEFGHI";
	std::vector <Mutation> mutations{{2, 5, "CGCGCG", "123456"}};
	int left = 0;
	bool right_align_ambig_dels = true;
	bool right_align_ambig_ins = true;

	std::vector <Mutation> shifted_mutations;
	shifted_mutations = shiftAmbigIndels(mutations,
										 seq,
										 qual,
										 left,
										 right_align_ambig_dels,
										 right_align_ambig_ins);
	std::vector <Mutation> exp_mutations{{3, 5, "G",    "2",    "_ambig"},
										 {4, 5, "CGCG", "3456", "_ambig"}};
	EXPECT_EQ(toString(exp_mutations), toString(shifted_mutations));
}

// Note: duplicated sequences are just stand-ins for Q-scores
TEST(Debug, ParseClassifyMutations) {
	std::string line = R"(M00236:2:000000000-A21YG:1:1106:15774:10066	0	136	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	32 33 "TCTTTC" "TCTTTC" "" 32 34 "T" "T" "" 82 84 "C" "C" "" 84 86 "A" "A" "" 114 118 "GA" "GA" "")";
	std::string line_segfault = R"(M00236:2:000000000-A21YG:1:1106:15774:10066	7	136	GGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	GGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	115 116 "C" "C" "")";
	std::string line_bad_insert = R"(M00236:2:000000000-A21YG:1:1106:15774:10066	7	136	GGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	GGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	25 29 "CCCC" "CCCC" "" 68 71 "G" "G" "")";
	std::string line_bad_insert2 = R"(M00236:2:000000000-A21YG:1:1106:15774:10066	11	136	CAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	CAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	46 50 "AAAA" "AAAA" "" 49 51 "A" "A" "" 88 90 "G" "G" "" 98 101 "C" "C" "" 109 111 "T" "T" "")";
	std::string line_crosses_primer = R"(M00236:2:000000000-A21YG:1:1106:15774:10066	0	136	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	32 34 "T" "T" "" 34 35 "T" "T" "" 68 70 "A" "A" "" 93 107 "" "" "" 108 112 "" "" "" 112 114 "A" "A" "" 115 117 "C" "C" "")";
	std::string line_substr_error = R"(M00236:2:000000000-A21YG:1:1106:15774:10066	1313	1447	CTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAAC	CTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAAC	1313 1317 "TGCTGCCTCCCGTAGGAGTCTGC" "TGCTGCCTCCCGTAGGAGTCTGC" "")";
	//std::cout << "\n################# Running ParseClassifyMutations test" << "\n" << std::flush;
	EXPECT_NO_THROW(processMutations_wrapper(line));
	EXPECT_NO_THROW(processMutations_wrapper(line_segfault));
	EXPECT_NO_THROW(processMutations_wrapper(line_bad_insert));
	EXPECT_NO_THROW(processMutations_wrapper(line_bad_insert2));
	EXPECT_NO_THROW(processMutations_wrapper(line_crosses_primer));
	EXPECT_NO_THROW(processMutations_wrapper(line_substr_error));
}

TEST(Debug, ParseClassifyMutationsSegfault2) {
	std::string line_segfault2 = R"(M00236:2:000000000-A21YG:1:1106:15774:10066	0	159	AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGA	AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGA	-1 1 "C" "C" "" 0 4 "AAACTTTTAAAT" "AAACTTTTAAAT" "")";
	EXPECT_NO_THROW(processMutations_wrapper(line_segfault2));
}

TEST(ProcessMutations, QualityFiltering) {
// "H" = phred score 39
// "!" = phred score 0
// "#" = phred score 2

	/*        q-scores:      H!
     *          insert:      CA
     *                      /
     *        q-scores: H HHHHHH#HHHHHHH
     *    aligned read: A-TTGGCCTTGCCGTA
     *          target: AATTGGCCATGCCGTA
     */
	std::string line = R"(M00236:dummy:QualityFiltering	0	15	AATTGGCCATGCCGTA	H!HHHHHH#HHHHHHH	0 2 "" "" "" 3 4 "CA" "H!" "" 7 9 "T" "#" "")";
	std::string included_mutations;
	std::string depth;
	std::string count;
	std::string exp_mutations;
	std::string exp_depth;
	std::string exp_count;

	// min_qual 0
	boost::tie(included_mutations,
			   depth,
			   count) = processMutations_wrapper_min_qual(line, 0);
	exp_mutations = R"(0 2 "" "" "A-" 3 4 "CA" "H!" "multinuc_insertion" 7 9 "T" "#" "AT")";
	exp_depth = "1111111111111111";
	exp_count = "0101000010000000";
	EXPECT_EQ(exp_depth, depth);
	EXPECT_EQ(exp_count, count);
	EXPECT_EQ(exp_mutations, included_mutations);

	// min_qual 2
	boost::tie(included_mutations,
			   depth,
			   count)  = processMutations_wrapper_min_qual(line, 2);
	exp_mutations = R"(0 2 "" "" "A-" 7 9 "T" "#" "AT")";
	exp_depth = "1111011111111111";
	exp_count = "0100000010000000";
	EXPECT_EQ(exp_depth, depth);
	EXPECT_EQ(exp_count, count);
	EXPECT_EQ(exp_mutations, included_mutations);

	// min_qual 40
	boost::tie(included_mutations,
			   depth,
			   count)  = processMutations_wrapper_min_qual(line, 40);
	exp_mutations = "";
	exp_depth = "0000000000000000";
	exp_count = "0000000000000000";
	EXPECT_EQ(exp_depth, depth);
	EXPECT_EQ(exp_count, count);
	EXPECT_EQ(exp_mutations, included_mutations);
}

TEST(ProcessMutations, QualityFilteringNeighbors) {
// "H" = phred score 39
// "!" = phred score 0
// "#" = phred score 2

	/*        q-scores:      HH
     *          insert:      CA
     *                      /
     *        q-scores: ! !!!HH#H#HHHHHH
     *    aligned read: A-TTGGCCTTGCCGTA
     *          target: AATTGGCCATGCCGTA
     */
	std::string line = R"(M00236:2:000000000-A21YG:1:1106:15774:10066	0	15	AATTGGCCATGCCGTA	!!!!!HH#H#HHHHHH	0 2 "" "" "" 3 4 "CA" "HH" "" 7 9 "T" "H" "")";
	std::string included_mutations;
	std::string depth;
	std::string count;
	std::string exp_mutations;
	std::string exp_depth;
	std::string exp_count;

	// min_qual 0
	//std::cout << "min_qual 0" << std::endl;
	boost::tie(included_mutations,
			   depth,
			   count) = processMutations_wrapper_min_qual(line, 0);
	//std::cout << "classified_mutations: " << included_mutations << std::endl;
	exp_mutations = R"(0 2 "" "" "A-" 3 4 "CA" "HH" "multinuc_insertion" 7 9 "T" "H" "AT")";
	exp_depth = "1111111111111111";
	exp_count = "0101000010000000";
	EXPECT_EQ(exp_mutations, included_mutations);
	EXPECT_EQ(exp_depth, depth);
	EXPECT_EQ(exp_count, count);

	// min_qual 2
	//std::cout << "min_qual 2" << std::endl;
	boost::tie(included_mutations,
			   depth,
			   count) = processMutations_wrapper_min_qual(line, 2);
	//std::cout << "classified_mutations: " << included_mutations << std::endl;
	exp_mutations = R"(7 9 "T" "H" "AT")";
	exp_depth = "0000001111111111";
	exp_count = "0000000010000000";
	EXPECT_EQ(exp_mutations, included_mutations);
	EXPECT_EQ(exp_depth, depth);
	EXPECT_EQ(exp_count, count);

	// min_qual 40
	//std::cout << "min_qual 40" << std::endl;
	boost::tie(included_mutations,
			   depth,
			   count) = processMutations_wrapper_min_qual(line, 40);
	//std::cout << "classified_mutations: " << included_mutations << std::endl;
	exp_mutations = "";
	exp_depth = "0000000000000000";
	exp_count = "0000000000000000";
	EXPECT_EQ(exp_mutations, included_mutations);
	EXPECT_EQ(exp_depth, depth);
	EXPECT_EQ(exp_count, count);
}

// problematic example
// - 4th and 5th mutations don't have the correct sequence
//   after shifting and merging
// - works now
TEST(Debug, ParseClassifyMutationsBug1) {
	std::vector<std::string> fields = {
			"M01228:25:000000000-A1CW0:1:1101:21515:5726",
			"16",
			"16S_crw",
			"344",
			"44",
			"58M1D1M1D46M2I2M1D91M",
			"*",
			"0",
			"0",
			"ACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCCAGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGTCAATGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGCAGGGTG",
			"??AAAABADDDDDDDDGGGGGGIIIIIIIIIGIIHHHHHHHHIIIHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJIIIIIIHHHHHHIIIIIHHHHHIHHIIIIHGGGGGGDDDDDDDEABAAA???",
			"AS:i:348",
			"XN:i:0",
			"XM:i:3",
			"XO:i:4",
			"XG:i:5",
			"NM:i:8",
			"MD:Z:58^G1^G0T47^G0G83G6",
			"YT:Z:UU"
	};

	Read read = mutation_parser::parseSamFields(fields, 10, true);
	EXPECT_NO_THROW(processMutations_wrapper(read.serializeForTest()));
}

// another problematic example
// - left-most 2 mutations have wrong sequence after
//   shifting and merging
// - 1st mutation appears correctly parsed, but oddly is a deletion including mismatches
// - shifting mutations doesn't correctly account for presence of internal mismatches
// - works now
TEST(Debug, ParseClassifyMutationsBug2) {
	std::vector<std::string> fields = {
			"M01228:25:000000000-A1CW0:1:1101:21515:5726",
			"16",
			"16S_crw",
			"344",
			"255",
			"57M2D48M1I94M",
			"*",
			"0",
			"0",
			"ACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCCAGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGTCAATGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGCAGGGTG",
			"??AAAABADDDDDDDDGGGGGGIIIIIIIIIGIIHHHHHHHHIIIHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJIIIIIIHHHHHHIIIIIHHHHHIHHIIIIHGGGGGGDDDDDDDEABAAA???",
			"MD:Z:57^CG1G0T45A1G0G83G6"
	};

	Read read = mutation_parser::parseSamFields(fields, 10, true);
	EXPECT_NO_THROW(processMutations_wrapper(read.serializeForTest()));
}

// similar situation to those above (ambig del with 2 deleted nucs)
// but this case is parsed correctly
TEST(Debug, ParseClassifyMutationsOkay) {
	std::vector<std::string> fields = {
			"M01228:25:000000000-A1CW0:1:1101:24421:5736",
			"16",
			"16S_crw",
			"304",
			"44",
			"60M2D84M",
			"*",
			"0",
			"0",
			"TGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCAGGGAG",
			"EEBFJJJHIJIIIJJIJJHHJJJJJJJJJJJJJHJJJJJJJJJJJJJFFFFDD=FFHHHDBF@@,@BGDDHHD<E7HHFCDE5DC+EBFC,CEEC>5,/A9?FFHHHHHHHFHEECCEEECCCFFCFFBBDDDDDD?BB?????",
			"AS:i:274",
			"XN:i:0",
			"XM:i:1",
			"XO:i:1",
			"XG:i:2",
			"NM:i:3",
			"MD:Z:60^AT78G5",
			"YT:Z:UU"
	};

	Read read = mutation_parser::parseSamFields(fields, 10, true);
	EXPECT_NO_THROW(processMutations_wrapper(read.serializeForTest()));
}

// merged mutation contains a spurious matching nucleotide on left
// - now fixed (collapseMutations() now strips matches from both ends)
TEST(Debug, ParseClassifyMutationsBug3) {
	std::vector<std::string> fields = {
			"M01228:25:000000000-A1CW0:1:1101:23724:5718",
			"16",
			"16S_crw",
			"1094",
			"44",
			"21M1D3M1I128M4S",
			"*",
			"0",
			"0",
			"GTCCCGCAACGAGCGCAACCCTATCCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCAGGGC",
			"???AA@JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJFBAAAA?",
			"AS:i:292",
			"XN:i:0",
			"XM:i:0",
			"XO:i:2",
			"XG:i:2",
			"NM:i:2",
			"MD:Z:21^T131",
			"YT:Z:UU"
	};

	Read read = mutation_parser::parseSamFields(fields, 10, true);
	EXPECT_NO_THROW(processMutations_wrapper(read.serializeForTest()));
}

// mutation near sequence end apparently results in substr() out-of-bounds at some point
// - working now (doing range check on left bound, and catching out_of_bounds exception to handle right bound)
TEST(Debug, ParseClassifyMutationsSubstrBug) {
	std::vector<std::string> fields = {
			"M01228:25:000000000-A1CW0:1:1101:17244:8581",
			"16",
			"16S_crw",
			"1314",
			"40",
			"1M20I134M",
			"*",
			"0",
			"0",
			"CTGCTGCCTCCCGTAGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAAC",
			"?????JHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJHBBB?????",
			"AS:i:-25"
					"XN:i:0",
			"XM:i:0",
			"XO:i:1",
			"XG:i:20",
			"NM:i:20",
			"MD:Z:135",
			"YT:Z:UU"
	};

	Read read = mutation_parser::parseSamFields(fields, 10, true);
	EXPECT_NO_THROW(processMutations_wrapper(read.serializeForTest()));
}

// ambiguous insert somehow turning into a deletion
// - error happens in collapseMutations()
// - working now (now making sure not to go outside mutation span when stripping matching nucs on either end)
TEST(Debug, InsToDelMisclassify) {
	std::vector<std::string> fields = {
			"M01228:25:000000000-A1CW0:1:1101:14887:6174",
			"16",
			"16S_crw",
			"35",
			"42",
			"76M1I26M1I17M",
			"*",
			"0",
			"0",
			"GCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGGTGAGTAATGTCTGGGAATCTGCCTTGATGGAGGTGGATAAC",
			"JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJIIIIIIIIIIIIIIIIIIIHGGGGGGEEEDDEEEBBBAA???",
			"AS:i:-16",
			"XN:i:0",
			"XM:i:2",
			"XO:i:2",
			"XG:i:2",
			"NM:i:4",
			"MD:Z:96A14G7",
			"YT:Z:UU"
	};

	Read read = mutation_parser::parseSamFields(fields, 10, true);
	EXPECT_NO_THROW(processMutations_wrapper(read.serializeForTest()));
}

TEST(Debug, shortRead){
	std::string mutation_line = R"(shortread	226	245	TCCTGGTAACGTTTTTATCC	@C,CC?FCA,8CF9FGDG<)";
	EXPECT_NO_THROW(processMutations_wrapper_exclude_3prime(mutation_line, 21));
}

TEST(FindClosestPrimers, A) {
	std::string primers = R"(>RNA-A
CTGGGACTTCCGAGGCAAC CATCACCTAGGAGGACGTACA
14 32 209 229
TGGGAAGGAGAGCGTCGTTA CAGTTCCAGGTGTCCTGCTT
147 166 336 355
GTCTGGTGGTGGGTCGTAAG GACAGTCGCTCCGTGACAG
419 438 593 611)";
	std::string r = R"([read]	UNPAIRED	418	611	-	INCLUDED	-999	AGTCTGGTGGTGGGTCGTAAGTTTAGGAGGTGACTGCATCCTCCAGCATCTCAACTCCGTCTGTCTACTGTGTGAGACTTCGGCGGACCATTAGGAATGAGATCCGTGAGATCCTTCCATCTTCTTGAAGTCGCCTTTAGGGTGGCTGCGAGGTAGAGGGTTGGGGGTTGGTGGGCTGTCACGGAGCGACTGTC	B=59DDFFFFFFFFFFFFEFFF<FFFFFFGGGGFFG?FAEFFFFGFCGGFCFFDGGEDGGFGGGGFGGGGGGFGFGFDGGGGGGGGGGFGGGGGGGFFGGGGF>GGFF<GGFFFFFGFF@GFGCGGGEGFGGFFAGGGGGGGGGFDGEDGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGG	11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111			)";
	int exp_fw_index = 2;
	int exp_rv_index = 2;

	std::string tmp_filename = getTestFileDir().string() + "/tmp/tmp.primers.txt";
	std::ofstream s(tmp_filename);
	s << primers << std::flush;
	s.close();
	std::vector<PrimerPair> primer_pairs = loadPrimerPairs(tmp_filename);
    // on buildbox, when /home/test/shapemapper*/bin/test_mutation_parser
    // executed from within /home/test, get failures for every test that
    // requires writing/reading temp file. Either a path/CWD issue or a 
    // permissions problem. Actually, probably just the fact that the tmp/
    // directory does not exist - need to create it before any attempt to open ofstream
    // here get "ERROR: Input file /home/shapemapper/cpp-src/test/files/tmp/tmp.primers.txt not found"    

	Read read;
	read = parseDebugRead(r);
	int fw_primer_index, rv_primer_index;
	boost::tie(fw_primer_index, rv_primer_index) =
			findClosestPrimers(read.left, read.right, primer_pairs, 10);
	EXPECT_EQ(exp_fw_index, fw_primer_index);
	EXPECT_EQ(exp_rv_index, rv_primer_index);
}

TEST(Debug, PairedReadMemoryCorruption) {
	std::string sample_reads = R"(M00236:570:000000000-BBRFH:1:1101:18209:3377:R1	99	RNA-B	1614	42	61S78M2D11M1I49M	=	1927	574	CCCAGTGTCTGTGTGTGGGAATTGGTATCTTGCACCCGTGGGAGTCGGGACATATAAATATCGACTTGGTTCTGTGTCTACTGTTAGAAGTTGAGTGGGAAGCTGCAGGCCCGCCAGGACCACTGGGTCCCTAGAGAGCAGGCCTGCGTGTTTCCCAAACAGCCGTTGCCCCCGTGGCCGAGGTAGGTAATCCATATTGG	CCCCCGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGDGGG6EGFFFGGGGGGDCGF	AS:i:255	XS:i:106	XN:i:0	XM:i:2	XO:i:2	XG:i:3	NM:i:5	MD:Z:72C5^AA40T19	YS:i:385	YT:Z:CP
M00236:570:000000000-BBRFH:1:1101:18209:3377:R2	147	RNA-B	1927	42	200M	=	1614	-574	CAAAGCGGAGGCCCGGGCTGCCCGCCGTCCCCCTTCAGCTCCCCACGGACTGTACCAGGCAGCTGGGCCCCGCAGGACAGGCCGCAGCGGGTGGCCGGCTCTGTCCCGTCTTGGGGTTCTTGGTGTCCACGTCTTGTGGGCCGTGGGCTTCTACCTGCCCTTGGCCTGCAGTGCTTTGCTGGAGAAGGGACTCCCTAGAC	A<?1:1???9>2>;??FFFB>>>95B>E>4DC6GFC;)?4EGD934CFGGFGF<CGFGGGGFFGGGGC6EGFGCEF9GEDGECGECGGGGGGEGF7GFFGGGGGGFEGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGCCCCC	AS:i:385	XN:i:0	XM:i:4	XO:i:0	XG:i:0	NM:i:4	MD:Z:37T60T97G1T1	YS:i:255	YT:Z:CP)";
	std::string test_dir = getTestFileDir().string();
	std::ofstream dummy(test_dir+"/tmp/tmp_prmc.sam");
	dummy << sample_reads;
	dummy.close();

	std::string primers = R"(>RNA-A
CTGGGACTTCCGAGGCAAC CATCACCTAGGAGGACGTACA
14 32 209 229
TGGGAAGGAGAGCGTCGTTA CAGTTCCAGGTGTCCTGCTT
147 166 336 355
GTCTGGTGGTGGGTCGTAAG GACAGTCGCTCCGTGACAG
419 438 593 611)";
	std::string tmp_filename = test_dir + "/tmp/tmp.primers.txt";
	std::ofstream s(tmp_filename);
	s << primers << std::flush;
	s.close();

	EXPECT_NO_THROW(
		parseSAM(test_dir+"/tmp/tmp_prmc.sam",
				 test_dir+"/tmp/tmp_prmc.mut",
				 "", // debug_outname
				 tmp_filename, // primers_filename
				 800, // max_paired_fragment_length
				 10, // min_mapq
				 false, // right_align_ambig_dels,
				 false, //right_align_ambig_ins,
				 5, // max_internal_match,
				 30, //min_qual,
				 0, // exclude_3prime,
				 "", //mutation_type,
				 false, // variant_mode,
				 true, // trim_primers
				 true, // require_forward_primer_mapped
				 true, // require_reverse_primer_mapped
				 10, // max_primer_offset
				 false, // input_is_unpaired
				 false, // debug,
				 true //warn_on_no_mapped
		);
	);
}

TEST(Debug, Segfault3) {
    std::string primers = R"(>RNA-C
AGGGGGCAAGAGGAATTACG AGTGGAATCCCACCCCCTAA
90 109 270 289
GGCGGCATGCGTTCCT CCAGGCTGGATGCAGTTAAGG
10 25 182 202
CCCTTAACTGCATCCAGCCT ACAATACAAAGCAATTTCCTCAGA
181 200 356 379
TCTGAGGGAGAACAAGACCGA GAAGCACCAAACGTGACCAT
145 165 329 348
CTGCACACCTACTAGTCACCA AGCCCCACAGAACTATTGTAAA
243 263 412 433)";
	
	std::string tmp_filename = getTestFileDir().string() + "/tmp/tmp.primers.txt";
	std::ofstream s(tmp_filename);
	s << primers << std::flush;
	s.close();
	//std::vector<PrimerPair> primer_pairs = loadPrimerPairs(tmp_filename);

	std::string sample_reads = R"(M00236:570:000000000-BBRFH:1:1102:20009:1455	153	RNA-C	360	36	24M2S	=	360	0	GAGGAAATTGCTTTGTATTGTACAGG	EFEFECE<<C,,9CF<E99F9C<CAA	AS:i:44	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:22T1	YT:Z:UP
)";
	std::string test_dir = getTestFileDir().string();
	std::ofstream dummy(test_dir+"/tmp/tmp_s3.sam");
	dummy << sample_reads;
	dummy.close();

	EXPECT_NO_THROW(
		parseSAM(test_dir+"/tmp/tmp_s3.sam",
				 test_dir+"/tmp/tmp_s3.mut",
				 "", // debug_outname
				 tmp_filename, // primers_filename
				 800, // max_paired_fragment_length
				 10, // min_mapq
				 false, // right_align_ambig_dels,
				 false, //right_align_ambig_ins,
				 5, // max_internal_match,
				 30, //min_qual,
				 0, // exclude_3prime,
				 "", //mutation_type,
				 false, // variant_mode,
				 true, // trim_primers
				 true, // require_forward_primer_mapped
				 true, // require_reverse_primer_mapped
				 10, // max_primer_offset
				 false, // input_is_unpaired
				 false, // debug,
				 true //warn_on_no_mapped
		);
	);
}


TEST(Debug, VectorOutOfRange) {
	std::string sample_reads = R"(K00270:121:HWNJGBBXX:6:1114:26666:6484	163	chr11:65505800-65505900	1	42	106S40M	=	65505813	265	GTTCAAGAACTGTAATGCTGGGTGGGAACATGTAACTTGTAGACTGGAGAAGTATTAAAACCACAGCTAAGTAGCTCTATTATAATACTTATCCAGAGACTAAAACAAACTTAAACCAGTCAGTGGAGAAATAACATGTTCAAGAA	A<-AAJFJF<AJFJF-AFJFFJ<AJ<FJJ-F-FAA--FFJJJF<JJ77-FAAAF<-7-<FF--77F-FFJ7--<--7A-AAA-7-77-77---7A<-AAJJ7FFFA--7----7----F-)----7AJ-----7--7--7-7---7	AS:i:186	XS:i:104	XN:i:0	XM:i:11	XO:i:1	XG:i:1	NM:i:12	MD:Z:C13A25	YS:i:294	YT:Z:CP
K00270:121:HWNJGBBXX:6:1106:2767:34565	147	chr11:65505800-65505900	1	44	77S73M	=	65505485	-388	CCTGCAAATTGTTAACAGAAGGGTATTAAAACCACAGCTAAGTAGCTCTATTATAATACTTATCCAGTGACTAAAACCAACTTAAACCAGTAAGTGGAGAAATAACATGTTCAAGAACTGTAATGCTGGGTGGGAACATGTAACTTGTAG	AJFJF<<<JJJJJJAJJJ<JFAJJJJFAFA7JJFJFJJA-AJJJFFJJJFJJ<F7-JFFJJ-AAFJJJJJJJF-AAJFJJJJJJJJJJFJFJAJJJFFJJJ-JJAJFJJJJJFFJFAJJJJFFJJJJJJJJJJFJFAJJJJJJJJFFFAA	AS:i:300	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:73	YS:i:300	YT:Z:CP)";

	std::string test_dir = getTestFileDir().string();
	std::ofstream dummy(test_dir+"/tmp/tmp_voor.sam");
	dummy << sample_reads;
	dummy.close();

	//EXPECT_NO_THROW(
		parseSAM(test_dir+"/tmp/tmp_voor.sam",
				 test_dir+"/tmp/tmp_voor.mut",
				 "", // debug_outname
				 "", // primers_filename
				 800, // max_paired_fragment_length
				 10, // min_mapq
				 false, // right_align_ambig_dels,
				 false, //right_align_ambig_ins,
				 5, // max_internal_match,
				 30, //min_qual,
				 0, // exclude_3prime,
				 "", //mutation_type,
				 false, // variant_mode,
				 true, // trim_primers
				 true, // require_forward_primer_mapped
				 true, // require_reverse_primer_mapped
				 10, // max_primer_offset
				 true, // input_is_unpaired
				 true, // debug,
				 true //warn_on_no_mapped
		);
	//);
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  // set location of top-level shapemapper path so we can locate data files
  if (argc > 1){
      BASEPATH = argv[1];
  }
  return RUN_ALL_TESTS();
}

// FIXME: seems like debug_out is shared/interacting across tests run in parallel. 
//         - created as a member of the mutation:: namespace.
//         - shouldn't be a real problem since debug output not used for anything in tests,
//           but means that to inspect debug output file, need to run individual test with
//           --gtest_filter=Debug.Whatever
// FIXME: maybe break up into multiple files
// FIXME: add tests for stripPrimers()
// FIXME: add tests for primer pair mapping filtering

// TODO: test reconstructed target sequences from a number of reads mapped to the same location are identical? cover more cases than just the handful here

// bunch of exception handling conditions?

// Test serialization/deserialization

// TODO: do some random stack stops in this module and find rate limiting steps
