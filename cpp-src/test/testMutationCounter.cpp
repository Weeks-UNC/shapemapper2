/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2017 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#include "gtest/gtest.h"
#include "MutationParser.cpp" // to load from SAM file or parse a SAM file line into mutations
                              // - not really testing this here, just for convenience
#include "MutationCounter.h"


//namespace BF = boost::filesystem;
using namespace mutation;
using namespace mutation_parser;
using namespace mutation_parser::detail;

// TODO: specific tests for collapseMutations, classifyMutations, VariantCounter




void updateCounts_wrapper(const std::string &line) {
    MutationCounter mc;
    std::string read_id;
    int left_target_pos;
    int right_target_pos;
    std::string local_target_seq;
    std::string local_target_qual;
    std::vector <Mutation> adjusted_mutations;
    boost::tie(read_id,
               left_target_pos,
               right_target_pos,
               local_target_seq,
               local_target_qual,
               adjusted_mutations) = parseReadInfo(line);
    //std::cout << toString(adjusted_mutations) << std::endl;
    mc.updateRightBound(right_target_pos);
    std::vector<TaggedMutation> classified_mutations;
    std::vector<bool> local_effective_depth;
    std::vector<bool> local_effective_count;
    boost::tie(classified_mutations,
               local_effective_depth,
               local_effective_count) = mc.updateCounts(read_id,
                                                        adjusted_mutations,
                                                        local_target_seq,
                                                        local_target_qual,
                                                        left_target_pos,
                                                        1,       // direction
                                                        false,    // separate_ambig_counts
                                                        false,   // right_align_ambig_dels
                                                        false,   // right_align_ambig_ins
                                                        6,      // max_internal_match
                                                        0,      // min_qual
                                                        1,       // exclude_3prime
                                                        "",      // mutation_type (empty = count all mutations)
                                                        false);  // print debug info

    //std::cout << toString(classified_mutations) << std::endl << std::endl;
}

void updateCounts_wrapper_exclude_3prime(const std::string &line,
                                         const int exclude_3prime) {
    MutationCounter mc;
    std::string read_id;
    int left_target_pos;
    int right_target_pos;
    std::string local_target_seq;
    std::string local_target_qual;
    std::vector <Mutation> adjusted_mutations;
    boost::tie(read_id,
               left_target_pos,
               right_target_pos,
               local_target_seq,
               local_target_qual,
               adjusted_mutations) = parseReadInfo(line);
    //std::cout << toString(adjusted_mutations) << std::endl;
    mc.updateRightBound(right_target_pos);
    std::vector<TaggedMutation> classified_mutations;
    std::vector<bool> local_effective_depth;
    std::vector<bool> local_effective_count;
    boost::tie(classified_mutations,
               local_effective_depth,
               local_effective_count) = mc.updateCounts(read_id,
                                                        adjusted_mutations,
                                                        local_target_seq,
                                                        local_target_qual,
                                                        left_target_pos,
                                                        1,       // direction
                                                        false,    // separate_ambig_counts
                                                        false,   // right_align_ambig_dels
                                                        false,   // right_align_ambig_ins
                                                        6,      // max_internal_match
                                                        0,      // min_qual
                                                        exclude_3prime,       // exclude_3prime
                                                        "",      // mutation_type (empty = count all mutations)
                                                        true);  // print debug info

    //std::cout << toString(classified_mutations) << std::endl << std::endl;
}

boost::tuple<std::string, std::string, std::string>
updateCounts_wrapper_min_qual(const std::string &line, 
                              const int min_qual) {
    MutationCounter mc;
    std::string read_id;
    int left_target_pos;
    int right_target_pos;
    std::string local_target_seq;
    std::string local_target_qual;
    std::vector <Mutation> adjusted_mutations;
    boost::tie(read_id,
               left_target_pos,
               right_target_pos,
               local_target_seq,
               local_target_qual,
               adjusted_mutations) = parseReadInfo(line);
    //std::cout << toString(adjusted_mutations) << std::endl;
    mc.updateRightBound(right_target_pos);
    std::vector<TaggedMutation> classified_mutations;
    std::vector<bool> local_effective_depth;
    std::vector<bool> local_effective_count;
    boost::tie(classified_mutations,
               local_effective_depth,
               local_effective_count) = mc.updateCounts(read_id,
                                                        adjusted_mutations,
                                                        local_target_seq,
                                                        local_target_qual,
                                                        left_target_pos,
                                                        1,       // direction
                                                        false,    // separate_ambig_counts
                                                        false,   // right_align_ambig_dels
                                                        false,   // right_align_ambig_ins
                                                        0,      // max_internal_match
                                                        min_qual,      // min_qual
                                                        1,       // exclude_3prime
                                                        "",      // mutation_type (empty = count all mutations)
                                                        false);  // print debug info
    //std::cout << toString(classified_mutations) << std::endl << std::endl;
    std::string cm = toString(classified_mutations);
    std::string led = toString(local_effective_depth);
    std::string lec = toString(local_effective_count);
    return boost::make_tuple(cm, led, lec);
}


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

    std::string local_target_seq =  "AATTGGCCATGCCGTA";
    std::string local_target_qual = "HHHHHHHHHHHHHHHH";
    std::vector<TaggedMutation> mutations{};
    int left_target_pos = 0;
    int exclude_3prime = 1;
    int min_qual = 30;

    std::string exp_local_effective_depth = "111111111111111";
    std::string exp_local_effective_count = "000000000000000";
    std::vector<TaggedMutation> exp_included_mutations = {};
    std::vector<TaggedMutation> exp_excluded_mutations = {};

    std::vector<bool> local_effective_depth;
    std::vector<bool> local_effective_count;
    std::vector<TaggedMutation> included_mutations;
    std::vector<TaggedMutation> excluded_mutations;
    boost::tie(local_effective_depth,
               local_effective_count,
               included_mutations,
               excluded_mutations) = filterQscoresCountDepths(mutations,
                                                              local_target_seq,
                                                              local_target_qual,
                                                              left_target_pos,
                                                              exclude_3prime,
                                                              min_qual,
                                                              "",
                                                              false); // variant_mode off
    EXPECT_EQ(exp_local_effective_depth, toString(local_effective_depth));
    EXPECT_EQ(exp_local_effective_count, toString(local_effective_count));
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

    std::string local_target_seq =  "AATTGGCCATGCCGTA";
    std::string local_target_qual = "HHHHHHHH!!HHHHHH";
    std::vector<TaggedMutation> mutations{{7, 10, "", "", ""}};
    int left_target_pos = 0;
    int exclude_3prime = 1;
    int min_qual = 30;

    std::string exp_local_effective_depth = "111111111111111";
    std::string exp_local_effective_count = "000000000100000";
    std::vector<TaggedMutation> exp_included_mutations = {};
    std::vector<TaggedMutation> exp_excluded_mutations = {};

    std::vector<bool> local_effective_depth;
    std::vector<bool> local_effective_count;
    std::vector<TaggedMutation> included_mutations;
    std::vector<TaggedMutation> excluded_mutations;
    boost::tie(local_effective_depth,
               local_effective_count,
               included_mutations,
               excluded_mutations) = filterQscoresCountDepths(mutations,
                                                              local_target_seq,
                                                              local_target_qual,
                                                              left_target_pos,
                                                              exclude_3prime,
                                                              min_qual,
                                                              "",
                                                              true); // variant_mode on
    EXPECT_EQ(exp_local_effective_depth, toString(local_effective_depth));
    EXPECT_EQ(exp_local_effective_count, toString(local_effective_count));
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

    std::string local_target_seq =  "AATTGGCCATGCCGTA";
    std::string local_target_qual = "HHHHHHHH!!HHHHHH";
    std::vector<TaggedMutation> mutations{{7, 10, "", "", ""}};
    int left_target_pos = 0;
    int exclude_3prime = 1;
    int min_qual = 30;

    std::string exp_local_effective_depth = "111111110111111";
    std::string exp_local_effective_count = "000000000100000";
    std::vector<TaggedMutation> exp_included_mutations = {};
    std::vector<TaggedMutation> exp_excluded_mutations = {};

    std::vector<bool> local_effective_depth;
    std::vector<bool> local_effective_count;
    std::vector<TaggedMutation> included_mutations;
    std::vector<TaggedMutation> excluded_mutations;
    boost::tie(local_effective_depth,
               local_effective_count,
               included_mutations,
               excluded_mutations) = filterQscoresCountDepths(mutations,
                                                              local_target_seq,
                                                              local_target_qual,
                                                              left_target_pos,
                                                              exclude_3prime,
                                                              min_qual,
                                                              "",
                                                              false); // variant_mode off
    EXPECT_EQ(exp_local_effective_depth, toString(local_effective_depth));
    EXPECT_EQ(exp_local_effective_count, toString(local_effective_count));
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

    std::string local_target_seq =  "AATTGGCCATGCCGTA";
    std::string local_target_qual = "HHHHHHHH!HHHHHHH";
    std::vector<TaggedMutation> mutations{{6, 8, "G", "H", ""},
                                          {7, 9, "", "", ""}};
    int left_target_pos = 0;
    int exclude_3prime = 1;
    int min_qual = 30;

    std::string exp_local_effective_depth = "111111111111111";
    std::string exp_local_effective_count = "000000011000000";
    std::vector<TaggedMutation> exp_included_mutations = {};
    std::vector<TaggedMutation> exp_excluded_mutations = {};

    std::vector<bool> local_effective_depth;
    std::vector<bool> local_effective_count;
    std::vector<TaggedMutation> included_mutations;
    std::vector<TaggedMutation> excluded_mutations;
    boost::tie(local_effective_depth,
               local_effective_count,
               included_mutations,
               excluded_mutations) = filterQscoresCountDepths(mutations,
                                                              local_target_seq,
                                                              local_target_qual,
                                                              left_target_pos,
                                                              exclude_3prime,
                                                              min_qual,
                                                              "",
                                                              true); // variant_mode on
    EXPECT_EQ(exp_local_effective_depth, toString(local_effective_depth));
    EXPECT_EQ(exp_local_effective_count, toString(local_effective_count));
}

TEST(filterQscoresCountDepths, largeExclude3prime){
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

    std::string local_target_seq =  "AATTGGCCATGCCGTA";
    std::string local_target_qual = "HHHHHHHHHHHHHHHH";
    std::vector<TaggedMutation> mutations{};
    int left_target_pos = 0;
    int exclude_3prime = 30;
    int min_qual = 30;

    std::string exp_local_effective_depth = "";
    std::string exp_local_effective_count = "";
    std::vector<TaggedMutation> exp_included_mutations = {};
    std::vector<TaggedMutation> exp_excluded_mutations = {};

    std::vector<bool> local_effective_depth;
    std::vector<bool> local_effective_count;
    std::vector<TaggedMutation> included_mutations;
    std::vector<TaggedMutation> excluded_mutations;
    boost::tie(local_effective_depth,
               local_effective_count,
               included_mutations,
               excluded_mutations) = filterQscoresCountDepths(mutations,
                                                              local_target_seq,
                                                              local_target_qual,
                                                              left_target_pos,
                                                              exclude_3prime,
                                                              min_qual,
                                                              "",
                                                              false); // variant_mode off
    EXPECT_EQ(exp_local_effective_depth, toString(local_effective_depth));
    EXPECT_EQ(exp_local_effective_count, toString(local_effective_count));
}

TEST(filterQscoresCountDepths, exclude3primeWithMutations){
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

    std::string local_target_seq =  "AATTGGCCATGCCGTA";
    std::string local_target_qual = "HHHHHHHHHHHHHHHH";
    std::vector<TaggedMutation> mutations{{7, 9, "A", "G", "H"}};
    int left_target_pos = 0;
    int exclude_3prime = 17;
    int min_qual = 30;

    std::string exp_local_effective_depth = "";
    std::string exp_local_effective_count = "";
    std::vector<TaggedMutation> exp_included_mutations = {};
    std::vector<TaggedMutation> exp_excluded_mutations = {};

    std::vector<bool> local_effective_depth;
    std::vector<bool> local_effective_count;
    std::vector<TaggedMutation> included_mutations;
    std::vector<TaggedMutation> excluded_mutations;
    boost::tie(local_effective_depth,
               local_effective_count,
               included_mutations,
               excluded_mutations) = filterQscoresCountDepths(mutations,
                                                              local_target_seq,
                                                              local_target_qual,
                                                              left_target_pos,
                                                              exclude_3prime,
                                                              min_qual,
                                                              "",
                                                              false); // variant_mode off
    EXPECT_EQ(exp_local_effective_depth, toString(local_effective_depth));
    EXPECT_EQ(exp_local_effective_count, toString(local_effective_count));
}


// Note: duplicated sequences are just stand-ins for Q-scores
TEST(Debug, ParseClassifyMutations) {
    std::string line = R"(M00236:2:000000000-A21YG:1:1106:15774:10066 0 136 GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC 32 33 "TCTTTC" "TCTTTC" 32 34 "T" "T" 82 84 "C" "C" 84 86 "A" "A" 114 118 "GA" "GA")";
    std::string line_segfault = R"(M00236:2:000000000-A21YG:1:1106:15774:10066 7 136 GGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC GGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC 115 116 "C" "C")";
    std::string line_bad_insert = R"(M00236:2:000000000-A21YG:1:1106:15774:10066 7 136 GGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC GGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC 25 29 "CCCC" "CCCC" 68 71 "G" "G")";
    std::string line_bad_insert2 = R"(M00236:2:000000000-A21YG:1:1106:15774:10066 11 136 CAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC CAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC 46 50 "AAAA" "AAAA" 49 51 "A" "A" 88 90 "G" "G" 98 101 "C" "C" 109 111 "T" "T")";
    std::string line_crosses_primer = R"(M00236:2:000000000-A21YG:1:1106:15774:10066 0 136 GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC 32 34 "T" "T" 34 35 "T" "T" 68 70 "A" "A" 93 107 "" "" 108 112 "" "" 112 114 "A" "A" 115 117 "C" "C")";
    std::string line_substr_error = R"(M00236:2:000000000-A21YG:1:1106:15774:10066 1313 1447 CTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAAC CTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAAC 1313 1317 "TGCTGCCTCCCGTAGGAGTCTGC" "TGCTGCCTCCCGTAGGAGTCTGC")";
    EXPECT_NO_THROW(updateCounts_wrapper(line));
    EXPECT_NO_THROW(updateCounts_wrapper(line_segfault));
    EXPECT_NO_THROW(updateCounts_wrapper(line_bad_insert));
    EXPECT_NO_THROW(updateCounts_wrapper(line_bad_insert2));
    EXPECT_NO_THROW(updateCounts_wrapper(line_crosses_primer));
    EXPECT_NO_THROW(updateCounts_wrapper(line_substr_error));
}


TEST(Debug, ParseClassifyMutationsSegfault2) {
    std::string line_segfault2 = R"(M00236:2:000000000-A21YG:1:1106:15774:10066 0 159 AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGA AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGA -1 1 "C" "C" 0 4 "AAACTTTTAAAT" "AAACTTTTAAAT")";
    EXPECT_NO_THROW(updateCounts_wrapper(line_segfault2));
}

TEST(UpdateCounts, QualityFiltering) {
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
    std::string line = R"(M00236:dummy:QualityFiltering 0 15 AATTGGCCATGCCGTA H!HHHHHH#HHHHHHH 0 2 "" "" 3 4 "CA" "H!" 7 9 "T" "#")";
    std::string included_mutations;
    std::string local_effective_depth;
    std::string local_effective_count;
    std::string exp_mutations;
    std::string exp_depth;
    std::string exp_count;

    // min_qual 0
    boost::tie(included_mutations,
               local_effective_depth,
               local_effective_count) = updateCounts_wrapper_min_qual(line, 0);
    exp_mutations =      R"(0 2 "" "" "A-" 3 4 "CA" "H!" "multinuc_insertion" 7 9 "T" "#" "AT")";
    exp_depth = "111111111111111";
    exp_count = "010100001000000";
    EXPECT_EQ(exp_depth, local_effective_depth);
    EXPECT_EQ(exp_count, local_effective_count);
    EXPECT_EQ(exp_mutations, included_mutations);
    
    // min_qual 2
    boost::tie(included_mutations,
               local_effective_depth,
               local_effective_count)  = updateCounts_wrapper_min_qual(line, 2);
    exp_mutations = R"(0 2 "" "" "A-" 7 9 "T" "#" "AT")";
    exp_depth = "111101111111111";
    exp_count = "010000001000000";
    EXPECT_EQ(exp_depth, local_effective_depth);
    EXPECT_EQ(exp_count, local_effective_count);
    EXPECT_EQ(exp_mutations, included_mutations);

    // min_qual 40
    boost::tie(included_mutations,
               local_effective_depth,
               local_effective_count)  = updateCounts_wrapper_min_qual(line, 40);
    exp_mutations = "";
    exp_depth = "000000000000000";
    exp_count = "000000000000000";
    EXPECT_EQ(exp_depth, local_effective_depth);
    EXPECT_EQ(exp_count, local_effective_count);
    EXPECT_EQ(exp_mutations, included_mutations);

}

TEST(UpdateCounts, QualityFilteringNeighbors) {
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
    std::string line = R"(M00236:2:000000000-A21YG:1:1106:15774:10066 0 15 AATTGGCCATGCCGTA !!!!!HH#H#HHHHHH 0 2 "" "" 3 4 "CA" "HH" 7 9 "T" "H")";
    std::string included_mutations;
    std::string local_effective_depth;
    std::string local_effective_count;
    std::string exp_mutations;
    std::string exp_depth;
    std::string exp_count;

    // min_qual 0
    //std::cout << "min_qual 0" << std::endl;
    boost::tie(included_mutations,
               local_effective_depth,
               local_effective_count) = updateCounts_wrapper_min_qual(line, 0);
    //std::cout << "classified_mutations: " << included_mutations << std::endl;
    exp_mutations = R"(0 2 "" "" "A-" 3 4 "CA" "HH" "multinuc_insertion" 7 9 "T" "H" "AT")";
    exp_depth = "111111111111111";
    exp_count = "010100001000000";
    EXPECT_EQ(exp_mutations, included_mutations);
    EXPECT_EQ(exp_depth, local_effective_depth);
    EXPECT_EQ(exp_count, local_effective_count);

    // min_qual 2
    //std::cout << "min_qual 2" << std::endl;
    boost::tie(included_mutations,
               local_effective_depth,
               local_effective_count) = updateCounts_wrapper_min_qual(line, 2);
    //std::cout << "classified_mutations: " << included_mutations << std::endl;
    exp_mutations = R"(7 9 "T" "H" "AT")";
    exp_depth = "000000111111111";
    exp_count = "000000001000000";
    EXPECT_EQ(exp_mutations, included_mutations);
    EXPECT_EQ(exp_depth, local_effective_depth);
    EXPECT_EQ(exp_count, local_effective_count);

    // min_qual 40
    //std::cout << "min_qual 40" << std::endl;
    boost::tie(included_mutations,
               local_effective_depth,
               local_effective_count) = updateCounts_wrapper_min_qual(line, 40);
    //std::cout << "classified_mutations: " << included_mutations << std::endl;
    exp_mutations = "";
    exp_depth = "000000000000000";
    exp_count = "000000000000000";
    EXPECT_EQ(exp_mutations, included_mutations);
    EXPECT_EQ(exp_depth, local_effective_depth);
    EXPECT_EQ(exp_count, local_effective_count);
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

    std::string read_id;
    int left_target_pos;
    int right_target_pos;
    std::string local_target_seq;
    std::string local_target_qual;
    std::vector<Mutation> adjusted_mutations;

    boost::tie(read_id,
               left_target_pos,
               right_target_pos,
               local_target_seq,
               local_target_qual,
               adjusted_mutations) = mutation_parser::parseSamMutations(fields);


    std::string mutation_line = serializeReadInfo(read_id,
                                                  left_target_pos,
                                                  right_target_pos,
                                                  local_target_seq,
                                                  local_target_qual,
                                                  adjusted_mutations);
    EXPECT_NO_THROW(updateCounts_wrapper(mutation_line));

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

    std::string read_id;
    int left_target_pos;
    int right_target_pos;
    std::string local_target_seq;
    std::string local_target_qual;
    std::vector<Mutation> adjusted_mutations;

    boost::tie(read_id,
               left_target_pos,
               right_target_pos,
               local_target_seq,
               local_target_qual,
               adjusted_mutations) = mutation_parser::parseSamMutations(fields);


    std::string mutation_line = serializeReadInfo(read_id,
                                                  left_target_pos,
                                                  right_target_pos,
                                                  local_target_seq,
                                                  local_target_qual,
                                                  adjusted_mutations);
    EXPECT_NO_THROW(updateCounts_wrapper(mutation_line));

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

    std::string read_id;
    int left_target_pos;
    int right_target_pos;
    std::string local_target_seq;
    std::string local_target_qual;
    std::vector<Mutation> adjusted_mutations;

    boost::tie(read_id,
               left_target_pos,
               right_target_pos,
               local_target_seq,
               local_target_qual,
               adjusted_mutations) = mutation_parser::parseSamMutations(fields);


    std::string mutation_line = serializeReadInfo(read_id,
                                                  left_target_pos,
                                                  right_target_pos,
                                                  local_target_seq,
                                                  local_target_qual,
                                                  adjusted_mutations);
    EXPECT_NO_THROW(updateCounts_wrapper(mutation_line));

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

    std::string read_id;
    int left_target_pos;
    int right_target_pos;
    std::string local_target_seq;
    std::string local_target_qual;
    std::vector<Mutation> adjusted_mutations;

    boost::tie(read_id,
               left_target_pos,
               right_target_pos,
               local_target_seq,
               local_target_qual,
               adjusted_mutations) = mutation_parser::parseSamMutations(fields);


    std::string mutation_line = serializeReadInfo(read_id,
                                                  left_target_pos,
                                                  right_target_pos,
                                                  local_target_seq,
                                                  local_target_qual,
                                                  adjusted_mutations);
    EXPECT_NO_THROW(updateCounts_wrapper(mutation_line));
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

    std::string read_id;
    int left_target_pos;
    int right_target_pos;
    std::string local_target_seq;
    std::string local_target_qual;
    std::vector<Mutation> adjusted_mutations;

    boost::tie(read_id,
               left_target_pos,
               right_target_pos,
               local_target_seq,
               local_target_qual,
               adjusted_mutations) = mutation_parser::parseSamMutations(fields);


    std::string mutation_line = serializeReadInfo(read_id,
                                                  left_target_pos,
                                                  right_target_pos,
                                                  local_target_seq,
                                                  local_target_qual,
                                                  adjusted_mutations);
    EXPECT_NO_THROW(updateCounts_wrapper(mutation_line));

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

    std::string read_id;
    int left_target_pos;
    int right_target_pos;
    std::string local_target_seq;
    std::string local_target_qual;
    std::vector<Mutation> adjusted_mutations;

    boost::tie(read_id,
               left_target_pos,
               right_target_pos,
               local_target_seq,
               local_target_qual,
               adjusted_mutations) = mutation_parser::parseSamMutations(fields);


    std::string mutation_line = serializeReadInfo(read_id,
                                                  left_target_pos,
                                                  right_target_pos,
                                                  local_target_seq,
                                                  local_target_qual,
                                                  adjusted_mutations);
    EXPECT_NO_THROW(updateCounts_wrapper(mutation_line));

}

TEST(exclude3prime, shortRead){
    std::string mutation_line = "shortread 226 245 TCCTGGTAACGTTTTTATCC @C,CC?FCA,8CF9FGDG<;";
    /*--exclude_3prime 21
      --max_internal_match 5
      --min_qual 30
      --length 631 */
    EXPECT_NO_THROW(updateCounts_wrapper_exclude_3prime(mutation_line, 21));
}


TEST(shiftAmbigIndels, gapWithMMShiftLeft){
/*
 * replacement seq:    CC
 *                     ^^^^^^
 *          target: TGCCGCGCGTGTA
 */
    std::string local_target_seq =  "TGCCGCGCGTGTA";
    std::string local_target_qual = "ABCDEFGHIJKLM";
    std::vector<Mutation> mutations{{2, 9, "CC", "#!"}};
    int left_target_pos = 0;
    bool right_align_ambig_dels = false;
    bool right_align_ambig_ins = false;

    std::vector<TaggedMutation> shifted_mutations;
    shifted_mutations = shiftAmbigIndels(mutations,
                                         local_target_seq,
                                         local_target_qual,
                                         left_target_pos,
                                         right_align_ambig_dels,
                                         right_align_ambig_ins);
    std::vector<TaggedMutation> exp_mutations{{2, 7, "", "", "_ambig"},
                                              {7, 9, "C", "!", "_ambig"}};
    EXPECT_EQ(toString(exp_mutations), toString(shifted_mutations));
}


TEST(shiftAmbigIndels, gapWithMMShiftRight){
/*
 * replacement seq:    CC
 *                     ^^^^^^
 *          target: TGCCGCGCGTGTA
 */
    std::string local_target_seq =  "TGCCGCGCGTGTA";
    std::string local_target_qual = "ABCDEFGHIJKLM";
    std::vector<Mutation> mutations{{2, 9, "CC", "#!"}};
    int left_target_pos = 0;
    bool right_align_ambig_dels = true;
    bool right_align_ambig_ins = true;

    std::vector<TaggedMutation> shifted_mutations;
    shifted_mutations = shiftAmbigIndels(mutations,
                                         local_target_seq,
                                         local_target_qual,
                                         left_target_pos,
                                         right_align_ambig_dels,
                                         right_align_ambig_ins);
    std::vector<TaggedMutation> exp_mutations{{3, 5, "C", "!", "_ambig"},
                                              {4, 9, "", "", "_ambig"}};
    EXPECT_EQ(toString(exp_mutations), toString(shifted_mutations));
}

TEST(shiftAmbigIndels, insertWithMMShiftLeft){
/*
 * replacement seq:    CGCGCG
 *                     ^^
 *          target: TGCCCTGTA
 */
    std::string local_target_seq =  "TGCCCTGTA";
    std::string local_target_qual = "ABCDEFGHI";
    std::vector<Mutation> mutations{{2, 5, "CGCGCG", "123456"}};
    int left_target_pos = 0;
    bool right_align_ambig_dels = false;
    bool right_align_ambig_ins = false;

    std::vector<TaggedMutation> shifted_mutations;
    shifted_mutations = shiftAmbigIndels(mutations,
                                         local_target_seq,
                                         local_target_qual,
                                         left_target_pos,
                                         right_align_ambig_dels,
                                         right_align_ambig_ins);
    std::vector<TaggedMutation> exp_mutations{{2, 3, "CGCG", "1234", "_ambig"},
                                              {3, 5, "G", "6", "_ambig"}};
    EXPECT_EQ(toString(exp_mutations), toString(shifted_mutations));
}

TEST(shiftAmbigIndels, insertWithMMShiftRight){
/*
 * replacement seq:    CGCGCG
 *                     ^^
 *          target: TGCCCTGTA
 */
    std::string local_target_seq = "TGCCCTGTA";
    std::string local_target_qual = "ABCDEFGHI";
    std::vector <Mutation> mutations{{2, 5, "CGCGCG", "123456"}};
    int left_target_pos = 0;
    bool right_align_ambig_dels = true;
    bool right_align_ambig_ins = true;

    std::vector <TaggedMutation> shifted_mutations;
    shifted_mutations = shiftAmbigIndels(mutations,
                                         local_target_seq,
                                         local_target_qual,
                                         left_target_pos,
                                         right_align_ambig_dels,
                                         right_align_ambig_ins);
    std::vector <TaggedMutation> exp_mutations{{3, 5, "G",    "2",    "_ambig"},
                                               {4, 5, "CGCG", "3456", "_ambig"}};
    EXPECT_EQ(toString(exp_mutations), toString(shifted_mutations));
}



// TODO: move these Histogram tests to a separate file

TEST(Histogram, linearSimple){
    Histogram mut_per_read("Mutations per read",1,5,5);
    mut_per_read.count(2);
    mut_per_read.count(2);
    mut_per_read.count(4);
    mut_per_read.count(5);
    mut_per_read.count(7);
    EXPECT_EQ(mut_per_read.total_reads, 5);
    EXPECT_EQ(mut_per_read.printCountsRow(), "0\t2\t0\t1\t2");
}

TEST(Histogram, linearSimpleTable){
    Histogram mut_per_read("Mutations per read",0,2,3);
    mut_per_read.count(1);
    mut_per_read.count(1);
    mut_per_read.count(1);
    mut_per_read.count(4);
    std::string exp = "Mutations per read\n"
                      "--------------------\n"
                      "bin_left\tfrequency\n"
                      "0\t0.000000\n"
                      "1\t0.750000\n"
                      "2\t0.250000\n"
                      "--------------------\n";
    EXPECT_EQ(mut_per_read.total_reads, 4);
    EXPECT_EQ(mut_per_read.printFreqTable(), exp);
}

TEST(Histogram, linearBiggerBins) {
    Histogram mut_per_read("Mutations per read",1, 10, 5);
    mut_per_read.count(4);
    mut_per_read.count(4);
    mut_per_read.count(8);
    mut_per_read.count(10);
    mut_per_read.count(14);
    EXPECT_EQ(mut_per_read.total_reads, 5);
    EXPECT_EQ(mut_per_read.printCountsRow(), "0\t2\t0\t1\t2");
}
