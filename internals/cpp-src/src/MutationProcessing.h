/** @file
 * @brief Mutation class and several related methods. 
 *        Used by MutationParser and MutationCounter.
 */

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2018 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#ifndef SHAPEMAPPER_MUTATIONPROCESSING_H
#define SHAPEMAPPER_MUTATIONPROCESSING_H

#include <iostream>
#include <fstream>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/newline.hpp>

#include "Mutation.h"
#include "Read.h"
#include "PrimerPair.h"
#include "util.h"

namespace BF = boost::filesystem;
namespace BI = boost::iostreams;

// FIXME: clean up/unify various std::to_string, mutation::toString, util::toString funcs
// FIXME: many (all?) of these functions should be updated to use the Read class
//        for input and output

namespace mutation {
    std::ofstream debug_out;

    /**
     * @brief Return a copied vector of mutations with ambiguously aligned
     *        indels realigned left or right.
     */
    std::vector<Mutation> shiftAmbigIndels(const std::vector<Mutation> &mutations,
                                           const std::string &local_target_seq,
                                           const std::string &local_target_qual,
                                           const int left_target_pos,
                                           const bool right_align_ambig_dels,
                                           const bool right_align_ambig_ins) {

        std::vector<Mutation> adjusted_mutations;
        // FIXME: replace iterators with auto
        if (mutations.size() > 0) {
            for (std::vector<Mutation>::const_iterator it = mutations.begin();
                 it != mutations.end();
                 ++it) {
                int left = (*it).left;
                int right = (*it).right;
                std::string seq = (*it).seq;
                std::string qual = (*it).qual;
                std::string tag = "";
                bool ambig = true;
                if ((*it).isAmbiguous()) {
                    // note: a bit ugly as a result of handling multinuc ambiguous mutations
                    //      with "internal" mismatches. checking for mismatches and
                    //      pushing back additional mutations if needed after shifting.
                    if ((*it).isGap()) {
                        if (right_align_ambig_dels) {
                            left = (*it).left + (*it).seq.length();
                            seq = "";
                            qual = "";
                            // check for mismatches and create mutations if needed
                            for (int i = 0; i < (*it).seq.length(); ++i) {
                                int local_left = (*it).left - left_target_pos + 1 + i;
                                std::string c = (*it).seq.substr(i, 1);
                                std::string q = (*it).qual.substr(i, 1);
                                //std::cout << "c: " << c << std::endl;
                                if (c != local_target_seq.substr(local_left, 1)) {
                                    Mutation m((*it).left + i,
                                               (*it).left + i + 2,
                                               c, q, tag, ambig);
                                    adjusted_mutations.push_back(m);
                                }
                            }
                            Mutation m(left, right, seq, qual, tag, ambig);
                            adjusted_mutations.push_back(m);
                        } else {
                            right = (*it).right - (*it).seq.length();
                            seq = "";
                            qual = "";
                            Mutation m(left, right, seq, qual, tag, ambig);
                            adjusted_mutations.push_back(m);
                            // check for mismatches and create mutations if needed
                            for (int i = 0; i < (*it).seq.length(); ++i) {
                                int local_left = right - left_target_pos + i;
                                std::string c = (*it).seq.substr(i, 1);
                                std::string q = (*it).qual.substr(i, 1);
                                if (c != local_target_seq.substr(local_left, 1)) {
                                    Mutation m(right + i - 1,
                                               right + i + 1,
                                               c, q, tag, ambig);
                                    adjusted_mutations.push_back(m);
                                }
                            }
                        }
                    } else if ((*it).isInsert()) {
                        if (right_align_ambig_ins) {
                            int d = (*it).seq.length() - ((*it).right - (*it).left - 1);
                            left = (*it).right - 1;
                            // take right-most portion of seq
                            seq = (*it).seq.substr((*it).seq.length() - d);
                            qual = (*it).qual.substr((*it).qual.length() - d);
                            // check for mismatches and create mutations if needed
                            for (int i = 0; i < (*it).seq.length() - d; ++i) {
                                int local_left = (*it).left - left_target_pos + 1 + i;
                                std::string c = (*it).seq.substr(i, 1);
                                std::string q = (*it).qual.substr(i, 1);
                                if (c != local_target_seq.substr(local_left, 1)) {
                                    Mutation m((*it).left + i,
                                               (*it).left + i + 2,
                                               c, q, tag, ambig);
                                    adjusted_mutations.push_back(m);
                                }
                            }
                            Mutation m(left, right, seq, qual, tag, ambig);
                            adjusted_mutations.push_back(m);
                        } else {
                            int d = (*it).seq.length() - ((*it).right - (*it).left - 1);
                            right = (*it).left + 1;
                            // take left-most portion of seq
                            seq = (*it).seq.substr(0, d);
                            qual = (*it).qual.substr(0, d);
                            Mutation m(left, right, seq, qual, tag, ambig);
                            adjusted_mutations.push_back(m);
                            // check for mismatches and create mutations if needed
                            for (int i = 0; i < (*it).seq.length() - d; ++i) {
                                int local_left = (*it).left - left_target_pos + 1 + i;
                                std::string c = (*it).seq.substr(d + i, 1);
                                std::string q = (*it).qual.substr(d + i, 1);
                                if (c != local_target_seq.substr(local_left, 1)) {
                                    Mutation m((*it).left + i,
                                               (*it).left + i + 2,
                                               c, q, tag, ambig);
                                    adjusted_mutations.push_back(m);
                                }
                            }
                        }

                    }
                } else {
                    ambig = false;
                    Mutation m(left, right, seq, qual, tag, ambig);
                    adjusted_mutations.push_back(m);
                }


            }
        }
        return adjusted_mutations;
    }


    /**
     * @brief Return a copied vector of mutations with 3-prime mutations within
     *        exclude_3prime nucleotides removed.
     */
    boost::tuple<std::vector<Mutation>,
            std::vector<bool>>
    stripEnd(const Read &r,
             const std::vector<bool> &effective_depth_in,
             const int exclude_length,
             const int which_end,
             const bool debug) {
        if (debug) {
            std::cout << "exclude_length: " << exclude_length << "\n";
        }
        std::vector<Mutation> stripped_mutations;
        std::vector<bool> effective_depth;
        for (auto d : effective_depth_in) {
            effective_depth.push_back(d);
        }
        for (int i = 0; i < effective_depth.size(); i++) {
            if (i + exclude_length >= effective_depth.size()) {
                effective_depth[i] = 0;
            }
        }
        if (r.mutations.size() > 0) {
            if (which_end == RIGHT) {
                int max_right = r.left + r.seq.length() - exclude_length - 1;
                for (auto &m : r.mutations) {
                    if ((m.right - 1) <= max_right) {
                        stripped_mutations.push_back(Mutation(m));
                    }
                }
            } else if (which_end == LEFT) {
                int min_left = r.left + exclude_length;
                for (auto &m : r.mutations) {
                    if ((m.left + 1) >= min_left) {
                        stripped_mutations.push_back(Mutation(m));
                    }
                }
            }
        }
        return boost::make_tuple(stripped_mutations,
                                 effective_depth);
    }

    /**
     * @brief
     */
    boost::tuple<std::vector<Mutation>,
            std::vector<bool>>
    stripPrimers(const std::vector<Mutation> &mutations,
                 const int left,
                 const std::vector<bool> &depth_in,
                 const PrimerPair &primer_pair,
                 const bool debug) {
        if (debug) { std::cout << primer_pair.toString() << std::endl << std::flush; }

        int len = depth_in.size();
        int right = left + len - 1;
        std::vector<bool> depth;
        for (auto d : depth_in) {
            depth.push_back(d);
        }

        std::vector<Mutation> stripped_mutations;
        if (primer_pair.fw_left > -1) {
            if (primer_pair.fw_right >= left) {
                for (int i = 0;
                     i + left <= primer_pair.fw_right;
                     i++) {
                    try {
                        depth.at(i) = 0;
                    } catch (std::out_of_range &e) { }
                }
            }
            if (primer_pair.rv_left <= right) {
                for (int i = len - 1;
                     i + left >= primer_pair.rv_left;
                     i--) {
                    try {
                        depth.at(i) = 0;
                    } catch (std::out_of_range &e) { }
                }
            }
            for (auto &m : mutations) {
                if (m.right < primer_pair.rv_left and m.left > primer_pair.fw_right) {
                    stripped_mutations.push_back(Mutation(m));
                }
            }
        } else {
            for (auto &m : mutations) {
                stripped_mutations.push_back(Mutation(m));
            }
        }
        return boost::make_tuple(stripped_mutations,
                                 depth);
    }

    /**
     * @brief Return a copied vector of mutations with adjacent mutations combined.
     *        Assumes mutations are sorted left-to-right. Mutations with up to
     *        max_internal_match unchanged nucleotides between them will be combined.
     */
    std::vector<Mutation> collapseMutations(const std::vector<Mutation> &mutations,
                                            const int max_internal_match,
                                            const std::string &local_target_seq,
                                            const std::string &local_target_qual,
                                            const int left_target_pos) {
        std::vector<Mutation> collapsed_mutations;
        std::vector<Mutation> unmerged_mutations; // temp storage for ambiguous mutations (if any),
        // to avoid merging with other mutations
        if (mutations.size() > 0) {
            for (std::vector<Mutation>::const_iterator it = mutations.begin();
                 it != mutations.end();
                 ++it) {
                Mutation m(*it);
                // if a mutation is not really a mutation (N match),
                // it should not be merged with other mutations
                if (m.tag == "N_match") {
                    unmerged_mutations.push_back(m);
                } else if (collapsed_mutations.size() > 0 and
                           m.left - (collapsed_mutations.back().right - 1) <= max_internal_match) {
                    // mutation within allowed distance for adjacency
                    std::string seq_sub = local_target_seq.substr(collapsed_mutations.back().right - left_target_pos,
                                                                  m.left - collapsed_mutations.back().right + 1);
                    // ^ pick up nucs between mutations
                    if (seq_sub.find('_') != std::string::npos) {
                        // don't collapse mutations across the gap between R1 and R2 if present
                        // (indicated by '_' chars in sequence)
                        collapsed_mutations.push_back(m);
                    } else {
                        collapsed_mutations.back().seq += seq_sub;
                        collapsed_mutations.back().qual += local_target_qual.substr(
                                collapsed_mutations.back().right - left_target_pos,
                                m.left - collapsed_mutations.back().right + 1);
                        collapsed_mutations.back().right = m.right;
                        collapsed_mutations.back().seq += m.seq;
                        collapsed_mutations.back().qual += m.qual;
                        collapsed_mutations.back().tag = "";
                        if (m.ambig) {
                            collapsed_mutations.back().ambig = true;
                        }
                    }
                } else {
                    // non-adjacent mutation
                    collapsed_mutations.push_back(m);
                }
            }
            // iteratively strip any matches from mutation ends (these can appear as a result of
            // shifted indels). For left-aligned ambig gaps, this only affects the 5-prime end of
            // a mutation, so shouldn't affect reactivity profiles. Adjusting it anyway for
            // consistency between aligners.
            if (collapsed_mutations.size() > 0) {
                for (std::vector<Mutation>::iterator it = collapsed_mutations.begin();
                     it != collapsed_mutations.end();
                     ++it) {
                    // from left end
                    int new_left = (*it).left;
                    std::string new_seq = (*it).seq;
                    std::string new_qual = (*it).qual;
                    try {
                        for (int i = 0; i < (*it).seq.length(); ++i) {
                            std::string c = (*it).seq.substr(i, 1);
                            std::string q = (*it).qual.substr(i, 1);
                            if ((*it).left + 1 + i >=
                                (*it).right) { break; } // don't go outside of region mutation spans
                            int p = (*it).left + 1 + i - left_target_pos;
                            if (p < 0) { break; }
                            std::string r = local_target_seq.substr(p, 1);
                            if (c == r) {
                                ++new_left;
                                new_seq = new_seq.substr(1);
                                new_qual = new_qual.substr(1);
                            } else {
                                break;
                            }
                        }
                    }
                    catch (std::out_of_range &exception) {
                        break;
                    }
                    (*it).left = new_left;
                    (*it).seq = new_seq;
                    (*it).qual = new_qual;

                    // from right end
                    int new_right = (*it).right;
                    new_seq = (*it).seq;
                    new_qual = (*it).qual;
                    try {
                        for (int i = (*it).seq.length() - 1; i >= 0; --i) {
                            std::string c = (*it).seq.substr(i, 1);
                            std::string q = (*it).qual.substr(i, 1);
                            int d = (*it).seq.length() - i;
                            if ((*it).right - d <= (*it).left) { break; } // don't go outside of region mutation spans
                            int p = (*it).right - d - left_target_pos;
                            if (p < 0) { break; }
                            std::string r = local_target_seq.substr(p, 1);
                            if (c == r) {
                                --new_right;
                                new_seq = new_seq.substr(0, i);
                                new_qual = new_qual.substr(0, i);
                            } else {
                                break;
                            }
                        }
                    }
                    catch (std::out_of_range &exception) {
                        break;
                    }
                    (*it).right = new_right;
                    (*it).seq = new_seq;
                    (*it).qual = new_qual;
                }
            }

            // put unmerged (mismatched Ns) mutations back into list
            collapsed_mutations.insert(collapsed_mutations.end(),
                                       unmerged_mutations.begin(),
                                       unmerged_mutations.end());
            // make sure unmerged are back in order
            std::sort(collapsed_mutations.begin(), collapsed_mutations.end());
        }
        return collapsed_mutations;
    }

    /**
     * @brief Assign a mutation classification to each mutation in vector.
     *        Return a copy with new tags.
     */
    std::vector<Mutation> classifyMutations(std::vector<Mutation> &mutations,
                                            const std::string &local_target_seq,
                                            const std::string &local_target_qual,
                                            const int target_pos) {
        std::vector<Mutation> classified_mutations;
        if (mutations.size() > 0) {
            for (std::vector<Mutation>::iterator it = mutations.begin();
                 it != mutations.end();
                 ++it) {
                Mutation m((*it));
                if (m.tag.length() == 0) {
                    m.tag = m.classify(local_target_seq,
                                       target_pos);
                }
                classified_mutations.push_back(m);
            }
        }
        return classified_mutations;
    }

    /**
     * @brief Filter mutations containing or neighboring bad basecalls.
     *        Filter non-mutation positions for qscore-aware read depth calculation.
     *
     *        In variant mode, mutation spans will not be excluded from depth, and
     *        gaps will contribute to depth if the basecalls on either side are good.
     */
    boost::tuple<std::vector<bool>,
            std::vector<bool>,
            std::vector<Mutation>,
            std::vector<Mutation>>
    filterQscoresCountDepths(const std::vector<Mutation> &mutations,
                             const std::string &seq,
                             const std::string &qual,
                             const std::vector<bool> &effective_depth_in,
                             const int left,

                             const int min_qual,
                             const std::string &mutation_type,
                             const bool variant_mode) {

        std::vector<Mutation> included_mutations;
        std::vector<Mutation> excluded_mutations;
        int len = seq.length();
        len = std::max(len, 0);
        //std::cout << "len: " << len << std::endl;
        std::vector<bool> effective_depth;
        for (auto d : effective_depth_in) {
            effective_depth.push_back(d);
        }

        // store mutation vector indices by left and right unchanged position
        // so we don't have to iterate over all mutations for each target position
        std::vector<int> left_mut_indices(len, -1);
        std::vector<int> right_mut_indices(len, -1);
        std::vector<bool> in_mutation(len, 0);
        for (int i = 0; i < mutations.size(); ++i) {
            Mutation m(mutations[i]);
            try {
                left_mut_indices.at(m.left - left) = i;
            } catch (std::out_of_range &exception) { }
            try {
                right_mut_indices.at(m.right - left) = i;
            } catch (std::out_of_range &exception) { }
            for (int n = m.left + 1 - left;
                 n < m.right - left;
                 ++n) {
                try {
                    in_mutation[n] = 1;
                } catch (std::out_of_range &exception) { }
            }
        }

        // go one basecall at a time, examining neighboring basecalls, and include or exclude from counts and depths
        // in first pass, skip over nucs in mutations (but do consider them as neighbors)
        // NOTE: '~' chars in qual are used to indicate unmapped positions between
        //       mate pair reads R1 and R2
        for (int i = 0;
             i < len;
             ++i) {
            if (in_mutation[i]) {
                continue;
            }
            bool bad_basecall = false;
            // look at self
            if (qual[i] - 33 < min_qual or qual[i] == '~') {
                bad_basecall = true;
            }
            // look at neighbor on left
            if (not bad_basecall) {
                if (right_mut_indices[i] != -1) {
                    Mutation m(mutations[right_mut_indices[i]]);
                    if (m.seq.length() > 0) {
                        // look at right-most basecall in mutation
                        // - should this include all basecalls in mutation?
                        if (m.qual.back() - 33 < min_qual) {
                            bad_basecall = true;
                        }
                    } else {
                        // look at basecall on other side of gap
                        int n = m.left - left;
                        try {
                            if (qual.at(n) - 33 < min_qual) {
                                bad_basecall = true;
                            }
                        } catch (std::out_of_range &exception) { }
                    }
                } else {
                    // look at basecall on left (not in a mutation)
                    try {
                        if (qual.at(i - 1) - 33 < min_qual) {
                            bad_basecall = true;
                        }
                    } catch (std::out_of_range &exception) { }
                }
            }
            // look at neighbor on right
            if (not bad_basecall) {
                if (left_mut_indices[i] != -1) {
                    Mutation m(mutations[left_mut_indices[i]]);
                    if (m.seq.length() > 0) {
                        // look at left-most basecall in mutation
                        // - should this include all basecalls in mutation?
                        if (m.qual.front() - 33 < min_qual) {
                            bad_basecall = true;
                        }
                    } else {
                        // look at basecall on other side of gap
                        int n = m.right - left;
                        try {
                            if (qual.at(n) - 33 < min_qual) {
                                bad_basecall = true;
                            }
                        } catch (std::out_of_range &exception) { }
                    }
                } else {
                    // look at basecall on right (not in a mutation)
                    try {
                        if (qual.at(i + 1) - 33 < min_qual) {
                            bad_basecall = true;
                        }
                    } catch (std::out_of_range &exception) { }
                }
            }
            if (bad_basecall and effective_depth[i]) {
                effective_depth[i] = 0;
            }
        }

        // in second pass, update depths for only nucs within mutations, 
        // and filter mutations by q-scores (and mutation_type, if given)
        //
        // mutation_type possible values: mismatch gap insert gap_multi insert_multi complex
        for (int i = 0; i < mutations.size(); ++i) {
            Mutation m(mutations[i]);
            bool bad_mutation = false;

            // limit to specific mutation type if given
            // FIXME: this relies on the template class being Mutation
            std::vector<std::string> mismatch_tags = {"AT", "AG", "AC",
                                                      "TA", "TG", "TC",
                                                      "GA", "GT", "GC",
                                                      "CA", "CT", "CG",
                                                      "multinuc_mismatch"};
            std::vector<std::string> insert_tags = {"-A", "-T", "-G", "-C", "-N"};
            std::vector<std::string> gap_tags = {"A-", "T-", "G-", "C-"};

            if (mutation_type.length() > 0) {
                if (mutation_type == "mismatch") {
                    if (std::find(mismatch_tags.begin(), mismatch_tags.end(), m.tag) == mismatch_tags.end()) {
                        bad_mutation = true;
                    }
                } else if (mutation_type == "insert") {
                    if (std::find(insert_tags.begin(), insert_tags.end(), m.tag) == insert_tags.end()) {
                        bad_mutation = true;
                    }
                } else if (mutation_type == "insert_multi") {
                    if (m.tag != "multinuc_insertion") {
                        bad_mutation = true;
                    }
                } else if (mutation_type == "gap") {
                    if (std::find(gap_tags.begin(), gap_tags.end(), m.tag) == gap_tags.end()) {
                        bad_mutation = true;
                    }
                } else if (mutation_type == "gap_multi") {
                    if (m.tag != "multinuc_deletion") {
                        bad_mutation = true;
                    }
                } else if (mutation_type == "complex") {
                    if (m.tag != "complex_deletion" && m.tag != "complex_insertion") {
                        bad_mutation = true;
                    }
                }
            }

            // look at basecalls within mutation
            if (not bad_mutation) {
                for (auto c : m.qual) {
                    if (c - 33 < min_qual) {
                        bad_mutation = true;
                        break;
                    }
                }
            }
            // look at neighbor on left
            if (not bad_mutation) {
                try {
                    int k = right_mut_indices.at(m.left + 1 - left);
                    if (k != -1 and k != i) {
                        Mutation neighbor(mutations[k]);
                        if (neighbor.seq.length() > 0) {
                            // look at right-most basecall in mutation
                            if (neighbor.qual.back() - 33 < min_qual) {
                                bad_mutation = true;
                            }
                        } else {
                            // look at basecall on other side of gap
                            int n = neighbor.left - left;
                            if (qual.at(n) - 33 < min_qual) {
                                bad_mutation = true;
                            }
                        }
                    } else {
                        // no mutation on left, just look at basecall
                        int n = m.left - left;
                        if (qual.at(n) - 33 < min_qual) {
                            bad_mutation = true;
                        }
                    }
                } catch (std::out_of_range &exception) { }
            }
            // look at neighbor on right
            if (not bad_mutation) {
                try {
                    int k = left_mut_indices.at(m.right - 1 - left);
                    if (k != -1 and k != i) {
                        Mutation neighbor(mutations[k]);
                        if (neighbor.seq.length() > 0) {
                            // look at left-most basecall in mutation
                            if (neighbor.qual.front() - 33 < min_qual) {
                                bad_mutation = true;
                            }
                        } else {
                            // look at basecall on other side of gap
                            int n = neighbor.right - left;
                            if (qual.at(n) - 33 < min_qual) {
                                bad_mutation = true;
                            }
                        }
                    } else {
                        // no mutation on right, just look at basecall
                        int n = m.right - left;
                        if (qual.at(n) - 33 < min_qual) {
                            bad_mutation = true;
                        }
                    }
                } catch (std::out_of_range &exception) { }
            }

            if (bad_mutation) {
                // exclude entire covered region from depths
                for (int n = m.left + 1 - left;
                     n < m.right - left;
                     ++n) {
                    try {
                        effective_depth.at(n) = 0;
                    } catch (std::out_of_range &exception) { }
                }
                excluded_mutations.push_back(Mutation(m));
            } else {
                if (variant_mode) {
                    // include entire covered region (even in gap positions)
                    for (int n = m.left + 1 - left;
                         n <= m.right - 1 - left;
                         ++n) {
                        try {
                            effective_depth.at(n) = 1;
                        } catch (std::out_of_range &exception) { }
                    }
                } else {
                    // exclude covered region, include inferred adduct site
                    for (int n = m.left + 1 - left;
                         n < m.right - 1 - left;
                         ++n) {
                        try {
                            effective_depth.at(n) = 0;
                        } catch (std::out_of_range &exception) { }
                    }
                    try {
                        effective_depth.at(m.right - 1 - left) = 1;
                    } catch (std::out_of_range &exception) { }
                }
                included_mutations.push_back(Mutation(m));
            }
        }

        // also mark mutation locations in simplified array
        std::vector<bool> local_effective_count(len, 0);
        for (auto &m:included_mutations) {
            try {
                local_effective_count.at(m.right - 1 - left) = 1;
            } catch (std::out_of_range &exception) { }
        }

        return boost::make_tuple(effective_depth,
                                 local_effective_count,
                                 included_mutations,
                                 excluded_mutations);
    };

    // FIXME: replace with Read constructor?
    /**
     * @brief Used in mutation_counter to read input from a stream.
     */
    boost::tuple<std::string, int, int, int, int, std::vector<bool>, std::vector<bool>, std::vector<bool>, std::vector<Mutation>>
    parseProcessedMutations(const std::string &line) {
        using std::stoi;

        std::string read_type;
        std::string read_id;
        int mapping_category;
        int primer_pair; // negative indicates no associated primer pair
        int left;
        int right;
        std::vector<bool> mapping_depth;
        std::vector<bool> local_effective_depth;
        std::vector<bool> local_effective_count;
        std::vector<Mutation> mutations;

        std::vector<std::string> fields;
        std::string trimmed = boost::trim_right_copy_if(line, boost::is_any_of("\n\r"));
        boost::split(fields, trimmed, boost::is_any_of("\t"), boost::token_compress_off);

        if (fields.size() < 10) {
            throw std::runtime_error("Error: unable to parse incomplete line. "
                                             "Trimmed line: '" + trimmed + "'");
        }
        try {
            read_type = fields[0];
            read_id = fields[1];
            left = stoi(fields[2]);
            right = stoi(fields[3]);
            mapping_category = util::indexOf(mapping_categories, fields[4]);
            primer_pair = stoi(fields[5]);
        } catch (std::exception &err) {
            throw std::runtime_error("Error: line is incorrectly formatted (couldn't parse left or right position)."
                                             "Trimmed line: '" + trimmed + "'");
        }
        for (char c : fields[6]) {
            mapping_depth.push_back((c == '1'));
        }
        for (char c : fields[7]) {
            local_effective_depth.push_back((c == '1'));
        }
        for (char c : fields[8]) {
            local_effective_count.push_back((c == '1'));
        }

        mutations = stringToMutationVect(fields[9]);

        return boost::make_tuple(read_id,
                                 mapping_category,
                                 primer_pair,
                                 left,
                                 right,
                                 mapping_depth,
                                 local_effective_depth,
                                 local_effective_count,
                                 mutations);
    };


    Read
    mergeMatePairs(std::vector<Read> reads) {
        std::vector<Mutation> mutations;
        std::string seq;
        std::string qual;
        std::vector<bool> depth;
        std::vector<bool> mapped_depth;
        int left;
        int right;
        int length;

        Read r1(reads[0]);
        Read r2(reads[1]);

        // seq can be merged no problem, there should be no conflicts
        // for qual, just use the higher score in overlap
        // conflicting mutations are trickier
        
        left = std::min(r1.left, r2.left);
        right = std::max(r1.left + r1.seq.size() - 1,
                         r2.left + r2.seq.size() - 1);
        length = right - left + 1;
        seq.resize(length);
        qual.resize(length);

        Read simple_merged = mergeMatePairsSimple(reads);
        mapped_depth = simple_merged.mapped_depth;

        char one, two;
        for (int i = 0; i < length; i++) {
            seq.at(i) = '_';
            qual.at(i) = '~';
            try {
                one = r1.seq.at(i - r1.left + left);
                seq.at(i) = one;
            } catch (const std::out_of_range &exception) { }
            try {
                two = r2.seq.at(i - r2.left + left);
                seq.at(i) = two;
            } catch (const std::out_of_range &exception) { }

            one = '~';
            two = '~';
            try {
                one = r1.qual.at(i - r1.left + left);
            } catch (const std::out_of_range &exception) { }
            try {
                two = r2.qual.at(i - r2.left + left);
            } catch (const std::out_of_range &exception) { }
            if (one != '~' and two != '~') {
                if ((int) one >= (int) two) {
                    qual.at(i) = one;
                } else {
                    qual.at(i) = two;
                }
            }
            else if (two == '~' and one != '~') {
                qual.at(i) = one;
            }
            else if (one == '~' and two != '~') {
                qual.at(i) = two;
            }
        }


        // now need to merge (possibly conflicting) mutations
        // index R1 mutations by left-most nuc
        std::vector<std::vector<Mutation>> indexed_r1_muts, indexed_r2_muts;
        indexed_r1_muts.resize(length);
        indexed_r2_muts.resize(length);
        for (auto &m : r1.mutations) {
            Mutation mc(m);
            indexed_r1_muts.at(m.left - left).push_back(mc);
        }
        // index R2 mutations by left-most nuc
        for (auto &m : r2.mutations) {
            Mutation mc(m);
            indexed_r2_muts.at(m.left - left).push_back(mc);
        }
        // group all overlapping mutations
        std::vector<MutationGroup> mutation_groups;
        MutationGroup mutation_group;
        //Mutation m1, m2;
        std::vector<Mutation> v1, v2;
        for (int i = 0; i < length; i++) {
            v1 = indexed_r1_muts[i];
            v2 = indexed_r2_muts[i];
            for (auto &m1 : v1) {
                if (mutation_group.r1_mutations.size() == 0 and
                    mutation_group.r2_mutations.size() == 0) {
                    mutation_group.r1_mutations.push_back(m1);
                    mutation_group.left = m1.left;
                    mutation_group.right = m1.right;
                } else {
                    // check for overlap
                    if (m1.left < mutation_group.right) {
                        mutation_group.r1_mutations.push_back(m1);
                        mutation_group.right = std::max(m1.right, mutation_group.right);
                    } else {
                        mutation_groups.push_back(mutation_group);
                        mutation_group = MutationGroup();
                        mutation_group.r1_mutations.push_back(m1);
                        mutation_group.left = m1.left;
                        mutation_group.right = m1.right;
                    }
                }
            }
            for (auto &m2 : v2) {
                if (mutation_group.r1_mutations.size() == 0 and
                    mutation_group.r2_mutations.size() == 0) {
                    mutation_group.r2_mutations.push_back(m2);
                    mutation_group.left = m2.left;
                    mutation_group.right = m2.right;
                } else {
                    // check for overlap
                    if (m2.left < mutation_group.right) {
                        mutation_group.r2_mutations.push_back(m2);
                        mutation_group.right = std::max(m2.right, mutation_group.right);
                    } else {
                        mutation_groups.push_back(mutation_group);
                        mutation_group = MutationGroup();
                        mutation_group.r2_mutations.push_back(m2);
                        mutation_group.left = m2.left;
                        mutation_group.right = m2.right;
                    }
                }
            }
        }
        // handle final group
        if (mutation_group.r1_mutations.size() > 0 or
            mutation_group.r2_mutations.size() > 0) {
            mutation_groups.push_back(mutation_group);
        }

        r1.depth.resize(r1.seq.size(), 1);
        r2.depth.resize(r2.seq.size(), 1);

        for (auto &mg : mutation_groups) {

            // resolve overlaps by choosing the higher quality read mutations
            float mean_qual_r1, mean_qual_r2;
            int num = 0;
            int denom = 0;
            if (mg.r1_mutations.size() > 0) {
                for (std::vector<Mutation>::const_iterator m = mg.r1_mutations.begin();
                     m != mg.r1_mutations.end();
                     ++m) {
                    for (int i = 0; i < (*m).qual.size(); i++) {
                        num += (int) (*m).qual.at(i);
                        denom++;
                    }
                    // also consider basecalls neighboring the mutation
                    try {
                        num += (int) r1.qual.at((*m).left - r1.left);
                        denom++;
                    } catch (const std::out_of_range &exception) { }
                    try {
                        num += (int) r1.qual.at((*m).right - r1.left);
                        denom++;
                    } catch (const std::out_of_range &exception) { }
                }
            } else {
                // mutations in r2, but not in r1.
                // consider basecalls in r1 over same region as mutation group in r2, if r1 covers
                // left and right
                int lindex = mg.left - r1.left;
                int rindex = mg.right - r1.left;
                if (lindex >= 0 and lindex < r1.qual.size() and
                    rindex >= 0 and rindex < r1.qual.size()) {
                    for (int p = mg.left; p <= mg.right; p++) {
                        num += (int) r1.qual.at(p - r1.left);
                        denom++;
                    }
                }
            }
            if (denom > 0.0) {
                mean_qual_r1 = float(num) / float(denom);
            } else {
                mean_qual_r1 = 0.0;
            }
            num = 0;
            denom = 0;
            if (mg.r2_mutations.size() > 0) {
                for (std::vector<Mutation>::const_iterator m = mg.r2_mutations.begin();
                     m != mg.r2_mutations.end();
                     ++m) {
                    for (int i = 0; i < (*m).qual.size(); i++) {
                        num += (int) (*m).qual.at(i);
                        denom++;
                    }
                    // also consider basecalls neighboring the mutation
                    try {
                        num += (int) r2.qual.at((*m).left - r2.left);
                        denom++;
                    } catch (const std::out_of_range &exception) { }
                    try {
                        num += (int) r2.qual.at((*m).right - r2.left);
                        denom++;
                    } catch (const std::out_of_range &exception) { }
                }
            } else {
                // mutations in r1, but not in r2.
                // consider basecalls in r2 over same region as mutation group in r1, if r2 covers
                // left and right
                int lindex = mg.left - r2.left;
                int rindex = mg.right - r2.left;
                if (lindex >= 0 and lindex < r2.qual.size() and
                    rindex >= 0 and rindex < r2.qual.size()) {
                    for (int p = mg.left; p <= mg.right; p++) {
                        num += (int) r2.qual.at(p - r2.left);
                        denom++;
                    }
                }
            }
            if (denom > 0.0) {
                mean_qual_r2 = float(num) / float(denom);
            } else {
                mean_qual_r2 = 0.0;
            }

            int selected_read = READ1;
            if (mean_qual_r2 > mean_qual_r1) {
                selected_read = READ2;
            }
            if (selected_read == READ1) {
                for (int i = 0; i < mg.r1_mutations.size(); i++) {
                    mutations.push_back(mg.r1_mutations[i]);
                }
                // remove excluded mutations from r2 effective depth
                if (mg.r2_mutations.size() > 0) {
                    for (int n = minLeft(mg.r2_mutations) + 1;
                         n < maxRight(mg.r2_mutations);
                         n++) {
                        try {
                            r2.depth.at(n - r2.left) = 0;
                        } catch (std::out_of_range &ex) { }
                    }
                }

            } else {
                for (int i = 0; i < mg.r2_mutations.size(); i++) {
                    mutations.push_back(mg.r2_mutations[i]);
                }
                // remove excluded mutations from r1 effective depth
                if (mg.r1_mutations.size() > 0) {
                    for (int n = minLeft(mg.r1_mutations) + 1;
                         n < maxRight(mg.r1_mutations);
                         n++) {
                        try {
                            r1.depth.at(n - r1.left) = 0;
                        } catch (std::out_of_range &ex) { }
                    }
                }

            }


        }

        // combine effective depths from r1 and r2
        depth.resize(seq.size(), 0);
        for (int i = 0; i < depth.size(); i++) {
            bool d1, d2;
            d1 = 0;
            d2 = 0;
            int loc = left + i;
            int r1_index = loc - r1.left;
            int r2_index = loc - r2.left;
            try {
                d1 = r1.depth.at(r1_index);
            } catch (std::out_of_range &ex) { }
            try {
                d2 = r2.depth.at(r2_index);
            } catch (std::out_of_range &ex) { }
            bool m = 0;
            if (d1 or d2) { m = 1; }
            depth[i] = m;
        }


        return Read(left,
                    right,
                    seq)
                .setReadType(PAIRED)
                .setStrand(FORWARD)
                .setId(r1.id)
                .setMappingCategory(r1.mapping_category)
                .setMappedDepth(mapped_depth)
                .setPrimerPair(r1.primer_pair)
                .setQual(qual)
                .setMutations(mutations)
                .setDepth(depth);
    }


    Read
    trimRightEnd(const Read &read,
                 const int exclude_3prime,
                 const bool debug = false) {
        /*
         * for stripping random RT primer sites
         */

        Read trimmed{};

        int read_type;
        std::string read_id;
        int left, right;
        std::vector<bool> depth;
        std::vector<Mutation> stripped_mutations;
        std::string seq;
        std::string qual;

        if (read.depth.size() > 0) {
            for (auto b : read.depth) {
                depth.push_back(b);
            }
        } else if (read.mapped_depth.size() > 0) {
            for (auto b : read.mapped_depth) {
                depth.push_back(b);
            }
        } else {
            depth.resize(read.seq.size(), 1);
        }

        if (debug_out) {
            debug_out << read << std::flush;
        }

        if (read.read_type == MERGED or
            read.read_type == PAIRED) {
            boost::tie(stripped_mutations,
                       depth) = stripEnd(read,
                                         depth,
                                         exclude_3prime,
                                         RIGHT,
                                         debug);
            if (debug_out) {
                debug_out << "trimmed " << exclude_3prime << " nts from "
                        "right end of read (for handling random primers)\n"
                << std::flush;
            }
        }
        else if (read.read_type == PAIRED_R1 or
                 read.read_type == UNPAIRED_R1 or
                 read.read_type == UNPAIRED) {
            if (read.strand == REVERSE) {
                int end_to_strip = RIGHT;
                // NOTE: this will not handle cases where the R1 read extends into/past the RT priming site,
                //       since we can't be sure where the random priming site is actually located
                boost::tie(stripped_mutations,
                           depth) = stripEnd(read,
                                             depth,
                                             exclude_3prime,
                                             end_to_strip,
                                             debug);
                if (debug_out) {
                    debug_out << "trimmed " << exclude_3prime << " nts from "
                            "right end of read (for handling random primers)\n"
                    << std::flush;
                }
            } else {
                stripped_mutations = read.mutations;
                if (debug_out) {
                    debug_out << "didn't trim right end of forward read\n"
                    << std::flush;
                }
            }
        } else if (read.read_type == PAIRED_R2 or
                   read.read_type == UNPAIRED_R2) {
            if (read.strand == REVERSE) {
                int end_to_strip = RIGHT;
                boost::tie(stripped_mutations,
                           depth) = stripEnd(read,
                                             depth,
                                             exclude_3prime,
                                             end_to_strip,
                                             debug);
                if (debug_out) {
                    debug_out << "trimmed " << exclude_3prime << " nts from "
                            "right end of read (for handling random primers)\n"
                    << std::flush;
                }
            } else {
                stripped_mutations = read.mutations;
                if (debug_out) {
                    debug_out << "didn't trim right end of forward read\n"
                    << std::flush;
                }
            }
        }
        seq = read.seq;
        qual = read.qual;
        left = read.left;
        right = left + seq.size() - 1;
        read_type = read.read_type;
        read_id = read.id;

        trimmed.setLeft(left)
                .setRight(right)
                .setReadType(read_type)
                .setId(read_id)
                .setSeq(seq)
                .setQual(qual)
                .setMutations(stripped_mutations)
                .setDepth(depth)
                .setPrimerPair(read.primer_pair);

        if (read.mapped_depth.size() > 0) {
            trimmed.setMappedDepth(read.mapped_depth);
        }

        if (debug_out) {
            debug_out << trimmed << std::flush;
        }

        return trimmed;
    }


    /**
     *
     * @direction process multinuc mutations from right or left-most nuc.
     *            >0 for right-to-left, <=0 for left-to-right. WARNING: left-to-right
     *            direction is not functional yet.
     * @brief
     *          handles:
     *          end trimming to account for random primers or directed primers,
                    - update to handle R1/R2/paired/merged differently and
                      intelligently handle multiple primer pairs
     *          shifting ambiguously aligned muts,
     *          collapsing muts to infer adduct location,
     *          classifying muts,
     *          post-hoc merging mate pairs *new*,
     *          applying quality filters and counting local depth,

     *            Return: tuple containing vector of "good" mutations and vector of bools
     *                    indicating "good" basecalls over the local_target_seq region, excluding
     *                    some range of nucs on right end according to exclude_3prime. Now using Read
     *                    class
     */
    Read
    processMutations(const std::vector<Read> reads,
                     const int direction,
                     const bool right_align_ambig_dels,
                     const bool right_align_ambig_ins,
                     const int max_internal_match,
                     const int min_qual,
                     const int exclude_3prime,
                     const std::string &mutation_type,
                     const bool variant_mode,
                     const bool trim_amplicon_primers,
                     const PrimerPair primer_pair,
                     const bool debug) {
        // FIXME: clarify/document coordinate system used here (Read.left vs. ScanningDeque.target_pos, etc., 0/1-based)
        /* example coordinates:
         *
         *  RNA ("genome" coords)
         *  0                                                        100
         *  |_________________________________________________________|
         *  Deque bounds:   ------------------------------
         *                  |    Read bounds: --------------------
         *   Deque.target_pos                 |
         *                            Read.left
         */

        /*if (debug_out) {
            debug_out << "running processMutations()\n";
            debug_out << std::flush;
        }*/

        Read read{};

        if (reads.size() == 2) {
            if (debug_out) {
                debug_out << reads[0];
                debug_out << reads[1];
                debug_out << std::flush;
            }

            // merge mate pair mutations before any trimming
            std::vector<Mutation> mutations;
            read = mergeMatePairs(reads);

            if (debug_out) {
                debug_out << "merged mate pairs\n"
                          << read
                          << std::flush;
            }
            if (debug) {
                std::cout << "merged mate pairs\n" << std::flush;
            }
        } else {
            read = reads[0];
            read.setDepth(read.mapped_depth);
        }

        if (trim_amplicon_primers) {
            // amplicon primer trim
            read.stripPrimers(primer_pair);
            if (debug_out) {
                debug_out << "trimmed amplicon primer sites\n"
                          << read
                          << std::flush;
            }
        } else {
            // random RT primer trim
            read.trimRightEnd(exclude_3prime);
        }

        

        std::vector<Mutation> classified_mutations;
        if (not variant_mode) {
            // don't collapse nearby mutations or adjust ambiguous placements if only
            // doing sequence variant detection, since these operations make it more difficult
            // to estimate the frequency of SNPs
            read.shiftAmbigIndels(right_align_ambig_dels,
                                  right_align_ambig_ins);
            if (debug_out) {
                debug_out << "shifted ambiguously aligned mutations\n"
                        << read
                        << std::flush;
            }

            // TODO: options for transcriptomes with reverse strand annotations (collapse to left edge instead of right)?
            read.collapseMutations(max_internal_match);
            if (debug_out) {
                debug_out << "collapsed nearby mutations\n"
                 << read
                 << std::flush;
            }
        }

        read.classifyMutations();
        if (debug_out) {
            debug_out << "classified mutations\n"
             << read
             << std::flush;
        }
        if (debug) {
            std::cout << "classified mutations\n" << std::flush;
        }


        read.filterQscoresCountDepths(min_qual,
                                      mutation_type,
                                      variant_mode);
        if (debug) { std::cout << read << std::flush; }
        if (debug_out) {
            debug_out << "filtered Q-scores and inferred adduct locations\n"
                      << read
                      << std::flush;
        }

        return read;
    }
}


// FIXME: should probably move these wrappers to be closer to the main Read class definitions

Read&
Read::trimRightEnd(const int exclude_3prime){
    *this = mutation::trimRightEnd(*this,
                                   exclude_3prime,
                                   false);
    return *this;
}

Read&
Read::stripPrimers(const PrimerPair& primer_pair) {
    std::vector<Mutation> stripped_mutations;
    std::vector<bool> depth;
    boost::tie(stripped_mutations,
               depth) = mutation::stripPrimers(this->mutations,
                                               this->left,
                                               this->depth,
                                               primer_pair,
                                               false);
    (*this).setMutations(stripped_mutations)
           .setDepth(depth);
    return *this;
}

Read&
Read::shiftAmbigIndels(const bool right_align_ambig_dels,
                       const bool right_align_ambig_ins) {
    std::vector<Mutation> shifted_mutations =
            mutation::shiftAmbigIndels(this->mutations,
                                       this->seq,
                                       this->qual,
                                       this->left,
                                       right_align_ambig_dels,
                                       right_align_ambig_ins);
    (*this).setMutations(shifted_mutations);
    return *this;
}

Read&
Read::collapseMutations(const int max_internal_match) {
    std::vector<Mutation> collapsed_mutations =
            mutation::collapseMutations(this->mutations,
                                        max_internal_match,
                                        this->seq,
                                        this->qual,
                                        this->left);
    (*this).setMutations(collapsed_mutations);
    return *this;
}

Read&
Read::classifyMutations() {
    std::vector<Mutation> classified_mutations =
            mutation::classifyMutations(this->mutations,
                                        this->seq,
                                        this->qual,
                                        this->left);
    (*this).setMutations(classified_mutations);
    return *this;
}

Read&
Read::filterQscoresCountDepths(const int min_qual,
                               const std::string &mutation_type,
                               const bool variant_mode) {
    std::vector<bool> depth;
    std::vector<bool> count;
    std::vector<Mutation> included_mutations;
    std::vector<Mutation> excluded_mutations;
    boost::tie(depth,
               count,
               included_mutations,
               excluded_mutations) =
            mutation::filterQscoresCountDepths(this->mutations,
                                               this->seq,
                                               this->qual,
                                               this->depth,
                                               this->left,
                                               min_qual,
                                               mutation_type,
                                               variant_mode);
    (*this).setMutations(included_mutations)
           .setDepth(depth)
           .setCount(count);
    return *this;
}

#endif //SHAPEMAPPER_MUTATIONPROCESSING_H
