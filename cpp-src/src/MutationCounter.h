/** @file
 * @brief Classes for counting sequencing depth, sequence variants, and
 *        reverse transcription mutations.
 */

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2017 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#ifndef SHAPEMAPPER_MUTATIONCOUNTER_H
#define SHAPEMAPPER_MUTATIONCOUNTER_H

#include <map>
#include <deque>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>

#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/newline.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "Mutation.h"
#include "Histogram.h"

namespace BF = boost::filesystem;
namespace BI = boost::iostreams;

using namespace mutation;

/**
 * @brief Output columns for mutation counts
 */
std::vector<std::string> column_names = mutation_classes;


/**
 * @brief Deque (list with fast append/delete from both ends) representing
 *        any form of data covering a range over an alignment target.
 */
template<typename T>
class ScanningDeque {
public:
    int target_pos = 0;
    /** Leftmost target position (0-based) covered by this region */
    std::deque<T> deq;

    /**
     * @brief Print values over a given range (in local deque coords,
     *        0-based inclusive; NOT alignment target coords).
     */
    // FIXME: remove zero-assignment or document reason for including
    virtual std::string printValues(const int, const int) = 0;

    /**
     * @brief Print data up to new_target_left (exclusive), then
     *        delete those data and update target_pos to new_target_left
     */
    std::string updateLeftBound(const int new_target_left) {
        std::string s;
        if (new_target_left > target_pos) {
            int n_to_drop = new_target_left - target_pos;
            s = printValues(0, n_to_drop - 1);
            // note: erase() right bound is exclusive
            deq.erase(deq.begin(), deq.begin() + n_to_drop);
            target_pos += n_to_drop;
        }
        return s;
    };

    /**
     * @brief Resize deque on right to accommodate given alignment target location.
     */
    void updateRightBound(const int new_target_right) {
        int current_right = target_pos + deq.size() - 1;
        if (new_target_right > current_right) {
            //std::cout << "Resized deque from " << deq.size();
            deq.resize(new_target_right - target_pos + 1);
            //std::cout << " to " << deq.size() << "." << std::endl;
        }
    }

};


struct VariantRow {
    int depth;
    std::map<Mutation, int> counts;
};

/**
 * @brief Sequence variants and variant counts covering a range of
 *        target alignment positions. Variants are indexed by leftmost
 *        unchanged alignment target nucleotide.
 */
class VariantCounter : public ScanningDeque<VariantRow> {
public:
    std::string printValues(const int left_inclusive,
                            const int right_inclusive) {
        std::string o;
        for (int i = left_inclusive; i <= right_inclusive; ++i) {
            o += std::to_string(deq[i].depth) + ' ';
            for (auto const &k : deq[i].counts) {
                Mutation m = k.first;
                int count = k.second;
                o += '(' + std::to_string(m.left) + '-';
                o += std::to_string(m.right) + ", ";
                o += "\"" + m.seq + "\", ";
                o += std::to_string(count) + ")";
            }
            o += '\n';
        }
        return o;
    }

    std::string printAllValues() {
        return printValues(0, deq.size() - 1);
    }

    void updateCounts(const std::vector<Mutation> &mutations,
                      const std::string &local_target_seq,
                      const std::string &local_target_qual,
                      const int left_target_pos,
                      const int min_qual,
                      const int exclude_3prime,
                      const std::string &mutation_type) {
        // mutations are grouped by leftmost unchanged alignment target nuc
        std::vector<Mutation> stripped_mutations = strip_3prime(mutations,
                                                                local_target_seq,
                                                                local_target_qual,
                                                                left_target_pos,
                                                                exclude_3prime);

        // convert to TaggedMutations because filterQscoresCountDepths function expects tag field, even if blank
        // FIXME: reorganize things so this is not an issue
        std::vector<TaggedMutation> tagged_mutations;
        for (auto &m : stripped_mutations) {
            tagged_mutations.push_back(TaggedMutation(m, ""));
        }

        std::vector<bool> local_effective_depth;
        std::vector<bool> local_effective_count;
        std::vector<TaggedMutation> included_mutations;
        std::vector<TaggedMutation> excluded_mutations;
        boost::tie(local_effective_depth,
                   local_effective_count,
                   included_mutations,
                   excluded_mutations) = filterQscoresCountDepths(tagged_mutations,
                                                                  local_target_seq,
                                                                  local_target_qual,
                                                                  left_target_pos,
                                                                  exclude_3prime,
                                                                  min_qual,
                                                                  mutation_type,
                                                                  true); // variant_mode on

        // update depths in deque coordinates
        for (int i = 0; i < local_effective_depth.size(); ++i) {
            // convert index to global coordinates, then to deque coords
            int n = i + left_target_pos - target_pos;
            if (local_effective_depth[i]) {
                //if (true) {
                try {
                    ++deq.at(n).depth;
                } catch (std::out_of_range &exception) { }
            }
        }

        // update variant counts
        for (auto &m : included_mutations) {
            m.qual = ""; // blank out quality scores. 
            // otherwise the same variant seqs with different q-scores
            // would be added to counts separately
            try {
                ++deq.at(m.left - target_pos).counts[m];
            }
            catch (std::out_of_range &exception) { ;
            }
        }
    }
};


/**
 * @brief Mutation counts and calculated read depths covering a range of
 *        target alignment positions. Mutations are indexed by rightmost
 *        changed alignment target nucleotide (will add options for
 *        reverse strand in the future, for mapping to transcriptomes).
 */
class MutationCounter : public ScanningDeque<std::map<std::string, int>> {
public:

    Histogram read_lengths;
    Histogram mutations_per_read;

    MutationCounter() : read_lengths("Read lengths", 0, 1000, 21),
                        mutations_per_read("Mutations per read", 0, 20, 21) { }

    std::string printValues(const int left_inclusive,
                            const int right_inclusive) {
        std::string o;
        for (int i = left_inclusive; i <= right_inclusive; ++i) {
            for (std::vector<std::string>::const_iterator it = column_names.begin();
                 it != column_names.end();
                 ++it) {
                o += std::to_string(deq[i][*it]);
                if (it + 1 != column_names.end()) {
                    o += '\t';
                }
            }
            o += '\n';
        }
        return o;
    }

    std::string printHeader() {
        std::string o;
        for (std::vector<std::string>::const_iterator it = column_names.begin();
             it != column_names.end();
             ++it) {
            o += *it;
            if (it + 1 != column_names.end()) {
                o += '\t';
            }
        }
        o += '\n';
        return o;
    }

    std::string printAllValues() {
        return printValues(0, deq.size() - 1);
    }

    std::string printHistograms() {
        return read_lengths.printFreqTable("range")+"\n"+mutations_per_read.printFreqTable();
    }

    /**
     *
     * @direction Count multinuc mutations from right or left-most nuc.
     *            >0 for right-to-left, <=0 for left-to-right. WARNING: left-to-right
     *            direction is not functional yet.
     *
     *            Return: tuple containing vector of "good" mutations and vector of bools
     *                    indicating "good" basecalls over the local_target_seq region, excluding
     *                    some range of nucs on right end according to exclude_3prime
     */
    boost::tuple<std::vector<TaggedMutation>,
                 std::vector<bool>,
                 std::vector<bool>>
    updateCounts(const std::string &read_id,
                 const std::vector<Mutation> &mutations,
                 const std::string &local_target_seq,
                 const std::string &local_target_qual,
                 const int left_target_pos,
                 const int direction,
                 const bool separate_ambig_counts,
                 const bool right_align_ambig_dels,
                 const bool right_align_ambig_ins,
                 const int max_internal_match,
                 const int min_qual,
                 const int exclude_3prime,
                 const std::string &mutation_type,
                 const bool debug) {
        // FIXME: clarify/document coordinate system used here (left_target_pos vs. ScanningDeque.target_pos, etc., 0/1-based)
        /* example coordinates:
         *
         *  RNA ("genome" coords)
         *  0                                                        100
         *  |_________________________________________________________|
         *  Deque bounds:   ------------------------------
         *                  |    Read bounds: --------------------
         *             target_pos             |
         *                             left_target_pos
         */

        if (debug){ std::cout << std::endl << "read_id: " << read_id << std::endl; }
        if (debug){ std::cout << "local_target_seq:  " << local_target_seq << std::endl; }
        if (debug){ std::cout << "local_target_qual: " << local_target_qual << std::endl; }
        if (debug){ std::cout << "mutations: " << toString(mutations) << std::endl; }

        std::vector<Mutation> stripped_mutations = strip_3prime(mutations,
                                                                local_target_seq,
                                                                local_target_qual,
                                                                left_target_pos,
                                                                exclude_3prime);
        if (debug) { std::cout << "stripped_mutations: " << toString(stripped_mutations) << std::endl; }

        std::vector<TaggedMutation> shifted_mutations = shiftAmbigIndels(stripped_mutations,
                                                                         local_target_seq,
                                                                         local_target_qual,
                                                                         left_target_pos,
                                                                         right_align_ambig_dels,
                                                                         right_align_ambig_ins);
        if (debug) { std::cout << "shifted_mutations: " << toString(shifted_mutations) << std::endl; }

        // TODO: options for transcriptomes with reverse strand annotations (collapse to left edge instead of right)
        std::vector<TaggedMutation> collapsed_mutations = collapseMutations(shifted_mutations,
                                                                            max_internal_match,
                                                                            local_target_seq,
                                                                            local_target_qual,
                                                                            left_target_pos);
        if (debug) { std::cout << "collapsed_mutations: " << toString(collapsed_mutations) << std::endl; }

        std::vector<TaggedMutation> classified_mutations = classifyMutations(collapsed_mutations,
                                                                             local_target_seq,
                                                                             local_target_qual,
                                                                             left_target_pos);
        if (debug) { std::cout << "classified_mutations: " << toString(classified_mutations) << std::endl; }

        std::vector<bool> local_effective_depth;
        std::vector<bool> local_effective_count;
        std::vector<TaggedMutation> included_mutations;
        std::vector<TaggedMutation> excluded_mutations;
        boost::tie(local_effective_depth,
                   local_effective_count,
                   included_mutations,
                   excluded_mutations) = filterQscoresCountDepths(classified_mutations,
                                                                  local_target_seq,
                                                                  local_target_qual,
                                                                  left_target_pos,
                                                                  exclude_3prime,
                                                                  min_qual,
                                                                  mutation_type,
                                                                  false); // variant_mode off

        // update histograms
        int len = local_target_seq.length() - exclude_3prime;
        len = std::max(len,0);
        read_lengths.count(len);
        mutations_per_read.count(included_mutations.size());

        if (debug) {
            std::string s = "local_effective_depth: ";
            for (auto b : local_effective_depth) {
                s += to_string(b);
            }
            std::cout << s << std::endl;
        }

        if (debug){
            std::string s = "local_effective_count: ";
            for (auto b : local_effective_count){
                s += to_string(b);
            }
            std::cout << s << std::endl;
        }
        if (debug){ std::cout << "excluded_mutations: " << toString(excluded_mutations) << std::endl; }
        if (debug){ std::cout << "included_mutations: " << toString(included_mutations) << std::endl; }

        //std::vector<int> occluded_depth(local_target_seq.length());

        // update mutation counts
        for (auto mut : included_mutations) {
            std::string s = mut.tag;

            if (separate_ambig_counts) {
                if (mut.ambig) {
                    s += "_ambig";
                }
            }

            if (direction > 0) {
                try {
                    // update mutation counts
                    ++deq.at(mut.right - target_pos - 1)[s];
                }
                catch (std::out_of_range &exception) { }
                //std::cout << "updated deque count" << std::endl;
            } else { // direction <= 0
                try {
                    // update mutation counts
                    ++deq.at(mut.left - target_pos + 1)[s];
                }
                catch (std::out_of_range &exception) { }
            }
        }


        // update depths in deque coordinates
        for (int n = left_target_pos - target_pos;
             n < left_target_pos + local_target_seq.length() - target_pos;
             ++n) {
            try {
                ++deq.at(n)["read_depth"];
            } catch (std::out_of_range &exception) { }
        }
        for (int i = 0; i < local_effective_depth.size(); ++i) {
            // convert index to global coordinates, then to deque coords
            int n = i + left_target_pos - target_pos;
            if (local_effective_depth[i]) {
                try {
                    ++deq.at(n)["effective_depth"];
                } catch (std::out_of_range &exception) { }
            }
        }

        return boost::make_tuple(included_mutations,
                                 local_effective_depth,
                                 local_effective_count);
    }

    /**
     * Default (no-debug)
     */
    boost::tuple<std::vector<TaggedMutation>,
                 std::vector<bool>,
                 std::vector<bool>>
    updateCounts(const std::string &read_id,
                 const std::vector<Mutation> &mutations,
                 const std::string &local_target_seq,
                 const std::string &local_target_qual,
                 const int left_target_pos,
                 const int direction,
                 const bool separate_ambig_counts,
                 const bool right_align_ambig_dels,
                 const bool right_align_ambig_ins,
                 const int max_internal_match,
                 const int min_qual,
                 const std::string &mutation_type,
                 const int exclude_3prime) {
        std::vector<TaggedMutation> included_mutations;
        std::vector<bool> local_effective_depth;
        std::vector<bool> local_effective_count;
        boost::tie(included_mutations,
                   local_effective_depth,
                   local_effective_count) = updateCounts(read_id,
                                                         mutations,
                                                         local_target_seq,
                                                         local_target_qual,
                                                         left_target_pos,
                                                         direction,
                                                         separate_ambig_counts,
                                                         right_align_ambig_dels,
                                                         right_align_ambig_ins,
                                                         max_internal_match,
                                                         min_qual,
                                                         exclude_3prime,
                                                         mutation_type,
                                                         false);
        return boost::make_tuple(included_mutations,
                                 local_effective_depth,
                                 local_effective_count);
    }

};


#endif //SHAPEMAPPER_MUTATIONCOUNTER_H
