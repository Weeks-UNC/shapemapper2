/** @file
 * @brief Classes for counting sequencing depth, sequence variants, and
 *        reverse transcription mutations.
 */

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2018 Steven Busan.                     *
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

#include "MutationProcessing.h"
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
     * @brief Resize deque on right to accommodate given alignment target location (0-based).
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

    void
    updateCounts(const std::vector<Mutation> &mutations,
                 const std::vector<bool> &local_effective_depth,
                 const std::vector<bool> &local_effective_count,
                 const int left_target_pos) {
        // mutations are in order by leftmost unchanged alignment target nuc

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

        // FIXME: make sure *_ambig tagged mutations aren't making it through to here

        // update variant counts
        for (auto &m : mutations) {
            // copy mutation object
            Mutation mut_copy(m);
            mut_copy.qual = ""; // blank out quality scores.
                                // otherwise the same variant seqs with different q-scores
                                // would be added to counts separately
            try {
                ++deq.at(mut_copy.left - target_pos).counts[mut_copy];
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
     *
     */
    void
    updateCounts(const std::vector<Mutation> &mutations,
                 const int mapping_category,
                 const int primer_pair,
                 const std::vector<bool> &mapped_depth,
                 const std::vector<bool> &local_effective_depth,
                 const std::vector<bool> &local_effective_count,
                 const int left_target_pos,
                 const bool separate_ambig_counts,
                 const bool debug = false) {

        // update histograms
        int len = local_effective_depth.size(); // FIXME: this should probably at least account for read trimming,
                                                // if not post-alignment basecall quality filtering or other effective read depth stuff
        len = std::max(len,0);
        read_lengths.count(len);
        mutations_per_read.count(mutations.size());

        // update mutation counts
        for (auto mut : mutations) {
            std::string s = mut.tag;

            // strip _ambig from mutation tag if we don't want to separate these
            if (not separate_ambig_counts) {
                int loc = s.find("_ambig");
                if (loc != std::string::npos){
                    s = s.substr(0, loc);
                }
            }

            try {
                // update mutation counts
                ++deq.at(mut.right - target_pos - 1)[s];
            }
            catch (std::out_of_range &exception) { }
            //std::cout << "updated deque count" << std::endl;

        }

        if (mapping_category == INCLUDED) {
            // update depths in deque coordinates
            for (int n = left_target_pos - target_pos;
                 n < left_target_pos + local_effective_depth.size() - target_pos;
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
        }

        // also update mapping depth columns broken down by amplicon primer pair (if used),
        // or by the reason particular reads were excluded from further processing
        std::string h = "mapped_depth";
        if (primer_pair >= 0) {
            h = "primer_pair_"+std::to_string(primer_pair+1)+"_mapped_depth";
        }
        if (mapping_category == OFF_TARGET) {
            h = "off_target_mapped_depth";
        } else if (mapping_category == LOW_MAPQ) {
            h = "low_mapq_mapped_depth";
        }

        for (int i = 0; i < mapped_depth.size(); ++i) {
            // convert index to global coordinates, then to deque coords
            int n = i + left_target_pos - target_pos;
            if (mapped_depth[i]) {
                try {
                    ++deq.at(n)[h];
                } catch (std::out_of_range &exception) { }
            }
        }

        return;
    }


};


#endif //SHAPEMAPPER_MUTATIONCOUNTER_H
