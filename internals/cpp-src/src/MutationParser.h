/** @file
 * @brief Parse alignment from a BAM file.
 */

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2018 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#ifndef SHAPEMAPPER_MUTATIONPARSER_H
#define SHAPEMAPPER_MUTATIONPARSER_H

#include <bitset>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>

#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/newline.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/erase.hpp>

#include "MutationProcessing.h"

// TODO: switch over from Boost tuple to C++11 std::tuple?

namespace BF = boost::filesystem;
namespace BI = boost::iostreams;
namespace BA = boost::algorithm;

namespace mutation_parser { namespace detail {

    using namespace mutation;

    // MD operations
    enum {
        MD_DELETION, MD_MATCH, MD_MISMATCH, MD_NONE
    };
    // char types for MD tag parsing
    enum {
        CHAR_NUMERIC, CHAR_NONNUMERIC, CHAR_NONE
    };
    // for readability during ambiguous alignment handling
    enum {
        MUT_DELETION, MUT_INSERTION, MUT_MISMATCH
    };

    // These are just for debugging convenience
    const std::string MD_op_names[] = {"deletion", "match", "mismatch"};

    /** MD ops (given as a length, type, and sequence (if any):
     * MD_DELETION - deletion from reference, seq is deleted reference.
     * MD_MISMATCH - mismatch, seq is reference.
     * MD_MATCH - match, no seq provided.
     */
    struct MdOp {
        int op;
        int length;
        std::string seq;
    };

    /** See SAM spec
     */
    struct CigarOp {
        char op;
        uint32_t length;
    };

    int charType(char c) {
        if (std::isdigit(c)) {
            return CHAR_NUMERIC;
        } else {
            return CHAR_NONNUMERIC;
        }
    }

    std::string toString(MdOp m) {
        // for debugging and testing
        std::string o;
        if (m.op == MD_DELETION) {
            o += "del ";
        } else if (m.op == MD_MATCH) {
            o += "match ";
        } else if (m.op == MD_MISMATCH) {
            o += "mismatch ";
        }
        o += "length ";
        o += std::to_string(m.length);
        o += ", seq \"" + m.seq + "\"";
        return o;
    }

    std::string toString(std::vector<MdOp> m) {
        std::string o;
        for (int i = 0; i < m.size(); ++i) {
            o += toString(m[i]);
            if (i + 1 != m.size()) {
                o += ". "; // add separator except after last element
            }
        }
        return o;
    }

    /**
     * @brief Split MD tag into fields by character type.
     *        Example: "15A^GC10T30" -> {"15","A","^GC","10","T","30"}
     */
    std::vector<std::string>
    splitMDtag(const std::string &tag_contents) {
        std::vector<std::string> fields;
        int current_char_type = CHAR_NONE;
        std::string current_field;
        for (int i = 0; i < tag_contents.length(); ++i) {
            char c = tag_contents[i];
            int char_type = charType(c);
            if (char_type == current_char_type) {
                current_field.push_back(c);
            } else {
                if (current_field.length() > 0) {
                    fields.push_back(current_field);
                }
                current_field = std::string(1, c);
                current_char_type = char_type;
            }
        }
        fields.push_back(current_field);
        return fields;
    }

    /**
      * @brief Parse an MD tag into operations (match, mismatch, and deletion)
      *
      * @param tag_contents String from SAM/BAM MD tag. Typically from
      *                     BamTools::BamAlignment.GetTag("MD", returned_contents)
      * @return Vector of MdOp (holds op type, length, and sequence (if any)).
      */
    std::vector<MdOp>
    parseMDtag(const std::string &tag_contents) {

        std::vector<std::string> fields = splitMDtag(tag_contents);

        //  convert strings to match lengths as appropriate,
        //  skip mismatches of zero-length,
        //  identify deletions
        std::vector<MdOp> ops;
        for (int i = 0; i < fields.size(); ++i) {
            char c = fields[i][0];
            int char_type = charType(c);
            if (char_type == CHAR_NUMERIC) {
                int len = std::stoi(fields[i]);
                if (len != 0) {
                    ops.push_back(MdOp{MD_MATCH, len, ""});
                }
            } else if (c == '^') {
                ops.push_back(MdOp{MD_DELETION,
                                   int(fields[i].length()) - 1,
                                   fields[i].substr(1)});
            } else {
                // must be a mismatch, since deletions are handled above
                ops.push_back(MdOp{MD_MISMATCH,
                                   int(fields[i].length()),
                                   fields[i]});
            }
        }
        return ops;
    };

    /**
      * @brief  Parse a CIGAR string into operations. Only needed if parsing a SAM file.
      *
      * @param  cigar_string SAM CIGAR string. 6th field in an alignment line.
      * @return Vector of CigarOp.
      */
    std::vector<CigarOp>
    parseCIGAR(const std::string &cigar_string) {
        std::vector<CigarOp> cigar_data;
        std::vector<std::string> fields;
        std::string trimmed = boost::trim_copy(cigar_string);
        //boost::char_separator<char> ops("MIDNSHPX=");
        //boost::tokenizer<boost::char_separator<char>> tokens(trimmed, ops);
        // having trouble getting a pre-existing tokenizer to preserve delimiters,
        // so just doing it myself
        char c = cigar_string[0];
        fields.push_back(std::string(1, c));
        bool currently_numeric = false;
        if (std::isdigit(c)) {
            currently_numeric = true;
        }
        for (int i = 1; i < cigar_string.length(); ++i) {
            c = cigar_string[i];
            bool is_numeric = std::isdigit(c);
            if (is_numeric == currently_numeric) {
                fields.back() = fields.back() + std::string(1, c);
            } else {
                fields.push_back(std::string(1, c));
            }
            currently_numeric = is_numeric;
        }

        for (int i = 0; i < fields.size(); i = i + 2) {
            char op;
            uint32_t len;
            if (i + 1 >= fields.size()) {
                throw std::runtime_error("Error: CIGAR string incorrectly formatted");
            }
            try {
                op = fields[i + 1][0];
                len = std::stoi(fields[i]);
            } catch (std::exception &err) {
                throw std::runtime_error("Error: CIGAR string incorrectly formatted");
            }
            cigar_data.push_back(CigarOp{op, len});
        }
        return cigar_data;
    }

    /**
     * @brief Parse 2nd field of SAM line into bit fields
     */
    std::bitset<12>
    flagsToBits(const std::string &flags) {
        std::bitset<12> bit_flags = std::bitset<12>(std::stoi(flags));
        return bit_flags;
    }

    /**
      * @brief Determine the right-most mapped position from an alignment.
      *        Only needed if parsing a SAM file.
      *
      * @param left_target_pos  Left-most alignment target position (0-based).
      * @param cigar_data       CIGAR operations.
      * @return Right-most mapped position in alignment (0-based target coordinate).
      */
    int
    calcRightTargetPos(const int left_target_pos,
                       const std::vector<CigarOp> &cigar_data) {
        // right-most mapped position is the sum of alignment matches and deletions,
        // (inserts do not count)

        // TODO: there should be some way to avoid doing a bunch of if statements
        // boost::is_any_of("MDNP=X");

        int right_target_pos = left_target_pos;
        for (int i = 0; i < cigar_data.size(); ++i) {
            char c = cigar_data[i].op;
            if (c == 'M' or
                c == 'D' or
                c == 'N' or
                c == 'P' or
                c == '=' or
                c == 'X') {
                right_target_pos = right_target_pos + cigar_data[i].length;
            }
        }
        --right_target_pos;
        return right_target_pos;
    }


    /**
     * @brief   Find an extended SAM tag and get contents. Similar to
     *          BamTools::BamAlignment.GetTag(), but only supports string
     *          contents for now.
     * @return  Stores tag contents in destination. Returns false if not found.
     */
    bool
    getSamTag(const std::vector<std::string> &fields,
              const std::string &tag,
              std::string &destination) {
        bool found = false;
        for (int i = 11; i < fields.size(); ++i) {
            if (BA::starts_with(fields[i], tag)) {
                found = true;
                // strip tag and format specifier
                destination = fields[i].substr(5, std::string::npos);
            }
        }
        return found;
    }


    /**
     * @brief Combine info from CIGAR, MD tag, and query read to locate mutations
     *        in alignment target coordinates.
     *
     * @param pos         Left-most alignment position in target coordinates (0-based),
     *                    typically from BamTools::BamAlignment.Position()-1
     * @param query_bases Sequence read, typically from BamTools::BamAlignment.QueryBases
     * @param cigar_data  Parsed CIGAR string, typically from BamTools::BamAlignment.CigarData
     * @param md_data     Parsed MD tag, typically from mutation_parser::detail::parseMdtag()
     * @param [reconstruct_target]       Option to reconstruct the alignment target sequence
     *                                   over the alignment region.
     * @param [reconstruct_aligned_read] Option to reconstruct the aligned read sequence
     *                                   (excluding inserts) over the alignment region.
     *
     * @return Returns a tuple containing parsed mutations, reconstructed target
     *         sequence and quality scores (if requested), and reconstructed aligned read (if requested).
     */
    boost::tuple<std::vector<Mutation>, std::string, std::string, std::string, std::string>
    locateMutations(const int pos,
                    const std::string &query_bases,
                    const std::string &query_qual,
                    const std::vector<CigarOp> &cigar_data,
                    const std::vector<MdOp> &md_data,
                    bool reconstruct_target = true,
                    bool reconstruct_aligned_read = true) {
        // TODO: call variants that result in soft-clipped regions? May not really be needed.

        std::vector<Mutation> mutations;
        int ts = pos; // alignment target sequence index, i.e. "genomic" coordinate (0-based)
        int qs = 0; // query sequence index
        int mo = 0; // MD op index
        int co = 0; // CIGAR op index
        std::string target_seq; // reconstructed alignment target sequence over aligned region
        std::string target_qual;
        std::string aligned_seq; // read sequence over aligned region (includes matches, mismatches, and deletions, not insertions)
        std::string aligned_qual;

        /* CIGAR ops (given as a length and type):
         * M - seq match or mismatch
         * I - insertion to reference
         * D - deletion from reference
         * N - skipped region from reference
         * S - soft-clipped region in query sequence
         * H - hard-clipped region (not in reported query sequence)
         * P - padding (silent deletion from padded reference) - I think this means that the reference has gaps
         * = - sequence match
         * X - sequence mismatch
         *
         */

        /* MD ops (given as a length, type, and sequence (if any):
         * MD_DELETION - deletion from reference, seq is deleted reference
         * MD_MISMATCH - mismatch, seq is reference
         * MD_MATCH - match, no seq provided
         */

        // MD deletions should exactly match CIGAR deletions
        // MD mismatches and matches should exactly match CIGAR matches unless there
        //  are intervening inserts

        bool in_match = false;
        int remaining_co_length = 0;
        int temp_md_mo = -1;

        int tmp_op = -1;
        int tmp_length = 0;
        std::string tmp_target_seq;
        std::string tmp_target_qual;
        std::string tmp_query_seq;
        std::string tmp_query_qual;

        std::string s;
        std::string q;

        std::string assoc_md[cigar_data.size()];
        for (co = 0; co < cigar_data.size(); ++co) {
            char c_type = cigar_data[co].op;
            int c_length = cigar_data[co].length;


            switch (c_type) {
                case ('M'): {
                    // alignment match (can be sequence mismatch and/or match)
                    if (md_data[mo].op != MD_MATCH and md_data[mo].op != MD_MISMATCH) {
                        throw std::runtime_error(
                                "Error: MD tag does not match CIGAR string at alignment match operator ('M').");
                    }

                    // iterate through MD ops until summed MD op lengths account for cigar op length
                    if (!in_match) {
                        in_match = true;
                        remaining_co_length = 0;
                        tmp_op = md_data[mo].op;
                        tmp_length = md_data[mo].length;
                        tmp_target_seq = md_data[mo].seq;
                        tmp_target_qual = query_qual.substr(qs, tmp_length);
                        tmp_query_seq = query_bases.substr(qs, tmp_length);
                        tmp_query_qual = query_qual.substr(qs, tmp_length);
                        temp_md_mo = mo;
                    }
                    remaining_co_length += c_length;
                    while (mo < md_data.size() and
                           (md_data[mo].op == MD_MATCH or md_data[mo].op == MD_MISMATCH) and
                           remaining_co_length > 0) {

                        if (mo != temp_md_mo) {
                            tmp_op = md_data[mo].op;
                            tmp_length = md_data[mo].length;
                            tmp_target_seq = md_data[mo].seq;
                            tmp_target_qual = query_qual.substr(qs, tmp_length);
                            tmp_query_seq = query_bases.substr(qs, tmp_length);
                            tmp_query_qual = query_qual.substr(qs, tmp_length);
                            temp_md_mo = mo;
                        }

                        int overlap_length;
                        if (tmp_length > remaining_co_length) {
                            // MD op is longer than current CIGAR match op length,
                            // so take a chunk out of temp MD op and keep iterating over CIGAR ops
                            overlap_length = remaining_co_length;
                        } else {
                            // MD op is shorter than current CIGAR match op
                            // so move to next MD op
                            overlap_length = tmp_length;
                            ++mo;
                        }

                        if (tmp_op == MD_MATCH) {
                            tmp_length -= overlap_length;
                            assoc_md[co].append(MD_op_names[tmp_op] + ' ' + std::to_string(overlap_length) + ' ');
                            if (reconstruct_target or reconstruct_aligned_read) {
                                s = query_bases.substr(qs, overlap_length);
                                q = query_qual.substr(qs, overlap_length);
                            }
                            if (reconstruct_target) {
                                target_seq.append(s);
                                target_qual.append(q);
                            }
                            if (reconstruct_aligned_read) {
                                aligned_seq.append(s);
                                aligned_qual.append(q);
                            }
                        } else {
                            // MISMATCH
                            std::string target_seq_overlap = tmp_target_seq.substr(0, overlap_length);
                            std::string target_seq_remain = tmp_target_seq.substr(overlap_length);
                            std::string target_qual_overlap = tmp_target_qual.substr(0, overlap_length);
                            std::string target_qual_remain = tmp_target_qual.substr(overlap_length);
                            std::string query_seq_overlap = tmp_query_seq.substr(0, overlap_length);
                            std::string query_seq_remain = tmp_query_seq.substr(overlap_length);
                            std::string query_qual_overlap = tmp_query_qual.substr(0, overlap_length);
                            std::string query_qual_remain = tmp_query_qual.substr(overlap_length);
                            tmp_target_seq = target_seq_remain;
                            tmp_target_qual = target_qual_remain;
                            tmp_query_seq = query_seq_remain;
                            tmp_query_qual = query_qual_remain;
                            tmp_length = tmp_target_seq.length();
                            if (reconstruct_target) {
                                target_seq.append(target_seq_overlap);
                                target_qual.append(target_qual_overlap);
                            }
                            if (reconstruct_aligned_read) {
                                aligned_seq.append(query_seq_overlap);
                                aligned_qual.append(query_qual_overlap);
                            }
                            assoc_md[co].append(
                                    MD_op_names[tmp_op] + ' ' + std::to_string(overlap_length) + ' ' +
                                    query_seq_overlap +
                                    ' ');
                            mutations.push_back(Mutation(ts - 1,
                                                         ts + overlap_length,
                                                         query_seq_overlap,
                                                         query_qual_overlap));
                        }
                        ts += overlap_length;
                        qs += overlap_length;
                        remaining_co_length -= overlap_length;
                    }

                    if (remaining_co_length == 0 and tmp_length == 0) {
                        in_match = false; // reached a point where cigar matches and md tag matches/mismatches sum to the same length
                    }

                    break;
                }
                case ('I'): {
                    // insertion
                    mutations.push_back(Mutation(ts - 1,
                                                 ts,
                                                 query_bases.substr(qs, c_length),
                                                 query_qual.substr(qs, c_length)));
                    qs += c_length;
                    break;
                }
                case ('D'): {
                    // deletion
                    if (md_data[mo].op != MD_DELETION or md_data[mo].length != c_length) {
                        throw std::runtime_error(
                                "Error: MD tag does not match CIGAR string at deletion operator ('D').");
                    }

                    mutations.push_back(Mutation(ts - 1,
                                                 ts + c_length,
                                                 "",
                                                 ""));
                    if (reconstruct_target) {
                        target_seq.append(md_data[mo].seq);
                        for (int i = 0; i < md_data[mo].length; i++) {
                            target_qual.append("!");
                        }
                    }
                    if (reconstruct_aligned_read) {
                        s = std::string(md_data[mo].length, '-');
                        aligned_seq.append(s);
                        for (int i = 0; i < md_data[mo].length; i++) {
                            aligned_qual.append("!");
                        }
                    }
                    ts += c_length;
                    assoc_md[co].append(MD_op_names[md_data[mo].op] + ' ' + std::to_string(md_data[mo].length) + ' ' +
                                        md_data[mo].seq + ' ');
                    ++mo;
                    break;
                }
                case ('N'): {
                    // skipped region from reference sequence
                    if (reconstruct_target or reconstruct_aligned_read) {
                        s = std::string(c_length, '~');
                    }
                    if (reconstruct_target) {
                        target_seq.append(s);
                        for (int i = 0; i < md_data[mo].length; i++) {
                            target_qual.append("!");
                        }
                    }
                    if (reconstruct_aligned_read) {
                        aligned_seq.append(s);
                        for (int i = 0; i < md_data[mo].length; i++) {
                            aligned_qual.append("!");
                        }
                    }
                    // TODO: treat skipped ref regions as deletions? CIGAR "N"
                    ts += 1;
                    break;
                }
                case ('S'): {
                    // soft-clipped region
                    qs += c_length;
                    break;
                }
                case ('H'): {
                    // hard-clipped region (not present in reference sequence)
                    break;
                }
                case ('P'): {
                    // padding (i.e. '-' chars in reference sequence)
                    if (reconstruct_target or reconstruct_aligned_read) {
                        s = std::string(c_length, '-');
                    }
                    if (reconstruct_target) {
                        target_seq.append(s);
                        target_qual.append(q);
                    }
                    if (reconstruct_aligned_read) {
                        aligned_seq.append(s);
                        aligned_qual.append(q);
                    }
                    ts += c_length;
                    break;
                }
                case ('='): {
                    // explicit sequence match
                    if (md_data[mo].op != MD_MATCH or md_data[mo].length != c_length) {
                        throw std::runtime_error(
                                "Error: MD tag does not match CIGAR string at explicit match operator ('=').");
                    }
                    if (reconstruct_target or reconstruct_aligned_read) {
                        s = query_bases.substr(qs, c_length);
                        q = query_qual.substr(qs, c_length);
                    }
                    if (reconstruct_target) {
                        target_seq.append(s);
                        target_qual.append(q);
                    }
                    if (reconstruct_aligned_read) {
                        aligned_seq.append(s);
                        aligned_qual.append(q);
                    }
                    ++mo;
                    break;
                }
                case ('X'): {
                    // explicit sequence mismatch
                    if (md_data[mo].op != MD_MISMATCH or md_data[mo].length != c_length) {
                        throw std::runtime_error(
                                "Error: MD tag does not match CIGAR string at explicit mismatch operator ('X').");
                    }
                    mutations.push_back(Mutation(ts - 1,
                                                 ts + c_length,
                                                 md_data[mo].seq,
                                                 query_qual.substr(qs, c_length)));
                    if (reconstruct_target) {
                        target_seq.append(md_data[mo].seq);
                        target_qual.append(query_qual.substr(qs, c_length));
                    }
                    if (reconstruct_aligned_read) {
                        aligned_seq.append(query_bases.substr(qs, c_length));
                        aligned_qual.append(query_qual.substr(qs, c_length));
                    }
                    qs += c_length;
                    ts += c_length;
                    ++mo;
                    break;
                }
                default: {
                    throw std::runtime_error("Error: Malformed CIGAR string.");
                }
            }

        }

        // DEBUG
        /*for (int i = 0; i < cigar_data.size(); ++i) {
            char l = '(';
            char r = ')';
            if (i == co) {
                l = '<';
                r = '>';
            }
            std::cout << l << cigar_data[i].length << cigar_data[i].op << r << "-[";
            std::cout << assoc_md[i] << "] ";
        }
        //std::cout <<'\n'<< target;
        std::cout << '\n';*/

        return boost::make_tuple(mutations, target_seq, target_qual, aligned_seq, aligned_qual);
    }

    /**
     * @brief Locate mutations without reconstructing alignment target sequence or
     *        aligned read sequence.
     */
    std::vector<Mutation>
    locateMutationsNoReconstruct(const int pos,
                                 const std::string &query_bases,
                                 const std::string &query_qual,
                                 const std::vector<CigarOp> &cigar_data,
                                 const std::vector<MdOp> &md_data) {

        return locateMutations(pos, query_bases, query_qual, cigar_data, md_data).get<0>();
    }


    void printBoolVector(std::vector<bool> v) {
        std::string p;
        for (int i = 0; i < v.size(); ++i) {
            p += v[i] ? '1' : '0';
        }
        p += '\n';
        std::cout << p;
    }

    // FIXME: add docstring
    void slideIndel(const std::string &local_target_seq,
                    const std::string &local_target_qual,
                    const std::string &aligned_seq,
                    const std::string &aligned_qual,
                    const std::vector<bool> &has_insert_left_of,
                    const int local_left,
                    const int local_right,
                    std::string &mut_seq,
                    std::string &mut_qual,
                    const int mut_type,
                    const int dir,
                    std::vector<Mutation> &adjusted_mutations,
                    std::vector<std::vector<int>> &appended_target_indices,
                    const int mut_index) {

        int offset = 0;
        while (true) {
            ++offset;
            int offset_left = local_left + offset * dir;
            int offset_right = local_right + offset * dir;

            // stop if out of bounds
            // TODO: handle insert one off end?
            // TODO: loop getting a bit large - break into smaller functions
            try {
                local_target_seq.at(offset_left);
                local_target_seq.at(offset_right);
            } catch (std::out_of_range) {
                break;
            }

            // stop sliding if we run into a gap
            if (mut_type == MUT_DELETION) {
                if ((dir == 1 and aligned_seq[offset_right] == '-') or
                    (dir == -1 and aligned_seq[offset_left] == '-')) {
                    break;
                }
            } else {
                // MUT_INSERTION
                if ((dir == 1 and aligned_seq[offset_left] == '-') or
                    (dir == -1 and aligned_seq[offset_right] == '-')) {
                    break;
                }
            }

            // stop sliding if we run into an insertion
            if (mut_type == MUT_DELETION) {
                if ((dir == 1 and has_insert_left_of[offset_right]) or
                    (dir == -1 and offset_left - 1 > 0 and has_insert_left_of[offset_left - 1])) {
                    break;
                }
            } else {
                // MUT_INSERTION
                if (dir == 1 and has_insert_left_of[offset_right]) {
                    break;
                }
            }


            // mut_seq should drop one char on left,
            // pick up one char from alignment target seq on right if del, left if ins
            // (or in the other direction analogously)
            // note: storing target index of picked up nuc to allow later joining of indels that picked up the same nuc(s)
            char dropped_seq;
            char dropped_qual;
            char from_target_seq;
            char from_target_qual;
            char from_aligned_seq;
            char from_aligned_qual;
            int from_target_index;
            if (dir == 1) {
                dropped_seq = mut_seq[0];
                dropped_qual = mut_qual[0];
                if (mut_type == MUT_DELETION) {
                    from_target_seq = local_target_seq[offset_right];
                    from_target_qual = local_target_qual[offset_right];
                    from_target_index = offset_right;
                    from_aligned_seq = aligned_seq[offset_right];
                    from_aligned_qual = aligned_qual[offset_right];
                } else {
                    // MUT_INSERTION
                    from_target_seq = local_target_seq[offset_left];
                    from_target_qual = local_target_qual[offset_left];
                    from_target_index = offset_left;
                    from_aligned_seq = aligned_seq[offset_left];
                    from_aligned_qual = aligned_qual[offset_left];
                }
                mut_seq = mut_seq.substr(1) + from_target_seq;
                mut_qual = mut_qual.substr(1) + from_target_qual;
                //appended_target_indices[mut_index].push_back(from_target_index);
            } else {
                //dir==-1
                dropped_seq = mut_seq[mut_seq.length() - 1];
                dropped_qual = mut_qual[mut_qual.length() - 1];
                if (mut_type == MUT_DELETION) {
                    from_target_seq = local_target_seq[offset_left];
                    from_target_qual = local_target_qual[offset_left];
                    from_target_index = offset_left;
                    from_aligned_seq = aligned_seq[offset_left];
                    from_aligned_qual = aligned_qual[offset_left];
                } else {
                    // MUT_INSERTION
                    from_target_seq = local_target_seq[offset_right];
                    from_target_qual = local_target_qual[offset_right];
                    from_target_index = offset_right;
                    from_aligned_seq = aligned_seq[offset_right];
                    from_aligned_qual = aligned_qual[offset_right];
                }
                mut_seq = from_target_seq + mut_seq.substr(0, mut_seq.length() - 1);
                mut_qual = from_target_qual + mut_qual.substr(0, mut_qual.length() - 1);
                //appended_target_indices[mut_index].insert(appended_target_indices[mut_index].begin(),
                //                                                          from_target_index);
            }

            if (dropped_seq == from_target_seq) {
                if (dir == 1) {
                    // update replacement seq
                    adjusted_mutations[mut_index].seq += from_aligned_seq;
                    adjusted_mutations[mut_index].qual += from_aligned_qual;
                    appended_target_indices[mut_index].push_back(from_target_index);
                    ++adjusted_mutations[mut_index].right;
                } else {
                    //dir==-1
                    adjusted_mutations[mut_index].seq = from_aligned_seq + adjusted_mutations[mut_index].seq;
                    adjusted_mutations[mut_index].qual = from_aligned_qual + adjusted_mutations[mut_index].qual;
                    appended_target_indices[mut_index].insert(appended_target_indices[mut_index].begin(),
                                                              from_target_index);
                    --adjusted_mutations[mut_index].left;
                }
            } else {
                break;
            }
        }
    }

    std::vector<Mutation>
    identifyAmbiguousMutations(const int pos,
                               const std::string &local_target_seq,
                               const std::string &local_target_qual,
                               const std::string &aligned_seq,
                               const std::string &aligned_qual,
                               const std::vector<Mutation> &mutations) {

        if (local_target_seq.length() != aligned_seq.length()) {
            throw std::runtime_error("Error: target sequence and aligned sequence lengths do not match.");
        }

        // store insert locations in linear array for marginally faster collision checking
        // when "sliding" gaps and inserts
        std::vector<bool> has_insert_left_of(aligned_seq.length() + 1, false);
        {
            int k = -1;
            for (std::vector<Mutation>::const_iterator it = mutations.begin();
                 it != mutations.end();
                 ++it) {
                ++k;
                if (!(*it).isSimpleInsert()) {
                    continue;
                }

                has_insert_left_of[(*it).right - pos] = true;
            }
        }

        std::vector<Mutation> adjusted_mutations(mutations.size());
        std::vector<std::vector<int>> appended_target_indices(
                mutations.size()); // helper array to hold temporary information about
        // nuc indices appended to ambiguous adjacent indels
        {
            int k = -1;
            for (std::vector<Mutation>::const_iterator it = mutations.begin();
                 it != mutations.end();
                 ++it) {
                ++k;
                adjusted_mutations[k] = Mutation((*it));
                if (not ((*it).isSimpleInsert() or (*it).isSimpleGap())) {
                    // MUT_MISMATCH
                    // skip these now, but before returning need to remove any mismatches
                    // that were merged with ambiguous gaps or inserts
                    continue;
                }

                int d[] = {1, -1}; // directions to "slide" indel (right and left, respectively)
                for (int i = 0; i < 2; ++i) {
                    int dir = d[i];

                    std::string mut_seq;
                    std::string mut_qual;
                    int local_left;
                    int local_right;
                    int mut_type;
                    if ((*it).isSimpleInsert()) {
                        mut_type = MUT_INSERTION;
                        local_left = (*it).left - pos; // 0-base target char left of insert
                        local_right = (*it).right - pos; // 0-based target char right of insert
                        mut_seq = (*it).seq;
                        mut_qual = (*it).qual;
                    } else if ((*it).isSimpleGap()) {
                        mut_type = MUT_DELETION;
                        local_left = (*it).left - pos + 1; // 0-based leftmost deleted char in target
                        local_right = (*it).right - pos - 1; // 0-based rightmost deleted char in target (inclusive)
                        mut_seq = local_target_seq.substr(local_left, local_right - local_left + 1);
                        mut_qual = local_target_qual.substr(local_left, local_right - local_left + 1);
                    }

                    //std::cout << "initial mut_seq: " << mut_seq << std::endl;
                    slideIndel(local_target_seq,
                               local_target_qual,
                               aligned_seq,
                               aligned_qual,
                               has_insert_left_of,
                               local_left,
                               local_right,
                               mut_seq,
                               mut_qual,
                               mut_type,
                               dir,
                               adjusted_mutations,
                               appended_target_indices,
                               k);

                }
            }
        }

        /*
            for (int i=0; i < adjusted_mutations.size(); ++i) {
                std::string s;
                s = std::to_string(adjusted_mutations[i].left) + ' ';
                s += std::to_string(adjusted_mutations[i].right) + " \"";
                s += adjusted_mutations[i].seq + "\" ";
                for (int j=0; j < appended_target_indices[i].size(); ++j) {
                    s += std::to_string(appended_target_indices[i][j]) + ' ';
                }
                std::cout << s << " | " << std::flush;
            }
            std::cout << '\n';
            */

        // store ambiguous gap and insert locations in linear form
        // for marginally faster removal of merged mismatches
        std::vector<bool> indel_covered(aligned_seq.length(), false);
        {
            for (std::vector<Mutation>::iterator it = adjusted_mutations.begin();
                 it != adjusted_mutations.end();
                 ++it) {
                if ((*it).isGapOrInsert()) {
                    for (int i = (*it).left + 1 - pos; i < (*it).right - pos; ++i) {
                        indel_covered[i] = true;
                    }
                } else {
                    continue;
                }
            }
        }


        // create a new vector with no duplicated mutations (since mismatches can be picked up
        // by adjacent ambiguously-aligned gaps or inserts)
        // - chain together ambiguously-aligned mutations that picked up the same adjacent nuc(s)
        std::vector<Mutation> merged_removed;
        for (int i = 0;
             i < adjusted_mutations.size();
             ++i) {
            if (adjusted_mutations[i].isGapOrInsert()) {
                // check mutation on left for shared nucs picked up during indel sliding
                // - merge mutations and eliminate shared nucs if needed
                // - (This is only needed to take care of rare corner cases in
                //    bad alignments)
                if (appended_target_indices[i].size() < 1 ||
                    merged_removed.size() < 1 ||
                    appended_target_indices[i - 1].size() < 1) {
                    merged_removed.push_back(adjusted_mutations[i]);
                    continue;
                } else {
                    bool do_merge = false;
                    // eliminate shared nucs
                    std::string tmp_seq = adjusted_mutations[i].seq;
                    std::string tmp_qual = adjusted_mutations[i].qual;
                    for (int k = 0; k < appended_target_indices[i].size(); ++k) {
                        // check if nuc shared with mutation on left
                        // if appended_target_indices[i][k] in appended_target_indices[i-1]:
                        if (std::find(appended_target_indices[i - 1].begin(),
                                      appended_target_indices[i - 1].end(),
                                      appended_target_indices[i][k]) != appended_target_indices[i - 1].end()) {
                            tmp_seq = tmp_seq.substr(1, tmp_seq.length());
                            tmp_qual = tmp_qual.substr(1, tmp_qual.length());
                            do_merge = true;
                        } else {
                            break;
                        }
                    }
                    // append replacement seq and update rightmost unchanged nuc index
                    if (do_merge) {
                        merged_removed.back().seq += tmp_seq;
                        merged_removed.back().qual += tmp_qual;
                        merged_removed.back().right = adjusted_mutations[i].right;
                    } else {
                        merged_removed.push_back(adjusted_mutations[i]);
                        continue;
                    }
                }
            } else {
                // is a mismatch
                bool previously_merged = false;
                for (int k = adjusted_mutations[i].left + 1 - pos;
                     k < adjusted_mutations[i].right - pos; ++k) {
                    if (indel_covered[k]) {
                        previously_merged = true;
                        break;
                    }
                }
                if (!previously_merged) {
                    merged_removed.push_back(adjusted_mutations[i]);
                }
            }
        }

        /*if (merged_removed.size() != adjusted_mutations.size()){
                std::cout << "Duplicated mutations in ambig handling func" << std::endl;
            std::cout << local_target << std::endl;
            std::cout << aligned << std::endl;
                for (int i=0; i < adjusted_mutations.size(); ++i) {
                    std::string s;
                    s = std::to_string(adjusted_mutations[i].left) + ' ';
                    s += std::to_string(adjusted_mutations[i].right) + " \"";
                    s += adjusted_mutations[i].seq + "\" ";
                    for (int j=0; j < appended_target_indices[i].size(); ++j) {
                        s += std::to_string(appended_target_indices[i][j]) + ' ';
                    }
                    std::cout << s << " | " << std::flush;
                }
                std::cout << '\n';
                for (int i=0; i < merged_removed.size(); ++i) {
                    std::string s;
                    s = std::to_string(merged_removed[i].left) + ' ';
                    s += std::to_string(merged_removed[i].right) + " \"";
                    s += merged_removed[i].seq + "\" ";
                    //for (int j=0; j < appended_target_indices[i].size(); ++j) {
                    //    s += std::to_string(appended_target_indices[i][j]) + ' ';
                    //}
                    std::cout << s << " | " << std::flush;
                }
                std::cout << '\n';

            }*/

        adjusted_mutations = merged_removed;

        return adjusted_mutations;
    }


    /**
     * @brief Parse a single BAM/SAM alignment.
     *
     * @param left_target_pos  Leftmost alignment target position covered by
     *                         a read (0-based target coord), typically from
     *                         BamTools::BamAlignment.Position
     * @param right_target_pos Rightmost alignment target position covered,
     *                         typically from
     *                         BamTools::BamAlignment.GetEndPosition()-1
     * @param query_bases      Read sequence, typically from
     *                         BamTools::BamAlignment.QueryBases
     * @param cigar_data       Parsed CIGAR string, typically from
     *                         BamTools::BamAlignment.cigar_data
     * @param md_tag_contents  MD tag, typically from
     *                         BamTools::BamAlignment.GetTag("MD", ...)
     *
     * @return Returns vector containing: Left-most position in alignment target
     *         sequence (0-based), right-most position (0-based), reconstructed
     *         target sequence over alignment region, and mutations with respect
     *         to the target sequence. Mutations are also adjusted so that
     *         ambiguously aligned indels include information about alternate
     *         valid placements (see identifyAmbiguousMutations()). Now using
     *         Read class.
     */
    Read
    parseMutations(const int left_target_pos,
                   const int right_target_pos,
                   const std::string &query_bases,
                   const std::string &query_qual,
                   const std::vector<CigarOp> &cigar_data,
                   const std::string &md_tag_contents) {

        std::vector<detail::MdOp> parsed_MD = parseMDtag(md_tag_contents);

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
                                                         parsed_MD,
                                                         true);

        Read read;
        read = Read(left_target_pos,
                    right_target_pos,
                    local_target_seq)
                .setQual(local_target_qual)
                .setMutations(mutations);

        if (debug_out) {
            debug_out << "parsed mutations from SAM read\n";
            debug_out << read;
            debug_out << std::flush;
        }

        //std::cout << local_target << std::endl;
        //std::cout << aligned_query << std::endl;

        // identify ambiguously aligned indels locations and adjust mutations
        // to include that information
        // TODO: option not to perform this calculation? this might slow things down
        //       for long and/or noisy reads
        std::vector<Mutation> adjusted_mutations = identifyAmbiguousMutations(left_target_pos,
                                                                              local_target_seq,
                                                                              local_target_qual,
                                                                              aligned_query_seq,
                                                                              aligned_query_qual,
                                                                              mutations);
        read.setMutations(adjusted_mutations);

        if (debug_out) {
            debug_out << "identified ambiguously aligned mutations\n";
            debug_out << read;
            debug_out << std::flush;
        }

        // note: returning left_target_pos and right_target_pos is redundant, but convenient
        return read;
    };




}

}


#endif //SHAPEMAPPER_MUTATIONPARSER_H
