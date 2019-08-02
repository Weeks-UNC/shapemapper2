/** @file
 * @brief Main interface functions for parsing SAM file or single alignments
 *        into mutations and reconstructed target sequences.
 */

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2018 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#include "MutationParser.h"


namespace mutation_parser {
    using namespace mutation;

    /**
     * @brief Given a single SAM alignment entry (assumed mapped and split into fields) -
     *        parse mutations, reconstruct target sequence over aligned region,
     *        and merge ambiguously aligned indels with alternate placements.
     *
     * @param fields Single SAM alignment, already split by tab character into vector
     * @return Returns the information needed to count RT mutations and depths or
     * count variants: Left-most position in alignment target sequence (0-based),
     * right-most position (0-based), reconstructed target sequence over this range,
     * and mutations present in this read. Also return read id for debugging purposes.
     * Now in a Read instance.
     */
    Read
    parseSamFields(const std::vector <std::string> &fields,
                   const int min_mapq,
                   const bool input_is_unpaired) {
        // TODO: add tests for reads mapped in the reverse sense

        if (fields.size() < 11) {
            throw std::runtime_error("Error: unable to parse incomplete line.");
        }

        int left_target_pos;
        try {
            left_target_pos = stoi(fields[3]) - 1;
        } catch (std::exception &err) {
            throw std::runtime_error("Error: line is incorrectly formatted (couldn't parse mapped location).");
        }

        std::string read_id = fields[0];
        std::string query_bases = fields[9];
        std::string query_qual = fields[10];
        std::string cigar_string = fields[5];

        int strand = FORWARD;
        std::bitset<12> flags = detail::flagsToBits(fields[1]);
        // 5th bit indicates read mapped to reverse strand
        if (flags[4]) {
            strand = REVERSE;
        }

        int read_type = MERGED;
        if (input_is_unpaired) {
            read_type = UNPAIRED;
            // command line flag "input_is_unpaired" distinguishes unpaired reads from mixed paired/unmerged paired reads
            // - STAR reports read 1's ID for both reads, so can't depend on *:R2 read id to parse
        } else if (flags[6]) { read_type = PAIRED_R1; }
        else if (flags[7]) { read_type = PAIRED_R2; }

        std::vector <detail::CigarOp> cigar_data = detail::parseCIGAR(cigar_string);

        int right_target_pos = detail::calcRightTargetPos(left_target_pos,
                                                          cigar_data);

        const std::string tag = "MD";
        std::string md_tag_contents;
        if (!detail::getSamTag(fields, tag, md_tag_contents)) {
            throw std::runtime_error("Error: no MD tag in alignment.");
        }

        if (debug_out) {
            debug_out << read_id << "\n"
            << std::flush;
        }


        int mapping_category = INCLUDED;
        if (fields[2] == "*") {
            if (debug_out) {
                debug_out << "read is unmapped\n" << std::flush;
            }
            mapping_category = UNMAPPED;
        }
        else if (std::stoi(fields[4]) < min_mapq) {
            // exclude reads below mapping quality threshold
            if (debug_out) {
                debug_out << "read has low mapping quality ("
                << fields[4] << ")\n"
                << std::flush;
            }
            mapping_category = LOW_MAPQ;
        }

        Read r{}; // call default constructor
        if (mapping_category == INCLUDED) {
            // only do these parsing steps for mapped reads with a decent MAPQ
            r = detail::parseMutations(left_target_pos,
                                       right_target_pos,
                                       query_bases,
                                       query_qual,
                                       cigar_data,
                                       md_tag_contents);
        } else if (mapping_category == LOW_MAPQ) {
            r.setLeft(left_target_pos)
             .setRight(right_target_pos);
        }
        r.setId(read_id)
         .setStrand(strand)
         .setReadType(read_type)
         .setMappingCategory(mapping_category);

        r.mapped_depth.resize(r.right - r.left + 1, 1);

        if (debug_out) {
            debug_out << "identified read type and mapped depth\n"
            << r
            << std::flush;
        }

        return r;
    }

    Read
    parseSamLine(const std::string &line,
                 const int min_mapq,
                 const bool input_is_unpaired) {

        // split into fields
        std::vector <std::string> fields;
        std::string trimmed = boost::trim_copy(line);
        boost::split(fields, trimmed, boost::is_any_of("\t"), boost::token_compress_off);

        return parseSamFields(fields, min_mapq, input_is_unpaired);
    }

    boost::tuple<int, int>
    findClosestPrimers(const int left,
                       const int right,
                       const std::vector <PrimerPair> &primer_pairs,
                       const int max_primer_offset) {
        // identify nearest primers from mapped read end locations

        int fw_primer_index = NO_ASSOCIATED_PRIMER_PAIR;
        int rv_primer_index = NO_ASSOCIATED_PRIMER_PAIR;
        std::vector<int> fw_primer_indices; // store indices of forward primers within max_primer_offset
        std::vector<int> fw_dists;
        std::vector<int> rv_primer_indices;
        std::vector<int> rv_dists;

        int fw_dist;
        int rv_dist;
        for (int i = 0; i < primer_pairs.size(); i++) {
            fw_dist = abs(left - primer_pairs[i].fw_left);
            if (fw_dist <= max_primer_offset) {
                fw_primer_indices.push_back(i);
                fw_dists.push_back(fw_dist);
            }
            rv_dist = abs(right - primer_pairs[i].rv_right);
            if (rv_dist <= max_primer_offset) {
                rv_primer_indices.push_back(i);
                rv_dists.push_back(rv_dist);
            }
        }
        int min_dist = 999999;
        for (int i = 0; i < fw_dists.size(); i++) {
            if (fw_dists[i] < min_dist) {
                fw_primer_index = fw_primer_indices[i];
                min_dist = fw_dists[i];
            }
        }
        min_dist = 999999;
        for (int i = 0; i < rv_dists.size(); i++) {
            if (rv_dists[i] < min_dist) {
                rv_primer_index = rv_primer_indices[i];
                min_dist = rv_dists[i];
            }
        }

        if (debug_out) {
            debug_out << "identified nearest primers to mapped location, out of "
            << primer_pairs.size() << " total pairs\n"
            << std::flush;
        }

        if (debug_out) { debug_out << "fw_primer_index: " << fw_primer_index << "\n"; }
        if (debug_out) { debug_out << "rv_primer_index: " << rv_primer_index << "\n"; }
        return boost::make_tuple(fw_primer_index, rv_primer_index);
    }

    int
    findOverlappingPrimers(const int left,
                           const int right,
                           const std::vector <PrimerPair> &primer_pairs) {
        /* find primer pair overlapping read ends*/
        int fw_primer_index = NO_ASSOCIATED_PRIMER_PAIR;
        int rv_primer_index = NO_ASSOCIATED_PRIMER_PAIR;
        for (int i = 0; i < primer_pairs.size(); i++) {
            if (left >= primer_pairs[i].fw_left and
                left <= primer_pairs[i].fw_right) {
                fw_primer_index = i;
            }
            if (right >= primer_pairs[i].rv_left and
                right <= primer_pairs[i].rv_right) {
                rv_primer_index = i;
            }
        }
        int primer_index = std::max(fw_primer_index, rv_primer_index);
        return primer_index;
    }

    bool
    isConcordant(const std::vector <Read> &reads,
                 const int max_paired_fragment_length) {
        // Recompute whether concordant alignment, since STAR has no --maxins param,
        // and may generate pairs mapped very far apart. Low mapping quality may also cause
        // us to reject one or the other reads.
        // FIXME: test
        bool concordant = true;
        if (reads[R1].mapping_category == UNMAPPED or reads[R2].mapping_category == UNMAPPED) {
            concordant = false;
        } else {
            // expect mapping to opposite strands
            if (reads[R1].strand == reads[R2].strand) {
                concordant = false;
            } else {
                // check fragment size
                int fragment_length = std::max(reads[R1].right,
                                               reads[R2].right)
                                      - std::min(reads[R1].left,
                                                 reads[R2].left);
                if (fragment_length > max_paired_fragment_length) {
                    concordant = false;
                }
                // check for dovetailing
                if (reads[R1].strand == FORWARD) {
                    if (reads[R2].left < reads[R1].left and
                        reads[R1].right > reads[R2].right) {
                        concordant = false;
                    }
                } else {
                    if (reads[R1].left < reads[R2].left and
                        reads[R2].right > reads[R1].right) {
                        concordant = false;
                    }
                }
            }
        }
        return concordant;
    }

    // FIXME: could probably unify and simplify these functions for read
    //        mapping location requirements
    bool
    isOffTargetUnpairedRead(const int fw_primer_index,
                            const int rv_primer_index,
                            const bool require_forward_primer_mapped,
                            const bool require_reverse_primer_mapped,
                            const int max_primer_offset) {
        bool off_target = false;
        // skip read if mapped locations don't match primer pairs
        if (require_forward_primer_mapped and fw_primer_index < 0) {
            if (debug_out) {
                debug_out << "skipped read because mapped end not within +/- "
                << max_primer_offset << " nts (inclusive) of required amplicon forward primer\n"
                << std::flush;
            }
            off_target = true;
        }
        else if (require_reverse_primer_mapped and rv_primer_index < 0) {
            if (debug_out) {
                debug_out << "skipped read because mapped end not within +/- "
                << max_primer_offset << " nts (inclusive) of required amplicon reverse primer\n"
                << std::flush;
            }
            off_target = true;
        }
        else if (require_forward_primer_mapped and
                 require_reverse_primer_mapped) {
            if (fw_primer_index != rv_primer_index) {
                if (debug_out) {
                    debug_out << "skipped read because mapped ends not near a matched amplicon primer pair\n"
                    << std::flush;
                }
                off_target = true;
            }
        }
        return off_target;
    }

    bool
    isOffTargetPairedRead(const int fw_primer_index,
                          const int rv_primer_index,
                          const bool require_forward_primer_mapped,
                          const bool require_reverse_primer_mapped,
                          const int max_primer_offset) {
        bool off_target = false;
        if (require_forward_primer_mapped and
            require_reverse_primer_mapped) {
            if (fw_primer_index != rv_primer_index) {
                if (debug_out) {
                    debug_out << "skipped read because mapped ends not near a matched amplicon primer pair\n"
                              << std::flush;
                }
                off_target = true;
            }
        } else if (require_forward_primer_mapped and fw_primer_index < 0) {
            if (debug_out) {
                debug_out << "skipped read because mapped end not within +/- "
                << max_primer_offset << " nts (inclusive) of required amplicon forward primer\n"
                << std::flush;
            }
            off_target = true;
        } else if (require_reverse_primer_mapped and rv_primer_index < 0) {
            if (debug_out) {
                debug_out << "skipped read because mapped end not within +/- "
                << max_primer_offset << " nts (inclusive) of required amplicon reverse primer\n"
                << std::flush;
            }
            off_target = true;
        }
        return off_target;
    }

    std::string
    parseUnpairedRead(const std::string &line,
                      const int min_mapq,
                      const bool right_align_ambig_dels,
                      const bool right_align_ambig_ins,
                      const int max_internal_match,
                      const int min_qual,
                      const int exclude_3prime,
                      const std::string &mutation_type,
                      const bool variant_mode,
                      const std::vector <PrimerPair> &primer_pairs,
                      const bool trim_primers,
                      const bool require_forward_primer_mapped,
                      const bool require_reverse_primer_mapped,
                      const int max_primer_offset,
                      const bool debug) {
        // parse single read (merged, unpaired, or paired but missing mate)

        if (debug_out) {
            debug_out << "[separator] -------------------------------------------------------------------------------\n"
            << std::flush;
        }


        if (debug) {
            std::cout << "in parseUnpairedReads()" << "\n" << std::flush;
        }


        Read read = parseSamLine(line, min_mapq, false);

        // skip unmapped reads
        if (read.mapping_category == UNMAPPED) {
            return "";
        }

        // skip reads below mapping quality threshold
        if (read.mapping_category == LOW_MAPQ) {
            return read.serializeMutations();
        }

        int fw_primer_index, rv_primer_index;
        boost::tie(fw_primer_index,
                   rv_primer_index) =
                findClosestPrimers(read.left,
                                   read.right,
                                   primer_pairs,
                                   max_primer_offset);

        bool off_target;
        if ((read.read_type == UNPAIRED) or (read.read_type == MERGED)) {
            off_target = isOffTargetUnpairedRead(fw_primer_index,
                                                  rv_primer_index,
                                                  require_forward_primer_mapped,
                                                  require_reverse_primer_mapped,
                                                  max_primer_offset);
        } else {
            // do not require both primer sites if this is, for example, a good
            // R1 with a missing R2 because of run quality
            bool rfp = require_forward_primer_mapped;
            bool rrp = require_reverse_primer_mapped;
            if (read.strand == FORWARD) {
                rrp = false;
            } else {
                rfp = false;
            }
            bool off_target = isOffTargetPairedRead(fw_primer_index,
                                                    rv_primer_index,
                                                    rfp,
                                                    rrp,
                                                    max_primer_offset);
        }

        /*if (debug_out and off_target) {
            debug_out << "OFF_TARGET\n" << std::flush;
        }*/
        if (off_target) {
            read.setMappingCategory(OFF_TARGET);
            return read.serializeMutations();
        }

        int primer_index = std::max(fw_primer_index,
                                    rv_primer_index);
        read.setPrimerPair(primer_index);

        // if no primer pair found within max_primer_offset, relax the search
        // and just find any primer that overlaps either end for trimming purposes
        // FIXME: this probably doesn't correctly handle multiple overlapping primers
        if (primer_index == NO_ASSOCIATED_PRIMER_PAIR) {
            primer_index = findOverlappingPrimers(read.left,
                                                  read.right,
                                                  primer_pairs);
        }

        PrimerPair primer_pair;
        if (primer_index >= 0) {
            primer_pair = PrimerPair(primer_pairs.at(primer_index));
        }


        //if (debug_out) {
        //    debug_out << line << "\n";
        //}

        std::vector <Read> read_list;
        read_list.push_back(read);
        Read processed_read =
                processMutations(read_list,
                                 FORWARD,
                                 right_align_ambig_dels,
                                 right_align_ambig_ins,
                                 max_internal_match,
                                 min_qual,
                                 exclude_3prime,
                                 mutation_type,
                                 variant_mode,
                                 trim_primers,
                                 primer_pair,
                                 debug);

        //if (debug_out) { debug_out << "[separator]\n"; }

        processed_read.setReadType(read.read_type)
                .setMappingCategory(read.mapping_category)
                .setMappedDepth(read.mapped_depth)
                .setPrimerPair(read.primer_pair);
        // Note: using primer pair matched within +/- max_primer_offset,
        //       not relaxed overlap search only used for trimming in rare cases
        std::string s = processed_read.serializeMutations();
        if (debug) {
            std::cout << "processed_read.serializeMutations(): " << s << "\n" << std::flush;
        }
        return s;
    }

    std::string
    parsePairedReads(const std::vector <std::string> &lines,
                     const int max_paired_fragment_length,
                     const int min_mapq,
                     const bool right_align_ambig_dels,
                     const bool right_align_ambig_ins,
                     const int max_internal_match,
                     const int min_qual,
                     const int exclude_3prime,
                     const std::string &mutation_type,
                     const bool variant_mode,
                     const std::vector <PrimerPair> &primer_pairs,
                     const bool trim_primers,
                     const bool require_forward_primer_mapped,
                     const bool require_reverse_primer_mapped,
                     const int max_primer_offset,
                     const bool debug) {
        if (debug_out) {
            debug_out << "[separator] ##############################################################################\n"
            << std::flush;
        }

        if (debug) {
            std::cout << "in parsePairedReads()" << "\n" << std::flush;
        }

        // parse group of two reads R1 and R2
        std::string s = "";
        bool concordant = true;
        std::vector <Read> reads;

        for (int i = 0; i < 2; i++) {
            Read read = parseSamLine(lines[i],
                                     min_mapq,
                                     true);
            reads.push_back(read);
        }

        if (debug) {
            for (auto r: reads) {
                std::cout << r << std::flush;
            }
        }

        if (reads[R1].mapping_category == UNMAPPED and
            reads[R2].mapping_category == UNMAPPED) {
            return "";
        }

        int fw_read_index = R1;
        int rv_read_index = R2;
        if (reads[R1].strand == REVERSE and reads[R2].strand == FORWARD) {
            fw_read_index = R2;
            rv_read_index = R1;
        }

        if (debug) {
            std::cout << "fw_read_index: " << fw_read_index << "\n" << std::flush;
            std::cout << "rv_read_index: " << rv_read_index << "\n" << std::flush;
        }


        concordant = isConcordant(reads,
                                  max_paired_fragment_length);

        if (debug) {
            std::cout << "concordant: " << concordant << "\n" << std::flush;
        }

        // if both reads are LOW_MAPQ, make a simple merged read and return
        // if one read is INCLUDED but the other read is LOW_MAPQ or UNMAPPED,
        // set concordant to False to process reads separately
        if (reads[R1].mapping_category == LOW_MAPQ and reads[R2].mapping_category == LOW_MAPQ) {
            if (debug) {
                std::cout << "both reads are LOW_MAPQ\n" << std::flush;
            }
            Read simple_merged = mergeMatePairsSimple(reads);
            simple_merged.setMappingCategory(LOW_MAPQ);
            if (debug) {
                std::cout << "returned from margeMatePairsSimple()\n" << std::flush;
            }
            return simple_merged.serializeMutations();
        } else {
            int included_count = 0;
            for (int i = 0; i < 2; i++) {
                if (reads[i].mapping_category == INCLUDED) {
                    included_count++;
                }
            }
            if (included_count != 2) {
                concordant = false;
            }
        }

        if (debug) {
            std::cout << "concordant: " << concordant << "\n" << std::flush;
        }

        //if (debug) { std::cout << "reads.size(): " << reads.size() << "\n"; }

        std::vector <PrimerPair> matching_primer_pairs;
        // dummy primer pair if no match found or no primers provided
        PrimerPair primer_pair;
        int primer_index = NO_ASSOCIATED_PRIMER_PAIR;
        matching_primer_pairs.push_back(primer_pair);

        // WIP: clean up, move to function (lot of duplicated code)

        matching_primer_pairs.resize(0);

        if (concordant) {
            // identify primer pair and filter reads if mapping conditions not met

            int fw_primer_index, rv_primer_index;
            boost::tie(fw_primer_index,
                       rv_primer_index) =
                    findClosestPrimers(reads[fw_read_index].left,
                                       reads[rv_read_index].right,
                                       primer_pairs,
                                       max_primer_offset);
            bool off_target = isOffTargetPairedRead(fw_primer_index,
                                                    rv_primer_index,
                                                    require_forward_primer_mapped,
                                                    require_reverse_primer_mapped,
                                                    max_primer_offset);
            /*if (debug_out and off_target) {
                debug_out << "OFF_TARGET\n" << std::flush;
            }*/

            if (off_target) {
                Read simple_merged = mergeMatePairsSimple(reads);
                simple_merged.setMappingCategory(OFF_TARGET);
                return simple_merged.serializeMutations();
            }

            int primer_index = std::max(fw_primer_index,
                                        rv_primer_index);
            reads[R1].setPrimerPair(primer_index);
            reads[R2].setPrimerPair(primer_index);

            if (debug) {
                std::cout << "set both reads primer_index to: " << primer_index << "\n" << std::flush;
            }

            if (primer_index == NO_ASSOCIATED_PRIMER_PAIR) {
                primer_index = findOverlappingPrimers(reads[fw_read_index].left,
                                                      reads[rv_read_index].right,
                                                      primer_pairs);
                if (debug) {
                    std::cout << "relaxed primer_index: " << primer_index << "\n" << std::flush;
                }
            }
            PrimerPair primer_pair{};
            if (primer_index >= 0) {
                primer_pair = PrimerPair(primer_pairs.at(primer_index));
            }

            matching_primer_pairs.push_back(primer_pair);

            if (debug) {
                std::cout << "primer_index: " << primer_index << "\n" << std::flush;
            }

        } else {
            // 1 or 2 discordantly mapped reads
            std::vector <Read> new_reads;
            for (int i = 0; i < reads.size(); i++) {
                int fw_primer_index, rv_primer_index;
                boost::tie(fw_primer_index,
                           rv_primer_index) =
                        findClosestPrimers(reads[i].left,
                                           reads[i].right,
                                           primer_pairs,
                                           max_primer_offset);
                // do not require both primer sites if this is, for example, a good
                // R1 with a missing R2 because of run quality
                bool rfp = require_forward_primer_mapped;
                bool rrp = require_reverse_primer_mapped;
                if (reads[i].strand == FORWARD) {
                    rrp = false;
                } else {
                    rfp = false;
                }
                bool off_target = isOffTargetPairedRead(fw_primer_index,
                                                        rv_primer_index,
                                                        rfp,
                                                        rrp,
                                                        max_primer_offset);

                /*if (debug_out and off_target) {
                    debug_out << "OFF_TARGET\n" << std::flush;
                }*/

                if (off_target) {
                    reads[i].setMappingCategory(OFF_TARGET);
                }

                int primer_index = std::max(fw_primer_index,
                                            rv_primer_index);
                reads[i].setPrimerPair(primer_index);

                // if no primer pair found within max_primer_offset, relax the search
                // and just find any primer that overlaps either end for trimming purposes
                if (primer_index == NO_ASSOCIATED_PRIMER_PAIR) {
                    primer_index = findOverlappingPrimers(reads[i].left,
                                                          reads[i].right,
                                                          primer_pairs);
                }

                PrimerPair primer_pair;
                if (primer_index >= 0) {
                    primer_pair = PrimerPair(primer_pairs.at(primer_index));
                }
                matching_primer_pairs.push_back(primer_pair);
                new_reads.push_back(reads[i]);
            }
            reads = new_reads;
        }

        // bit ugly, since R1/R2 alone need separate calls
        if (concordant) {
            Read processed_read =
                    processMutations(reads,
                                     FORWARD,
                                     right_align_ambig_dels,
                                     right_align_ambig_ins,
                                     max_internal_match,
                                     min_qual,
                                     exclude_3prime,
                                     mutation_type,
                                     variant_mode,
                                     trim_primers,
                                     matching_primer_pairs[0],
                                     debug);

            //if (debug_out) { debug_out << "[separator]\n"; }
            processed_read.setReadType(PAIRED);
            s += processed_read.serializeMutations();

        } else {
            // paired reads mapped separately
            for (int i = 0; i < reads.size(); i++) {
                // update read type to reflect discordantly mapped if needed
                if (reads[i].read_type == PAIRED_R1) {
                    reads[i].read_type = UNPAIRED_R1;
                } else if (reads[i].read_type == PAIRED_R2) {
                    reads[i].read_type = UNPAIRED_R2;
                }

                if (reads[i].mapping_category == UNMAPPED) {
                    continue;
                }

                if (reads[i].mapping_category != INCLUDED) {
                    s += reads[i].serializeMutations();
                    continue;
                }

                std::vector <Read> read_list;
                read_list.push_back(reads[i]);
                Read processed_read =
                        processMutations(read_list,
                                         FORWARD,
                                         right_align_ambig_dels,
                                         right_align_ambig_ins,
                                         max_internal_match,
                                         min_qual,
                                         exclude_3prime,
                                         mutation_type,
                                         variant_mode,
                                         trim_primers,
                                         matching_primer_pairs[i],
                                         debug);

                s += processed_read.serializeMutations();
            }
        }
        return s;
    }

    // FIXME: just make a class to pass these names around
    boost::tuple<bool, bool, bool, bool, bool, bool, bool, bool>
    getReadMappingProperties(std::string line) {
        /*bool paired,
             concordant,
             unmapped,
             mate_unmapped,
             reverse_strand,
             mate_reverse_strand,
             first_in_pair,
             second_in_pair;*/

        // split into fields
        std::vector <std::string> fields;
        std::string trimmed = boost::trim_copy(line);
        boost::split(fields, trimmed, boost::is_any_of("\t"), boost::token_compress_off);

        if (fields.size() < 11) {
            throw std::runtime_error("Error: unable to parse incomplete line.");
        }

        /*
        The SAM FLAGS field, the second field in a SAM record, has multiple bits that describe the paired-end nature
        of the read and alignment. The first (least significant) bit (1 in decimal, 0x1 in hexadecimal) is set if the
        read is part of a pair. The second bit (2 in decimal, 0x2 in hexadecimal) is set if the read is part of a pair
        that aligned in a paired-end fashion. The fourth bit (8 in decimal, 0x8 in hexadecimal) is set if the read is
        part of a pair and the other mate in the pair had at least one valid alignment. The sixth bit (32 in decimal,
        0x20 in hexadecimal) is set if the read is part of a pair and the other mate in the pair aligned to the
        Crick strand (or, equivalently, if the reverse complement of the other mate aligned to the Watson strand).
        The seventh bit (64 in decimal, 0x40 in hexadecimal) is set if the read is mate 1 in a pair. The eighth bit
        (128 in decimal, 0x80 in hexadecimal) is set if the read is mate 2 in a pair.
        */

        std::bitset<12> flags = detail::flagsToBits(fields[1]);

        return boost::make_tuple(flags[0],
                                 flags[1],
                                 flags[2],
                                 flags[3],
                                 flags[4],
                                 flags[5],
                                 flags[6],
                                 flags[7]);
    }


    /**
     * @brief Parse mutations and reconstruct target sequences for mapped reads in
     *        a given SAM file. Process mutations to ... FIXME: complete description
     *
     * @param paired          Parse SAM file as pairs of R1,R2 reads
     * @param min_mapq
     *                        Minimum mapping quality to include read

     * @param right_align_ambig_dels
     *                        Realign ambiguously placed deletions to their rightmost
     *                        valid position if true, leftmost if false.
     * @param right_align_ambig_ins
     *                        Realign ambiguously placed insertions to their rightmost
     *                        valid position if true, leftmost if false.
     * @param max_internal_match
     *                        Combine nearby mutations separated by up to this many
     *                        unchanged reference nucleotides.
     * @param min_qual        Minimum basecall quality Phred score to allow in a mutation
     *                        before excluding from counting.
     * @param exclude_3prime  Exclude any mutations overlapping the region within this
     *                        many nucleotides of the right end of a read.
     *
     */
    void parseSAM(const std::string &filename,
                  const std::string &outname,
                  const std::string &debug_outname,
                  const std::string &primers_filename,
                  const int max_paired_fragment_length,
                  const unsigned int min_mapq,
                  const bool right_align_ambig_dels,
                  const bool right_align_ambig_ins,
                  const int max_internal_match,
                  const int min_qual,
                  const int exclude_3prime,
                  const std::string mutation_type,
                  const bool variant_mode,
                  const bool trim_primers,
                  const bool require_forward_primer_mapped,
                  const bool require_reverse_primer_mapped,
                  const int max_primer_offset,
                  const bool input_is_unpaired,
                  const bool debug,
                  const bool warn_on_no_mapped = false) {

        std::vector <PrimerPair> primer_pairs;
        if (primers_filename != "") {
            primer_pairs = loadPrimerPairs(primers_filename);
        }

        try {
            int file_size = BF::file_size(filename);
            if (file_size == 0) {
                throw std::runtime_error("ERROR: Input file " + filename + " is empty.");
            }
        } catch (BF::filesystem_error &e) {
            // handled below
        }

        std::ifstream file_in(filename, std::ios_base::in | std::ios_base::binary);
        if (!file_in) {
            // Do additional checks to see if the file exists or if it's a permissions issue
            if (!(BF::is_regular_file(filename))) {
                throw std::runtime_error("ERROR: Input file " + filename + " not found.");
            }
            throw std::runtime_error("ERROR: Could not open input file " + filename +
                                     " - unknown error.\nCheck file and folder permissions.");
        }

        BI::filtering_istream in;
        // universal newline support filter
        in.push(BI::newline_filter(BI::newline::posix));
        if (BF::extension(BF::path(filename)) == ".gz") {
            // decompress gzip if file looks compressed
            in.push(BI::gzip_decompressor());
        }
        in.push(file_in);

        std::ofstream file_out(outname, std::ios_base::out | std::ios_base::binary);
        if (!file_out) {
            throw std::runtime_error(
                    "ERROR: Could not open output file " + outname + "\nCheck file and folder permissions.");
        }

        BI::filtering_ostream out;
        if (BF::extension(BF::path(outname)) == ".gz") {
            // compress using gzip if requested
            out.push(BI::gzip_compressor());
        }
        out.push(file_out);

        // init debug_out if filename provided
        if (debug_outname.size() > 0) {
            mutation::debug_out.open(debug_outname, std::ios_base::out | std::ios_base::binary);
            if (!mutation::debug_out) {
                throw std::runtime_error(
                        "ERROR: Could not open debug output file " + debug_outname +
                        "\nCheck file and folder permissions.");
            }
        }

        std::string line;
        size_t c = 0;

        std::vector <std::string> lines; // store R1 as first element and R2 as second element if present
        while (std::getline(in, line)) {
            // skip headers
            if (line.length() < 1 or line[0] == '@') {
                continue;
            }

            bool paired,
                 concordant,
                 unmapped,
                 mate_unmapped,
                 reverse_strand,
                 mate_reverse_strand,
                 first_in_pair,
                 second_in_pair;
            boost::tie(paired,
                       concordant,
                       unmapped,
                       mate_unmapped,
                       reverse_strand,
                       mate_reverse_strand,
                       first_in_pair,
                       second_in_pair) = getReadMappingProperties(line);

            lines.push_back(line);

            if ((not input_is_unpaired) and
                (not mate_unmapped) and
                concordant and
                lines.size() == 2) {

                if (debug) {
                    for (auto &l : lines) {
                        std::cout << l << std::endl << std::flush;
                    }
                }


                /*const std::vector<std::string> &lines,
                     const int max_paired_fragment_length,
                     const int min_mapq,
                     const bool right_align_ambig_dels,
                     const bool right_align_ambig_ins,
                     const int max_internal_match,
                     const int min_qual,
                     const bool exclude_3prime,
                     const std::string &mutation_type,
                     const bool variant_mode,
                      const std::vector<PrimerPair> &primer_pairs,
                      const bool trim_primers,
                      const bool require_forward_primer_mapped,
                      const bool require_reverse_primer_mapped,
                      const int max_primer_offset,
                     const bool debug*/
                std::string s = parsePairedReads(lines,
                                                 max_paired_fragment_length,
                                                 min_mapq,
                                                 right_align_ambig_dels,
                                                 right_align_ambig_ins,
                                                 max_internal_match,
                                                 min_qual,
                                                 exclude_3prime,
                                                 mutation_type,
                                                 variant_mode,
                                                 primer_pairs,
                                                 trim_primers,
                                                 require_forward_primer_mapped,
                                                 require_reverse_primer_mapped,
                                                 max_primer_offset,
                                                 debug);

                out << s;
                out << std::flush; // FIXME: remove if possible?
                c++; // FIXME: don't increment for unmapped reads
                lines.clear();
            } else {
                if (debug) {
                    std::cout << lines[0] << std::endl << std::flush;
                }
                if (debug) { std::cout << "exclude_3prime just before parseUnpairedRead: " << exclude_3prime << "\n"; }

                /*const std::string &line,
                      const int min_mapq,
                      const bool right_align_ambig_dels,
                      const bool right_align_ambig_ins,
                      const int max_internal_match,
                      const int min_qual,
                      const bool exclude_3prime,
                      const std::string &mutation_type,
                      const bool variant_mode,
                      const std::vector<PrimerPair> &primer_pairs,
                      const bool trim_primers,
                      const bool require_forward_primer_mapped,
                      const bool require_reverse_primer_mapped,
                      const int max_primer_offset,
                      const bool debug*/
                std::string s = parseUnpairedRead(lines[0],
                                                  min_mapq,
                                                  right_align_ambig_dels,
                                                  right_align_ambig_ins,
                                                  max_internal_match,
                                                  min_qual,
                                                  exclude_3prime,
                                                  mutation_type,
                                                  variant_mode,
                                                  primer_pairs,
                                                  trim_primers,
                                                  require_forward_primer_mapped,
                                                  require_reverse_primer_mapped,
                                                  max_primer_offset,
                                                  debug);

                out << s;
                out << std::flush; // FIXME: remove if possible
                c++;
                lines.clear();
            }
        }

        out << std::flush;

        if (c < 1) {
            if (warn_on_no_mapped) {
                std::cout << "WARNING: Input file " + filename + " contains no mapped reads." << std::endl;
            } else {
                throw std::runtime_error("ERROR: Input file " + filename + " contains no mapped reads.");
            }
        }

    }


}
