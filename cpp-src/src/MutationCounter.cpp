/** @file
 * @brief Main interface function for counting mutations and/or depths
 *        from a mutations file.
 */

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2017 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#include "MutationCounter.h"

#include <memory>


// TODO: test scanning (sorted) file output
// TODO: option for max variant insert length to store (poor quality reads could potentially use up a fair amount of memory otherwise)

namespace mutation_counter {

    /**
     * @brief Count sequencing depth, variants, and/or mutations from parsed mutations.
     *
     * @param filename        Input file (mutations and reconstructed target
     *                        sequences). Typically from mutation_parser::parseBam()
     * @param seq_len         Length of reference sequence. Ignored if 0. If provided,
     *                        output files are guaranteed to have this many lines, even
     *                        if there are regions with no read coverage.
     * @param depth_out       Sequencing depth output file. Ignored if
     *                        zero-length.
     * @param variant_out     Sequence variants and counts output file. Ignored
     *                        if zero-length.
     * @param count_out       Mutation counts output file. Ignored if zero-length.
     * @param classified_out  Mutation classification debug output file. Ignored if
     *                        zero-length.
     * @param input_is_sorted Set true if mutations came from a BAM/SAM file
     *                        sorted by left index. If true, less memory will
     *                        be used, as counts will be periodically written
     *                        to file instead of all at once at the end of
     *                        execution.
     * @param separate_ambig_counts
     *                        Count mutations derived from ambiguously-aligned mutations
     *                        in separate columns with the additional header "_ambig".
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
     */
    void countSelected(const std::string &filename,
                       const int seq_len,
                       const std::string &variant_out,
                       const std::string &count_out,
                       const std::string &classified_out,
                       const bool input_is_sorted,
                       const bool separate_ambig_counts,
                       const bool right_align_ambig_dels,
                       const bool right_align_ambig_ins,
                       const int max_internal_match,
                       const int min_qual,
                       const int exclude_3prime,
                       const std::string &mutation_type,
                       bool warn_on_no_mapped = false) {
        // set selected output columns (this is a global from MutationCounter.h)
        std::vector<std::string> new_column_names;
        for (std::vector<std::string>::iterator it = column_names.begin();
             it != column_names.end();
             ++it) {
            new_column_names.push_back((*it));
            if (separate_ambig_counts){
                new_column_names.push_back((*it)+"_ambig");
            }
        }
        new_column_names.push_back("read_depth");
        new_column_names.push_back("effective_depth");
        column_names = new_column_names;

        try {
            int file_size = BF::file_size(filename);
            if (file_size == 0) {
                if (warn_on_no_mapped) {
                    std::cout << "WARNING: Input file " + filename + " is empty." << std::endl;
                } else {
                    throw std::runtime_error("ERROR: Input file " + filename + " is empty.");
                }
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

        std::vector<bool> so;
        std::vector <std::string> out_names = {variant_out,
                                               count_out,
                                               classified_out};
        std::vector <std::shared_ptr<std::ofstream>> out_files;
        for (int i = 0; i < out_names.size(); ++i) {
            if (out_names[i].length() == 0) {
                out_files.emplace_back(new std::ofstream());
                so.push_back(false);
                continue;
            }
            //std::ofstream *file_out = new std::ofstream(out_names[i], std::ios_base::out | std::ios_base::binary);
            std::shared_ptr <std::ofstream> file_out(
                    new std::ofstream(out_names[i], std::ios_base::out | std::ios_base::binary));

            if (!(*file_out)) {
                throw std::runtime_error(
                        "ERROR: Could not open output file " + out_names[i] + "\nCheck file and folder permissions.");
            }
            std::unique_ptr <BI::filtering_ostream> out(new BI::filtering_ostream);
            if (BF::extension(BF::path(out_names[i])) == ".gz") {
                // compress using gzip if requested
                out->push(BI::gzip_compressor());
            }
            out->push((*file_out));
            out_files.push_back(file_out);
            so.push_back(true);
        }

        VariantCounter vc;
        MutationCounter mc;

        if (so[1]) { *out_files[1] << mc.printHeader(); }

        std::string line;
        size_t c = 0;

        std::string read_id;
        int left_target_pos;
        int right_target_pos;
        std::string local_target_seq;
        std::string local_target_qual;
        std::vector <Mutation> adjusted_mutations;

        // output intermediate per-read mutation info (mostly for debugging)
        std::vector <TaggedMutation> classified_mutations;
        std::vector <bool> local_effective_depth;
        std::vector <bool> local_effective_count;

        int count = 0;
        while (std::getline(in, line)) {
            count += 1;
	        //std::cout << "line " << count << ": " << line << std::endl;
            boost::tie(read_id,
                       left_target_pos,
                       right_target_pos,
                       local_target_seq,
                       local_target_qual,
                       adjusted_mutations) = parseReadInfo(line);
            if (so[0]) { vc.updateRightBound(right_target_pos); }
            if (so[1]) { mc.updateRightBound(right_target_pos); }

            if (input_is_sorted) {
                // use less memory if the input file is sorted. Print and remove counts to
                // the left of the current read
                if (so[0]) { *out_files[0] << vc.updateLeftBound(left_target_pos); }
                if (so[1]) { *out_files[1] << mc.updateLeftBound(left_target_pos); }
            }

            if (so[0]) {
                vc.updateCounts(adjusted_mutations,
                                local_target_seq,
                                local_target_qual,
                                left_target_pos,
                                min_qual,
                                exclude_3prime,
                                mutation_type);
            }
            if (so[1]) {
                boost::tie(classified_mutations,
                           local_effective_depth,
                           local_effective_count) =
                        mc.updateCounts(read_id,
                                        adjusted_mutations,
                                        local_target_seq,
                                        local_target_qual,
                                        left_target_pos,
                                        1,
                                        separate_ambig_counts,
                                        right_align_ambig_dels,
                                        right_align_ambig_ins,
                                        max_internal_match,
                                        min_qual,
                                        mutation_type,
                                        exclude_3prime);
            }

            // output a local effective mutation count string too (bit redundant but would be
            // simpler to parse for applications that don't care about the specific details of
            // each mutation)
            // NOTE: left and right are 0-based inclusive
            if (so[2]) {
                *out_files[2] << read_id << ' '
                              << left_target_pos << ' '
                              << right_target_pos - exclude_3prime << ' '
                              << toString(local_effective_depth) << ' '
                              << toString(local_effective_count) << ' '
                              << toString(classified_mutations) << std::endl;
            }

            c++;
            /*if (c%100000==0){
                std::cout << "Counted "<<c<<" reads" <<std::endl;
            }*/
	        //std::cout << "made it to end of loop" << std::endl;
        }
        if (c < 1) {
            if (warn_on_no_mapped) {
                std::cout << "WARNING: Input file " + filename + " contains no reads." << std::endl;
            } else {
                throw std::runtime_error("ERROR: Input file " + filename + " contains no reads.");
            }
        }

        // update right end to last nucleotide in sequence
        // (ensures that output files will have the correct
        //  number of lines even if there are no reads that cover
        //  the 3-prime end of the sequence).
        if (seq_len>0) {
            // (bounds are given in 0-based coordinates)
            if (so[0]) { vc.updateRightBound(seq_len-1); }
            if (so[1]) { mc.updateRightBound(seq_len-1); }
        }

        // print any remaining values in each deque
        if (so[0]) { *out_files[0] << vc.printAllValues(); }
        if (so[1]) { *out_files[1] << mc.printAllValues(); }

        // write read length and mutations per read histograms
        // TODO: make separate output files for these tables
        std::cout << mc.printHistograms();
    }


}
