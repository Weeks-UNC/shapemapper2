/** @file
 * @brief Main interface function for counting mutations and/or depths
 *        from a mutations file.
 */

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2018 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#include "MutationCounter.h"

#include <memory>
#include <boost/iostreams/device/file.hpp>


// TODO: test scanning (sorted) file output
// TODO: option for max variant insert length to store (poor quality reads could potentially use up a fair amount of memory otherwise)

// FIXME: move mutation collapsing and quality filtering functions to MutationParser, so
//        output of that stage is actually usable for other applications

namespace mutation_counter {

    /**
     * @brief Count sequencing depth, variants, and/or mutations from parsed mutations.
     *
     * @param filename        Input file (mutations and reconstructed target
     *                        sequences). Typically from mutation_parser::parseSam()
     * @param seq_len         Length of reference sequence. Ignored if 0. If provided,
     *                        output files are guaranteed to have this many lines, even
     *                        if there are regions with no read coverage.
     * @param primer_pairs    Number of amplicon primer pairs (if any) used to filter
     *                        mapped reads.
     * @param depth_out       Sequencing depth output file. Ignored if
     *                        zero-length.
     * @param variant_out     Sequence variants and counts output file. Ignored
     *                        if zero-length.
     * @param count_out       Mutation counts output file. Ignored if zero-length.
     * @param input_is_sorted Set true if mutations came from a BAM/SAM file
     *                        sorted by left index. If true, less memory will
     *                        be used, as counts will be periodically written
     *                        to file instead of all at once at the end of
     *                        execution.
     * @param separate_ambig_counts
     *                        Count mutations derived from ambiguously-aligned mutations
     *                        in separate columns with the additional header "_ambig".
     */
    void countSelected(const std::vector<std::string> &filenames,
                       const int seq_len,
                       const int primer_pairs,
                       const std::string &variant_out,
                       const std::string &count_out,
                       const bool hist,
                       const bool input_is_sorted,
                       const bool separate_ambig_counts,
                       bool debug,
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
        new_column_names.push_back("off_target_mapped_depth");
        new_column_names.push_back("low_mapq_mapped_depth");
        if (primer_pairs > 0) {
            for (int i=1; i<=primer_pairs; i++){
                new_column_names.push_back("primer_pair_"+std::to_string(i)+"_mapped_depth");
            }
        } else {
            new_column_names.push_back("mapped_depth");
        }
        column_names = new_column_names;

        // ------------------------------------------------------------------------------
        // open input files, do some checks, set up streams

        // this garbage is needed because streams don't support move
        std::vector<std::unique_ptr<BI::filtering_istream>> files;

        for (auto & filename : filenames) {
            // again, this junk is required because streams don't support move
            files.emplace_back(
                    std::unique_ptr<BI::filtering_istream>(new BI::filtering_istream()));

            files.back()->push(BI::newline_filter(BI::newline::posix));
            try {
                int file_size = BF::file_size(filename);
                if (file_size == 0) {
                    if (warn_on_no_mapped) {
                        std::cout << "WARNING: Input file " + filename + " is empty." << std::endl;
                    } else {
                        throw std::runtime_error("ERROR: Input file " + filename + " is empty.");
                    }
                }
            } catch (BF::filesystem_error &e) { }


            BI::file_source fs(filename);
            if (not fs.is_open()) {
                // Do additional checks to see if the file exists or if it's a permissions issue
                if (!(BF::is_regular_file(filename))) {
                    throw std::runtime_error("ERROR: Input file " + filename + " not found.");
                }
                throw std::runtime_error("ERROR: Could not open input file " + filename +
                                         " - unknown error.\nCheck file and folder permissions.");
            }
            if (BF::extension(BF::path(filename)) == ".gz") {
                // decompress gzip if file looks compressed
                files.back()->push(BI::gzip_decompressor());
            }
            files.back()->push(fs);

        }

        // ------------------------------------------------------------------------------
        // open output files, set up streams

        std::vector<bool> so;
        std::vector <std::string> out_names = {variant_out,
                                               count_out};
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
        
        // ---------------------------------------------------------------------------
    
        VariantCounter vc;
        MutationCounter mc;

        if (so[1]) { *out_files[1] << mc.printHeader(); }

        std::string line;

        std::string read_id;
        int mapping_category;
        int primer_pair;
        int left_target_pos;
        int right_target_pos;
        std::vector<bool> mapping_depth;
        std::vector<bool> local_effective_depth;
        std::vector<bool> local_effective_count;
        std::vector<Mutation> processed_mutations;

        size_t count = 0;
        std::vector<std::string> input_lines(files.size());
        std::vector<bool> eof(files.size());
        std::vector<std::string> non_empty_lines;
        while(true) {
            for (int i=0; i<files.size(); i++) {
                eof[i] = !std::getline(*(files[i]), input_lines[i]);
            }

            // terminate if we reach the end of all files
            if (std::all_of(eof.begin(), eof.end(), [](bool v){return v;})) {
                break;
            }
            // interleave all inputs by line if present
            for (int i=0; i<files.size(); i++) {
                if (not eof[i]) { non_empty_lines.push_back(input_lines[i]); }
            }

            for (auto & line : non_empty_lines) {
                // process single line

                count += 1;

	            if(debug) { std::cout << "line " << count << ": " << line << std::endl << std::flush; }

                boost::tie(read_id,
                           mapping_category,
                           primer_pair,
                           left_target_pos,
                           right_target_pos,
                           mapping_depth,
                           local_effective_depth,
                           local_effective_count,
                           processed_mutations) = parseProcessedMutations(line);

                if(debug) { std::cout << "read_id: " << read_id << std::endl << std::flush; }
                if(debug) { std::cout << "mapping_category: " << mapping_category << "\n" <<std::flush; }
                if(debug) { std::cout << "primer_pair: " << std::to_string(primer_pair) << "\n" << std::flush; }
                if(debug) { std::cout << "left_target_pos: " << left_target_pos << std::endl << std::flush; }
                if(debug) { std::cout << "right_target_pos: " << right_target_pos << std::endl << std::flush; }
                if(debug) { std::cout << "mapping_depth: " << util::toString(mapping_depth) << "\n" << std::flush; }
                if(debug) { std::cout << "local_effective_depth: " << util::toString(local_effective_depth) << std::endl << std::flush; }
                if(debug) { std::cout << "local_effective_count: " << util::toString(local_effective_count) << std::endl << std::flush; }
                if(debug) { std::cout << "mutations: " << toString(processed_mutations) << std::endl << std::flush; }

                if (so[0]) { vc.updateRightBound(right_target_pos); }
                if (so[1]) { mc.updateRightBound(right_target_pos); }

                if (input_is_sorted) {
                    // use less memory if the input file is sorted. Print and remove counts to
                    // the left of the current read
                    if (so[0]) { *out_files[0] << vc.updateLeftBound(left_target_pos); }
                    if (so[1]) { *out_files[1] << mc.updateLeftBound(left_target_pos); }
                }

                if (so[0]) {
                    vc.updateCounts(processed_mutations,
                                    local_effective_depth,
                                    local_effective_count,
                                    left_target_pos);
                }
                if (so[1]) {
                    mc.updateCounts(processed_mutations,
                                    mapping_category,
                                    primer_pair,
                                    mapping_depth,
                                    local_effective_depth,
                                    local_effective_count,
                                    left_target_pos,
                                    separate_ambig_counts,
                                    debug);
                }

                /*if (count%100000==0){
                    std::cout << "Counted "<<count<<" reads" <<std::endl;
                }*/
            }
            non_empty_lines.resize(0);
            
        }
        if (count < 1) {
            if (warn_on_no_mapped) {
                std::cout << "WARNING: No reads were found in the input files." << std::endl;
            } else {
                throw std::runtime_error("ERROR: Input files contained no reads.");
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
        if (hist) {
            std::cout << mc.printHistograms();
        }
    }


}
