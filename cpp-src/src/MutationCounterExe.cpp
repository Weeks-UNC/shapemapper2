/** @file
 * @brief Count sequencing depth, sequence variants, and/or reverse
 *        transcription mutations. Commandline executable.
 */

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2017 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#include <boost/program_options.hpp>

#include "MutationCounter.cpp"

namespace po = boost::program_options;


int main(int argc, char *argv[]) {
    try {

        po::options_description desc("Usage");
        desc.add_options()
                ("help,h", "print usage message")

                ("in,i", po::value<std::string>()->required(), "input file path (parsed mutations)")

                ("length,n", po::value<int>()->default_value(0),
                 "length of reference sequence. If provided, output files are guaranteed to have this many lines even if there are regions of no read coverage.")

                ("variant_out,v", po::value<std::string>()->default_value(""),
                 "sequence variant counts output file path")

                ("count_out,c", po::value<std::string>()->default_value(""), "mutation counts output file path")

                ("classified_out,l", po::value<std::string>()->default_value(""),
                 "mutations and classifications debug output file path")

                ("exclude_3prime", po::value<int>()->default_value(0),
                 "exclude mutations occurring within this many nucleotides of 3-prime end of read")
                 // TODO: add option to exclude 5-prime as well (could be useful for people pooling multiple targeted primer pairs)
                 // TODO: add explicit filters for mapping requirements? e.g. left end within specific target nuc range, right end etc.

                ("input_is_sorted,s", po::bool_switch()->default_value(false),
                 "use less memory if scanning along mutations from a BAM file sorted by leftmost alignment position")

                ("right_align_ambig_dels", po::bool_switch()->default_value(false),
                 "realign ambiguously aligned deletions to right end (not recommended), otherwise realign left")

                ("right_align_ambig_ins", po::bool_switch()->default_value(false),
                 "realign ambiguously aligned insertions to right end (not recommended), otherwise realign left")

                ("max_internal_match", po::value<int>()->default_value(7),
                 "allow up to N unchanged reference sequence nucs between merged mutations")

                ("min_qual", po::value<int>()->default_value(30),
                 "Exclude mutations that contain or are adjacent to any basecalls with Phred quality scores below this value. This filter is also applied to the calculation of the effective read depth.")

                ("separate_ambig_counts", po::bool_switch()->default_value(false),
                 "output ambiguously aligned derived mutation counts in separate columns")

                ("use_only_mutation_type", po::value<std::string>()->default_value(""),
                 "use only mutations from a specific mutation class (not recommended). Possible values: mismatch gap insert gap_multi insert_multi complex")

                ("warn_on_no_mapped,w",
                  po::bool_switch()->default_value(false),
                  "exit with warning instead of error if no mapped reads present in input")
                ;



        po::variables_map vm;

        try {

            po::store(po::parse_command_line(argc, argv, desc), vm);

            if (vm.count("help") or argc == 1) {
                std::cout << desc << std::endl;
                return 0;
            }
            po::notify(vm);
        }
        catch (const po::error &e) {
            std::cerr << "ERROR: " << e.what() << "\n" << std::endl;
            std::cerr << desc << std::endl;
            return 1;
        }
        catch (const std::exception &e) {
            std::cerr << "ERROR: " << e.what() << "\n" << std::endl;
            return 1;
        }

        std::string in = vm["in"].as<std::string>();
        std::string variant_out = vm["variant_out"].as<std::string>();
        std::string count_out = vm["count_out"].as<std::string>();
        std::string classified_out = vm["classified_out"].as<std::string>();
        bool input_is_sorted = vm["input_is_sorted"].as<bool>();
        bool separate_ambig_counts = vm["separate_ambig_counts"].as<bool>();
        bool right_align_ambig_dels = vm["right_align_ambig_dels"].as<bool>();
        bool right_align_ambig_ins = vm["right_align_ambig_ins"].as<bool>();
        int max_internal_match = vm["max_internal_match"].as<int>();
        int min_qual = vm["min_qual"].as<int>();
        int exclude_3prime = vm["exclude_3prime"].as<int>();
        int seq_len = vm["length"].as<int>();
        std::string mutation_type = vm["use_only_mutation_type"].as<std::string>();
        bool warn_on_no_mapped = vm["warn_on_no_mapped"].as<bool>();

        if ((variant_out + count_out + classified_out).length() == 0) {
            std::cerr << "ERROR: must include at least one output file.\n";
            std::cout << desc << std::endl;
            return 1;
        }

        std::cout << "Attempting to count from parsed mutations file " << in;
        if (input_is_sorted) {
            std::cout << " (sorted)";
        } else {
            std::cout << " (unsorted)";
        }
        if (seq_len > 0){
            std::cout << " with reference sequence length " << seq_len << '\n';
        }
        std::cout << " and write\n";
        if (variant_out.length() != 0) {
            std::cout << "\tsequence variants and counts to " << variant_out << '\n';
        }
        if (count_out.length() != 0) {
            std::cout << "\treverse transcription mutation counts to " << count_out << '\n';
        }
        if (classified_out.length() != 0) {
            std::cout << "\tmutation classification debug info to " << classified_out << '\n';
        }

        if (right_align_ambig_dels){
            std::cout << "\ttreating ambiguously aligned deletions as right-aligned\n";
        } else {
            std::cout << "\ttreating ambiguously aligned deletions as left-aligned\n";
        }
        if (right_align_ambig_ins) {
            std::cout << "\ttreating ambiguously aligned insertions as right-aligned\n";
        } else {
            std::cout << "\ttreating ambiguously aligned insertions as left-aligned\n";
        }

        if (mutation_type.length() != 0) {
            std::cout << "\tusing only mutations of the type: " << mutation_type << '\n';
        }

        if (exclude_3prime > 0){
            std::cout << "\texcluding mutations within " << exclude_3prime << " nucleotides of read 3-prime end\n";
        }
        std::cout << "\tmerging adjacent mutations within " << max_internal_match << " nucleotides of each other\n";
        std::cout << "\texcluding mutations with any basecall q-scores below " << min_qual << "\n";

        if (separate_ambig_counts) {
            std::cout << "\toutputing ambiguous mutation counts in separate columns\n";
        }

        if (warn_on_no_mapped) {
            std::cout << "\twarning (not exiting with error) if no mapped reads in input\n";
        }
        

        mutation_counter::countSelected(in,
                                        seq_len,
                                        variant_out,
                                        count_out,
                                        classified_out,
                                        input_is_sorted,
                                        separate_ambig_counts,
                                        right_align_ambig_dels,
                                        right_align_ambig_ins,
                                        max_internal_match,
                                        min_qual,
                                        exclude_3prime,
                                        mutation_type,
                                        warn_on_no_mapped);

        std::cout << "... Successfully counted mutations." << std::endl;
    }
    catch (const BF::filesystem_error &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Unknown error." << std::endl;
        return 1;
    }
    return 0; //SUCCESS
}
