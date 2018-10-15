/** @file
 * @brief Count sequencing depth, sequence variants, and/or reverse
 *        transcription mutations. Commandline executable.
 */

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2018 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#include <boost/program_options.hpp>

#include "MutationCounter.cpp"

namespace po = boost::program_options;


int main(int argc, char *argv[]) {
    try {
        std::vector<std::string> in;
        int length;
        int primer_pairs;
        std::string variant_out;
        std::string count_out;
        bool hist;
        bool input_is_sorted;
        bool separate_ambig_counts;
        bool debug;
        bool warn_on_no_mapped;

        po::options_description desc("Usage");
        desc.add_options()
            ("help,h", "print usage message")

            ("in,i", po::value<std::vector<std::string>>(&in)->multitoken(), "input file path(s) (parsed mutations). comma-separated")

            ("length,n", po::value<int>(&length)->default_value(0),
             "length of reference sequence. If provided, output files are guaranteed to have this many lines even if there are regions of no read coverage.")

            ("n_primer_pairs,p", po::value<int>(&primer_pairs)->default_value(0),
             "number of primer pairs (if any) previously used for read mapping location filtering. If provided, read mapping depth columns will be split up by amplicon.")

            ("variant_out,v", po::value<std::string>(&variant_out)->default_value(""),
             "sequence variant counts output file path")

            ("count_out,c", po::value<std::string>(&count_out)->default_value(""), "mutation counts output file path")

            ("hist,h", po::bool_switch(&hist)->default_value(false), "output read length and mutation frequency histogram tables")

            ("input_is_sorted,s", po::bool_switch(&input_is_sorted)->default_value(false),
             "use less memory if scanning along mutations from a BAM file sorted by leftmost alignment position")

            ("separate_ambig_counts", po::bool_switch(&separate_ambig_counts)->default_value(false),
             "output ambiguously aligned derived mutation counts in separate columns")

            ("debug", po::bool_switch(&debug)->default_value(false),"debug info")

            ("warn_on_no_mapped,w",
              po::bool_switch(&warn_on_no_mapped)->default_value(false),
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

        if (in.size() < 1) {
            std::cerr << "ERROR: must provide at least one input file.\n";
            std::cerr << desc << std::endl;
            return 1; 
        }

        if ((variant_out + count_out ).length() == 0) {
            std::cerr << "ERROR: must include at least one output file.\n";
            std::cerr << desc << std::endl;
            return 1;
        }

        std::cout << "Attempting to count from parsed mutations file(s):";
        for (auto &filename : in) {
            std::cout << " " << filename << "\n";
        }
        if (input_is_sorted) {
            std::cout << "\t(sorted)\n";
        } else {
            std::cout << "\t(unsorted)\n";
        }
        if (length > 0){
            std::cout << " with reference sequence length " << length << '\n';
        }
        if (primer_pairs > 0) {
            std::cout << " with " << primer_pairs << " amplicon primer pairs\n";
        }
        std::cout << " and write\n";
        if (variant_out.length() != 0) {
            std::cout << "\tsequence variants and counts to " << variant_out << '\n';
        }
        if (count_out.length() != 0) {
            std::cout << "\treverse transcription mutation counts to " << count_out << '\n';
        }

        if (hist) {
            std::cout << "\tprinting read length and mutation frequency histogram tables\n";
        }

        if (separate_ambig_counts) {
            std::cout << "\toutputting ambiguous mutation counts in separate columns\n";
        }

        if (warn_on_no_mapped) {
            std::cout << "\twarning (not exiting with error) if no mapped reads in input\n";
        }
        std::cout << std::flush;

        // FIXME: move separate_ambig_counts option to MutationParser?

        mutation_counter::countSelected(in,
                                        length,
                                        primer_pairs,
                                        variant_out,
                                        count_out,
                                        hist,
                                        input_is_sorted,
                                        separate_ambig_counts,
                                        debug,
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
