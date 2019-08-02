/** @file
 * @brief Parse mapped SAM alignments into mutations and
 *        reconstructed alignment target sequences. Commandline executable.
 */

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2018 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "MutationParser.cpp"


namespace po = boost::program_options;
namespace BA = boost::algorithm;

int main(int argc, char *argv[]) {
    try {
        std::string in;
        std::string out;
        std::string debug_out;
        std::string primers;
        bool input_is_unpaired;
        int max_paired_fragment_length;
        int min_mapq;
        int exclude_3prime;
        bool right_align_ambig_dels;
        bool right_align_ambig_ins;
        int max_internal_match;
        int min_qual;
        std::string use_only_mutation_type;
        bool warn_on_no_mapped;
        bool variant_mode;
        bool trim_primers;
        bool require_forward_primer_mapped;
        bool require_reverse_primer_mapped;
        int max_primer_offset;
        bool debug;

        po::options_description desc("Usage");
        desc.add_options()
            ("help,h", "print usage message")

            ("in,i", po::value<std::string>(&in)->required(), "SAM input file path")

            ("out,o", po::value<std::string>(&out)->required(), "parsed mutations output file path")

            ("debug_out,d", po::value<std::string>(&debug_out)->default_value(""), "intermediate debug info file path")

            ("max_paired_fragment_length",
             po::value<int>(&max_paired_fragment_length)->default_value(800),
             "analogous to bowtie2's --maxins param. Paired reads mapping to a fragment size above this threshold will "
             "be treated as separate reads.")

            ("min_mapq,m",
             po::value<int>(&min_mapq)->default_value(30),
             "minimum reported mapping quality to allow")

            ("exclude_3prime", po::value<int>(&exclude_3prime)->default_value(0),
             "exclude mutations occurring within this many nucleotides of 3-prime end of read")

            // new flag distinguishing merged from unpaired input
            ("input_is_unpaired", po::bool_switch(&input_is_unpaired)->default_value(false),
            "specify that reads are unpaired (as opposed to paired and/or unmerged paired reads)")

            /*
            params: not sure if simpler to pass and parse input file or define CLI args
                primers,p: filename with primers and locations to parse
               trim_primers
               require_forward_primer_mapped
               require_reverse_primer_mapped
               max_primer_offset
            */

            ("primers", po::value<std::string>(&primers)->default_value(""),
            "")
            ("trim_primers", po::bool_switch(&trim_primers)->default_value(false),
            "")
            ("require_forward_primer_mapped", po::bool_switch(&require_forward_primer_mapped)->default_value(false),
            "")
            ("require_reverse_primer_mapped", po::bool_switch(&require_reverse_primer_mapped)->default_value(false),
            "")
            ("max_primer_offset", po::value<int>(&max_primer_offset)->default_value(0),
            "")

            ("right_align_ambig_dels", po::bool_switch(&right_align_ambig_dels)->default_value(false),
             "realign ambiguously aligned deletions to right end (not recommended), otherwise realign left")

            ("right_align_ambig_ins", po::bool_switch(&right_align_ambig_ins)->default_value(false),
             "realign ambiguously aligned insertions to right end (not recommended), otherwise realign left")

            ("max_internal_match", po::value<int>(&max_internal_match)->default_value(7),
             "allow up to N unchanged reference sequence nucs between merged mutations")

            ("min_qual", po::value<int>(&min_qual)->default_value(30),
             "Exclude mutations that contain or are adjacent to any basecalls with Phred quality scores below this value. This filter is also applied to the calculation of the effective read depth.")

            ("use_only_mutation_type", po::value<std::string>(&use_only_mutation_type)->default_value(""),
             "use only mutations from a specific mutation class (not recommended). Possible values: mismatch gap insert gap_multi insert_multi complex")

            ("variant_mode,v",  po::bool_switch(&variant_mode)->default_value(false),
            "If true, nearby mutation merging and ambiguous mutation realignment steps will not be performed. Used by shapemapper to simplify sequence variant detection, i.e. SNP calling.")

            ("debug", po::bool_switch(&debug)->default_value(false),
            "print debugging information")

            ("warn_on_no_mapped,w",
              po::bool_switch(&warn_on_no_mapped)->default_value(false),
              "exit with warning instead of error if no mapped reads present in input")
             ;

        po::variables_map vm;

        try {

            po::store(po::parse_command_line(argc, argv, desc), vm);

            if (vm.count("help") or argc == 1) {
                std::cout << desc << std::endl;
                return 0; //SUCCESS
            }
            po::notify(vm);
        }
        catch (const po::error &e) {
            std::cerr << "ERROR: " << e.what() << "\n" << std::endl;
            std::cerr << desc << std::endl;
            return 1; //FAILURE
        }
        catch (const std::exception &e) {
            std::cerr << "ERROR: " << e.what() << "\n" << std::endl;
            return 1;
        }

        // check extension to determine if this is a SAM file
        std::string file_type = "";
        std::string infile = BA::to_lower_copy(in);
        if (BA::ends_with(infile, ".sam") or
            BA::ends_with(infile, ".sam.gz")) {
            file_type = "SAM";
        } else {
            std::cerr << "Unable to determine file type of "
            << in << std::endl
            << "Recognized extensions are .sam, and .sam.gz"
            << " (capitalization not important)."
            << std::endl;
            return 1;
        }


        // don't allow negative parameters (won't get checked by base method, since its
        // params are unsigned
        if (min_mapq < 0) {
            throw std::invalid_argument("ERROR: min_mapq must be positive.");
        }

        std::cout << "Attempting to parse " << file_type << " file "
        << in
        << " and write to "
        << out;
        std::cout << std::endl;
        if (debug_out.size() > 0) {
            std::cout << "\n\twriting debug intermediate info to " << debug_out << std::endl;
        }
        std::cout << "\n\tusing min_mapq=" << min_mapq << "." << std::endl;

        std::cout << "\ttreating input reads as ";
        if (input_is_unpaired) {std::cout << "unpaired reads\n"; }
        else { std::cout << "merged and/or paired reads\n"; }

        if (not input_is_unpaired) { std::cout << "\ttreating paired reads mapping to a max fragment size of "
         << max_paired_fragment_length << " as a single read\n"; }

        if (require_forward_primer_mapped) {
            std::cout << "\trequiring read mapping to expected forward primer location within "
            << max_primer_offset << " nucleotides\n";
        }

        if (require_reverse_primer_mapped) {
            std::cout << "\trequiring read mapping to expected reverse primer location within "
            << max_primer_offset << " nucleotides\n";
        }

        if (trim_primers) {
            std::cout << "\ttrimming amplicon primers provided in " << primers << "\n";
        } else {
            std::cout << "\ttrimming " << exclude_3prime << "from right end of reads (to account for random primer)\n";
        }

        // FIXME: describe more fully and indicate interactions with right_align/other params
        std::cout << "\tsequence variant mode is ";
        if (variant_mode) { std::cout << "on\n"; }
        else { std::cout << "off\n"; }

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

        if (exclude_3prime > 0){
            std::cout << "\texcluding mutations within " << exclude_3prime << " nucleotides of read 3-prime end\n";
        }
        std::cout << "\tmerging adjacent mutations within " << max_internal_match << " nucleotides of each other\n";
        std::cout << "\texcluding mutations with any basecall q-scores below " << min_qual << "\n";

        if (use_only_mutation_type.length() != 0) {
            std::cout << "\tusing only mutations of the type: " << use_only_mutation_type << '\n';
        }

        std::cout << std::flush;

        if (file_type=="SAM"){
            /*const std::string &filename,
                  const std::string &outname,
                  const std::string &primers_filename,
                  const bool paired,
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
                  const bool debug,
                  const bool warn_on_no_mapped = false*/
            mutation_parser::parseSAM(in,
                                      out,
                                      debug_out,
                                      primers,
                                      max_paired_fragment_length,
                                      min_mapq,
                                      right_align_ambig_dels,
                                      right_align_ambig_ins,
                                      max_internal_match,
                                      min_qual,
                                      exclude_3prime,
                                      use_only_mutation_type,
                                      variant_mode,
                                      trim_primers,
                                      require_forward_primer_mapped,
                                      require_reverse_primer_mapped,
                                      max_primer_offset,
                                      input_is_unpaired,
                                      debug,
                                      warn_on_no_mapped);
        }

        std::cout << "... Successfully parsed mutations from file." << std::endl;
    }
    catch (const BF::filesystem_error &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return 1; //FAILURE
    }
    catch (...) {
        std::cerr << "Unknown error." << std::endl;
        return 1;
    }
    return 0; //SUCCESS
}

