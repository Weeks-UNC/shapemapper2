/** @file
 * @brief Parse mapped BAM alignments into mutations and
 *        reconstructed alignment target sequences. Commandline executable.
 */

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2017 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "MutationParser.cpp"


namespace po = boost::program_options;
namespace BA = boost::algorithm;

int main(int argc, char *argv[]) {
    try {

        po::options_description desc("Usage");
        desc.add_options()
                ("help,h", "print usage message")

                ("in,i", po::value<std::string>()->required(), "SAM or BAM input file path")

                ("out,o", po::value<std::string>()->required(), "parsed mutations output file path")

                ("min_mapq,m",
                 po::value<int>()->default_value(DEFAULT_MIN_MAPQ, std::to_string(DEFAULT_MIN_MAPQ)),
                 "minimum reported mapping quality to allow")

                 ("warn_on_no_mapped,w",
                  po::bool_switch()->default_value(false),
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

        // check extension to determine whether this is a SAM or BAM file
        std::string file_type = "";
        std::string infile = vm["in"].as<std::string>();
        BA::to_lower(infile);
        if (BA::ends_with(infile, ".sam") or
            BA::ends_with(infile, ".sam.gz")) {
            file_type = "SAM";
        } else if (BA::ends_with(infile, ".bam")) {
            file_type = "BAM";
        } else {
            std::cerr << "Unable to determine file type of "
            << vm["in"].as<std::string>() << std::endl
            << "Recognized extensions are .sam, .bam, and .sam.gz"
            << " (capitalization not important)."
            << std::endl;
            return 1;
        }

        std::cout << "Attempting to parse " << file_type << " file "
        << vm["in"].as<std::string>()
        << " and write to "
        << vm["out"].as<std::string>();
        std::cout << std::endl;
        std::cout << "\n... using min_mapq=" << vm["min_mapq"].as<int>() << "." << std::endl;

        // don't allow negative parameters (won't get checked by base method, since its
        // params are unsigned
        if (vm["min_mapq"].as<int>() < 0) {
            throw std::invalid_argument("ERROR: min_mapq must be positive.");
        }


        if (file_type=="BAM") {
            mutation_parser::parseBAM(vm["in"].as<std::string>(),
                                      vm["out"].as<std::string>(),
                                      vm["min_mapq"].as<int>());
        } else if (file_type=="SAM"){
            mutation_parser::parseSAM(vm["in"].as<std::string>(),
                                      vm["out"].as<std::string>(),
                                      vm["min_mapq"].as<int>(),
                                      vm["warn_on_no_mapped"].as<bool>());
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

