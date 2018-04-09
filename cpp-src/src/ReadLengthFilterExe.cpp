/** @file
 * @brief Filter FASTQ reads to remove reads below some length. Commandline executable.
 */

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2017 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#include <boost/program_options.hpp>

#include "ReadLengthFilter.cpp"


namespace po = boost::program_options;

int main(int argc, char *argv[]) {
    try {

        po::options_description desc("Usage");
        desc.add_options()
                ("help,h", "print usage message")

                ("in,i", po::value<std::string>()->required(), "FASTQ input file path")

                ("out,o", po::value<std::string>()->required(), "filtered FASTQ output file path")

                ("min_length,l",
                 po::value<int>()->default_value(DEFAULT_MIN_LENGTH, std::to_string(DEFAULT_MIN_LENGTH)),
                 "minimum read length to allow");

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
        std::cout << "Attempting to filter fastq file "
        << vm["in"].as<std::string>()
        << " and write to "
        << vm["out"].as<std::string>();
        std::cout << "\n... using min_length=" << vm["min_length"].as<int>() << "." << std::endl;

        // don't allow negative parameters (won't get checked by base method, since its
        // params are unsigned
        if (vm["min_length"].as<int>() < 0) {
            throw std::invalid_argument("ERROR: min_length must be positive.");
        }

        read_filter::filterFastq(vm["in"].as<std::string>(),
                                 vm["out"].as<std::string>(),
                                 vm["min_length"].as<int>());

        std::cout << "... Successfully filtered fastq file." << std::endl;
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