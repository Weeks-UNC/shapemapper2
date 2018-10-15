/** @file
 * @brief Trim FASTQ reads by windowed mean phred score. Commandline executable.
 */

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2018 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#include <boost/program_options.hpp>

#include "ReadTrimmer.cpp"


namespace po = boost::program_options;

int main(int argc, char *argv[]) {
    try {
        std::string in;
        std::string out;
        int window_size;
        int min_phred;
        int min_length;

        po::options_description desc("Usage");
        desc.add_options()
                ("help,h", "print usage message")

                ("in,i", po::value<std::string>(&in)->required(), "FASTQ input file path")

                ("out,o", po::value<std::string>(&out)->required(), "trimmed FASTQ output file path")

                ("window_size,w",
                 po::value<int>(&window_size)->default_value(DEFAULT_WINDOW_SIZE, std::to_string(DEFAULT_WINDOW_SIZE)),
                 "size of window to average quality scores")

                ("min_phred,p",
                 po::value<int>(&min_phred)->default_value(DEFAULT_MIN_PHRED, std::to_string(DEFAULT_MIN_PHRED)),
                 "minimum average quality score to allow")

                ("min_length,l",
                 po::value<int>(&min_length)->default_value(DEFAULT_MIN_LENGTH, std::to_string(DEFAULT_MIN_LENGTH)),
                 "minimum trimmed read length to allow");

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
        std::cout << "Attempting to trim fastq file "
        << in
        << " and write to "
        << out;
        std::cout << "\n... Using params: "
        << "window_size=" << window_size
        << ", min_phred=" << min_phred
        << ", min_length=" << min_length << "." << std::endl;

        // don't allow negative parameters (won't get checked by base method, since its
        // params are unsigned
        if (window_size < 1) {
            throw std::invalid_argument("ERROR: window_size must be at least 1.");
        }
        if (min_phred < 0) {
            throw std::invalid_argument("ERROR: min_phred must be positive.");
        }
        if (min_length < 0) {
            throw std::invalid_argument("ERROR: min_length must be positive.");
        }

        read_trimmer::trimFastq(in,
                                out,
                                window_size,
                                min_phred,
                                min_length);


        std::cout << "... Successfully trimmed fastq file." << std::endl;
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