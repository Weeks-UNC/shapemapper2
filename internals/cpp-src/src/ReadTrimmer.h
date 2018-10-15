/** @file
 * @brief Trim FASTQ reads by windowed mean phred score. Utility functions.
 */

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2018 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#ifndef SHAPEMAPPER_READTRIMMER_H
#define SHAPEMAPPER_READTRIMMER_H

#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>

#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/newline.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace BF = boost::filesystem;
namespace BI = boost::iostreams;

// "Utility" function implementations
namespace read_trimmer { namespace detail {

    int charToPhred(char c) {
        int ascii_value = (int) c;
        int phred_score = ascii_value - 33;
        // Do some sanity checks, and throw exception if needed.
        if (phred_score < 0 or phred_score > 93) {
            throw std::invalid_argument(
                    "Phred score of " + std::to_string(phred_score) + " for char " + std::to_string(c) +
                    " is out of expected range 0-93. Input quality score lines may have unprintable characters.");
        }
        return phred_score;
    }

    /**
     * Return the leftmost index of the first window (scanning left-to-right) with mean
     * phred score below min_phred.
     * Returns -1 if all positions pass
     */
    int locateLowQualityWindow(const std::string &phred_scores,
                               unsigned int window_size,
                               unsigned int min_phred) {
        std::vector<unsigned int> scores(phred_scores.length());
        for (int i = 0; i < phred_scores.length(); ++i) {
            try {
                scores[i] = charToPhred(phred_scores[i]);
            } catch (std::invalid_argument) {
                throw std::invalid_argument(
                        "ERROR: Phred score string contains whitespace or non-printable characters. Check line endings.");
            }
        }
        int low_qual_index = -1;
        for (int i = 0; i < scores.size() - window_size; i++) {
            int sum = 0;
            for (int j = 0; j < window_size; j++) {
                sum += scores[i + j];
            }
            float mean = (float) sum / (float) window_size;
            if (mean < min_phred) {
                low_qual_index = i;
                break;
            }
        }
        return low_qual_index;

    }


    /**
     * Count lines in file. (probably only works correctly for platform-specific line-endings.)
     */
    // TODO: move to a utility namespace outside read_trimmer
    // Used in tests
    int
    countLines(std::string filename) {
        std::ifstream f(filename);
        int c = 0;
        std::string line;
        while (std::getline(f, line)) {
            c++;
        }
        return c;
    }

}}

#endif //SHAPEMAPPER_READTRIMMER_H
