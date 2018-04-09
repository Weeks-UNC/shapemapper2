/** @file
 * @brief Filter FASTQ reads to remove reads below some length. Primary interface functions.
 */

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2017 Steven Busan.                     *
 *----------------------------------------------------------------------*/

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

#define DEFAULT_MIN_LENGTH 2

namespace BF = boost::filesystem;
namespace BI = boost::iostreams;

namespace read_filter {

    /**
     * @brief Open a FASTQ file, filter reads by length, and write passing reads to
     * new file. Input and/or output may be gzip compressed, specified
     * with a ".gz" file extension.
     */
    void
    filterFastq(std::string filename,
              std::string outname,
              unsigned int min_length = DEFAULT_MIN_LENGTH) {

        // ifstream constructor doesn't seem to throw exceptions, just returns null

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
            else {
                // Can't figure out how to check read permissions using boost,
                // Not worth any more time.
            }
            throw std::runtime_error("ERROR: Could not open input file " + filename +
                                     " - unknown error.\nCheck file and folder permissions.");
        }

        // TODO: there may be some performance hit to use iostreams. Consider setting std::ios_base::sync_with_stdio, other stuff
        BI::filtering_istream in;
        // universal newline support filter
        in.push(BI::newline_filter(BI::newline::posix));
        if (BF::extension(BF::path(filename)) == ".gz") {
            // decompress gzip if file looks compressed
            in.push(BI::gzip_decompressor());
        }
        in.push(file_in);

        //BF::ofstream file_out(outname);
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

        std::string line;
        std::string block[4];

        unsigned int k = 0;
        size_t c = 0;
        while (std::getline(in, line)) {
            block[k] = line;
            k++;
            if (k > 3) {
                if (block[0].substr(0, 1) != "@" or block[2] != "+") {
                    throw std::runtime_error("ERROR: Input file " + filename + " does not appear FASTQ formatted.");
                }

                if (block[1].length() >= min_length) {
                    // write reads meeting length threshold to output file
                    for (int n = 0; n < 4; n++) {
                        out << block[n] << '\n';
                    }
                }
                c++;
                k = 0;
            }

        }
        out << std::flush;
        if (c < 1) {
            throw std::runtime_error("ERROR: Input file " + filename + " contains no reads.");
        }
    }

}



































