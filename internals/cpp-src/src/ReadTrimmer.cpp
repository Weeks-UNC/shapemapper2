/** @file
 * @brief Trim FASTQ reads by windowed mean phred score. Primary interface functions.
 */

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2018 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#include "ReadTrimmer.h"


// This is the simplest way to get doxygen docs to show actual values in function
// call signature.
#define DEFAULT_WINDOW_SIZE 5
#define DEFAULT_MIN_PHRED 20
#define DEFAULT_MIN_LENGTH 25

namespace read_trimmer {

    /**
     * @brief Base method for read trimming by windowed phred score.
     *
     * @param read          Basecalls
     * @param phred_scores  Basecall quality scores (in PHRED+33 ASCII encoding)
     * @param [min_phred]   Minimum windowed average phred score to keep
     * @param [min_length]  Minimum read length after trimming to keep
     * @param [window_size] Window size in basecalls
     *
     * @return Returns a tuple containing the trimmed read and trimmed phred scores.
     *         If trimmed read is shorter than min_length, <"N","!"> will be returned.
     */
    boost::tuple<std::string, std::string>
    trimRead(const std::string &read,
             const std::string &phred_scores,
             unsigned int window_size = DEFAULT_WINDOW_SIZE,
             unsigned int min_phred = DEFAULT_MIN_PHRED,
             unsigned int min_length = DEFAULT_MIN_LENGTH) {

        // Check sensible params
        if (min_length < window_size) {
            throw std::invalid_argument("ERROR: Read trimming min_length cannot be less than window_size.");
        } else if (read.length() != phred_scores.length()) {
            throw std::invalid_argument("ERROR: Read length does not match phred scores length.");
        }
        // if read is shorter than window, go ahead and return
        if (read.length() < window_size or read.length() == 0) {
            return boost::make_tuple("N", "!");
        }

        int low_qual_index = detail::locateLowQualityWindow(phred_scores, window_size, min_phred);
        if (low_qual_index != -1) {
            // we have identified a low-quality window
            std::string trimmed_read = read.substr(0, low_qual_index);
            std::string trimmed_phred = phred_scores.substr(0, low_qual_index);
            if (trimmed_read.length() < min_length) {
                return boost::make_tuple("N", "!");
            } else {
                return boost::make_tuple(trimmed_read, trimmed_phred);
            }

        } else {
            // otherwise return untrimmed copy
            std::string r = read;
            std::string p = phred_scores;
            return boost::make_tuple(r, p);
        }

    }


    /**
     * @brief Open a FASTQ file, trim reads, and write trimmed reads to
     * new file. Input and/or output may be gzip compressed, specified
     * with a ".gz" file extension.
     */
    void
    trimFastq(std::string filename,
              std::string outname,
              unsigned int window_size = DEFAULT_WINDOW_SIZE,
              unsigned int min_phred = DEFAULT_MIN_PHRED,
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

        BI::filtering_istream in;
        // universal newline support filter
        in.push(BI::newline_filter(BI::newline::posix));
        if (BF::extension(BF::path(filename)) == ".gz") {
            // decompress gzip if file looks compressed
            in.push(BI::gzip_decompressor());
        }
        in.push(file_in);

        // create path to output file if needed
        BF::path outpath(outname);
        if (outpath.has_parent_path() and
            not (BF::exists(outpath.parent_path()))){
                BF::create_directories(outpath.parent_path());
        }

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
            
                try {
                    boost::tie(block[1], block[3]) = read_trimmer::trimRead(block[1],
                                                                            block[3],
                                                                            window_size,
                                                                            min_phred,
                                                                            min_length);
                } catch (const std::invalid_argument &e) {
                    std::cout << "Error at line " << c*4+1 << " in input file " << filename << ":" << std::endl;  
                    throw e;                    
                    //std::cout << e.what() << std::endl;
                }
                
                // write trimmed reads to output file
                for (int n = 0; n < 4; n++) {
                    out << block[n] << '\n';
                }
                c++;
                k = 0;
            }

        }
        out << std::flush;
        if (c < 1) {
            throw std::runtime_error("ERROR: Input file " + filename + " contains no reads.");
        }

        //throw std::runtime_error("ERROR: intentionally triggered exception for testing");
    }


    /* TODO: move gzip functionality to stream wrappers in detail namespace?
     */

}



































