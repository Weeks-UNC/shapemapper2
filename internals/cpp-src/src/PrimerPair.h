/** @file
 *  @brief PrimerPair class for amplicon primer pair read filtering
 */

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2018 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#ifndef SHAPEMAPPER_PRIMERPAIR_H
#define SHAPEMAPPER_PRIMERPAIR_H

#include <iostream>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/newline.hpp>

namespace BF = boost::filesystem;
namespace BI = boost::iostreams;


class PrimerPair {
public:
    int fw_left, fw_right, rv_left, rv_right;

    PrimerPair() {
        fw_left = -999;
        fw_right = -999;
        rv_left = -999;
        rv_right = -999;
    };

    PrimerPair(const PrimerPair &p) {
        fw_left = p.fw_left;
        fw_right = p.fw_right;
        rv_left = p.rv_left;
        rv_right = p.rv_right;
    }

    PrimerPair(const std::string line) {
        std::vector<std::string> fields;
        std::string trimmed = boost::trim_copy(line);
        boost::split(fields, trimmed, boost::is_space(), boost::token_compress_on);

        if (fields.size() < 4) {
            throw std::runtime_error("Error: unable to parse incomplete line in primer file.");
        }
        try {
            fw_left = stoi(fields[0]);
            fw_right = stoi(fields[1]);
            rv_left = stoi(fields[2]);
            rv_right = stoi(fields[3]);
        } catch (std::exception &err) {
            throw std::runtime_error("Error: line is incorrectly formatted (couldn't parse primer locations).");
        }
    }

    std::string toString() const {
        using std::to_string;
        std::string s = "";
        s += "fw_left: " + to_string(fw_left) + "\n";
        s += "fw_right: " + to_string(fw_right) + "\n";
        s += "rv_left: " + to_string(rv_left) + "\n";
        s += "rv_right: " + to_string(rv_right) + "\n";
        return s;
    }
};

std::vector<PrimerPair>
loadPrimerPairs(std::string filename) {
    std::ifstream file_in(filename, std::ios_base::in | std::ios_base::binary);
    if (!file_in) {
        // Do additional checks to see if the file exists or if it's a permissions issue
        if (!(BF::is_regular_file(filename))) {
            throw std::runtime_error("ERROR: Input file " + filename + " not found.");
        }
        throw std::runtime_error("ERROR: Could not open input file " + filename +
                                 " - unknown error.\nCheck file and folder permissions.");
    }
    BI::filtering_istream in;
    // universal newline support filter
    in.push(BI::newline_filter(BI::newline::posix));
    in.push(file_in);

    std::vector<PrimerPair> primer_pairs;
    std::string line;
    while (std::getline(in, line)) {
        std::string trimmed = boost::trim_copy(line);
        // skip blank lines, RNA name, and primer sequences
        if (trimmed.length()<1 or
            trimmed[0]=='>' or
            std::isalpha(trimmed[0])){
            continue;
        }
        primer_pairs.push_back(PrimerPair(trimmed));
    }
    return primer_pairs;
}

#endif //SHAPEMAPPER_PRIMERPAIR_H