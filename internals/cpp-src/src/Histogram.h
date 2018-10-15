/** @file
 * @brief Utility class for simple histogram calculation.
 */

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2018 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#ifndef SHAPEMAPPER_HISTOGRAM_H
#define SHAPEMAPPER_HISTOGRAM_H

#include <math.h>

/**
 * @brief Simple histogram table for storing approximate distribution of read
 *        lengths or mutations per read.
 */
class Histogram {
public:
    std::string title;
    int total_reads = 0;
    std::vector<int> bin_lefts;
    std::vector<int> counts;
    int start;
    int end;
    std::string scale;

    /**
     * @brief Set up N bins from start to end, inclusive, including a bin for
     *        values greater than end. (N = total bins)
     *
     * @param title
     * @param start
     * @param end
     * @param total_bins
     * @param scale
     */
    Histogram(const std::string title,
              const int start,
              const int end,
              const int total_bins,
              const std::string scale = "linear") {
        // TODO: implement log scale for read length distribution
        this->title = title;
        this->start = start;
        this->end = end;
        this->scale = scale;
        if (scale == "linear") {
            for (int i = 0; i < total_bins; i++) {
                int left = start + floor(i*(float(end - start) / float(total_bins - 1)));
                bin_lefts.push_back(left);
                counts.push_back(0);
            }
        } else if (scale == "log") {

        }
    }

    /**
     * @brief Find bin and update counts
     */
    void count(const int value) {
        total_reads++;
        if (scale == "linear") {
            int total_bins = bin_lefts.size();
            int i = std::min(total_bins - 1, int(floor(float(value - start) /
                                                      (float(end - start) / float(total_bins - 1)))));
            counts.at(i)++;
        } else if (scale == "log") {

        }
    }

    std::string printCountsRow() {
        std::string o;
        for (auto &count:counts) {
            o += std::to_string(count);
            // add tabs after all elements except the last one
            if (&count != &counts.back()) o += "\t";
        }
        return o;
    }

    std::string printFreqTable(std::string bin_format="simple") {
        std::string o;
        o += title + "\n";
        o += "--------------------\n";
        std::string bin_header = "bin_left";
        if (bin_format=="range") {
            bin_header = "bin_range";
        }
        o += bin_header + "\tfrequency\n";
        for (int i = 0; i < counts.size(); i++) {
            std::string bin_label = std::to_string(bin_lefts.at(i));
            if (bin_format=="range"){
                if (i < counts.size()-1) {
                    bin_label = "[" + std::to_string(bin_lefts.at(i)) + "," + std::to_string(bin_lefts.at(i + 1)-1) + "]";
                } else {
                    bin_label = ">="+std::to_string(bin_lefts.at(i));
                }
            }
            o += bin_label + "\t";
            o += std::to_string(float(counts.at(i))/float(total_reads))+"\n";
        }
        o += "--------------------\n";
        return o;
    }
};

#endif //SHAPEMAPPER_HISTOGRAM_H