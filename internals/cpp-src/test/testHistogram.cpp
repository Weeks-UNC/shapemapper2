/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2018 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#include "gtest/gtest.h"
#include "Histogram.h"


TEST(Histogram, linearSimple){
Histogram mut_per_read("Mutations per read",1,5,5);
mut_per_read.count(2);
mut_per_read.count(2);
mut_per_read.count(4);
mut_per_read.count(5);
mut_per_read.count(7);
EXPECT_EQ(mut_per_read.total_reads, 5);
EXPECT_EQ(mut_per_read.printCountsRow(), "0\t2\t0\t1\t2");
}

TEST(Histogram, linearSimpleTable){
Histogram mut_per_read("Mutations per read",0,2,3);
mut_per_read.count(1);
mut_per_read.count(1);
mut_per_read.count(1);
mut_per_read.count(4);
std::string exp = "Mutations per read\n"
        "--------------------\n"
        "bin_left\tfrequency\n"
        "0\t0.000000\n"
        "1\t0.750000\n"
        "2\t0.250000\n"
        "--------------------\n";
EXPECT_EQ(mut_per_read.total_reads, 4);
EXPECT_EQ(mut_per_read.printFreqTable(), exp);
}

TEST(Histogram, linearBiggerBins) {
Histogram mut_per_read("Mutations per read",1, 10, 5);
mut_per_read.count(4);
mut_per_read.count(4);
mut_per_read.count(8);
mut_per_read.count(10);
mut_per_read.count(14);
EXPECT_EQ(mut_per_read.total_reads, 5);
EXPECT_EQ(mut_per_read.printCountsRow(), "0\t2\t0\t1\t2");
}


