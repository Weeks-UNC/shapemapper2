/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2017 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#include "gtest/gtest.h"
#include "ReadTrimmer.cpp" // has a useful line counting function, will eventually move to util namespace
#include "ReadLengthFilter.cpp"

namespace BF = boost::filesystem;

std::string FILEPATH = __FILE__;
std::string BASEPATH = "";

BF::path getTestFileDir() {
    BF::path filedir;
    if (BASEPATH == "") {
        // data file location if explicit path to main shapemapper directory
        // not provided by test runner (assumes software location has not
        // changed since compilation)
        filedir = BF::path(FILEPATH).parent_path() / "files";
    } else {
        filedir = BF::path(BASEPATH) / "cpp-src" / "test" / "files";
    }
    return filedir;
}

std::string getTestFilePath() {
    return (getTestFileDir() / "contains_short_read.fastq").string();
}

TEST(FileHandling, FastqFileFilter) {
    std::string file_in = getTestFilePath();
    std::string file_out = (getTestFileDir() / "tmp" / "filtered.fastq").string();
    EXPECT_NO_THROW(read_filter::filterFastq(file_in, file_out));
    EXPECT_EQ(read_trimmer::detail::countLines(file_in),
              read_trimmer::detail::countLines(file_out)+4);
}


int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc >= 2) {
        BASEPATH = argv[1];
    }
    return RUN_ALL_TESTS();
}