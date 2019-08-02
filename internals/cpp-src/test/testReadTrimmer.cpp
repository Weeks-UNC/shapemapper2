/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2018 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#include "gtest/gtest.h"
#include "ReadTrimmer.cpp"

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
        filedir = BF::path(BASEPATH) / "internals" / "cpp-src" / "test" / "files";
        //std::cout << "Using BASEPATH " << BASEPATH << " passed by argument" << std::endl;
    }
    BF::create_directory(filedir / "tmp");
    return filedir;
}

std::string getTestFilePath() {
    return (getTestFileDir() / "3_R1.fastq").string();
}

std::string getTestGzipFilePath() {
    return (getTestFileDir() / "3_R1.fastq.gz").string();
}

std::string getTestCRFilePath() {
    return (getTestFileDir() / "CR_line_ending.fastq").string();
}


TEST(ReadTrimmerTest, charToPhred) {
    char phred_char = '@';
    int result = read_trimmer::detail::charToPhred(phred_char);
    EXPECT_EQ(31, result);
}

TEST(ReadTrimmerTest, ErrorOnCharToPhredBounds) {
    char phred_char = 32; // Should produce a phred score of -1
    EXPECT_THROW(read_trimmer::detail::charToPhred(phred_char), std::invalid_argument);
}

TEST(ReadTrimmerTest, LeavesHighQualityRead) {
    const std::string read("ATGCATGCATGCATGCATGC");
    const std::string phred("~~~~~~~~~~~~~~~~~~~~");
    unsigned int window_size = 1;
    unsigned int min_phred = 30;
    unsigned int min_length = 5;
    std::string trimmed_read;
    std::string trimmed_phred;
    boost::tie(trimmed_read, trimmed_phred) = read_trimmer::trimRead(read, phred, window_size, min_phred, min_length);
    EXPECT_EQ(read, trimmed_read);
    EXPECT_EQ(phred, trimmed_phred);
}

TEST(ReadTrimmerTest, EliminatesLowQualityRead) {
    const std::string read("ATGCATGCATGCATGCATGC");
    const std::string phred("!!!!!!!!!!!!!!!!!!!!");
    unsigned int window_size = 1;
    unsigned int min_phred = 30;
    unsigned int min_length = 5;
    std::string trimmed_read;
    std::string trimmed_phred;
    boost::tie(trimmed_read, trimmed_phred) = read_trimmer::trimRead(read, phred, window_size, min_phred, min_length);
    EXPECT_EQ("N", trimmed_read);
    EXPECT_EQ("!", trimmed_phred);
}

TEST(ReadTrimmerTest, EliminatesShortRead) {
    const std::string read("ATGCATGCATGCATGCATGC");
    const std::string phred("~~~~!!!!!!!!!!!!!!!!");
    unsigned int window_size = 1;
    unsigned int min_phred = 30;
    unsigned int min_length = 5;
    std::string trimmed_read;
    std::string trimmed_phred;
    boost::tie(trimmed_read, trimmed_phred) = read_trimmer::trimRead(read, phred, window_size, min_phred, min_length);
    EXPECT_EQ("N", trimmed_read);
    EXPECT_EQ("!", trimmed_phred);
}

TEST(ReadTrimmerTest, TrimsModerateQualityRead) {
    const std::string read("ATGCATGCATGCATGCATGC");
    const std::string phred("~~~~~~~~~~~~!!!!!!!!");
    unsigned int window_size = 1;
    unsigned int min_phred = 30;
    unsigned int min_length = 5;
    std::string trimmed_read;
    std::string trimmed_phred;
    boost::tie(trimmed_read, trimmed_phred) = read_trimmer::trimRead(read, phred, window_size, min_phred, min_length);
    EXPECT_EQ("ATGCATGCATGC", trimmed_read);
    EXPECT_EQ("~~~~~~~~~~~~", trimmed_phred);
}

TEST(ReadTrimmerTest, HandlesWindow) {
    // + 10
    // ! 0
    // Should cut off first window that includes a '!' given a min_phred of 10
    const std::string read("ATGCATGCATGCATGCATGC");
    const std::string phred("++++++++++++!!!!!!!!");
    unsigned int window_size = 2;
    unsigned int min_phred = 10;
    unsigned int min_length = 5;
    std::string trimmed_read;
    std::string trimmed_phred;
    boost::tie(trimmed_read, trimmed_phred) = read_trimmer::trimRead(read, phred, window_size, min_phred, min_length);
    EXPECT_EQ("ATGCATGCATG", trimmed_read);
    EXPECT_EQ("+++++++++++", trimmed_phred);
}

TEST(ReadTrimmerTest, HandlesEmptyString) {
    const std::string read("");
    const std::string phred("");
    unsigned int window_size = 1;
    unsigned int min_phred = 30;
    unsigned int min_length = 5;
    std::string trimmed_read;
    std::string trimmed_phred;
    boost::tie(trimmed_read, trimmed_phred) = read_trimmer::trimRead(read, phred, window_size, min_phred, min_length);
    EXPECT_EQ("N", trimmed_read);
    EXPECT_EQ("!", trimmed_phred);
}

TEST(ReadTrimmerTest, ErrorOnWindowSizeLargerThanMinLength) {
    const std::string read("ATGC");
    const std::string phred("~~~~");
    unsigned int window_size = 10;
    unsigned int min_phred = 30;
    unsigned int min_length = 5;
    std::string trimmed_read;
    std::string trimmed_phred;
    EXPECT_THROW(read_trimmer::trimRead(read, phred, window_size, min_phred, min_length), std::invalid_argument);
}

TEST(ReadTrimmerTest, ErrorOnNewlineChar) {
    const std::string read("ATGCATGC\n");
    const std::string phred("~~~~~~~~\n");
    unsigned int window_size = 10;
    unsigned int min_phred = 30;
    unsigned int min_length = 5;
    std::string trimmed_read;
    std::string trimmed_phred;
    EXPECT_THROW(read_trimmer::trimRead(read, phred, window_size, min_phred, min_length), std::invalid_argument);
}

TEST(ReadTrimmerTest, ErrorOnMismatchedLengths) {
    const std::string read("ATGCATGC");
    const std::string phred("~~~~~~~~~");
    unsigned int window_size = 10;
    unsigned int min_phred = 30;
    unsigned int min_length = 5;
    std::string trimmed_read;
    std::string trimmed_phred;
    EXPECT_THROW(read_trimmer::trimRead(read, phred, window_size, min_phred, min_length), std::invalid_argument);
}


TEST(FileHandling, FastqFileTrim) {
    std::string file_in = getTestFilePath();
    std::string file_out = (getTestFileDir() / "tmp" / "trimmed.fastq").string();
    EXPECT_NO_THROW(read_trimmer::trimFastq(file_in, file_out));
    EXPECT_EQ(read_trimmer::detail::countLines(file_in),
              read_trimmer::detail::countLines(file_out));
}

TEST(FileHandling, HandlesLineEndingCR) {
    std::string file_in = getTestCRFilePath();
    std::string file_out = (getTestFileDir() / "tmp" / "trimmed_from_CR.fastq").string();
    EXPECT_NO_THROW(read_trimmer::trimFastq(file_in, file_out));
    EXPECT_EQ(read_trimmer::detail::countLines(getTestFilePath()),
              read_trimmer::detail::countLines(file_out));
}

TEST(FileHandling, GzipInputFastqFileTrim) {
    std::string file_in = getTestGzipFilePath();
    std::string file_out = (getTestFileDir() / "tmp" / "trimmed_from_gz.fastq").string();
    EXPECT_NO_THROW(read_trimmer::trimFastq(file_in, file_out));
    EXPECT_EQ(read_trimmer::detail::countLines(getTestFilePath()),
              read_trimmer::detail::countLines(file_out));
}

TEST(FileHandling, GzipOutputFastqFileTrim) {
    std::string file_in = getTestFilePath();
    std::string file_out = (getTestFileDir() / "tmp" / "trimmed.fastq.gz").string();
    EXPECT_NO_THROW(read_trimmer::trimFastq(file_in, file_out));
    // probably need additional check here
}

// TODO: check for correct exception using gmock's StartsWith or HasSubstr matchers
TEST(FileHandling, ExceptionOnInputNotFound) {
    std::string file_in = (getTestFileDir() / "does_not_exist.fastq").string();
    std::string file_out = (getTestFileDir() / "tmp" / "trash.fastq").string();
    EXPECT_THROW(read_trimmer::trimFastq(file_in, file_out), std::runtime_error);
}

TEST(FileHandling, ExceptionOnMalformedFastq) {
    std::string file_in = (getTestFileDir() / "malformed.fastq").string();
    std::string file_out = (getTestFileDir() / "tmp" / "trash.fastq").string();
    EXPECT_THROW(read_trimmer::trimFastq(file_in, file_out), std::runtime_error);
}

TEST(FileHandling, ExceptionOnNearlyEmptyFastq) {
    std::string file_in = (getTestFileDir() / "nearly_empty.fastq").string();
    std::string file_out = (getTestFileDir() / "tmp" / "trash.fastq").string();
    EXPECT_THROW(read_trimmer::trimFastq(file_in, file_out), std::runtime_error);
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  // set location of top-level shapemapper path so we can locate data files
  if (argc > 1){
      BASEPATH = argv[1];
  }
  return RUN_ALL_TESTS();
}
