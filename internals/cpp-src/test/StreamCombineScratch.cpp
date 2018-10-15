#include <memory>
#include <map>
#include <deque>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>

#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/newline.hpp>
#include <boost/iostreams/device/file.hpp>

namespace BF = boost::filesystem;
namespace BI = boost::iostreams;

std::string FILEPATH = __FILE__;

BF::path getTestFileDir() {
    BF::path PARENTPATH = BF::path(FILEPATH).parent_path();
    BF::path FILEDIR = PARENTPATH / "files";
    return FILEDIR;
}

namespace BF = boost::filesystem;
namespace BI = boost::iostreams;

int main(int argc, char *argv[]) {
    std::vector<std::string> filenames;
    filenames.push_back(getTestFileDir().string()+"/unpaired_mutations.mut");
    filenames.push_back(getTestFileDir().string()+"/paired_mutations.mut");

    // this garbage is needed because streams don't support move
    std::vector<std::unique_ptr<BI::filtering_istream>> files;

    for (auto & filename : filenames) {
        // again, this junk is required because streams don't support move
        files.emplace_back(
                std::unique_ptr<BI::filtering_istream>(new BI::filtering_istream()));

        files.back()->push(BI::newline_filter(BI::newline::posix));
        //std::ifstream file_in(filenames[0], std::ios_base::in | std::ios_base::binary);
        BI::file_source fs(filename);
        if (not fs.is_open()) {
            throw std::runtime_error("ERROR: Could not open input file " + filename);
        }
        files.back()->push(fs);

    }

    std::vector<std::string> input_lines(files.size());
    std::vector<bool> eof(files.size());
    std::vector<std::string> non_empty_lines;
    while(true) {
        for (int i=0; i<files.size(); i++) {
            eof[i] = !std::getline(*(files[i]), input_lines[i]);
        }

        // terminate if we reach the end of all files
        if (std::all_of(eof.begin(), eof.end(), [](bool v){return v;})) {
            break;
        }
        for (int i=0; i<files.size(); i++) {
            if (not eof[i]) { non_empty_lines.push_back(input_lines[i]); }
        }

        for (auto & line : non_empty_lines) {
            // process line
            std::cout << line << std::endl;
        }
        non_empty_lines.resize(0);
    }

    return 0;
}
