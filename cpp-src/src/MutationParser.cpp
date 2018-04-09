/** @file
 * @brief Main interface functions for parsing BAM file or single alignments
 *        into mutations and reconstructed target sequences.
 */

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2017 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#include "MutationParser.h"

#define DEFAULT_MIN_MAPQ 30

namespace mutation_parser {
    using namespace mutation;

    /**
     * @brief Given a single BAM alignment entry (assumed mapped) -
     *        parse mutations, reconstruct target sequence over aligned region,
     *        and merge ambiguously aligned indels with alternate placements.
     *
     * @param al Single alignment
     * @return Returns the information needed to count RT mutations and depths or
     * count variants: Left-most position in alignment target sequence (0-based),
     * right-most position (0-based), reconstructed target sequence over this range,
     * and mutations present in this read.
     */
    boost::tuple<int, int, std::string, std::string, std::vector<Mutation>>
    parseMutations(const BAM::BamAlignment &al) {
        const std::string tag = "MD";
        std::string md_tag_contents;
        if (!al.GetTag(tag, md_tag_contents)) {
            throw std::runtime_error("Error: no MD tag in alignment.");
        }

        return detail::parseMutations(al.Position,
                                      al.GetEndPosition()-1, // fixed: this function returns a half-open interval value 
                                      al.QueryBases,
                                      al.Qualities,
                                      al.CigarData,
                                      md_tag_contents);

    }

    /**
     * @brief Given a single SAM alignment entry (assumed mapped and split into fields) -
     *        parse mutations, reconstruct target sequence over aligned region,
     *        and merge ambiguously aligned indels with alternate placements.
     *
     * @param fields Single SAM alignment, already split by tab character into vector
     * @return Returns the information needed to count RT mutations and depths or
     * count variants: Left-most position in alignment target sequence (0-based),
     * right-most position (0-based), reconstructed target sequence over this range,
     * and mutations present in this read. Also return read id for debugging purposes.
     */
    boost::tuple<std::string, int, int, std::string, std::string, std::vector<Mutation>>
    parseSamMutations(const std::vector<std::string> &fields) {
        // TODO: test whether reads mapped in the reverse sense are parsed correctly
        int left_target_pos;
        try {
            left_target_pos = stoi(fields[3])-1;
        } catch (std::exception &err) {
            throw std::runtime_error("Error: line is incorrectly formatted (couldn't parse mapped location).");
        }

        std::string read_id = fields[0];
        std::string query_bases = fields[9];
        std::string query_qual = fields[10];
        std::string cigar_string = fields[5];

        std::vector<BamTools::CigarOp> cigar_data = detail::parseCIGAR(cigar_string);

        int right_target_pos = detail::calcRightTargetPos(left_target_pos,
                                                          cigar_data);

        const std::string tag = "MD";
        std::string md_tag_contents;
        if (!detail::getSamTag(fields, tag, md_tag_contents)) {
            throw std::runtime_error("Error: no MD tag in alignment.");
        }

        std::string local_target_seq;
        std::string local_target_qual;
        std::vector <Mutation> adjusted_mutations;
        boost::tie(left_target_pos,
                   right_target_pos,
                   local_target_seq,
                   local_target_qual,
                   adjusted_mutations) = detail::parseMutations(left_target_pos,
                                                                right_target_pos,
                                                                query_bases,
                                                                query_qual,
                                                                cigar_data,
                                                                md_tag_contents);

        return boost::make_tuple(read_id,
                                 left_target_pos,
                                 right_target_pos,
                                 local_target_seq,
                                 local_target_qual,
                                 adjusted_mutations);

    }

    /**
     * @brief Parse mutations and reconstruct target sequences for mapped reads in
     *        a given BAM file. Write mutations and local target sequences to file.
     */
    void parseBAM(const std::string &filename,
                  const std::string &outname,
                  unsigned int min_mapq = DEFAULT_MIN_MAPQ) {

        // FIXME: support warn_on_no_mapped flag here too

        BAM::BamReader reader;
        if (!reader.Open(filename)) {
            std::string error = "Error: Could not open " + filename + " for BAM reading.";
            std::cout << error << std::endl;
            throw std::runtime_error(error);
        }

        // TODO: support '-' arg for writing to stdout instead of file, see http://stackoverflow.com/questions/366955/obtain-a-stdostream-either-from-stdcout-or-stdofstreamfile
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

        BAM::BamAlignment al;

        while (reader.GetNextAlignment(al)) {

            // skip unmapped reads
            if (al.RefID == -1) {
                //std::cout << "Skipping unmapped read." << std::endl;
                continue;
            }

            // skip reads below mapping quality threshold
            if (al.MapQuality < min_mapq){
                continue;
            }

            std::string read_id = al.Name;

            int left_target_pos;
            int right_target_pos;
            std::string local_target_seq;
            std::string local_target_qual;
            std::vector<Mutation> adjusted_mutations;

            boost::tie(left_target_pos,
                       right_target_pos,
                       local_target_seq,
                       local_target_qual,
                       adjusted_mutations) = mutation_parser::parseMutations(al);


            std::string s = serializeReadInfo(read_id,
                                              left_target_pos,
                                              right_target_pos,
                                              local_target_seq,
                                              local_target_qual,
                                              adjusted_mutations);

            out << s;
        }

        out << std::flush;
    }

    // WIP: SAM file support (for cases like bowtie2 output). Could use lib or just roll my own quick since it's tab-delimited.

    /**
     * @brief Parse mutations and reconstruct target sequences for mapped reads in
     *        a given SAM file. Write mutations and local target sequences to file.
     */
    void parseSAM(const std::string &filename,
                  const std::string &outname,
                  unsigned int min_mapq = DEFAULT_MIN_MAPQ,
                  bool warn_on_no_mapped = false) {

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
        size_t c = 0;

        while (std::getline(in, line)) {
            // skip headers
            if (line.length()<1 or line[0]=='@'){
                continue;
            }

            // split into fields
            std::vector<std::string> fields;
            std::string trimmed = boost::trim_copy(line);
            boost::split(fields, trimmed, boost::is_any_of("\t"), boost::token_compress_off);

            if (fields.size() < 11){
                throw std::runtime_error("Error: unable to parse incomplete line.");
            }

            // skip unmapped reads
            if (fields[2]=="*"){
                continue;
            }

            // skip reads below mapping quality threshold
            if (std::stoi(fields[4])<min_mapq){
                continue;
            }

            std::string read_id;
            int left_target_pos;
            int right_target_pos;
            std::string local_target_seq;
            std::string local_target_qual;
            std::vector<Mutation> adjusted_mutations;

            boost::tie(read_id,
                       left_target_pos,
                       right_target_pos,
                       local_target_seq,
                       local_target_qual,
                       adjusted_mutations) = mutation_parser::parseSamMutations(fields);

            std::string s = serializeReadInfo(read_id,
                                              left_target_pos,
                                              right_target_pos,
                                              local_target_seq,
                                              local_target_qual,
                                              adjusted_mutations);

            out << s;

            c++;
            /*if (c%100000==0){
                std::cout << "Counted "<<c<<" reads" <<std::endl;
            }*/
        }

        out << std::flush;

        if (c < 1) {
            if (warn_on_no_mapped) {
                std::cout << "WARNING: Input file " + filename + " contains no mapped reads." << std::endl;
            } else {
                throw std::runtime_error("ERROR: Input file " + filename + " contains no mapped reads.");
            }
        }

    }


}
