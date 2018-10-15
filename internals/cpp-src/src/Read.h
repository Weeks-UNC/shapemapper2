/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2018 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#ifndef SHAPEMAPPER_READ_H
#define SHAPEMAPPER_READ_H

#include <algorithm>

#include "Mutation.h"
#include "MutationProcessing.h"
#include "PrimerPair.h"
#include "util.h"

namespace M = mutation;

enum {
    FORWARD, REVERSE, UNSPECIFIED_STRAND
};

const std::vector<std::string> strands = {
    "FORWARD", "REVERSE", "UNSPECIFIED_STRAND"
};

enum {
    READ1, READ2
};

enum {
    R1, R2
};

/*
 * generally, these indicate read type *at the time of alignment*
 * - so MERGED indicates a read that was a result of mate pair
 *   merging before alignment
 * - PAIRED indicates a read that was aligned as one of two mate pairs
 *   and merged after alignment
 * - however, UNPAIRED_R1/R2 indicates a read that was provided to the
 *   aligner as a pair, but aligned discordantly or has a missing mate
 * 
 */


enum {
    PAIRED_R1,
    PAIRED_R2,
    UNPAIRED_R1,
    UNPAIRED_R2,
    UNPAIRED,
    MERGED,
    PAIRED,
    UNSPECIFIED_READ_TYPE
};


const std::vector<std::string> read_types = {
    "PAIRED_R1",
    "PAIRED_R2",
    "UNPAIRED_R1",
    "UNPAIRED_R2",
    "UNPAIRED",
    "MERGED",
    "PAIRED",
    "UNSPECIFIED_READ_TYPE"
};

enum {
    INCLUDED,
    LOW_MAPQ,
    OFF_TARGET,
    UNMAPPED,
};

std::vector<std::string> mapping_categories = {
    "INCLUDED",
    "LOW_MAPQ",
    "OFF_TARGET",
    "UNMAPPED",
};

enum {
    NO_ASSOCIATED_PRIMER_PAIR = -999,
};

enum {
    RIGHT, LEFT
};

// FIXME: document coordinate conventions used here
class Read {
public:
    int left;
    int right;
    int strand;
    int read_type;
    int mapping_category;
    /* ^ for plotting purposes indicating:
        included, low mapping quality, or off-target
        */
    int primer_pair; // negative indicates no associated primer pair
    std::string id;
    std::string seq;
    std::string qual;
    std::vector<bool> mapped_depth; // simple end-to-end read depths. excludes missing regions between mate pairs. does not exclude primer regions.
    std::vector<bool> depth;
    std::vector<bool> count;
    std::vector<Mutation> mutations;

    Read();
    Read(const int left_,
          const int right_,
          const std::string &seq_);
    Read(const std::string &serialized);
    Read& setLeft(const int left_);
    Read& setRight(const int right_);
    Read& setStrand(const int strand_);
    Read& setMutations(const std::vector<Mutation> &mutations_);
    Read& setReadType(const int read_type_);
    Read& setMappingCategory(const int mapping_category_);
    Read& setMappingCategory(const std::string &mapping_category_);
    Read& setPrimerPair(const int primer_pair_);
    Read& setSeq(const std::string &seq_);
    Read& setQual(const std::string &qual_);
    Read& setMappedDepth(const std::vector<bool> &mapped_depth_);
    Read& setDepth(const std::vector<bool> &depth_);
    Read& setCount(const std::vector<bool> &count_);
    Read& setId(const std::string &id_);
    Read& trimRightEnd(const int exclude_3prime);
    Read& stripPrimers(const PrimerPair& primer_pair);
    Read& shiftAmbigIndels(const bool right_align_ambig_dels,
                           const bool right_align_ambig_ins);
    Read& collapseMutations(const int max_internal_match);
    Read& classifyMutations();
    Read& filterQscoresCountDepths(const int min_qual,
                                   const std::string &mutation_type,
                                   const bool variant_mode);
    std::string toString() const;
    std::string serializeMutations() const;
    std::string serializeForTest() const;
};

// FIXME: simplify this to something like https://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Named_Parameter
// - don't use all this inline garbage?

Read::Read()
        : left (0)
        , right (0)
        , strand (UNSPECIFIED_STRAND)
        , mapping_category (INCLUDED)
        , primer_pair (NO_ASSOCIATED_PRIMER_PAIR)
        , read_type (UNSPECIFIED_READ_TYPE)
        , id ("")
        , seq ("")
        , qual ("")
        , mutations {}
        , depth {}
        , count {}
        , mapped_depth {}
{ }

Read::Read(const int left_, // 0-based
           const int right_, // 0-based
           const std::string &seq_)
        : left (left_)
        , right (right_)
        , strand (UNSPECIFIED_STRAND)
        , mapping_category (INCLUDED)
        , primer_pair (NO_ASSOCIATED_PRIMER_PAIR)
        , seq (seq_)
        , id ("")
        , mutations {}
        , read_type (UNSPECIFIED_READ_TYPE)
        , qual ("")
        , depth {}
        , count {}
        , mapped_depth {}
{ }

// used in testMutationParser to read debug outputs back in
Read::Read(const std::string &serialized){
    std::vector<std::string> fields;
    std::string trimmed = boost::trim_right_copy_if(serialized, boost::is_any_of("\n\r"));
    boost::split(fields, trimmed, boost::is_any_of("\t"));
    read_type = UNSPECIFIED_READ_TYPE;
    primer_pair = NO_ASSOCIATED_PRIMER_PAIR;
    mapping_category = INCLUDED;
    if (fields.at(1) != "") {
        read_type = util::indexOf(read_types, fields.at(1));
    }
    left = stoi(fields.at(2));
    right = stoi(fields.at(3));
    strand = UNSPECIFIED_STRAND;
    if (fields.at(4) == "+"){ strand = FORWARD; }
    if (fields.at(4) == "-"){ strand = REVERSE; }
    mapping_category = util::indexOf(mapping_categories, fields.at(5));
    primer_pair = std::stoi(fields.at(6));
    seq = fields.at(7);
    qual = fields.at(8);
    try {
        mapped_depth = util::stringToBoolVect(fields.at(9));
        depth = util::stringToBoolVect(fields.at(10));
        count = util::stringToBoolVect(fields.at(11));
        mutations = stringToMutationVect(fields.at(12));
    } catch (const std::out_of_range& ex) { }
}

Read&
Read::setStrand(const int strand_){
    strand = strand_;
    return *this;
}

Read&
Read::setLeft(const int left_){
    left = left_;
    return *this;
}

Read&
Read::setRight(const int right_){
    right = right_;
    return *this;
}

Read&
Read::setPrimerPair(const int primer_pair_){
    primer_pair = primer_pair_;
    return *this;
}

Read&
Read::setMutations(const std::vector<Mutation> &mutations_){
    mutations.clear();
    if (mutations_.size() > 0) {
        for (auto &m : mutations_) {
            mutations.push_back(Mutation(m));
        }
    }
    return *this;
}

Read&
Read::setReadType(const int read_type_){
    read_type = read_type_;
    return *this;
}

Read&
Read::setId(const std::string &id_){
    id = id_;
    return *this;
}

Read&
Read::setMappingCategory(const int mapping_category_){
    mapping_category = mapping_category_;
    return *this;
}

Read&
Read::setMappingCategory(const std::string &mapping_category_){
    mapping_category = util::indexOf(mapping_categories, mapping_category_);
    return *this;
}

Read&
Read::setSeq(const std::string &seq_){
    seq = seq_;
    return *this;
}

Read&
Read::setQual(const std::string &qual_){
    qual = qual_;
    return *this;
}

Read&
Read::setMappedDepth(const std::vector<bool> &mapped_depth_){
    mapped_depth.clear();
    if (mapped_depth_.size() > 0) {
        for (auto b : mapped_depth_) {
            mapped_depth.push_back(b);
        }
    }
    return *this;

}

Read&
Read::setDepth(const std::vector<bool> &depth_){
    depth.clear();
    if (depth_.size() > 0) {
        for (auto b : depth_) {
            depth.push_back(b);
        }
    }
    return *this;
}

Read&
Read::setCount(const std::vector<bool> &count_){
    count.clear();
    if (count_.size() > 0) {
        for (auto b : count_) {
            count.push_back(b);
        }
    }
    return *this;
}

/*
 * @brief For debugging and mutation rendering
 */
std::string
Read::toString() const {
    std::string s = "";
    s += "[read]\t";
    s += read_types[read_type] + '\t';
    s += std::to_string(left) + '\t';
    s += std::to_string(right) + '\t';
    std::string strand_str = "N/A";
    if (strand == FORWARD) { strand_str = "+"; }
    if (strand == REVERSE) { strand_str = "-"; }
    s += strand_str + '\t';
    s += mapping_categories.at(mapping_category) + '\t';
    s += std::to_string(primer_pair) + '\t';
    s += seq + '\t';
    s += qual + '\t';
    s += util::toString(mapped_depth) + '\t';
    s += util::toString(depth) + '\t';
    s += util::toString(count) + '\t';
    s += mutation::toString(mutations);
    s += "\n";
    return s;
}

/*
 * @brief For mutation_parser file output
 */
std::string
Read::serializeMutations() const {
    using std::to_string;
    std::string sep = "\t";
    std::string o = read_types.at(read_type) + sep +
                    id + sep +
                    to_string(left) + sep +
                    to_string(right) + sep +
                    mapping_categories.at(mapping_category) + sep +
                    to_string(primer_pair) + sep +
                    util::toString(mapped_depth) + sep +
                    util::toString(depth) + sep +
                    util::toString(count) + sep +
                    mutation::toString(mutations) + "\n";
    return o;
}

/*
 * @brief For testing purposes (can probably remove)
 */
std::string
Read::serializeForTest() const {
    using std::to_string;
    std::string o =
                    id + "\t" +
                    to_string(left) + "\t" +
                    to_string(right) + "\t" +
                    seq + "\t" +
                    qual + "\t" +
                    mutation::toString(mutations) + "\n";
    return o;
}

// overload insertion operator to allow simple Read "pipe" to ofstream
std::ostream & operator << (std::ostream &out, const Read &r)
{
    out << r.toString();
    return out;
}




// end of Read member functions
//-----------------------------------------------------------------------------------------


// FIXME: unify or clarify various deserialization funcs. have several funcs with slight variations for historical reasons


/*
 * @brief For debugging and mutation rendering
 */
Read
parseDebugRead(const std::string &line) {
    // FIXME: move this to the Read(string) constructor and update tests with newer format
    //[read]	UNPAIRED	418	611	-	INCLUDED	-999	AGTCTGGTGGTGGGTCGTAAGTTTAGGAGGTGACTGCATCCTCCAGCATCTCAACTCCGTCTGTCTACTGTGTGAGACTTCGGCGGACCATTAGGAATGAGATCCGTGAGATCCTTCCATCTTCTTGAAGTCGCCTTTAGGGTGGCTGCGAGGTAGAGGGTTGGGGGTTGGTGGGCTGTCACGGAGCGACTGTC	B=59DDFFFFFFFFFFFFEFFF<FFFFFFGGGGFFG?FAEFFFFGFCGGFCFFDGGEDGGFGGGGFGGGGGGFGFGFDGGGGGGGGGGFGGGGGGGFFGGGGF>GGFF<GGFFFFFGFF@GFGCGGGEGFGGFFAGGGGGGGGGFDGEDGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGG	11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
    std::vector<std::string> fields;
    std::string trimmed = boost::trim_right_copy_if(line, boost::is_any_of("\n\r"));
    boost::split(fields, trimmed, boost::is_any_of("\t"));
    int read_type = UNSPECIFIED_READ_TYPE;
    if (fields.at(1) != "") {
        read_type = util::indexOf(read_types, fields.at(1));
    }
    int left = stoi(fields.at(2));
    int right = stoi(fields.at(3));
    int strand = UNSPECIFIED_STRAND;
    if (fields.at(4) == "+"){ strand = FORWARD; }
    if (fields.at(4) == "-"){ strand = REVERSE; }
    int mapping_category = util::indexOf(mapping_categories, fields.at(5));
    int primer_pair = stoi(fields.at(6));
    std::string seq = fields.at(7);
    std::string qual = fields.at(8);
    std::vector<bool> mapped_depth;
    std::vector<bool> depth;
    std::vector<bool> count;
    std::vector<Mutation> mutations;
    try {
        mapped_depth =  util::stringToBoolVect(fields.at(9));
        depth = util::stringToBoolVect(fields.at(10));
        count = util::stringToBoolVect(fields.at(11));
        mutations = stringToMutationVect(fields.at(12));
    } catch (const std::out_of_range& ex) { }
    Read r = Read();
    r.setLeft(left)
            .setRight(right)
            .setStrand(strand)
            .setReadType(read_type)
            .setSeq(seq)
            .setQual(qual)
            .setDepth(depth)
            .setMutations(mutations)
            .setCount(count)
            .setMappingCategory(mapping_category)
            .setPrimerPair(primer_pair)
            .setMappedDepth(mapped_depth);
    return r;
}


Read
parseTestRead(const std::string &line) {
    /*
     * parse debug outputs
     */
    // example format
    //M00236:2:000000000-A21YG:1:1106:15774:10066	0	136	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	GGCCTTCGGGCCAAGGACTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGTATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATCCAAATCGGGCTTCGGTCCGGTTC	32 33 "TCTTTC" "TCTTTC" "" 32 34 "T" "T" "" 82 84 "C" "C" "" 84 86 "A" "A" "" 114 118 "GA" "GA" ""
    std::vector<std::string> fields;
    std::string trimmed = boost::trim_right_copy_if(line, boost::is_any_of("\n\r"));
    boost::split(fields, trimmed, boost::is_any_of("\t"));
    std::string id = fields.at(0);
    int read_type = UNPAIRED;
    int left = stoi(fields.at(1));
    int right = stoi(fields.at(2));
    int strand = FORWARD;
    std::string seq = fields.at(3);
    std::string qual = fields.at(4);
    std::vector<Mutation> mutations = {};
    try {
        mutations = stringToMutationVect(fields.at(5));
    } catch (const std::out_of_range& ex) { }
    Read r(left, right, seq);
    r.setMutations(mutations)
            .setQual(qual)
            .setReadType(read_type)
            .setId(id)
            .setStrand(strand);
    return r;
}

Read
mergeMatePairsSimple(const std::vector<Read> &reads,
                     const bool debug = false){
    /*
     * Simple merge of two mate paired reads, without considering sequence or mutations.
     * Used to allow read depth counting for off-target and low mapping quality reads.
     */
    Read simple_merged{}; // explicitly call default constructor
    std::vector<bool> mapped_depth;

    int fw_read_index = R1;
    int rv_read_index = R2;
    if (reads[R1].strand == REVERSE and reads[R2].strand == FORWARD) {
        fw_read_index = R2;
        rv_read_index = R1;
    }

    int left = std::min(reads.at(fw_read_index).left, reads.at(rv_read_index).left);
    int right = std::max(reads.at(rv_read_index).right, reads.at(fw_read_index).right);
    mapped_depth.resize(right - left + 1, 0);
    if (debug) {
        std::cout << "mapped_depth.size() after resize: " << mapped_depth.size() << "\n" << std::flush;
    }
    int l1 = reads.at(fw_read_index).right - reads.at(fw_read_index).left + 1;
    int l2 = reads.at(rv_read_index).right - reads.at(rv_read_index).left + 1;
    std::fill(mapped_depth.begin(), mapped_depth.begin() + l1, 1);
    if (debug) {   
        std::cout << "filled fw read depth\n" << std::flush;
    }
    std::fill(mapped_depth.end() - l2, mapped_depth.end(), 1);
    if (debug) {
        std::cout << "filled fw read depth\n" << std::flush;
    }
    simple_merged.setReadType(PAIRED)
                 .setId(reads[0].id)
                 .setLeft(left)
                 .setRight(right)
                 .setMappedDepth(mapped_depth);
    return simple_merged;
}




#endif //SHAPEMAPPER_READ_H
