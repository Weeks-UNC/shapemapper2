/** @file
 * @brief Mutation class and several related methods. 
 *        Used by MutationParser and MutationCounter.
 */

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2017 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#ifndef SHAPEMAPPER_MUTATION_H_H
#define SHAPEMAPPER_MUTATION_H_H

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>

namespace mutation {
    const std::vector<std::string> mutation_classes = {
            "A-", "T-", "G-", "C-",
            "-A", "-T", "-G", "-C", "-N",
            "AT", "AG", "AC",
            "TA", "TG", "TC",
            "GA", "GT", "GC",
            "CA", "CT", "CG",
            "multinuc_deletion",
            "multinuc_insertion",
            "multinuc_mismatch",
            "complex_deletion",
            "complex_insertion",
    };

    // FIXME: move to a utility header
    std::string toString(const std::vector<bool> &bools){
        std::string s;
        for (auto b : bools){
            s += std::to_string(b);
        }
        return s;
    }

    /**
     * @brief Stores deviations from an alignment target between two nucleotide
     *        positions.
     */
    class Mutation {
        // TODO: refactor functions from MutationParser to be member funcs within Mutation class?
    public:

        /** @brief Leftmost unchanged alignment target nucleotide (0-based) */
        int left;
        /** @brief Rightmost unchanged alignment target nucleotide (0-based) */
        int right;
        /** @brief Read sequence replacing alignment target sequence between left and right (exclusive) */
        std::string seq;
        /** @brief Basecall quality scores (ASCII encoded Phred scores) for nucs in seq */
        std::string qual;

        Mutation() { };

        Mutation(const int l, const int r, const std::string &s, const std::string &q) {
            left = l;
            right = r;
            seq = s;
            qual = q;
        }

        Mutation(const Mutation &m) {
            left = m.left;
            right = m.right;
            seq = m.seq;
            qual = m.qual;
        }

        // need a "<" operator so std::map knows how to compare objects
        bool operator<(const Mutation &a) const {
            if (left < a.left) {
                return true;
            }
            if (a.left < left) {
                return false;
            }
            if (right < a.right) {
                return true;
            }
            if (a.right < right) {
                return false;
            }
            if (seq < a.seq) {
                return true;
            }
            if (a.seq < seq) {
                return false;
            }
            if (qual < a.qual) {
                return true;
            }
            if (a.qual < qual) {
                return false;
            }
            return false;
        }

        /*bool operator==(const Mutation &a) const{
            return (left==a.left and right==a.right and seq==a.seq);
        }*/

        std::string toString() const {
            using std::to_string;
            return to_string(left) + ' ' + to_string(right) + ' ' + '"' + seq + '"' + ' ' + '"' + qual + '"';
        }

        bool isSimpleInsert() const {
            return right - left == 1;
        }

        bool isSimpleGap() const {
            return seq.length() == 0;
        }

        bool isGapOrInsert() const {
            return seq.length() != right - left - 1;
        }

        bool isGap() const {
            return seq.length() < right - left - 1;
        }

        bool isInsert() const {
            return seq.length() > right - left - 1;
        }

        /**
         * @brief Check if a mutation was previously detected as ambiguously aligned.
         *        Assumes adjustAmbiguousMutations() was previously run.
         */
        bool isAmbiguous() const {
            int d = right - left - 1;
            return ((d > seq.length() and seq.length() > 0) or
                    (d < seq.length() and d > 0));
        }

        /**
         * @brief Classify a mutation. Examples: "AG" A in target, G in read.
         *        "-C" insert of C. "A-" deletion of A. Others: "multinuc_deletion",
         *        "multinuc_insertion", "multinuc_mismatch", "complex_insertion/complex_deletion" (other multinuc mutation)
         */
        std::string classify(const std::string &local_target,
                             const int target_pos) const {
            int d = right - left - 1;
            std::string s;
            try {
                if (d == 1 and seq.length() == 0) {
                    // deletion of a single nucleotide
                    int i = left + 1 - target_pos;
                    s = local_target.substr(left + 1 - target_pos, 1) + '-';
                } else if (d == 0 and seq.length() == 1) {
                    // insertion of a single nucleotide
                    s = '-' + seq;
                } else if (d == 1 and seq.length() == 1) {
                    // single-nucleotide mismatch
                    if (seq == "N"){
                        s = "N_match";
                    } else {
                        s = local_target[left + 1 - target_pos] + seq;
                    }
                } else if (d > 1 and seq.length() == 0) {
                    s = "multinuc_deletion";
                } else if (d == 0 and seq.length() > 1) {
                    s = "multinuc_insertion";
                } else if (d == seq.length()) {
                    s = "multinuc_mismatch";
                } else if (seq.length() < d) {
                    s = "complex_deletion";
                } else if (seq.length() > d) {
                    s = "complex_insertion";
                } else {
                    throw std::runtime_error("Error: Unknown error. Unable to classify mutation.");
                }
            } catch (std::out_of_range &err) {
                throw std::runtime_error(
                        "Error: Unable to classify mutation. Mutation location falls outside local target sequence.");
            }
            return s;
        }

    };

    std::string toString(const std::vector<Mutation> &m) {
        if (m.size() < 1) {
            return "";
        }
        std::string o;
        for (std::vector<Mutation>::const_iterator it = m.begin();
             it != m.end();
             ++it) {
            o += it->toString();
            if (it + 1 != m.end()) {
                o += ' '; // add separator except after last item
            }
        }
        return o;
    }

    /**
    * @brief Stores mutations with additional information
    */
    class TaggedMutation : public Mutation {

    public:

        /** @brief Mutation tag (usually mutation classification)*/
        std::string tag;
        /** @brief Whether a mutation is or is derived from an ambiguous alignment */
        bool ambig;

        TaggedMutation(const TaggedMutation &m) {
            left = m.left;
            right = m.right;
            seq = m.seq;
            qual = m.qual;
            tag = m.tag;
            ambig = m.ambig;
        }

        TaggedMutation(const int l,
                       const int r,
                       const std::string &s,
                       const std::string &q,
                       const std::string &t) {
            left = l;
            right = r;
            seq = s;
            qual = q;
            tag = t;
            ambig = false;
        }

        TaggedMutation(const int l,
                       const int r,
                       const std::string &s,
                       const std::string &q,
                       const std::string &t,
                       const bool a) {
            left = l;
            right = r;
            seq = s;
            qual = q;
            tag = t;
            ambig = a;
        }

        TaggedMutation(const Mutation &m,
                       const std::string &t) {
            left = m.left;
            right = m.right;
            seq = m.seq;
            qual = m.qual;
            tag = t;
            ambig = false;
        }

        TaggedMutation(const Mutation &m,
                       const std::string &t,
                       const bool a) {
            left = m.left;
            right = m.right;
            seq = m.seq;
            qual = m.qual;
            tag = t;
            ambig = a;
        }

        std::string toString() const {
            using std::to_string;
            std::string s = to_string(left) + ' ' + to_string(right) + " \"";
            s += seq + "\" \"" + qual + "\" \"" + tag;
            if (ambig){
                s += "_ambig";
            }
            s += '"';
            return s;
        }
    };

    std::string toString(const std::vector<TaggedMutation> &m) {
        if (m.size() < 1) {
            return "";
        }
        std::string o;
        for (std::vector<TaggedMutation>::const_iterator it = m.begin();
             it != m.end();
             ++it) {
            o += it->toString();
            if (it + 1 != m.end()) {
                o += ' '; // add separator except after last item
            }
        }
        return o;
    }

    /**
     * @brief Return a copied vector of mutations with previously identified ambiguously aligned
     *        indels realigned left or right.
     */
    std::vector<TaggedMutation> shiftAmbigIndels(const std::vector<Mutation> &mutations,
                                                 const std::string &local_target_seq,
                                                 const std::string &local_target_qual,
                                                 const int left_target_pos,
                                                 const bool right_align_ambig_dels,
                                                 const bool right_align_ambig_ins) {

        std::vector<TaggedMutation> adjusted_mutations;
        if (mutations.size() > 0) {
            for (std::vector<Mutation>::const_iterator it = mutations.begin();
                 it != mutations.end();
                 ++it) {
                int left = (*it).left;
                int right = (*it).right;
                std::string seq = (*it).seq;
                std::string qual = (*it).qual;
                std::string tag = "";
                bool ambig = true;
                if ((*it).isAmbiguous()){
                    // note: a bit ugly as a result of handling multinuc ambiguous mutations
                    //      with "internal" mismatches. checking for mismatches and
                    //      pushing back additional mutations if needed after shifting.
                    if ((*it).isGap()) {
                        if (right_align_ambig_dels){
                            left = (*it).left + (*it).seq.length();
                            seq = "";
                            qual = "";
                            // check for mismatches and create mutations if needed
                            for (int i=0; i<(*it).seq.length(); ++i){
                                int local_left = (*it).left - left_target_pos + 1 + i;
                                std::string c = (*it).seq.substr(i,1);
                                std::string q = (*it).qual.substr(i,1);
                                //std::cout << "c: " << c << std::endl;
                                if ( c != local_target_seq.substr(local_left, 1)){
                                    TaggedMutation m((*it).left + i,
                                                     (*it).left + i + 2,
                                                     c, q, tag, ambig);
                                    adjusted_mutations.push_back(m);
                                }
                            }
                            TaggedMutation m(left, right, seq, qual, tag, ambig);
                            adjusted_mutations.push_back(m);
                        } else {
                            right = (*it).right - (*it).seq.length();
                            seq = "";
                            qual = "";
                            TaggedMutation m(left, right, seq, qual, tag, ambig);
                            adjusted_mutations.push_back(m);
                            // check for mismatches and create mutations if needed
                            for (int i=0; i<(*it).seq.length(); ++i){
                                int local_left = right - left_target_pos + i;
                                std::string c = (*it).seq.substr(i,1);
                                std::string q = (*it).qual.substr(i,1);
                                if ( c != local_target_seq.substr(local_left, 1)){
                                    TaggedMutation m(right + i - 1,
                                                     right + i + 1,
                                                     c, q, tag, ambig);
                                    adjusted_mutations.push_back(m);
                                }
                            }
                        }
                    } else if ((*it).isInsert()) {
                        if (right_align_ambig_ins){
                            int d = (*it).seq.length() - ((*it).right - (*it).left - 1);
                            left = (*it).right-1;
                            // take right-most portion of seq
                            seq = (*it).seq.substr((*it).seq.length()-d);
                            qual = (*it).qual.substr((*it).qual.length()-d);
                            // check for mismatches and create mutations if needed
                            for (int i=0; i<(*it).seq.length()-d; ++i){
                                int local_left = (*it).left - left_target_pos + 1 + i;
                                std::string c = (*it).seq.substr(i,1);
                                std::string q = (*it).qual.substr(i,1);
                                if ( c != local_target_seq.substr(local_left, 1)){
                                    TaggedMutation m((*it).left+i,
                                                     (*it).left+i+2,
                                                      c, q, tag, ambig);
                                    adjusted_mutations.push_back(m);
                                }
                            }
                            TaggedMutation m(left, right, seq, qual, tag, ambig);
                            adjusted_mutations.push_back(m);
                        } else {
                            int d = (*it).seq.length() - ((*it).right - (*it).left - 1);
                            right = (*it).left+1;
                            // take left-most portion of seq
                            seq = (*it).seq.substr(0, d);
                            qual = (*it).qual.substr(0, d);
                            TaggedMutation m(left, right, seq, qual, tag, ambig);
                            adjusted_mutations.push_back(m);
                            // check for mismatches and create mutations if needed
                            for (int i=0; i<(*it).seq.length()-d; ++i) {
                                int local_left = (*it).left - left_target_pos + 1 + i;
                                std::string c = (*it).seq.substr(d + i, 1);
                                std::string q = (*it).qual.substr(d + i, 1);
                                if (c != local_target_seq.substr(local_left, 1)) {
                                    TaggedMutation m((*it).left + i,
                                                     (*it).left + i + 2,
                                                     c, q, tag, ambig);
                                    adjusted_mutations.push_back(m);
                                }
                            }
                        }

                    }
                } else {
                    ambig = false;
                    TaggedMutation m(left, right, seq, qual, tag, ambig);
                    adjusted_mutations.push_back(m);
                }


            }
        }
        return adjusted_mutations;
    }

    /**
     * @brief Return a copied vector of mutations with 3-prime mutations within
     *        exclude_3prime nucleotides removed.
     */
    template<typename T>
    std::vector<T> strip_3prime(const std::vector<T> &mutations,
                                const std::string &local_target_seq,
                                const std::string &local_target_qual,
                                const int left_target_pos,
                                const int exclude_3prime) {
        std::vector<T> stripped_mutations;
        if (mutations.size() > 0) {
            int max_right = left_target_pos + local_target_seq.length() - exclude_3prime - 1;
            for (std::vector<Mutation>::const_iterator it = mutations.begin();
                 it != mutations.end();
                 ++it) {
                 if (((*it).right-1) <= max_right){
                     stripped_mutations.push_back(T((*it)));
                 }
            }
        }
        return stripped_mutations;
    }

    /**
     * @brief Return a copied vector of mutations with adjacent mutations combined.
     *        Assumes mutations are sorted left-to-right. Mutations with up to
     *        max_internal_match unchanged nucleotides between them will be combined.
     */
    std::vector<TaggedMutation> collapseMutations(const std::vector<TaggedMutation> &mutations,
                                                  const int max_internal_match,
                                                  const std::string &local_target_seq,
                                                  const std::string &local_target_qual,
                                                  const int left_target_pos) {
        std::vector<TaggedMutation> collapsed_mutations;
        std::vector<TaggedMutation> unmerged_mutations; // temp storage for ambiguous mutations (if any),
                                                        // to avoid merging with other mutations
        if (mutations.size() > 0) {
            for (std::vector<TaggedMutation>::const_iterator it = mutations.begin();
                 it != mutations.end();
                 ++it) {
                TaggedMutation m(*it);
                // if a mutation is not really a mutation (N match),
                // it should not be merged with other mutations
                if (m.tag=="N_match") {
                    unmerged_mutations.push_back(m);
                } else if (collapsed_mutations.size() > 0 and
                           collapsed_mutations.back().right - 1 == m.left) {
                    // adjacent mutation
                    collapsed_mutations.back().right = m.right;
                    collapsed_mutations.back().seq += m.seq;
                    collapsed_mutations.back().qual += m.qual;
                    collapsed_mutations.back().tag = "";
                    if (m.ambig) {
                        collapsed_mutations.back().ambig = true;
                    }
                } else if (collapsed_mutations.size() > 0  and
                           m.left - (collapsed_mutations.back().right - 1) <= max_internal_match){
                    // mutation within allowed distance for adjacency
                    collapsed_mutations.back().seq += local_target_seq.substr(collapsed_mutations.back().right-left_target_pos,
                                                                              m.left-collapsed_mutations.back().right+1) ;
                                                   // ^ pick up nucs between mutations
                    collapsed_mutations.back().qual += local_target_qual.substr(collapsed_mutations.back().right-left_target_pos,
                                                                                m.left-collapsed_mutations.back().right+1) ;
                    collapsed_mutations.back().right = m.right;
                    collapsed_mutations.back().seq += m.seq;
                    collapsed_mutations.back().qual += m.qual;
                    collapsed_mutations.back().tag = "";
                    if (m.ambig) {
                        collapsed_mutations.back().ambig = true;
                    }
                } else {
                    // non-adjacent mutation
                    collapsed_mutations.push_back(m);
                }
            }
            // iteratively strip any matches from mutation ends (these can appear as a result of
            // shifted indels). For left-aligned ambig gaps, this only affects the 5-prime end of
            // a mutation, so shouldn't affect reactivity profiles. Adjusting it anyway for
            // consistency between aligners.
            if (collapsed_mutations.size() > 0){
                for (std::vector<TaggedMutation>::iterator it = collapsed_mutations.begin();
                it != collapsed_mutations.end();
                ++it){
                    // from left end
                    int new_left = (*it).left;
                    std::string new_seq = (*it).seq;
                    std::string new_qual = (*it).qual;
                    try{
                        for (int i = 0; i < (*it).seq.length(); ++i){
                            std::string c = (*it).seq.substr(i,1);
                            std::string q = (*it).qual.substr(i,1);
                            if ((*it).left + 1 + i >= (*it).right) { break; } // don't go outside of region mutation spans
                            int p = (*it).left + 1 + i - left_target_pos;
                            if (p<0) { break; }
                            std::string r = local_target_seq.substr(p, 1);
                            if (c==r){
                                ++new_left;
                                new_seq = new_seq.substr(1);
                                new_qual = new_qual.substr(1);
                            } else {
                                break;
                            }
                        }
                    }
                    catch( std::out_of_range& exception){
                        break;
                    }
                    (*it).left = new_left;
                    (*it).seq = new_seq;
                    (*it).qual = new_qual;

                    // from right end
                    int new_right = (*it).right;
                    new_seq = (*it).seq;
                    new_qual = (*it).qual;
                    try{
                        for (int i = (*it).seq.length()-1; i >= 0; --i){
                            std::string c = (*it).seq.substr(i,1);
                            std::string q = (*it).qual.substr(i,1);
                            int d = (*it).seq.length() - i;
                            if ((*it).right - d <= (*it).left) { break; } // don't go outside of region mutation spans
                            int p = (*it).right - d - left_target_pos;
                            if (p<0){ break; }
                            std::string r = local_target_seq.substr(p, 1);
                            if (c==r){
                                --new_right;
                                new_seq = new_seq.substr(0,i);
                                new_qual = new_qual.substr(0,i);
                            } else {
                                break;
                            }
                        }
                    }
                    catch( std::out_of_range& exception){
                        break;
                    }
                    (*it).right = new_right;
                    (*it).seq = new_seq;
                    (*it).qual = new_qual;
                }
            }

            // put unmerged (mismatched Ns) mutations back into list
            collapsed_mutations.insert(collapsed_mutations.end(),
                                       unmerged_mutations.begin(),
                                       unmerged_mutations.end());
            // make sure unmerged are back in order
            std::sort(collapsed_mutations.begin(), collapsed_mutations.end());
        }
        return collapsed_mutations;
    }

    /**
     * @brief Assign a mutation classification to each mutation in vector.
     *        Return a copy with new tags.
     */
    std::vector<TaggedMutation> classifyMutations(std::vector<TaggedMutation> &mutations,
                                                  const std::string &local_target_seq,
                                                  const std::string &local_target_qual,
                                                  const int target_pos) {
        std::vector<TaggedMutation> classified_mutations;
        if (mutations.size() > 0) {
            for (std::vector<TaggedMutation>::iterator it = mutations.begin();
                 it != mutations.end();
                 ++it) {
                TaggedMutation m((*it));
                if (m.tag.length() == 0){
                    m.tag = m.classify(local_target_seq,
                                       target_pos);
                }
                classified_mutations.push_back(m);
            }
        }
        return classified_mutations;
    }

    /**
     * @brief Filter mutations containing or neighboring bad basecalls.
     *        Filter non-mutation positions for qscore-aware read depth calculation.
     *
     *        In variant mode, mutation spans will not be excluded from depth, and
     *        gaps will contribute to depth if the basecalls on either side are good.
     */
    boost::tuple<std::vector<bool>,
                 std::vector<bool>,
                 std::vector<TaggedMutation>,
                 std::vector<TaggedMutation>>
    filterQscoresCountDepths(const std::vector<TaggedMutation> &mutations,
                             const std::string &local_target_seq,
                             const std::string &local_target_qual,
                             const int left_target_pos,
                             const int exclude_3prime,
                             const int min_qual,
                             const std::string &mutation_type,
                             const bool variant_mode) {
        
        std::vector<TaggedMutation> included_mutations;
        std::vector<TaggedMutation> excluded_mutations;
        int len = local_target_seq.length() - exclude_3prime;
        len = std::max(len, 0);
        //std::cout << "len: " << len << std::endl;
        std::vector<bool> local_effective_depth(len, 0);

        // store mutation vector indices by left and right unchanged position
        // so we don't have to iterate over all mutations for each target position
        std::vector<int> left_mut_indices(len,-1);
        std::vector<int> right_mut_indices(len,-1);
        std::vector<bool> in_mutation(len,0);
        for (int i=0; i<mutations.size(); ++i){
            TaggedMutation m(mutations[i]);
            try {
                left_mut_indices.at(m.left-left_target_pos) = i;
            } catch (std::out_of_range &exception) { }
            try {
                right_mut_indices.at(m.right - left_target_pos) = i;
            } catch (std::out_of_range &exception) { }
            for (int n = m.left + 1 - left_target_pos;
                 n < m.right - left_target_pos;
                 ++n) {
                try {
                    in_mutation[n] = 1;
                } catch (std::out_of_range &exception) { }
            }
        }

        // go one basecall at a time, examining neighboring basecalls, and include or exclude from counts and depths
        // in first pass, skip over nucs in mutations (but do consider them as neighbors)
        for (int i=0;
             i < len;
             ++i){
            if (in_mutation[i]){
                continue;
            }
            bool bad_basecall = false;
            // look at self
            if (local_target_qual[i]-33 < min_qual){
                bad_basecall = true;
            }
            // look at neighbor on left
            if (not bad_basecall){
                if (right_mut_indices[i] != -1){
                    TaggedMutation m(mutations[right_mut_indices[i]]);
                    if (m.seq.length() > 0){
                        // look at right-most basecall in mutation
                        // - should this include all basecalls in mutation?
                        if (m.qual.back()-33 < min_qual){
                            bad_basecall = true;
                        }
                    } else {
                        // look at basecall on other side of gap
                        int n = m.left - left_target_pos;
                        try {
                            if (local_target_qual.at(n)-33 < min_qual){
                                bad_basecall = true;
                            }
                        } catch (std::out_of_range &exception) { }
                    }
                } else {
                    // look at basecall on left (not in a mutation)
                    try {
                        if (local_target_qual.at(i-1)-33 < min_qual) {
                            bad_basecall = true;
                        }
                    } catch (std::out_of_range &exception) { }
                }
            }
            // look at neighbor on right
            if (not bad_basecall){
                if (left_mut_indices[i] != -1){
                    TaggedMutation m(mutations[left_mut_indices[i]]);
                    if (m.seq.length() > 0){
                        // look at left-most basecall in mutation
                        // - should this include all basecalls in mutation?
                        if (m.qual.front()-33 < min_qual){
                            bad_basecall = true;
                        }
                    } else {
                        // look at basecall on other side of gap
                        int n = m.right - left_target_pos;
                        try {
                            if (local_target_qual.at(n)-33 < min_qual){
                                bad_basecall = true;
                            }
                        } catch (std::out_of_range &exception) { }
                    }
                } else {
                    // look at basecall on right (not in a mutation)
                    try {
                        if (local_target_qual.at(i+1)-33 < min_qual) {
                            bad_basecall = true;
                        }
                    } catch (std::out_of_range &exception) { }
                }
            }
            if (not bad_basecall){
                local_effective_depth[i] = 1;
            }
        }

        // in second pass, update depths for only nucs within mutations, 
        // and filter mutations by q-scores (and mutation_type, if given)
        //
        // mutation_type possible values: mismatch gap insert gap_multi insert_multi complex
        for (int i=0; i<mutations.size(); ++i) {
            TaggedMutation m(mutations[i]);
            bool bad_mutation = false;

            // limit to specific mutation type if given
            // FIXME: this relies on the template class being TaggedMutation
            std::vector<std::string> mismatch_tags = {"AT", "AG", "AC",
                                                      "TA", "TG", "TC",
                                                      "GA", "GT", "GC",
                                                      "CA", "CT", "CG", 
                                                      "multinuc_mismatch"};
            std::vector<std::string> insert_tags = {"-A", "-T", "-G", "-C", "-N"};
            std::vector<std::string> gap_tags = {"A-", "T-", "G-", "C-"};

            if (mutation_type.length() > 0) {
                if (mutation_type == "mismatch") {
                    if (std::find(mismatch_tags.begin(), mismatch_tags.end(), m.tag) == mismatch_tags.end()) {
                        bad_mutation = true;
                    } 
                } else if (mutation_type == "insert") {
                    if (std::find(insert_tags.begin(), insert_tags.end(), m.tag) == insert_tags.end()) {
                        bad_mutation = true;
                    } 
                } else if (mutation_type == "insert_multi") {
                    if (m.tag != "multinuc_insertion") {
                        bad_mutation = true;
                    }
                } else if (mutation_type == "gap") {
                    if (std::find(gap_tags.begin(), gap_tags.end(), m.tag) == gap_tags.end()) {
                        bad_mutation = true;
                    } 
                } else if (mutation_type == "gap_multi") {
                    if (m.tag != "multinuc_deletion") {
                        bad_mutation = true;
                    }
                } else if (mutation_type == "complex") {
                    if (m.tag != "complex_deletion" && m.tag != "complex_insertion") {
                        bad_mutation = true;
                    }
                }
            }

            // look at basecalls within mutation
            if (not bad_mutation) {
                for (auto c : m.qual){
                    if (c-33 < min_qual){
                        bad_mutation = true;
                        break;
                    }
                }
            }
            // look at neighbor on left
            if (not bad_mutation){
                try {
                    int k = right_mut_indices.at(m.left + 1 - left_target_pos);
                    if ( k != -1 and k != i){
                        TaggedMutation neighbor(mutations[k]);
                        if (neighbor.seq.length() > 0){
                            // look at right-most basecall in mutation
                            if (neighbor.qual.back()-33 < min_qual){
                                bad_mutation = true;
                            }
                        } else {
                            // look at basecall on other side of gap
                            int n = neighbor.left - left_target_pos;
                            if (local_target_qual.at(n)-33 < min_qual){
                                bad_mutation = true;
                            }
                        }
                    } else {
                        // no mutation on left, just look at basecall
                        int n = m.left-left_target_pos;
                        if (local_target_qual.at(n)-33 < min_qual){
                            bad_mutation = true;
                        }
                    }
                } catch (std::out_of_range &exception) { }
            }
            // look at neighbor on right
            if (not bad_mutation){
                try{
                    int k = left_mut_indices.at(m.right - 1 - left_target_pos);
                    if (k != -1 and k != i){
                        TaggedMutation neighbor(mutations[k]);
                        if (neighbor.seq.length() > 0){
                            // look at left-most basecall in mutation
                            if (neighbor.qual.front()-33 < min_qual){
                                bad_mutation = true;
                            }
                        } else {
                            // look at basecall on other side of gap
                            int n = neighbor.right - left_target_pos;
                            if (local_target_qual.at(n)-33 < min_qual){
                                bad_mutation = true;
                            }
                        }
                    } else {
                        // no mutation on right, just look at basecall
                        int n = m.right-left_target_pos;
                        if (local_target_qual.at(n)-33 < min_qual){
                           bad_mutation = true;
                        }
                    }
                } catch (std::out_of_range &exception) { }
            }

            if (bad_mutation){
                // exclude entire covered region from depths
                for (int n = m.left + 1 - left_target_pos;
                     n < m.right - left_target_pos;
                     ++n) {
                    try {
                        local_effective_depth.at(n) = 0;
                    } catch (std::out_of_range &exception) { }
                }
                excluded_mutations.push_back(TaggedMutation(m));
            } else {
                if (variant_mode){
                    // include entire covered region (even in gap positions)
                    for (int n = m.left + 1 - left_target_pos;
                         n <= m.right - 1 - left_target_pos;
                         ++n) {
                        try {
                            local_effective_depth.at(n) = 1;
                        } catch (std::out_of_range &exception) { }
                    }
                } else {
                    // exclude covered region, include inferred adduct site
                    for (int n = m.left + 1 - left_target_pos;
                         n < m.right - 1 - left_target_pos;
                         ++n) {
                        try {
                            local_effective_depth.at(n) = 0;
                        } catch (std::out_of_range &exception) { }
                    }
                    try {
                        local_effective_depth.at(m.right-1-left_target_pos) = 1;
                    } catch (std::out_of_range &exception) { }
                }
                included_mutations.push_back(TaggedMutation(m));
            }
        }

        // also mark mutation locations in simplified array
        // (unused by downstream shapemapper components, but simpler to parse
        //  for single-molecule applications)
        std::vector<bool> local_effective_count(len, 0);
        for (auto &m:included_mutations) {
            try {
                local_effective_count.at(m.right - 1 - left_target_pos) = 1;
            } catch (std::out_of_range &exception) { }
        }

        return boost::make_tuple(local_effective_depth,
                                 local_effective_count,
                                 included_mutations,
                                 excluded_mutations);
    };


    /**
     * @brief Used in mutation_parser to write output to a stream.
     */
    std::string serializeReadInfo(const std::string &read_id,
                                  const int left,
                                  const int right,
                                  const std::string &local_target_seq,
                                  const std::string &local_target_qual,
                                  const std::vector<Mutation> &mutations) {
        using std::to_string;
        std::string o = read_id + " " +
                        to_string(left) + " " +
                        to_string(right) + " " +
                        local_target_seq + " " +
                        local_target_qual + " " +
                        toString(mutations) + "\n";
        return o;
    }

    /**
     * @brief Used in mutation_counter to read input from a stream.
     */
    boost::tuple<std::string, int, int, std::string, std::string, std::vector<Mutation>>
    parseReadInfo(const std::string &line) {
        using std::stoi;

        std::string read_id;
        int left;
        int right;
        std::string local_target_seq;
        std::string local_target_qual;
        std::vector<Mutation> mutations;

        std::vector<std::string> fields;
        std::string trimmed = boost::trim_copy(line);
        boost::split(fields, trimmed, boost::is_space(), boost::token_compress_on);

        if (fields.size() < 5) {
            throw std::runtime_error("Error: unable to parse incomplete line.");
        }
        try {
            read_id = fields[0];
            left = stoi(fields[1]);
            right = stoi(fields[2]);
        } catch (std::exception &err) {
            throw std::runtime_error("Error: line is incorrectly formatted (couldn't parse left or right position).");
        }
        local_target_seq = fields[3];
        local_target_qual = fields[4];

        if ((fields.size()-5) % 4 != 0) {
            throw std::runtime_error(
                    "Error: unable to parse mutation from incomplete line. " + std::to_string(fields.size()) +
                    " fields in line.");
        }
        for (int i = 5; i <= fields.size() - 4; i += 4) {
            try {
                Mutation m(stoi(fields[i]),
                           stoi(fields[i + 1]),
                           fields[i + 2].substr(1, fields[i + 2].length() - 2),
                           fields[i + 3].substr(1, fields[i + 2].length() - 2));
                mutations.push_back(m);
            } catch (std::exception &err) {
                throw std::runtime_error(
                        "Error: line is incorrectly formatted (couldn't parse mutation left or right bounds).");
            }
        }

        return boost::make_tuple(read_id, left, right, local_target_seq, local_target_qual, mutations);
    };


}

#endif //SHAPEMAPPER_MUTATION_H_H
