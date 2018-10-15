/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2018 Steven Busan.                     *
 *----------------------------------------------------------------------*/


#ifndef SHAPEMAPPER_MUTATION_H
#define SHAPEMAPPER_MUTATION_H

#include <iostream>
#include <fstream>
#include <exception>


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


/**
* @brief Stores deviations from an alignment target between two nucleotide
*        positions.
*/
class Mutation {
public:

    /** @brief Leftmost unchanged alignment target nucleotide (0-based) */
    int left;
    /** @brief Rightmost unchanged alignment target nucleotide (0-based) */
    int right;
    /** @brief Read sequence replacing alignment target sequence between left and right (exclusive) */
    std::string seq;
    /** @brief Basecall quality scores (ASCII encoded Phred scores) for nucs in seq */
    std::string qual;
    /** @brief Mutation tag (usually mutation classification)*/
    std::string tag;
    /** @brief Whether a mutation is or is derived from an ambiguous alignment */
    bool ambig;

    Mutation() {
        left = -999;
        right = -999;
        seq = "";
        qual = "";
        tag = "";
        ambig = false;
    };

    Mutation(const int l,
             const int r,
             const std::string &s,
             const std::string &q) {
        left = l;
        right = r;
        seq = s;
        qual = q;
        tag = "";
        ambig = false;
    }

    Mutation(const int l,
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

    Mutation(const int l,
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

    Mutation(const Mutation &m) {
        left = m.left;
        right = m.right;
        seq = m.seq;
        qual = m.qual;
        tag = m.tag;
        ambig = m.ambig;
    }

    Mutation(const Mutation &m,
             const std::string &t) {
        left = m.left;
        right = m.right;
        seq = m.seq;
        qual = m.qual;
        tag = t;
        ambig = false;
    }

    Mutation(const Mutation &m,
             const std::string &t,
             const bool a) {
        left = m.left;
        right = m.right;
        seq = m.seq;
        qual = m.qual;
        tag = t;
        ambig = a;
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


    std::string toString() const {
        using std::to_string;
        std::string s = to_string(left) + ' ' + to_string(right) + " \"";
        s += seq + "\" \"" + qual + "\" \"" + tag;
        if (ambig) {
            s += "_ambig";
        }
        s += '"';
        return s;
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
     *        Assumes identifyAmbiguousMutations() was previously run.
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
                if (seq == "N") {
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

bool operator==(const Mutation &a, const Mutation &b) {
    return (a.left == b.left and a.right == b.right and a.seq == b.seq);
}

// FIXME: need namespace to resolve function signature for some reason. otherwise
//        gets confused with DebugRead.toString()
namespace mutation {
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
}

int minLeft(const std::vector<Mutation> &vm) {
    if (vm.size() == 0) { return -9999; }
    int min_left = vm[0].left;
    for (int i = 1; i < vm.size(); i++) {
        min_left = std::min(vm[i].left, min_left);
    }
    return min_left;
}

int maxRight(const std::vector<Mutation> &vm) {
    if (vm.size() == 0) { return -9999; }
    int max_right = vm[0].right;
    for (int i = 1; i < vm.size(); i++) {
        max_right = std::max(vm[i].right, max_right);
    }
    return max_right;
}


class MutationGroup {
public:
    int left;
    int right;
    std::vector<Mutation> r1_mutations;
    std::vector<Mutation> r2_mutations;

    MutationGroup() { };
};

struct FieldsSizeException : public std::exception
{
	const char * what () const throw ()
    {
    	return "ERROR";
    }
};

std::vector<Mutation>
fieldsToMutationVect(const std::vector<std::string> &fields,
                     int start_index=0) {
    std::vector<Mutation> mutations;
    if (fields.size() % 5 != 0) {
        /*throw std::runtime_error(
            "Error: unable to read mutations from incomplete line. " + std::to_string(fields.size()) +
            " field(s) in right-most column.");*/
        throw FieldsSizeException();
    }
    for (int i = start_index; i <= fields.size() - 5; i += 5) {
        try {
            Mutation m(stoi(fields[i]),
                             stoi(fields[i + 1]),
                             fields[i + 2].substr(1, fields[i + 2].length() - 2),
                             fields[i + 3].substr(1, fields[i + 3].length() - 2),
                             fields[i + 4].substr(1, fields[i + 4].length() - 2)
            );
            mutations.push_back(m);
        } catch (std::exception &err) {
            throw std::runtime_error(
                    "Error: line is incorrectly formatted (couldn't read mutation left or right bounds).");
        }
    }
    return mutations;
}

std::vector<Mutation>
stringToMutationVect(const std::string &s){
    std::vector<std::string> fields;
    std::string trimmed = boost::trim_copy(s);
    boost::split(fields, trimmed, boost::is_space(), boost::token_compress_off);
    std::vector<Mutation> mutations;
    if (s.size() > 0){
        try {
            mutations = fieldsToMutationVect(fields);
            return mutations;
        } catch (FieldsSizeException &err) {
            throw std::runtime_error(
                "Error: unable to read mutations from incomplete line. " + std::to_string(fields.size()) +
                " field(s) in right-most column. Right-most column: \"" + s + "\"");
        }
    }
    return mutations;
}


#endif //SHAPEMAPPER_MUTATION_H
