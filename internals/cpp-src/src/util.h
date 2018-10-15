

/*-----------------------------------------------------------------------
 * This file is a part of ShapeMapper, and is licensed under the terms  *
 * of the MIT license. Copyright 2018 Steven Busan.                     *
 *----------------------------------------------------------------------*/

#ifndef SHAPEMAPPER_UTIL_H
#define SHAPEMAPPER_UTIL_H

#include <iostream>
#include <fstream>


namespace util {
    std::string toString(const std::vector<bool> &bools) {
        std::string s;
        for (auto b : bools) {
            s += std::to_string(b);
        }
        return s;
    }


    bool endsWith(const std::string &s,
                  const std::string &end) {
        if (end.size() > s.size()) {
            return false;
        } else {
            return std::equal(end.rbegin(), end.rend(), s.rbegin());
        }
    }


    int indexOf(const std::vector<std::string> &vect,
                const std::string &s) {
        for (int i = 0; i < vect.size(); i++) {
            if (vect[i] == s) {
                return i;
            }
        }
        throw std::runtime_error("Error: string not found in vector");
    }

    std::vector<bool>
    stringToBoolVect(const std::string &s) {
        std::vector<bool> vect;
        bool b;
        for (auto c : s) {
            b = false;
            if (c == '1') { b = true; }
            vect.push_back(b);
        }
        return vect;
    }
}

#endif //SHAPEMAPPER_UTIL_H
