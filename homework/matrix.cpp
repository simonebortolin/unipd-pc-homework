//
// Created by simon on 01/06/2022.
//

#include <algorithm>
#include <limits>
#include <sstream>
#include "matrix.h"

namespace pc {


    std::vector<std::string> split (const std::string& s, const std::string& delimiter) {
        size_t pos_start = 0, pos_end, delim_len = delimiter.length();
        std::string token;
        std::vector<std::string> res;

        while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos) {
            token = s.substr (pos_start, pos_end - pos_start);
            pos_start = pos_end + delim_len;
            res.push_back (token);
        }

        res.push_back (s.substr (pos_start));
        return res;
    }

    template <>
    double convertTo<double>(const std::string &str) {
        if(str == "i") return std::numeric_limits<double>::infinity();
        if(str == "inf") return std::numeric_limits<double>::infinity();
        return std::stod(str);
    }

    template <>
    int convertTo<int>(const std::string &str) {
        return std::stoi(str);
    }

} // pc