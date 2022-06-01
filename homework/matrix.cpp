//
// Created by simon on 01/06/2022.
//

#include <algorithm>
#include <limits>
#include "matrix.h"

namespace pc {
    matrix::matrix(int size) {
        create(size);
    }
    matrix::~matrix(){
        destroy();
    }

    void matrix::destroy() {
        for(int i =0 ; i < _size; i++) {
            delete[] _matrix[i];
        }
        delete[] _matrix;
    }

    void matrix::create(int size) {
        _size = size;
        _matrix = new double*[size];
        for(int i =0 ; i< size; i++) {
            _matrix[i] = new double[size];
        }
    }

    int  matrix::getSize() const {
        return _size;
    }

    double matrix::get(int i, int j) const {
        if(i < _size && j < _size)
            return _matrix[i][j];
    }

    void matrix::set(int i, int j, double w){
        if(i < _size && j < _size)
            _matrix[i][j] = w;
    }

    void matrix::changeSize(int newSize) {
        if(newSize != _size) {
            destroy();
            create(newSize);
        }
    }

    double *matrix::get(int i) {
        if(i<_size)
            return _matrix[i];
        return nullptr;
    }

    std::ostream& operator<<(std::ostream& os, const matrix &mt) {
        os << mt.getSize() << ";";
        for(int i =0; i<mt.getSize();i++) {
            for(int j = 0; j< mt.getSize();j++) {
                if(j == mt.getSize() -1) {
                    os << mt.get(i,j) << ";";
                } else {
                    os << mt.get(i,j) << ",";
                }
            }
        }
        return os;
    }

    std::istream &operator>>(std::istream &is, matrix& mt) {
        std::string line;
        std::getline(is, line, ';');
        int size = std::stoi(line);

        mt.changeSize(size);
        for(int i = 0; i<size; i++) {
            std::getline(is, line, ';');

            std::vector<std::string> result = split(line, ",");
            std::vector<double> doubleVector(result.size());
            std::transform(result.begin(), result.end(), doubleVector.begin(),
                      [](std::string const& val) {
                          if(val == "i") return std::numeric_limits<double>::infinity();
                          return std::stod(val);});
            std::copy(doubleVector.begin(), doubleVector.end(), mt.get(i));
        }


        return is;
    }

    std::vector<std::string> split (std::string s, std::string delimiter) {
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

} // pc