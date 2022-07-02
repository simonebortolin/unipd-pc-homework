//
// Created by simon on 01/06/2022.
//

#ifndef HOMEWORK_MATRIX_H
#define HOMEWORK_MATRIX_H

#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include <cstring>
#include <sstream>
#include "basic_matrix.h"

namespace pc {
    template <class T>
    class matrix : public basic_matrix<T> {
    private:
        void destroy();
    public:
        // constructor
        explicit matrix(int);
        matrix();

        // copy and move constructor
        matrix(const matrix<T>&);
        matrix(matrix<T> &&m);

        // copy and move operator
        matrix<T>& operator=(const matrix<T>&);
        matrix<T> &operator=(matrix<T> &&m);

        // operator == and !=
        bool operator==(matrix<T> &);
        bool operator!=(matrix<T> &);

        // destructor
        ~matrix();
    };

    template <class T>
    std::ostream &operator<<(std::ostream &os, const matrix<T>& m);

    template <class T>
    std::istream &operator>>(std::istream &is, matrix<T>& mt);

    std::vector<std::string> split(const std::string& s, const std::string& delimiter);

    template <class T>
    bool matrix<T>::operator==(matrix<T> &b) {
        if(b._size != this->_size) return false;
        for(int i =0; i<this->_size; i++) {
            for(int j =0; j<this->_size; j++) {
                if(b.get(i,j) != this->get(i,j)) {
                    return false;
                }
            }
        }
        return true;
    }

    template <class T>
    bool matrix<T>::operator!=(matrix<T> &b) {
        return !((*this) == b);
    }

    template <class T>
    T convertTo (const std::string &str);

    template <class T>
    void matrix<T>::destroy() {
        delete[] this->_matrix;
        this->_matrix = nullptr;
    }

    template <class T>
    matrix<T>::matrix(int size) :  basic_matrix<T>(new T[size*size], size) { }

    template <class T>
    matrix<T>::matrix(const matrix<T> &m) : basic_matrix<T>(new T[m._size*m._size], m._size) {
        for(int i =0; i<m._size; i++) {
            for(int j =0; j<m._size; j++) {
                this->set(i,j,m._matrix[i*this->getSize() + j]);
            }
        }
    }

    template <class T>
    std::istream &operator>>(std::istream &is, matrix<T>& m) {

        std::string line;
        std::getline(is, line, ';');
        int size = std::stoi(line);

        if(m.getSize() != size)
            m = matrix<T>(size);

        for(int i = 0; i<size; i++) {
            std::getline(is, line, ';');

            std::vector<std::string> result = split(line, ",");
            std::vector<double> doubleVector(result.size());
            std::transform(result.begin(), result.end(), doubleVector.begin(),
                           [](std::string const& val) { return convertTo<T>(val); });
            std::copy(doubleVector.begin(), doubleVector.end(), m.get(i));
        }


        return is;
    }

    template <class T>
    std::ostream& operator<<(std::ostream& os, const matrix<T>& m) {
        std::stringstream ss;
        ss << m.getSize() << ";";
        for(int i =0; i < m.getSize(); i++) {
            for(int j = 0; j < m.getSize(); j++) {
                if(j == m.getSize() - 1) {
                    ss << m.get(i, j) << ";";
                } else {
                    ss << m.get(i, j) << ",";
                }
            }
        }
        os << ss.str();
        return os;
    }

    template<class T>
    matrix<T> &matrix<T>::operator=(const matrix<T> &m) {
        if(&m != this) {
            if(this->getSize() != m._size) {
                destroy();
                (*this) = matrix<T>(m._size);
            }
            std::copy(m._matrix, m._matrix + (m._size * m._size), this->_matrix);
        }
        return *this;
    }

    template<class T>
    matrix<T>::matrix() : basic_matrix<T>(nullptr, 0, nullptr) { }

    template<class T>
    matrix<T>::matrix(matrix<T> &&m) : basic_matrix<T>(m._matrix, m._size) {
        m._size = 0;
        m._matrix = nullptr;
    }

    template<class T>
    matrix<T> &matrix<T>::operator=(matrix<T> &&m) {
        if(&m != this) {
            destroy();
            this->_matrix = m._matrix;
            this->_size = m._size;
            m._matrix = nullptr;
            m._size = 0;
        }
        return *this;
    }

    template<class T>
    matrix<T>::~matrix() {
        destroy();
    }
} // pc

#endif //HOMEWORK_MATRIX_H