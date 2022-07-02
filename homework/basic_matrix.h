//
// Created by simone on 7/2/22.
//

#ifndef HOMEWORK_BASIC_MATRIX_H
#define HOMEWORK_BASIC_MATRIX_H

#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include <cstring>
#include <sstream>

namespace pc {
    template <class T>
    class basic_matrix {
    protected:
        T * _matrix;
        int _size;
    public:
        int getSize() const;
        T get(int, int) const;
        T* get(int);
        void set(int, int, T);
        T* operator[](int);
        T* begin();
        T* end();
        void fill(T t);
        basic_matrix(T*, int);
        basic_matrix();
    };

    template <class T>
    void basic_matrix<T>::fill(T t) {

        std::fill (_matrix, _matrix + _size*_size, t);
    }

    template <class T>
    int basic_matrix<T>::getSize() const {
        return _size;
    }

    template <class T>
    T basic_matrix<T>::get(int i, int j) const {
        if(i < _size && j < _size && i >= 0 && j >= 0)
            return _matrix[i*_size + j];
        throw std::invalid_argument( "invalid i or j" );
    }

    template <class T>
    void basic_matrix<T>::set(int i, int j, T w){
        if(i < _size && j < _size && i >= 0 && j >= 0)
            _matrix[i*_size + j] = w;
        else {
            throw std::invalid_argument( "invalid i or j" );
        }
    }

    template <class T>
    T* basic_matrix<T>::get(int i) {
        if(i<_size && i >= 0) {
            return &_matrix[i*_size];
        }
        return nullptr;
    }

    template<class T>
    T *basic_matrix<T>::operator[](int i) {
        return get(i);
    }

    template<class T>
    T *basic_matrix<T>::begin() {
        return _matrix;
    }

    template<class T>
    T *basic_matrix<T>::end() {
        return _matrix + _size * _size;
    }

    template<class T>
    basic_matrix<T>::basic_matrix(T * matrix, int size) {
        _matrix = matrix;
        _size = size;
    }

    template<class T>
    basic_matrix<T>::basic_matrix() {
        _size = 0;
        _matrix = nullptr;
    }
} // pc


#endif //HOMEWORK_BASIC_MATRIX_H
