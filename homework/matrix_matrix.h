//
// Created by simone on 7/2/22.
//

#ifndef HOMEWORK_MATRIX_MATRIX_H
#define HOMEWORK_MATRIX_MATRIX_H


#include "basic_matrix.h"

namespace pc {
    template<class T>
    class matrix_matrix : public basic_matrix<basic_matrix<T>> {
        T* _allocator;
        int _step;
    public:
        matrix_matrix(int size, int step);
        matrix_matrix(matrix_matrix<T> &&);
        matrix_matrix<T> operator=(matrix_matrix<T> &&);
        ~matrix_matrix();

    };

    template<class T>
    matrix_matrix<T>::matrix_matrix(int size, int step) : basic_matrix<basic_matrix<T>>(new basic_matrix<T>[size * size], size) {
       // size_t size_of_buffer = size * size;
	//std::cerr << size_of_buffer * sizeof(T) /1024 /1024 << "MB" << std::endl;
	    _allocator = new T[size*size * step*step];
        _step = step;
        for(int i = 0; i < size; i++) {
            for(int j = 0; j < size; j++) {
                this->set(i,j, basic_matrix<T>(_allocator + (i * step * step * size + j * step * step), _step));
            }
        }
    }

    template<class T>
    matrix_matrix<T>::~matrix_matrix() {
        delete[] _allocator;
        delete[] this->_matrix;
        _allocator = nullptr;
        this->_matrix = nullptr;
    }

    template<class T>
    matrix_matrix<T>::matrix_matrix(matrix_matrix<T> &&m) {
        this->_matrix = m._matrix;
        this->_size = m._size;
        _allocator = m._allocator;
        _step = m._step;
        m._matrix = nullptr;
        m._size = 0;
        m._allocator = nullptr;
        m._step  = 0;
    }

    template<class T>
    matrix_matrix<T> matrix_matrix<T>::operator=(matrix_matrix<T> &&m) {
        if(&m != this) {
            delete[] _allocator;
            delete[] this->_matrix;
            this->_matrix = m._matrix;
            this->_size = m._size;
            _allocator = m._allocator;
            _step = m._step;
            m._matrix = nullptr;
            m._size = 0;
            m._allocator = nullptr;
            m._step  = 0;
        }
        return *this;
    }


}


#endif //HOMEWORK_MATRIX_MATRIX_H
