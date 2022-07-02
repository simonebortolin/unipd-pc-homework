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
        void (* _destroy)(T* &matrix);
    public:
        int getSize() const;
        T get(int, int) const;
        T* get(int);
        void set(int, int, T);
        bool operator==(basic_matrix<T> &);
        bool operator!=(basic_matrix<T> &);
        T* operator[](int);
        T* begin();
        T* end();
        void fill(T t);
        basic_matrix(basic_matrix<T>&&) noexcept;
        basic_matrix<T>& operator=(basic_matrix<T>&&);
    protected:
        basic_matrix(T*, int, void (*)(T*&));
        ~basic_matrix();
    };

    template <class T>
    std::ostream &operator<<(std::ostream &os, const basic_matrix<T>& mt);

    template <class T>
    std::istream &operator>>(std::istream &is, basic_matrix<T>& mt);

    template <class T>
    bool basic_matrix<T>::operator==(basic_matrix<T> &b) {
        if(b._size != _size) return false;
        for(int i =0; i<_size; i++) {
            for(int j =0; j<_size; j++) {
                if(b.get(i,j) != get(i,j)) {
                    return false;
                }
            }
        }
        return true;
    }

    template <class T>
    bool basic_matrix<T>::operator!=(basic_matrix<T> &b) {
        return !((*this) == b);
    }

    template <class T>
    T convertTo (const std::string &str);

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

    template <class T>
    std::ostream& operator<<(std::ostream& os, const basic_matrix<T>& mt) {
        std::stringstream ss;
        ss << mt.getSize() << ";";
        for(int i =0; i<mt.getSize();i++) {
            for(int j = 0; j< mt.getSize();j++) {
                if(j == mt.getSize() -1) {
                    ss << mt.get(i,j) << ";";
                } else {
                    ss << mt.get(i,j) << ",";
                }
            }
        }
        os << ss.str();
        return os;
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
    basic_matrix<T>::basic_matrix(T * matrix, int size, void (* destroy)(T* &matrix)) {
        _matrix = matrix;
        _size = size;
        _destroy = destroy;
    }

    template<class T>
    basic_matrix<T>::~basic_matrix() {
        if(_destroy != nullptr)  _destroy(_matrix);
    }

    template<class T>
    basic_matrix<T>::basic_matrix(basic_matrix<T> &&m)  noexcept {
        _size = m._size;
        _matrix = m._matrix;
        _destroy = m._destroy;
        m._size = 0;
        m._matrix = nullptr;
        m._destroy = nullptr;
    }

    template<class T>
    basic_matrix<T> &basic_matrix<T>::operator=(basic_matrix<T> &&m) {
        if(&m != this) {
            if(_destroy != nullptr) _destroy(_matrix);
            _matrix = m._matrix;
            _size = m._size;
            _destroy = m._destroy;
            m._matrix = nullptr;
            m._size = 0;
            m._destroy = nullptr;
        }
        return *this;
    }
} // pc


#endif //HOMEWORK_BASIC_MATRIX_H
