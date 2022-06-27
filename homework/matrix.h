//
// Created by simon on 01/06/2022.
//

#ifndef HOMEWORK_MATRIX_H
#define HOMEWORK_MATRIX_H

#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>

namespace pc {
    template <class T>
    class matrix {
    private:
        T * _matrix;
        int _size;
        void create(int);
        void destroy();
    public:
        matrix(const matrix&);
        matrix(matrix&&);
        explicit matrix(int);
        int getSize() const;
        T get(int, int ) const;
        T* get(int);
        void set(int, int, T );
        void changeSize(int);
        bool operator==(matrix<T> &);
        bool operator!=(matrix<T> &);
        matrix<T>& operator=(matrix&&);
        T* operator[](int);

        ~matrix();

    };

    template <class T>
    std::ostream &operator<<(std::ostream &os, const matrix<T>& mt);

    template <class T>
    std::istream &operator>>(std::istream &is, matrix<T>& mt);

    template <class T>
    bool matrix<T>::operator==(matrix<T> &b) {
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
    bool matrix<T>::operator!=(matrix<T> &b) {
        return !((*this) == b);
    }

    std::vector<std::string> split(const std::string& s, const std::string& delimiter);

    template <class T>
    T convertTo (const std::string &str);

    template <class T>
    matrix<T>::matrix(int size) {
        _size = size;
        create(size);
    }

    template <class T>
    matrix<T>::~matrix(){
        destroy();
    }

    template <class T>
    void matrix<T>::destroy() {
        delete[] _matrix;
        _matrix = nullptr;
    }

    template <class T>
    void matrix<T>::create(int size) {
        _size = size;
        _matrix = new T[size*size];
    }

    template <class T>
    int matrix<T>::getSize() const {
        return _size;
    }

    template <class T>
    T matrix<T>::get(int i, int j) const {
        if(i < _size && j < _size && i >= 0 && j >= 0)
            return _matrix[i*_size + j];

        throw std::invalid_argument( "invalid i or j" );
    }

    template <typename T>
    void matrix<T>::set(int i, int j, T w){
        if(i < _size && j < _size && i >= 0 && j >= 0)
            _matrix[i*_size + j] = w;
        else throw std::invalid_argument( "invalid i or j" );
    }

    template <class T>
    void matrix<T>::changeSize(int newSize) {
        if(newSize != _size) {
            destroy();
            create(newSize);
        }
    }

    template <class T>
    T* matrix<T>::get(int i) {
        if(i<_size && i >= 0) {
            return &_matrix[i*_size];
        }
        return nullptr;
    }

    template<class T>
    matrix<T>::matrix(const matrix &m) {
        _size =m._size;
        _matrix = new T[_size];
        for(int i =0; i<m._size; i++) {
            for(int j =0; j<m._size; j++) {
                set(i,j,m._matrix[i*_size + j]);
            }
        }
    }

    template<class T>
    matrix<T>::matrix(matrix &&m) {
        _size = m._size;
        _matrix = m._matrix;
        m._size = 0;
        m._matrix = nullptr;
    }

    template <class T>
    std::ostream& operator<<(std::ostream& os, const matrix<T>& mt) {
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

    template <class T>
    std::istream &operator>>(std::istream &is, matrix<T>& mt) {

        std::string line;
        std::getline(is, line, ';');
        int size = std::stoi(line);

        mt.changeSize(size);
        for(int i = 0; i<size; i++) {
            std::getline(is, line, ';');

            std::vector<std::string> result = split(line, ",");
            std::vector<double> doubleVector(result.size());
            std::transform(result.begin(), result.end(), doubleVector.begin(),
                           [](std::string const& val) { return convertTo<T>(val); });
            std::copy(doubleVector.begin(), doubleVector.end(), mt.get(i));
        }


        return is;
    }

    template<class T>
    matrix<T> &matrix<T>::operator=(matrix &&m) {
        if(&m != this) {
            delete[] _matrix;
            _matrix = m._matrix;
            _size = m._size;
            m._matrix = nullptr;
            m._size = 0;
        }
        //return m;
    }

    template<class T>
    T *matrix<T>::operator[](int i) {
        return get(i);
    }


} // pc

#endif //HOMEWORK_MATRIX_H