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
            T ** _matrix;
            int _size;
            void create(int size);
            void destroy();
        public:
            explicit matrix(int size);
            int getSize() const;
            T get(int i, int j) const;
            T* get(int i);
            void set(int i, int j, T w);
            void changeSize(int newSize);
            bool operator==(matrix<T> &);
            bool operator!=(matrix<T> &);

        ~matrix();
    };

    template <class T>
    bool matrix<T>::operator==(matrix<T> &b) {
        if(b._size != _size) return false;
        for(int i =0; i<_size; i++) {
            for(int j =0; j<_size; j++) {
                if(b.get(i,j) != get(i,j)) {
                    std::cout << i << " " << j << std::endl;
                    //return false;
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
    std::ostream& operator<<(std::ostream& os, const matrix<T>& mt);

    template <class T>
    std::istream &operator>>(std::istream &is, matrix<T>& mt);

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
        for(int i =0 ; i < _size; i++) {
            delete[] _matrix[i];
        }
        delete[] _matrix;
        _matrix = nullptr;
    }

    template <class T>
    void matrix<T>::create(int size) {
        _size = size;
        _matrix = new T*[size];
        for(int i =0 ; i< size; i++) {
            _matrix[i] = new T[size];
        }
    }

    template <class T>
    int matrix<T>::getSize() const {
        return _size;
    }

    template <class T>
    T matrix<T>::get(int i, int j) const {
        if(i < _size && j < _size)
            return _matrix[i][j];
    }

    template <typename T>
    void matrix<T>::set(int i, int j, T w){
        if(i < _size && j < _size)
            _matrix[i][j] = w;
    }

    template <class T>
    void matrix<T>::changeSize(int newSize) {
        if(newSize != _size) {
            destroy();
            create(newSize);
        }
    }

    template <class T>
    T *matrix<T>::get(int i) {
        if(i<_size)
            return _matrix[i];
        return nullptr;
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


} // pc

#endif //HOMEWORK_MATRIX_H
