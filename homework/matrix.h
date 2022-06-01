//
// Created by simon on 01/06/2022.
//

#ifndef HOMEWORK_MATRIX_H
#define HOMEWORK_MATRIX_H

#include <iostream>
#include <vector>

namespace pc {

    class matrix {
        private:
            double ** _matrix;
            int _size;
            void create(int size);
            void destroy();
        public:
            explicit matrix(int size);
            int getSize() const;
            double get(int i, int j) const;
            double* get(int i);
            void set(int i, int j, double w);
            void changeSize(int newSize);



        ~matrix();
    };

    std::ostream& operator<<(std::ostream& os, const matrix& mt);
    std::istream &operator>>(std::istream &is, matrix& mt);


    std::vector<std::string> split(std::string s, std::string delimiter);

} // pc

#endif //HOMEWORK_MATRIX_H
