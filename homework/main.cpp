#include <iostream>
#include "matrix.h"
#include "sequential.h"

int main(int argc , char ** argv) {
    int s;
    std::cin >> s;

    pc::matrix<double> cseq = pc::matrix<double>(0);
    std::cin >> cseq;
    pc::matrix<double> csq = cseq;

    pc::matrix<int> out = floydWarhsallSeq(cseq);
    pc::matrix<int> outc = floydWarhsallSquared(csq,s);

    if(cseq != csq) {
        std::cout << "not ok dist" << std::endl;
    } else {
        std::cout << csq << std::endl; //distances matrix: min costs
        std::cout << outc << std::endl; //predecessors matrix: min path
    }

    return 0;

}


