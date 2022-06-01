#include <iostream>
#include "matrix.h"

pc::matrix floydWarhsallSeq(pc::matrix &dist);

int main() {
    pc::matrix c = pc::matrix(0);
    std::cin >> c;

    pc::matrix out = floydWarhsallSeq(c);
    std::cout << c << std::endl;
    std::cout << out << std::endl;

    return 0;

}


/*
 * Input
 * 5;0,3,8,i,-4;i,0,i,1,7;i,4,0,i,i;2,i,-5,0,i;i,i,i,6,0;
 * Output
 * 5;0,1,-3,2,-4;3,0,-4,1,-1;7,4,0,5,3;2,-1,-5,0,-2;8,5,1,6,0;
 * 5;0,2,3,4,0;3,1,3,1,0;3,2,2,1,0;3,2,3,3,0;3,2,3,4,4;
 */

pc::matrix floydWarhsallSeq(pc::matrix &dist){

    double tr, pr;
    int d = dist.getSize();
    pc::matrix pred = pc::matrix(d);

    for(int i=0; i<d;i++){
        for(int j=0; j<d;j++){
            pred.set(i,j,i);
        }
    }

    for(int h=0; h<d;h++){
        for(int i=0; i<d;i++){
            for(int j=0; j<d;j++){
                tr = dist.get(i,h) + dist.get(h,j);
                if(  tr < dist.get(i,j) ){
                    dist.set(i,j,tr);
                    pr = pred.get(h,j);
                    pred.set(i,j,pr);
                }
            }
        }
        for(int i=0;i<d;i++ ){
            if(dist.get(i,i) < 0) {
                std:: cout<< "error";
                return pred ; // TODO
            }
        }
    }
    return pred;
}
