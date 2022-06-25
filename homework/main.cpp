#include <iostream>
#include "matrix.h"

pc::matrix<int> floydWarhsallSeq(pc::matrix<double> &dist);
pc::matrix<int> floydWarhsallSquared(pc::matrix<double> &dist, int s);

int main() {
    int s;
    std::cin >> s;

    pc::matrix<double> c = pc::matrix<double>(0);
    std::cin >> c;

    pc::matrix<int> out = floydWarhsallSquared(c,s);
    std::cout << c << std::endl; //distances matrix: min costs
    std::cout << out << std::endl; //predecessors matrix: min path



    return 0;

}


/*
 * Input
 * 5;0,3,8,i,-4;i,0,i,1,7;i,4,0,i,i;2,i,-5,0,i;i,i,i,6,0;
 * Output
 * 5;0,1,-3,2,-4;3,0,-4,1,-1;7,4,0,5,3;2,-1,-5,0,-2;8,5,1,6,0;
 * 5;0,2,3,4,0;3,1,3,1,0;3,2,2,1,0;3,2,3,3,0;3,2,3,4,4;
 */

/*
 * Input
 * 6;0,3,8,i,-4,i;i,0,i,1,7,i;i,4,0,i,i,2;2,i,-5,0,i,i;i,i,i,6,0,i;i,4,i,i,i,0;
 * Output
 * 6;0,1,-3,2,-4,-1;3,0,-4,1,-1,-2;7,4,0,5,3,2;2,-1,-5,0,-2,-3;8,5,1,6,0,3;7,4,0,5,3,0; <---
 * 6;0,2,3,4,0,2;3,1,3,1,0,2;3,2,2,1,0,2;3,2,3,3,0,2;3,2,3,4,4,2;3,5,3,1,0,5;
 */

/*
 * Input
 * 10;0,3,4,i,i,i,i,i,i,i;i,0,2,i,7,i,i,i,i,i;i,i,0,3,i,i,i,i,i,i;i,i,i,0,i,i,i,i,5,i;i,i,i,2,0,5,i,i,i,i;i,i,i,i,i,0,1,i,i,4;i,i,i,i,i,i,0,2,i,i;i,i,i,i,i,i,i,0,i,8;i,i,i,i,i,i,2,i,0,i;2,i,i,i,i,i,i,i,i,0;
 * Output
 * 10;0,3,4,7,10,15,14,16,12,19;18,0,2,5,7,12,12,14,10,16;22,25,0,3,32,37,10,12,8,20;19,22,23,0,29,34,7,9,5,17;11,14,15,2,0,5,6,8,7,9;6,9,10,13,16,0,1,3,18,4;12,15,16,19,22,27,0,2,24,10;10,13,14,17,20,25,24,0,22,8;14,17,18,21,24,29,2,4,0,12;2,5,6,9,12,17,16,18,14,0;
 * 10;0,0,0,2,1,4,8,6,3,5;9,1,1,2,1,4,8,6,3,5;9,0,2,2,1,4,8,6,3,7;9,0,0,3,1,4,8,6,3,7;9,0,0,4,4,4,5,6,3,5;9,0,0,2,1,5,5,6,3,5;9,0,0,2,1,4,6,6,3,7;9,0,0,2,1,4,8,7,3,7;9,0,0,2,1,4,8,6,8,7;9,0,0,2,1,4,8,6,3,9;
 */

pc::matrix<int> floydWarhsallSeq(pc::matrix<double> &dist){

    double tr;
    int pr;
    int d = dist.getSize();
    pc::matrix<int> pred = pc::matrix<int>(d);

    // initialised predecessors
    for(int i=0; i<d;i++){
        for(int j=0; j<d;j++){
            pred.set(i,j,i);
        }
    }

    //computation of min path costs
    for(int h=0; h<d;h++){
        for(int i=0; i<d;i++){
            for(int j=0; j<d;j++){

                //path comparison between current i->j and i->h->j
                tr = dist.get(i,h) + dist.get(h,j);

                if(  tr < dist.get(i,j) ){

                    //update of min path distance and predecessors
                    dist.set(i,j,tr);
                    pr = pred.get(h,j);
                    pred.set(i,j,pr);

                }
            }
        }

        //possible negative cycle
        for(int i=0;i<d;i++ ){
            if(dist.get(i,i) < 0) {
                std:: cout<< "error";
                return pred ; // TODO
            }
        }

        std::cout << dist << std::endl;
    }
    return pred;
}

int pippo = 0;
// s: submatrix size
// dx,dy: start index of submatrix d
// ax,ay: start index of submatrix a
// bx,by: start index of submatrix b

//0,2,0
void floyd(pc::matrix<double> &dist, pc::matrix<int> &pred, int dx, int dy, int ax, int ay, int bx, int by, int s) {
    double tr;
    int pr;
    pippo++;
    for(int h=0; h<s;h++){
        for(int i=0; i<s;i++){
            for(int j=0; j<s;j++){

                //path comparison between current i->j and i->h->j
                tr = dist.get(i+ax,h+ay) + dist.get(h+bx,j+by);

                if(  tr < dist.get(i+dx,j+dy) ){

                    //update of min path distance and predecessors
                    dist.set(i+dx,j+dy,tr);
                    pr = pred.get(h+bx,j+by);
                    pred.set(i+dx,j+dy,pr);

                }
            }
        }

    }
}

pc::matrix<int> floydWarhsallSquared(pc::matrix<double> &dist, int s){


    double tr;
    int pr;
    int n = dist.getSize();

    if(((int)((double)n/(double)s))%2 == 0) {
        std:: cout<< n/s;
        std:: cout<< "work only with a odd division of matrix";
        return pc::matrix<int>(0);
    }

    pc::matrix<int> pred = pc::matrix<int>(n);


    // initialised predecessors
    for(int i=0; i<n;i++){
        for(int j=0; j<n;j++){
            pred.set(i,j,i);
        }
    }

    //std::cout << dist << std::endl;

    //computation of min path costs
    for(int h=0; h<n;h+=s){
        //verde scuro self dependent
        floyd(dist, pred, h,h,h,h,h,h,s); // first


        for(int k=1; k<=(n/s)/2;k++){ // second - ....
            std:: cout << "h: " << h <<  " k: " << k << std::endl;
            //verde chiaro - riga h
            int j = (h-k*s+n) % n;
            floyd(dist,pred,h,j,h,h,h,j,s);

            j = (h+k*s+n) % n;
            floyd(dist,pred,h,j,h,h,h,j,s);

            // verde chiaro - colonna h
            int i = (h-k*s+n) % n;
            floyd(dist, pred, i, h, i, h, h, h,s);

            i = (h+k*s+n) % n;
            floyd(dist, pred, i, h, i, h, h, h,s);

            // bianco - resto

            for(int l = k-(n/s)/2; l<=(n/s)/2+k; l++) {
                int m = (l+k)*s;
                i = (m+n) % n;
                if(i == h) continue;
                j = (h-k*s+n) % n;
                floyd(dist,pred,i,j,i,h,h,j,s);
                j = (h+k*s+n) % n;
                floyd(dist,pred,i,j,i,h,h,j,s);
            }
            for(int l = k-(n/s)/2+s; l<=(n/s)/2+k-s; l++) {
                int m = (l+k)*s;
                i = (m+n) % n;
                if(i == h) continue;
                j = (h-k*s+n) % n;
                floyd(dist,pred,i,j,i,h,h,j,s);
                j = (h+k*s+n) % n;
                floyd(dist,pred,i,j,i,h,h,j,s);
            }
            std::cout << dist << std::endl;


        }


        //possible negative cycle
        for(int i=0;i<n;i++ ){
            if(dist.get(i,i) < 0) {
                std:: cout<< "error";
                return pred ; // TODO
            }
        }
    }
    return pred;
}


