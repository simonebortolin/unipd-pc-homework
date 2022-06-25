#include <iostream>
#include "matrix.h"

pc::matrix<int> floydWarhsallSeq(pc::matrix<double> &dist);
pc::matrix<int> floydWarhsallSquared(pc::matrix<double> &dist, int s);

int main() {
    int s;
    std::cin >> s;

    pc::matrix<double> c = pc::matrix<double>(0);
    std::cin >> c;
    pc::matrix<double> cc = c;

    pc::matrix<int> out = floydWarhsallSeq(c);
    pc::matrix<int> outc = floydWarhsallSquared(cc,s);

    if(c != cc) {
        std::cout << "not ok dist" << std::endl;
    }
    else if(out != outc)  {
        std::cout << "not ok pred" << std::endl;
    } else {

        std::cout << c << std::endl; //distances matrix: min costs
        std::cout << out << std::endl; //predecessors matrix: min path
    }


    return 0;

}

<<<<<<< HEAD
=======

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

>>>>>>> 0d280cf (revert)
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

        //std::cout << dist << std::endl;
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


<<<<<<< HEAD
        for(int k=1; k<=(n/s)/2;k++){ // second - ....
=======
        for(int k=0; k<(n/s)/2;k++){ // second - ....
            std:: cout << "h: " << h <<  " k: " << k << std::endl;
>>>>>>> 0d280cf (revert)
            //verde chiaro - riga h
            int j = (h-s+n) % n;
            floyd(dist,pred,h,j,h,h,h,j,s);

            j = (h+s+n) % n;
            floyd(dist,pred,h,j,h,h,h,j,s);

            // verde chiaro - colonna h
            int i = (h-s+n) % n;
            floyd(dist, pred, i, h, i, h, h, h,s);

            i = (h+s+n) % n;
            floyd(dist, pred, i, h, i, h, h, h,s);

            // bianco - resto

            for(int l = k-(n/s)/2; l<=(n/s)/2+k; l++) {
                int m = (l+k)*s;
                i = (m+n) % n;
                if(i == h) continue;
                j = (h-s+n) % n;
                floyd(dist,pred,i,j,i,h,h,j,s);
                j = (h+s+n) % n;
                floyd(dist,pred,i,j,i,h,h,j,s);
            }
            for(int l = k-(n/s)/2+s; l<=(n/s)/2+k-s; l++) {
                int m = (l+k)*s;
                i = (m+n) % n;
                if(i == h) continue;
                j = (h-s+n) % n;
                floyd(dist,pred,i,j,i,h,h,j,s);
                j = (h+s+n) % n;
                floyd(dist,pred,i,j,i,h,h,j,s);
            }
        }


        //possible negative cycle
        for(int i=0;i<n;i++ ){
            if(dist.get(i,i) < 0) {
                throw std::invalid_argument( "negative cycle" );
            }
        }
    }
    return pred;
}