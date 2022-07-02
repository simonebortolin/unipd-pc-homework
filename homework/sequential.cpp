//
// Created by simon on 25/06/2022.
//
#include "matrix.h"
#include "sequential.h"

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
                throw std::invalid_argument( "negative cycle" );
            }
        }
        //std::cout << dist << std::endl;
    }
    return pred;
}

// s: submatrix size
// dx,dy: start index of submatrix d
// ax,ay: start index of submatrix a
// bx,by: start index of submatrix b
void floyd(pc::matrix<double> &dist, pc::matrix<int> &pred, int dx, int dy, int ax, int ay, int bx, int by, int s) {
    double tr;
    int pr;
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
    int n = dist.getSize();

    pc::matrix<int> pred = pc::matrix<int>(n);


    if(((int)((double)n/(double)s))%2 == 0) {
        std::cout << "ERROR! work only with a odd division of matrix";
        return pc::matrix<int>(0);
    }


    // initialised predecessors
    for(int i=0; i<n;i++){
        for(int j=0; j<n;j++){
            pred.set(i,j,i);
        }
    }

    //computation of min path costs
    for(int h=0; h<n;h+=s){
        // dark green self dependent
        floyd(dist, pred, h,h,h,h,h,h,s); // first


        for(int k=1; k<=(n/s)/2;k++){ // second - ....
            // light green - row h
            int j = (h-k*s+n) % n;
            floyd(dist,pred,h,j,h,h,h,j,s);

            j = (h+k*s+n) % n;
            floyd(dist,pred,h,j,h,h,h,j,s);

            // light green - row h
            int i = (h-k*s+n) % n;
            floyd(dist, pred, i, h, i, h, h, h,s);

            i = (h+k*s+n) % n;
            floyd(dist, pred, i, h, i, h, h, h,s);

            // white - rest of cells

            for(int c = k - (n / s) / 2; c <= (n / s) / 2 + k; c++) {
                i = (((c + k) * s) + n) % n;
                if(i == h) continue;
                j = (h-k*s+n) % n;
                floyd(dist,pred,i,j,i,h,h,j,s);
                j = (h+k*s+n) % n;
                floyd(dist,pred,i,j,i,h,h,j,s);
            }
            for(int r = k - (n / s) / 2 + s; r <= (n / s) / 2 + k - s; r++) {
                j = (((r + k) * s) + n) % n;
                if(j == h) continue;
                i = (h-k*s+n) % n;
                floyd(dist,pred,i,j,i,h,h,j,s);
                i = (h+k*s+n) % n;
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
