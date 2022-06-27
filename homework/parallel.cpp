#include <iostream>
#include <mpi.h>
#include "matrix.h"
#include "sequential.h"

#define TAG 21

void floydWarhsallSquaredParallel(pc::matrix<double> &dist, pc::matrix<int> &pred, int n, int s, int thread_size, int thread_rank);

void whiteFor(int s, int n, int h, int k, pc::matrix<double> &dist, pc::matrix<int> &pred, pc::matrix<int> &predFilter);

void fill(pc::matrix<int> &m, int dx, int dy, int s);

template <class T>
void plusOfSquareToLinear(pc::matrix<T> &matrix, T* dest,  int s, int n, int h, int k);

template <class T>
void plusOfSquareFromLinear(pc::matrix<T> &matrix, T* dest,  int s, int n, int h, int k);

template <class T>
void submatrixToLinear(pc::matrix<T> &m, T* dest, int dx, int dy, int s);

template <class T>
void linearToSubatrix(pc::matrix<T> &m, T* source, int dx, int dy, int s);

int main(int argc , char ** argv) {
    int thread_rank , thread_size ;
    MPI_Init (& argc , & argv );
    MPI_Comm_rank ( MPI_COMM_WORLD , & thread_rank );
    MPI_Comm_size ( MPI_COMM_WORLD , & thread_size );

    pc::matrix<double> dist = pc::matrix<double>(0);
    pc::matrix<int> pred = pc::matrix<int>(0);
    int s = 0;
    int n = 0;

    if(thread_rank == 0) {
        std::cin >> s;

        std::cin >> dist;

        pred = pc::matrix<int>(dist.getSize());

        n = dist.getSize();
    }


    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Bcast( &n, 1, MPI_INT, 0, MPI_COMM_WORLD );
    MPI_Bcast( &s, 1, MPI_INT, 0, MPI_COMM_WORLD );

    if(thread_rank != 0) {
        dist = pc::matrix<double>(n);
        pred = pc::matrix<int>(n);
        dist.fill(std::numeric_limits<double>::infinity());
    }

    floydWarhsallSquaredParallel(dist, pred, n, s,  thread_size, thread_rank);

    MPI_Barrier( MPI_COMM_WORLD );


    if(thread_rank == 0) {
        std::cout << dist << std::endl; //distances matrix: min costs
        std::cout << pred << std::endl; //predecessors matrix: min path
    }

    MPI_Finalize ();
    return 0;
}

void floydWarhsallSquaredParallel(pc::matrix<double> &dist, pc::matrix<int> &pred, int n, int s,  int thread_size, int thread_rank){
    if(((int)((double)n/(double)s))%2 == 0) {
        if(thread_rank == 0) {
            std::cout << "ERROR! work only with a odd division of matrix" << std::endl;
        }
        return;
    }

    if(thread_rank == 0) {

        // initialised predecessors
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
               pred.set(i, j, i);
            }
        }

        //computation of min path costs
        for (int h = 0; h < n; h += s) {
            pc::matrix<int> pf = pc::matrix<int>(n);
            pf.fill(0);

            // dark green self dependent
            floyd(dist, pred, h, h, h, h, h, h, s); // first
            fill(pf, h,h, s);



            for (int k = 1; k <= (n / s) / 2; k++) { // second - ....
                // light green - row h
                int j = (h - k * s + n) % n;
                floyd(dist, pred, h, j, h, h, h, j, s);
                fill(pf, h,j, s);

                j = (h + k * s + n) % n;
                floyd(dist, pred, h, j, h, h, h, j, s);
                fill(pf, h,j, s);

                // light green - row h
                int i = (h - k * s + n) % n;
                floyd(dist, pred, i, h, i, h, h, h, s);
                fill(pf, i,h, s);

                i = (h + k * s + n) % n;
                floyd(dist, pred, i, h, i, h, h, h, s);
                fill(pf, i,h, s);

                // white - rest of cells

                if (thread_size == 1) {
                    whiteFor(s, n, h, k, dist, pred, pf);
                } else {
                    int size = 12*(k+1)*s*s;
                    double * linearSquareDist = new double[size];
                    int * linearSquarePred = new int[size];
                    plusOfSquareToLinear<double>(dist, linearSquareDist, s, n, h, k);
                    plusOfSquareToLinear<int>(pred, linearSquarePred, s, n, h, k);

                    MPI_Send(&k, 1, MPI_INT, (k-1) % (thread_size -1) +1, TAG, MPI_COMM_WORLD);
                    MPI_Send(linearSquareDist, size, MPI_DOUBLE, (k-1) % (thread_size -1) +1, TAG, MPI_COMM_WORLD);
                    MPI_Send(linearSquarePred, size, MPI_INT, (k-1) % (thread_size -1) +1, TAG, MPI_COMM_WORLD);
                    delete[] linearSquareDist;
                    delete[] linearSquarePred;

                }
            }

            if (thread_size != 1) {
                MPI_Barrier(MPI_COMM_WORLD);
                pc::matrix<double> merge = pc::matrix<double>(n);
                pc::matrix<int> mergePred = pc::matrix<int>(n);

                for(int *p = pred.begin(), * ppf = pf.begin(); p < pred.end() && ppf < pf.end(); p++, ppf++) {
                    *p = *p**ppf;
                }

                MPI_Reduce(dist.begin(), merge.begin(), n*n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
                MPI_Reduce(pred.begin(), mergePred.begin(), n*n, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

                dist = merge;
                pred = mergePred;
            }
            //possible negative cycle
            for (int i = 0; i < n; i++) {
                if (dist.get(i, i) < 0) {
                    MPI_Abort(MPI_COMM_WORLD, MPI_ERR_UNKNOWN);

                    throw std::invalid_argument("negative cycle");
                }
            }
        }
    }
    else {
        int kk = (n / s) / 2;
        for (int h = 0; h < n; h += s) {
            pc::matrix<int> pf = pc::matrix<int>(n);
            pf.fill(0);

            for (int i = 0; i < kk / (thread_size -1); i++) {
                MPI_Status status;
                int k = 0;
                MPI_Recv(&k, 1, MPI_INT, MPI_ANY_SOURCE, TAG, MPI_COMM_WORLD, &status);
                int size = 12*(k+1)*s*s;
                auto * linearSquareDist = new double[size];
                int * linearSquarePred = new int[size];
                MPI_Recv(linearSquareDist, size, MPI_DOUBLE, MPI_ANY_SOURCE, TAG, MPI_COMM_WORLD, &status);
                MPI_Recv(linearSquarePred, size, MPI_INT, MPI_ANY_SOURCE, TAG, MPI_COMM_WORLD, &status);

                plusOfSquareFromLinear(dist, linearSquareDist, s, n, h, k);
                plusOfSquareFromLinear(pred, linearSquarePred, s, n, h, k);

                whiteFor(s, n, h, k, dist, pred, pf);

                delete[] linearSquareDist;
                delete[] linearSquarePred;
            }

            // if we have to do one more execution
            if((thread_rank -1) < (kk % (thread_size -1))) {
                //printf("cpu %d is ready to recive last \n", thread_rank);

                MPI_Status status;
                int k = 0;
                MPI_Recv(&k, 1, MPI_INT, MPI_ANY_SOURCE, TAG, MPI_COMM_WORLD, &status);
                int size = 12*(k+1)*s*s;
                auto * linearSquareDist = new double[size];
                int * linearSquarePred = new int[size];
                MPI_Recv(linearSquareDist, size, MPI_DOUBLE, MPI_ANY_SOURCE, TAG, MPI_COMM_WORLD, &status);
                MPI_Recv(linearSquarePred, size, MPI_INT, MPI_ANY_SOURCE, TAG, MPI_COMM_WORLD, &status);
                //printf("cpu %d recive packet with data s: %d n: %d h: %d k: %d \n", thread_rank, s, n, h, k);

                //std::cout << dist << std::endl;

                plusOfSquareFromLinear(dist, linearSquareDist, s, n, h, k);
                plusOfSquareFromLinear(pred, linearSquarePred, s, n, h, k);

                whiteFor(s, n, h, k, dist, pred, pf);

                //std::cout << dist << std::endl;

                //printf("cpu %d end of calc \n", thread_rank);

                delete[] linearSquareDist;
                delete[] linearSquarePred;
            }

            for(int *p = pred.begin(), * ppf = pf.begin(); p < pred.end() && ppf < pf.end(); p++, ppf++) {
                *p = *p**ppf;
            }

            MPI_Barrier ( MPI_COMM_WORLD );
            MPI_Reduce(dist.begin(), dist.begin(), n*n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
            MPI_Reduce(pred.begin(), pred.begin(), n*n, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        }
    }
}

// O(s^2)
void fill(pc::matrix<int> &m, int dx, int dy, int s) {
    for(int i=0; i<s;i++) {
        for (int j = 0; j < s; j++) {
            m[dx + i][dy + j] = 1;
        }
    }
}

// O(s^3* n/2s) = O(s^2*n)
void whiteFor(int s, int n, int h, int k, pc::matrix<double> &dist, pc::matrix<int> &pred, pc::matrix<int> &predFilter) {
    int i,j;
    for(int c = k - (n / s) / 2; c <= (n / s) / 2 + k; c++) {
        i = (((c + k) * s) + n) % n;
        if(i == h) continue;
        j = (h-k*s+n) % n;
        floyd(dist,pred,i,j,i,h,h,j,s);
        fill(predFilter, i,j, s);
        j = (h+k*s+n) % n;
        floyd(dist,pred,i,j,i,h,h,j,s);
        fill(predFilter, i,j, s);
    }
    for(int r = k - (n / s) / 2 + s; r <= (n / s) / 2 + k - s; r++) {
        j = (((r + k) * s) + n) % n;
        if(j == h) continue;
        i = (h-k*s+n) % n;
        floyd(dist,pred,i,j,i,h,h,j,s);
        fill(predFilter, i,j, s);
        i = (h+k*s+n) % n;
        floyd(dist,pred,i,j,i,h,h,j,s);
        fill(predFilter, i,j, s);
    }
}

template <class T> // O(s^2 * n/2s) = O(sn)
void plusOfSquareToLinear(pc::matrix<T> &matrix, T* dest,  int s, int n, int h, int k){
    int a = 0;
    int i,j;
    for(int l = k - (n / s) / 2; l <= (n / s) / 2 + k; l++) {
        int m = (l + k) * s;
        i = (m + n) % n;
        j = (h - k * s + n) % n;
        submatrixToLinear<T>(matrix, dest + a, i, j, s);
        a+= s*s;
        j = (h + n) % n;
        submatrixToLinear<T>(matrix, dest + a, i, j, s);
        a+= s*s;
        j = (h+k*s+n) % n;
        submatrixToLinear<T>(matrix, dest + a, i, j, s);
        a+= s*s;
    }
    for(int l = k-(n/s)/2+s; l<=(n/s)/2+k-s; l++) {
        int m = (l + k) * s;
        j = (m + n) % n;
        i = (h - k * s + n) % n;
        submatrixToLinear<T>(matrix, dest + a, i, j, s);
        a+= s*s;
        i = (h + n) % n;
        submatrixToLinear<T>(matrix, dest + a, i, j, s);
        a+= s*s;
        i = (h+k*s+n) % n;
        submatrixToLinear<T>(matrix, dest + a, i, j, s);
        a+= s*s;
    }
}

template <class T> // O(s^2 * n/2s) = O(sn)
void plusOfSquareFromLinear(pc::matrix<T> &matrix, T* dest,  int s, int n, int h, int k) {
    int a = 0;
    int i,j;
    for(int l = k - (n / s) / 2; l <= (n / s) / 2 + k; l++) {
        int m = (l + k) * s;
        i = (m + n) % n;
        j = (h - k * s + n) % n;
        linearToSubatrix<T>(matrix, dest + a, i, j, s);
        a+= s*s;
        j = (h + n) % n;
        linearToSubatrix<T>(matrix, dest + a, i, j, s);
        a+= s*s;
        j = (h+k*s+n) % n;
        linearToSubatrix<T>(matrix, dest + a, i, j, s);
        a+= s*s;
    }
    for(int l = k-(n/s)/2+s; l<=(n/s)/2+k-s; l++) {
        int m = (l + k) * s;
        j = (m + n) % n;
        i = (h - k * s + n) % n;
        linearToSubatrix<T>(matrix, dest + a, i, j, s);
        a+= s*s;
        i = (h + n) % n;
        linearToSubatrix<T>(matrix, dest + a, i, j, s);
        a+= s*s;
        i = (h+k*s+n) % n;
        linearToSubatrix<T>(matrix, dest + a, i, j, s);
        a+= s*s;
    }
}

template <class T> // O(n^2)
void submatrixToLinear(pc::matrix<T> &m, T* dest, int dx, int dy, int s) {
    for(int i = 0, k=0; i< s; i++, k+=s) {
        std::copy(&m[i+dx][dy], &m[i+dx][dy+s],&dest[k]);
    }
}

template <class T> // O(n^2)
void linearToSubatrix(pc::matrix<T> &m, T* source, int dx, int dy, int s) {
    for(int i = 0, k =0; i< s; i++, k+=s) {
        std::copy(&source[k], &source[k+s], &m[i+dx][dy]);

    }
}