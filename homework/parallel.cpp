#include <iostream>
#include <mpi.h>
#include "matrix.h"
#include "sequential.h"

#define TAG 21

void floydWarhsallSquaredParallel(pc::matrix<double> &dist, pc::matrix<int> &pred, int s, int thread_size, int thread_rank);

void whiteFor(int [], pc::matrix<double> &dist, pc::matrix<int> &pred);

int main(int argc , char ** argv) {
    int thread_rank , thread_size ;
    MPI_Init (& argc , & argv );
    MPI_Comm_rank ( MPI_COMM_WORLD , & thread_rank );
    MPI_Comm_size ( MPI_COMM_WORLD , & thread_size );
    printf (" Hello world from process %d of %d\n", thread_rank , thread_size );

    pc::matrix<double> dist = pc::matrix<double>(0);
    int s;

    if(thread_rank == 0) {
        std::cin >> s;

        std::cin >> dist;
    }

    pc::matrix<int> pred = pc::matrix<int>(dist.getSize());

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Bcast( &s, 1, MPI_INT, 0, MPI_COMM_WORLD );
    MPI_Bcast( &dist, sizeof(dist), MPI_BYTE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &pred, sizeof(pred), MPI_BYTE, 0, MPI_COMM_WORLD );

    floydWarhsallSquaredParallel(dist, pred, s,  thread_size, thread_rank);

    if(thread_rank == 0) {
        std::cout << dist << std::endl; //distances matrix: min costs
        std::cout << pred << std::endl; //predecessors matrix: min path
    }

    MPI_Finalize ();
    return 0;
}

void floydWarhsallSquaredParallel(pc::matrix<double> &dist, pc::matrix<int> &pred, int s,  int thread_size, int thread_rank){
    int n = dist.getSize();

    if(((int)((double)n/(double)s))%2 == 0) {
        if(thread_rank == 0) {
            std::cout << "ERROR! work only with a odd division of matrix";
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
            // dark green self dependent
            floyd(dist, pred, h, h, h, h, h, h, s); // first


            for (int k = 1; k <= (n / s) / 2; k++) { // second - ....
                // light green - row h
                int j = (h - k * s + n) % n;
                floyd(dist, pred, h, j, h, h, h, j, s);

                j = (h + k * s + n) % n;
                floyd(dist, pred, h, j, h, h, h, j, s);

                // light green - row h
                int i = (h - k * s + n) % n;
                floyd(dist, pred, i, h, i, h, h, h, s);

                i = (h + k * s + n) % n;
                floyd(dist, pred, i, h, i, h, h, h, s);

                // white - rest of cells

                int data[] = {s, n, h, k};
                if (thread_size == 1) {
                    whiteFor(data, dist, pred);
                } else {
                    MPI_Send(data, 4, MPI_INT, k % (thread_size -1) +1, TAG, MPI_COMM_WORLD);

                }
            }

            if (thread_size != 1) {
                MPI_Barrier(MPI_COMM_WORLD);
            }
            //possible negative cycle
            for (int i = 0; i < n; i++) {
                if (dist.get(i, i) < 0) {
                    throw std::invalid_argument("negative cycle");
                }
            }
        }
    }
    else {
        int kk = (n / s) / 2;
        for (int h = 0; h < n; h += s) {
            for (int i = 0; i < kk / (thread_size -1); i++) {
                int data[4];
                MPI_Status status;
                MPI_Recv(data, 4, MPI_INT, MPI_ANY_SOURCE, TAG, MPI_COMM_WORLD, &status);
                whiteFor(data, dist, pred);

            }
            if(thread_rank - 1 < kk % (thread_size -1)) {
                int data[4];
                MPI_Status status;
                MPI_Recv(data, 4, MPI_INT, MPI_ANY_SOURCE, TAG, MPI_COMM_WORLD, &status);
                whiteFor(data, dist, pred);
            }
            MPI_Barrier ( MPI_COMM_WORLD );
        }
    }
}

void whiteFor(int args[], pc::matrix<double> &dist, pc::matrix<int> &pred) {
    int s = args[0];
    int n = args[1];
    int h = args[2];
    int k = args[3];
    int i,j;
    for(int l = k - (n / s) / 2; l <= (n / s) / 2 + k; l++) {
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
}
