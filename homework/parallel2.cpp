#include <iostream>
#include <mpi.h>
#include "matrix.h"
#include "sequential.h"

#define TAG 21

void floyd(pc::matrix<double> &d, pc::matrix<double> &a, pc::matrix<double> &b, pc::matrix<int> &preda, pc::matrix<int> &predb);
template <class T>
void to_blocked_matrix(pc::matrix<T> & from, pc::matrix<pc::matrix<T>> & to, int s);
template <class T>
void from_blocked_matrix(pc::matrix<pc::matrix<T>> & from, pc::matrix<T>  & to, int s);
int gi = 0;

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
    if(thread_rank == 0) {
        if(((int)((double)n/(double)s))%2 == 0) {
            std::cout << "ERROR! work only with a odd division of matrix" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
            return -1;
        }

        // initialised predecessors
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                pred.set(i, j, j);
            }
        }
    }
    else {
        dist = pc::matrix<double>(dist.getSize());
        pred = pc::matrix<int>(dist.getSize());
    }
    MPI_Barrier( MPI_COMM_WORLD );

    MPI_Bcast( dist.begin(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( pred.begin(), 1, MPI_INT, 0, MPI_COMM_WORLD );

    double t0 , t1 , time;

    t0 = MPI_Wtime();

    pc::matrix< pc::matrix<double>> d = pc::matrix< pc::matrix<double>>(n/s);
    pc::matrix< pc::matrix<int>> p = pc::matrix< pc::matrix<int>>(n/s);
    to_blocked_matrix(dist, d, s);
    to_blocked_matrix(pred, p, s);
    std::cout << p << std::endl;


    //computation of min path costs
    for (int h = 0; h < n/s; h ++) {
        if(thread_rank == 0) {
            // dark green self dependent
            floyd(d[h][h], d[h][h], d[h][h], p[h][h], p[h][h]); // first
        }
        MPI_Barrier( MPI_COMM_WORLD );
        MPI_Bcast(d[h][h].begin(), s*s, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(p[h][h].begin(), s*s, MPI_INT, 0, MPI_COMM_WORLD);
        int ns = (n / s);
        for(int k = 1; k <= ns / 2; k++) {
            int j = (h - 1 + ns) % ns;
            if(thread_rank == (k*s) % thread_size) floyd(d[h][j], d[h][h], d[h][j], p[h][h], p[h][j]);
            MPI_Bcast(d[h][j].begin(), s*s, MPI_DOUBLE, (k*s) % thread_size, MPI_COMM_WORLD);
            MPI_Bcast(p[h][j].begin(), s*s, MPI_INT, (k*s) % thread_size, MPI_COMM_WORLD);
            j = (h + 1 + ns) % ns;
            if(thread_rank == (k*s+1) % thread_size) floyd(d[h][j], d[h][h], d[h][j], p[h][j], p[h][j]);
            MPI_Bcast(d[h][j].begin(), s*s, MPI_DOUBLE, (k*s+1) % thread_size, MPI_COMM_WORLD);
            MPI_Bcast(p[h][j].begin(), s*s, MPI_INT, (k*s+1) % thread_size, MPI_COMM_WORLD);
            int i = (h - 1 + ns) % ns;
            if(thread_rank == (k*s+2) % thread_size) floyd(d[i][h], d[i][h], d[h][h], p[i][h], p[h][h]);
            MPI_Bcast(d[i][h].begin(), s*s, MPI_DOUBLE, (k*s+2) % thread_size, MPI_COMM_WORLD);
            MPI_Bcast(p[i][h].begin(), s*s, MPI_INT, (k*s+2) % thread_size, MPI_COMM_WORLD);
            i = (h + 1 + ns) % ns;
            if(thread_rank == (k*s+3) % thread_size) floyd(d[i][h], d[i][h], d[h][h], p[i][h], p[h][h]);
            MPI_Bcast(d[i][h].begin(), s*s, MPI_DOUBLE, (k*s+3) % thread_size, MPI_COMM_WORLD);
            MPI_Bcast(p[i][h].begin(), s*s, MPI_INT, (k*s+3) % thread_size, MPI_COMM_WORLD);

            for(int jj = h-k; jj <= h+k; jj++) {
                j = (jj+ns) % ns;
                if(j==h) continue;
                i = (h - k + ns) % ns;
                if(thread_rank == 0) floyd(d[i][j], d[i][h], d[h][j], p[i][j], p[h][j]);
                MPI_Bcast(d[h][j].begin(), s*s, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                MPI_Bcast(p[h][j].begin(), s*s, MPI_INT, 0, MPI_COMM_WORLD);
                i = (h + k + ns) % ns;
                if(thread_rank == 0) floyd(d[i][j], d[i][h], d[h][j], p[i][j], p[h][j]);
                MPI_Bcast(d[h][j].begin(), s*s, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                MPI_Bcast(p[h][j].begin(), s*s, MPI_INT, 0, MPI_COMM_WORLD);
            }
            for(int ii = h-k +1; ii < k+h; ii++) {
                i = (ii+ns) % ns;
                if(i==h) continue;
                j = (h - k + ns) % ns;
                if(thread_rank == 0) floyd(d[i][j], d[i][h], d[h][j], p[i][j], p[h][j]);
                MPI_Bcast(d[h][j].begin(), s*s, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                MPI_Bcast(p[h][j].begin(), s*s, MPI_INT, 0, MPI_COMM_WORLD);
                j = (h + k + ns) % ns;
                if(thread_rank == 0) floyd(d[i][j], d[i][h], d[h][j], p[i][j], p[h][j]);
                MPI_Bcast(d[h][j].begin(), s*s, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                MPI_Bcast(p[h][j].begin(), s*s, MPI_INT, 0, MPI_COMM_WORLD);
            }

        }
        MPI_Barrier( MPI_COMM_WORLD );

        from_blocked_matrix(d,dist,s);
        std::cout << dist << std::endl;
    }

    t1 = MPI_Wtime();
    if(thread_rank == 0) {
        from_blocked_matrix(d,dist,s);
        std::cout << dist << std::endl;
    }

}

template <class T>
void to_blocked_matrix(pc::matrix<T> & from, pc::matrix<pc::matrix<T>> & to, int s) {
    int n = from.getSize();
    for(int dx = 0, x=0; dx < n; dx+= s, x++) {
        for(int dy = 0, y=0; dy < n; dy+= s, y++) {
            to.set(x, y, pc::matrix<T>(s));
            for(int k=0; k<s;k++) {
                std::copy(&from[dx+k][dy], &from[dx+k][dy] +s, to[x][y].begin() + k*s);
            }
        }
    }

}

template <class T>
void from_blocked_matrix(pc::matrix<pc::matrix<T>> & from, pc::matrix<T>  & to, int s) {
    int n = from.getSize();
    for(int dx = 0, x=0; dx < n; dx+= s, x++) {
        for(int dy = 0, y=0; dy < n; dy+= s, y++) {
            for(int k=0; k<s;k++) {
               std::copy(from[x][y].begin() + k*s, from[x][y].begin() + k*s + s, &to[dx+k][dy]);
            }
        }
    }
}
void floyd(pc::matrix<double> &d, pc::matrix<double> &a, pc::matrix<double> &b, pc::matrix<int> &preda, pc::matrix<int> &predb) {
    gi++;
    int  s = d.getSize();
    double tr;
    int pr;
    for(int h=0; h<s;h++){
        for(int i=0; i<s;i++){
            for(int j=0; j<s;j++){

                //path comparison between current i->j and i->h->j
                tr = a.get(i,h) + b.get(h,j);

                if(  tr < d.get(i,j) ){

                    //update of min path distance and predecessors
                    d.set(i,j,tr);
                    pr = predb.get(h,j);
                    preda.set(i,j,pr);

                }
            }
        }

    }
}