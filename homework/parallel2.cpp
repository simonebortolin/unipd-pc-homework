#include <iostream>
#include <mpi.h>
#include "matrix.h"
#include "matrix_matrix.h"

#define TAG 21

void floyd(pc::basic_matrix<float> &dist, pc::basic_matrix<float> &a, pc::basic_matrix<float> &b, pc::basic_matrix<int> &preda, pc::basic_matrix<int> &predb);
template <class T>
void to_blocked_matrix(pc::matrix<T> & from, pc::matrix_matrix<T> & to, int s);
template <class T>
void to_linear_matrix(pc::matrix_matrix<T> & from, pc::matrix<T>  & to, int s);

int main(int argc , char ** argv) {
    int thread_rank , thread_size ;
    MPI_Init (& argc , & argv );
    MPI_Comm_rank ( MPI_COMM_WORLD , & thread_rank );
    MPI_Comm_size ( MPI_COMM_WORLD , & thread_size );

    pc::matrix<float> dist = pc::matrix<float>(0);
    pc::matrix<int> pred = pc::matrix<int>(0);
    int s = 0;
    int n = 0;

    int b[2];
    int position = 0;

    if(thread_rank == 0) {
        std::cin >> s;

        std::cin >> dist;

        pred = pc::matrix<int>(dist.getSize());

        n = dist.getSize();

        MPI_Pack(&n, 1, MPI_INT, b, sizeof(b), &position, MPI_COMM_WORLD);
        MPI_Pack(&s, 1, MPI_INT, b, sizeof(b), &position, MPI_COMM_WORLD);
    }


    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Bcast(&b, sizeof(b), MPI_PACKED, 0, MPI_COMM_WORLD );
    if(thread_rank != 0) {
        MPI_Unpack(b, sizeof(b), &position, &n, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack(b, sizeof(b), &position, &s, 1, MPI_INT, MPI_COMM_WORLD);
    }

    if(thread_rank == 0) {
        if(((int)((float)n/(float)s))%2 == 0) {
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
        dist = pc::matrix<float>(n);
        pred = pc::matrix<int>(n);
    }
    MPI_Barrier( MPI_COMM_WORLD );



    MPI_Bcast( dist.begin(), n*n, MPI_FLOAT, 0, MPI_COMM_WORLD );
    MPI_Bcast( pred.begin(), n*n, MPI_INT, 0, MPI_COMM_WORLD );

    double t0 , t1 , time;


    pc::matrix_matrix<float> d = pc::matrix_matrix<float>(n , s);
    pc::matrix_matrix<int> p = pc::matrix_matrix<int>(n , s);
    to_blocked_matrix(dist, d, s);
    to_blocked_matrix(pred, p, s);

    t0 = MPI_Wtime();


    //computation of min path costs
    for (int h = 0; h < n/s; h ++) {
        int size_of_buffer = s*s*(sizeof(int)+sizeof(float));
        char * buffer = new char[size_of_buffer];
        position = 0;
        if(thread_rank == 0) {
            // dark green self dependent
            floyd(d[h][h], d[h][h], d[h][h], p[h][h], p[h][h]); // first
            MPI_Pack(d[h][h].begin(), s*s, MPI_FLOAT, buffer,  size_of_buffer, &position, MPI_COMM_WORLD);
            MPI_Pack(p[h][h].begin(), s*s, MPI_INT, buffer,  size_of_buffer, &position, MPI_COMM_WORLD);
        }
        MPI_Bcast(buffer, size_of_buffer, MPI_PACKED, 0, MPI_COMM_WORLD);
        if(thread_rank != 0) {
            MPI_Unpack(buffer, size_of_buffer, &position, d[h][h].begin(), s*s, MPI_FLOAT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, size_of_buffer, &position, p[h][h].begin(), s*s, MPI_INT, MPI_COMM_WORLD);
        }
        delete[] buffer;
        buffer = nullptr;

        int ns = (n / s);
        for(int k = 1; k <= ns / 2; k++) {
            size_of_buffer = 4*s*s*(sizeof(int)+sizeof(float));
            buffer = new char[size_of_buffer];
            position = 0;
            if(thread_rank == 0) {
                int j = (h - k + ns) % ns;
                floyd(d[h][j], d[h][h], d[h][j], p[h][h], p[h][j]);
                MPI_Pack(d[h][j].begin(), s*s, MPI_FLOAT, buffer,  size_of_buffer, &position, MPI_COMM_WORLD);
                MPI_Pack(p[h][j].begin(), s*s, MPI_INT, buffer,  size_of_buffer, &position, MPI_COMM_WORLD);
                j = (h + k + ns) % ns;
                floyd(d[h][j], d[h][h], d[h][j], p[h][h], p[h][j]);
                MPI_Pack(d[h][j].begin(), s*s, MPI_FLOAT, buffer,  size_of_buffer, &position, MPI_COMM_WORLD);
                MPI_Pack(p[h][j].begin(), s*s, MPI_INT, buffer,  size_of_buffer, &position, MPI_COMM_WORLD);
                int i = (h - k + ns) % ns;
                floyd(d[i][h], d[i][h], d[h][h], p[i][h], p[h][h]);
                MPI_Pack(d[i][h].begin(), s*s, MPI_FLOAT, buffer,  size_of_buffer, &position, MPI_COMM_WORLD);
                MPI_Pack(p[i][h].begin(), s*s, MPI_INT, buffer,  size_of_buffer, &position, MPI_COMM_WORLD);
                i = (h + k + ns) % ns;
                floyd(d[i][h], d[i][h], d[h][h], p[i][h], p[h][h]);
                MPI_Pack(d[i][h].begin(), s*s, MPI_FLOAT, buffer,  size_of_buffer, &position, MPI_COMM_WORLD);
                MPI_Pack(p[i][h].begin(), s*s, MPI_INT, buffer,  size_of_buffer, &position, MPI_COMM_WORLD);
            }
            MPI_Barrier( MPI_COMM_WORLD );
            MPI_Bcast(buffer, size_of_buffer, MPI_PACKED, 0, MPI_COMM_WORLD);
            if(thread_rank != 0) {
                int j = (h - k + ns) % ns;
                MPI_Unpack(buffer, size_of_buffer, &position, d[h][j].begin(), s * s, MPI_FLOAT, MPI_COMM_WORLD);
                MPI_Unpack(buffer, size_of_buffer, &position, p[h][j].begin(), s * s, MPI_INT,  MPI_COMM_WORLD);
                j = (h + k + ns) % ns;
                MPI_Unpack(buffer, size_of_buffer, &position, d[h][j].begin(), s * s, MPI_FLOAT, MPI_COMM_WORLD);
                MPI_Unpack(buffer, size_of_buffer, &position, p[h][j].begin(), s * s, MPI_INT, MPI_COMM_WORLD);
                int i = (h - k + ns) % ns;
                MPI_Unpack(buffer, size_of_buffer, &position, d[i][h].begin(), s * s, MPI_FLOAT, MPI_COMM_WORLD);
                MPI_Unpack(buffer, size_of_buffer, &position, p[i][h].begin(), s * s, MPI_INT, MPI_COMM_WORLD);
                i = (h + k + ns) % ns;
                MPI_Unpack(buffer, size_of_buffer, &position, d[i][h].begin(), s * s, MPI_FLOAT, MPI_COMM_WORLD);
                MPI_Unpack(buffer, size_of_buffer, &position, p[i][h].begin(), s * s, MPI_INT, MPI_COMM_WORLD);
            }

            delete[] buffer;
            buffer = nullptr;
            if(thread_rank == k % (thread_size -1) +1) {
                size_of_buffer = 4*s*s*(sizeof(int)+sizeof(float))*k*k;
                buffer = new char[size_of_buffer];
                position = 0;
                for (int jj = h - k; jj <= h + k; jj++) {
                    int j = (jj + ns) % ns;
                    if (j == h) continue;
                    int i = (h - k + ns) % ns;
                    floyd(d[i][j], d[i][h], d[h][j], p[i][j], p[h][j]);
                    MPI_Pack(d[i][j].begin(), s*s, MPI_FLOAT, buffer,  size_of_buffer, &position, MPI_COMM_WORLD);
                    MPI_Pack(p[i][j].begin(), s*s, MPI_INT, buffer,  size_of_buffer, &position, MPI_COMM_WORLD);
                    i = (h + k + ns) % ns;
                    floyd(d[i][j], d[i][h], d[h][j], p[i][j], p[h][j]);
                    MPI_Pack(d[i][j].begin(), s*s, MPI_FLOAT, buffer,  size_of_buffer, &position, MPI_COMM_WORLD);
                    MPI_Pack(p[i][j].begin(), s*s, MPI_INT, buffer,  size_of_buffer, &position, MPI_COMM_WORLD);
                }
                for (int ii = h - k + 1; ii < k + h; ii++) {
                    int i = (ii + ns) % ns;
                    if (i == h) continue;
                    int j = (h - k + ns) % ns;
                    floyd(d[i][j], d[i][h], d[h][j], p[i][j], p[h][j]);
                    MPI_Pack(d[i][j].begin(), s*s, MPI_FLOAT, buffer,  size_of_buffer, &position, MPI_COMM_WORLD);
                    MPI_Pack(p[i][j].begin(), s*s, MPI_INT, buffer,  size_of_buffer, &position, MPI_COMM_WORLD);
                    j = (h + k + ns) % ns;
                    floyd(d[i][j], d[i][h], d[h][j], p[i][j], p[h][j]);
                    MPI_Pack(d[i][j].begin(), s*s, MPI_FLOAT, buffer,  size_of_buffer, &position, MPI_COMM_WORLD);
                    MPI_Pack(p[i][j].begin(), s*s, MPI_INT, buffer,  size_of_buffer, &position, MPI_COMM_WORLD);
                }
                MPI_Send(buffer, position, MPI_PACKED, 0, TAG, MPI_COMM_WORLD);

                delete[] buffer;
                buffer = nullptr;
            }
            if(thread_rank == 0) {
                size_of_buffer = 4*s*s*(sizeof(int)+sizeof(float))*k*k;
                buffer = new char[size_of_buffer];
                position = 0;
                MPI_Status status;
                MPI_Recv(buffer, size_of_buffer, MPI_PACKED, k % (thread_size -1) +1, TAG, MPI_COMM_WORLD, &status);
                for (int jj = h - k; jj <= h + k; jj++) {
                    int j = (jj + ns) % ns;
                    if (j == h) continue;
                    int i = (h - k + ns) % ns;
                    MPI_Unpack(buffer, size_of_buffer, &position,d[i][j].begin(), s * s, MPI_FLOAT, MPI_COMM_WORLD);
                    MPI_Unpack(buffer, size_of_buffer, &position,p[i][j].begin(), s * s, MPI_INT, MPI_COMM_WORLD);
                    i = (h + k + ns) % ns;
                    MPI_Unpack(buffer, size_of_buffer, &position,d[i][j].begin(), s * s, MPI_FLOAT, MPI_COMM_WORLD);
                    MPI_Unpack(buffer, size_of_buffer, &position,p[i][j].begin(), s * s, MPI_INT, MPI_COMM_WORLD);
                }
                for (int ii = h - k + 1; ii < k + h; ii++) {
                    int i = (ii + ns) % ns;
                    if (i == h) continue;
                    int j = (h - k + ns) % ns;
                    MPI_Unpack(buffer, size_of_buffer, &position,d[i][j].begin(), s * s, MPI_FLOAT, MPI_COMM_WORLD);
                    MPI_Unpack(buffer, size_of_buffer, &position,p[i][j].begin(), s * s, MPI_INT, MPI_COMM_WORLD);
                    j = (h + k + ns) % ns;
                    MPI_Unpack(buffer, size_of_buffer, &position,d[i][j].begin(), s * s, MPI_FLOAT, MPI_COMM_WORLD);
                    MPI_Unpack(buffer, size_of_buffer, &position,p[i][j].begin(), s * s, MPI_INT, MPI_COMM_WORLD);
                }

                delete[] buffer;
                buffer = nullptr;
            }
        }
        MPI_Barrier( MPI_COMM_WORLD );
    }

    t1 = MPI_Wtime();
    time = 1.e6*( t1 - t0 );
    if(thread_rank==0) {
        printf ("par2sq,%d,%f\n", thread_size, time );
        if(n <= 20) {
            to_linear_matrix(d, dist, s);
            std::cout << dist << std::endl;
        }
    }
    MPI_Finalize ();
    return 0;
}

template <class T>
void to_blocked_matrix(pc::matrix<T> & from, pc::matrix_matrix<T> & to, int s) {
    int n = from.getSize();
    for(int dx = 0, x=0; dx < n; dx+= s, x++) {
        for(int dy = 0, y=0; dy < n; dy+= s, y++) {
            for(int k=0; k<s;k++) {
                for(int h=0; h< s; h++) {
                    to[x][y][k][h] = from[dx+k][dy+h];
                }
            }
        }
    }

}

template <class T>
void to_linear_matrix(pc::matrix_matrix<T> & from, pc::matrix<T>  & to, int s) {
    int n = to.getSize();
    for(int dx = 0, x=0; dx < n; dx+= s, x++) {
        for(int dy = 0, y=0; dy < n; dy+= s, y++) {
            for(int k=0; k<s;k++) {
                for(int h=0; h< s; h++) {
                    to[dx+k][dy+h] = from[dx/s][dy/s][k][h];
                }
            }
        }
    }
}
void floyd(pc::basic_matrix<float> &dist, pc::basic_matrix<float> &a, pc::basic_matrix<float> &b, pc::basic_matrix<int> &preda, pc::basic_matrix<int> &predb) {
    int  s = dist.getSize();
    float tr;
    int pr;
    for(int h=0; h<s;h++){
        for(int i=0; i<s;i++){
            for(int j=0; j<s;j++){

                //path comparison between current i->j and i->h->j
                tr = a.get(i,h) + b.get(h,j);

                if(tr < dist.get(i, j) ){

                    //update of min path distance and predecessors
                    dist.set(i, j, tr);
                    pr = predb.get(h,j);
                    preda.set(i,j,pr);

                }
            }
        }

    }
}
