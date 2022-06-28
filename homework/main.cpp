#include <iostream>
#include <mpi.h>
#include "matrix.h"
#include "sequential.h"

int main(int argc , char ** argv) {
    int s;
    std::cin >> s;

    pc::matrix<double> cseq = pc::matrix<double>(0);
    std::cin >> cseq;
    pc::matrix<double> csq = cseq;

    double t0 , t1 , time;
    t0 = MPI_Wtime ();

    pc::matrix<int> out = floydWarhsallSeq(cseq);
    pc::matrix<int> outc = floydWarhsallSquared(csq,s);

    t1 = MPI_Wtime ();
    time = 1.e6 * ( t1 - t0 );
    printf (" That took %f seconds \n ", time );

    if(cseq != csq) {
        printf ("not ok - time: %d", time );
    } else {
        printf ("pass - time: %d", time );
    }

    return 0;

}


