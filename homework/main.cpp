#include <iostream>
#include <mpi.h>
#include "matrix.h"
#include "sequential.h"

int main(int argc , char ** argv) {
    int s;
    std::cin >> s;

    pc::matrix<float> cseq = pc::matrix<float>(0);
    std::cin >> cseq;
    pc::matrix<float> csq = cseq;

    double t0 , t1 ,t2, time;
    t0 = MPI_Wtime ();

    pc::matrix<int> out = floydWarhsallSeq(cseq);

    t1 = MPI_Wtime();
    pc::matrix<int> outc = floydWarhsallSquared(csq, s);

    t2 = MPI_Wtime ();
    
    time = 1.e6*( t1 - t0 );
    printf ("seq,,%f\n", time );

    time = 1.e6*( t2 - t1 );
    printf ("sq,,%f\n", time );

    if(cseq != csq) {
        printf ("not ok" );
    } else {
	if(  cseq.getSize()<= 20) std::cout << cseq << std::endl;
    }

    return 0;

}


