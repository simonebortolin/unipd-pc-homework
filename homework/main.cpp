#include <iostream>
#include <mpi.h>
#include <fstream>
#include "matrix.h"
#include "sequential.h"

int main(int argc , char ** argv) {
    int s;
    std::fstream fin(argv[1], std::ios::in);
    if(!fin){
        MPI_Abort(MPI_COMM_WORLD, 99);
    }
    fin >> s;

    pc::matrix<float> cseq = pc::matrix<float>(0);
    fin >> cseq;
    pc::matrix<float> csq = cseq;

    double t0 , t1 ,t2;
    t0 = MPI_Wtime ();

    pc::matrix<int> out = floydWarhsallSeq(cseq);

    t1 = MPI_Wtime();
    pc::matrix<int> outc = floydWarhsallSquared(csq, s);

    t2 = MPI_Wtime ();
    
    printf ("seq,,%f\n", t1 - t0 );

    printf ("sq,,%f\n", t2-t1 );

    if(cseq != csq) {
        printf ("not ok" );
    } else {
	    if(  cseq.getSize()<= 20) std::cout << cseq << std::endl;
    }

    return 0;
}


