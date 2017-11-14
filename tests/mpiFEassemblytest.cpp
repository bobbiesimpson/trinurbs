
#include <iostream>

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FECrsGraph.h"

int main(int argc, char *argv[]) {
    
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif
    
    //
    int NumGlobalElements = 10;
    Epetra_Map Map(NumGlobalElements,0,Comm);
    
    Epetra_FECrsGraph graph;
    
    // create a diagonal FE crs matrix (one nonzero per row)
    Epetra_FECrsMatrix K(Copy,Map,1);
    
    //
    // Note 2: We fill the matrix using 'InsertGlobalValues'. An
    // alternative approach that would be more efficient for large
    // matrices in most cases would be to first create and fill a
    // graph (Epetra_FECrsGraph), then construct the matrix with the
    // graph (after calling graph.FillComplete) and fill the matrix
    // using the method 'SumIntoGlobalValues'.
    //
//    if( Comm.MyPID() == 0 ) {
//        for( int i=0 ; i<NumGlobalElements ; ++i ) {
//            int index = i;
//            double value = 1.0*i;
//            A.InsertGlobalValues(1,&index,&value);
//        }
//    }
//    
//    A.GlobalAssemble();
//    
//    std::cout << A;
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    
    return(EXIT_SUCCESS);
    
}
