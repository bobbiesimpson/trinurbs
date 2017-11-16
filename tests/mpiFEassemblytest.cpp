
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
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseVector.h"

int main(int argc, char *argv[]) {
    
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif
    
    int nel = 3;
    int global_dof_n = nel + 2;
    
    Epetra_Map element_map(nel, 0, Comm);
    
    int num_local_els = element_map.NumMyElements();
    int my_elements[num_local_els];
    element_map.MyGlobalElements(my_elements);
    
    std::vector<std::vector<int>> conn;
    
    for(int ilocal = 0; ilocal < element_map.NumMyElements(); ++ilocal)
    {
        int iel = my_elements[ilocal];
        
        std::vector<int> v{iel, iel+1, iel+2};
        conn.push_back(v);
    }
    
    
    // Create a standard map
    Epetra_Map map(global_dof_n,0,Comm);
    
    int num_my_rows = map.NumMyElements();
    
    int my_global_rows[num_my_rows];
    map.MyGlobalElements(my_global_rows);

    
    // now determine the number of nonzeros in each row of the matrix
    std::vector<int> num_nonzeros(num_my_rows, 0);
    
    std::vector<std::vector<int>> global_col_indices;
//    
//    for(int ilocal = 0; ilocal < num_my_rows; ++ilocal)
//    {
//        int i = my_global_rows[ilocal];
//        
//        std::vector<int> nonzero_cols;
//        
//        // loop over all the element connectivities
//        for(auto col_conn : conn)
//        {
//            std::sort(col_conn.begin(), col_conn.end());
//            
//            // discard any where the row index is less than or greater than the min or max col index
//            if(i < col_conn.front() || i > col_conn.back())
//                continue;
//            
//            auto find = std::find(col_conn.begin(), col_conn.end(), i);
//            
//            // if found, add this connectivity to the current col list
//            if(find != col_conn.end())
//                nonzero_cols.insert(nonzero_cols.end(), col_conn.begin(), col_conn.end());
//        }
//        
//        // determine the unique set of col indices for each row
//        std::sort(nonzero_cols.begin(), nonzero_cols.end());
//        auto last = std::unique(nonzero_cols.begin(), nonzero_cols.end());
//        nonzero_cols.erase(last, nonzero_cols.end());
//        global_col_indices.push_back(nonzero_cols);
//        
//        // the size of this unique set is the number of non-zeros for this row
//        num_nonzeros[ilocal] = nonzero_cols.size();
//    }
//    
    // Create the appropriate graph from the vector of nonzero entries
    //Epetra_FECrsGraph graph(Copy, map, num_nonzeros.data(), true);
    
    Epetra_FECrsGraph graph(Copy, map, 3);

    
//    for(int ilocal_el = 0; ilocal_el < num_local_els; ++ilocal_el)
//    {
//        const auto& el_conn = conn[ilocal_el];
//        
//        std::cout << "inserting global indices for local element: " << ilocal_el << "\n";
//        for(const auto& index : el_conn)
//            std::cout << index << "\t";
//        std::cout << "\n";
//        
//        graph.InsertGlobalIndices(el_conn.size(), el_conn.data(), el_conn.size(), el_conn.data());
//        
//    }
    
//    // Fill in the relevant column indices
//    for(int i = 0; i < global_col_indices.size(); ++i)
//        graph.InsertGlobalIndices(i, global_col_indices[i].size(), global_col_indices[i].data());
    
    graph.FillComplete();

//    if(0 == Comm.MyPID())
//        std::cout << graph << "\n";
    
    // now build the Stiffness matrix
    Epetra_FECrsMatrix K(Copy, map, 3);
    
    // lop over elements on this process
    for(int ilocal = 0; ilocal < num_local_els; ++ilocal)
    {
//        int iel = my_elements[ilocal];
        
        auto& el_conn = conn[ilocal];
        //std::cout << "on element: " << my_elements[ilocal] << " ";
//        for(const auto& index : el_conn)
//            std::cout << index << "\t";
//        std::cout << "\n";
        
        Epetra_SerialDenseMatrix k_local(el_conn.size(), el_conn.size());
        Epetra_IntSerialDenseVector scatter(Copy, el_conn.data(), el_conn.size());
        
        for(int irow = 0; irow < el_conn.size(); ++irow)
            for(int icol = 0; icol < el_conn.size(); ++icol)
                k_local(irow, icol) = 1.0;
        K.InsertGlobalValues(scatter, k_local);
    }

    K.GlobalAssemble();
    
    std::cout << "done!!\n";

    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    
    return(EXIT_SUCCESS);
    
}
