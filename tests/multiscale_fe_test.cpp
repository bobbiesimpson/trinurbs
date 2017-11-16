#include <iostream>

#include "Forest.h"
#include "Geometry.h"
#include "MultiscaleForest.h"
#include "base.h"
#include "OutputVTK.h"
#include "IElemIntegrate.h"
#include "PeriodicForest.h"

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"

using namespace trinurbs;

// function that we use for projection
double myfunc(double x,
              double y,
              double z)
{
    return x;
    //return std::sqrt((x-0.5) * (x-0.5) + (y-0.5) * (y-0.5) +(z-0.5) * (z-0.5) );
}

uint refinementInput(const char* c)
{
    uint refine = 0;
    const uint max_refine = 10;
    
    auto input = std::atoi(c);
    
    if(input < 0)
        std::cout << "Cannot supplied negative refinement integer. Carrying on with no h-refinement";
    if(input > max_refine) {
        std::cout << "Truncating refinement to " << max_refine << " levels.";
        refine = max_refine;
    }
    else {
        std::cout << "Applying " << input << " levels of h-refinement\n";
        refine = input;
    }
    return refine;
}

int main(int argc, char* argv[])
{
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif
    
    //
    // Read in input files and generate IGA discretisation
    //
    std::ifstream ifs(argv[1]);
    if(!ifs)
        trinurbs::error("Cannot open macro file for reading\n");
    
    trinurbs::Geometry macro_geom;
    if(!macro_geom.loadHBSFile(ifs))
        trinurbs::error("Error loading macro geometry file");
    
    std::ifstream ifs2(argv[2]);
    if(!ifs2)
        trinurbs::error("Cannot open micro geometry file for reading\n");
    
    trinurbs::Geometry micro_geom;
    if(!micro_geom.loadHBSFile(ifs2))
        trinurbs::error("Error loading micro geometry file\n");
    
    MultiscaleForest multiscaleforest(macro_geom, micro_geom);
    
    uint macro_refine = 0;
    uint micro_refine = 0;
    
    if(argc > 3)
        macro_refine = refinementInput(argv[3]);
    if(argc > 4)
        micro_refine = refinementInput(argv[4]);
    
    PeriodicForest p_forest(micro_geom);
    
    multiscaleforest.hrefine(macro_refine, micro_refine);
    
    std::cout << "Multiscale geometry with "
    << multiscaleforest.macroForest().elemN() << " elements (macro), "
    << multiscaleforest.microForest().elemN() << " elements (micro), "
    << multiscaleforest.elemN() << " elements total\n\n";
    
    std::cout << "Process: " << Comm.MyPID() << "\n";
    
    long long ndof = multiscaleforest.globalDofN();
    
    // Map for mpi communication of matrix and vector terms
    Epetra_Map map(ndof, 0, Comm);
    
    // Construct sparse matrix
    Epetra_FECrsMatrix K(Copy, map, 5);
    
    // Construct RHS vector
    Epetra_FEVector f(map);
    
    // Set up element map for MPI assembly routines
    long long nel = multiscaleforest.elemN();
    Epetra_Map element_map(nel, 0, Comm);
    
    // Get data structure corresponding to local elements to this process
    int num_local_els = element_map.NumMyElements();
    long long* my_elements = element_map.MyGlobalElements64();

    // Now loop over 'local' elements
    for(auto ilocal = 0; ilocal < num_local_els; ++ilocal)
    {
        // global element index
        auto iel = my_elements[ilocal];
        
        const auto bel = multiscaleforest.multiscaleBezierElement(iel);
        
        const auto unsigned_conn = bel->globalBasisIVec();
        
        std::vector<int> conn;
        for(const auto& uval : unsigned_conn)
            conn.push_back(static_cast<int>(uval));
        
        // local stiffness matrix terms
        //Epetra_SerialDenseMatrix k_local(conn.size(), conn.size());
        Epetra_SerialDenseMatrix k_local(2,2);
        k_local(0,0) = 1.0; k_local(1,1) = 1.0;
        
        // local force vector terms
        Epetra_SerialDenseVector f_local(conn.size());
        //Epetra_IntSerialDenseVector scatter(Copy, conn.data(), conn.size());
        Epetra_IntSerialDenseVector scatter(2);
        scatter(0) = 0; scatter(1) = 1;
        
        
        // Perform quadrature
//        for(trinurbs::IElemIntegrate igpt(bel->integrationOrder()); !igpt.isDone(); ++igpt)
//        {
//            const auto gpt = igpt.get();
//            const auto basis = bel->basis(gpt.xi, gpt.eta, gpt.zeta);
//            const auto jdet = bel->jacDet(gpt.xi, gpt.eta, gpt.zeta);
//            const auto phys_coord = bel->eval(gpt.xi, gpt.eta, gpt.zeta);
//            const double funcval = myfunc(phys_coord[0], phys_coord[1], phys_coord[2]);
//            
//            for(int itest = 0; itest < conn.size(); ++itest)
//            {
//                const auto gtest_i = conn[itest];
//                for(int itrial = 0; itrial < conn.size(); ++itrial)
//                    k_local(itest, itrial) += basis[itest] * basis[itrial] * jdet * igpt.getWeight();
//                
//
//                // force vector assembly
//                f_local(itest) += basis[itest] * funcval * igpt.getWeight() * jdet;
//            }
//        }
//        std::cout << k_local << "\n";
        
        K.InsertGlobalValues(scatter, k_local);
//        f.SumIntoGlobalValues(scatter, f_local);

    }
    
    K.GlobalAssemble();
    f.GlobalAssemble();

    if(0 == Comm.MyPID())
    {
        std::cout << K << "\n";
        std::cout << f << "\n";
    }
    
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    
    return EXIT_SUCCESS;
}