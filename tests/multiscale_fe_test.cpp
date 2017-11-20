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
#include "Epetra_Import.h"

#include "AztecOO.h"
#include "Amesos.h"

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
    
   int ndof = multiscaleforest.globalDofN();
    
    // Map for mpi communication of matrix and vector terms
    Epetra_Map map(ndof, 0, Comm);
    
    // Construct sparse matrix
    Epetra_FECrsMatrix K(Copy, map, 5);
    
    // Construct RHS vector
    Epetra_FEVector f(map);
    
    // Set up element map for MPI assembly routines
    int nel = multiscaleforest.elemN();
    Epetra_Map element_map(nel, 0, Comm);
    
    // Get data structure corresponding to local elements to this process
    int num_local_els = element_map.NumMyElements();
    int* my_elements = element_map.MyGlobalElements();

    if(0 == Comm.MyPID())
        std::cout << "Assembly progress....\n";
    
    std::vector<int> my_element_indices;
    
    const double barwidth = 70.0;
    
    // Now loop over 'local' elements
    for(auto ilocal = 0; ilocal < num_local_els; ++ilocal)
    {
        const double progress = static_cast<double>(ilocal) / num_local_els;
        
        if(0 == Comm.MyPID())
        {
            std::cout << "[";
            int pos = barwidth * progress;
            for (int i = 0; i < barwidth; ++i) {
                if (i < pos) std::cout << "=";
                else if (i == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << int(progress * 100.0) << " %\r";
            std::cout.flush();
        }
        
        // global element index
        auto iel = my_elements[ilocal];
        
        const auto bel = multiscaleforest.multiscaleBezierElement(iel);
        
        const auto unsigned_conn = bel->globalBasisIVec();
        
        std::vector<int> conn;
        for(const auto& uval : unsigned_conn)
            conn.push_back(static_cast<int>(uval));
        
        // add element indices to global vector
        my_element_indices.insert(my_element_indices.end(), conn.begin(), conn.end());
        
        // local stiffness matrix terms
        Epetra_SerialDenseMatrix k_local(conn.size(), conn.size());
        
        // local force vector terms
        Epetra_SerialDenseVector f_local(conn.size());
        Epetra_IntSerialDenseVector scatter(Copy, conn.data(), conn.size());
        
        // Perform quadrature
        for(trinurbs::IElemIntegrate igpt(bel->integrationOrder()); !igpt.isDone(); ++igpt)
        {
            const auto gpt = igpt.get();
            const auto basis = bel->basis(gpt.xi, gpt.eta, gpt.zeta);
            const auto jdet = bel->jacDet(gpt.xi, gpt.eta, gpt.zeta);
            const auto phys_coord = bel->eval(gpt.xi, gpt.eta, gpt.zeta);
            const double funcval = myfunc(phys_coord[0], phys_coord[1], phys_coord[2]);
            
            for(int itest = 0; itest < conn.size(); ++itest)
            {
                for(int itrial = 0; itrial < conn.size(); ++itrial)
                    k_local(itest, itrial) += basis[itest] * basis[itrial] * jdet * igpt.getWeight();

                // force vector assembly
                f_local(itest) += basis[itest] * funcval * igpt.getWeight() * jdet;
            }
        }
        
        K.InsertGlobalValues(scatter, k_local);
        f.SumIntoGlobalValues(scatter, f_local);

    }
    
    K.GlobalAssemble();
    f.GlobalAssemble();
    
    if(0 == Comm.MyPID())
        std::cout << "Starting solver....\n";

    // Now generate the solution to Kx=f
    Epetra_Vector x(map);

    Epetra_LinearProblem problem(&K, &x, &f);
    
    // Use an interative CG solver
    AztecOO solver(problem);

    solver.SetAztecOption(AZ_precond, AZ_Jacobi);
    solver.Iterate(100, 1.0E-6);
    
    //
    // Output
    //
    
    // First determine the set of node indices for the elements local to this process
    std::sort(my_element_indices.begin(), my_element_indices.end());
    auto last = std::unique(my_element_indices.begin(), my_element_indices.end());
    my_element_indices.erase(last, my_element_indices.end());
    
    // And now generate a map for this set
    Epetra_Map target_map(-1, my_element_indices.size(), my_element_indices.data(), 0, Comm);
    
    Epetra_Import import(target_map, map);
    
    // Create a vector with this new local set and populate it with the relevant values
    Epetra_Vector output_vec(target_map);
    output_vec.Import(x, import, Epetra_CombineMode::Insert);
    
    // Create a local map (used by OutputVTK class)
    std::map<int, double> local_soln_map;
    for(int i = 0; i < target_map.NumMyElements(); ++i)
        local_soln_map[target_map.GID(i)] = output_vec[i];
    
    const uint ngridpts = 2;
    std::string fname = "multiscale" + std::to_string(Comm.MyPID());
    OutputVTK output(fname, ngridpts);
    output.outputNodalField(multiscaleforest, "testdata", local_soln_map, std::vector<int>(my_elements, my_elements + num_local_els));

    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    
    return EXIT_SUCCESS;
}