#include <iostream>

#include "Forest.h"
#include "Geometry.h"
#include "MultiscaleForest.h"
#include "base.h"
#include "OutputVTK.h"
#include "IElemIntegrate.h"
#include "PeriodicForest.h"

#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <vector>

using namespace trinurbs;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> Triplet;

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
    
    double v = 0.0;
    
    std::cout << "Now performing simple L2 projection test....\n";
    
    std::vector<Triplet> coefficients;      // coefficients for sparse matrix
    const uint ndof = multiscaleforest.globalDofN();
    SpMat K(ndof, ndof);                    // K matrix
    Eigen::VectorXd f(ndof);                // force vector
    
    for(uint i = 0; i < ndof; ++i)
        f(i) = 0.0;
    
    for(uint ielem = 0; ielem < multiscaleforest.elemN(); ++ielem)
    {
        const auto bel = multiscaleforest.multiscaleBezierElement(ielem);
        const auto conn = bel->globalBasisIVec();
        
        // local stiffness matrix
        std::vector<std::vector<double>> submatrix;
        for(size_t i = 0; i < conn.size(); ++i)
            submatrix.push_back(std::vector<double>(conn.size(), 0.0));
        
        for(trinurbs::IElemIntegrate igpt(bel->integrationOrder()); !igpt.isDone(); ++igpt)
        {
            const auto gpt = igpt.get();
            const auto basis = bel->basis(gpt.xi, gpt.eta, gpt.zeta);
            const auto jdet = bel->jacDet(gpt.xi, gpt.eta, gpt.zeta);
            const auto phys_coord = bel->eval(gpt.xi, gpt.eta, gpt.zeta);
            const double funcval = myfunc(phys_coord[0], phys_coord[1], phys_coord[2]);
            
            for(size_t itest = 0; itest < conn.size(); ++itest)
            {
                const auto gtest_i = conn[itest];
                for(size_t itrial = 0; itrial < conn.size(); ++itrial)
                    submatrix[itest][itrial] += basis[itest] * basis[itrial] * jdet * igpt.getWeight();
                
                // force vector assembly
                f(gtest_i) += basis[itest] * funcval * igpt.getWeight() * jdet;
            }
            v += jdet * igpt.getWeight();
        }
        
        for(size_t itest = 0; itest < conn.size(); ++itest)
        {
            const auto gtest_i = conn[itest];
            for(size_t itrial = 0; itrial < conn.size(); ++itrial)
            {
                const auto gtrial_i = conn[itrial];
                //                std::cout << submatrix[itest][itrial] << "\n";
                coefficients.push_back(Triplet(gtest_i, gtrial_i, submatrix[itest][itrial]));
            }
        }
    }
    
    K.setFromTriplets(coefficients.begin(), coefficients.end());
    
//    std::cout << K << "\n";
//    std::cout << f << "\n";
    
    Eigen::SimplicialCholesky<SpMat> chol(K);  // performs a Cholesky factorization of A
    Eigen::VectorXd x_soln = chol.solve(f);
    
    std::cout << "SOLUTION!!!!!\n";
    
    std::cout << x_soln << "\n";
    
    std::vector<double> soln(multiscaleforest.globalDofN(), 0.0);
    for(size_t i = 0; i < soln.size(); ++i)
        soln[i] = x_soln(i);

    std::cout << "volume of multiscale geometry = " << v << "\n";
    
    const uint ngridpts = 2;
    OutputVTK output("multiscale", ngridpts);
    output.outputNodalField(multiscaleforest, "testdata", soln);
    
    output.outputForestGeometry(multiscaleforest.microForest());
    
    std::cout << "done!\n";
    
    return EXIT_SUCCESS;
}