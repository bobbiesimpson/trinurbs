#include <iostream>
#include <cstdlib>
#include "base.h"
#include "InputDataStructures.h"
#include "BSplineSpace.h"
#include "NURBSCommon.h"
#include "Point.h"
#include "Point4D.h"
#include "Geometry.h"
#include "IElemIntegrate.h"
#include "IParentSample.h"
#include "OutputVTK.h"

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

int main(int argc, char* argv[])
{
    std::ifstream ifs(argv[1]);
    if(!ifs)
        trinurbs::error("Cannot open file for reading\n");
    
    trinurbs::Geometry g;
    if(!g.loadHBSFile(ifs))
        trinurbs::error("Can't open geometry file");
    
    uint refine = 0;
    const uint max_refine = 10;
    if(argc > 2) {
        auto input = std::atoi(argv[2]);
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
    }

    // create forest and refine if necessary
    trinurbs::Forest forest(g);
    forest.hrefine(refine);
    std::cout << "Created forest with " << forest.elemN() << " elements\n";
    
    // now do a simple L2 projection
    
    // first set up data structures
    //
    double v = 0.0;                         // volume
    std::vector<Triplet> coefficients;      // coefficients for sparse matrix
    const uint ndof = forest.globalDofN();
    SpMat K(ndof, ndof);                    // K matrix
    Eigen::VectorXd f(ndof);                // force vector
    
    for(uint i = 0; i < ndof; ++i)
        f(i) = 0.0;
    
    for(uint ielem = 0; ielem < forest.elemN(); ++ielem)
    {
        const auto bel = forest.bezierElement(ielem);
        const auto conn = bel->globalBasisIVec();
        
//        for(const auto& i : conn)
//            std::cout << i << "\t";
//        std::cout << "\n";
        
        // local stiffness matrix
        std::vector<std::vector<double>> submatrix;
        for(size_t i = 0; i < conn.size(); ++i)
            submatrix.push_back(std::vector<double>(conn.size(), 0.0));
        
        for(trinurbs::IElemIntegrate igpt(bel->integrationOrder()); !igpt.isDone(); ++igpt)
        {
            const auto gpt = igpt.get();
            const auto basis = bel->basis(gpt.xi, gpt.eta, gpt.zeta);
            const auto jdet = bel->jacDet(gpt);
            const auto phys_coord = bel->eval(gpt);
            const double funcval = myfunc(phys_coord[0], phys_coord[1], phys_coord[2]);
            
//            std::cout << bel->eval(gpt) << "\n";
            //std::cout << bel->tangent(gpt.xi, gpt.eta, gpt.zeta, trinurbs::W) << "\n";
//            std::cout << jdet << "\n";
            
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
    
    Eigen::SimplicialCholesky<SpMat> chol(K);  // performs a Cholesky factorization of A
    Eigen::VectorXd x_soln = chol.solve(f);
    
    std::vector<double> soln(forest.globalDofN(), 0.0);
    for(size_t i = 0; i < soln.size(); ++i)
        soln[i] = x_soln(i);
    
    std::cout << "solution vector:\n\n" << soln << "\n";
    
    const uint ngridpts = 30;
    trinurbs::OutputVTK output("trivariate", ngridpts);
    output.outputNodalField(forest, "trial-field", soln);
    
    return EXIT_SUCCESS;
}