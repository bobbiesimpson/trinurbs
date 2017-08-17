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

#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <vector>

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> Triplet;

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
    
    trinurbs::Forest forest(g);
    forest.hrefine(refine);
    
    std::cout << "forest with " << forest.elemN() << " elements\n";
    
    double v = 0.0;
    for(uint ielem = 0; ielem < forest.elemN(); ++ielem)
    {
        const auto bel = forest.bezierElement(ielem);
        
        for(trinurbs::IElemIntegrate igpt(bel->integrationOrder()); !igpt.isDone(); ++igpt)
        {
            const auto gpt = igpt.get();
            const auto basis = bel->basis(gpt.xi, gpt.eta, gpt.zeta);
            const auto ders = bel->basisDers(gpt.xi, gpt.eta, gpt.zeta, trinurbs::DU);
            
            const auto jdet = bel->jacDet(gpt);
//            std::cout << bel->eval(gpt) << "\n";
            //std::cout << bel->tangent(gpt.xi, gpt.eta, gpt.zeta, trinurbs::W) << "\n";
//            std::cout << jdet << "\n";
            
            v += jdet * igpt.getWeight();
        }
    }
    
    std::cout << "volume = " << v << "\n";
    
    
    
//    forest.hrefine(1);
    
    return EXIT_SUCCESS;
}