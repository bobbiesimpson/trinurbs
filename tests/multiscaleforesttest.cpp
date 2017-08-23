#include <iostream>

#include "Forest.h"
#include "Geometry.h"
#include "MultiscaleForest.h"
#include "base.h"
#include "OutputVTK.h"
#include "IElemIntegrate.h"

using namespace trinurbs;

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
    
    multiscaleforest.hrefineMacro(macro_refine);
    multiscaleforest.hrefineMicro(micro_refine);
    
    std::cout << "Multiscale geometry with "
              << multiscaleforest.macroForest().elemN() << " elements (macro), "
              << multiscaleforest.microForest().elemN() << " elements (micro), "
              << multiscaleforest.elemN() << " elements total\n\n";
    
    double v = 0.0;
    
    for(size_t iel = 0; iel < multiscaleforest.elemN(); ++iel)
    {
        const auto el = multiscaleforest.multiscaleBezierElement(iel);
        for(IElemIntegrate igpt(el->integrationOrder()); !igpt.isDone(); ++igpt)
        {
            const auto gpt = igpt.get();
            v += el->jacDet(gpt.xi, gpt.eta, gpt.zeta) * igpt.getWeight();
        }
    }
    
    std::cout << "volume of multiscale geometry = " << v << "\n";
    
    const uint ngridpts = 2;
    OutputVTK output("micro_", ngridpts);

    std::cout << "outputing micro geometry....\n";
    output.outputForestGeometry(multiscaleforest.microForest());
    
    output.setFilename("macro_");
    output.outputForestGeometry(multiscaleforest.macroForest());
    
    output.setFilename("multiscale_");
    std::cout << "outputting multiscale geometry....\n";
    output.outputMultiscaleForestGeometry(multiscaleforest);
    
    
    std::cout << "done!\n";
    
    return EXIT_SUCCESS;
}