#include <iostream>

#include "Forest.h"
#include "Geometry.h"
#include "MultiscaleForest.h"
#include "base.h"
#include "OutputVTK.h"

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
    
    std::cout << "Writing output vtk files.....\n";
    const uint ngridpts = 2;
    OutputVTK output("trivariate_ms", ngridpts);
    //output.outputForestGeometry(multiscaleforest.macroForest());
    output.outputForestGeometry(multiscaleforest.microForest());
    
    output.outputMultiscaleForestGeometry(multiscaleforest);
    
    
    return EXIT_SUCCESS;
}