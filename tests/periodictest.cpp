#include "base.h"
#include "Geometry.h"
#include "PeriodicForest.h"
#include "OutputVTK.h"
#include "Point3D.h"

#include <iostream>
#include <map>


using namespace trinurbs;

int main(int argc, char* argv[])
{
    std::ifstream ifs(argv[1]);
    if(!ifs)
        trinurbs::error("Cannot open file for reading\n");
    
    trinurbs::Geometry g;
    if(!g.loadHBSFile(ifs))
        trinurbs::error("Can't open geometry file");
    
    PeriodicForest pforest(g);
    //pforest.hrefine(2);
    
    const uint npoints = 10;
    OutputVTK output("periodic", npoints);
    output.outputForestGeometry(pforest);
    
    
    return EXIT_SUCCESS;
}
    