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

int main(int argc, char* argv[])
{
    std::ifstream ifs(argv[1]);
    if(!ifs)
        trinurbs::error("Cannot open file for reading\n");
    
    trinurbs::Geometry g;
    if(!g.loadHBSFile(ifs))
        trinurbs::error("Can't open geometry file");
    
    trinurbs::Forest forest(g);
    forest.hrefine(1);
    
    return EXIT_SUCCESS;
}