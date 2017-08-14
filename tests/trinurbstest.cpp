#include <iostream>
#include <cstdlib>
#include "base.h"
#include "InputDataStructures.h"
#include "BSplineSpace.h"
#include "NURBSCommon.h"
#include "Point.h"
#include "Point4D.h"
#include "Geometry.h"

using namespace trinurbs;

int main(int argc, char* argv[])
{
    trinurbs::Geometry g;
    trinurbs::Point4D p(0.0, 1.0, 0.0, 1.0);
    
    for(uint iface = 0; iface < NFACES; ++iface)
        std::cout << localBasisIVec(faceType(iface), 3, 3, 3) << "\n";
    
    return EXIT_SUCCESS;
}