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
    
    std::vector<std::vector<uint>> m{{1,2,3,4}, {5,6,7,8}};


    std::vector<uint> igvec{ 12, 15, 0, 56};
    
    auto result = std::min_element(std::begin(igvec), std::end(igvec));
    auto index = std::distance(std::begin(igvec), result);
    
    // perform cyclic rotation such that index defining origin
    // is the first element
    std::rotate(igvec.begin(), igvec.begin() + index, igvec.end());
    
    std::cout << igvec << "\n";
    
    return EXIT_SUCCESS;
}