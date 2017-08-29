#include "base.h"
#include "Geometry.h"
#include "PeriodicForest.h"

#include <iostream>


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
    //pforest.hrefine(3);
    
//    uint volume_count = 0;
//    
//    for(uint iel = 0; iel < pforest.elemN(); ++iel)
//    {
//        std::cout << "ELEMENT: " << iel << "\n\n";
//        const auto el = pforest.bezierElement(iel);
//        bool surface = false;
//        
//        for(uint iface = 0; iface < NFACES; ++iface)
//        {
//            if(el->liesOnParentElFace(faceType(iface)))
//            {
//                surface = true;
//                std::cout << "Element lies on geometry face: " << iface << "\n";
//            }
//            
//        }
//        
//        for(uint iedge = 0; iedge < NEDGES; ++iedge)
//        {
//            if(el->liesOnParentElEdge(edgeType(iedge)))
//            {
//                surface = true;
//                std::cout << "Element lies on geometry edge: " << iedge << "\n";
//            }
//        }
//        
//        for(uint ivertex = 0; ivertex < NVERTICES; ++ivertex)
//        {
//            if(el->liesOnParentElVertex(vertexType(ivertex)))
//            {
//                surface = true;
//                std::cout << "Element lies on geometry vertex: " << ivertex << "\n";
//            }
//        }
//        
//        if(!surface)
//            ++volume_count;
//        std::cout << "\n\n";
//    }
//    
//    std::cout << "There are " << volume_count << " interior elements\n";
//    
    
    return EXIT_SUCCESS;
}
    