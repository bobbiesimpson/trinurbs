#include "MultiscaleForest.h"
#include "MultiscaleBezierNodalElement.h"

namespace trinurbs {
    
    size_t MultiscaleForest::elemN() const
    {
        return mGlobalToLocalPairElMap.size();
    }
    
    void MultiscaleForest::init()
    {
        clearData();
        
        // construct mapping from global (micro) element index to its (macro,micro) index pair.
        uint icurrent = 0;
        for(uint imacro = 0; imacro < macroForest().elemN(); ++imacro)
            for(uint imicro = 0; imicro < microForest().elemN(); ++imicro)
            {
                mGlobalToLocalPairElMap.insert(std::make_pair(icurrent, std::make_tuple(imacro, imicro)));
                ++icurrent;
            }
        
        // construct global nodal connectivty
        
        // general idea is to loop over each macro element and search each of the
        // neighbouring macro elements for previusly assigned node indices.
        //
        // We do a search over neighbouring 'geometry' elements at the micro scale
        // making sure we check any elements connected to vertices, edge or faces
        // of the periodic cell domain.  We use the present macro element as a reference
        // for determining the order in which we traverse node indices over elements
        // in the micro scale.
        
        
        
        
        
        
        
        // do some preprocessing of period (micro) cell connectivity
//        const auto micro_forest = microForest();
//        
//        // loop over micro elements
//        for(uint iel = 0; iel < micro_forest.elemN(); ++iel)
//        {
//            for(
//        }
        
        
        
    }
    
}