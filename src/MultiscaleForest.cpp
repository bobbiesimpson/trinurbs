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
        
        // construct global nodal connectivty (big!)
        
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