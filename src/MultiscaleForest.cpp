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
        
        uint icurrent = 0;
        for(uint imacro = 0; imacro < macroForest().elemN(); ++imacro)
            for(uint imicro = 0; imicro < microForest().elemN(); ++imicro)
            {
                mGlobalToLocalPairElMap.insert(std::make_pair(icurrent, std::make_tuple(imacro, imicro)));
                ++icurrent;
            }
    }
    
}