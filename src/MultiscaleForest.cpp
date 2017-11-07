#include "MultiscaleForest.h"
#include "MultiscaleBezierNodalElement.h"

#include <unordered_map>

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
        
        const auto& macro_forest = macroForest();
        const auto& micro_forest = microForest();
        
        // Cache of previously stored global connectivity
        // (Possibly dangerous since we are using essentially a vector of doubles
        //  as a key for searching. When is a Point3D equal/not equal?)
        std::map<Point3D, uint> pt_cache;
        
        // Temporary global connectivity map for each macro element
        std::map<uint, std::vector<uint>> temp_gconn_map;
        
        uint current_gdof = 0;
        
        for(uint imacro = 0; imacro < macro_forest.elemN(); ++imacro)
        {
            std::cout << "Element: " << imacro << "\n\n";
            // Macro element
            const auto macro_el = macro_forest.bezierElement(imacro);
            
            // Connectivity array for this macro element which maps
            // a (global) periodic node index to a global (macro) node index
            std::vector<uint> pconn(micro_forest.globalDofN());
            
            // loop over all surface (Greville) points of micro element
            for(const auto& local_index : micro_forest.surfaceGrevilleISet())
            {
                // Get the Greville point in (micro) physical space
                const Point3D gpt = micro_forest.grevillePt(local_index);
                
                //std::cout << gpt << "\n";
                
                // And map to physical space using the macro element
                const Point3D x = macro_el->eval(gpt[0], gpt[1], gpt[2]);
                
                std::cout << x << "\n";
                
                auto find = pt_cache.find(x);
                if(find != pt_cache.end())
                {
                    // Global index has been assigned for this point.
                    // Set the relevant entry for this
                    pconn[local_index] = find->second;
                }
                else
                {
                    pconn[local_index] = current_gdof;
                    pt_cache[x] = current_gdof;
                    ++current_gdof;
                }
            }
            
            // And loop over remaining interior points which we know are not
            // connected to any other macro elements
            for(const auto& local_index : micro_forest.interiorGrevilleISet())
            {
                pconn[local_index] = current_gdof;
                ++current_gdof;
            }
            
            // Assign global connectivity for this macro element
            temp_gconn_map[imacro] = pconn;
            
        }
        
        // And now set member data to fully complete global connectivity map
        mPeriodicCellNodalConnectivity = temp_gconn_map;
        
        
        

        
        
    }
    
}