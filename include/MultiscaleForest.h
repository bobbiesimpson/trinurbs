#ifndef TRINURBS_MULTISCALE_FOREST_H
#define TRINURBS_MULTISCALE_FOREST_H

#include "base.h"
#include "Forest.h"
#include "PeriodicForest.h"
#include "Geometry.h"
#include "MultiscaleBezierNodalElement.h"

#include <map>
#include <memory>

namespace trinurbs {

    /// Forward declarations
    class Forest;
    class Geometry;
    
    class MultiscaleForest
    {
        
    public:
        
        /// Default constructor
        //MultiscaleForest() {}
        
        /// Construct with given geometries
        MultiscaleForest(Geometry& macrogeom,
                         Geometry& microgeom)
        :
            mMacroForest(macrogeom),
            mMicroForest(microgeom)
        {
            init();
        }
        
        // hrefine macro forest
        void hrefineMacro(const uint refine)
        {
            if(0 == refine) return;
            
            macroForest().hrefine(refine);
            init();
        }
        
        /// h refine micro forest
        void hrefineMicro(const uint refine)
        {
            if(0 == refine) return;
            
            microForest().hrefine(refine);
            init();
        }
        
//        /// Macroscale geometry setter
//        void setMacroscaleGeometry(const Geometry& g)
//        { mMacroForest = Forest(g); }
//        
//        /// Microscale geometry setter
//        void setMicroscaleGeometry(const Geometry& g)
//        {
//            mMicroForest = Forest(g);
//        }
        
        /// Default destructor
        virtual ~MultiscaleForest() = default;
        
        /// Default copy constructor
        MultiscaleForest(const MultiscaleForest& f) = default;
        
        /// Default copy assignment
        MultiscaleForest& operator=(const MultiscaleForest& f) = default;
        
        /// Default Move constructor
        MultiscaleForest(MultiscaleForest&& f) = default;
        
        /// Default move assignment
        MultiscaleForest& operator=(MultiscaleForest&& f) = default;
        
        /// Get the multiscale element. We never store it explicitly in this class, but
        /// return just a pointer which the caller can manage.
        std::unique_ptr<MultiscaleBezierNodalElement> multiscaleBezierElement(const uint i) const
        {
            const auto find = mGlobalToLocalPairElMap.find(i);
            if(find == mGlobalToLocalPairElMap.end())
                error("Bad multiscale element index requested.");
            
            return make_unique<MultiscaleBezierNodalElement>(macroForest().bezierElement(std::get<0>(find->second)),
                                                             microForest().bezierElement(std::get<1>(find->second)));
        }
        
        /// Total number of elements for entire microscale geometry
        /// (i.e. # elements microscale * # elements macroscale)
        size_t elemN() const;
        
        /// macro forest getter
        const Forest& macroForest() const
        {
            return mMacroForest;
        }
        
        /// microforest getter
        const Forest& microForest() const
        {
            return mMicroForest;
        }
        
        /// non const macro forest getter
        Forest& macroForest()
        {
            return mMacroForest;
        }
        
        /// non const microforest getter
        Forest& microForest()
        {
            return mMicroForest;
        }
        
    protected:
        
    private:
        
        /// initiate data structures
        void init();
        
        /// clear data structures. Doesn't affect forest structures
        void clearData()
        {
            mGlobalToLocalPairElMap.clear();
            mPeriodicCellNodalConnectivity.clear();
        }
        
        /// Parameterisation of macroscale geometry
        Forest mMacroForest;
        
        /// Parameterisation of microscale geometry
        PeriodicForest mMicroForest;
        
        /// Nodal connectivity from periodic (global) cell node
        /// index to global index.  We can then easily
        /// retrieve the mapping from a micro element node
        /// index to the global node index.
        std::map<uint, std::vector<uint>> mPeriodicCellNodalConnectivity;
        
        /// map from global element index to (macro, micro) element pairing
        std::map<uint, std::tuple<uint, uint>> mGlobalToLocalPairElMap;
    };
}
#endif