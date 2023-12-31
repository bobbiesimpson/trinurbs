#ifndef TRINURBS_PERIODIC_FOREST_H
#define TRINURBS_PERIODIC_FOREST_H

#include "Forest.h"
#include "base.h"
#include "Geometry.h"

#include <map>
#include <vector>
#include <tuple>
#include <set>

namespace trinurbs {
    
    
    /// A represnetation of a periodic forest used to represent a unit cell
    /// for multiscale geometry.
    ///
    /// There are several functions will are useful for generating the connectivity
    /// of a multiscale geometry such as determining which elements are
    /// 'opposite' to one another, node indices along edges, faces, vertices etc.
    ///
    class PeriodicForest : public Forest
    {
    public:
        
        /// Default constructor
        PeriodicForest()
        :
            Forest()
        {}
        
        /// Constructor with geometry refernece
        PeriodicForest(Geometry& g)
        :
            Forest(g)
        {
            g.normaliseToParentInterval();
            initAnalysisData();
        }
        
        /// Use default destructor
        virtual ~PeriodicForest() = default;
        
        /// Copy constructor
        PeriodicForest(const PeriodicForest& f)
        :
            Forest(f)
        {}
        
        /// Assignment operator
        PeriodicForest& operator=(const PeriodicForest& f)
        {
            PeriodicForest temp(f);
            *this = std::move(temp);
            return *this;
        }
        
        
        /// Move constructor. Simply move all the data across
        PeriodicForest(PeriodicForest&& f) = default;
        
        /// Move assignment operator. As before, move all the data across.
        PeriodicForest& operator=(PeriodicForest&& f) = default;
        
        /// Override hrefinement
        virtual void hrefine(const uint n) override
        {
            Forest::hrefine(n);
            initAnalysisData();
        }
        
        /// Greville point accessor
        Point3D grevillePt(const uint i) const
        {
            return mGrevilleAbscissaMap.at(i);
        }
        
        /// Total number of Greville abscissa
        size_t grevillePtN() const
        {
            return mGrevilleAbscissaMap.size();
        }
        
        /// Accessor for set of surface Greville
        /// abscissa indices
        const std::set<uint>& surfaceGrevilleISet() const
        {
            return mSurfaceGrevilleISet;
        }
        
        /// Accessor for set of interior Greville
        /// abscissa indices
        const std::set<uint>& interiorGrevilleISet() const
        {
            return mInteriorGrevilleISet;
        }
        
    protected:
        
    private:
        
        /// Initialise analysis data structures (called after refinement).
        void initAnalysisData();
        
        /// Clear all analysis data
        void clearAnalysisData()
        {
            mGrevilleAbscissaMap.clear();
            mSurfaceGrevilleISet.clear();
            mInteriorGrevilleISet.clear();
        }

        /// Vector of greville abscissa for this forest. The order
        /// corresponds directly to the global nodal connectivity
        /// of this forest.
        std::map<uint, Point3D> mGrevilleAbscissaMap;
        
        /// Index of Greville parametric coordinates that lie
        /// on the unit cell surface.
        std::set<uint> mSurfaceGrevilleISet;
        
        /// Index of Greville parametric coordinates that lie
        /// in the cell interior (not on the surface).
        std::set<uint> mInteriorGrevilleISet;
        
    };
    
    /// Do all the given control points lie on the given face
    /// assuming a cell domain of [-1,1]^3?
    bool allCPtsLieOnPeriodicCellFace(const std::vector<Point3D>& cpts,
                                      const Face cellface);
    
    /// Do all the given points lie on the given cell edge?
    bool allCPtsLieOnPeriodicCellEdge(const std::vector<Point3D>& cpts,
                                      const Edge celledge);
    
    // does the given point lie on the specified cell vertex? If yes,
    // return true and the associated cell vertex.
    std::pair<bool, Vertex> liesOnPeriodicCellVertex(const Point3D& p);
    
    
}
#endif