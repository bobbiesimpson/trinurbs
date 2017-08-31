#ifndef TRINURBS_PERIODIC_FOREST_H
#define TRINURBS_PERIODIC_FOREST_H

#include "Forest.h"
#include "base.h"
#include "Geometry.h"

#include <map>
#include <vector>
#include <tuple>

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
            initGeometry();
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
        
    protected:
        
    private:

        
        typedef std::vector<std::pair<const NAnalysisElement*, Face>> ElemFacePairVec;
        
        typedef std::vector<std::pair<const NAnalysisElement*, Edge>> ElemEdgePairVec;
        
        typedef std::vector<std::pair<const NAnalysisElement*, Vertex>> ElemVertexPairVec;

        
        /// Initialise data structures (after construction or refinement)
        void initGeometry();
        
        /// Initialise analysis data structures (called after refinement).
        void initAnalysisData();
        
        /// Clear all geometry data
        void clearGeometryData()
        {
            mFaceGeomEls.clear();
            mEdgeGeomEls.clear();
            mVertexGeomEls.clear();
        }

        /// Vectors of (geometry element, local face) pairs for each
        /// face of the periodic cell.  These are elements that intersect
        /// the periodic cell domain and are independent of refinement.
        std::map<Face, ElemFacePairVec> mFaceGeomEls;
        
        /// Vector of geometry elements connected to each edge of the
        /// periodic cell.
        std::map<Edge, ElemEdgePairVec> mEdgeGeomEls;
        
        /// Vector of geometry elements connected to each vertex of the
        /// periodic cell.
        std::map<Vertex, ElemVertexPairVec> mVertexGeomEls;
        
        
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