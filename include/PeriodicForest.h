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
        
    protected:
        
    private:
        
        /// Typedef of map which maps a macro face permutations to a vector of
        /// local elements (with local face indices and face permutations)
        typedef std::map<FacePermutation, std::vector<std::tuple<uint, Face, FacePermutation>>> FacePermutationMap;
        
        /// Initialise data structures (after construction or refinement)
        void initGeometry();
        
        /// Vector of maps of elements that lie on each macro face with
        /// associcated permutations.  The size of this vector must equal 6.
        std::vector<FacePermutationMap> mFaceElMaps;
        
        /// Vectors of (geometry element, local face) pairs for each
        /// face of the periodic cell.  These are elements that intersect
        /// the periodic cell domain and are independent of refinement.
        std::vector<std::vector<std::pair<const NAnalysisElement*, Face>>> mFaceGeomEls;
        
        
    };
    
    /// Do all the given control points lie on the given face
    /// assuming a cell domain of [-1,1]^3?
    bool allCPtsLieOnPeriodicCellFace(const std::vector<Point3D>& cpts,
                                      const Face f);
    
    
}
#endif