#ifndef TRINURBS_GEOMETRY_H
#define TRINURBS_GEOMETRY_H

#include <iostream>
#include <ostream>
#include <vector>
#include <memory>
#include <cassert>

#include "base.h"
#include "Forest.h"
#include "Point3D.h"
#include "Point4D.h"

namespace trinurbs
{
    /// A class that represents a surface geometry description.
    /// Comprises of a Forest which contains the relevant B-spline spaces
    /// and a set of control points that are interpolated using the given
    /// set of B-spline basis functions.
    
    /// The thinking behind this class is to maintain a forest class hierarchy that contains
    /// no geometry information. However, every forest must maintain a link to a primal
    /// geometry data structure to evaluate quantities such as normals, jacobians, etc.
    /// which is provided by this class.
    
    class Geometry {
        
    public:
        
        /// Default constructor
        Geometry() {}
        
        Geometry(const Forest& f,
                 const std::vector<Point4D>& cpts,
                 bool b = false)
        :
            mPrimalForest(f),
            mCPts(cpts),
            mFlipNormals(b) {}
        
        /// Clear all the data for this geometry object
        void clear()
        {
            mPrimalForest.clear();
            mCPts.empty();
        }
        /// Const accessor for primal forest
        const Forest& primalForest() const
        {return mPrimalForest; }
        
        /// Const accessor of control points
        const Point4D& controlPt(const uint i) const
        {
            assert(i < mCPts.size());
            return mCPts[i];
        }
        
        /// number of control points
        uint controlPtN() const { return mCPts.size(); }
        
        /// Number of B-spline spaces
        uint spaceN() const
        { return primalForest().spaceN(); }
        
        /// Degrees for each parametric diredtion of this space
        UIntVec degree(const uint ispace) const;
        
        /// Jacobian determinant
        double jacDet(const double u,
                      const double v,
                      const double w,
                      const uint ispace) const;
        
        /// Jacobian
        DoubleVecVec jacob(const double u,
                           const double v,
                           const double w,
                           const uint ispace) const;
        
        /// Normal
//        Point3D normal(const double u,
//                       const double v,
//                       const double w,
//                       const uint ispace) const;
        
        /// Tangent
        Point3D tangent(const double u,
                        const double v,
                        const double w,
                        const uint ispace,
                        const ParamDir dir) const;
        
        /// Interpolate the surface
        Point3D eval(const double u,
                     const double v,
                     const double w,
                     const uint ispace) const;
//
//        /// Write to vtu file.
//        void writeVTKOutput(const std::string& file,
//                            const uint nsample = DEFAULT_NGRID_PTS) const;
//        
//        /// Load from a file stream
        bool load(std::istream& ist);

        /// Load from python script output in HBS format.
        bool loadHBSFile(std::istream& ist);

        /// Print to a file stream
        void print(std::ostream& ost) const;
//
//        /// set flipped normals flag
//        void flipNormals(bool b) { mFlipNormals = b; }
//        
//        /// normal flipped getter
//        bool normalsFlipped() const { return mFlipNormals; }
//        
//        /// rescale geometry by a given factor.
//        void rescale(const double sf);
//        
//        /// Rotate the geometry around a given axis by an angle (in radians)
//        void rotate(const nurbs::CartesianComponent comp,
//                    const double rotation);
//        
//        /// Translate geometry by given vector
//        void translate(const nurbs::Point3D& p);
//        
        /// normalise to [-1,1]^3
        void normaliseToParentInterval();
        
    protected:
        
        /// Get point given a set of local indices for each parmetric direction
        const Point4D& point(const uint sp,
                             const uint i,
                             const uint j,
                             const uint k) const
        { return mCPts.at(primalForest().globalI(sp,i,j,k)); }
        
        /// Non-const point accessor
        Point4D& point(const uint i)
        {
            return mCPts[i];
        }
        
        /// Get a point given a 'global' space index. Use a row major numbering system
        const Point4D& point(const uint sp, const uint i) const
        { return mCPts.at(primalForest().globalI(sp,i)); }
        
    private:
        
        /// The forest which is used to discretise the geometry.
        Forest mPrimalForest;
        
        /// The set of control points which are interpolated by the basis functions
        /// defined by the forest.
        std::vector<Point4D> mCPts;
        
        /// flag to indicate if normals are flipped
        bool mFlipNormals;
        
        /// Overload input operator
        friend std::istream& operator>>(std::istream& ist, Geometry& g);
        
        /// Overload outut operator
        friend std::ostream& operator<<(std::ostream& ost, const Geometry& g);
    };
}

#endif