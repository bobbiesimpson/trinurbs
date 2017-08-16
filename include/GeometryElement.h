#ifndef TRINURBS_GEOMETRY_ELEMENT_H
#define TRINURBS_GEOMETRY_ELEMENT_H

#include <vector>
#include <cassert>
#include <algorithm>
#include <map>
#include <utility>
#include <mutex>
#include <limits>
#include <tuple>
#include <boost/icl/continuous_interval.hpp>

#include "base.h"
#include "Point3D.h"
#include "IElemIntegrate.h"

namespace trinurbs
{
    /// The base element class. Can be instantianted (non-abstract class)
    /// Represents a geometry element providing the neccessary interface
    /// for derived classes.
    
    class Geometry;
    
    class GeometryElement {
        
    public:
        
        /// Construct with a geometry object, space index
        /// and knot coordinates which define the element bounds
        GeometryElement(const Geometry& g,
                        const uint ispace,
                        const DoublePairVec& knots);
        
        /// Get the const geometry object
        const Geometry& geometry() const {return mGeom;}
        
        /// Space index accessor
        uint spaceI() const { return mSpaceI; }
        
        /// Evaluate the physical point given a parent coordinate.
        virtual Point3D eval(const double xi,
                             const double eta,
                             const double zeta) const;
        
        /// Same as above but passing a 2d guass point
        Point3D eval(const GPt3D& gp) const { return eval(gp.xi, gp.eta, gp.zeta); }
        
        /// Get the physical coordinate given a vertex in the parent space
        Point3D evalVertex(const Vertex v) const;
        
        /// Get jacobian determinant from parent to physical space
        virtual double jacDet(const double xi,
                              const double eta,
                              const double zeta) const;
        
        /// Same as above but passing a 2d guass point
        double jacDet(const GPt3D& gp) const {return jacDet(gp.xi, gp.eta, gp.zeta); }
        
        /// Get jacobian determinant from parent to physical space
        virtual double jacDet(const double xi,
                              const double eta,
                              const double zeta,
                              const Point3D& t1,
                              const Point3D& t2) const
        {
            return cross(t1,t2).length() * jacDetParam(xi,eta,zeta);
        }
        
        /// Get jacobian determinant from parent to physical space
        virtual double jacDet(const GPt3D& gp,
                              const Point3D& t1,
                              const Point3D& t2) const
        {
            return cross(t1,t2).length() * jacDetParam(gp.xi, gp.eta, gp.zeta);
        }
        
        /// Get jacobian
        virtual DoubleVecVec jacob(const double xi,
                                   const double eta,
                                   const double zeta) const;
        
        /// Evaluate jacobian with given tangent vectors
        virtual DoubleVecVec jacob(const double xi,
                                   const double eta,
                                   const double zeta,
                                   const Point3D& t1,
                                   const Point3D& t2) const
        {
            throw std::runtime_error("Jacob function not implemented.");
        }
        
        /// Same as above but passing a 2d guass point
        DoubleVecVec jacob(const GPt3D& gp) const {return jacob(gp.xi, gp.eta, gp.zeta); }
        
        /// Evaluate jacobian inverse
        virtual DoubleVecVec jacobInv(const double xi,
                                      const double eta,
                                      const double zeta,
                                      const Point3D& t1,
                                      const Point3D& t2) const
        {
            throw std::runtime_error("Jacob inverse function not implemented.");
        }
        
        /// TODO: we will need to generare the normal for faces where loads
        /// are applied.
        
        /// Get the normal
//        virtual Point3D normal(const double xi,
//                               const double eta,
//                               const double zeta) const;
        
        /// Same as above but passing a 2d guass point
//        Point3D normal(const GPt3D& gp) const {return normal(gp.xi, gp.eta, gp.zeta); }
        
        /// Get the tangent
        virtual Point3D tangent(const double xi,
                                const double eta,
                                const double zeta,
                                const ParamDir d) const;
        
        /// Same as above but passing a 2d guass point
        Point3D tangent(const GPt3D& gp, const ParamDir d) const {return tangent(gp.xi, gp.eta, gp.zeta, d); }
        
        /// Get the degree for all parametric directions.
        UIntVec geometryDegree() const;
        
        /// Wrapper for other function
        std::pair<bool, GPt3D> containsParamPt(const ParamCoord& p) const
        {
            return containsParamPt(p.u, p.v, p.w);
        }
        
        /// Does this element contain the given parametric point? If so, return true
        /// with local parent coordinate
        std::pair<bool, GPt3D> containsParamPt(const double u,
                                               const double v,
                                               const double w) const;
        
        /// Is this element a subelement of this element
        bool contains(const GeometryElement& e) const;
        
        /// Get approximate element size
//        double size() const
//        {
//            if(!mSize.first)
//            {
//                std::lock_guard<std::mutex> lock(*mMutex);
//                const double size = std::max(dist(eval(-1.0, -1.0), eval(1.0, 1.0)),
//                                             dist(eval(1.0, -1.0), eval(-1.0, 1.0)));
//                mSize = std::make_pair(true, size);
//                return size;
//            }
//            else
//                return mSize.second;
//            
//        }
        
//        double aspectRatio() const
//        {
//            return dist(eval(-1.0, 0.0), eval(1.0, 0.0)) / dist(eval(0.0, -1.0), eval(0.0, 1.0));
//        }
        
        /// Get the approximate min coordinate of this element
        /// based on the four vertex points
//        Point3D approxLowerBound() const
//        {
//            const double max = std::numeric_limits<double>::max();
//            Point3D currentLB(max, max, max);
//            std::vector<std::pair<double, double>> vpts
//            {
//                {-1.0,-1.0},
//                {1.0,-1.0},
//                {1.0,1.0},
//                {-1.0,1.0},
//                {0.0, 0.0},
//                {0.0, -1.0},
//                {1.0, 0.0},
//                {0.0, 1.0},
//                {-1.0, 0.0}
//            };
//            for(const auto& p : vpts)
//                currentLB = min(currentLB, eval(p.first, p.second));
//            return currentLB;
//        }
        
        /// Get the approximate max coordinate of this element
        /// based on the four vertex points
//        Point3D approxUpperBound() const
//        {
//            const double min = std::numeric_limits<double>::lowest();
//            Point3D currentUB(min, min, min);
//            std::vector<std::pair<double, double>> vpts
//            {
//                {-1.0,-1.0},
//                {1.0,-1.0},
//                {1.0,1.0},
//                {-1.0,1.0},
//                {0.0, 0.0},
//                {0.0, -1.0},
//                {1.0, 0.0},
//                {0.0, 1.0},
//                {-1.0, 0.0}
//            };
//            for(const auto& p : vpts)
//                currentUB = max(currentUB, eval(p.first, p.second));
//            return currentUB;
//        }
        
        /// parametric lower bound
        ParamCoord lowerBound() const
        {
            return ParamCoord(lowerBound(U), lowerBound(V), lowerBound(W));
        }
        
        /// parametric upper bound
        ParamCoord upperBound() const
        {
            return ParamCoord(upperBound(U), upperBound(V), upperBound(W));
        }
        
        /// Print function that may be overridden as required.
        virtual void print(std::ostream& ost) const;
        
        /// Jacobian from parent space to parametric space
        DoubleVecVec jacobParam(const double xi,
                                const double eta,
                                const double zeta) const
        {
            return {
                { (upperBound(U) - lowerBound(U)) / 2.0, 0.0, 0.0},
                { 0.0, (upperBound(V) - lowerBound(V)) / 2.0, 0.0},
                { 0.0, 0.0, (upperBound(W) - lowerBound(W))}
            };
        }
        
        /// Jacobian determinant for mapping from parent to parametric
        /// space
        double jacDetParam(const double xi,
                           const double eta,
                           const double zeta) const
        {
            return (upperBound(U) - lowerBound(U)) / 2.0
            * (upperBound(V) - lowerBound(V)) / 2.0
            * (upperBound(W) - lowerBound(W)) / 2.0;
        }
        
        /// Transform parent coordinate to parametric coordinate
        ParamCoord getParamCoord(const double xi,
                                 const double eta,
                                 const double zeta) const
        {
            ParamCoord p;
            p.u = xi * (upperBound(U) - lowerBound(U)) / 2.0 + (upperBound(U) + lowerBound(U)) / 2.0;
            p.v = eta * (upperBound(V) - lowerBound(V)) / 2.0 + (upperBound(V) + lowerBound(V)) / 2.0;
            p.w = zeta * (upperBound(W) - lowerBound(W)) / 2.0 + (upperBound(W) + lowerBound(W)) / 2.0;
            return p;
        }
        
        /// Wrapper for above function
        ParamCoord getParamCoord(const GPt3D& g) const
        {
            return getParamCoord(g.xi, g.eta, g.zeta);
        }
        
        /// Get the parent coordinate given a paramatric coordinate
        GPt3D getParentCoord(const double u,
                             const double v,
                             const double w) const
        {
            return GPt3D((2.0 * u - (upperBound(U) + lowerBound(U))) / (upperBound(U) - lowerBound(U)),
                         (2.0 * v - (upperBound(V) + lowerBound(V))) / (upperBound(V) - lowerBound(V)),
                         (2.0 * w - (upperBound(W) + lowerBound(W))) / (upperBound(W) - lowerBound(W)));
        }
        
        /// Get the parent coordinate given a paramatric coordinate
        GPt3D parentCoord(const ParamCoord& paramcoord) const
        {
            return getParentCoord(paramcoord.u, paramcoord.v, paramcoord.w);
        }
        
        /// Get lower bound (knot coordinate) for specified parametric direction
        double lowerBound(const ParamDir d) const
        {
            assert(d < mKnotIntervals.size());
            return mKnotIntervals[d].lower();
        }
        
        /// Get upper bound (knot coordinate) for specified parametric direction
        double upperBound(const ParamDir d) const
        {
            assert(d < mKnotIntervals.size());
            return mKnotIntervals[d].upper();
        }
        
        /// Get the parametric interval for this element for the given parametric direction
        const boost::icl::continuous_interval<double>& boostInterval(const ParamDir dir) const
        {
            assert(mKnotIntervals.size() > 1);
            return mKnotIntervals[dir];
        }
        
        /// Does this element contain any degenerate edges?
//        bool degenerate() const
//        {
//            const double tol = 1.e-9;
//            std::vector<Point3D> vcoords
//            {
//                evalVertex(Vertex::VERTEX0),
//                evalVertex(Vertex::VERTEX1),
//                evalVertex(Vertex::VERTEX2),
//                evalVertex(Vertex::VERTEX3)
//            };
//            for(size_t ivertex = 0; ivertex < vcoords.size(); ++ivertex)
//            {
//                for(size_t iinner = 0; iinner < vcoords.size(); ++iinner)
//                {
//                    if(ivertex == iinner)
//                        continue;
//                    if(dist(vcoords[ivertex], vcoords[iinner]) < tol)
//                        return true;
//                }
//            }
//            return false;
//        }
        
//        std::pair<bool, Edge> degenerateEdge() const
//        {
//            if(!degenerate())
//                return std::make_pair(false, Edge::EDGE0);
//            
//            const double tol = 1.e-9;
//            std::vector<Point3D> vcoords
//            {
//                evalVertex(Vertex::VERTEX0),
//                evalVertex(Vertex::VERTEX1),
//                evalVertex(Vertex::VERTEX2),
//                evalVertex(Vertex::VERTEX3)
//            };
//            
//            if(dist(vcoords[0], vcoords[1])< tol)
//                return std::make_pair(true, Edge::EDGE0);
//            else if(dist(vcoords[2], vcoords[3])< tol)
//                return std::make_pair(true, Edge::EDGE1);
//            else if(dist(vcoords[1], vcoords[3])< tol)
//                return std::make_pair(true, Edge::EDGE3);
//            else if(dist(vcoords[2], vcoords[0])< tol)
//                return std::make_pair(true, Edge::EDGE2);
//            else
//            {
//                throw std::runtime_error("Bad degenerate edge call.");
//                return std::make_pair(false, Edge::EDGE0);
//            }
//        }
        
    protected:
        
    private:
        
        /// Reference to the forest where this element belongs. Non-owning pointer
        const Geometry& mGeom;
        
        /// The space that this element lies in
        const uint mSpaceI;
        
        /// Mutex to protected cached data
        std::shared_ptr<std::mutex> mMutex;
        
        void insertCachedPoint(const std::tuple<double, double, double>& p, const Point3D& pt) const
        {
            std::lock_guard<std::mutex> lock(*mMutex);
            mPointCache.insert(std::make_pair(p,pt));
        }
        
        void insertCachedJDet(const std::tuple<double, double, double>& p, const double jdet) const
        {
            std::lock_guard<std::mutex> lock(*mMutex);
            mJDetCache.insert(std::make_pair(p,jdet));
        }
        
        void insertCachedJacob(const std::tuple<double, double, double>& p, const DoubleVecVec& jacob) const
        {
            std::lock_guard<std::mutex> lock(*mMutex);
            mJacobCache.insert(std::make_pair(p,jacob));
        }
        
        /// The non-zero knot interval that defines this element
        std::vector<boost::icl::continuous_interval<double>> mKnotIntervals;
        
        /// Cache evaluation of physical points
        mutable std::map<std::tuple<double, double, double>, Point3D> mPointCache;
        
        /// Cache jacobian determinant
        mutable std::map<std::tuple<double, double, double>, double> mJDetCache;
        
        /// Cache jacobian matrix
        mutable std::map<std::tuple<double, double, double>, DoubleVecVec> mJacobCache;
        
        /// Size cache
        mutable std::pair<bool, double> mSize;
        
        /// Overload output operator
        friend std::ostream& operator<<(std::ostream& ost, const GeometryElement& e)
        { e.print(ost); return ost; }
        
    };
    
    /// Get distance between two elements
    double dist(const GeometryElement& e1, const GeometryElement& e);
    
    
}

#endif