#ifndef TRINURBS_NURBS_SPACE_H
#define TRINURBS_NURBS_SPACE_H

#include <vector>
#include <iostream>
#include <ostream>
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <map>
#include <tuple>
#include <boost/icl/continuous_interval.hpp>
#include <boost/icl/closed_interval.hpp>

#include "base.h"
#include "InputDataStructures.h"

namespace trinurbs
{
    
    /// A representation of the set of basis functions that
    /// constitute a b-spline space. We deliberately separate this
    /// from the geometry information since the spaces
    /// require no knowledge of the geometry
    
    class BSplineSpace {
        
    public:
        
        /// Default constructor
        BSplineSpace() = default;
        
        /// Constructor with knot vector and degree vector
        BSplineSpace(const DoubleVecVec& knotvec,
                     const UIntVec& d,
                     std::string n = "no_name")
        :   mKnotVecs(knotvec),
            mDegrees(d),
            mName(n)
        {
            if( knotvec.size() != 3 )
                error( "Only R^3 nurbs spaces are can be used." );
            init();
        }
        
        /// Alternative constructor. Simply delegates work to former constructor
        BSplineSpace(const DoubleVecVec& knotvec,
                     const uint p,
                     const uint q,
                     const uint r,
                     std::string n = "no_name")
        : BSplineSpace(knotvec,{p,q,r}, n) {}
        
        /// Return the total number of basis functions
        uint basisFuncN() const
        {
            uint n = 1;
            for(const auto& v : mNumBasisVec)
                n *= v;
            return n;
        }
        
        /// basis function number getter
        uint basisFuncN(const ParamDir dir) const
        {
            assert(dir < mNumBasisVec.size());
            return mNumBasisVec[dir];
        }
        
        /// Get a vector of degrees
        UIntVec degree() const { return mDegrees; }
        
        /// degree getter
        inline uint degree(const ParamDir dir) const
        {
            assert(dir < mDegrees.size());
            return mDegrees[dir];
        }
        
        /// Const name accessor
        std::string name() const { return mName; }
        
        /// Non-const name accessor
        std::string& name() { return mName; }
        
        /// Get knot vector for specified parametric direction
        inline DoubleVec knotVec(const ParamDir dir) const
        {
            assert(dir < paramDimN());
            return mKnotVecs[dir];
        }
        
        /// knot coordinate getter
        inline double knot(const uint i, const ParamDir dir) const
        {
            assert(dir < mKnotVecs.size());
            assert(i < mKnotVecs[dir].size());
            return mKnotVecs[dir][i];
        }
        
        /// Unique knot vector getter
        std::vector<double> uniqueKnotVec(const ParamDir dir) const
        {
            return mUniqueKnotVecs[dir];
        }
        
        /// unique knot coordinate getter
        inline double uniqueKnot(const uint i, const ParamDir dir) const
        {
            assert(dir < mUniqueKnotVecs.size());
            assert(i < mUniqueKnotVecs[dir].size());
            return mUniqueKnotVecs[dir][i];
        }
        
        /// unique knot number getter
        inline uint uniqueKnotN(const ParamDir dir) const
        {
            assert(dir < mUniqueKnotVecs.size());
            return mUniqueKnotVecs[dir].size();
        }
        
        /// return number of times knot value is repeated
        inline uint knotRepeats(const uint i, const ParamDir dir) const
        {
            assert(i < uniqueKnotN(dir));
            const auto kv = knotVec(dir);
            return std::count(kv.begin(), kv.end(), uniqueKnot(i, dir));
        }
        
        /// Get the range for the given parametric direction
        std::pair<double, double> range(const ParamDir dir) const
        {
            return std::make_pair(mUniqueKnotVecs[dir].front(), mUniqueKnotVecs[dir].back());
        }
        
        /// Greville abscissa not yet implemented
        
//        /// Get a 'global' greville abscissa index given parametric indices
//        uint globalGrevillePtI(const uint i,
//                               const uint j) const
//        {
//            assert(i < grevilleAbscissaPtN(ParamDir::S));
//            assert(j < grevilleAbscissaPtN(ParamDir::T));
//            return j * grevilleAbscissaPtN(ParamDir::S) + i;
//            
//        }
//        /// Number of greville abscissa points in given direction
//        inline uint grevilleAbscissaPtN(const ParamDir dir) const { return basisFuncN(dir); }
//        
//        /// Number of greville abscissa points on this space
//        inline uint grevilleAbscissaPtN() const { return basisFuncN(); }
//        
//        /// Greville abscissa point getter
//        double grevilleAbscissaPt(const uint i, const ParamDir dir) const
//        {
//            assert(i < grevilleAbscissaPtN(dir));
//            double xi = 0.0;
//            for(uint k = i + 1; k < i + degree(dir) + 1; ++k)
//                xi += knot(k, dir) / degree(dir);
//            return xi;
//        }
//        
//        /// Return the greville abscissa vector for a given parametric direction
//        std::vector<double> grevilleAbscissa(const ParamDir dir) const
//        {
//            std::vector<double> v;
//            for(uint i = 0; i < grevilleAbscissaPtN(dir); ++i)
//                v.push_back(grevilleAbscissaPt(i, dir));
//            return v;
//        }
//        
//        /// Get a two dimensionl Greville point
//        GPt2D grevilleAbscissaPt(const uint i) const
//        {
//            assert(i < grevilleAbscissaPtN());
//            const uint sindex = i % grevilleAbscissaPtN(S);
//            const uint tindex = i / grevilleAbscissaPtN(S);
//            return GPt2D(grevilleAbscissaPt(sindex, S),
//                         grevilleAbscissaPt(tindex, T));
//        }
//        
//        /// REturn the Greville abscissa  points along the given edge
//        std::vector<GPt2D> grevilleAbscissaPts(const Edge e) const
//        {
//            std::vector<GPt2D> v;
//            switch (e) {
//                case Edge::EDGE0:
//                    for(uint i = 0; i < grevilleAbscissaPtN(S); ++i)
//                        v.push_back(GPt2D(grevilleAbscissaPt(i, S), range(T).first));
//                    break;
//                case Edge::EDGE1:
//                    for(uint i = 0; i < grevilleAbscissaPtN(S); ++i)
//                        v.push_back(GPt2D(grevilleAbscissaPt(i, S), range(T).second));
//                    break;
//                case Edge::EDGE2:
//                    for(uint j = 0; j < grevilleAbscissaPtN(T); ++j)
//                        v.push_back(GPt2D(range(S).first, grevilleAbscissaPt(j, T)));
//                    break;
//                case Edge::EDGE3:
//                    for(uint j = 0; j < grevilleAbscissaPtN(T); ++j)
//                        v.push_back(GPt2D(range(S).second, grevilleAbscissaPt(j, T)));
//                    break;
//            }
//            return v;
//        }
        
        // Number of non-zero knot spans in a given direction
        inline uint nonzeroKnotSpanN(const ParamDir dir) const
        {
            return mUniqueKnotVecs[dir].size() - 1;
        }
        
        /// Total number of non-zero knot spans (i.e. elements)
        inline uint nonzeroKnotSpanN() const
        {
            uint c = 1;
            for(const auto& v : mUniqueKnotVecs)
                c *= v.size() - 1;
            return c;
        }
        
        /// We are dealing with spaces in R^3
        inline uint paramDimN() const { return 3; }
        
        /// Get b-spline basis function at given parametric coordinate.
        /// Prefer to use other member function however for efficiency.
        inline DoubleVecVec tensorBasis(const double u,
                                        const double v,
                                        const double w) const
        {return tensorBasis(u,v,w,span(u,v,w));}
        
        /// Get b-spline basis with precalculated knot span index
        DoubleVecVec tensorBasis(const double u,
                                 const double v,
                                 const double w,
                                 const UIntVec& span) const;
        
        DoubleVec basis(const double u,
                        const double v,
                        const double w) const
        {return basis(u,v,w,span(u,v,w));}
        
        /// Get basis functions (not in tensor product format)
        DoubleVec basis(const double u,
                        const double v,
                        const double w,
                        const UIntVec& span) const
        {
            DoubleVecVec bvec = tensorBasis(u,v,w,span);
            DoubleVec vec;
            for(const auto& b_k : bvec[W])
                for(const auto& b_j : bvec[V])
                    for(const auto& b_i : bvec[U])
                        vec.emplace_back(b_i * b_j * b_k);
            return vec;
        }
        
        /// Get b-spline basis function derivatives at given parametric coordinate.
        /// Prefer to use other member function however for efficiency.
        DoubleVecVec tensorBasisDers(const double u,
                                     const double v,
                                     const double w,
                                     const DerivOrder der = D1) const
        {return tensorBasisDers(u,v,w,span(u,v,w),der);}
        
        /// Get b-spline basis with precalculate knot span index
        DoubleVecVec tensorBasisDers(const double u,
                                     const double v,
                                     const double w,
                                     const UIntVec& span,
                                     const DerivOrder der = D1) const;
        
        /// Get basis function derivatives with tensor product multiplied out
        DoubleVec basisDers(const double u,
                            const double v,
                            const double w,
                            const DerivOrder der = D1) const
        {return basisDers(u,v,w,span(u,v,w),der);}
        
        /// Get basis function derivatives with tensor product multiplied out
        DoubleVec basisDers(const double u,
                            const double v,
                            const double w,
                            const UIntVec& span,
                            const DerivOrder der = D1) const
        {
            DoubleVecVec bvec = tensorBasisDers(u,v,w,span,der);
            DoubleVec vec;
            for(const auto& b_k : bvec[W])
                for(const auto& b_j : bvec[V])
                    for(const auto& b_i : bvec[U])
                        vec.emplace_back(b_i * b_j * b_k);
            return vec;
        }
        
        /// Get basis function derivatives with specified derivative
        DoubleVec basisDers(const double u,
                            const double v,
                            const double w,
                            const DerivType dtype) const
        { return basisDers(u,v,w,span(u,v,w), dtype); }
        
        /// Get basis function derivatives
        DoubleVec basisDers(const double u,
                            const double v,
                            const double w,
                            const UIntVec& span,
                            const DerivType dtype) const
        {
            const DoubleVecVec tbasis = tensorBasis(u,v,w, span);
            const DoubleVecVec tbasis_d = tensorBasisDers(u,v,w, span);
            DoubleVec dvec;
            
            if(DU == dtype) {
                for(const auto& b_k : tbasis[W])
                    for(const auto& b_j : tbasis[V])
                        for(const auto& bd_i : tbasis_d[U])
                            dvec.push_back(bd_i * b_j * b_k);
            }
            else if(DV == dtype) {
                
                for(const auto& b_k : tbasis[W])
                    for(const auto& bd_j : tbasis_d[V])
                        for(const auto& b_i : tbasis[U])
                            dvec.push_back(b_i * bd_j * b_k);
            }
            else if(DW == dtype) {
                
                for(const auto& bd_k : tbasis_d[W])
                    for(const auto& b_j : tbasis[V])
                        for(const auto& b_i : tbasis[U])
                            dvec.push_back(b_i * b_j * bd_k);
            }
            else
                throw std::runtime_error("Bad derivative type specified.");
            return dvec;
        }
        
        /// Return the non-zero local basis function indices for this parametric coordinate
        /// Returns a set of vectors corresponding to the indices in each parametric direction
        UIntVecVec localBasisIVec(const double u,
                                   const double v,
                                   const double w) const;
        
        /// Get the non-zero basis function indices using a row major numbering system.
        /// These are global in the sense that we return a single vector of indices rather
        /// than a tensor product format.
        UIntVec globalBasisIVec(const double u,
                                 const double v,
                                 const double w) const
        {
            UIntVec gb_vec;
            const uint nb_u = basisFuncN(U);
            const uint nb_v = basisFuncN(V);
            
            auto localvec = localBasisIVec(u,v,w);  // vector of basis functions in each direction
            
            for(const auto& w_index : localvec[W])
                for(const auto& v_index : localvec[V])
                    for(const auto& u_index : localvec[U])
                        gb_vec.push_back(w_index * nb_u * nb_v + v_index * nb_u + u_index);
            return gb_vec;
        }
        
        /// Get the knot indices for given parametric coordinate
        UIntVec span(const double u,
                     const double v,
                     const double w) const;
        
        /// Get local element indices (i,j,k) along each parametric direction
        /// given a 'global' element index.
        std::tuple<uint, uint, uint> localIndices(const uint iel) const;
        
        /// Extraction operator getter
        const std::vector<std::vector<double>>& extractionOperator(const uint iel,
                                                                   const ParamDir dir) const
        {
            switch(dir) {
                case ParamDir::U:
                    return mUExtractionOperators.at(iel);
                    break;
                case ParamDir::V:
                    return mVExtractionOperators.at(iel);
                    break;
                case ParamDir::W:
                    return mWExtractionOperators.at(iel);
                    break;
                default:
                    throw std::runtime_error("Bad parametric direction in Bezier element");
            }
        }
        
        /// Degree reduce in given parametric direction
        void degreeReduce(const ParamDir dir);
        
        /// Degree reduce space in both parametric directions
        void degreeReduce()
        {
            degreeReduce(U);
            degreeReduce(V);
            degreeReduce(W);
            
        }
        
        /// Raise degree in all parametric directions
        void degreeElevate()
        {
            degreeElevate(U);
            degreeElevate(V);
            degreeElevate(W);
        }
        
        /// Degree elevate in given parametric direction
        void degreeElevate(const ParamDir dir);
        
        /// Apply h-refinement (knot insertion) n times
        void hrefine(const uint n = 1);
        
        /// Apply graded h-refinement (knot insertion) n elements
        void graded_hrefine(const uint n, const double coeff);
        
        /// Load from an input stream
        void load(std::istream& ist);
        
        /// Remove internal knots
        void removeInternalKnots()
        {
            mKnotVecs[U].clear();
            for(size_t i = 0; i < mDegrees[U] + 1; ++i)
                mKnotVecs[U].push_back(0.0);
            for(size_t i = 0; i < mDegrees[U] + 1; ++i)
                mKnotVecs[U].push_back(1.0);
            
            mKnotVecs[V].clear();
            for(size_t i = 0; i < mDegrees[V] + 1; ++i)
                mKnotVecs[V].push_back(0.0);
            for(size_t i = 0; i < mDegrees[V] + 1; ++i)
                mKnotVecs[V].push_back(1.0);
            
            mKnotVecs[W].clear();
            for(size_t i = 0; i < mDegrees[W] + 1; ++i)
                mKnotVecs[W].push_back(0.0);
            for(size_t i = 0; i < mDegrees[W] + 1; ++i)
                mKnotVecs[W].push_back(1.0);
            
            init();
        }
        
        /// Insert knots such that the patch becomes a Bezier patch
        void convertToBezier()
        {
            for(uint iparam = 0; iparam < paramDimN(); ++iparam)
            {
                auto& kvec = mKnotVecs[iparam];
                const auto& unique_kvec = mUniqueKnotVecs[iparam];
                
                kvec.clear();
                for(size_t i = 0; i < unique_kvec.size(); ++i)
                {
                    const double kval = unique_kvec[i];
                    
                    if(i == 0 || i == unique_kvec.size() - 1)
                    {
                        for(size_t j = 0; j < mDegrees[iparam] + 1; ++j)
                            kvec.push_back(kval);
                    }
                    else
                    {
                        const uint mult = mDegrees[iparam] > 0 ? mDegrees[iparam]: 1;
                        for(size_t j = 0; j < mult; ++j)
                            kvec.push_back(kval);
                    }
                }
            }
            init();
        }
        
        
    private:
        
        /// recalculate unique knots and intervals
        void init();
        
        /// Clear all data
        void clear()
        {
            mKnotVecs.clear();
            mUniqueKnotVecs.clear();
            mDegrees.clear();
            mNumBasisVec.clear();
            mIntervalVec.clear();
            mUExtractionOperators.clear();
            mVExtractionOperators.clear();
            mWExtractionOperators.clear();
        }
        
        /// Is this a valid parametric coordinate? I.e. does it lie in the parametric space?
        bool validCoord(const double u,
                        const double v,
                        const double w) const
        {
            DoubleVec p{u,v,w};
            for(uint i = 0; i < paramDimN(); ++i)
                if(!boost::icl::contains(mIntervalVec[i], p[i]))
                    return false;
            return true;
        }
        
        /// Set the extraction operator for the given element index
        /// and parametric direction
        void setExtractionOperator(const uint iel,
                                   const ParamDir dir,
                                   const std::vector<std::vector<double>>& op)
        {
            switch(dir) {
                    
                case ParamDir::U:
                    mUExtractionOperators[iel] = op;
                    break;
                case ParamDir::V:
                    mVExtractionOperators[iel] = op;
                    break;
                case ParamDir::W:
                    mWExtractionOperators[iel] = op;
                    break;
                default:
                    throw std::runtime_error("Bad parametric direction in Bspline space");
            }
        }
        
        /// Global knot vectors
        DoubleVecVec mKnotVecs;
        
        /// Unique global knot vectors (define elements with non-zero
        /// parametric area
        DoubleVecVec mUniqueKnotVecs;
        
        /// Vector of degrees for each direction
        UIntVec mDegrees;
        
        /// Vector containing number of basis function in each direction
        UIntVec mNumBasisVec;
        
        /// Vector of interval that specify bounds of this space
        std::vector<Interval> mIntervalVec;
        
        /// The name of this space
        std::string mName;
        
        /// U-direction element extraction operators
        std::map<uint, std::vector<std::vector<double>>> mUExtractionOperators;
        
        /// V-direction element extraction operators
        std::map<uint, std::vector<std::vector<double>>> mVExtractionOperators;
        
        /// W- direction element extraction operators
        std::map<uint, std::vector<std::vector<double>>> mWExtractionOperators;
        
        /// Print function
        void printData(std::ostream& ost) const;
        
        /// overload output operator
        friend std::ostream& operator<<(std::ostream& ost, const BSplineSpace& s)
        {
            s.printData(ost);
            return ost;
        }
        
        /// Overload input operator
        friend std::istream& operator>>(std::istream& ist, BSplineSpace& s)
        {
            s.load(ist);
            return ist;
        }
        
    };
    
    /// helper functions
    
    /// Get the the (i,j,k) tuple given a 'space' element index and
    /// number of elements in the u- and v- parametric directions
    std::tuple<uint, uint, uint> localElementSpaceIndices(const uint ielem,
                                                          const uint nel_u,
                                                          const uint nel_v);
    
    // Get local (space) index for given vertex and # basis functions in
    // each parametric direction
    uint localBasisI(const Vertex v,
                     const uint nb_u,
                     const uint nb_v,
                     const uint nb_w);
    
    /// Get local (space) index for given vertex and space
    uint localBasisI(const Vertex v,
                     const BSplineSpace& s);
    
    /// Get local (space) basis indices for given edge
    UIntVec localBasisIVec(const Edge e,
                           const uint nb_u,
                           const uint nb_v,
                           const uint nv_w);
    
    /// Get local (space) indicies for a given edge and space
    UIntVec localBasisIVec(const Edge e, const BSplineSpace& s);
    
    /// Get ordered local vertex indices for the given edge and space
    std::pair<uint, uint> localBasisIPair(const Edge e,
                                          const BSplineSpace& s);
    
    /// Given a face enum and space, return the set of ordered
    /// local vertex indices. A positive order is defined by working anti-
    /// clockwise around the face, starting at the origin as defined
    /// by the local (u,v) face coord system.
    std::tuple<uint, uint, uint, uint> localBasisITuple(const Face f,
                                                        const BSplineSpace& s);
    
    /// Return the set of local basis function indices for a given face
    /// using the coordinate system as defined in the comments in base.h
    UIntVecVec localBasisIVec(const Face f, const BSplineSpace& s);
    
    /// Return the set of local basis function indices for a given face
    /// using the coordinate system as defined in the comments in base.h
    UIntVecVec localBasisIVec(const Face f,
                              const uint nb_u,
                              const uint nb_v,
                              const uint nb_w);
    
    /// Get the basis functions on each of the edges on this face
    /// Note: basis indices are repeated at vertices
    UIntVecVec localBasisOnEdgesIVec(const Face f,
                                     const BSplineSpace& s);
    
    /// Get the basis functions on each of the edges on this face
    /// Note: basis indices are repeated at vertices
    UIntVecVec localBasisOnEdgesIVec(const Face f,
                                     const uint nb_u,
                                     const uint nb_v,
                                     const uint nb_w);
    
    
    
    /// Reorder the local indices for a given face such that they align
    /// with the tuple that specifies the global vertex indices for the face
    void reorderLocalFaceIndices(std::vector<std::vector<uint>>& lindices,
                                 const std::tuple<uint, uint, uint, uint>& igvert);
    
    

    
    
}

#endif
