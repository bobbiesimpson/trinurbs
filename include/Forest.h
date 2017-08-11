#ifndef TRINURBS_FOREST_H
#define TRINURBS_FOREST_H

#include <vector>
#include <map>
#include <iostream>
#include <memory>
#include <cassert>
#include <utility>
#include <string>
#include <fstream>
#include <algorithm>

#include "base.h"
#include "BSplineSpace.h"

namespace trinurbs
{
    /// A forest class is responsible for the set of Bspline spaces, associated
    /// degrees of freedom (with D components) and the connectivity between
    /// the trees (b-spline patches).
    
    /// Forward declarations
    class Geometry;
    
    class Forest {
        
    public:
        
        /// Default constructor. No spaces created.
        Forest()
        :   mGeom(nullptr),
            mGlobalDofN(std::make_pair(false, 0))
        {}
        
        /// Construct with a geometry object. The geometry pointer may be null
        /// if this forest represents a primal forest.
        Forest(const Geometry& g);
        
        /// Construct with vector of bspline spaces and connectivity.
        Forest(const std::vector<BSplineSpace>& s_vec,
               const std::vector<UIntVec>& conn)
        :   mGeom(nullptr),
            mSpaces(s_vec),
            mGlobalDofN(std::make_pair(false, 0))
        {
            if(s_vec.size() != conn.size())
                error("Invalid data passed in to Forest constructor");
            for(uint s = 0; s < spaceN(); ++s) {
                mNodalConn.insert(std::make_pair(s, conn[s]));
                mSpaceMap.insert(std::make_pair("space" + std::to_string(s), s));
            }
            initEdgeConn();
            initFaceConn();
        }
        
        /// Use default destructor
        virtual ~Forest() = default;
        
        /// Copy constructor
        Forest(const Forest& f);
        
        /// Assignment operator
        Forest& operator=(const Forest& f);
        
        /// Move constructor. Simply move all the data across
        Forest(Forest&& f) = default;
        
        /// Move assignment operator. As before, move all the data across.
        Forest& operator=(Forest&& f) = default;
        
        /// clear all data
        void clear();
        
        /// Get the geometry
        const Geometry* geometry() const
        {
            //assert(mGeom != nullptr);
            return mGeom;
        }
        
        /// Geometry setter
        void setGeometry(const Geometry* g) { mGeom = g; }
        
        /// Const b-spline space accessor
        const BSplineSpace& space(const uint i) const
        {
            assert(i < mSpaces.size());
            
            return mSpaces.at(i);
        }
        
        /// Non-const B-spline space getter
        const BSplineSpace& space(const uint i)
        {
            assert(i < mSpaces.size());
            
            return mSpaces[i];
        }
        
        /// Const connectivity vector accessor
        const UIntVec& connectivityVec(const uint i) const
        {
            return mNodalConn.at(i);
        }
        
        /// Total number of elements in forest
        uint elemN() const
        {
            uint n = 0;
            for(uint i = 0; i < spaceN(); ++i)
                n += space(i).nonzeroKnotSpanN();
            return n;
        }
        
        /// Get the number of basis functions for a given space
        uint basisFuncN(const uint ispace) const
        {
            return space(ispace).basisFuncN();
        }
        
        ///
        /// TODO
        ///
        
//        /// Get element at global index
//        const NAnalysisElement* element(const uint i) const;
//        
//        /// Get an element given a space index and local element index
//        const NAnalysisElement* element(const uint ispace,
//                                        const uint ielem) const
//        { return element(globalElI(ispace, ielem)); }
//        
//        /// Get element at global index
//        const NAnalysisElement* bezierElement(const uint i) const;
//        
//        /// Get an element given a space index and local element index
//        const NAnalysisElement* bezierElement(const uint ispace,
//                                              const uint ielem) const
//        { return bezierElement(globalElI(ispace, ielem)); }
        
        /// END TODO
        
        /// Get the number of elements on the given space
        const uint elemN(const uint ispace) const
        {
            return space(ispace).nonzeroKnotSpanN();
        }
        
        /// Get a global element index given a local (space) index numbering
        uint globalElI(const uint ispace,
                       uint i,
                       uint j,
                       uint k) const
        {
            // check indices are valid
            assert(i < space(ispace).nonzeroKnotSpanN(ParamDir::U));
            assert(j < space(ispace).nonzeroKnotSpanN(ParamDir::V));
            assert(k < space(ispace).nonzeroKnotSpanN(ParamDir::W));
            
            const uint nel_u = space(ispace).nonzeroKnotSpanN(ParamDir::U);
            const uint nel_v = space(ispace).nonzeroKnotSpanN(ParamDir::V);
            
            return globalElI(ispace, nel_u * nel_v * k + nel_u * j + i);
        }
        
        /// Get a global element index given an space index and local element index
        uint globalElI(const uint ispace,
                       const uint ielem) const
        {
            assert(ielem < space(ispace).nonzeroKnotSpanN());
            
            uint glb_index = 0;
            for(uint i = 0; i < ispace; ++i)
                glb_index += space(i).nonzeroKnotSpanN();
            return glb_index + ielem;
        }
        
        /// Number of B-spline spaces
        inline uint spaceN() const { return mSpaces.size(); }
        
        /// Get the global index for this space and local basis function index
        uint globalI(const uint sp, const uint i) const
        {
            assert(sp < spaceN() && i < space(sp).basisFuncN());
            
            return mNodalConn.at(sp)[i];
        }
        
        /// Is the given edge degenerate (i.e. equal connectivity along edge)
        bool degenerateEdge(const uint ispace,
                            const uint iedge) const;
        
        /// Get a global index given indices for each parametric direction
        uint globalI(const uint sp,
                     const uint i,
                     const uint j,
                     const uint k) const
        {
            const uint nu = space(sp).basisFuncN(U);
            const uint nv = space(sp).basisFuncN(V);
            
            const uint index = k * nu * nv + j * nu + i;
            return globalI(sp, index);
        }
        
        /// Given a set of local indices for the given space, return a vector of
        /// global nodal indices.
        UIntVec globalIVec(const uint sp, const UIntVec& lvec) const
        {
            assert(sp < spaceN());
            
            UIntVec rvec; // return vector
            for(const auto& v : lvec) {
                assert(v < space(sp).basisFuncN());
                rvec.push_back(globalI(sp, v));
            }
            return rvec;
        }
        
        /// Get global edge index
        uint globalEdgeI(const uint sp, const uint e) const
        {
            assert(sp < spaceN());
            assert(e < NEDGES);
            
            return mEdgeConn.at(sp)[e];
        }
        
        /// Get set of all global edge indices for this space
        std::vector<uint> globalEdgeIVec(const uint sp) const
        {
            assert(sp < spaceN());
            
            return mEdgeConn.at(sp);
        }
        
        /// Get global face index
        uint globalFaceI(const uint ispace, const uint iface)
        {
            assert(ispace < spaceN());
            assert(iface < NFACES);
            
            return mFaceConn.at(ispace)[iface];
        }
        
        /// An advancement over connectedEdges() in that this function can determine
        /// if there are multiple edges connected between the two spaces. This happens
        /// for instance in the torus topology.
//        std::pair<bool, std::map<Edge, Edge>> multiConnectedEdges(const uint ispace1,
//                                                                  const uint ispace2) const;
        
        /// Return the (local) connected edges between the given spaces
        /// Assumes only one connected edge therefore deprecated.
//        bool connectedEdges(const uint ispace1,
//                            const uint ispace2,
//                            Edge& edge1,
//                            Edge& edge2) const;
//        
//        /// Return the set of connected edges between two given spaces
//        /// and global edge index.
//        bool connectedEdges(const uint ispace1,
//                            const uint ispace2,
//                            const uint iedge,
//                            Edge& edge1,
//                            Edge& edge2) const;
        
        /// Number of edges in this forest
        uint edgeN() const
        {
            return mEdgeSpaceMap.size();
        }
        
        /// Return the (local) connected vertices between two spaces
//        bool connectedVertices(const uint ispace1,
//                               const uint ispace2,
//                               Vertex& vert1,
//                               Vertex& vert2) const;
        
        /// Is this edge positive or negative?
        Sign globalEdgeDir(const uint ispace, const Edge e) const;
        
        /// Get edge direction using integer index
        Sign globalEdgeDir(const uint ispace, const uint e) const
        { return globalEdgeDir(ispace, edgeType(e)); }
        
        /// get a global node index given a space index and local vertex index
        uint globalVertexI(const uint ispace, const Vertex v) const;
        
        /// Get vertex type using integer index
        uint globalVertexI(const uint ispace, const uint v) const
        { return globalVertexI(ispace, vertexType(v)); }
        
        /// Get the global vertex indices given a space and local edge index
        std::pair<uint, uint> globalVertexPairI(const uint ispace, const Edge e) const;
        
        /// Get vertex pair using integer indices.
        std::pair<uint, uint> globalVertexPairI(const uint ispace, const uint e) const
        { return globalVertexPairI(ispace, edgeType(e)); }
        
        /// Get the knot interval pairs for this element
        DoublePairVec knotIntervals(const uint ispace,
                                    const uint iel) const
        {
            assert(ispace < spaceN() && iel < space(ispace).nonzeroKnotSpanN());
            
            const BSplineSpace& bspace = space(ispace);
            DoublePairVec p_vec;
            
            const uint nel_u = bspace.uniqueKnotN(U) - 1;
            const uint nel_v = bspace.uniqueKnotN(V) - 1;
            
            /// Need to check this works during debugging
            const auto tuple = localElementSpaceIndices(iel, nel_u, nel_v);
            
            const uint i_index = std::get<0>(tuple);;
            const uint j_index = std::get<1>(tuple);;
            const uint k_index = std::get<2>(tuple);
            
            p_vec.emplace_back(std::make_pair(bspace.uniqueKnot(i_index, U), bspace.uniqueKnot(i_index + 1, U)));
            p_vec.emplace_back(std::make_pair(bspace.uniqueKnot(j_index, V), bspace.uniqueKnot(j_index + 1, V)));
            p_vec.emplace_back(std::make_pair(bspace.uniqueKnot(k_index, W), bspace.uniqueKnot(k_index + 1, W)));
            
            return p_vec;
        }
        
        /// Get the global dof for this forest
        uint globalDofN() const
        {
            if(!mGlobalDofN.first) {
                uint max = 0;
                for(uint ispace = 0; ispace < spaceN(); ++ispace) {
                    const BSplineSpace& s = space(ispace);
                    for(uint ibasis = 0; ibasis < s.basisFuncN(); ++ibasis) {
                        const uint g_index = globalI(ispace, ibasis);
                        if(g_index > max)
                            max = g_index;
                    }
                }
                mGlobalDofN = std::make_pair(true, max + 1);
            }
            return mGlobalDofN.second;
        }
        
        /// Get the vector of connected spaces given a global edge index
        std::vector<uint> connectedSpacesOnEdge(const uint iedge) const
        {
            return mEdgeSpaceMap.at(iedge);
        }
        
        /// Get the vector of connected spaces given a global vertex
        std::vector<uint> connectedSpacesOnVertex(const uint ivert) const
        {
            return mCVertexSpaceMap.at(ivert);
        }
        
        /// Get space indices connected to given face index
        std::vector<uint> connectedSpacesOnFace(const uint iface) const
        {
            return mFaceSpaceMap.at(iface);
        }
        
        /// Get space indices that are vertex adjacent (not edge adjcent) to the given space
        /// and vertex
//        std::vector<uint> vertexConnectedSpaces(const uint ispace,
//                                                const Vertex v) const
//        {
//            std::vector<uint> e_vec; // spaces that are edge adjacent which need to be excluded
//            switch (v) {
//                case Vertex::VERTEX0:
//                    e_vec = {globalEdgeI(ispace, 0), globalEdgeI(ispace, 2)};
//                    break;
//                case Vertex::VERTEX1:
//                    e_vec = {globalEdgeI(ispace, 0), globalEdgeI(ispace, 3)};
//                    break;
//                case Vertex::VERTEX2:
//                    e_vec = {globalEdgeI(ispace, 2), globalEdgeI(ispace, 1)};
//                    break;
//                case Vertex::VERTEX3:
//                    e_vec = {globalEdgeI(ispace, 1), globalEdgeI(ispace, 3)};
//                    break;
//            }
//            
//            
//            std::vector<uint> s_vec; // vector of spaces connected at edges (including ispace)
//            for(const auto& e : e_vec)
//                for(const auto& s : connectedSpacesOnEdge(e))
//                    s_vec.push_back(s);
//            std::sort(s_vec.begin(), s_vec.end());
//            auto l = std::unique(s_vec.begin(), s_vec.end()); s_vec.erase(l, s_vec.end());
//            auto full_ivec = connectedSpacesOnVertex(globalVertexI(ispace, v));
//            std::sort(full_ivec.begin(), full_ivec.end());
//            std::vector<uint> rvec;
//            std::set_difference(s_vec.begin(), s_vec.end(), full_ivec.begin(), full_ivec.end(), std::inserter(rvec, rvec.begin()));
//            return rvec;
//        }
        
        /// Apply uniform hrefinement with n levels of refinement
        void hrefine(const uint n)
        {
            if(0 == n)
                return;
            for(auto& s : mSpaces)
                s.hrefine(n);
            initNodalConn();
            clearElementData();
        }
        
        /// Apply degree reduction n times
        void degreeReduce(const uint n)
        {
            if(0 == n) {
                std::cerr << "Degree reduction of 0 requested: carrying on\n";
                return;
            }
            for(uint i = 0; i < n; ++i) {
                for(auto& s : mSpaces) {
                    const auto dvec = s.degree();
                    if(dvec[0] == 0 || dvec[1] == 0 || dvec[2]) {
                        std::cerr << "Cannot degree reduce beyond C^0";
                        break;
                    }
                    s.degreeReduce();
                }
            }
            initNodalConn();
            clearElementData();
        }
        
        /// Apply degree elevation n times
        void degreeElevate(const uint n)
        {
            if(0 == n) {
                std::cerr << "Degree elevate of 0 requested: carrying on\n";
                return;
            }
            const uint max_degree = 20;
            
            for(uint i = 0; i < n; ++i) {
                for(auto& s : mSpaces) {
                    const auto dvec = s.degree();
                    if(dvec[0] >= max_degree || dvec[1] >= max_degree || dvec[2] >= max_degree) {
                        std::cerr << "Reached max degree p=20";
                        break;
                    }
                    s.degreeElevate();
                }
            }
            initNodalConn();
            clearElementData();
        }
        
        /// Load from an input stream
        void load(std::istream& ist);
        
        /// Print this Forest to the given output stream
        void print(std::ostream& ost) const;
        
    protected:
        
        /// Add a space to the Forest
        void addSpace(const BSplineSpace& s)
        {
            auto i = mSpaceMap.find(s.name());
            if(i != mSpaceMap.end())
                error("Cannot add two spaces with identical names");
            mSpaces.push_back(s);
            mSpaceMap.insert(std::make_pair(s.name(), mSpaces.size() - 1));
        }
        
        void addConnectivityVec(const std::string& s, const std::vector<uint>& v)
        {
            auto p = validSpace(s);
            if(!p.first)
                error("Bad connectivity vector specified. Cannot add to forest.");
            mNodalConn.insert(std::make_pair(p.second, v));
        }
        
        /// Is there a space in this forest identified by the given string?
        std::pair<bool, uint> validSpace(const std::string& s)
        {
            auto i = mSpaceMap.find(s);
            if(i == mSpaceMap.end())
                return std::make_pair(false, 0);
            return std::make_pair(true, i->second);
        }
        
    private:
        
        /// This is reuquired to construct the nodal connectivity after refinement
        // Assumes edge connectivity exists
        void initNodalConn();
        
        /// Construct edge connectivity
        void initEdgeConn();
        
        /// Construct face connectivity
        void initFaceConn();
        
        /// Clear all data related to elements (used after applying refinement)
        void clearElementData()
        {
            /// TODO
            
//            mElems.clear();
//            mBezierElems.clear();
            mElemIndexMap.clear();
            mGlobalDofN = std::make_pair(false, 0);
        }
        
        /// Pointer to the geometry object associated with this forest.
        const Geometry* mGeom;
        
        /// Vector of rational b-spline patches denoted as 'trees'.
        std::vector<BSplineSpace> mSpaces;
        
        /// A mapping that identifies valid spaces and returns their index
        std::map<std::string, uint> mSpaceMap;
        
        /// Connectivity of trees
        std::map<uint, std::vector<uint> > mNodalConn;
        
        /// Edge connectivity of trees
        std::map<uint, std::vector<uint>> mEdgeConn;
        
        /// Face connectivity of trees
        std::map<uint, std::vector<uint>> mFaceConn;
        
        /// Map from edge index to connected spaces
        std::map<uint, std::vector<uint>> mEdgeSpaceMap;
        
        /// Map from corner vertex indices to connected space indices
        std::map<uint, std::vector<uint>> mCVertexSpaceMap;
        
        /// Map from face index to connected spaces
        std::map<uint, std::vector<uint>> mFaceSpaceMap;
        
        /// Map from an edge index to 'sub'-edges that may arise due
        /// to T-junctions
//        std::map<uint, std::vector<uint>> mSubEdgeMap;
        
        /// Map from a sub-edge to its parent edge index
//        std::map<uint, uint> mSuperEdgeMap;
        
        /// A mapping from a global element index (over all elements over all spaces) to
        /// a space index and element index local to the particular space
        mutable std::map<uint, std::pair<uint, uint>> mElemIndexMap;
        
        /// TODO
        
//        /// Map from a (space, local element index) pair to an element instance
//        mutable std::map<uint, std::unique_ptr<NAnalysisElement>> mElems;
//        
//        /// Map from a (space, local element index) pair to a bezier element instance
//        mutable std::map<uint, std::unique_ptr<NAnalysisElement>> mBezierElems;
        
        /// Cached global dof
        mutable std::pair<bool, uint> mGlobalDofN;
        
        /// Overload input operator
        friend std::istream& operator>>(std::istream& ist, Forest& f);
        
        /// Overload output operator
        friend std::ostream& operator<<(std::ostream& ost, const Forest& f);
        
    };
}

#endif
