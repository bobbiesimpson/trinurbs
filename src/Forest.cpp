#include "Forest.h"
#include "Geometry.h"

namespace trinurbs
{
    Forest::Forest(const Geometry& g) : mGeom(&g)
    {
        const Forest& f = g.primalForest();
        mSpaces = f.mSpaces;
        mSpaceMap = f.mSpaceMap;
        mNodalConn = f.mNodalConn;
        mGlobalDofN = f.mGlobalDofN;
        
        initEdgeConn();
        initFaceConn();
    }
    
    Forest::Forest(const Forest& f)
    {
        clear();
        mGeom = f.mGeom;
        mSpaces = f.mSpaces;
        mSpaceMap = f.mSpaceMap;
        mNodalConn = f.mNodalConn;
        mEdgeConn = f.mEdgeConn;
        mFaceConn = f.mFaceConn;
        mEdgeSpaceMap = f.mEdgeSpaceMap;
        mCVertexSpaceMap = f.mCVertexSpaceMap;
        mFaceSpaceMap = f.mFaceSpaceMap;
        mElemIndexMap = f.mElemIndexMap;
        
        /// TODO
        
//        for(const auto& e : f.mElems)
//            mElems.insert(std::make_pair(e.first, e.second->copy()));
//        for(const auto& e : f.mBezierElems)
//            mBezierElems.insert(std::make_pair(e.first, e.second->copy()));
        mGlobalDofN = f.mGlobalDofN;
    }
    
    /// Assignment operator
    Forest& Forest::operator=(const Forest& f)
    {
        if(this == &f) // check for self assignment
            return *this;
        clear();
        mGeom = f.mGeom;
        mSpaces = f.mSpaces;
        mSpaceMap = f.mSpaceMap;
        mNodalConn = f.mNodalConn;
        mEdgeConn = f.mEdgeConn;
        mFaceConn = f.mFaceConn;
        mEdgeSpaceMap = f.mEdgeSpaceMap;
        mCVertexSpaceMap = f.mCVertexSpaceMap;
        mFaceSpaceMap = f.mFaceSpaceMap;
        mElemIndexMap = f.mElemIndexMap;
        
        /// TODO
//        for(const auto& e : f.mElems)
//            mElems.insert(std::make_pair(e.first, e.second->copy()));
//        for(const auto& e : f.mBezierElems)
//            mBezierElems.insert(std::make_pair(e.first, e.second->copy()));
        
        mGlobalDofN = f.mGlobalDofN;
        return *this;
    }
    
    /// clear all data
    void Forest::clear()
    {
        mSpaces.clear();
        mSpaceMap.clear();
        mNodalConn.clear();
        mEdgeConn.clear();
        mFaceConn.clear();
        mEdgeSpaceMap.clear();
        mCVertexSpaceMap.clear();
        mFaceSpaceMap.clear();
        mElemIndexMap.clear();
        
        /// TODO
//        mElems.clear();
//        mBezierElems.clear();
        mGlobalDofN = std::make_pair(false, 0);
    }
    
    Sign Forest::globalEdgeDir(const uint ispace, const Edge e) const
    {
        auto v_pair = globalVertexPairI(ispace, e);
        if(v_pair.second > v_pair.first)
            return Sign::POSITIVE; // increasing global vertex index == positive edge
        else return Sign::NEGATIVE; // negative otherwise
    }
    
    uint Forest::globalVertexI(const uint ispace, const Vertex v) const
    {
        const BSplineSpace& s = space(ispace);
        const uint nbasis_u = s.basisFuncN(U);
        const uint nbasis_v = s.basisFuncN(V);
        const uint nbasis_w = s.basisFuncN(W);
        
        // vertex is specified in terms of (i,j,k) indices in parametric space
        switch(v) {
            case Vertex::VERTEX0: return globalI(ispace,0,0,0);
            case Vertex::VERTEX1: return globalI(ispace, nbasis_u - 1,0,0);
            case Vertex::VERTEX2: return globalI(ispace,0,0,nbasis_w - 1);
            case Vertex::VERTEX3: return globalI(ispace, nbasis_u - 1,0, nbasis_w - 1);
            case Vertex::VERTEX4: return globalI(ispace,0,nbasis_v - 1,0);
            case Vertex::VERTEX5: return globalI(ispace, nbasis_u - 1, nbasis_v -1,0);
            case Vertex::VERTEX6: return globalI(ispace,0, nbasis_v - 1, nbasis_w - 1);
            case Vertex::VERTEX7: return globalI(ispace, nbasis_u - 1,nbasis_v - 1, nbasis_w - 1);
        }
    }
    
    std::pair<uint, uint> Forest::globalVertexPairI(const uint ispace, const Edge e) const
    {
        
        switch(e) {
            case Edge::EDGE0:
                return std::make_pair(globalVertexI(ispace, Vertex::VERTEX0),
                                      globalVertexI(ispace, Vertex::VERTEX1));
            case Edge::EDGE1:
                return std::make_pair(globalVertexI(ispace, Vertex::VERTEX2),
                                      globalVertexI(ispace, Vertex::VERTEX3));
            case Edge::EDGE2:
                return std::make_pair(globalVertexI(ispace, Vertex::VERTEX0),
                                      globalVertexI(ispace, Vertex::VERTEX2));
            case Edge::EDGE3:
                return std::make_pair(globalVertexI(ispace, Vertex::VERTEX1),
                                      globalVertexI(ispace, Vertex::VERTEX3));
            case Edge::EDGE4:
                return std::make_pair(globalVertexI(ispace, Vertex::VERTEX4),
                                      globalVertexI(ispace, Vertex::VERTEX5));
            case Edge::EDGE5:
                return std::make_pair(globalVertexI(ispace, Vertex::VERTEX6),
                                      globalVertexI(ispace, Vertex::VERTEX7));
            case Edge::EDGE6:
                return std::make_pair(globalVertexI(ispace, Vertex::VERTEX4),
                                      globalVertexI(ispace, Vertex::VERTEX6));
            case Edge::EDGE7:
                return std::make_pair(globalVertexI(ispace, Vertex::VERTEX5),
                                      globalVertexI(ispace, Vertex::VERTEX7));
            case Edge::EDGE8:
                return std::make_pair(globalVertexI(ispace, Vertex::VERTEX0),
                                      globalVertexI(ispace, Vertex::VERTEX4));
            case Edge::EDGE9:
                return std::make_pair(globalVertexI(ispace, Vertex::VERTEX1),
                                      globalVertexI(ispace, Vertex::VERTEX5));
            case Edge::EDGE10:
                return std::make_pair(globalVertexI(ispace, Vertex::VERTEX2),
                                      globalVertexI(ispace, Vertex::VERTEX6));
            case Edge::EDGE11:
                return std::make_pair(globalVertexI(ispace, Vertex::VERTEX3),
                                      globalVertexI(ispace, Vertex::VERTEX7));
        }
    }
    
    void Forest::load(std::istream& ist)
    {
        while(true) {
            BSplineSpace s;
            if(!(ist >> s))
                break;
            std::cout << "found space\n";
            addSpace(s);
        }
        endOfLoop(ist, '}', "Could not read B-spline space\n");
        
        while(true) {
            ConnVecInput cv;
            if(!(ist >> cv))
                break;
            addConnectivityVec(cv.name, cv.data);
        }
        endOfLoop(ist, '}', "Could not read connectivity vectors\n");
        initEdgeConn();
        initFaceConn();
    }
    
    void Forest::print(std::ostream& ost) const
    {
        ost << "Forest: " << spaceN() << " trunks\n";
        for(const auto& n : mSpaceMap) {
            ost << mSpaces[n.second] << "\n";
        }
        ost << "Nodal connectivity:\n";
        for(const auto& c : mNodalConn)
            ost << c.first << ": " << c.second << "\n";
    }
    
    void Forest::initNodalConn()
    {
        // The purpose of this function is to generate the nodal connectivity
        // of Bspline spaces for any level of refinement. The general algorithm
        // is to use the connectivity of the geometry to guide this process.
        
        const Forest& primal = geometry()->primalForest();  // reference to geometry (primal forest)
        
        std::map<uint, std::vector<uint>> temp_conn;        // temporary nodal connectivity - will be copied at end
        std::map<uint, std::vector<uint>> edge_map;         // map from global edge index to global node indices
        std::map<uint, std::vector<std::vector<uint>>> face_map;  // map from global face index to global node indicies
        std::map<uint, uint> vx_map;                        // map from geometry vertex index to new vertex index
        
        uint current_index = 0;                             // the current global node index
        
        // Now loop over Bspline spaces
        for(uint ispace = 0; ispace < spaceN(); ++ispace) {
            
            const BSplineSpace& s = space(ispace);
            std::vector<int> gnode_vec(s.basisFuncN(), -1); // global node indices for this space: -1 means unassigned
            
            // First fill in vertex global indices if previously assigned
            for(uint ivertex = 0; ivertex < NVERTICES; ++ivertex)  {
                
                const Vertex vertex = vertexType(ivertex);
                const uint geom_ivertex = primal.globalVertexI(ispace, vertex);
                auto find = vx_map.find(geom_ivertex);
                if(find != vx_map.end())
                    gnode_vec[localBasisI(vertex, s)] = find->second;
                else {
                    const uint new_index = current_index++;
                    vx_map[geom_ivertex] = new_index;
                    gnode_vec[localBasisI(vertex, s)] = new_index;
                }
            }
            
            // Now fill in edge nodal indices if previously assigned
            for(uint iedge = 0; iedge < NEDGES; ++iedge) {
                
                const Edge edge = edgeType(iedge);
                
                // get global edge index
                const uint global_iedge = globalEdgeI(ispace, iedge);
                auto find = edge_map.find(global_iedge);
                
                // Get ordered local vertex indices for this edge
                auto vpair = localBasisIPair(edge, s); 
                
                assert(vpair.first != vpair.second);
                
                // Get the sign of the edge according to positively increasing global indices
                const Sign sign = (gnode_vec[vpair.second] > gnode_vec[vpair.first]) ? Sign::POSITIVE : Sign::NEGATIVE;
                
                // Get the local (space) indices along this edge
                auto local_indices = localBasisIVec(edge, s);
                
                // order local indices in the positive direction
                if(Sign::NEGATIVE == sign)
                    std::reverse(local_indices.begin(), local_indices.end());
                
                // edge indices have already been assigned so let's assign them to this space
                if(find != edge_map.end()) {
                    const auto glb_indices = find->second;  // the vector of global edge indices
                    
                    assert(local_indices.size() == glb_indices.size());
                    
                    for(uint i = 0; i < local_indices.size(); ++i) {
                        
                        // check to make sure it hasn't previously been assigned (e.g. vertex)
                        if(gnode_vec[local_indices[i]] != -1)
                            continue;
                        else
                            gnode_vec[local_indices[i]] = glb_indices[i];
                    }
                }
                else { // the edge node indices have not been assigned. Let's do it now.
                    
                    std::vector<uint> gvec;
                    
                    for(const auto& i : local_indices) {
                        
                        // Be careful to consider already assigned vertex index
                        if(gnode_vec[i] != -1)
                            gvec.push_back(gnode_vec[i]);
                        else {
                            const uint gindex = current_index++;
                            gvec.push_back(gindex);
                            gnode_vec[i] = gindex;
                        }
                    }
                    
                    // and store the new global edge connectivity
                    edge_map[global_iedge] = gvec;
                }
            }
            
            // Now fill in face nodal indices if previously assigned
            for(uint iface = 0; iface < NFACES; ++iface)
            {
                const Face face = faceType(iface);
                
                // get the global face index
                const uint global_iface = globalFaceI(ispace, iface);
                auto find = face_map.find(global_iface);
                
                // We define a global face origin by the lowest global
                // index in the set that defines the face.
                // We can then define a 2D coord system where
                // the +ve u-axis is defined by the edge from the origin to
                // the lowest global vertex index (out of the possible 2).
                
                // Given this definition, we need to 'rotate' and/or
                // 'flip' the local indices depending on the orientation
                // of the present face to this global face coord system.
                
                // get local indices for this face ordered using the
                // coord system as detailed in base.h
                auto local_indices = localBasisIVec(face, s);
                
                // now re-order the local indices to match that of the
                // 'global' indices
                
                
                
            }
            
            
            // now fill up unassigned indices
            std::vector<uint> unsigned_gvec;
            for(auto& i : gnode_vec) {
                if(i != -1) {
                    unsigned_gvec.push_back(static_cast<uint>(i));
                    continue;
                }
                i = current_index++;
                unsigned_gvec.push_back(i);
            }
            temp_conn[ispace] = unsigned_gvec;
        }
        mNodalConn = temp_conn;
    }
    
    void Forest::initEdgeConn()
    {
    
    }
    
    void Forest::initFaceConn()
    {
    
    }
    
    std::istream& operator>>(std::istream& ist, Forest& f)
    {
        f.load(ist);
        return ist;
    }
    
    std::ostream& operator<<(std::ostream& ost, const Forest& f)
    {
        f.print(ost);
        return ost;
    }
    
    

}