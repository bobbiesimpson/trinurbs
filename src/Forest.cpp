#include "Forest.h"
#include "Geometry.h"
#include "BezierNodalElement.h"

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
        for(const auto& e : f.mBezierElems)
            mBezierElems.insert(std::make_pair(e.first, e.second->copy()));
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
        for(const auto& e : f.mBezierElems)
            mBezierElems.insert(std::make_pair(e.first, e.second->copy()));
        
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
        mBezierElems.clear();
        mGlobalDofN = std::make_pair(false, 0);
    }
    
    const NAnalysisElement* Forest::bezierElement(const uint i) const
    {
        // First search cache for element and return if found
        auto e = mBezierElems.find(i);
        if(e != mBezierElems.end())
            return e->second.get();
        
        // otherwise we need to construct the relevant Bezier element
        //
        // First find local index of element (in space)
        uint start_i = 0;
        for(uint s = 0; s < spaceN(); ++s) {
            const uint el_n = space(s).nonzeroKnotSpanN();
            if((i - start_i) > (el_n - 1)) {
                start_i += el_n;
                continue;
            }
            const uint local_i = i - start_i;
            
            // create the element
            auto r = mBezierElems.insert(std::make_pair(i, make_unique<BezierNodalElement>(this,
                                                                                           s,
                                                                                           local_i)));
            if(!r.second)
                error("Failed attempt to create element");
            
            // create mapping from global element index to space and local element index
            auto result = mElemIndexMap.insert(std::make_pair(i, std::make_pair(s, local_i)));
            if(!result.second)
                 error("Failure inserting bezier element index mapping.");
            
            auto el = mBezierElems[i].get();
            
            // finally, create a referece to the parent element in the primal forest
            
            // search for parent element
            for(uint ielem = 0; ielem < geometry()->primalForest().space(s).nonzeroKnotSpanN(); ++ielem) {
                const auto pel = geometry()->primalForest().bezierElement(s, ielem);
                if(pel->contains(*el))
                    el->setParent(pel);
            }
            assert(el->parent() != nullptr);
            return el;
        }
        error("Failed to create Bezier element"); return nullptr;
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
    
    std::tuple<uint, uint, uint, uint> Forest::globalFaceVertexITuple(const uint ispace,
                                                                      const Face f) const
    {
        const auto local_ivertices = localBasisITuple(f,space(ispace));
        return std::make_tuple(globalVertexI(ispace, std::get<0>(local_ivertices)),
                               globalVertexI(ispace, std::get<1>(local_ivertices)),
                               globalVertexI(ispace, std::get<2>(local_ivertices)),
                               globalVertexI(ispace, std::get<3>(local_ivertices)));
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
            
            //
            // REQUIRES TESTING!!!!!!
            //
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
                // the next lowest global vertex index. The v-axis
                // is then defined from the origin to the index with the
                // next least difference.
                
                // Given this definition, we need to 'rotate' and/or
                // 'flip' the local indices depending on the orientation
                // of the present face to this global face coord system.
                
                // get local indices for this face ordered using the
                // coord system as detailed in base.h
                auto ilocal_basis_vec = localBasisIVec(face, s);
                
                // get the local indices of the vertices that define this
                // face (returned as a tuple)
                auto ilocal_verts = localBasisITuple(face, s);
                
                // and fill up the tuple which defines the
                auto iglobal_verts = std::make_tuple(gnode_vec[std::get<0>(ilocal_verts)],
                                                     gnode_vec[std::get<1>(ilocal_verts)],
                                                     gnode_vec[std::get<2>(ilocal_verts)],
                                                     gnode_vec[std::get<3>(ilocal_verts)]);
                
                // now re-order the local indices to match that of the
                // 'global' indices
                reorderLocalFaceIndices(ilocal_basis_vec, iglobal_verts);
                
                if(find != face_map.end())
                {
                    // We reach here if indices have already been assigned to this
                    // face
                    const auto iglobal_basis_vec = find->second;
                    
                    assert(ilocal_basis_vec.size() == iglobal_basis_vec.size());
                    assert(ilocal_basis_vec[0].size() == iglobal_basis_vec[0].size());
                    
                    for(size_t i = 0; i < ilocal_basis_vec.size(); ++i)
                    {
                        for(size_t j = 0; j < ilocal_basis_vec[0].size(); ++j)
                        {
                            const auto lindex = ilocal_basis_vec[i][j];
                            
                            // if previously assigned, silently carry on
                            if(gnode_vec[lindex] != -1)
                                continue;
                            else
                                gnode_vec[lindex] = iglobal_basis_vec[i][j];
                        }
                    }
                }
                else
                {
                    // face vertices have not been assigned. Let's do it now.
                    std::vector<std::vector<uint>> gmatrix;
                    
                    for(const auto& local_row : ilocal_basis_vec)
                    {
                        gmatrix.push_back(std::vector<uint>());
                        auto& current = gmatrix.back();

                        for(const auto& i : local_row)
                        {
                            if(gnode_vec[i] != -1)
                                current.push_back(gnode_vec[i]);
                            else
                            {
                                const uint gindex = current_index++;
                                current.push_back(gindex);
                                gnode_vec[i] = gindex;
                            }
                        }
                    }
                    face_map[global_iface] = gmatrix;
                }
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
        // clear entries before initialising
        mEdgeConn.clear();
        mEdgeSpaceMap.clear();
        mCVertexSpaceMap.clear();
        
        std::map<std::vector<uint>, uint> edge_map; // mapping from edge vertices to global edge index
        std::vector<Edge> e_types{
            Edge::EDGE0, Edge::EDGE1, Edge::EDGE2, Edge::EDGE3,
            Edge::EDGE4, Edge::EDGE5, Edge::EDGE6, Edge::EDGE7,
            Edge::EDGE8, Edge::EDGE9, Edge::EDGE10, Edge::EDGE11
        };
        
        uint edge_count = 0;
        for(uint ispace = 0; ispace < spaceN(); ++ispace)
        {
            std::vector<uint> edge_conn; // edge connectivity for this space
            std::vector<uint> c_vertices;
            
            for(const auto& e : e_types)
            {
                auto v_pair = globalVertexPairI(ispace, e);
                auto e_conn = globalIVec(ispace, localBasisIVec(e, space(ispace)));
                std::sort(e_conn.begin(), e_conn.end());
                
                c_vertices.push_back(v_pair.first);
                c_vertices.push_back(v_pair.second);
                
                if(v_pair.first > v_pair.second) // make sure ordering is +ve
                    v_pair = std::make_pair(v_pair.second, v_pair.first);
                
                auto search = edge_map.find(e_conn);
                
                if(search != edge_map.end())
                    edge_conn.push_back(search->second);
                else
                {
                    edge_map[e_conn] = edge_count;
                    edge_conn.push_back(edge_count);
                    ++edge_count;
                }
            }
            std::sort(c_vertices.begin(), c_vertices.end());
            auto last = std::unique(c_vertices.begin(), c_vertices.end());
            c_vertices.erase(last, c_vertices.end());
            for(const auto& c : c_vertices)
                mCVertexSpaceMap[c].push_back(ispace);
            for(const auto& e : edge_conn)
                mEdgeSpaceMap[e].push_back(ispace);
            mEdgeConn[ispace] = edge_conn;
        }
    }
    
    void Forest::initFaceConn()
    {
        // clear appropriate entries
        mFaceConn.clear();
        mFaceSpaceMap.clear();
        
        // map from basis indices around edge of face to a global face index
        std::map<std::vector<uint>, uint> face_map;
        
        // vector of all the face types
        std::vector<Face> face_types{
            Face::FACE0, Face::FACE1, Face::FACE2, Face::FACE3,
            Face::FACE4, Face::FACE5
        };
        
        uint face_count = 0;
        for(uint ispace = 0; ispace < spaceN(); ++ispace)
        {
            // let's generate a vector of global basis indices
            // for all edges on this face
            
            // face connectivity for this space
            std::vector<uint> face_conn;
            
            // loop over faces
            for(const auto& f : face_types)
            {
                // get vectors of basis indices along all edges of this face
                auto igbasis_edgevecs = localBasisOnEdgesIVec(f, space(ispace));
                
                // now condense in a single vector
                std::vector<uint> edge_basis;
                for(const auto& vec : igbasis_edgevecs)
                    for(const auto& e : vec)
                        edge_basis.push_back(e);
                
                // sort and remove duplicates
                std::sort(edge_basis.begin(), edge_basis.end());
                edge_basis.erase(std::unique(edge_basis.begin(), edge_basis.end()), edge_basis.end());
                
                auto search = face_map.find(edge_basis);
                
                // if the face has been encountered previously, add the space
                if(search != face_map.end())
                    face_conn.push_back(search->second);
                
                // otherwise we create a new face index
                else
                {
                    face_map[face_conn] = face_count;
                    face_conn.push_back(face_count);
                    ++face_count;
                }
            }
            
            // finally generate the space map of face and add
            // the face connectivity for this space to the global
            // member variable
            for(const auto& f : face_conn)
                mFaceSpaceMap[f].push_back(ispace);
            mFaceConn[ispace] = face_conn;
        }
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