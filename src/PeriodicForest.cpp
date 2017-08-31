#include "PeriodicForest.h"
#include "Geometry.h"

#include <vector>

namespace trinurbs
{
    void PeriodicForest::initGeometry()
    {
        // generate data structures that allow node indices to be
        // determined across vertices, faces and edges.
        
        // e.g. given two 'macro' elements connected by a vertex,
        // edge or face, what is the correspondence between node
        // indices along that connected vertex, edge or face?
        // This is the fundamental requirement for generating
        // connectivity at the fine-scale.
        
        // We can loop over elements and populate sets of node indices
        // that are connected to each vertex, edge and face.
        
        // We need to construct the vectors of local elements and faces
        // which are connected to each macro cell face. We also need to
        // generate equivalent vectors for the various macro cell face permutations
        // (e.g. flips, rotations) which are required during construction
        // of the global connectivity
        
        // Clear all geometry data first
        clearGeometryData();
        
        // First determine the sets of geometry elements with faces on each macro
        // face
        const Geometry* geom = geometry();
        const Forest& primal = geom->primalForest();

        // Determine sets of geometry elements that lie on each face,
        // edge and vertex of the macro (cell).
        
        for(uint iel = 0; iel < primal.elemN(); ++iel)
        {
            const auto geom_el = primal.bezierElement(iel);  // geometry element
            
            // loop over geometry element vertices
            for(uint ielvertex = 0; ielvertex < NVERTICES; ++ielvertex)
            {
                const Vertex elvertex = vertexType(ielvertex);
            
                // perhaps not the most robust approach, but we simply evaluate the coordinate
                // at each vertex of the element and check for equality with the cell vertices.
                const auto pair = liesOnPeriodicCellVertex(geom_el->eval(parentIntervalVertex(elvertex)));
                
                // this geometry vertex lies on the cell vertex. Add it to the map.
                if(pair.first)
                    mVertexGeomEls[pair.second].push_back(std::make_pair(geom_el, elvertex));
            }
            
            // loop over cell edges
            for(uint icelledge = 0; icelledge < NEDGES; ++icelledge)
            {
                const Edge celledge = edgeType(icelledge);
                
                // loop over edges of the geometry element
                for(uint ieledge = 0; ieledge < NEDGES; ++ieledge)
                {
                    const Edge eledge = edgeType(ieledge);
                    const auto econn = geom_el->globalBasisIVec(eledge);
                    
                    // populate vector of control points on this element face
                    std::vector<Point3D> cpts;
                    for(const auto& ipt : econn)
                        cpts.push_back(geom->controlPt(ipt).asCartesian());
                    
                    if(allCPtsLieOnPeriodicCellEdge(cpts, celledge))
                        mEdgeGeomEls[celledge].push_back(std::make_pair(geom_el, eledge));
                }
            }
            
            // loop over cell faces
            for(uint icellface = 0; icellface < NFACES; ++icellface)
            {
                const Face cellface = faceType(icellface);
                
                // loop over geometry element faces
                for(uint ielface = 0; ielface < NFACES; ++ielface)
                {
                    const Face elface = faceType(ielface);
                    const auto fconn = geom_el->globalBasisIVec(elface);
                    
                    // populate vector of control points on this element face
                    std::vector<Point3D> cpts;
                    for(const auto& ipt : fconn)
                        cpts.push_back(geom->controlPt(ipt).asCartesian());
                    
                    // if all control points on element face lie on cell face
                    // then this element must lie on the cell face
                    if(allCPtsLieOnPeriodicCellFace(cpts, cellface))
                        mFaceGeomEls[cellface].push_back(std::make_pair(geom_el,elface));
                }
            }
        }
        
    }
    
    void PeriodicForest::initAnalysisData()
    {
        
        // Construct vector of 'child' elements for each geometry element
        // that lies on a face/edge/vertex of the cell.
        
        // The order in which these are stored is important. We order them
        // according to face and edge coordinate systems described in base.h
        
        // Edge elements are ordered in terms of positive edges
        
        // Face elements are ordered in a local u-v face coordinate system
        // first in terms of positively increasing u values and then v value.
        
        // Vertex element have no order since we assume only one analysis
        // element at each vertex.
        
        // map from geometry element to child analysis elements
        // could be optimised to reduce to set of geometry elements
        // that only lie on faces/edges/vertices of cell geometry.
        
        clearAnalysisData();
        
        // TODO: fill up children maps
        
        std::map<const NAnalysisElement*, std::vector<const NAnalysisElement*>> geom2child_map;
        
        // loop over (refined) analysis elements and populate child map
        for(uint iel = 0; iel < elemN(); ++iel)
        {
            auto el = bezierElement(iel);
            geom2child_map[el->parent()].push_back(el);
        }
        
        // loop over cell (macro) geometry faces
//        for(uint imacroface = 0; imacroface < NFACES; ++imacroface)
//        {
//            // get vector of geometry/face pairs for this face
//            const auto& geomface_vec = geomElFacePairVec(imacroface);
//            
//            // loop over geoemtry face pairs
//            for(const auto& gfacepair : geomface_vec)
//            {
//                const auto& g_pel = gfacepair.first; // pointer to geometry el
//                const auto& face = gfacepair.second; // local face of geom el
//                
//                // loop over children of this geometry element
//                for(const auto& pel : geom2child_map[g_pel])
//                {
//                    // temp code just to check map
//                    if(pel->liesOnParentElFace(face))
//                        mFaceElMaps[imacroface][FacePermutation::NOPERMUTATION].push_back(std::make_tuple(pel,face,FacePermutation::NOPERMUTATION));
//                }
//                
//            }
//        }
        
        
        std::cout << "Finished setting up analysis data\n";
    }
    
    bool allCPtsLieOnPeriodicCellFace(const std::vector<Point3D>& cpts,
                                      const Face f)
    {
        uint comp = 0;
        double val = 0.0;
        
        switch(f)
        {
            
            case Face::FACE0:
                comp = 1;   // v-component
                val = -1.0; // equals -1.0 on this face
                break;
            case Face::FACE1:
                comp = 1;   // v-component
                val = 1.0; // equals -1.0 on this face
                break;
            case Face::FACE2:
                comp = 0;   // u-component
                val = -1.0; // equals -1.0 on this face
                break;
            case Face::FACE3:
                comp = 0;   // u-component
                val = 1.0; // equals -1.0 on this face
                break;
            case Face::FACE4:
                comp = 2;   // w-component
                val = -1.0; // equals -1.0 on this face
                break;
            case Face::FACE5:
                comp = 2;   // w-component
                val = 1.0; // equals -1.0 on this face
                break;
        }
        
        for(const auto& p : cpts)
            if(!logically_equal(p[comp], val))
                return false;
        return true;
    }
    
    bool allCPtsLieOnPeriodicCellEdge(const std::vector<Point3D>& cpts,
                                      const Edge celledge)
    {

        std::vector<std::pair<uint,double>> cvec{};
        
        switch(celledge)
        {
            case Edge::EDGE0:
                cvec.push_back(std::make_pair(1, -1.0));
                cvec.push_back(std::make_pair(2, -1.0));
                break;
            case Edge::EDGE1:
                cvec.push_back(std::make_pair(1, -1.0));
                cvec.push_back(std::make_pair(2, 1.0));
                break;
            case Edge::EDGE2:
                cvec.push_back(std::make_pair(0, -1.0));
                cvec.push_back(std::make_pair(1, -1.0));
                break;
            case Edge::EDGE3:
                cvec.push_back(std::make_pair(0, 1.0));
                cvec.push_back(std::make_pair(1, -1.0));
                break;
            case Edge::EDGE4:
                cvec.push_back(std::make_pair(1, 1.0));
                cvec.push_back(std::make_pair(2, -1.0));
                break;
            case Edge::EDGE5:
                cvec.push_back(std::make_pair(1, 1.0));
                cvec.push_back(std::make_pair(2, 1.0));
                break;
            case Edge::EDGE6:
                cvec.push_back(std::make_pair(0, -1.0));
                cvec.push_back(std::make_pair(1, 1.0));
                break;
            case Edge::EDGE7:
                cvec.push_back(std::make_pair(0, 1.0));
                cvec.push_back(std::make_pair(1, 1.0));
                break;
            case Edge::EDGE8:
                cvec.push_back(std::make_pair(0, -1.0));
                cvec.push_back(std::make_pair(2, -1.0));
                break;
            case Edge::EDGE9:
                cvec.push_back(std::make_pair(0, 1.0));
                cvec.push_back(std::make_pair(2, -1.0));
                break;
            case Edge::EDGE10:
                cvec.push_back(std::make_pair(0, -1.0));
                cvec.push_back(std::make_pair(2, 1.0));
                break;
            case Edge::EDGE11:
                cvec.push_back(std::make_pair(0, 1.0));
                cvec.push_back(std::make_pair(2, 1.0));
                break;
        }
        
        
        for(const auto& p : cpts)
            for(const auto& pair : cvec)
                if(!logically_equal(p[pair.first], pair.second))
                    return false;
        return true;
    }
    
    
    std::pair<bool, Vertex> liesOnPeriodicCellVertex(const Point3D& p)
    {
        for(uint icellvtx = 0; icellvtx < NVERTICES; ++icellvtx)
        {
            const auto vtx = vertexType(icellvtx);
            const auto gpt = parentIntervalVertex(vtx);
            const Point3D ptvtx(gpt.xi, gpt.eta,gpt.zeta);
            if(ptvtx == p)
                return std::make_pair(true,vtx);
        }
        return std::make_pair(false, Vertex::VERTEX0);
    }
    
}