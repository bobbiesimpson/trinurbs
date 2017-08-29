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
        
        // first rescale geometry to biunit interval [-1,1]^3
        
        // set up maps which hold vectors of elements
        mFaceElMaps.clear();
        mFaceElMaps.resize(NFACES);
        
        mFaceGeomEls.clear();
        mFaceGeomEls.resize(NFACES);
        
        // First determine the sets of geometry elements with faces on each macro
        // face
        const Geometry* geom = geometry();
        const Forest& primal = geom->primalForest();
        
        std::vector<Interval> cell_intervals;
        for(uint i = 0; i < 3; ++i)
            cell_intervals.push_back(boost::icl::construct<boost::icl::continuous_interval<double>>(-1.0,1.0));

        // Determine set of geometry elements and local (geometry element) faces
        // that lie on each face of the periodic domain. These elements (and any
        // subelements thereof) will be used to connect neighbouring periodic
        // cells
        
        for(uint iel = 0; iel < primal.elemN(); ++iel)
        {
            const auto geom_el = primal.bezierElement(iel);  // geometry element
            const auto gconn = geom_el->globalBasisIVec();
            
            // if all cpts of geometry element lie in plane of face
            // then element intersects this macro cell face
            
            // loop over cell faces
            for(uint iface = 0; iface < NFACES; ++iface)
            {
                const Face face = faceType(iface);
                std::vector<Point3D> cpts;
                
                for(const auto& ipt : gconn)
                    cpts.push_back(geom->controlPt(ipt).asCartesian());
                
                if(allCPtsLieOnPeriodicCellFace(cpts, face))
                    mFaceGeomEls[iface].push_back(std::make_pair(geom_el,face));
            }
            
            
        }
        
        // Construct vector of 'child' elements for each geometry element
        // that lies on a face/edge/vertex of the cell
        
        
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
    
}