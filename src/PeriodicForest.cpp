#include "PeriodicForest.h"
#include "Geometry.h"

#include <vector>
#include <set>

namespace trinurbs
{
//    void PeriodicForest::initGeometry()
//    {
//
//    }
    
    void PeriodicForest::initAnalysisData()
    {
        clearAnalysisData();
        
        // Generate the Greville abscissa points for the Forest by looping
        // over all elements.
        
        // Points that lie on the surface are flagged by adding their indices
        // to the relevant surface point vector. Likewise all interior points
        // are added to the interior point vector.
        
        /// temp cache for Greville points for Forest
        std::map<uint, Point3D> greville_pt_cache;
        std::set<uint> greville_surface_cache;
        std::set<uint> greville_interior_cache;
        
        // First populate the cache with all the Greville points (in physical space)
        // By looping over each of the B-spline spaces of this Forest
        for(uint ispace = 0; ispace < spaceN(); ++ispace)
        {
            const auto& analysis_space = space(ispace);
            const auto& geometry = this->geometry();
            
            // loop over Greville abscissa
            for(uint igpt = 0; igpt < analysis_space.grevilleAbscissaPtN(); ++igpt)
            {
                const uint g_index = globalI(ispace, igpt);
                const ParamCoord gpt_param = analysis_space.grevilleAbscissaPt(igpt);
                const Point3D gpt_phys = geometry->eval(gpt_param.u, gpt_param.v, gpt_param.w, ispace);
                greville_pt_cache[g_index] = gpt_phys;
            }
        }
        
        // Now determine which points lie on the surface of the periodic cell and those in the interior
        for(const auto& mpair : greville_pt_cache)
        {
            const auto& index = mpair.first;
            const auto& gpt = mpair.second;
            
            if(logically_equal(std::abs(gpt[0]), 1.0) || logically_equal(std::abs(gpt[1]), 1.0) || logically_equal(std::abs(gpt[2]), 1.0))
                greville_surface_cache.insert(index);
            else
                greville_interior_cache.insert(index);
        }
        
//        std::cout << "surface points\n\n";
//        for(const auto& pt : greville_surface_cache)
//            std::cout << greville_pt_cache[pt] << "\n";
//        
//        std::cout << "interior points\n\n";
//        for(const auto& pt : greville_interior_cache)
//            std::cout << greville_pt_cache[pt] << "\n";

        // Now copy cache to member data
        mGrevilleAbscissaMap = greville_pt_cache;
        mSurfaceGrevilleISet = greville_surface_cache;
        mInteriorGrevilleISet = greville_interior_cache;
        
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