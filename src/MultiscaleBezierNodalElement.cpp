#include "MultiscaleBezierNodalElement.h"
#include "BezierNodalElement.h"
#include "MultiscaleForest.h"

namespace trinurbs {
    
    Point3D MultiscaleBezierNodalElement::eval(const double xi,
                                               const double eta,
                                               const double zeta) const
    {
        // get point in 'cell' space in [-1,1]^3
        const auto pt = microBezierElement()->eval(xi, eta, zeta);
        
        return macroBezierElement()->GeometryElement::eval(pt[0], pt[1], pt[2]);
        
    }
    
    /// get determinant of jacobian
    double MultiscaleBezierNodalElement::jacDet(const double xi,
                                                const double eta,
                                                const double zeta) const
    {
        const auto micro_jdet = microBezierElement()->jacDet(xi, eta, zeta);
        const auto pt = microBezierElement()->eval(xi, eta, zeta);
        const auto macro_jdet = macroBezierElement()->jacDet(pt[0], pt[1], pt[2]);
        return micro_jdet * macro_jdet;
    }
    
    void MultiscaleBezierNodalElement::init()
    {
        if(nullptr == microBezierElement())
            return;
        
        // mapping from local basis index to PeriodicForest global index
        const auto& micro_conn = microBezierElement()->globalBasisIVec();
        
        // mapping from PeriodicForest global index to global index in the multiscale forest
        const auto& macro_conn = multiscaleForest()->globalBasisIVec(macroElI());
        
        std::vector<uint> temp_conn(micro_conn.size());
        
        for(size_t ilocal = 0; ilocal < micro_conn.size(); ++ilocal)
            temp_conn[ilocal] = macro_conn[micro_conn[ilocal]];
        
        // set member data
        mMicroGlobalBasisIVec = temp_conn;
    }
    
    
    
    
    
    
}