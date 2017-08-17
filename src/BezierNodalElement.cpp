#include "BezierNodalElement.h"
#include "Point4D.h"
#include "Geometry.h"

namespace trinurbs
{
    
    Point3D BezierNodalElement::eval(const double xi,
                                     const double eta,
                                     const double zeta) const
    {
        
        // We get a reference to the parent element, convert to its local
        // parent coordinates and evaluate using the basis that belongs to the
        // geometry (often much coarser than the analysis basis)
        
        const auto pel = parent();
        const auto p_gpt = transformToParentElParentCoord(GPt3D(xi,eta,zeta));
        const auto& b = pel->basis(p_gpt.xi, p_gpt.eta, p_gpt.zeta);
        const auto gvec = pel->globalBasisIVec();
        
        Point4D x;
        for(uint i = 0; i < gvec.size(); ++i)
            x += geometry().controlPt(gvec[i]) * b[i];
        return x.asCartesian();
    }
    
    Point3D BezierNodalElement::tangent(const double xi,
                                        const double eta,
                                        const double zeta,
                                        const ParamDir dir) const
    {
        // efficiency is crucial here for fast BE analysis
        const auto pel = parent();
        const auto ppt = transformToParentElParentCoord(GPt3D(xi,eta,zeta));
        const auto gvec = pel->globalBasisIVec();
        
        const auto& basis = pel->basis(ppt.xi, ppt.eta, ppt.zeta);
        //        std::cout << ppt << "\n";
        //        std::cout << basis << "\n";
        
        double w_interp = 0.0;
        for(uint i = 0; i < gvec.size(); ++i)
            w_interp += geometry().controlPt(gvec[i]).getWeight() * basis[i];
        
        // get non-rational bspline basis function parametric derivatives
        DerivType dtype = (U == dir) ? DU : (V == dir ? DV: DW);
        
        const auto& basis_der = pel->basisDers(ppt.xi, ppt.eta, ppt.zeta, dtype);
        
        //        std::cout << basis_der << "\n";
        
        // calculate weight function derivative
        double w_der = 0.0;
        for(uint i = 0; i < gvec.size(); ++i)
            w_der += geometry().controlPt(gvec[i]).getWeight() * basis_der[i];
        
        Point3D result;
        for(uint i = 0; i < gvec.size(); ++i)
        {
            const auto term = geometry().controlPt(gvec[i]).getWeight() *
            (1.0 / w_interp * basis_der[i] -
             1.0 / (w_interp * w_interp) * w_der * basis[i]);
            result += geometry().controlPt(gvec[i]).asCartesian() * term;
        }
        return result;
    }
}