#include "MultiscaleBezierNodalElement.h"
#include "BezierNodalElement.h"

namespace trinurbs {
    
    Point3D MultiscaleBezierNodalElement::eval(const double xi,
                                               const double eta,
                                               const double zeta) const
    {
        // get point in 'cell' space in [-1,1]^3
        const auto pt = microBezierElement()->eval(xi, eta, zeta);
        return macroBezierElement()->GeometryElement::eval(pt[0], pt[1], pt[2]);
        
    }
    
    
}