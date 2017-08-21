#ifndef TRINURBS_MULTISCALE_BEZIER_NODAL_ELEMENT_H
#define TRINURBS_MULTISCALE_BEZIER_NODAL_ELEMENT_H

#include "base.h"
#include "Point3D.h"
#include "AnalysisElement.h"

namespace trinurbs {
    
    /// The representation of an element which lives in the microscale geometry
    /// It is simply defined by the relevant element in the microscale
    /// geometry and macroscale geometry
    
    class MultiscaleBezierNodalElement
    {
        
    public:
        
        // Construct with pointers to macro and micro element
        MultiscaleBezierNodalElement(const NAnalysisElement* macroel = nullptr,
                                     const NAnalysisElement* microel = nullptr)
        :
            mpMacroScaleElement(macroel),
            mpMicroScaleElement(microel)
        {}
        
        /// evaluate physical coordinate given parent coordinate
        Point3D eval(const double xi,
                     const double eta,
                     const double zeta) const;
        
        // macro element getter
        const NAnalysisElement* macroBezierElement() const
        { return mpMacroScaleElement; }
        
        /// micro element getter
        const NAnalysisElement* microBezierElement() const
        { return mpMicroScaleElement; }
        
    protected:
        
    private:
        
        
        const NAnalysisElement* mpMacroScaleElement;
        
        const NAnalysisElement* mpMicroScaleElement;
        
    };
    
}

#endif