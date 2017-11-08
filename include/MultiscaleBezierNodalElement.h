#ifndef TRINURBS_MULTISCALE_BEZIER_NODAL_ELEMENT_H
#define TRINURBS_MULTISCALE_BEZIER_NODAL_ELEMENT_H

#include "base.h"
#include "Point3D.h"
#include "AnalysisElement.h"


namespace trinurbs {
    
    /// Forward declaration
    class MultiscaleForest;
    
    /// The representation of an element which lives in the microscale geometry
    /// It is simply defined by the relevant element in the microscale
    /// geometry and macroscale geometry
    
    class MultiscaleBezierNodalElement
    {
        
    public:
        
        // Construct with pointers to macro and micro element
        MultiscaleBezierNodalElement(const MultiscaleForest* multiscaleforest,
                                     const uint imacro,
                                     const uint imicro,
                                     const NAnalysisElement* macroel = nullptr,
                                     const NAnalysisElement* microel = nullptr)
        :
            mpMultiscaleForest(multiscaleforest),
            mMacroElI(imacro),
            mMicroElI(imicro),
            mpMacroScaleElement(macroel),
            mpMicroScaleElement(microel)
        {
            init();
        }
        
        /// evaluate physical coordinate given parent coordinate
        Point3D eval(const double xi,
                     const double eta,
                     const double zeta) const;
        
        /// get determinant of jacobian
        double jacDet(const double xi,
                      const double eta,
                      const double zeta) const;
        
        /// NUmber of basis functions
        uint basisFuncN() const
        {
            return microBezierElement()->basisFuncN();
        }
        
        /// Get basis functions in microscale element
        std::vector<double> basis(const double xi,
                                  const double eta,
                                  const double zeta) const
        {
            return microBezierElement()->basis(xi, eta, zeta);
        }
        
        /// Basis function derivatives of micro scale element
        std::vector<double> basisDers(const double xi,
                                      const double eta,
                                      const double zeta,
                                      const DerivType dtype) const
        {
            return microBezierElement()->basisDers(xi, eta, zeta, dtype);
        }
        
        /// Degree getter
        uint degree(const ParamDir dir, const uint comp = 0) const
        {
            return microBezierElement()->degree(dir, comp);
        }
        
        // macro element getter
        const NAnalysisElement* macroBezierElement() const
        { return mpMacroScaleElement; }
        
        /// micro element getter
        const NAnalysisElement* microBezierElement() const
        { return mpMicroScaleElement; }
        
        /// Get integratino order of element (micro element)
        UIntVec integrationOrder(const uint offset = 0) const
        {
            return microBezierElement()->integrationOrder(offset);
        }
        
        /// Multiscale forest getter
        const MultiscaleForest* multiscaleForest() const
        {
            return mpMultiscaleForest;
        }
        
        /// Macro element getter
        const uint macroElI() const
        {
            return mMacroElI;
        }
        
        /// Micro element getter
        const uint microElI() const
        {
            return mMicroElI;
        }
        
        /// Get mapping from local micro element basis index
        /// to global index in the multiscale forest
        std::vector<uint> globalBasisIVec() const
        {
            return mMicroGlobalBasisIVec;
        }
        
    protected:
        
    private:
        
        /// Initialise data structures, specifically global connectivity
        void init();
        
        /// Pointer to multiscale forest
        const MultiscaleForest* mpMultiscaleForest;
        
        /// Macro element index
        const uint mMacroElI;
        
        /// Micro element index
        const uint mMicroElI;
        
        /// Pointer to macro scale element
        const NAnalysisElement* mpMacroScaleElement;
        
        /// Pointer to micro scale element
        const NAnalysisElement* mpMicroScaleElement;
        
        /// Connectivity array of the micro element
        std::vector<uint> mMicroGlobalBasisIVec;
        
    };
    
}

#endif