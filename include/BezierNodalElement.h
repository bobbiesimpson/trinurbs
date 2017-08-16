#ifndef TRINURBS_BEZIER_NODAL_ELEMENT_H
#define TRINURBS_BEZIER_NODAL_ELEMENT_H

#include <mutex>

#include "AnalysisElement.h"
#include "BSplineSpace.h"
#include "Forest.h"
#include "base.h"
#include "NURBSCommon.h"

namespace trinurbs
{
    
    ///
    /// A reprenestation of a Bezier element that
    /// provides fast evaluation of (rational) B-spline
    /// basis functions and derivatives for analysis.
    ///
    class BezierNodalElement : public NAnalysisElement
    {
        
    public:
        
        /// Constructor
        BezierNodalElement(const Forest* f,
                           const uint ispace,
                           const uint iel,
                           const BezierNodalElement* pel = nullptr)
        :
            NAnalysisElement(*f->geometry(),
                             ispace,
                             f->knotIntervals(ispace,iel)),
            mForest(f),
            mSpace(&f->space(ispace)),
            mElemI(iel),
            mpParentEl(pel)
        {
            // initialise connectivity
            ParamCoord c = getParamCoord(lowerBound(U),
                                         lowerBound(V),
                                         lowerBound(W));
            mLocalBasisIVec = space()->globalBasisIVec(c.u,
                                                       c.v,
                                                       c.w);
            
            for(const auto& i : mLocalBasisIVec)
                mGlobalBasisIVec.emplace_back(forest()->globalI(spaceI(),i));
            
        }
        
        /// Virtual copy function
        std::unique_ptr<AnalysisElement> copy() const override
        {
            return make_unique<BezierNodalElement>(*this);
        }
        
        /// Scalar basis. Only one component.
        uint componentN() const override { return 1; }
        
        /// Number of non-zero basis functions over this element.
        uint basisFuncN() const override { return localBasisIVec().size(); }
        
        /// Evaluate physical coordinates using Bezier extraction
        virtual Point3D eval(const double xi,
                             const double eta,
                             const double zeta) const override;
        
        /// Evaluate tangent vector using Bezier extraction
        virtual Point3D tangent(const double xi,
                                const double eta,
                                const double zeta,
                                const ParamDir dir) const override;
        
        /// TODO: We will probably need to write a normal() function
        /// which computes the normal on a given face of the element
        /// used for computing force vector terms
        
        /// Override normal evaluation
//        virtual Point3D normal(const double u, const double v) const override
//        {
//            return cross(tangent(u,v,S), tangent(u,v,T)).normalise();
//        }
        
        /// Override jacobian determinant evaluation
        virtual double jacDet(const double xi,
                              const double eta,
                              const double zeta) const override
        {
            DoubleVecVec jtemp({{t1.asVec()}, {t2.asVec()}, {t3.asVec()}});
            transpose(jtemp);
            
            return det3x3(jtemp) * jacDetParam(gp.xi,gp.eta,gp.zeta);
        }
        
        /// Override jacobian evaluation
        virtual DoubleVecVec jacob(const double xi,
                                   const double eta,
                                   const double zeta) const override
        {
            DoubleVecVec jacob_param 

            jacob_param.push_back(tangent(xi,eta,zeta,U).asVec());
            jacob_param.push_back(tangent(xi,eta,zeta,V).asVec());
            jacob_param.push_back(tangent(xi,eta,zeta,W).asVec());
            
            transpose(jacob_param);
            
            const auto jacob_parent = jacobParam(xi,eta,zeta);
            
            DoubleVecVec r{
                { 0.0, 0.0, 0.0 },
                { 0.0, 0.0, 0.0 },
                { 0.0, 0.0, 0.0 }
            };
            
            for(uint i = 0; i < 3; ++i)
                for(uint j = 0; j < 3; ++j)
                    for(uint k = 0; k < 3; ++k)
                        r[i][j] += jacob_param[i][k] * jacob_parent[k][j];
            return r;
        }
        
        /// Return local basis function indices that are non-zero over the
        /// relevant Bspline space
        UIntVec localBasisIVec() const override
        {
            return mLocalBasisIVec;
        }
        
        /// Return the global (to the forest) basis function indices
        UIntVec globalBasisIVec() const override
        {
            return mGlobalBasisVecI;
        }
        
        /// REturn the basis function values
        DoubleVec basis(const double xi,
                        const double eta,
                        const double zeta) const override
        {
            auto indices = space()->localIndices(localElementI());
            
            // Using a reference makes a huge difference to speed
            const auto& op_u = space()->extractionOperator(indices.first, U);
            const auto& op_v = space()->extractionOperator(indices.second, V);
            const auto& op_w = space()->extractionOperator(indices.second, W);
            
            const auto b_u = nurbshelper::bernsteinPolynomial(xi, degree(U));
            const auto b_v = nurbshelper::bernsteinPolynomial(eta, degree(V));
            const auto b_w = nurbshelper::bernsteinPolynomial(zeta, degree(W));
            
            const uint n_u = op_u.size();
            std::vector<double> basis_u(n_u, 0.0);
            for(uint i = 0; i < n_u; ++i)
                for(uint j = 0; j < op_u[0].size(); ++j)
                    basis_u[i] += op_u[i][j] * b_u[j];
            
            const uint n_v = op_v.size();
            std::vector<double> basis_v(n_v, 0.0);
            for(uint i = 0; i < n_v; ++i)
                for(uint j = 0; j < op_v[0].size(); ++j)
                    basis_v[i] += op_v[i][j] * b_v[j];
            
            const uint n_w = op_w.size();
            std::vector<double> basis_w(n_w, 0.0);
            for(uint i = 0; i < n_w; ++i)
                for(uint j = 0; j < op_w[0].size(); ++j)
                    basis_w[i] += op_w[i][j] * b_w[j];
            
            std::vector<double> final;
            for(uint k = 0; k < basis_w.size(); ++k)
                for(uint j = 0; j < basis_v.size(); ++j)
                    for(uint i = 0; i < basis_u.size(); ++i)
                        final.push_back(basis_u[i] * basis_v[j] * basis_w[k]);
            return final;
        }
        
        /// TODO: finsihed coverting to trivariate up to here
        /// on 16/8/2017
        /// Get the local basis derivatives either in S or T direction
        DoubleVec localBasisDers(const double xi,
                                 const double eta,
                                 const double zeta,
                                 const DerivType dtype) const override
        {
            auto indices = space()->localIndices(localElementI());
            
            // Using a reference makes a huge difference to speed
            const auto& op_u = space()->extractionOperator(indices.first, S);
            const auto& op_v = space()->extractionOperator(indices.second, T);
            
            // Get Bertnein basis in each parametric direction taking account
            // of derivatives.
            DoubleVec bu;
            DoubleVec bv;
            if(DerivType::DS == dtype) {
                bu = nurbshelper::bernsteinPolynomialDeriv(u, degree(S));
                bv = nurbshelper::bernsteinPolynomial(v, degree(T));
            }
            else {
                bu = nurbshelper::bernsteinPolynomial(u, degree(S));
                bv = nurbshelper::bernsteinPolynomialDeriv(v, degree(T));
            }
            
            const uint n_u = op_u.size();
            std::vector<double> basis_u(n_u, 0.0);
            for(uint i = 0; i < n_u; ++i)
                for(uint j = 0; j < op_u[0].size(); ++j)
                    basis_u[i] += op_u[i][j] * bu[j];
            
            const uint n_v = op_v.size();
            std::vector<double> basis_v(n_v, 0.0);
            for(uint i = 0; i < n_v; ++i)
                for(uint j = 0; j < op_v[0].size(); ++j)
                    basis_v[i] += op_v[i][j] * bv[j];
            
            std::vector<double> final;
            const auto param_j = jacobParam(u, v);
            const double jterm = (DerivType::DS == dtype) ? 1.0/param_j[0][0] : 1.0 / param_j[1][1];
            
            for(uint j = 0; j < basis_v.size(); ++j)
                for(uint i = 0; i < basis_u.size(); ++i)
                    final.push_back(basis_u[i] * basis_v[j] * jterm);
                    
                    return final;
        }
        
        /// Basis function degrees
        uint degree(const ParamDir dir, const uint comp = 0) const override
        {
            return space()->degree(dir);
        }
        
        /// Basis function degrees
        UIntVec degree(const uint comp = 0) const override
        {
            return space()->degree();
        }
        
        /// TODO
        void print(std::ostream& ost) const override
        {
            GeometryElement::print(ost);
        }
        
        /// Local element index getter
        uint localElementI() const { return mElemI; }
        
        /// Forest accessor
        const Forest* forest() const {return mForest;}
        
        /// Accessor for space
        const BSplineSpace* space() const{ return mSpace; }
        
        
    protected:
        
    private:
        
        /// Reference to forest
        const Forest* mForest;
        
        /// Reference to space this element belongs to
        const BSplineSpace* mSpace;
        
        /// Element index (local to space)
        const uint mElemI;
        
        /// Reference to 'parent' element which exists
        /// in the primal forest
        const BezierNodalElement* mpParentEl;
        
        /// Cache the knot span indices for evaluating basis functions
        UIntVec mLocalBasisIVec;
        
        /// A vector containing the global basis function indices. (i.e. the global forest basis index)
        UIntVec mGlobalBasisIVec;
    };
}

#endif