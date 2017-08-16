#include "GeometryElement.h"
#include "Geometry.h"

namespace trinurbs
{
    
    GeometryElement::GeometryElement(const Geometry& g,
                                     const uint ispace,
                                     const DoublePairVec& knots)
    :
        mGeom(g),
        mSpaceI(ispace),
        mMutex(std::make_shared<std::mutex>()),
        mSize(std::make_pair(false, 0.0))
    {
        assert(knots.size() >= 3);
        for(const auto& kp : knots)
            mKnotIntervals.emplace_back(boost::icl::construct<boost::icl::continuous_interval<double>>(kp.first,
                                                                                                       kp.second,
                                                                                                       boost::icl::interval_bounds::closed()));
    }
    
    /// Get jacobian determinant from parent to physical space
    double GeometryElement::jacDet(const double xi,
                                   const double eta,
                                   const double zeta) const
    {
        auto it = mJDetCache.find(std::make_tuple(xi,eta,zeta));
        if(it != mJDetCache.end())
            return it->second;
        ParamCoord c = getParamCoord(xi,eta,zeta);
        const double jdet =  geometry().jacDet(c.u, c.v, c.w, spaceI()) * jacDetParam(xi,eta,zeta);
        insertCachedJDet(std::make_tuple(xi,eta,zeta), jdet);
        return jdet;
    }
    
    Point3D GeometryElement::eval(const double xi,
                                  const double eta,
                                  const double zeta) const
    {
        auto it = mPointCache.find(std::make_tuple(xi,eta,zeta));
        if(it != mPointCache.end())
            return it->second;
        
        ParamCoord c = getParamCoord(xi,eta,zeta);
        const auto pt = geometry().eval(c.u,c.v,c.w, spaceI());
        insertCachedPoint(std::make_tuple(xi,eta,zeta), pt);
        return pt;
    }
    
    Point3D GeometryElement::evalVertex(const Vertex v) const
    {
        switch (v) {
            case Vertex::VERTEX0:
                return eval(-1.0, -1.0, -1.0);
                break;
            case Vertex::VERTEX1:
                return eval(1.0,-1.0, -1.0);
                break;
            case Vertex::VERTEX2:
                return eval(-1.0, -1.0, 1.0);
                break;
            case Vertex::VERTEX3:
                return eval(1.0, -1.0, 1.0);
                break;
            case Vertex::VERTEX4:
                return eval(-1.0, 1.0, -1.0);
                break;
            case Vertex::VERTEX5:
                return eval(1.0, 1.0, -1.0);
                break;
            case Vertex::VERTEX6:
                return eval(-1.0, 1.0, 1.0);
                break;
            case Vertex::VERTEX7:
                return eval(1.0, 1.0, 1.0);
                break;
            default:
                error("Bad vertex enum");
                break;
        }
    }
    
    /// Get jacobian
    DoubleVecVec GeometryElement::jacob(const double xi,
                                        const double eta,
                                        const double zeta) const
    {
        auto it = mJacobCache.find(std::make_tuple(xi,eta,zeta));
        if(it != mJacobCache.end())
            return it->second;
        
        ParamCoord c = getParamCoord(xi,eta,zeta);
        
        const auto jacob_param = geometry().jacob(c.u, c.v, c.w, spaceI());
        const auto jacob_parent = jacobParam(xi,eta,zeta);
        
        DoubleVecVec r{
            { 0.0, 0.0, 0.0 },
            { 0.0, 0.0, 0.0 },
            { 0.0, 0.0, 0.0 }
        };
        
        // multiply jacobian matrices above to generate final jacobian
        for(uint i = 0; i < 3; ++i)
            for(uint j = 0; j < 3; ++j)
                for(uint k = 0; k < 3; ++k)
                    r[i][j] += jacob_param[i][k] * jacob_parent[k][j];
        
        insertCachedJacob(std::make_tuple(xi,eta,zeta), r);
        return r;
    }
    
    /// Get the normal
//    Point3D GeometryElement::normal(const double xi,
//                                    const double eta,
//                                    const double zeta) const
//    {
//        ParamCoord c = paramCoord(u,v);
//        return geometry().normal(c.s, c.t, spaceI());
//    }
    
    /// Get the tangent
    Point3D GeometryElement::tangent(const double xi,
                                     const double eta,
                                     const double zeta,
                                     const ParamDir d) const
    {
        ParamCoord c = getParamCoord(xi,eta,zeta);
        return geometry().tangent(c.u, c.v, c.w, spaceI(), d);
    }
    
    UIntVec GeometryElement::geometryDegree() const
    { return geometry().degree(spaceI()); }
    
    std::pair<bool, GPt3D> GeometryElement::containsParamPt(const double u,
                                                            const double v,
                                                            const double w) const
    {
        if(boost::icl::contains(boostInterval(U),u) && boost::icl::contains(boostInterval(V), v) && boost::icl::contains(boostInterval(W), w))
            return std::make_pair(true, getParentCoord(u,v,w));
        else
            return std::make_pair(false, GPt3D());
    }
    
    bool GeometryElement::contains(const GeometryElement& e) const
    {
        if(containsParamPt(e.lowerBound()).first && containsParamPt(e.upperBound()).first)
            return true;
        return false;
    }
    
    void GeometryElement::print(std::ostream& ost) const
    {
        ost << "Geometry: " << &mGeom << "\n";
        ost << "Space index: " << mSpaceI << "\n";
        ost << "Knot intervals: ";
        for(const auto& i : mKnotIntervals)
            ost << "(" << i.lower() << ", " << i.upper() << ")  ";
        
    }
    
    double dist(const GeometryElement& e1, const GeometryElement& e2)
    {
        return dist(e1.eval(0.0,0.0,0.0), e2.eval(0.0, 0.0, 0.0));
    }
    
}