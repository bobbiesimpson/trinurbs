//
//  Geometry.cpp
//  trinurbs
//
//  Created by Robert Simpson on 14/07/2014.
//
//

#include <ios>
#include <memory>
#include <string>
#include <numeric>

#include "Geometry.h"
#include "BSplineSpace.h"
#include "Forest.h"
#include "NURBSCommon.h"
#include "base.h"

//#include "vtkSmartPointer.h"
//#include "vtkUnstructuredGrid.h"
//#include "vtkDoubleArray.h"
//#include "vtkPoints.h"
//#include "vtkXMLUnstructuredGridWriter.h"
//#include "vtkQuad.h"
//#include "vtkPointData.h"

namespace trinurbs
{
    UIntVec Geometry::degree(const uint ispace) const
    {
        return primalForest().space(ispace).degree();
    }

    double Geometry::jacDet(const double u,
                            const double v,
                            const double w,
                            const uint ispace) const
    {
        const auto j = jacob(u,v,w,ispace);
        return det3x3(j);
    }
    
    /// Jacobian
    /// [ dx/du dx/dv dx/dw
    ///   dy/du dy/dv dy/dw
    ///   dz/du dz/dv dz/dw ]
    DoubleVecVec Geometry::jacob(const double u,
                                 const double v,
                                 const double w,
                                 const uint ispace) const
    {
        // first generate the tranpose of the jacobian
        // from tangent vectors
        auto jtemp = DoubleVecVec{
            tangent(u,v,w,ispace,U).asVec(),
            tangent(u,v,w,ispace,V).asVec(),
            tangent(u,v,w,ispace,W).asVec()
        };
        
        // transpose
        transpose(jtemp);
        
        return jtemp;
    }
    
    /// Normal
//    Point3D Geometry::normal(const double u,
//                             const double v,
//                             const double w,
//                             const uint ispace) const
//    {
//        Point3D n = cross(tangent(s,t,sp,S), tangent(s,t,sp,T));
//        if(normalsFlipped())
//            n *= -1;
//        return n / n.length();
//    }
    
    /// 'Tangent' vector at a given parametric coordinate
    Point3D Geometry::tangent(const double u,
                              const double v,
                              const double w,
                              const uint ispace,
                              const ParamDir dir) const
    {
        /// The two 'other' parametric directions
        ParamDir dir2, dir3;
        
        switch(dir)
        {
            case U:
                dir2 = V;
                dir3 = W;
                break;
            case V:
                dir2 = U;
                dir3 = W;
                break;
            case W:
                dir2 = U;
                dir3 = V;
                break;
            default:
                error("Bad parametric direction in Geometry::Tangent()");
        }
        
        Point4D temp;
        
        // Get basis functions, derivatives, indices and span
        const BSplineSpace& b_space = primalForest().space(ispace);
        UIntVecVec indices = b_space.localBasisIVec(u,v,w);
        UIntVec span = b_space.span(u,v,w);
        DoubleVecVec basis = b_space.tensorBasis(u,v,w,span);
        DoubleVecVec ders = b_space.tensorBasisDers(u,v,w,span);
        
        Point3D a, ader;
        double wt = 0.0, wder = 0.0;
        
        // now loop over tensor product structure of space
        for(uint i = 0; i < b_space.degree(U) + 1; ++i )
        {
            for(uint j = 0; j < b_space.degree(V) + 1; ++j )
            {
                for(uint k = 0; k < b_space.degree(W) + 1; ++k)
                {
                    const uint der_i = (dir == U) ? i : ( (dir == V) ? j : k);
                    const uint nder2_i = (dir == U) ? j : ( (dir == V) ? i : i);
                    const uint nder3_i = (dir == U) ? k : ( (dir == V) ? k : j);
                    
                    const Point4D& p = point(ispace, indices[U][i], indices[V][j], indices[W][k]);
                    
                    a += p.asUnweighted() * basis[U][i] * basis[V][j] * basis[W][k];
                    wt += p.getWeight() * basis[U][i] * basis[V][j] * basis[W][k];
                    
                    ader += p.asUnweighted() * ders[dir][der_i] * basis[dir2][nder2_i] * basis[dir3][nder3_i];
                    wder += p.getWeight() * ders[dir][der_i] * basis[dir2][nder2_i] * basis[dir3][nder3_i];
                }
            }
        }
        //        std::cout << "w = " << w << "\n\n";
        //        std::cout << "wder = " << wder << "\n\n";
        
        return nurbshelper::getNonRationalDeriv({ a, ader }, { w, wder });
    }
    
    /// Interpolate the surface
    Point3D Geometry::eval(const double u,
                           const double v,
                           const double w,
                           const uint ispace) const
    {
        const BSplineSpace& b_space = primalForest().space(ispace);
        UIntVec indices = b_space.globalBasisIVec(u,v,w);
        DoubleVec basis = b_space.basis(u,v,w);
        
        assert(basis.size() == indices.size());
        
        Point4D p;
        for(uint i = 0; i < basis.size(); ++i)
            p += point(ispace,indices[i]) * basis[i];
        return p.asCartesian();
    }
//
//    /// Write to vtu file.
//    void Geometry::writeVTKOutput(const std::string& file,
//                                  const uint nsample) const
//    {
//        /// TODO
//    }
//    
    bool Geometry::load(std::istream& ist)
    {
        char ch;
        if(ist >> ch && ch != '{') {
            ist.unget();
            ist.clear(std::ios_base::failbit);
            return false;
        }
        clear();
        Forest f_load;
        if(!(ist >> f_load))
            error("Cannot load forest into geometry data structure");
        mPrimalForest = f_load;
        
        if(ist >> ch && ch != '{') { // We might be reading something else
            ist.unget();
            ist.clear(std::ios_base::failbit);
            return false;
        }
        mCPts.clear();
        Point4D p;
        while(ist >> p)
            mCPts.push_back(p);
        endOfLoop(ist, '}', "Bad control point set");
        return true;
    }

    bool Geometry::loadHBSFile(std::istream& ist)
    {
        if(!ist)
            error("Cannot load HBS file: bad istream");
        clear();
        
        uint d; // get the dimension
        std::string s;
        std::getline(ist, s);
        std::istringstream ifs(s);
        ifs >> d;
        if(d != 3)
            error("Currently only parametric dimension = 3 is implemented ");
        
        uint npts;
        std::getline(ist, s);
        ifs.clear(); ifs.str(s);
        ifs >> npts;
        
        uint ntrunk;
        std::getline(ist, s);
        ifs.clear(); ifs.str(s);
        ifs >> ntrunk;
        
        std::vector<Point4D> cp_vec;
        for(uint i = 0; i < npts; ++i) {
            std::getline(ist, s);
            ifs.clear(); ifs.str(s);
            double x,y,z,w;
            ifs >> x >> y >> z >> w;
            cp_vec.emplace_back(x,y,z,w);
        }
        
        std::vector<BSplineSpace> s_vec;
        std::vector<UIntVec> conn_vec;
        for(uint i = 0; i < ntrunk; ++i) {
            UIntVec degree;
            std::getline(ist, s);
            ifs.clear(); ifs.str(s);
            for(uint c = 0; c < d; ++c) {
                uint v;
                ifs >> v;
                degree.push_back(v);
            }
            
            UIntVec point_vec;
            std::getline(ist, s);
            ifs.clear(); ifs.str(s);
            for(uint c = 0; c < d; ++c) {
                uint v;
                ifs >> v;
                point_vec.push_back(v);
            }
            uint n = 1;
            std::for_each(point_vec.begin(), point_vec.end(), [&n](uint& val){ n *= val; });
            
            std::vector<DoubleVec> knot_vec(degree.size());
            for(uint c = 0; c < point_vec.size(); ++c) {
                std::getline(ist, s);
                ifs.clear(); ifs.str(s);
                double kval;
                while(ifs >> kval)
                    knot_vec[c].push_back(kval);
            }
            s_vec.emplace_back(knot_vec, degree);
            
            UIntVec conn;
            for(uint c = 0; c < n; ++c) {
                std::getline(ist, s);
                ifs.clear(); ifs.str(s);
                uint val;
                ifs >> val;
                conn.push_back(val);
            }
            conn_vec.push_back(conn);
        }
        mPrimalForest = Forest(s_vec, conn_vec);
        mPrimalForest.setGeometry(this);
        mCPts = cp_vec;
        return true;
    }
//
//    
    void Geometry::print(std::ostream& ost) const
    {
        ost << "- Geometry object -\n\n";
        ost << "Primal forest: \n\n" << mPrimalForest << "\n";
        ost << "Control point set: \n\n";
        for(const auto& p : mCPts)
            std::cout << p << "\n";
        std::cout << "\n- end geometry object -\n";
    }
//
//    void Geometry::rescale(const double sf)
//    {
//        for(uint ipoint = 0; ipoint < controlPtN(); ++ipoint) {
//            auto& pt = point(ipoint);
//            for(uint i = 0; i < 3; ++i)
//                pt[i] *= sf;
//        }
//    }
//    
//    void Geometry::rotate(const nurbs::CartesianComponent comp,
//                          const double theta)
//    {
//        // make a copy of all control points
//        auto cvec_copy = mCPts;
//        
//        // compute rotation matrix depending on given componet
//        DoubleVecVec rotmat;
//        switch (comp) {
//            case nurbs::CartesianComponent::X:
//                rotmat = DoubleVecVec
//            {
//                {1.0, 0.0, 0.0},
//                {0.0, std::cos(theta), -std::sin(theta)},
//                {0.0, std::sin(theta), std::cos(theta)}
//            };
//                break;
//            case nurbs::CartesianComponent::Y:
//                rotmat = DoubleVecVec
//            {
//                {std::cos(theta),0.0, std::sin(theta)},
//                {0.0, 1.0, 0.0},
//                {-std::sin(theta), 0.0, std::cos(theta)}
//            };
//                break;
//            case nurbs::CartesianComponent::Z:
//                rotmat = DoubleVecVec
//            {
//                {std::cos(theta), -std::sin(theta), 0.0},
//                {std::sin(theta), std::cos(theta), 0.0},
//                {0.0, 0.0, 1.0}
//            };
//                break;
//        }
//        
//        // loop over all points and apply rotation matrix to each
//        for(auto& pt : cvec_copy)
//        {
//            const auto pt_copy = pt;
//            Point4D newpt;
//            newpt.setWeight(pt.getWeight());
//            for(uint i = 0; i < 3; ++i)
//            {
//                for(uint j = 0; j < 3; ++j)
//                    newpt[i] += rotmat[i][j] * pt_copy[j];
//            }
//            pt = newpt;
//        }
//        
//        // finally copy rotated points to member variable
//        mCPts = cvec_copy;
//    }
//    
//    void Geometry::translate(const nurbs::Point3D& p)
//    {
//        auto cpvec_copy = mCPts;
//        
//        for(auto& cpt : cpvec_copy)
//        {
//            cpt[0] += p[0];
//            cpt[1] += p[1];
//            cpt[2] += p[2];
//        }
//        mCPts = cpvec_copy;
//    }
//    
//    void Geometry::normalise()
//    {
//        const double max = std::numeric_limits<double>::max();
//        const double min = std::numeric_limits<double>::lowest();
//        Point3D minpt(max, max, max);
//        Point3D maxpt(min, min, min);
//        
//        for(unsigned ipoint = 0; ipoint < controlPtN(); ++ipoint)
//        {
//            const auto& cpt = point(ipoint);
//            for(unsigned i = 0; i < 3; ++i)
//            {
//                if(cpt[i] < minpt[i])
//                    minpt[i] = cpt[i];
//                if(cpt[i] > maxpt[i])
//                    maxpt[i] = cpt[i];
//            }
//        }
//        const auto diff = maxpt - minpt;
//        std::vector<double> diffvec = diff.asVec();
//        auto result = std::max_element(diffvec.begin(), diffvec.end());
//        const double scalefactor = 1.0 / (*result);
//        rescale(scalefactor);
//        
//    }
    
    std::istream& operator>>(std::istream& ist, Geometry& g)
    {
        g.load(ist);
        return ist;
    }
    
    /// Overload outut operator
    std::ostream& operator<<(std::ostream& ost, const Geometry& g)
    {
        g.print(ost);
        return ost;
    }
}