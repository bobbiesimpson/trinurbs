#include <algorithm>
#include <stdexcept>
#include <numeric>

#include "BSplineSpace.h"
#include "base.h"
#include "NURBSCommon.h"

namespace trinurbs
{
    
    DoubleVecVec BSplineSpace::tensorBasis(const double u,
                                           const double v,
                                           const double w,
                                           const UIntVec& span) const
    {
        return {nurbshelper::getBsplineBasis(u, span[U], knotVec(U), degree(U)),
                nurbshelper::getBsplineBasis(v, span[V], knotVec(V), degree(V)),
                nurbshelper::getBsplineBasis(w, span[W], knotVec(W), degree(W))
        };
    }
    
    DoubleVecVec BSplineSpace::tensorBasisDers(const double u,
                                               const double v,
                                               const double w,
                                               const UIntVec& span,
                                               const DerivOrder der) const
    {
        return {nurbshelper::getBsplineBasisDers(u, span[U], knotVec(U), degree(U), der).at(der),
                nurbshelper::getBsplineBasisDers(v, span[V], knotVec(V), degree(V), der).at(der),
                nurbshelper::getBsplineBasisDers(w, span[W], knotVec(W), degree(W), der).at(der)};
    }
    
    UIntVecVec BSplineSpace::localBasisFuncI(const double u,
                                             const double v,
                                             const double w) const
    {
        return {nurbshelper::getBasisFnIndices(u, knotVec(U), degree(U)),
                nurbshelper::getBasisFnIndices(v, knotVec(V), degree(V)),
                nurbshelper::getBasisFnIndices(w, knotVec(W), degree(W))
        };
    }
    
    UIntVec BSplineSpace::span(const double u,
                               const double v,
                               const double w) const
    {
        return {nurbshelper::getKnotSpan(u, knotVec(U), degree(U)),
                nurbshelper::getKnotSpan(v, knotVec(V), degree(V)),
                nurbshelper::getKnotSpan(w, knotVec(W), degree(W))};
    }
    
    std::tuple<uint, uint, uint> BSplineSpace::localIndices(const uint iel) const
    {
        const uint nu = uniqueKnotN(U) - 1;
        const uint nv = uniqueKnotN(V) - 1;
        
        return localElementSpaceIndices(iel, nu, nv);
    }
    
    
    void BSplineSpace::degreeReduce(const ParamDir dir)
    {
        if(0 == mDegrees[dir])
            throw std::runtime_error("Cannot degree reduce degree 0 spline\n");
        auto& knotvec = mKnotVecs[dir];
        std::vector<double> kv_reduced;
        std::copy(knotvec.begin() + 1, knotvec.end() - 1, std::back_inserter(kv_reduced));
        knotvec = kv_reduced;
        
        //        auto it = knotvec.begin();
        //        for(uint i = 0; i < uniqueKnotN(dir); ++i) {
        //            const double knot = uniqueKnot(i, dir);
        //            it = std::find(it, knotvec.end(), knot);
        //            knotvec.erase(it); // only remove first and last knots
        //        }
        mDegrees[dir] -= 1; // reduce degree
        init();
    }
    
    void BSplineSpace::degreeElevate(const ParamDir dir)
    {
        auto& knotvec = mKnotVecs[dir];
        const auto lastknot = knotvec.back();
        const auto firstknot = knotvec.front();
        knotvec.push_back(lastknot);
        knotvec.push_back(firstknot);
        
        //        for(uint i = 0; i < uniqueKnotN(dir); ++i) {
        //            const double knot = uniqueKnot(i, dir);
        //            knotvec.push_back(knot);
        //        }
        mDegrees[dir] += 1; // elevate degree
        std::sort(knotvec.begin(), knotvec.end());
        init();
    }
    
    void BSplineSpace::hrefine(const uint n)
    {
        for(auto& kv : mKnotVecs)
            kv = nurbshelper::uniformKnotInsertion(kv,n);
        init();
    }
    
    void BSplineSpace::graded_hrefine(const uint n, const double coeff)
    {
        for(auto& kv : mKnotVecs)
            kv = nurbshelper::gradedKnotInsertion(kv, n, coeff);
        init();
    }
    
    void BSplineSpace::load(std::istream& ist)
    {
        clear();
        char ch;
        if(ist >> ch && ch != '{') {
            ist.unget();
            ist.clear(std::ios_base::failbit);
            return;
        }
        std::string s;
        ist >> s;
        if(!ist )
            error("Cannot read b-spline space");
        mName = s;
        InputVec<uint> d_vec;
        if(!(ist >> d_vec))
            error("Bad degree vector");
        mDegrees = d_vec.data;
        if(ist >> ch && ch != '{')
            error("Bad knot vector reading");
        while(true) {
            InputVec<double> k_vec;
            if(!(ist >> k_vec))
                break;
            mKnotVecs.push_back(k_vec.data);
        }
        endOfLoop(ist, '}', "Bad knot vector.");
        init(); // recalculate unique knots + intervals
    }
    
    void BSplineSpace::init()
    {
        mUniqueKnotVecs.clear(); // clear data
        mIntervalVec.clear();
        mNumBasisVec.clear();
        mUExtractionOperators.clear();
        mVExtractionOperators.clear();
        
        for(uint i = 0; i < paramDimN(); ++i) {
            mIntervalVec.push_back(boost::icl::construct<Interval>
                                   (mKnotVecs[i].front(), mKnotVecs[i].back(),
                                    boost::icl::interval_bounds::closed()));
            Interval& last = mIntervalVec.back();
            if(essentiallyEqual(last.upper(), last.lower(), TOL))
                error("Cannot prescribe a b-spline space with non-zero parametric area.");
            mNumBasisVec.push_back(mKnotVecs[i].size() - mDegrees[i] - 1);
        }
        
        // construct unique knot vectors
        for(auto kv : mKnotVecs) {
            auto last = std::unique(kv.begin(), kv.end());
            mUniqueKnotVecs.emplace_back( DoubleVec( kv.begin(), last ) );
        }
        
        typedef std::vector<std::vector<double>> Matrix;
        
        // construct element extraction operators
        for(uint iparam = 0; iparam < paramDimN(); ++iparam)
        {
            const ParamDir dir = ParamDirType(iparam);
            const uint p = degree(dir);
            const auto kvec = knotVec(dir);
            const auto uniqueknotvec = uniqueKnotVec(dir);
            
            size_t currentindex = 0;
            size_t insertedknot_n = 0;
            
            Matrix global_cmatrix;  // the global extraction operator
            auto kvec_copy = kvec;
            
            const auto n = kvec_copy.size() - p - 1;
            
            // keep going until we have C^0 continuity at all knots
            while(currentindex < uniqueknotvec.size())
            {
                const auto kval = uniqueknotvec[currentindex];
                const auto count = std::count_if(kvec_copy.begin(), kvec_copy.end(), [&](double k) {return essentiallyEqual(k, kval, 1.0e-5);});
                auto requiredknots = p - count;
                
                // perform required knot insertion
                while(requiredknots > 0)
                {
                    const unsigned k = nurbshelper::getKnotSpan(kval, kvec_copy, p);
                    //                    std::cout << requiredknots << " knots required at " << kval << " at knot span " <<  k << "\n";
                    
                    // initialise Cmatrix
                    Matrix Cmatrix;
                    for(size_t i = 0; i < n + insertedknot_n; ++i)
                        Cmatrix.push_back(std::vector<double>(n + insertedknot_n + 1, 0.0));
                    
                    for(unsigned a = 0; a < n + insertedknot_n + 1; ++a)
                    {
                        // First compute alpha
                        double alpha;
                        if(a <= k - p - 1)
                            alpha = 1.0;
                        else if(a >= k - p && a <= k - 1)
                        {
                            double numer = kval - kvec_copy[a];
                            alpha = numer / (kvec_copy[a + p] - kvec_copy[a]);
                        }
                        else
                            alpha = 0.0;
                        
                        // now populate the matrix for the current knot insertion
                        if(a == 0)
                            Cmatrix[a][a] = alpha;
                        else if(a == n + insertedknot_n)
                            Cmatrix[a-1][a] = 1 - alpha;
                        else
                        {
                            Cmatrix[a][a] = alpha;
                            Cmatrix[a-1][a] = 1- alpha;
                        }
                    }
                    
                    if(insertedknot_n == 0)
                        global_cmatrix = Cmatrix;
                    else
                    {
                        // now update the global matrix
                        auto global_cmatrix_copy = global_cmatrix;
                        
                        global_cmatrix.clear();
                        // resize global matrix
                        for(size_t i = 0; i < n; ++i)
                            global_cmatrix.push_back(std::vector<double>(n + insertedknot_n + 1, 0.0));
                        
                        for(size_t i = 0; i < n; ++i)
                            for(size_t j = 0; j < n + insertedknot_n + 1; ++j)
                                for(size_t k = 0; k < n + insertedknot_n; ++k)
                                    global_cmatrix[i][j] += global_cmatrix_copy[i][k] * Cmatrix[k][j];
                    }
                    // insert the new knot into the knot vector
                    kvec_copy.insert(kvec_copy.begin() + k , kval);
                    requiredknots--;
                    ++insertedknot_n;
                }
                ++currentindex;
            }
            
            // in the case that no knots have been inserted, the knot vector must already
            // be in bezier form and we simply set the identity matrix for each extraction
            // operator. Otherwise, use the global Cmatrix to construct element extraciton
            // operators
            
            // We've inserted knots
            if(insertedknot_n != 0)
            {
                // now assign the element extraction operators
                for(uint iel = 0; iel < uniqueknotvec.size() - 1; ++iel)
                {
                    const auto rows = nurbshelper::getBasisFnIndices(uniqueknotvec[iel], kvec, p);
                    UIntVec cols(p+1);
                    std::iota(cols.begin(), cols.end(), p * iel);
                    
                    Matrix extraction_op;
                    for(size_t r = 0; r < rows.size(); ++r)
                        extraction_op.push_back(std::vector<double>(cols.size(), 0.0));
                    
                    for(size_t i = 0; i < rows.size(); ++i)
                        for(size_t j = 0; j < cols.size(); ++j)
                            extraction_op[i][j] = global_cmatrix[rows[i]][cols[j]];
                    
                    setExtractionOperator(iel, dir, extraction_op);
                }
            }
            // no knots inserted
            else
            {
                for(size_t iel = 0; iel < uniqueknotvec.size() - 1; ++iel)
                {
                    Matrix I(p + 1, std::vector<double>(p+1, 0.0));
                    for(uint i = 0; i < p + 1; ++i)
                        I[i][i] = 1.0;
                    setExtractionOperator(iel, dir, I);
                }
            }
        }
        
    }
    
    void BSplineSpace::printData(std::ostream& ost) const
    {
        ost << "--------- NURBS SPACE DATA -------------\n\n";
        for(uint i = 0; i < paramDimN(); ++i) {
            ost << "Parameteric direction: " << i << "\n";
            ost << "Degree: " << mDegrees[i] << "\n";
            ost << "Knot vector: ";
            printVector(mKnotVecs[i], ost);
            ost << "\n";
            ost << "Unique knot vector: ";
            printVector(mUniqueKnotVecs[i], ost);
            ost << "\n---------\n";
            
        }
        ost << "\n--------------------------------------\n";				
    }
    
    std::tuple<uint, uint, uint> localElementSpaceIndices(const uint ielem,
                                                          const uint nel_u,
                                                          const uint nel_v)
    {
        const uint nuv = nel_u * nel_v;
        const uint k = ielem / nuv;
        
        return std::make_tuple(ielem % nel_u, (ielem - k * nuv) / nel_u, k);
    }
    
    uint localBasisI(const Vertex v,
                     const uint nb_u,
                     const uint nb_v,
                     const uint nb_w)
    {
        switch(v) {
            case Vertex::VERTEX0: return 0;
            case Vertex::VERTEX1: return nb_u - 1;
            case Vertex::VERTEX2: return nb_u * nb_v * (nb_w - 1);
            case Vertex::VERTEX3: return nb_u * nb_v * (nb_w - 1) + nb_u - 1;
            case Vertex::VERTEX4: return nb_u * (nb_v - 1);
            case Vertex::VERTEX5: return nb_u * nb_v - 1;
            case Vertex::VERTEX6: return nb_u * nb_v * nb_w - nb_u;
            case Vertex::VERTEX7: return nb_u * nb_v * nb_w - 1;
        }
    }
    
    uint localBasisI(const Vertex v,
                     const BSplineSpace& s)
    { return localBasisI(v, s.basisFuncN(U),
                            s.basisFuncN(V),
                            s.basisFuncN(W)); }
    
    UIntVec localBasisIVec(const Edge e,
                           const uint nb_u,
                           const uint nb_v,
                           const uint nb_w)
    {
        std::vector<uint> l_vec;
        
        switch(e) {
            case Edge::EDGE0:
                for(uint i = 0; i < nb_u; ++i)
                    l_vec.push_back(i);
                break;
            case Edge::EDGE1:
                for(uint i = 0; i < nb_u; ++i)
                    l_vec.push_back(i + nb_u * nb_v * (nb_w - 1));
                break;
            case Edge::EDGE2:
                for(uint i = 0; i < nb_w; ++i)
                    l_vec.push_back(i * nb_u * nb_v);
                break;
            case Edge::EDGE3:
                for(uint i = 0; i < nb_w; ++i)
                    l_vec.push_back(i * nb_u * nb_v + (nb_u - 1));
                break;
            case Edge::EDGE4:
                for(uint i = 0; i < nb_u; ++i)
                    l_vec.push_back(i + nb_u * (nb_v - 1));
                break;
            case Edge::EDGE5:
                for(uint i = 0; i < nb_u; ++i)
                    l_vec.push_back(i + nb_u * nb_v * nb_w - nb_u);
                break;
            case Edge::EDGE6:
                for(uint i = 0; i < nb_w; ++i)
                    l_vec.push_back((i + 1) * nb_u * nb_v - nb_u);
                break;
            case Edge::EDGE7:
                for(uint i = 0; i < nb_w; ++i)
                    l_vec.push_back((i + 1) * nb_u * nb_v - 1);
                break;
            case Edge::EDGE8:
                for(uint i = 0; i < nb_v; ++i)
                    l_vec.push_back(i * nb_u);
                break;
            case Edge::EDGE9:
                for(uint i = 0; i < nb_v; ++i)
                    l_vec.push_back(i * nb_u + (nb_u - 1));
                break;
            case Edge::EDGE10:
                for(uint i = 0; i < nb_v; ++i)
                    l_vec.push_back(i * nb_u + nb_u * nb_v * (nb_w - 1));
                break;
            case Edge::EDGE11:
                for(uint i = 0; i < nb_v; ++i)
                    l_vec.push_back(i * nb_u + nb_u * nb_v * (nb_w - 1) + (nb_u - 1));
                break;
        }
        
        return l_vec;
    }
    
    UIntVec localBasisIVec(const Edge e, const BSplineSpace& s)
    {
        return localBasisIVec(e, s.basisFuncN(U),
                                 s.basisFuncN(V),
                                 s.basisFuncN(W));
    }
    
    std::pair<uint, uint> localBasisIPair(const Edge e,
                                          const BSplineSpace& s)
    {
        switch (e) {
            case Edge::EDGE0:
                return std::make_pair(localBasisI(Vertex::VERTEX0, s),
                                      localBasisI(Vertex::VERTEX1, s));
                break;
            case Edge::EDGE1:
                return std::make_pair(localBasisI(Vertex::VERTEX2, s),
                                      localBasisI(Vertex::VERTEX3, s));
                break;
            case Edge::EDGE2:
                return std::make_pair(localBasisI(Vertex::VERTEX0, s),
                                      localBasisI(Vertex::VERTEX2, s));
                break;
            case Edge::EDGE3:
                return std::make_pair(localBasisI(Vertex::VERTEX1, s),
                                      localBasisI(Vertex::VERTEX3, s));
                break;
            case Edge::EDGE4:
                return std::make_pair(localBasisI(Vertex::VERTEX4, s),
                                      localBasisI(Vertex::VERTEX5, s));
                break;
            case Edge::EDGE5:
                return std::make_pair(localBasisI(Vertex::VERTEX6, s),
                                      localBasisI(Vertex::VERTEX7, s));
                break;
            case Edge::EDGE6:
                return std::make_pair(localBasisI(Vertex::VERTEX4, s),
                                      localBasisI(Vertex::VERTEX6, s));
                break;
            case Edge::EDGE7:
                return std::make_pair(localBasisI(Vertex::VERTEX5, s),
                                      localBasisI(Vertex::VERTEX7, s));
                break;
            case Edge::EDGE8:
                return std::make_pair(localBasisI(Vertex::VERTEX0, s),
                                      localBasisI(Vertex::VERTEX4, s));
                break;
            case Edge::EDGE9:
                return std::make_pair(localBasisI(Vertex::VERTEX1, s),
                                      localBasisI(Vertex::VERTEX5, s));
                break;
            case Edge::EDGE10:
                return std::make_pair(localBasisI(Vertex::VERTEX2, s),
                                      localBasisI(Vertex::VERTEX6, s));
                break;
            case Edge::EDGE11:
                return std::make_pair(localBasisI(Vertex::VERTEX3, s),
                                      localBasisI(Vertex::VERTEX7, s));
                break;
        }
    }
    
    std::tuple<uint, uint, uint, uint> localBasisITuple(const Face f,
                                                        const BSplineSpace& s)
    {
        switch(f)
        {
            case Face::FACE0:
                return std::make_tuple(localBasisI(Vertex::VERTEX0, s),
                                       localBasisI(Vertex::VERTEX1, s),
                                       localBasisI(Vertex::VERTEX3, s),
                                       localBasisI(Vertex::VERTEX2, s));
                break;
            case Face::FACE1:
                return std::make_tuple(localBasisI(Vertex::VERTEX4, s),
                                       localBasisI(Vertex::VERTEX5, s),
                                       localBasisI(Vertex::VERTEX7, s),
                                       localBasisI(Vertex::VERTEX6, s));
                break;
            case Face::FACE2:
                return std::make_tuple(localBasisI(Vertex::VERTEX0, s),
                                       localBasisI(Vertex::VERTEX2, s),
                                       localBasisI(Vertex::VERTEX6, s),
                                       localBasisI(Vertex::VERTEX4, s));
                break;
            case Face::FACE3:
                return std::make_tuple(localBasisI(Vertex::VERTEX1, s),
                                       localBasisI(Vertex::VERTEX3, s),
                                       localBasisI(Vertex::VERTEX7, s),
                                       localBasisI(Vertex::VERTEX5, s));
                break;
            case Face::FACE4:
                return std::make_tuple(localBasisI(Vertex::VERTEX0, s),
                                       localBasisI(Vertex::VERTEX1, s),
                                       localBasisI(Vertex::VERTEX5, s),
                                       localBasisI(Vertex::VERTEX4, s));
                break;
            case Face::FACE5:
                return std::make_tuple(localBasisI(Vertex::VERTEX2, s),
                                       localBasisI(Vertex::VERTEX3, s),
                                       localBasisI(Vertex::VERTEX7, s),
                                       localBasisI(Vertex::VERTEX6, s));
                break;

        }
    }
    
    UIntVecVec localBasisIVec(const Face f, const BSplineSpace& s)
    {
        return localBasisIVec(f, s.basisFuncN(U),
                                 s.basisFuncN(V),
                                 s.basisFuncN(W));
    }
    
    /// Return the set of local basis function indices for a given face
    /// using the coordinate system as defined in the comments in base.h
    UIntVecVec localBasisIVec(const Face f,
                              const uint nb_u,
                              const uint nb_v,
                              const uint nb_w)
    {
        UIntVecVec lvec;
        switch(f)
        {
            case Face::FACE0:
                // Origin: 0
                // u axis: 0-1
                // v axis 0-2
                for(uint iw = 0; iw < nb_w; ++iw)
                {
                    UIntVec temp;
                    for(uint iu = 0; iu < nb_u; ++iu)
                        temp.push_back(iu + nb_u * nb_v * iw);
                    lvec.push_back(temp);
                }
                break;
            case Face::FACE1:
                // Origin: 4
                // u axis: 4-5
                // v axis 4-6
                for(uint iw = 0; iw < nb_w; ++iw)
                {
                    UIntVec temp;
                    for(uint iu = 0; iu < nb_u; ++iu)
                        temp.push_back(iu + nb_u * (nb_v - 1) + nb_u * nb_v * iw);
                    lvec.push_back(temp);
                }
                break;
            case Face::FACE2:
                // Origin: 0
                // u axis: 0-2
                // v axis 0-4
                for(uint iv = 0; iv < nb_v; ++iv)
                {
                    UIntVec temp;
                    for(uint iw = 0; iw < nb_w; ++iw)
                        temp.push_back(iw * nb_u * nb_v + iv * nb_u);
                    lvec.push_back(temp);
                }
                break;
            case Face::FACE3:
                // Origin: 1
                // u axis: 1-3
                // v axis 1-5
                for(uint iv = 0; iv < nb_v; ++iv)
                {
                    UIntVec temp;
                    for(uint iw = 0; iw < nb_w; ++iw)
                        temp.push_back(iw * nb_u * nb_v + iv * nb_u + (nb_u - 1));
                    lvec.push_back(temp);
                }
                break;
            case Face::FACE4:
                // Origin: 0
                // u axis: 0-1
                // v axis 0-4
                for(uint iv = 0; iv < nb_v; ++iv)
                {
                    UIntVec temp;
                    for(uint iu = 0; iu < nb_u; ++iu)
                        temp.push_back(iv * nb_u + iu);
                    lvec.push_back(temp);
                }
                break;
            case Face::FACE5:
                // Origin: 2
                // u axis: 2-3
                // v axis 2-6
                for(uint iv = 0; iv < nb_v; ++iv)
                {
                    UIntVec temp;
                    for(uint iu = 0; iu < nb_u; ++iu)
                        temp.push_back(iv * nb_u + iu + (nb_w - 1) * nb_u * nb_v);
                    lvec.push_back(temp);
                }
                break;
        }
        return lvec;
    }
    
    void reorderLocalFaceIndices(std::vector<std::vector<uint>>& lindices,
                                 const Face f,
                                 const std::tuple<uint, uint, uint, uint>& igvert)
    {
        // general algorithm is to first consider the set igvert as cyclic
        // and to rotate the face until the lowest index is the first index.
        // The number of rotations required to do this dicates how many
        // rotations to apply to the lindices matrix
        
        // See this page for how transpose/rotate matrix
        // https://stackoverflow.com/questions/42519/how-do-you-rotate-a-two-dimensional-array
        
        // And finally, if required apply a transpose to ensure u-axis
        // is aligned along the bottom edge
        
        
        
        
        
    }
    
    
}
