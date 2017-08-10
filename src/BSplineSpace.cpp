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
    
}
