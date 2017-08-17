#ifndef TRINURBS_IPARENT_SAMPLE_H
#define TRINURBS_IPARENT_SAMPLE_H

#include "base.h"

namespace trinurbs
{
    struct ElementBounds
    {
        ElementBounds(double l_u, double u_u,
                double l_v, double u_v,
                double l_w, double u_w)
        :
            lower_u(l_u),
            upper_u(u_u),
            lower_v(l_v),
            upper_v(u_v),
            lower_w(l_w),
            upper_w(u_w)
        {}
        
        double lower_u;
        double upper_u;
        double lower_v;
        double upper_v;
        double lower_w;
        double upper_w;
        
    };
    
    /// An iterator class that iterates over a set of evenly spaced points over
    /// an element (knot span)
    class IParentSample
    {
    public:
        
        explicit IParentSample(const ElementBounds& e,
                               const uint npts)
        :
            mCurrentEl(e),
            mKnotLinePts(npts),
            mCurrentIndex(0)
        {
            mSamplePtN = npts * npts * npts;
            mUIncrement = ( mCurrentEl.upper_u - mCurrentEl.lower_u ) / static_cast< double >( npts - 1 );
            mVIncrement = ( mCurrentEl.upper_v - mCurrentEl.lower_v ) / static_cast< double >( npts - 1 );
            mWIncrement = ( mCurrentEl.upper_w - mCurrentEl.lower_w ) / static_cast< double >( npts - 1 );
        }
        
        /// Default constuctors
        IParentSample(const uint npts,
                      const double lu = -1.0,
                      const double uu = 1.0,
                      const double lv = -1.0,
                      const double uv = 1.0,
                      const double lw = -1.0,
                      const double uw = 1.0)
        :
            IParentSample(ElementBounds(lu,uu,lv,uv,lw,uw), npts)
        {}
        
        bool isDone() const
        {
            return mCurrentIndex == mSamplePtN;
        }
        
        /// prefix increment operator
        IParentSample& operator++()
        {
            ++mCurrentIndex;
            return *this;
        }
        
        /// postfix increment operator
        IParentSample operator++( int i )
        {
            IParentSample temp( *this );
            ++(*this);
            return temp;
        }
        
        ParamCoord getCurrentPt() const
        {
            const uint nuv = pointN() * pointN();
            
            const uint k = currentIndex() / nuv;
            const uint i = currentIndex() % pointN();
            const uint j = (currentIndex() - k * nuv) / pointN();

            /// TODO: this needs testing 
            return ParamCoord(element().lower_u + sampleIncrement(U) * i,
                              element().lower_v + sampleIncrement(V) * j,
                              element().lower_w + sampleIncrement(W) * k);
        }
        
        /// index getter
        uint currentIndex() const
        {return mCurrentIndex; }
        
    protected:
        
    private:
        
        /// number of points getter (same for all directions)
        uint pointN() const
        { return mKnotLinePts; }
        
        /// element getter
        const ElementBounds& element() const { return mCurrentEl; }
        
        /// sample point incrmeent getter
        double sampleIncrement(ParamDir dir) const
        {
            switch(dir)
            {
                case ParamDir::U:
                    return mUIncrement;

                case ParamDir::V:
                    return mVIncrement;

                case ParamDir::W:
                    return mWIncrement;
                default:
                    error("Bad parametric direction");
            }
        }
        
        /// The current element we are iterating over
        const ElementBounds mCurrentEl;
        
        /// The number of points along each knot line
        const uint mKnotLinePts;
        
        /// the total number of sample points for the element
        uint mSamplePtN;
        
        /// the current point index
        uint mCurrentIndex;
        
        /// The increment between points in the u direction
        double mUIncrement;
        
        /// The increment between points in the v direction
        double mVIncrement;
        
        /// The increment between points in the w direction
        double mWIncrement;
        
        
    };
}
#endif