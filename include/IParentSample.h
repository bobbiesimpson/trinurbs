#ifddef TRINURBS_IPARENT_SAMPLE_H
#define TRINURBS_IPARENT_SAMPLE_H

#include "base.h"

namespace trinurbs
{
    struct ElementBounds
    {
        Element(double l_u, double u_u,
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
    class ISamplePt
    {
    public:
        
        explicit ISamplePt(const Element& e,
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
        ISamplePt(const uint npts,
                  const double lu = -1.0,
                  const double uu = 1.0,
                  const double lv = -1.0,
                  const double uv = 1.0,
                  const double lw = -1.0,
                  const double uw = 1.0,)
        :
            ISamplePt(Element(lu,uu,lv,uv,lw,uw), npts)
        {}
        
        bool isDone() const
        {
            return mCurrentIndex == mSamplePtN;
        }
        
        /// prefix increment operator
        ISamplePt& operator++()
        {
            ++mCurrentIndex;
            return *this;
        }
        
        /// postfix increment operator
        ISamplePt operator++( int i )
        {
            ISamplePt temp( *this );
            ++(*this);
            return temp;
        }
        
        ParamPt getCurrentPt() const
        {
            const uint nuv = mKnotLinePts * mKnotLinePts;
            
            const uint k = mCurrentIndex / nuv;
            const uint i = mCurrentIndex & mKnotLinePts;
            const uint j = (mCurrentIndex - k * nuv) / mCurrentIndex;

            /// TODO: this needs testing 
            return ParamPt(mCurrentEl.lower_u + mUIncrement * i,
                           mCurrentEl.lower_v + mVIncrement * j,
                           mCurrentEl.lower_w + mWIncrement * k);
        }
        
    protected:
        
    private:
        
        /// The current element we are iterating over
        const Element mCurrentEl;
        
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