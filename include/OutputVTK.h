#ifndef TRINURBS_OUTPUTVTK_H
#define TRINURBS_OUTPUTVTK_H

#include <string>
#include <vector>
#include <mutex>

#include "base.h"
#include "IParentSample.h"
#include "Point3D.h"

namespace trinurbs
{
    /// Forward declarations
    class Forest;
    class MultiscaleForest;
    
    /// Class responsible for output to VTK format
    class OutputVTK {
        
    public:
        
        /// Default constructor
        OutputVTK() : OutputVTK("unnamed_output") {}
        
        /// Construct with filename
        OutputVTK(const std::string& f,
                  const uint nsample = DEFAULT_NGRID_PTS)
        :
            mFilename(f),
            mSamplePtN(nsample) {}
        
        /// output the geometry of a forest to VTK
        void outputForestGeometry(const Forest& f) const;
        
        /// output geometry of multiscale forest
        void outputMultiscaleForestGeometry(const MultiscaleForest& f) const;
        
        /// Write a complex nodal field to a vtu file
        /// if nlocaldof > 1 then we assume the soln is ordered
        /// as e.g. [u^1_x u^1_y u^1_z u^2_x u^2_y ...]
        void outputNodalField(const Forest& f,
                              const std::string& fieldname,
                              const std::vector<double>& soln,
                              const uint nlocaldof = 1) const;
        
        /// Sample point number setter
        void setSamplePtN(const uint n)
        {
            if(n < 1)
                error("Must specifiy non-zero sample points for output");
            mSamplePtN = n;
        }
        
        /// Filename getter
        const std::string& filename() const
        {
            return mFilename;
        }
        
        /// Filename setter
        void setFilename(const std::string& newfile)
        {
            mFilename = newfile;
        }
        
        /// Sample point number getter
        uint samplePtN() const { return mSamplePtN; }
        
    private:
        
        /// Filename we write to
        std::string mFilename;
        
        /// Number of sample points
        uint mSamplePtN;
        
        /// Get the mutex
        std::mutex& mutex() const
        {
            return mMutex;
        }
        
        /// Mutex for this class
        mutable std::mutex mMutex;
    };
    
}

#endif