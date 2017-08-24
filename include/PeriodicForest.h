#ifndef TRINURBS_PERIODIC_FOREST_H
#define TRINURBS_PERIODIC_FOREST_H

#include "Forest.h"

namespace trinurbs {
    
    /// A represnetation of a periodic forest used to represent a unit cell
    /// for multiscale geometry.
    ///
    /// There are several functions will are useful for generating the connectivity
    /// of a multiscale geometry such as determining which elements are
    /// 'opposite' to one another, node indices along edges, faces, vertices etc.
    ///
    class PeriodicForest : public Forest
    {
    public:
        
        /// Default constructor
        PeriodicForest()
        :
            Forest()
        {}
        
        /// Constructor with geometry refernece
        PeriodicForest(const Geometry& g)
        :
            Forest(g)
        {}
        
        /// Use default destructor
        virtual ~PeriodicForest() = default;
        
        /// Copy constructor
        PeriodicForest(const PeriodicForest& f)
        :
            Forest(f)
        {
            // TODO: specifiy how any member data is copied
        }
        
        /// Assignment operator
        PeriodicForest& operator=(const PeriodicForest& f)
        {
            PeriodicForest temp(f);
            *this = std::move(temp);
            return *this;
        }
        
        
        /// Move constructor. Simply move all the data across
        PeriodicForest(PeriodicForest&& f) = default;
        
        /// Move assignment operator. As before, move all the data across.
        PeriodicForest& operator=(PeriodicForest&& f) = default;
        
    protected:
        
        
    private:
        
        
        /// Initialise data structures (after construction or refinement)
        void init();
        
        
    };
    
    
}
#endif