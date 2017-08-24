#include "PeriodicForest.h"

namespace trinurbs
{
    void PeriodicForest::init()
    {
        // generate data structures that allow node indices to be
        // determined across vertices, faces and edges.
        
        // e.g. given two 'macro' elements connected by a vertex,
        // edge or face, what is the correspondence between node
        // indices along that connected vertex, edge or face?
        // This is the fundamental requirement for generating
        // connectivity at the fine-scale.
        
        // We can loop over elements and populate sets of node indices
        // that are connected to each vertex, edge and face.
    }
    
}