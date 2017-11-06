#ifndef TRINURBS_ANALYSIS_ELEMENT_H
#define TRINURBS_ANALYSIS_ELEMENT_H

#include <memory>
#include <cassert>
#include <algorithm>
#include <map>
#include <stdexcept>

#include <boost/math/special_functions.hpp>

#include "GeometryElement.h"
#include "base.h"

namespace trinurbs
{
    
    /// An abstract class that specifies the interface for an element
    /// used for analysis. The two main types of elements that will
    /// be subclasses from this elements with nodal and vector basis functions.
    
    /// The element is parameterised by T which specifies the format
    /// that a basis function is returned.
    ///
    /// E.g. Nodal basis functions T = double
    /// Vector basis functions T = DoubleVec
    
    template<typename T>
    class AnalysisElement : public GeometryElement {
        
    public:
        
        /// Get the basis function return type as specified by
        /// the template parameter.
        typedef T basisReturnType;
        
        /// A copy function that must be implemented by dervied classes
        virtual std::unique_ptr<AnalysisElement> copy() const = 0;
        
        /// Number of components. (= 1 for nodal basis, >1 for vector basis)
        virtual uint componentN() const = 0;
        
        /// Number of non-zero basis functions
        virtual uint basisFuncN() const = 0;
        
        /// Get basis function indices local to this parameter space
        virtual UIntVec localBasisIVec() const = 0;
        
        /// Get the global basis function indices. These are global
        /// in the sense they relate to dofs in the forest.
        virtual UIntVec globalBasisIVec() const = 0;
        
        /// Get the global basis function indices for the given face.
        virtual UIntVec globalBasisIVec(const Face f) const = 0;
        
        /// Get the global basis function indices for the given edge.
        virtual UIntVec globalBasisIVec(const Edge e) const = 0;
        
        /// Get the signed global basis function index vector. Any degenerate
        /// indices will be signified by -1
        virtual IntVec signedGlobalBasisIVec() const
        {
            throw std::runtime_error("signed basis function index generator not implemented yet");
        }
        
        /// Return non-zero basis function set for the given parent
        /// coordinate
        virtual std::vector<T> basis(const double xi,
                                     const double eta,
                                     const double zeta) const = 0;
        
        /// Evauate basis with given tangent vectors
        virtual std::vector<T> basis(const double xi,
                                     const double eta,
                                     const double zeta,
                                     const Point3D& t1,
                                     const Point3D& t2,
                                     const Point3D& t3) const
        {
            throw std::runtime_error("Basis function not implemented.");
        }
        
        
        
        /// Return basis without Piola tranform.
        /// Default behaviour is to return the above basis() function
        virtual std::vector<T> localBasis(const double xi,
                                          const double eta,
                                          const double zeta) const
        { return basis(xi,eta,zeta); }
        
        /// Get the basis derivatives
        virtual std::vector<T> basisDers(const double xi,
                                         const double eta,
                                         const double zeta,
                                         const DerivType dtype) const = 0;
        
        /// Get the degree for given parametric direction and component
        virtual uint degree(const ParamDir dir, const uint comp = 0) const = 0;
        
        /// Get the degree vector for a specified component
        virtual UIntVec degree(const uint comp = 0) const = 0;
        
        /// Get the integration order for a specified component with an
        /// offset applied to the usual p + 1 rule.
        UIntVec integrationOrder(const int offset = 0,
                                 const uint comp = 0) const
        {
            assert(comp < componentN());
            UIntVec d = degree(comp);
            std::for_each(d.begin(), d.end(), [&offset](uint& v){ v += 1 + offset; });
            return d;
        }
        
        /// Return the integration order (equal in both parameteric directions)
        /// and equal to the maximum degree + 1
        UIntVec equalIntegrationOrder(const int offset = 0) const
        {
            uint max = 0;
            for(uint icomp = 0; icomp < componentN(); ++icomp) {
                const auto int_order = integrationOrder(offset, icomp);
                const auto d_max = std::max_element(int_order.begin(), int_order.end());
                max = (*d_max > max) ? *d_max : max;
            }
            return UIntVec(componentN(), max);
        }
        
        /// Add a vertex connected element
//        void setVertexConnectedEl(const Vertex v, const AnalysisElement* el)
//        { mVertexEls[v] = el; }
        
        /// Add an edge connected element
//        void setEdgeConnectedEl(const Edge e, const AnalysisElement* el)
//        { mEdgeEls[e] = el; }
        
        /// Get vertex connected element at specified vertex
//        const AnalysisElement* getVertexConnectedEl(const Vertex v) const
//        {
//            auto it = mVertexEls.find(v);
//            if(it != mVertexEls.end())
//                return it->second;
//            else
//                return nullptr;
//        }
        
        /// Get edge connected element at specified vertex
//        const AnalysisElement* getEdgeConnectedEl(const Edge e) const
//        {
//            auto it = mEdgeEls.find(e);
//            if(it != mEdgeEls.end())
//                return it->second;
//            else
//                return nullptr;
//        }
        
        /// Print function as used by output operator
        virtual void print(std::ostream& ost) const = 0;
        
        /// parent element getter
        const AnalysisElement<double>* parent() const { return mpParent; }
        
        /// parent setter
        void setParent(const AnalysisElement<double>* p) { mpParent = p; }
        
        /// Given a parent coordinate for this element in [-1,1]
        /// transform this to the parent coordinate in the parent element
        /// also in [-1,1]
        GPt3D transformToParentElParentCoord(const GPt3D& gp) const
        {
            ParamCoord c = getParamCoord(gp);
            return parent()->getParentCoord(c.u, c.v, c.w);
        }
        
        /// Does the given face on this analysis element lie on the parent
        /// (geometry) face?
        bool liesOnParentElFace(const Face f) const
        {
            assert(parent());
            
            switch(f) {
                case Face::FACE0:
                    return logically_equal(lowerBound(ParamDir::V), parent()->lowerBound(ParamDir::V));
                case Face::FACE1:
                    return logically_equal(upperBound(ParamDir::V), parent()->upperBound(ParamDir::V));
                case Face::FACE2:
                    return logically_equal(lowerBound(ParamDir::U), parent()->lowerBound(ParamDir::U));
                case Face::FACE3:
                    return logically_equal(upperBound(ParamDir::U), parent()->upperBound(ParamDir::U));
                case Face::FACE4:
                    return logically_equal(lowerBound(ParamDir::W), parent()->lowerBound(ParamDir::W));
                case Face::FACE5:
                    return logically_equal(upperBound(ParamDir::W), parent()->upperBound(ParamDir::W));
            }
        }
        
        /// Does this element lie on the given parent edge?
        bool liesOnParentElEdge(const Edge e) const
        {
            switch(e)
            {
                case Edge::EDGE0:
                    return liesOnParentElFace(Face::FACE0) && liesOnParentElFace(Face::FACE4);
                case Edge::EDGE1:
                    return liesOnParentElFace(Face::FACE0) && liesOnParentElFace(Face::FACE5);
                case Edge::EDGE2:
                    return liesOnParentElFace(Face::FACE0) && liesOnParentElFace(Face::FACE2);
                case Edge::EDGE3:
                    return liesOnParentElFace(Face::FACE0) && liesOnParentElFace(Face::FACE3);
                case Edge::EDGE4:
                    return liesOnParentElFace(Face::FACE1) && liesOnParentElFace(Face::FACE4);
                case Edge::EDGE5:
                    return liesOnParentElFace(Face::FACE1) && liesOnParentElFace(Face::FACE5);
                case Edge::EDGE6:
                    return liesOnParentElFace(Face::FACE1) && liesOnParentElFace(Face::FACE2);
                case Edge::EDGE7:
                    return liesOnParentElFace(Face::FACE1) && liesOnParentElFace(Face::FACE3);
                case Edge::EDGE8:
                    return liesOnParentElFace(Face::FACE2) && liesOnParentElFace(Face::FACE4);
                case Edge::EDGE9:
                    return liesOnParentElFace(Face::FACE3) && liesOnParentElFace(Face::FACE4);
                case Edge::EDGE10:
                    return liesOnParentElFace(Face::FACE2) && liesOnParentElFace(Face::FACE5);
                case Edge::EDGE11:
                    return liesOnParentElFace(Face::FACE3) && liesOnParentElFace(Face::FACE5);
            }
        }
        
        /// Does this element lie on the given parent vertex?
        bool liesOnParentElVertex(const Vertex v) const
        {
            switch(v)
            {
                case Vertex::VERTEX0:
                    return liesOnParentElFace(Face::FACE0) && liesOnParentElFace(Face::FACE2) && liesOnParentElFace(Face::FACE4);
                case Vertex::VERTEX1:
                    return liesOnParentElFace(Face::FACE0) && liesOnParentElFace(Face::FACE3) && liesOnParentElFace(Face::FACE4);
                case Vertex::VERTEX2:
                    return liesOnParentElFace(Face::FACE0) && liesOnParentElFace(Face::FACE2) && liesOnParentElFace(Face::FACE5);
                case Vertex::VERTEX3:
                    return liesOnParentElFace(Face::FACE0) && liesOnParentElFace(Face::FACE3) && liesOnParentElFace(Face::FACE5);
                case Vertex::VERTEX4:
                    return liesOnParentElFace(Face::FACE1) && liesOnParentElFace(Face::FACE2) && liesOnParentElFace(Face::FACE4);
                case Vertex::VERTEX5:
                    return liesOnParentElFace(Face::FACE1) && liesOnParentElFace(Face::FACE3) && liesOnParentElFace(Face::FACE4);
                case Vertex::VERTEX6:
                    return liesOnParentElFace(Face::FACE1) && liesOnParentElFace(Face::FACE2) && liesOnParentElFace(Face::FACE5);
                case Vertex::VERTEX7:
                    return liesOnParentElFace(Face::FACE1) && liesOnParentElFace(Face::FACE3) && liesOnParentElFace(Face::FACE5);
                    
            }
        }
        
        
    protected:
        
        /// Constructor which simply calls BaseElement constructor
        AnalysisElement(const Geometry& g,
                        const uint ispace,
                        const DoublePairVec& knots)
        :
            GeometryElement(g, ispace, knots),
            mpParent(nullptr)
        {}
        
    private:
        
        /// Reference to elements connected to vertices of this element.
        /// Assumption that no more than one element is connected to each vertex
//        std::map<Vertex, const AnalysisElement*> mVertexEls;
        
        /// Reference to elements connected to edges of this element.
//        std::map<Edge, const AnalysisElement*> mEdgeEls;
        
        /// Reference to geometry space this element lives in
        //const BSplineSpace* mpGeomSpace;
        
        /// Pointer to parent element that lives in the primal forest
        const AnalysisElement<double>* mpParent;
        
        /// Overload output operator
        friend std::ostream& operator<<(std::ostream& ost, const AnalysisElement& e)
        { e.print(ost); return ost; }
    };
    
    /// Are the given elements connected at an edge?
    /// This implementation does not treat a degenerate shared edge
    /// as a valid case
//    template<typename T>
//    bool edgeConnected(const AnalysisElement<T>& e1,
//                       const AnalysisElement<T>& e2,
//                       Edge& edg1,
//                       Edge& edg2)
//    {
//        
//        for(uint edge1 = 0; edge1 < NEDGES; ++edge1) {
//            const auto adj_e1 = e1.getEdgeConnectedEl(edgeType(edge1));
//            if(nullptr == adj_e1) continue;
//            if(adj_e1 == &e2) { // The elements are edge connected
//                for(uint edge2 = 0; edge2 < NEDGES; ++edge2) {
//                    const auto adj_e2 = e2.getEdgeConnectedEl(edgeType(edge2));
//                    if(nullptr == adj_e2) return false; // is this correct?
//                    if(adj_e2 == &e1) {
//                        edg1 = edgeType(edge1); edg2 = edgeType(edge2);
//                        return true;
//                    }
//                }
//                throw std::runtime_error("Bad edge connectivity for edge connected elements.");
//            }
//        }
//        return false;
//    }
    
    /// Are the given elements connected at a vertex
//    template<typename T>
//    bool vertexConnected(const AnalysisElement<T>& e1,
//                         const AnalysisElement<T>& e2,
//                         Vertex& vtx1,
//                         Vertex& vtx2)
//    {
//        for(uint v1 = 0; v1 < NVERTICES; ++v1) {
//            const auto adj_e1 = e1.getVertexConnectedEl(vertexType(v1));
//            if(nullptr == adj_e1) continue;
//            if(adj_e1 == &e2) { // the elements are vertex connected
//                for(uint v2 = 0; v2 < NVERTICES; ++v2) {
//                    const auto adj_e2 = e2.getVertexConnectedEl(vertexType(v2));
//                    if(nullptr == adj_e2) continue;
//                    if(adj_e2 == &e1) {
//                        vtx1 = vertexType(v1); vtx2 = vertexType(v2);
//                        return true;
//                    }
//                }
//                throw std::runtime_error("Bad vertex connectivity for vertex connected elements.");
//            }
//        }
//        return false;
//    }
//    
//    template<typename T>
//    bool connectedAtDegeneratePt(const AnalysisElement<T>& e1,
//                                 const AnalysisElement<T>& e2)
//    {
//        const auto e1_pair = e1.degenerateEdge();
//        const auto e2_pair = e2.degenerateEdge();
//        
//        if(!e1_pair.first || !e2_pair.first)
//            return false;
//        const double tol = 1.e-9;
//        
//        Point3D p1;
//        switch(e1_pair.second)
//        {
//            case Edge::EDGE0:
//                p1 = e1.evalVertex(Vertex::VERTEX0);
//                break;
//            case Edge::EDGE1:
//                p1 = e1.evalVertex(Vertex::VERTEX3);
//                break;
//            case Edge::EDGE2:
//                p1 = e1.evalVertex(Vertex::VERTEX2);
//                break;
//            case Edge::EDGE3:
//                p1 = e1.evalVertex(Vertex::VERTEX1);
//                break;
//        }
//        
//        Point3D p2;
//        switch(e2_pair.second)
//        {
//            case Edge::EDGE0:
//                p2 = e2.evalVertex(Vertex::VERTEX0);
//                break;
//            case Edge::EDGE1:
//                p2 = e2.evalVertex(Vertex::VERTEX3);
//                break;
//            case Edge::EDGE2:
//                p2 = e2.evalVertex(Vertex::VERTEX2);
//                break;
//            case Edge::EDGE3:
//                p2 = e2.evalVertex(Vertex::VERTEX1);
//                break;
//        }
//        if(nurbs::dist(p1, p2) < tol)
//            return true;
//        else
//            return false;
//    }
    
    /// Typedefs for nodal and vector analysis elements
    typedef AnalysisElement<double> NAnalysisElement;
    typedef AnalysisElement<DoubleVec> VAnalysisElement;
    
}

#endif


