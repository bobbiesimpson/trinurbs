#ifndef TRINURBS_BASE_H
#define TRINURBS_BASE_H

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <limits>
#include <iterator>
#include <memory>
#include <complex>
#include <tuple>

#include <boost/icl/continuous_interval.hpp>

namespace trinurbs
{
    
    /// error function (wraps a runtime exception)
    void error( const std::string& msg );
    
    /// Useful typedefs
    typedef unsigned int uint;
    typedef std::vector< double > DoubleVec;
    typedef std::vector< std::vector< double > > DoubleVecVec;
    typedef std::vector<int> IntVec;
    typedef std::vector< uint > UIntVec;
    typedef std::vector<UIntVec> UIntVecVec;
    typedef std::complex< double > ComplexDouble;
    typedef std::vector<std::complex<double>> ComplexDoubleVec;
    typedef std::pair<double, double> DoublePair;
    typedef std::vector<std::pair<double, double>> DoublePairVec;
    
    typedef boost::icl::continuous_interval<double> Interval;
    
    /// Tolerance that we use throughout the code
    const double TOL = std::numeric_limits< double >::epsilon();
    
    /// Default increment for NURBS curve interpolation
    const double INC = 0.1;
    
    /// Default number of grid points for vtk output
    const uint DEFAULT_NGRID_PTS = 12;
    
    /// Default number of guass points
    const uint DEFAULT_NGPS = 4;
    
    /// Default number of elements per quadtree cell
    const uint DEFAULT_QT_ELN = 1;
    
    /// Default quadtree tolerance
    const double DEFAULT_QT_TOL = 1.0e-3;
    
    /// pi = 3.141...
    const double PI = atan( 1.0 ) * 4.0;
    
    /// Specifies the order of a B-spline derivative
    enum DerivOrder
    {
        D1 = 1,
        D2
    };
    
    /// Enumeration of derivative type
    enum DerivType {
        DU = 1,
        DV = 2,
        DW = 3
    };
    
    /// Specifies a parametric direction
    enum ParamDir
    {
        U = 0,
        V,
        W
        
    };
    
    /// Cast uint to ParamDir type
    ParamDir ParamDirType(const uint d);
    
    /// A simple struct for holding a gauss point
    struct GPt3D
    {
        GPt3D(const double xi_in = 0.0,
              const double eta_in = 0.0,
              const double zeta_in = 0.0)
        :
        xi(xi_in),
        eta(eta_in),
        zeta(zeta_in) {}
        
        /// Access component values
        double get(const uint i) const
        {
            if(i == 0)
                return xi;
            else if(i == 1)
                return eta;
            else if(i == 2)
                return zeta;
            else {
                error("Bad index specified for GaussPt struct");
                return 0.0;
            }
        }
        
        /// Components
        double xi;
        double eta;
        double zeta;
        
    };
    
    /// Overload comparison operator
//    bool operator<(const GPt3D& g1, const GPt3D& g2);
    
    /// Overload output operator
    std::ostream& operator<<(std::ostream& ost, const GPt3D& gpt);
    
    /// Representing each of the carteisan components
    enum CartesianComponent
    {
        X = 0,
        Y = 1,
        Z = 2
    };
    
    /// Boundary condition type enumeration
    enum BCType
    {
        DIRICHLET = 0,
        NEUMANN
    };
    
    /// A simple struct for representing a parametric coordinate
    typedef struct {
        double u;
        double v;
        double w;
    } ParamCoord;
    
    /// Vertex enumeration
    ///    6 ------7
    ///   /|      /|
    ///  / |     / |
    /// 2 ----- 3  |
    /// |  4 ---|- 5
    /// | /     | /
    /// |/      |/
    /// 0 ------1
    
    /// Positive u: 0 -> 1
    /// Positive v: 0 -> 4
    /// Positive w: 0 -> 2
    
    enum class Vertex {
        VERTEX0 = 0,
        VERTEX1 = 1,
        VERTEX2 = 2,
        VERTEX3 = 3,
        VERTEX4 = 4,
        VERTEX5 = 5,
        VERTEX6 = 6,
        VERTEX7 = 7
    };
    
    /// Global const for # vertices in each cell
    const std::size_t NVERTICES = 8;
    
    /// Overload vertex output operator
    std::ostream& operator<<(std::ostream& ost, Vertex v);
    
    /// Get vertex type given an unsigned integer
    Vertex vertexType(const uint v);
    
    /// Edge enumeration
    ///      ---5----
    ///    /|       /|
    ///  10 |      11|
    ///  /  6     /  7
    ///  ---|-1---   |
    /// |   |--4-|---|
    /// |  /     |   /
    /// 2 8      3  9
    /// |/       | /
    ///  ---0----
    
    /// Parametric coord system same as defined for Vertex numbering above

    enum class Edge {
        EDGE0 = 0,
        EDGE1 = 1,
        EDGE2 = 2,
        EDGE3 = 3,
        EDGE4 = 4,
        EDGE5 = 5,
        EDGE6 = 6,
        EDGE7 = 7,
        EDGE8 = 8,
        EDGE9 = 9,
        EDGE10 = 10,
        EDGE11 = 11
    };
    
    /// Global const for # edges in each cell
    const std::size_t NEDGES = 12;
    
    /// FACE enumeration
    ///     -------
    ///   /|  5   /|
    ///  / |   1 / |
    ///  -------   |
    /// |2  ----|3-|
    /// | / 0   | /
    /// |/   4  |/
    ///  -------
    
    /// Face coord system:
    
    /// Given a face index, we define a face coord system (u,v) where the
    /// origin is defined by the vertex with the lowest index and the
    /// u-axis is defined from this vertex to the vertex with
    /// the next lowest index.  This coord system is used to define
    /// the ordering of local basis function indices when returned
    /// from functions like localBasisIVec(Face f, BSplineSpace& s)
    
    /// Parametric coord system same as defined for Vertex numbering above
    enum class Face {
        FACE0 = 0,
        FACE1,
        FACE2,
        FACE3,
        FACE4,
        FACE5
    };
    
    /// Global const for # edges in each cell
    const std::size_t NFACES = 6;
    
    /// Cast uint into face enum
    Face faceType(const uint f);
    
    /// Overload output operator
    std::ostream& operator<<(std::ostream& ost, Edge e);
    
    /// Get edge type given an unsigned integer
    Edge edgeType(const uint e);
    
    /// Sign enumeration. Used for edges and vector basis
    enum class Sign {
        POSITIVE,
        NEGATIVE
    };
    
    /// Cast Sign to double
    double asDouble(Sign s);
    
    /// Overload output operator for Sign
    std::ostream& operator<<(std::ostream& ost, Sign s);
    
    /// Vector basis direction
    enum class ContinuityType {
        NORMAL,
        TANGENT
    };
    
    /// Enum for rotation angles in increments of 90 degrees
    /// Counter clockwise positive
    enum class RotationAngle {
        NINETYDEGREES,
        ONEEIGHTYDEGREES,
        TWOSEVENTYDEGREES
    };
    
    /// overload output operator for BCType
    std::ostream& operator<<( std::ostream& ost, const BCType& bc );
    
    /// helper functions
    template< class T >
    bool approximatelyEqual( T a, T b, T epsilon )
    {
        return fabs( a - b ) <= ( ( fabs( a ) < fabs( b ) ? fabs( b ) : fabs( a ) ) * epsilon );
    }
    
    template< class T >
    bool essentiallyEqual( T a, T b, T epsilon )
    {
        return fabs( a - b ) <= ( ( fabs( a ) > fabs( b ) ? fabs( b ) : fabs( a ) ) * epsilon );
    }
    
    template< class T >
    bool lessThanOrEqual( T a, T b, T epsilon )
    {
        return ( a < b ) || essentiallyEqual( a, b, epsilon );
    }
    
    template< class T >
    bool greaterThanOrEqual( T a, T b, T epsilon )
    {
        return ( a > b ) || essentiallyEqual( a, b, epsilon );
    }
    
    /// print a vector
    template< typename T >
    void printVector( const std::vector< T >& vec, std::ostream& ost )
    {
        ost << "{";
        if(vec.size() != 0) {
            std::copy( vec.begin(), vec.end() - 1, std::ostream_iterator< T >( ost, "," ) );
            ost << vec.back();
        }
        ost <<  "}";
    }
    
    /// overload output operator for a generic vector
    template<typename T>
    std::ostream& operator<<(std::ostream& ost, const std::vector<T>& v)
    {
        ost << "{";
        for(typename std::vector<T>::const_iterator it = v.begin(); it != v.end(); ++it) {
            ost << *it;
            if(it != v.end() - 1)
                ost << ", ";
        }
        ost << "}";
        return ost;
    }
    
    template< class T >
    void printVal( const std::string& n, T val, std::ostream& ost )
    {
        ost << n << " = " << val << "\n";
    }
    
    /// the factory function for creating unique_ptr instances (due to be implemented in c++14)
    /// Use as e.g. std::unique_ptr< Widget > pw = make_unique< Widget >();
    template< typename T, typename ...Args >
    std::unique_ptr< T > make_unique( Args&& ...args )
    {
        return std::unique_ptr< T >( new T( std::forward< Args >(args)... ) );
    }
    
    /// check the input stream if we have reached the specified end char. Set the stream
    /// to good if equal, otherwise throw an exception
    void endOfLoop( std::istream& ist, const char delim, const std::string& msg );
    
    /// Does the given parent coordinate lie in the standard [-1,1] x [-1,1]
    /// interval?
    bool validParentCoord(const double u, const double v);
    
    /// N Evenly spaced values defined in the interval [a,b]
    DoubleVec range(const double a, const double b, const uint N);
    
    /// Return a vector that splits 'mem' into equal parts.
    /// If a remainder exists, the last term is modified accordingly.
    std::vector<long int> bounds(long int parts, long int mem);
    
    /// Transpose a matrix (vector of vectors)
    template <typename T>
    void transpose(std::vector<std::vector<T>>& m)
    {
        auto copy = m;
        const size_t nrows = m.size();
        const size_t ncols = m.front().size();
        
        m.clear();
        
        m.resize(ncols);
        for(size_t i = 0; i < ncols; ++i)
            m[i].resize(nrows);
        
        for(size_t i = 0; i < ncols; ++i)
            for(size_t j = 0; j < nrows; ++j)
                m[i][j] = copy[j][i];
        
    }
    
    template<typename T>
    void reverseColumns(std::vector<std::vector<T>>& m)
    {
        auto copy = m;
        const size_t nrows = m.size();
        const size_t ncols = m.front().size();
        
        for(size_t j = 0; j < ncols; ++j)
            for(size_t i = 0; i < nrows; ++i)
                m[i][j] = copy[nrows-i-1][j];
    }
    
    template<typename T>
    void reverseRows(std::vector<std::vector<T>>& m)
    {
        auto copy = m;
        const size_t nrows = m.size();
        const size_t ncols = m.front().size();
        
        for(size_t i = 0; i < nrows; ++i)
            for(size_t j = 0; j < ncols; ++j)
                m[i][j] = copy[i][ncols-j-1];
    }
    
    /// Given a matrix (vector of vectors), rotation angle
    /// and (optional) transpose flag, modify the matrix entries
    /// accordingly.
    template<typename T>
    void rotateMatrixEntries(std::vector<std::vector<T>>& m,
                             const RotationAngle rangle,
                             const bool transposeapplied = false)
    {
        // now apply each of the three angle cases
        switch(rangle)
        {
            case RotationAngle::NINETYDEGREES:
                transpose(m);
                reverseColumns(m);
                break;
            case RotationAngle::ONEEIGHTYDEGREES:
                reverseRows(m);
                reverseColumns(m);
                break;
            case RotationAngle::TWOSEVENTYDEGREES:
                transpose(m);
                reverseRows(m);
                break;
        }
        if(transposeapplied)
            transpose(m);
    }
}

// code for allowing for tuple types as keys in unordered_maps

namespace std{
    namespace
    {
        
        // Code from boost
        // Reciprocal of the golden ratio helps spread entropy
        //     and handles duplicates.
        // See Mike Seymour in magic-numbers-in-boosthash-combine:
        //     http://stackoverflow.com/questions/4948780
        
        template <class T>
        inline void hash_combine(std::size_t& seed, T const& v)
        {
            seed ^= hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        
        // Recursive template code derived from Matthieu M.
        template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
        struct HashValueImpl
        {
            static void apply(size_t& seed, Tuple const& tuple)
            {
                HashValueImpl<Tuple, Index-1>::apply(seed, tuple);
                hash_combine(seed, get<Index>(tuple));
            }
        };
        
        template <class Tuple>
        struct HashValueImpl<Tuple,0>
        {
            static void apply(size_t& seed, Tuple const& tuple)
            {
                hash_combine(seed, get<0>(tuple));
            }
        };
    }
    
    template <typename ... TT>
    struct hash<std::tuple<TT...>>
    {
        size_t
        operator()(std::tuple<TT...> const& tt) const
        {
            size_t seed = 0;
            HashValueImpl<std::tuple<TT...> >::apply(seed, tt);
            return seed;
        }
        
    };
}

#endif
