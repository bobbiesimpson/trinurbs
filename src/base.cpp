#include "base.h"
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <cmath>

#include <boost/math/special_functions.hpp>

namespace trinurbs
{
    void error( const std::string& msg )
    {
        throw std::runtime_error( "Error: " + msg );
    }
//    
//    bool operator<(const GPt3D& g1, const GPt3D& g2)
//    {
//        if(g1.xi != g2.xi)
//            return g1.s < g2.s;
//        return g1.t < g2.t;
//    }
    
    std::ostream& operator<<(std::ostream& ost, const GPt3D& gpt)
    {
        ost << gpt.xi << "\t" << gpt.eta << "\t" << gpt.zeta;
        return ost;
    }
    
    ParamDir ParamDirType(const uint d)
    {
        assert(d < 3);
        if(d == 0)
            return U;
        else if(d == 1)
            return V;
        else
            return W;
    }
    
    Vertex vertexType(const uint v)
    {
        switch(v)
        {
            case 0: return Vertex::VERTEX0;
            case 1: return Vertex::VERTEX1;
            case 2: return Vertex::VERTEX2;
            case 3: return Vertex::VERTEX3;
            case 4: return Vertex::VERTEX4;
            case 5: return Vertex::VERTEX5;
            case 6: return Vertex::VERTEX6;
            case 7: return Vertex::VERTEX7;
            default: throw std::runtime_error("Bad vertex specified");
        }
    }
    
    Face faceType(const uint f)
    {
        switch(f)
        {
            case 0: return Face::FACE0;
            case 1: return Face::FACE1;
            case 2: return Face::FACE2;
            case 3: return Face::FACE3;
            case 4: return Face::FACE4;
            case 5: return Face::FACE5;
            default: throw std::runtime_error("Bad face index specified");
                
        }
    }
    
    Edge edgeType(const uint e)
    {
        switch (e) {
            case 0: return Edge::EDGE0;
            case 1: return Edge::EDGE1;
            case 2: return Edge::EDGE2;
            case 3: return Edge::EDGE3;
            case 4: return Edge::EDGE4;
            case 5: return Edge::EDGE5;
            case 6: return Edge::EDGE6;
            case 7: return Edge::EDGE7;
            case 8: return Edge::EDGE8;
            case 9: return Edge::EDGE9;
            case 10: return Edge::EDGE10;
            case 11: return Edge::EDGE11;
            default: throw std::runtime_error("Bad edge specified");
        }
    }
    
    double asDouble(Sign s)
    {
        switch(s)
        {
            case Sign::POSITIVE: return 1.0;
            case Sign::NEGATIVE: return -1.0;
            default: throw std::runtime_error("Bad sign specified");
        }
    }
    
    std::ostream& operator<<(std::ostream& ost, Vertex v)
    {
        switch (v) {
            case Vertex::VERTEX0:
                ost << "VERTEX0";
                break;
            case Vertex::VERTEX1:
                ost << "VERTEX1";
                break;
            case Vertex::VERTEX2:
                ost << "VERTEX2";
                break;
            case Vertex::VERTEX3:
                ost << "VERTEX3";
                break;
            case Vertex::VERTEX4:
                ost << "VERTEX4";
                break;
            case Vertex::VERTEX5:
                ost << "VERTEX5";
                break;
            case Vertex::VERTEX6:
                ost << "VERTEX6";
                break;
            case Vertex::VERTEX7:
                ost << "VERTEX7";
                break;
        }
        return ost;
    }
    
    std::ostream& operator<<(std::ostream& ost, Edge e)
    {
        switch (e) {
            case Edge::EDGE0:
                ost << "EDGE0";
                break;
            case Edge::EDGE1:
                ost << "EDGE1";
                break;
            case Edge::EDGE2:
                ost << "EDGE2";
                break;
            case Edge::EDGE3:
                ost << "EDGE3";
            case Edge::EDGE4:
                ost << "EDGE4";
                break;
            case Edge::EDGE5:
                ost << "EDGE5";
                break;
            case Edge::EDGE6:
                ost << "EDGE6";
                break;
            case Edge::EDGE7:
                ost << "EDGE7";
            case Edge::EDGE8:
                ost << "EDGE8";
                break;
            case Edge::EDGE9:
                ost << "EDGE9";
                break;
            case Edge::EDGE10:
                ost << "EDGE10";
                break;
            case Edge::EDGE11:
                ost << "EDGE11";
                break;
        }
        return ost;
    }
    
    std::ostream& operator<<(std::ostream& ost, Sign s)
    {
        if(Sign::POSITIVE == s) ost << "+";
        else ost << "-";
        return ost;
    }
    
    
    std::ostream& operator<<( std::ostream& ost, const BCType& bc )
    {
        if( bc == DIRICHLET ) ost << "Dirichlet";
        else ost << "Neumann";
        return ost;
    }
    
    void endOfLoop( std::istream& ist, const char delim, const std::string& msg )
    {
        if( ist.fail() )
        {
            ist.clear();
            char ch;
            if( ist >> ch && ch == delim )
                return;
            error( msg );
        }
    }
    
    bool validParentCoord(const double u, const double v)
    {
        Interval i = boost::icl::construct<boost::icl::continuous_interval<double>>(-1.0,1.0);
        if(boost::icl::contains(i,u) && boost::icl::contains(i, v))
            return true;
        return false;
    }
    
    DoubleVec range(const double a, const double b, const uint N)
    {
        assert(N > 0);
        if (1 == N) {
            return {(a + b) / 2.0};
        }
        else {
            DoubleVec rvec;
            const double h = (b - a) / (N - 1);
            for(uint i = 0; i < N; ++i)
                rvec.push_back(i * h + a);
            return rvec;
        }
    }
    
    std::vector<long int> bounds(long int parts, long int mem)
    {
        std::vector<long int>bnd;
        long int delta = mem / parts;
        long int reminder = mem % parts;
        long int N1 = 0, N2 = 0;
        bnd.push_back(N1);
        for (long int i = 0; i < parts; ++i) {
            N2 = N1 + delta;
            if (i == parts - 1)
                N2 += reminder;
            bnd.push_back(N2);
            N1 = N2;
        }
        return bnd;
    }
}
