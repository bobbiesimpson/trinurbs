#include "Forest.h"
#include "Geometry.h"

namespace trinurbs
{
    Forest::Forest(const Geometry& g) : mGeom(&g)
    {
        const Forest& f = g.primalForest();
        mSpaces = f.mSpaces;
        mSpaceMap = f.mSpaceMap;
        mNodalConn = f.mNodalConn;
        mGlobalDofN = f.mGlobalDofN;
        
        initEdgeConn();
        initFaceConn();
    }
    
    Forest::Forest(const Forest& f)
    {
        clear();
        mGeom = f.mGeom;
        mSpaces = f.mSpaces;
        mSpaceMap = f.mSpaceMap;
        mNodalConn = f.mNodalConn;
        mEdgeConn = f.mEdgeConn;
        mFaceConn = f.mFaceConn;
        mEdgeSpaceMap = f.mEdgeSpaceMap;
        mCVertexSpaceMap = f.mCVertexSpaceMap;
        mFaceSpaceMap = f.mFaceSpaceMap;
        mElemIndexMap = f.mElemIndexMap;
        
        /// TODO
        
//        for(const auto& e : f.mElems)
//            mElems.insert(std::make_pair(e.first, e.second->copy()));
//        for(const auto& e : f.mBezierElems)
//            mBezierElems.insert(std::make_pair(e.first, e.second->copy()));
        mGlobalDofN = f.mGlobalDofN;
    }
    
    /// Assignment operator
    Forest& Forest::operator=(const Forest& f)
    {
        if(this == &f) // check for self assignment
            return *this;
        clear();
        mGeom = f.mGeom;
        mSpaces = f.mSpaces;
        mSpaceMap = f.mSpaceMap;
        mNodalConn = f.mNodalConn;
        mEdgeConn = f.mEdgeConn;
        mFaceConn = f.mFaceConn;
        mEdgeSpaceMap = f.mEdgeSpaceMap;
        mCVertexSpaceMap = f.mCVertexSpaceMap;
        mFaceSpaceMap = f.mFaceSpaceMap;
        mElemIndexMap = f.mElemIndexMap;
        
        /// TODO
//        for(const auto& e : f.mElems)
//            mElems.insert(std::make_pair(e.first, e.second->copy()));
//        for(const auto& e : f.mBezierElems)
//            mBezierElems.insert(std::make_pair(e.first, e.second->copy()));
        
        mGlobalDofN = f.mGlobalDofN;
        return *this;
    }
    
    /// clear all data
    void Forest::clear()
    {
        mSpaces.clear();
        mSpaceMap.clear();
        mNodalConn.clear();
        mEdgeConn.clear();
        mFaceConn.clear();
        mEdgeSpaceMap.clear();
        mCVertexSpaceMap.clear();
        mFaceSpaceMap.clear();
        mElemIndexMap.clear();
        
        /// TODO
//        mElems.clear();
//        mBezierElems.clear();
        mGlobalDofN = std::make_pair(false, 0);
    }
    
    
    
    
    
    
    
    
    
    void Forest::load(std::istream& ist)
    {
        while(true) {
            BSplineSpace s;
            if(!(ist >> s))
                break;
            std::cout << "found space\n";
            addSpace(s);
        }
        endOfLoop(ist, '}', "Could not read B-spline space\n");
        
        while(true) {
            ConnVecInput cv;
            if(!(ist >> cv))
                break;
            addConnectivityVec(cv.name, cv.data);
        }
        endOfLoop(ist, '}', "Could not read connectivity vectors\n");
        initEdgeConn();
        initFaceConn();
    }
    
    void Forest::print(std::ostream& ost) const
    {
        ost << "Forest: " << spaceN() << " trunks\n";
        for(const auto& n : mSpaceMap) {
            ost << mSpaces[n.second] << "\n";
        }
        ost << "Nodal connectivity:\n";
        for(const auto& c : mNodalConn)
            ost << c.first << ": " << c.second << "\n";
    }

}