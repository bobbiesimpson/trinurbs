project(trinurbs)
cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)

# look for the boost package
find_package(Boost REQUIRED)
if( Boost_FOUND )
include_directories(${Boost_INCLUDE_DIRS})
add_definitions( "-DHAS_BOOST" )
endif()

# look for VTK package
find_package( VTK REQUIRED )
include(${VTK_USE_FILE})
include_directories( ${VTK_INCLUDE_DIRS} )
set(THIS_LINK_LIBS ${THIS_LINK_LIBS} ${VTK_LIBRARIES})

# add Eigen library include fileds
include_directories(SYSTEM ${fast-ibem_SOURCE_DIR}/libs)

set(LIB_INCLUDE_DIR ${trinurbs_SOURCE_DIR}/include)
include_directories(${LIB_INCLUDE_DIR})

set(LIB_HEADERS
    ${LIB_INCLUDE_DIR}/base.h
    ${LIB_INCLUDE_DIR}/AnalysisElement.h
    ${LIB_INCLUDE_DIR}/BezierNodalElement.h
    ${LIB_INCLUDE_DIR}/BSplineSpace.h
    ${LIB_INCLUDE_DIR}/Forest.h
    ${LIB_INCLUDE_DIR}/Geometry.h
    ${LIB_INCLUDE_DIR}/GeometryElement.h
    ${LIB_INCLUDE_DIR}/IBaseIntegrate.h
    ${LIB_INCLUDE_DIR}/IElemIntegrate.h
    ${LIB_INCLUDE_DIR}/InputDataStructures.h
    ${LIB_INCLUDE_DIR}/IParentSample.h
    ${LIB_INCLUDE_DIR}/MultiscaleBezierNodalElement.h
    ${LIB_INCLUDE_DIR}/MultiscaleForest.h
    ${LIB_INCLUDE_DIR}/NURBSCache.h
    ${LIB_INCLUDE_DIR}/NURBSCommon.h
    ${LIB_INCLUDE_DIR}/OutputVTK.h
    ${LIB_INCLUDE_DIR}/Point.h
    ${LIB_INCLUDE_DIR}/Point3D.h
    ${LIB_INCLUDE_DIR}/Point4D.h)

set(LIB_SOURCE_DIR ${trinurbs_SOURCE_DIR}/src/)
set(LIB_SOURCES
    ${LIB_SOURCE_DIR}/base.cpp
    ${LIB_SOURCE_DIR}/BezierNodalElement.cpp
    ${LIB_SOURCE_DIR}/BSplineSpace.cpp
    ${LIB_SOURCE_DIR}/Forest.cpp
    ${LIB_SOURCE_DIR}/Geometry.cpp
    ${LIB_SOURCE_DIR}/GeometryElement.cpp
    ${LIB_SOURCE_DIR}/IElemIntegrate.cpp
    ${LIB_SOURCE_DIR}/MultiscaleBezierNodalElement.cpp
    ${LIB_SOURCE_DIR}/MultiscaleForest.cpp
    ${LIB_SOURCE_DIR}/NURBSCommon.cpp
    ${LIB_SOURCE_DIR}/OutputVTK.cpp
    ${LIB_SOURCE_DIR}/Point.cpp
    ${LIB_SOURCE_DIR}/Point3D.cpp
    ${LIB_SOURCE_DIR}/Point4D.cpp
    ${LIB_HEADERS})

add_library(trinurbs ${LIB_SOURCES})

include(CTest)
add_subdirectory(tests)
enable_testing()


