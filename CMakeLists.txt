project(trinurbs)
cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)

# look for the boost package
find_package(Boost REQUIRED)
if( Boost_FOUND )
include_directories(${Boost_INCLUDE_DIRS})
add_definitions( "-DHAS_BOOST" )
endif()

set(LIB_INCLUDE_DIR ${trinurbs_SOURCE_DIR}/include)
include_directories(${LIB_INCLUDE_DIR})

set(LIB_HEADERS ${LIB_INCLUDE_DIR}/base.h)

set(LIB_SOURCE_DIR ${trinurbs_SOURCE_DIR}/src/)
set(LIB_SOURCES ${LIB_SOURCE_DIR}/base.cpp
    ${LIB_HEADERS})

add_library(trinurbs ${LIB_SOURCES})

include(CTest)
add_subdirectory(tests)
enable_testing()


