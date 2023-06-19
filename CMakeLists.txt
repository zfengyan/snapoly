# CMakeList.txt : CMake project for snpoly_cmake, include source and define
# project specific logic here.
#

cmake_minimum_required(VERSION 3.8)
project(snapoly)

# set cpp standard 17 - as #include<filesystem> requires minimum cpp 17
set(CMAKE_CXX_STANDARD 17)

# CGAL
find_package(CGAL)
if (CGAL_FOUND)
    include(${CGAL_USE_FILE})
    message(STATUS "CGAL found")
    message(STATUS "CGAL version -- " ${CGAL_VERSION})
else()
    message(SEND_ERROR "this code requires the CGAL library")
    return()
endif()

# gdal
find_package(GDAL)
if (GDAL_FOUND)
    message(STATUS "GDAL found")
    message(STATUS "GDAL version -- " ${GDAL_VERSION})
else()
    message(SEND_ERROR "this code requires the GDAL library")
    return()
endif()
include_directories(${GDAL_INCLUDE_DIR}) # include GDAL directories

# GEOS
find_package(GEOS REQUIRED)
if (GEOS_FOUND)
    message(STATUS "GEOS found")
    message(STATUS "GEOS version -- " ${GEOS_VERSION})
else()
    message(SEND_ERROR "this code requires the GEOS library")
    return()
endif()
include_directories(${GEOS_INCLUDE_DIR}) # include GEOS directories

# Add source to this project's executable.
add_executable (snapoly "snapoly.cpp" "Headers.h" "Enhanced_constrained_delaunay_triangulation_2.h" "Enhanced_constrained_delaunay_triangulation_2.cpp" "Snap_rounding_2.h" "Snap_rounding_2.cpp" "io.h" "io.cpp")

target_link_libraries(snapoly ${GDAL_LIBRARY} GEOS::geos) # link to GDAL and GEOS::geos




