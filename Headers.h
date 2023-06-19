#ifndef HEADER_FILE_H
#define HEADER_FILE_H

// standard library headers
#include <iostream> // for cout, etc.
#include <cstring> // for comparing two Cstrings
#include <limits> // for max and min values
#include <utility> // for std::pair
#include <queue> // for std::queue
#include <tuple> // for std::tuple
#include <algorithm> // for std::find
#include <chrono> // for Timer
#include <filesystem> // for checking whether the file has existed yet - cpp17 or higher

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h> // use exact_inexact kernel
#include <CGAL/Constrained_Delaunay_triangulation_2.h> // CDT
#include <CGAL/Triangulation_vertex_base_with_info_2.h> // for customized vertex
#include <CGAL/Triangulation_face_base_with_info_2.h> // for customized face
#include <CGAL/Polygon_2.h> // 2D polygon, for declaring 2d boundaries
#include <CGAL/Polygon_2_algorithms.h> // for 2D polygon related algorithms
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/squared_distance_2.h> //for 2D distance functions
#include <CGAL/Boolean_set_operations_2.h> // for polygon boolean operations
#include <CGAL/Segment_2.h> // for constructing segments

#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/CoordinateArraySequence.h>
#include <geos/geom/LineString.h>
#include <geos/geom/Polygon.h>
#include <geos/operation/polygonize/Polygonizer.h>

#include <ogrsf_frmts.h> // GDAL

#endif