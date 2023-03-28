#ifndef SNAPOLY_IO_H
#define SNAPOLY_IO_H

#include "Snap_rounding_2.h"

// Attention
// geos namespace may contain same class name like Point
// which will probably be conflict with CDT typedefs
using namespace geos::geom;
using namespace geos::operation::polygonize;
using GEOSPolygon = geos::geom::Polygon; // differentiate from CDTPolygon

class io {

public:

	// store the spatialReference information
	static OGRSpatialReference* m_spatialReference;

	// store the minX and minY of a layer for coordinates shifting - note that a dataset could have more than one layer
	// now just consider one layer in a dataset
	static double minX;
	static double minY;

public:

	/*
	* add OGRPolygon to polygons
	* the function accepts a const reference of OGRFeatureUniquePtr and based on which the associated information is added
	* and a OGRPolygon* pointer indicating which OGRPolygon is added
	* poFeature is used for adding the associated information
	*/
	static void add_OGRPolygon_to_polygons(
		const std::unique_ptr<OGRFeature, OGRFeatureUniquePtrDeleter>& poFeature,
		OGRPolygon* poOGRPolygon,
		vector<CDTPolygon>& polygons);


	/*
	* add OGRMultiPolygon to polygons
	* the function accepts a const reference of OGRFeatureUniquePtr and based on which the associated information is added
	* OGRMultiPolygon* indicates which OGRMultiPolygon is added
	*/
	static void add_OGRMultiPolygon_to_polygons(
		const std::unique_ptr<OGRFeature, OGRFeatureUniquePtrDeleter>& poFeature,
		OGRMultiPolygon* poOGRMultiPolygon,
		vector<CDTPolygon>& polygons);


	/*
	* reading from the file using GDAL and store the polygons (with attributes attached)
	* store in polygons vector
	*/
	static void add_polygons_from_input_file(const char* input_file, vector<CDTPolygon>& polygons);


	/*
	* export to GPKG file
	*/
	static void export_to_gpkg(const char* filename, const list<Constraint>& constraintsWithInfo);


	/*
	* using constraintsWithInfo to build output polygons
	* the output polygons will be stored into outputPolygonsVec
	*/
	static void build_polygons_from_constraints(
		list<Constraint>& constraintsWithInfo,
		vector<CDTPolygon>& resPolygonsVec);


	/*
	* output resPolygons directly
	* polygons: store the 2D polygon with hole(s) (type: Polygon)
	* currently attributes only contain id (std::string)
	*/
	static void export_to_gpkg(const char* filename, vector<CDTPolygon>& resPolygonsVec);


	/**
	* Output a constrained delaunay triangulation to a .gpkg(geopackage) file.
	* TODO: need to fix edge output issues - const auto& - solved using index-based
	* Output three layers: polygons(triangles), edges, vertices
	* @param filename: output file name, e.g. "example.gpkg"
	* @param cdt: constrained dealunay triangulation
	* @return bool: true for success otherwise false
	*/
	static void export_to_gpkg(const char* filename, CDT& cdt);



}; // class io



#endif // !SNAPOLY_IO_H

