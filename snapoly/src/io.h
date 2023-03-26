#ifndef SNAPOLY_IO_H
#define SNAPOLY_IO_H

#include "Snap_rounding_2.h"


namespace snapoly {

	namespace io {

		// Attention
		// geos namespace may contain same class name like Point
		// which will probably be conflict with CDT typedefs
		using namespace geos::geom;
		using namespace geos::operation::polygonize;
		using GEOSPolygon = geos::geom::Polygon; // differentiate from CDTPolygon


		/*
		* add OGRPolygon to polygons
		* the function accepts a const reference of OGRFeatureUniquePtr and based on which the associated information is added
		* and a OGRPolygon* pointer indicating which OGRPolygon is added
		* poFeature is used for adding the associated information
		*/
		void add_OGRPolygon_to_polygons(
			const std::unique_ptr<OGRFeature, OGRFeatureUniquePtrDeleter>& poFeature,
			OGRPolygon* poOGRPolygon,
			vector<CDTPolygon>& polygons);


		/*
		* add OGRMultiPolygon to polygons
		* the function accepts a const reference of OGRFeatureUniquePtr and based on which the associated information is added
		* OGRMultiPolygon* indicates which OGRMultiPolygon is added
		*/
		void add_OGRMultiPolygon_to_polygons(
			const std::unique_ptr<OGRFeature, OGRFeatureUniquePtrDeleter>& poFeature,
			OGRMultiPolygon* poOGRMultiPolygon,
			vector<CDTPolygon>& polygons);


		/*
		* reading from the file using GDAL and store the polygons (with attributes attached)
		* store in polygons vector
		*/
		void add_polygons_from_input_file(const char* input_file, vector<CDTPolygon>& polygons);


		/*
		* export to GPKG file
		*/
		void export_to_gpkg(const char* filename, const list<Constraint>& constraintsWithInfo);


		/*
		* using constraintsWithInfo to build output polygons
		* the output polygons will be stored into outputPolygonsVec
		*/
		void build_polygons_from_constraints(
			list<Constraint>& constraintsWithInfo, 
			vector<CDTPolygon>& resPolygonsVec);


		/*
		* output resPolygons directly
		* polygons: store the 2D polygon with hole(s) (type: Polygon)
		* currently attributes only contain id (std::string)
		*/
		void export_to_gpkg(const char* filename, vector<CDTPolygon>& resPolygonsVec);


	} // namespace io
} // namespace snapoly








#endif // !SNAPOLY_IO_H

