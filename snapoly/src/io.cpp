#include "pch.h"
#include "io.h"

// definition of static variables -------------------------------------------------------------------
// definition of static class member
// definition allows the linker to find the memory space for the io::spatialReference variable.
OGRSpatialReference* io::m_spatialReference = nullptr;

// min values of the current extent of a layer
double io::minX = std::numeric_limits<double>::infinity();
double io::minY = std::numeric_limits<double>::infinity();

// add OGRPolygon to polygons vector
void io::add_OGRPolygon_to_polygons(
	const std::unique_ptr<OGRFeature, OGRFeatureUniquePtrDeleter>& poFeature,
	OGRPolygon* poOGRPolygon,
	vector<CDTPolygon>& polygons)
{
	// nullptr check
	if (poFeature == nullptr) {
		std::cerr << "the pointer of OGRFeature is empty, please check "
			<< "add_OGRPolygon_to_polygons() function in Snap_rounding.hpp file \n";
		return;
	}
	if (poOGRPolygon == nullptr) {
		std::cerr << "the OGRPolygon* pointer is empty, please check "
			<< "add_OGRPolygon_to_polygons() function in Snap_rounding.hpp file \n";
		return;
	}


	// add OGRPolygon
	int NumInteriorRings = poOGRPolygon->getNumInteriorRings(); // number of interior rings (holes)

	// for each OGRPolygon, make a Polygon object
	CDTPolygon polygon(NumInteriorRings);

	// add outer ring
	// NB: getNumPoints() will return the first point twice, e.g. for a square, getNumPoints() will be 5
	// and the outerRing(0) is the same point with outerRing(getNumPoints()-1)
	OGRLinearRing* poOuterRing = poOGRPolygon->getExteriorRing();
	if (poOuterRing == nullptr) {
		cout << "no exterior ring fetched, please check!\n";
		return;
	}

	int NumOuterPoints = poOuterRing->getNumPoints() - 1; // first point is the same as the last point

	//for (int i = 0; i < NumOuterPoints; ++i) // degenerate cases: NumOfouterRingPoints < 3?
	//{
	//	polygon.outer_boundary().push_back(Point(outerRing->getX(i), outerRing->getY(i)));
	//}

	for (int i = NumOuterPoints; i > 0; --i)
	{
		polygon.outer_boundary().push_back(CDTPoint(poOuterRing->getX(i) - minX, poOuterRing->getY(i) - minY));
		//cout << poOuterRing->getX(i) << "-" << poOuterRing->getY(i) << '\n';
	}

	// check
	//for (auto iter = polygon.outer_boundary().vertices_begin(); iter != polygon.outer_boundary().vertices_end(); ++iter)
	//	cout << iter->x() << "--" << iter->y() << '\n';
	//cout << "num of holes: " << poOGRPolygon->getNumInteriorRings() << '\n';

	// add inner rings
	Polygon_2 hole;
	for (int i = 0; i < poOGRPolygon->getNumInteriorRings(); ++i)
	{
		OGRLinearRing* innerRing = poOGRPolygon->getInteriorRing(i);
		int NumInnerPoints = innerRing->getNumPoints() - 1; // first point is the same as the last point
		for (int j = NumInnerPoints; j > 0; --j) // degenerate cases: NumOfouterRingPoints < 3?
		{
			hole.push_back(CDTPoint(innerRing->getX(j) - minX, innerRing->getY(j) - minY));
		}
		polygon.holes().push_back(hole);
		//sr.print_polygon_2(hole);
		hole.clear();
	}

	// check
	//for (auto const& hole : polygon.holes())sr.print_polygon_2(hole);

	// now let's store the fields (attributes)

	// cout << "Field name: " << oFeature->GetDefnRef()->GetName() << '\n'; // this is the layer name
	for (auto&& oField : *poFeature) // rvalue reference?
	{
		// field name - currently not added
		// cout << oField.GetName() << '\n';

		// if field unset or not exist 
		if (oField.IsUnset()) {
			printf("(unset),");
			continue;
		}
		if (oField.IsNull()) {
			printf("(null),");
			continue;
		}

		// get the type
		//const char* id = nullptr;
		switch (oField.GetType())
		{
		case OFTInteger:
			printf("%d,", oField.GetInteger());
			break;
		case OFTInteger64:
			printf(CPL_FRMT_GIB ",", oField.GetInteger64());
			break;
		case OFTReal:
			printf("%.3f,", oField.GetDouble());
			break;
		case OFTString:
			// GetString() returns a C string - const char*
			polygon.id() = string(oField.GetString()); // MUST use std::string
			break;
		default:
			// Note: we use GetAsString() and not GetString(), since
			// the later assumes the field type to be OFTString while the
			// former will do a conversion from the original type to string.
			printf("%s,", oField.GetAsString());
			break;
		}
	} // end for: all fields of a feature

	// now add this polygon to the polygons vector
	polygons.push_back(polygon);
} // end of function

// add OGRPolygon to polygons vector
void io::add_OGRMultiPolygon_to_polygons(
	const std::unique_ptr<OGRFeature, OGRFeatureUniquePtrDeleter>& poFeature,
	OGRMultiPolygon* poOGRMultiPolygon,
	vector<CDTPolygon>& polygons)
{
	// nullptr check
	if (poFeature == nullptr) {
		std::cerr << "the pointer of OGRFeature is empty, please check "
			<< "add_OGRMultiPolygon_to_polygons() function in Snap_rounding.hpp file \n";
		return;
	}
	if (poOGRMultiPolygon == nullptr) {
		std::cerr << "the OGRMultiPolygon* pointer is empty, please check "
			<< "add_OGRMultiPolygon_to_polygons() function in Snap_rounding.hpp file \n";
		return;
	}


	// add each OGRPolygon og this OGRMultiPolygon
	for (int currentPolygon = 0; currentPolygon < poOGRMultiPolygon->getNumGeometries(); ++currentPolygon) {
		OGRPolygon* poOGRPolygon = static_cast<OGRPolygon*>(poOGRMultiPolygon->getGeometryRef(currentPolygon));
		add_OGRPolygon_to_polygons(poFeature, poOGRPolygon, polygons);
	} // all OGRPolygons in one OGRMultiPolygon


	// note that inside the for loop, poFeature->getGeometryRef() will return a OGRMultiPolygon
}

// add polygons from the input file
void io::add_polygons_from_input_file(const char* input_file, vector<CDTPolygon>& polygons)
{
	cout << "reading polygons ... \n";


	GDALAllRegister();


	GDALDataset* poDataset = static_cast<GDALDataset*>(GDALOpenEx(input_file, GDAL_OF_VECTOR, NULL, NULL, NULL)); // get dataset
	if (poDataset == nullptr) {
		std::cerr << "Error: couldn't open file\n";
		return;
	}


	// prompt
	cout << "\tPath: " << input_file << endl;
	cout << "\tType: " << poDataset->GetDriverName() << endl;
	cout << "\tNum of Layers: " << poDataset->GetLayerCount() << endl;


	// all layers
	for (OGRLayer* poLayer : poDataset->GetLayers()) {
		
		poLayer->ResetReading(); // ensure we are starting at the beginning of the layer

		// get SpatialRef - will be saved and used for creating output files - how to store it?
		OGRSpatialReference* tmp = poLayer->GetSpatialRef();
		if (tmp != nullptr) {
			m_spatialReference = tmp->Clone(); // store the CRS information
		}


		long long numberOfPolygons = poLayer->GetFeatureCount(true); // true some implementations will actually scan the entire layer once to count objects.
		cout << "\tnumber of polygons: " << numberOfPolygons << '\n';
		polygons.reserve(polygons.size() + numberOfPolygons);
		cout << "\tnumber of fields: " << poLayer->GetLayerDefn()->GetFieldCount() << '\n';

		// get the extent - data type: double
		OGREnvelope extentOfLayer;
		poLayer->GetExtent(&extentOfLayer);
		cout << "extent: \n";
		cout << "min X: " << extentOfLayer.MinX << '\t' << "max X: " << extentOfLayer.MaxX << '\n';
		cout << "min Y: " << extentOfLayer.MinY << '\t' << "max Y: " << extentOfLayer.MaxY << '\n';
		minX = extentOfLayer.MinX;
		minY = extentOfLayer.MinY;
		// get the extent


		// Save the polygons, field names and types
		// int countFeature = 0; // for Debug
		for (const auto& poFeature : poLayer) // traverse all the features in one layer
		{
			//++countFeature;
			//cout << countFeature << '\n'; // for Debug

			OGRGeometry* poGeometry = poFeature->GetGeometryRef();
			if (poGeometry == nullptr) {
				std::cerr << "Warning: couldn't get the geometry of feature: \t" << poFeature->GetFID() << '\n';
				continue;
			}


			// major types would be multipolygon and polygon
			switch (poGeometry->getGeometryType())
			{
			case wkbPolygon: { // if the feature's geometry type is polygon
				OGRPolygon* poOGRPolygon = static_cast<OGRPolygon*>(poFeature->GetGeometryRef());
				add_OGRPolygon_to_polygons(poFeature, poOGRPolygon, polygons); // after this function, poFeature must exist, thus we use a const reference
				break;
			} // case: Polygon


			case wkbMultiPolygon: { // if the feature's geometry type is multipolygon
				OGRMultiPolygon* poOGRMultiPolygon = static_cast<OGRMultiPolygon*>(poFeature->GetGeometryRef());
				add_OGRMultiPolygon_to_polygons(poFeature, poOGRMultiPolygon, polygons);
				break;
			} // case: Polygon


			default:
				break;
			} // default: other types						

		} // end for: all features of a layer

	} // end for: all layers

	//cout << "vertices size: " << sr.et.number_of_vertices() << '\n';

	GDALClose(poDataset); // close the dataset

	cout << "done\n";

}

// output boundaries
void io::export_to_gpkg(const char* filename, const list<Constraint>& constraintsWithInfo)
{
	GDALAllRegister();

	// get driver
	const char* out_driver_name = "GPKG"; // output as .gpkg file
	GDALDriver* out_driver = GetGDALDriverManager()->GetDriverByName(out_driver_name);
	if (out_driver == nullptr) {
		std::cerr << "Error: OGR driver not found\n";
		return;
	}

	/* if file has already existed, overwrite * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	bool file_exist = std::filesystem::exists(std::filesystem::path(filename));
	if (file_exist) {
		std::cout << "file " << filename << " has already existed, overwriting ...\n";
		if (out_driver->Delete(filename) != OGRERR_NONE) {
			std::cerr << "Error: couldn't overwrite file\n";
			return;
		}
	}
	else std::cout << "Writing " << out_driver_name << " file " << filename << "...\n";
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	// get datasource
	const char* out_name = filename;
	GDALDataset* out_dataset = out_driver->Create(out_name, 0, 0, 0, GDT_Unknown, nullptr);
	if (out_dataset == nullptr) {
		std::cerr << "Error: couldn't create file: " << out_name << '\n';
		return;
	}

	// add constraints ---------------------------------------------------------------------------------------------------------
	std::cout << "-- output gpkg, write constraints with info" << '\n';
	// get layer for Edges
	OGRLayer* out_layer_edges = out_dataset->CreateLayer("edges", m_spatialReference);
	if (out_layer_edges == nullptr) {
		std::cerr << "Error: couldn't create layer - edges." << '\n';
		return;
	}
	// field for layer - edges, the first id
	OGRFieldDefn out_field_edges("osm_id", OFTString);
	out_field_edges.SetWidth(32);
	if (out_layer_edges->CreateField(&out_field_edges) != OGRERR_NONE) {
		std::cerr << "Error: Creating type field failed - layer edges" << '\n';
		return;
	}
	// for storing other ids - currently only support for at most two ids of a constraint
	OGRFieldDefn out_field_edges_1("osm_multiple_id", OFTString);
	out_field_edges_1.SetWidth(32);
	if (out_layer_edges->CreateField(&out_field_edges_1) != OGRERR_NONE) {
		std::cerr << "Error: Creating type field failed - layer edges" << '\n';
		return;
	}
	// for stroing the id numbers
	OGRFieldDefn out_field_edges_2("osm_id_numbers", OFTInteger);
	out_field_edges_2.SetWidth(32);
	if (out_layer_edges->CreateField(&out_field_edges_2) != OGRERR_NONE) {
		std::cerr << "Error: Creating type field failed - layer edges" << '\n';
		return;
	}

	// add constraints with info to OGRLineString
	for (auto const& c : constraintsWithInfo) {

		// - create local feature for each triangle face and set attribute (if any)
		OGRFeature* ogr_feature = OGRFeature::CreateFeature(out_layer_edges->GetLayerDefn());

		// for each constraint, the idCollection has a default value "unkown"
		// if they are attached with an ID, then the id will be idCollection[1]
		// if they have multiple IDs, then the idCollection.size() > 2
		ogr_feature->SetField("osm_id", c.idCollection[1].c_str()); // set attribute

		// if the constraint contains other id
		// now only supports for at most two ids for a constraint
		if (c.idCollection.size() > 2) {
			ogr_feature->SetField("osm_multiple_id", c.idCollection[2].c_str()); // set attribute
		}

		// store the number(amount) of id(s) - minus 1 because the default value "unknown"
		ogr_feature->SetField("osm_id_numbers", (int)c.idCollection.size() - 1); // set attribute

		// - create local geometry object
		OGRLineString ogr_line;
		ogr_line.addPoint(c.p0.x() + minX, c.p0.y() + minY);
		ogr_line.addPoint(c.p1.x() + minX, c.p1.y() + minY);

		// - set geometry
		ogr_feature->SetGeometry(&ogr_line);

		// - create feature in the file
		if (out_layer_edges->CreateFeature(ogr_feature) != OGRERR_NONE) {
			std::cout << "Error: couldn't create feature." << '\n';
			return;
		}
		OGRFeature::DestroyFeature(ogr_feature); // - clean up the local feature

	} // for loop: all edges
	// add edges ---------------------------------------------------------------------------------------------------------

	// close dataset
	GDALClose(out_dataset);

	// clean up
	GDALDestroyDriverManager();

	cout << "file saved at: " << filename << '\n';

}

// build output polygons from constraints
void io::build_polygons_from_constraints(
	list<Constraint>& constraintsWithInfo,
	vector<CDTPolygon>& resPolygonsVec)
{
	cout << "building polygons from constraints ... \n";
	
	// build an unordered_map and use the id of constraints as key, the corresponding constraints as value
	// i.e.constraintsMap[id] = constraints having the id "id"
	unordered_map<string, vector<Constraint>> constraintsMap;

	// populate the map
	// one constraint may contain multiple more than one id indicating it is a common boundary
	for (auto const& constraint : constraintsWithInfo) {
		auto number_of_id = constraint.idCollection.size();
		if (number_of_id > 1) { // the idCollection[0] is "unknown"
			for (int current_id = 1; current_id < number_of_id; ++current_id) {
				const string& id = constraint.idCollection[current_id];
				constraintsMap[id].push_back(constraint);
			} // end for: all ids
		} // end if: constraint.idCollection.size() > 1
	} // end for: each constraint in the constraints with info list

	//Debug
	/*for (auto const& element : constraintsMap) {
		cout << "id: " << element.first << '\n';
		cout << "constraints number: " << element.second.size() << '\n';
		for (auto const& constraint : element.second) {
			cout << constraint.p0.x() + minX << '\n';
		}cout << '\n';
	}*/
	//Debug

	// build CoordinateSequence for constraints with the same id
	// Coordinate: (x, y), for each constraint two points -> two Coordinates objects
	int count = 0; // Debug
	for (auto& element : constraintsMap) {

		//++count; // Debug
		//cout << count << '\n'; // Debug

		vector<vector<Coordinate>> coordinates;

		for (auto const& cons : element.second) { // element.second: vector<Constraint> -> store the constraints with same id
			coordinates.emplace_back();
			coordinates.back().emplace_back(Coordinate(cons.p0.x(), cons.p0.y()));
			coordinates.back().emplace_back(Coordinate(cons.p1.x(), cons.p1.y()));
		}

		// each LineString represents for one constraint
		std::size_t numOfVerticesOfLineString = 2;

		// createLineString accepts pointer of CoordinateSequence
		// the ownership of CoordinateSequence object will be assumed by global_factory->createLineString()
		// thus we don't need to free the memory manually
		vector<CoordinateSequence*> coordinateSequences;

		for (int i = 0; i < coordinates.size(); ++i) {
			CoordinateSequence* c = new CoordinateArraySequence(numOfVerticesOfLineString, 2);
			c->setPoints(coordinates[i]);
			coordinateSequences.push_back(c); // make a copy of the pointer
		}
		
		// create and store the line strings
		GeometryFactory::Ptr global_factory = GeometryFactory::create(); // unique_ptr, don't need to delete manually

		vector<Geometry*> geoms; //= new vector<Geometry*>; // need to delete later

		for (int i = 0; i < coordinateSequences.size(); ++i) {
			LineString* ls = global_factory->createLineString(coordinateSequences[i]); // takes the ownership of cl
			geoms.push_back(ls); // shallow copy, pointing to the same memory
		}

		Polygonizer pgnizer; // one ognizer for a set of constraints having same id
		pgnizer.add(&geoms);

		/*Ownership of vector is transferred to caller, subsequent
		* calls will return NULL.
		* @return a collection of Polygons*/
		std::vector<std::unique_ptr<GEOSPolygon>> polys = pgnizer.getPolygons();

		// it is IMPORTANT to check whether the polys vector is populated
		if (!polys.size()) {
			cout << "the constraints with osm_id: " << element.first
				<< " are not closed, will not form any polygon \n";
			continue; // continue to the next for loop
		}

		//cout << "polygon vector size: " << polys.size() << '\t';
		//cout << "has dangles? " << pgnizer.hasDangles() << '\n';
		// check
		//if (polys.size() != 1 || pgnizer.hasDangles())
			//cout << element.first << '\n';
		// check

		// if polygon(s) is(are) built


		// if polys.size() > 1 --------------------------------------------------------------------------------------------------
		// <1> the polygon contains hole(s), polys[0] is the exterior 
		// <2> there are more than one polygon with the same tag, e.g. overlapping area
		//  ___________a_______________
		// |          /|               |   e.g. vertex a, c, b belong to A2
		// |         / |               |        vertex a, b belong to A1
		// |  A1   c/A1|      A2       |        face a-c-b belongs to A1 and A2, which makes it an overlapping area
		// |        \A2|               |		there can be two possibilites when tagging -> see examples in the thesis
		// |         \ |               |        when polygonizer form the polygon A2, there will be two polygons tagged as A2
		// |__________\|_______________|
		//             b
		// ----------------------------------------------------------------------------------------------------------------------



		// <1> the polygon contains hole(s), polys[0] is the exterior
		// when adding exterior ring, align with the reading functions: add from backward
		if (polys[0]->getNumInteriorRing()) { // if there is(are) hole(s)
			
			// build resPolygon
			CDTPolygon resPolygon;
			resPolygon.id() = element.first;
											  									  
			// exterior ring --------------------------------------------------------------------------------------------------
			std::unique_ptr<CoordinateSequence> exteriorCoordSeq = polys[0]->getExteriorRing()->getCoordinates();
			std::size_t numOfExteriorPoints = exteriorCoordSeq->getSize() - 1; // last point is the same as the first

			// add points of exterior
			for (std::size_t i = numOfExteriorPoints; i > 0; --i) {
				const Coordinate& coord = exteriorCoordSeq->getAt(i);
				resPolygon.outer_boundary().push_back(CDTPoint(coord.x, coord.y));
				//std::cout << "(" << coord.x << ", " << coord.y << ")" << std::endl;
			}
			// add exterior points --------------------------------------------------------------------------------------------


			// when adding interior rings, the sequence of interior rings must be oriented CW not CCW
			// add interior rings ---------------------------------------------------------------------------------------------
			for (std::size_t i = 0; i < polys[0]->getNumInteriorRing(); ++i) {
				std::unique_ptr<CoordinateSequence> holeCoordSeq = polys[0]->getInteriorRingN(i)->getCoordinates();
				
				// represent a hole
				Polygon_2 hole;
				std::size_t numOfInteriorPoints = holeCoordSeq->getSize() - 1; // last point is the same as the first
				for (std::size_t j = numOfInteriorPoints; j > 0; --j) {
					const Coordinate& coord = holeCoordSeq->getAt(j);
					hole.push_back(CDTPoint(coord.x, coord.y));
					//std::cout << "(" << coord.x << ", " << coord.y << ")" << std::endl;
				} // end for: all points of an interior ring

				// add the hole to the CDTPolygon
				resPolygon.holes().push_back(hole);

			}// end for: all interior rings
			// add interior points --------------------------------------------------------------------------------------------

			// add built CDTPolygon to the vec
			resPolygonsVec.push_back(resPolygon);

		} // end if: if there is(are) hole(s)
		else {
					
			// <2> there is(are) no hole(s), there can be one or more than one polygon with the same tag, e.g. overlapping area
	        // each exterior ring for each polygon
			// we union them into one polygon - will union introduce cascading effects?

			//vector<Polygon_2> polygons_with_same_info;
			//vector<Polygon_with_holes_2> unions; // receive the union results

			//for (int i = 0; i < polys.size(); ++i) { // for each poly in polys
			//	// only add exterior
			//	// build Polygon_2 - for each poly
			//	Polygon_2 polygon;
			//	
			//	// exterior ring
			//	std::unique_ptr<CoordinateSequence> exteriorCoordSeq = polys[i]->getExteriorRing()->getCoordinates();
			//	std::size_t numOfExteriorPoints = exteriorCoordSeq->getSize() - 1; // last point is the same as the first

			//	// add points of exterior
			//	for (std::size_t i = numOfExteriorPoints; i > 0; --i) {
			//		const Coordinate& coord = exteriorCoordSeq->getAt(i);
			//		polygon.push_back(CDTPoint(coord.x, coord.y));
			//		//std::cout << "(" << coord.x << ", " << coord.y << ")" << std::endl;
			//	}
			//	//add exterior points

			//	// add built CDTPolygon to the vec
			//	polygons_with_same_info.push_back(polygon);

			//} // end for: each poly in polys	
				

			for (int i = 0; i < polys.size(); ++i) {
				// only add exterior
				// build resPolygon - for each poly
				CDTPolygon resPolygon;
				resPolygon.id() = element.first;

				// exterior ring
				std::unique_ptr<CoordinateSequence> exteriorCoordSeq = polys[i]->getExteriorRing()->getCoordinates();
				std::size_t numOfExteriorPoints = exteriorCoordSeq->getSize() - 1; // last point is the same as the first

				// add points of exterior
				for (std::size_t i = numOfExteriorPoints; i > 0; --i) {
					const Coordinate& coord = exteriorCoordSeq->getAt(i);
					resPolygon.outer_boundary().push_back(CDTPoint(coord.x, coord.y));
					//std::cout << "(" << coord.x << ", " << coord.y << ")" << std::endl;
				}
				//add exterior points

				// add built CDTPolygon to the vec
				resPolygonsVec.push_back(resPolygon);

			} // end for: each poly in polys	

		} // end else: no hole(s) in the polygon 
		

		// clean up
		for (unsigned int i = 0; i < geoms.size(); i++) {
			delete geoms[i];
		}

		
	} // end for: constraintsMap

	cout << "done \n";

}

// export the res polygons to gpkg file
void io::export_to_gpkg(const char* filename, vector<CDTPolygon>& resPolygonsVec)
{
	GDALAllRegister();

	// get driver
	const char* out_driver_name = "GPKG"; // output as .gpkg file
	GDALDriver* out_driver = GetGDALDriverManager()->GetDriverByName(out_driver_name);
	if (out_driver == nullptr) {
		std::cerr << "Error: OGR driver not found\n";
		return;
	}

	/* if file has already existed, overwrite * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	bool file_exist = std::filesystem::exists(std::filesystem::path(filename));
	if (file_exist) {
		std::cout << "file " << filename << " has already existed, overwriting ...\n";
		if (out_driver->Delete(filename) != OGRERR_NONE) {
			std::cerr << "Error: couldn't overwrite file\n";
			return;
		}
	}
	else std::cout << "Writing " << out_driver_name << " file " << filename << "...\n";
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	// get datasource
	const char* out_name = filename;
	GDALDataset* out_dataset = out_driver->Create(out_name, 0, 0, 0, GDT_Unknown, nullptr);
	if (out_dataset == nullptr) {
		std::cerr << "Error: couldn't create file: " << out_name << '\n';
		return;
	}

	// add faces ------------------------------------------------------------------------------------------------
	std::cout << "-- output gpkg, write polygons" << '\n';
	// get layer for polygons
	OGRLayer* out_layer_polygons = out_dataset->CreateLayer("polygons", m_spatialReference);
	if (out_layer_polygons == nullptr) {
		std::cerr << "Error: couldn't create layer - polygons." << '\n';
		return;
	}
	// field for layer - polygons
	OGRFieldDefn out_field_polygons("osm_id", OFTString);
	out_field_polygons.SetWidth(32);
	if (out_layer_polygons->CreateField(&out_field_polygons) != OGRERR_NONE) {
		std::cerr << "Error: Creating type field failed - layer polygons" << '\n';
		return;
	}

	// add faces - traverse all polygons
	for (auto const& pgn : resPolygonsVec) { // must use reference here: auto const& or Face_handle&

		// - create local feature for each triangle face and set attribute (if any)
		OGRFeature* ogr_feature = OGRFeature::CreateFeature(out_layer_polygons->GetLayerDefn());
		ogr_feature->SetField("osm_id", pgn.id().c_str()); // set attribute

		// - create local geometry object
		OGRPolygon* ogr_polygon = new OGRPolygon; // using new keyword - need to be deleted by using `delete` later

		// add exterior ring
		OGRLinearRing* ogr_exterior_ring = new OGRLinearRing;
		for (auto iter = pgn.outer_boundary().vertices_begin(); iter != pgn.outer_boundary().vertices_end(); ++iter) { // add points of exterior ring
			ogr_exterior_ring->addPoint(CGAL::to_double(iter->x() + minX), CGAL::to_double(iter->y() + minY));
		}ogr_exterior_ring->closeRings(); // ogr ring must be closed
		ogr_polygon->addRingDirectly(ogr_exterior_ring); // assumes ownership, no need to use 'delete' key word on the created (OGR) ring

		// add interior rings if any
		if (pgn.has_holes()) {
			for (auto const& hole : pgn.holes()) { // add each interior ring (hole)				
				OGRLinearRing* ogr_interior_ring = new OGRLinearRing; // for each hole, there must be a corresponding OGRLinearRing
				for (auto iter = hole.vertices_begin(); iter != hole.vertices_end(); ++iter) {
					ogr_interior_ring->addPoint(CGAL::to_double(iter->x() + minX), CGAL::to_double(iter->y() + minY));
				}ogr_interior_ring->closeRings();
				ogr_polygon->addRingDirectly(ogr_interior_ring);
			}
		}

		// - set geometry
		ogr_feature->SetGeometryDirectly(ogr_polygon); // assumes ownership, no need to use 'delete' key word on the created (OGR) polygon

		// - create feature in the file
		if (out_layer_polygons->CreateFeature(ogr_feature) != OGRERR_NONE) {
			std::cout << "Error: couldn't create feature." << '\n';
			return;
		}

		OGRFeature::DestroyFeature(ogr_feature); // - clean up the local feature

	} // end for: polygons
	// add faces ------------------------------------------------------------------------------------------------



	// add vertices ------------------------------------------------------------------------------------------------------
	//std::cout << "-- output gpkg, write polygon vertices" << '\n';
	//// get layer for vertices
	//OGRLayer* out_layer_vertices = out_dataset->CreateLayer("vertices");
	//if (out_layer_vertices == nullptr) {
	//	std::cerr << "Error: couldn't create layer - vertices." << '\n';
	//	return;
	//}

	//// field for layer vertices - point type
	//OGRFieldDefn out_field_vertices_type("type", OFTString);
	//out_field_vertices_type.SetWidth(32);
	//if (out_layer_vertices->CreateField(&out_field_vertices_type) != OGRERR_NONE) {
	//	std::cerr << "Error: Creating type field failed - layer vertices type" << '\n';
	//	return;
	//}

	//// field for layer vertices - point x coordinates
	//OGRFieldDefn out_field_vertices_x("X", OFTReal);
	//out_field_vertices_x.SetPrecision(10); // create a real field with a precision of 10 decimal places - can store values with up to 10 digits after the decimal point.
	//if (out_layer_vertices->CreateField(&out_field_vertices_x) != OGRERR_NONE) {
	//	std::cerr << "Error: Creating type field failed - layer vertices x coordinates" << '\n';
	//	return;
	//}

	//// field for layer vertices - point y coordinates
	//OGRFieldDefn out_field_vertices_y("Y", OFTReal);
	//out_field_vertices_y.SetPrecision(10); // create a real field with a precision of 10 decimal places - can store values with up to 10 digits after the decimal point.
	//if (out_layer_vertices->CreateField(&out_field_vertices_y) != OGRERR_NONE) {
	//	std::cerr << "Error: Creating type field failed - layer vertices x coordinates" << '\n';
	//	return;
	//}

	//// add vertices
	//for (auto const& pgn : resPolygonsVec) {

	//	// add points of the exterior ring
	//	for (auto iter = pgn.outer_boundary().vertices_begin(); iter != pgn.outer_boundary().vertices_end(); ++iter) { // add points of exterior ring
	//		OGRFeature* ogr_feature = OGRFeature::CreateFeature(out_layer_vertices->GetLayerDefn()); // - create local feature and set attribute (if any)
	//		ogr_feature->SetField("type", "point 2D"); // set attribute
	//		ogr_feature->SetField("X", CGAL::to_double(iter->x())); // set attribute
	//		ogr_feature->SetField("Y", CGAL::to_double(iter->y())); // set attribute

	//		OGRPoint ogr_pt; // create local geometry object
	//		ogr_pt.setX(CGAL::to_double(iter->x()));
	//		ogr_pt.setY(CGAL::to_double(iter->y()));
	//		ogr_feature->SetGeometry(&ogr_pt); // set geometry

	//		// - create feature in the file
	//		if (out_layer_vertices->CreateFeature(ogr_feature) != OGRERR_NONE) {
	//			std::cout << "Error: couldn't create feature." << '\n';
	//			return;
	//		}
	//		OGRFeature::DestroyFeature(ogr_feature); // - clean up the local feature
	//	}

	//	// add points of the interior ring
	//	if (pgn.has_holes()) {
	//		for (auto const& hole : pgn.holes()) {
	//			for (auto iter = hole.vertices_begin(); iter != hole.vertices_end(); ++iter) {
	//				OGRFeature* ogr_feature = OGRFeature::CreateFeature(out_layer_vertices->GetLayerDefn()); // - create local feature and set attribute (if any)
	//				ogr_feature->SetField("type", "point 2D"); // set attribute
	//				ogr_feature->SetField("X", CGAL::to_double(iter->x())); // set attribute
	//				ogr_feature->SetField("Y", CGAL::to_double(iter->y())); // set attribute

	//				OGRPoint ogr_pt; // create local geometry object
	//				ogr_pt.setX(CGAL::to_double(iter->x()));
	//				ogr_pt.setY(CGAL::to_double(iter->y()));
	//				ogr_feature->SetGeometry(&ogr_pt); // set geometry

	//				// - create feature in the file
	//				if (out_layer_vertices->CreateFeature(ogr_feature) != OGRERR_NONE) {
	//					std::cout << "Error: couldn't create feature." << '\n';
	//					return;
	//				}
	//				OGRFeature::DestroyFeature(ogr_feature); // - clean up the local feature
	//			}
	//		}
	//	}

	//}
	// add vertices ------------------------------------------------------------------------------------------------------

	// close dataset
	GDALClose(out_dataset);

	// clean up
	GDALDestroyDriverManager();

	cout << "file saved at: " << filename << '\n';

}

// export the triangulation to gpkg file
void io::export_to_gpkg(const char* filename, CDT& cdt)
{
	GDALAllRegister();


	// get driver
	const char* out_driver_name = "GPKG"; // output as .gpkg file
	GDALDriver* out_driver = GetGDALDriverManager()->GetDriverByName(out_driver_name);
	if (out_driver == nullptr) {
		std::cerr << "Error: OGR driver not found\n";
		return;
	}


	/* if file has already existed, overwrite * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	bool file_exist = std::filesystem::exists(std::filesystem::path(filename));
	if (file_exist) {
		std::cout << "file " << filename << " has already existed, overwriting ...\n";
		if (out_driver->Delete(filename) != OGRERR_NONE) {
			std::cerr << "Error: couldn't overwrite file\n";
			return;
		}
	}
	else std::cout << "Writing " << out_driver_name << " file " << filename << "...\n";
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


	// get datasource
	const char* out_name = filename;
	GDALDataset* out_dataset = out_driver->Create(out_name, 0, 0, 0, GDT_Unknown, nullptr);
	if (out_dataset == nullptr) {
		std::cerr << "Error: couldn't create file: " << out_name << '\n';
		return;
	}



	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	* for polygons / edges / vertices
	* stored in different layers
	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */



	// add triangle faces ------------------------------------------------------------------------------------------------
	std::cout << "-- output gpkg, write polygons" << '\n';
	// get layer for polygons
	OGRLayer* out_layer_polygons = out_dataset->CreateLayer("polygons", m_spatialReference);
	if (out_layer_polygons == nullptr) {
		std::cerr << "Error: couldn't create layer - polygons." << '\n';
		return;
	}
	// field for layer - polygons
	OGRFieldDefn out_field_polygons("id", OFTString);
	out_field_polygons.SetWidth(32);
	if (out_layer_polygons->CreateField(&out_field_polygons) != OGRERR_NONE) {
		std::cerr << "Error: Creating type field failed - layer polygons" << '\n';
		return;
	}

	// add triangle faces
	//int count = 0;
	for (auto const& fh : cdt.finite_face_handles()) { // must use reference here: auto const& or Face_handle&

		//std::cout << "count: " << count++ << '\n';

		// if the current fh contains at least one constrained edge
		/*bool contains_constrained = false;
		for (int i = 0; i < 3; ++i) {
			CDT::Edge cdt_ege(fh, i);
			if (cdt.is_constrained(cdt_ege))
				contains_constrained = true;
		}*/

		// - create local feature for each triangle face and set attribute (if any)
		OGRFeature* ogr_feature = OGRFeature::CreateFeature(out_layer_polygons->GetLayerDefn());
		
		// if the finite face belongs to a polygon, then faceid_collection.size() = 2
		// if not faceid_collection.size() = 1
		if(fh->info().faceid_collection.size() > 1)
			ogr_feature->SetField("id", fh->info().faceid_collection[1].c_str()); // set attribute - id
		else
			ogr_feature->SetField("id", fh->info().faceid_collection[0].c_str());

		// - create local geometry object
		OGRPolygon* ogr_polygon = new OGRPolygon; // using new keyword
		OGRLinearRing* ogr_ring = new OGRLinearRing;
		//std::cout << "--check: construct ring" << '\n';
		// when exporting, add minX and minY to export the coordinates in the original extent
		ogr_ring->addPoint(CGAL::to_double(fh->vertex(0)->point().x() + minX), CGAL::to_double(fh->vertex(0)->point().y() + minY));
		ogr_ring->addPoint(CGAL::to_double(fh->vertex(1)->point().x() + minX), CGAL::to_double(fh->vertex(1)->point().y() + minY));
		ogr_ring->addPoint(CGAL::to_double(fh->vertex(2)->point().x() + minX), CGAL::to_double(fh->vertex(2)->point().y() + minY));
		ogr_ring->addPoint(CGAL::to_double(fh->vertex(0)->point().x() + minX), CGAL::to_double(fh->vertex(0)->point().y() + minY));
		ogr_ring->closeRings(); // ogr ring must be closed
		ogr_polygon->addRingDirectly(ogr_ring); // assumes ownership, no need to use 'delete' key word of created OGRLinearRing
		//std::cout << "--check: finish construct ring" << '\n';

		// - set geometry
		ogr_feature->SetGeometryDirectly(ogr_polygon); // assumes ownership, no need to use 'delete' key word of created OGRPolygon
		//std::cout << "--check: set geometry" << '\n';

		// - create feature in the file
		if (out_layer_polygons->CreateFeature(ogr_feature) != OGRERR_NONE) {
			std::cout << "Error: couldn't create feature." << '\n';
			return;
		}
		//std::cout << "--check: create feature" << '\n';
		OGRFeature::DestroyFeature(ogr_feature); // - clean up the local feature
	} // end for: finite face handles
	// add triangle faces ------------------------------------------------------------------------------------------------


	// add vertices ------------------------------------------------------------------------------------------------------
	std::cout << "-- output gpkg, write vertices" << '\n';
	// get layer for vertices
	OGRLayer* out_layer_vertices = out_dataset->CreateLayer("vertices", m_spatialReference);
	if (out_layer_vertices == nullptr) {
		std::cerr << "Error: couldn't create layer - vertices." << '\n';
		return;
	}
	// field for layer - vertices
	OGRFieldDefn out_field_vertices("type", OFTString);
	out_field_vertices.SetWidth(32);
	if (out_layer_vertices->CreateField(&out_field_vertices) != OGRERR_NONE) {
		std::cerr << "Error: Creating type field failed - layer vertices" << '\n';
		return;
	}

	// add coordinates
	// field for layer vertices - point x coordinates
	OGRFieldDefn out_field_vertices_x("X", OFTReal);
	out_field_vertices_x.SetPrecision(10); // create a real field with a precision of 10 decimal places - can store values with up to 10 digits after the decimal point.
	if (out_layer_vertices->CreateField(&out_field_vertices_x) != OGRERR_NONE) {
		std::cerr << "Error: Creating type field failed - layer vertices x coordinates" << '\n';
		return;
	}

	// field for layer vertices - point y coordinates
	OGRFieldDefn out_field_vertices_y("Y", OFTReal);
	out_field_vertices_y.SetPrecision(10); // create a real field with a precision of 10 decimal places - can store values with up to 10 digits after the decimal point.
	if (out_layer_vertices->CreateField(&out_field_vertices_y) != OGRERR_NONE) {
		std::cerr << "Error: Creating type field failed - layer vertices x coordinates" << '\n';
		return;
	}

	// add vertices
	for (auto const& vh : cdt.finite_vertex_handles()) {

		// - create local feature for each triangle face and set attribute (if any)
		OGRFeature* ogr_feature = OGRFeature::CreateFeature(out_layer_vertices->GetLayerDefn());
		ogr_feature->SetField("type", "point 2D"); // set attribute
		ogr_feature->SetField("X", CGAL::to_double(vh->point().x() + minX)); // set attribute
		ogr_feature->SetField("Y", CGAL::to_double(vh->point().y() + minY)); // set attribute

		// - create local geometry object
		OGRPoint ogr_pt;
		ogr_pt.setX(CGAL::to_double(vh->point().x() + minX));
		ogr_pt.setY(CGAL::to_double(vh->point().y() + minY));

		// - set geometry
		ogr_feature->SetGeometry(&ogr_pt);

		// - create feature in the file
		if (out_layer_vertices->CreateFeature(ogr_feature) != OGRERR_NONE) {
			std::cout << "Error: couldn't create feature." << '\n';
			return;
		}
		OGRFeature::DestroyFeature(ogr_feature); // - clean up the local feature
	}
	// add vertices ------------------------------------------------------------------------------------------------------


	// add edges ---------------------------------------------------------------------------------------------------------
	std::cout << "-- output gpkg, write edges" << '\n';
	// get layer for Edges
	OGRLayer* out_layer_edges = out_dataset->CreateLayer("edges", m_spatialReference);
	if (out_layer_edges == nullptr) {
		std::cerr << "Error: couldn't create layer - edges." << '\n';
		return;
	}
	// field for layer - edges
	OGRFieldDefn out_field_edges("type", OFTString);
	out_field_edges.SetWidth(32);
	if (out_layer_edges->CreateField(&out_field_edges) != OGRERR_NONE) {
		std::cerr << "Error: Creating type field failed - layer edges" << '\n';
		return;
	}

	// add edges
	std::vector<Edge_with_info> edges;
	Edge_with_info edge; // one declaration for multiple usages
	for (Face_handle fh : cdt.finite_face_handles()) {
		for (int i = 0; i != 3; ++i) {

			// mark the constrained status
			CDT::Edge cdt_ege(fh, i);
			if (cdt.is_constrained(cdt_ege)) {
				edge.is_constrained = true;
			}
			else {
				edge.is_constrained = false;
			}

			// get the custom edge with (x1, y1, x2, y2)
			int cw = CDT::cw(i);
			int ccw = CDT::ccw(i);
			edge.x1 = CGAL::to_double(fh->vertex(cw)->point().x() + minX);
			edge.y1 = CGAL::to_double(fh->vertex(cw)->point().y() + minY);
			edge.x2 = CGAL::to_double(fh->vertex(ccw)->point().x() + minX);
			edge.y2 = CGAL::to_double(fh->vertex(ccw)->point().y() + minY);
			if (std::find(edges.begin(), edges.end(), edge) == edges.end())
				edges.emplace_back(edge); // if not existed
		} // for each face - traverse its three vertices and get the edge opposite to vertex i
	} // for loop: all finite face handles

	// add edges to OGRLineString
	for (auto const& e : edges) { // change const auto& to index-based?

		// - create local feature for each triangle face and set attribute (if any)
		OGRFeature* ogr_feature = OGRFeature::CreateFeature(out_layer_edges->GetLayerDefn());
		const char* constrained_field = e.is_constrained ? "constrained" : "non constrained";
		ogr_feature->SetField("type", constrained_field); // set attribute

		// - create local geometry object
		OGRLineString ogr_line;
		ogr_line.addPoint(e.x1, e.y1);
		ogr_line.addPoint(e.x2, e.y2);

		// - set geometry
		ogr_feature->SetGeometry(&ogr_line);

		// - create feature in the file
		if (out_layer_edges->CreateFeature(ogr_feature) != OGRERR_NONE) {
			std::cout << "Error: couldn't create feature." << '\n';
			return;
		}
		OGRFeature::DestroyFeature(ogr_feature); // - clean up the local feature

	} // for loop: all edges
	// add edges ---------------------------------------------------------------------------------------------------------

	// close dataset
	GDALClose(out_dataset);

	// clean up
	GDALDestroyDriverManager();

	cout << "file saved at: " << filename << '\n';

}


