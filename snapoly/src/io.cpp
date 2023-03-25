#include "pch.h"
#include "io.h"

// add OGRPolygon to polygons vector
void snapoly::io::add_OGRPolygon_to_polygons(
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
		polygon.outer_boundary().push_back(CDTPoint(poOuterRing->getX(i), poOuterRing->getY(i)));
		//cout << poOuterRing->getX(i) << "-" << poOuterRing->getY(i) << '\n';
	}

	// check
	//for (auto iter = polygon.outer_boundary().vertices_begin(); iter != polygon.outer_boundary().vertices_end(); ++iter)
	//	cout << iter->x() << "--" << iter->y() << '\n';

	// add inner rings
	Polygon_2 hole;
	for (int i = 0; i < poOGRPolygon->getNumInteriorRings(); ++i)
	{
		OGRLinearRing* innerRing = poOGRPolygon->getInteriorRing(i);
		int NumInnerPoints = innerRing->getNumPoints() - 1; // first point is the same as the last point
		for (int j = NumInnerPoints; j > 0; --j) // degenerate cases: NumOfouterRingPoints < 3?
		{
			hole.push_back(CDTPoint(innerRing->getX(j), innerRing->getY(j)));
		}polygon.holes().push_back(hole);
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
void snapoly::io::add_OGRMultiPolygon_to_polygons(
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
void snapoly::io::add_polygons_from_input_file(const char* input_file, vector<CDTPolygon>& polygons)
{
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
			char* spatialReference = nullptr;
			tmp->exportToWkt(&spatialReference);
			cout << "SpatialRef: " << spatialReference << '\n';
		}


		long long numberOfPolygons = poLayer->GetFeatureCount(true); // true some implementations will actually scan the entire layer once to count objects.
		cout << "number of polygons: " << numberOfPolygons << '\n';
		polygons.reserve(polygons.size() + numberOfPolygons);
		cout << "number of fields: " << poLayer->GetLayerDefn()->GetFieldCount() << '\n';


		// Save the polygons, field names and types
		for (const auto& poFeature : poLayer) // traverse all the features in one layer
		{
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

}
