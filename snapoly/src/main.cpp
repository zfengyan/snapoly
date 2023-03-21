
#include "pch.h"
#include "Snap_rounding_2.h"

using namespace geos::geom;
using namespace geos::operation::polygonize;


/*
* if there are more than one polygons
* the output would be two polygons
* when we write to gpkg, we need to check the polygons' relationship
* i.e. which is the biggest polygon -> the biggest polygon should contain the other polygons
*/

int main()
{
	
	
	vector<vector<Coordinate>> coordinates;

	// line string 1
	coordinates.emplace_back();
	coordinates.back().emplace_back(Coordinate(3, 2));
	coordinates.back().emplace_back(Coordinate(3, 4));

	// line string 2
	coordinates.emplace_back();
	coordinates.back().emplace_back(Coordinate(5, 2));
	coordinates.back().emplace_back(Coordinate(5, 4));

	// line string 3
	coordinates.emplace_back();
	coordinates.back().emplace_back(Coordinate(3, 4));
	coordinates.back().emplace_back(Coordinate(5, 4));

	// line string 4
	coordinates.emplace_back();
	coordinates.back().emplace_back(Coordinate(3, 2));
	coordinates.back().emplace_back(Coordinate(5, 2));

	// hole 1 - line string 1
	coordinates.emplace_back();
	coordinates.back().emplace_back(Coordinate(4, 2.5));
	coordinates.back().emplace_back(Coordinate(3.5, 3));

	// hole 1 - line string 2
	coordinates.emplace_back();
	coordinates.back().emplace_back(Coordinate(4.5, 3));
	coordinates.back().emplace_back(Coordinate(4, 3.5));

	// hole 1 - line string 3
	coordinates.emplace_back();
	coordinates.back().emplace_back(Coordinate(4, 2.5));
	coordinates.back().emplace_back(Coordinate(4.5, 3));

	// hole 1 - line string 4
	coordinates.emplace_back();
	coordinates.back().emplace_back(Coordinate(3.5, 3));
	coordinates.back().emplace_back(Coordinate(4, 3.5));

	// dangling line string
	coordinates.emplace_back();
	coordinates.back().emplace_back(Coordinate(5, 2));
	coordinates.back().emplace_back(Coordinate(6, 1));

	// coordinate sequences
	std::size_t coordinates_size = 2;

	// createLineString accepts pointer of CoordinateSequence
	// the ownership of CoordinateSequence object will be assumed by global_factory->createLineString()
	// thus we don't need to free the memory manually
	vector<CoordinateSequence*> coordinateSequences;



	for (int i = 0; i < coordinates.size(); ++i) {
		CoordinateSequence* c = new CoordinateArraySequence(coordinates_size, 2);
		c->setPoints(coordinates[i]);
		coordinateSequences.push_back(c); // make a copy of the pointer
	}
	
	
	GeometryFactory::Ptr global_factory = GeometryFactory::create(); // unique_ptr, don't need to delete manually

	vector<Geometry*> geoms; //= new vector<Geometry*>; // need to delete later

	for (int i = 0; i < coordinateSequences.size(); ++i) {
		LineString* ls = global_factory->createLineString(coordinateSequences[i]); // takes the ownership of cl
		geoms.push_back(ls); // shallow copy, pointing to the same memory
	}

	////////////////////////////////////////////////////

	// by far, the memory we need to explicitly call `delete` is:
	// delete each pointer in geoms vector -> it holds the memory for each LineString
	// delete geoms itself -> it holds the memory for the dynamically allocated vector

	////////////////////////////////////////////////////

	Polygonizer plgnizer;
	plgnizer.add(&geoms);

	/*Ownership of vector is transferred to caller, subsequent
	* calls will return NULL.
	* @return a collection of Polygons*/
	std::vector<std::unique_ptr<Polygon>> polys = plgnizer.getPolygons();
	cout << "polygon vector size: " << polys.size();
	cout << "has dangles? " << plgnizer.hasDangles() << '\n';

	for (unsigned int i = 0; i < polys.size(); i++) {
		cout << "polygon: " << polys[i]->toString() << '\n';
		cout << "area: " << polys[i]->getArea() << '\n';


		// Print the coordinates

		// exterior ring
		//std::unique_ptr<CoordinateSequence> exteriorCoordSeq = polys[i]->getExteriorRing()->getCoordinates();
		//for (std::size_t m = 0; m < exteriorCoordSeq->getSize(); ++m) {
		//	const Coordinate& coord = exteriorCoordSeq->getAt(m);
		//	std::cout << "(" << coord.x << ", " << coord.y << ")" << std::endl;
		//}

		//// interior ring(s)
		//cout << "number of inner rings: " << polys[i]->getNumInteriorRing() << '\n';
		//for (std::size_t m = 0; m < polys[i]->getNumInteriorRing(); ++m) {
		//	std::unique_ptr<CoordinateSequence> holeCoordSeq = polys[i]->getInteriorRingN(m)->getCoordinates();

		//	// Print the coordinates of the i-th hole
		//	std::cout << "Hole " << i << ":" << std::endl;
		//	for (std::size_t j = 0; j < holeCoordSeq->getSize(); ++j) {
		//		const Coordinate& coord = holeCoordSeq->getAt(j);
		//		std::cout << "(" << coord.x << ", " << coord.y << ")" << std::endl;
		//	}
		//}
	}


	// clean up
	for (unsigned int i = 0; i < geoms.size(); i++) {
		delete geoms[i];
	}
	//delete geoms;


	return 0;
}