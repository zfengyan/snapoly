
#include "pch.h"
#include "io.h"


/*
* if there are more than one polygons
* the output would be two polygons
* when we write to gpkg, we need to check the polygons' relationship
* i.e. which is the biggest polygon -> the biggest polygon should contain the other polygons
*/

// Attention
// geos namespace may contain same class name like Point
// which will probably be conflict with CDT typedefs
using namespace geos::geom;
using namespace geos::operation::polygonize;



int main()
{
	Snap_rounding_2 sr_2;

	Constraint c1(Kernel::Point_2(3.0, 2.0), Kernel::Point_2(3.0, 4.0));
	Constraint c2(Kernel::Point_2(5.0, 2.0), Kernel::Point_2(5.0, 4.0));
	Constraint c3(Kernel::Point_2(3.0, 4.0), Kernel::Point_2(5.0, 4.0));
	Constraint c4(Kernel::Point_2(3.0, 2.0), Kernel::Point_2(5.0, 2.0)); 

	Constraint c5(Kernel::Point_2(4.0, 2.5), Kernel::Point_2(3.5, 3.0));
	Constraint c6(Kernel::Point_2(4.5, 3.0), Kernel::Point_2(4.0, 3.5));
	Constraint c7(Kernel::Point_2(4.0, 2.5), Kernel::Point_2(4.5, 3.0));
	Constraint c8(Kernel::Point_2(3.5, 3), Kernel::Point_2(4.0, 3.5));

	Constraint c9(Kernel::Point_2(5.0, 2.0), Kernel::Point_2(6.0, 1.0));

	sr_2.constraintsWithInfo.emplace_back(c1);
	sr_2.constraintsWithInfo.emplace_back(c2);
	sr_2.constraintsWithInfo.emplace_back(c3);
	sr_2.constraintsWithInfo.emplace_back(c4);
	sr_2.constraintsWithInfo.emplace_back(c5);
	sr_2.constraintsWithInfo.emplace_back(c6);
	sr_2.constraintsWithInfo.emplace_back(c7);
	sr_2.constraintsWithInfo.emplace_back(c8);
	sr_2.constraintsWithInfo.emplace_back(c9);


	//////////////////////////////////////////////////////////////////////

	// for each constraint, two points -> two Coordinates objects

	//////////////////////////////////////////////////////////////////////
	
	
	vector<vector<Coordinate>> coordinates;

	for (auto const& cons : sr_2.constraintsWithInfo) {
		coordinates.emplace_back();
		coordinates.back().emplace_back(Coordinate(cons.p0.x(), cons.p0.y()));
		coordinates.back().emplace_back(Coordinate(cons.p1.x(), cons.p1.y()));
	}

	// coordinate sequences
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



	Kernel::Point_2 p0(1, 1), p1(1, 3);
	cout << CGAL::midpoint(p0, p1);

	//Point p;

	//Edge e;


	return 0;
}