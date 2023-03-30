#ifndef SNAP_ROUNDING_H
#define SNAP_ROUNDING_H


#include "Enhanced_constrained_delaunay_triangulation_2.h"


/*
* Timer struct for counting the run time
* ATTENTION: Timer struct requires the C++ header file: #include<chrono>
* this is provided by MSVC, version: 14.29.30133
* and this may not suitable for gcc/g++ or other compilers like clang
* if it is not compiled using other compilers, feel free to comment the Timer struct
* 
* The code and usage is inspired from an outstanding C++ series:
* https://www.youtube.com/watch?v=oEx5vGNFrLk
*/
struct Timer //for counting the time
{
	std::chrono::time_point<std::chrono::steady_clock>start, end;
	std::chrono::duration<float>duration;

	Timer() //set default value
	{
		start = end = std::chrono::high_resolution_clock::now();
		duration = end - start;
	}

	~Timer()
	{
		end = std::chrono::high_resolution_clock::now();
		duration = end - start;

		std::cout << "Time: " << duration.count() * 0.01666666 << "min\n"; // 1/60
	}
};


// Constraint struct
struct Constraint {
	CDTPoint p0;
	CDTPoint p1;
	vector<string> idCollection; // for common boundary, there will be two constraints (same location) with different id

	Constraint()
		: p0(), p1(), idCollection()
	{
		idCollection.emplace_back("unknown");
	}
	Constraint(const CDTPoint& m, const CDTPoint& n)
		: p0(m), p1(n), idCollection()
	{
		idCollection.emplace_back("unknown");
	}

	// compare if two constraints are same by location
	bool operator==(const Constraint& rhs) const {
		if (
			(p0 == rhs.p0 && p1 == rhs.p1) ||
			(p0 == rhs.p1 && p1 == rhs.p0)) {
			return true;
		}
		else return false;
	}

	// compare if two constraints are not same by location
	bool operator!=(const Constraint& rhs) const {
		return !(*this == rhs);
	}

	// merge the id collections
	// the idCollection of rhs will be merged into the idCollection of *this
	// the idCollection of rhs will not be changed
	// e.g.
	// idCollection of *this: [A, B], idCollection of rhs: [B, C]
	// after the merge operation
	// idCollection of *this: [A, B, C]
	void merge_id_collection(const vector<string>& rhs_idCollection);


};


// class for hash function - unordered_set<Constraint>
class ConstraintHashFunction {
public:
	// id is returned as hash function
	size_t operator()(const Constraint& t) const
	{
		//return t.id;
		return 0;
		// return x_hash ^ (y_hash << 1)
	}
};


// CDTPolygon class
// represent for a 2D polygon with hole(s) - similar like Polygon_with_holes_2<Kernel>
class CDTPolygon
{
public:
	// easier usage
	typedef vector<Polygon_2>::const_iterator Hole_const_iterator; // const iterator
	typedef vector<Polygon_2>::iterator Hole_iterator; // iterator

	// constructors
	CDTPolygon(int NumInteriorRings, const string& id)
		: m_outerRing(), m_innerRings(), m_id(id)
	{
		m_innerRings.reserve(NumInteriorRings);
	}

	CDTPolygon(int NumInteriorRings)
		: m_outerRing(), m_innerRings(), m_id()
	{
		m_innerRings.reserve(NumInteriorRings);
	}

	CDTPolygon(const string& id)
		: m_outerRing(), m_innerRings(), m_id(id)
	{}

	CDTPolygon()
		: m_outerRing(), m_innerRings(), m_id()
	{}

	// access
	Polygon_2& outer_boundary() { return m_outerRing; }
	const Polygon_2& outer_boundary() const { return m_outerRing; }

	vector<Polygon_2>& holes() { return m_innerRings; }
	const vector<Polygon_2>& holes() const { return m_innerRings; }

	Hole_iterator holes_begin() { return m_innerRings.begin(); }
	const Hole_const_iterator holes_begin() const { return m_innerRings.begin(); }

	Hole_iterator holes_end() { return m_innerRings.end(); }
	const Hole_const_iterator holes_end() const { return m_innerRings.end(); }

	// vertices_begin() returns const iterator
	const Polygon_2::Vertex_const_iterator exterior_begin() const { return m_outerRing.vertices_begin(); }
	const Polygon_2::Vertex_const_iterator exterior_end() const { return m_outerRing.vertices_end(); }

	// id
	string& id() { return m_id; }
	const string& id() const { return m_id; }

	// predicate
	bool has_holes() const { return (m_innerRings.size() != 0); }

	// getters
	int number_of_holes() const { return static_cast<int>(m_innerRings.size()); }

	int number_of_exterior_points() const
	{
		int count = 0;
		for (auto iter = m_outerRing.vertices_begin(); iter != m_outerRing.vertices_end(); ++iter)++count;
		return count;
	}

	/*
	* check whether a point lies inside a CDTPolygon (possibly with hole(s))
	* Point type must be Traits::Point_2
	*
	* for each polygon
	* let's assume each polygon has at least one hole
	* so we need two variables: outerSide and innerSide to decide whether a point lies in the polygon
	* the condition is:
	* outerSide == CGAL::ON_BOUNDED_SIDE -> inside the exterior ring
	* innerSide == CGAL::ON_UNBOUNDED_SIDE -> outside the EACH interior ring
	*/
	static bool is_point_inside_polygon(const CDTPoint& pt, const CDTPolygon& pgn);


	/*
	* calculate the area of a polygon
	* if containing any holes:
	* area = areaOfExterior - areaOfHoles
	*/
	double area() const;

protected:
	Polygon_2 m_outerRing; // exterior ring of a polygon
	vector<Polygon_2> m_innerRings; // possible holes in the polygon
	string m_id; // fields of a polygon
};


// Snap_rounding_2 class
// for performing snap rounding
class Snap_rounding_2 {
public:

	// initialize according to the declaraion: m_tolerance(0.3),
	Snap_rounding_2():
		m_tolerance(0.3),
		m_squared_tolerance(0.09),
		m_et(),
		m_constraintsWithInfo(),
		m_polygons(),
		m_result_polygons()
	{
		cout << "default snap rounding tolerance is: " << m_tolerance << '\n';
	}


	/*
	* set tolerance
	*/
	void set_tolerance(double tolerance_param);

	
	/*
	* build the constrained (Delaunay) triangulation based on the polygons provided
	* polygons: a std::vector to hold all input polygons(assume to be valid)
	*/
	void insert_polygons_to_triangulation();

	
	/*
	* assume a startingFace is found
	* starting from the startingFace, using a BFS search to tag all the faces of the polygon
	* tagg the triangles until a constrained edge is encountered, which means the edge could be
	* a. the outer boundary of the polygon
	* b. the inner boundary of the polygon (hole(s))
	*
	* also tag the constraints, store the constraints with info to a vector
	*
	* when using emplace_back(), extra care needs to be taken
	* for example, constructed elements may not be in the same order when they were added into the vector
	*
	* @param:
	* @param - startingFace: where we start the tagging process
	* @param - refPgn: the reference polygon
	*/
	void add_tag_to_one_polygon(Face_handle& startingFace, const CDTPolygon& refPgn);


	/*
	* add tag to the triangulation
	* input: polygons with possible hole(s) and attributes attached (type: CDTPolygon)
	* output:
	* (1) a triangulation with face info (which polygon it belongs to)
	* (2) and constraints with tags stored in a vector
	*
	* Note that when tagging the triangles, we need to consider the polygon with hole(s):
	* ONLY if a point lies inside the outer boundary of a polygon but outside all of the holes
	* will this point is considered as "inside" the polygon with hole(s)
	* For example, if a point is inside the exterior ring but on the boundary of the interiors (holes)
	* that is not considered as "inside"? (or should we?)
	*
	* the precision is 1e-15 -> the highest precision for double type
	* since we are not using exact constructions, thus the calculated centroid has the precision of double
	*/
	void add_tag_to_triangulation();


	/*
	* remove dangles in the constrainsWithInfo list
	* there can be more than one dangle at the same location
	*    a   b           centroid
	*    /\  /\           / | \ 
	*   /  \/  \    ->   /  |  \
	*   |      |         |     |
	* after snap rounding a and b, there will be redundant constraint dangles in the middle (connecting to the centroid)
	*/
	void remove_dangles();


	/*
	* snap close vertex to vertex
	* update constraintsWithInfo first and then modify the triangulation
	* @param: edgeOfVertexToVertex - indicates the edge which connects two close vertices
	*/
	void snap_vertex_to_vertex(Edge& edgeOfVertexToVertex);


	/*
	* snap close vertex to boundary
	* update constraintsWithInfo first and then modify the triangulation
	* @param:
	* sliverFace: a sliver triangle which has one constrained base
	* oppositeVertexIndex: Edge(sliverFace, oppositeVertexIndex) indicates the constrained base
	*/
	void snap_vertex_to_boundary(Face_handle& sliverFace, int capturingVertexIndex);


	/*
	* snap from the minimum, maybe two close vertices or vertex - edge
	* depends on the length -> each time snap the minimum
	*/
	void snap_rounding();


	/*
	* measure the distortions by comparing the changed area
	* std::abs(area_of_result_polygon - area_of_original_polygon) / area_of_result_polygon
	*/
	double measure_distortions() const;


	/*
	* check the minimum distance under the given tolerance between:
	* <1> point to point
	* <2> point to boundary
	* return the minimum distance
	*/
	double minimum_distance() const;


	/*
	* find tolerance
	* two types of distances:
	* <1> point to point distance
	* <2> point to boundary distance
	* store the distances in a priority queue
	*/
	void find_tolerance(std::priority_queue<double>& lengthQueue);


	/*
	* get polygons vector
	*/
	vector<CDTPolygon>& polygons() { return m_polygons; }
	const vector<CDTPolygon>& polygons() const { return m_polygons; }


	/*
	* get boundaries
	*/
	list<Constraint>& constraintsWithInfo() { return m_constraintsWithInfo; }
	const list<Constraint>& constraintsWithInfo() const { return m_constraintsWithInfo; }


	/*
	* get the triangulation
	*/
	Enhanced_triangulation& triangulation() { return m_et; }
	const Enhanced_triangulation& triangulation() const { return m_et; }


	/*
	* get the result polygons 
	*/
	vector<CDTPolygon>& result_polygons() { return m_result_polygons; }
	const vector<CDTPolygon>& result_polygons() const { return m_result_polygons; }


	/*
	* get the tolerance
	* tolerance must be modified using set_tolerance() function
	*/
	const double& tolerance()const { return m_tolerance; }

protected:
	/*
	* the sequence of member declarations must be aligned with the constructors
	*/
	double m_tolerance; // snap rounding tolerance
	double m_squared_tolerance; // squared tolerance

	Enhanced_triangulation m_et; // the enhanced constrained Delaunay triangulation 

	list<Constraint> m_constraintsWithInfo; // store the constraints with the id attached, constraints at common boundary will have multiple ids

	vector<CDTPolygon> m_polygons; // store the OGRPolygons
	vector<CDTPolygon> m_result_polygons; // store the polygons recovered from the constraints with info list
};


namespace snapoly {
	namespace printer {
		void print(const Polygon_2& polygon_2);
	}
}




#endif //SNAP_ROUNDING_H