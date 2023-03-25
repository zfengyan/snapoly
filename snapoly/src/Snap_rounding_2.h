#ifndef SNAP_ROUNDING_H
#define SNAP_ROUNDING_H


#include "Enhanced_constrained_delaunay_triangulation_2.h"


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
	string& getID() { return m_id; }
	const string& getID() const { return m_id; }

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

protected:
	Polygon_2 m_outerRing; // exterior ring of a polygon
	vector<Polygon_2> m_innerRings; // possible holes in the polygon
	string m_id; // fields of a polygon
};


// Snap_rounding_2 class
// for performing snap rounding
class Snap_rounding_2 {

	list<Constraint> constraintsWithInfo;

	Snap_rounding_2()
		: constraintsWithInfo()
	{}

};




#endif //SNAP_ROUNDING_H