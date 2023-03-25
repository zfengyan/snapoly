#ifndef SNAP_ROUNDING_H
#define SNAP_ROUNDING_H

#include "Enhanced_constrained_delaunay_triangulation_2.h"

typedef CDT::Point Point;

// Constraint
struct Constraint {
	CDTPoint p0;
	Kernel::Point_2 p1;
	string id; // for common boundary, there will be two constraints (same location) with different id

	Constraint()
		: p0(), p1(), id("unkown")
	{}
	Constraint(const Kernel::Point_2& m, const Kernel::Point_2& n)
		: p0(m), p1(n), id("unkown")
	{}

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


struct Snap_rounding_2 {

	list<Constraint> constraintsWithInfo;

	Snap_rounding_2()
		: constraintsWithInfo()
	{}

};




#endif //SNAP_ROUNDING_H