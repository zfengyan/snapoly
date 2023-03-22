#include "pch.h"
#include "Enhanced_constrained_delaunay_triangulation_2.h"


/*
* compare if two edges are the same
* 
* @param rhs: an rhs object which will be compared
* @return bool: return true if same by location, otherwise false
*/
bool Edge_with_info::operator==(const Edge_with_info& rhs) const
{
	return (
		((std::abs(this->x1 - rhs.x1) < snapoly::constant::epsilon) &&
			(std::abs(this->y1 - rhs.y1) < snapoly::constant::epsilon) &&
			(std::abs(this->x2 - rhs.x2) < snapoly::constant::epsilon) &&
			(std::abs(this->y2 - rhs.y2) < snapoly::constant::epsilon)) ||
		((std::abs(this->x1 - rhs.x2) < snapoly::constant::epsilon) &&
			(std::abs(this->y1 - rhs.y2) < snapoly::constant::epsilon) &&
			(std::abs(this->x2 - rhs.x1) < snapoly::constant::epsilon) &&
			(std::abs(this->y2 - rhs.y1) < snapoly::constant::epsilon))
		);
}



/*
// basic compute functions for Enhanced_constrained_delaunay_triangulation_2 class
// ====================================================================================================================================================================
*/

Kernel::FT Enhanced_constrained_delaunay_triangulation_2::squared_length(const Face_handle& f, int i) const
{
	int ccw = f->ccw(i);
	int cw = f->cw(i);
	return CGAL::squared_distance(f->vertex(ccw)->point(), f->vertex(cw)->point());
}

Kernel::FT Enhanced_constrained_delaunay_triangulation_2::squared_length(const Edge& edge) const
{
	Face_handle incident_face = edge.first;
	int opposite_vertex = edge.second;
	return squared_length(incident_face, opposite_vertex);
}

Kernel::FT Enhanced_constrained_delaunay_triangulation_2::area(const Point& p, const Point& q, const Point& r) const
{
	Kernel::FT signed_area = CGAL::area(p, q, r);
	return (signed_area > snapoly::constant::epsilon ? signed_area : -signed_area); // return the positive value
}

Kernel::FT Enhanced_constrained_delaunay_triangulation_2::area(const Face_handle& face) const
{
	if (is_infinite(face)) {
		cout << "The area of an infinite face is infinite. \n";
		cout << "std::numeric_limits<double>::max() will be returned. \n";
		return snapoly::constant::DOUBLE_MAX;
	}
	Point pa = face->vertex(0)->point();
	Point pb = face->vertex(1)->point();
	Point pc = face->vertex(2)->point();

	Kernel::FT signed_area = CGAL::area(pa, pb, pc);
	return (signed_area > snapoly::constant::epsilon ? signed_area : -signed_area); // return the positive value
}

Kernel::FT Enhanced_constrained_delaunay_triangulation_2::squared_height(const Face_handle& face, int i) const
{	
	auto area_ = area(face);// get the area according to the face handle
	auto squared_length_ = squared_length(face, i); // get the edge ledngth of the edge

	// get the squared height -> 1/2 * h * l = s -> h = 2s / l -> h^2 = 4 * s^2 / l^2
	auto squared_height_ = squared_length_ > snapoly::constant::epsilon ? ((4 * area_ * area_) / squared_length_) : snapoly::constant::DOUBLE_MAX; // _squared_length must > 0

	return squared_height_;
}
// ====================================================================================================================================================================


