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
		((std::abs(this->x1 - rhs.x1) < snapoly::constants::EPSILON) &&
			(std::abs(this->y1 - rhs.y1) < snapoly::constants::EPSILON) &&
			(std::abs(this->x2 - rhs.x2) < snapoly::constants::EPSILON) &&
			(std::abs(this->y2 - rhs.y2) < snapoly::constants::EPSILON)) ||
		((std::abs(this->x1 - rhs.x2) < snapoly::constants::EPSILON) &&
			(std::abs(this->y1 - rhs.y2) < snapoly::constants::EPSILON) &&
			(std::abs(this->x2 - rhs.x1) < snapoly::constants::EPSILON) &&
			(std::abs(this->y2 - rhs.y1) < snapoly::constants::EPSILON))
		);
}



/*
// basic compute functions for Enhanced_constrained_delaunay_triangulation_2 class
// ====================================================================================================================================================================
*/

// numerical computations
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

Kernel::FT Enhanced_constrained_delaunay_triangulation_2::area(
	const CDTPoint& p, const CDTPoint& q, const CDTPoint& r) const
{
	Kernel::FT signed_area = CGAL::area(p, q, r);
	return (signed_area > snapoly::constants::EPSILON ? signed_area : -signed_area); // return the positive value
}

Kernel::FT Enhanced_constrained_delaunay_triangulation_2::area(const Face_handle& face) const
{
	if (is_infinite(face)) {
		cout << "The area of an infinite face is infinite. \n";
		cout << "std::numeric_limits<double>::max() will be returned. \n";
		return snapoly::constants::DOUBLE_MAX;
	}
	CDTPoint pa = face->vertex(0)->point();
	CDTPoint pb = face->vertex(1)->point();
	CDTPoint pc = face->vertex(2)->point();

	Kernel::FT signed_area = CGAL::area(pa, pb, pc);
	return (signed_area > snapoly::constants::EPSILON ? signed_area : -signed_area); // return the positive value
}

Kernel::FT Enhanced_constrained_delaunay_triangulation_2::squared_height(const Face_handle& face, int i) const
{	
	auto area_ = area(face);// get the area according to the face handle
	auto squared_length_ = squared_length(face, i); // get the edge ledngth of the edge

	// get the squared height -> 1/2 * h * l = s -> h = 2s / l -> h^2 = 4 * s^2 / l^2
	auto squared_height_ = squared_length_ > snapoly::constants::EPSILON ? ((4 * area_ * area_) / squared_length_) : snapoly::constants::DOUBLE_MAX; // _squared_length must > 0

	return squared_height_;
}


// face centroids
CDTPoint Enhanced_constrained_delaunay_triangulation_2::face_centroid(
	const CDTPoint& pa, const CDTPoint& pb, const CDTPoint& pc) const
{
	return CGAL::centroid(pa, pb, pc);
}

CDTPoint Enhanced_constrained_delaunay_triangulation_2::face_centroid(const Vertex_handle& va, const Vertex_handle& vb, const Vertex_handle& vc) const
{
	return CGAL::centroid(va->point(), vb->point(), vc->point());
}

CDTPoint Enhanced_constrained_delaunay_triangulation_2::face_centroid(const Face_handle& face) const
{
	if (is_infinite(face)) // if infinite face
	{
		std::cout << "the centroid of a face can not be calculated if the face is infinite"
			<< "an infinite CDTPoint will be returned" << '\n';
		return infinite_vertex()->point();
	}
	Vertex_handle va = face->vertex(0);
	Vertex_handle vb = face->vertex(1);
	Vertex_handle vc = face->vertex(2);
	return face_centroid(va, vb, vc);
}

// edge centroid
CDTPoint Enhanced_constrained_delaunay_triangulation_2::edge_centroid(const Edge& edge) const
{
	// get two vertices of this edge
	Face_handle incident_face = edge.first;
	int opposite_vertex = edge.second;
	int ccw = incident_face->ccw(opposite_vertex);
	int cw = incident_face->cw(opposite_vertex);
	Vertex_handle v_ccw = incident_face->vertex(ccw);
	Vertex_handle v_cw = incident_face->vertex(cw);

	return CGAL::midpoint(v_ccw->point(), v_cw->point());
}


// get the vertices of a edge
std::pair<Vertex_handle, Vertex_handle> Enhanced_constrained_delaunay_triangulation_2::vertices_of_edge(
	const Face_handle& incidentFace, int oppositeVertex) const
{
	int ccw = incidentFace->ccw(oppositeVertex);
	int cw = incidentFace->cw(oppositeVertex);
	Vertex_handle v_ccw = incidentFace->vertex(ccw);
	Vertex_handle v_cw = incidentFace->vertex(cw);

	return std::make_pair(v_ccw, v_cw);
}

// get the vertices of a edge
std::pair<Vertex_handle, Vertex_handle> Enhanced_constrained_delaunay_triangulation_2::vertices_of_edge(const Edge& edge) const
{
	// get incident face and opposite vertex of this edge
	Face_handle incidentFace = edge.first;
	int oppositeVertex = edge.second;
	return vertices_of_edge(incidentFace, oppositeVertex);
}

// ====================================================================================================================================================================






/*
// functions for supporting snap rounding, detailed descriptions please refer to Enhanced_constrained_delaunay_triangulation_2.h file
// ====================================================================================================================================================================
*/

// get constrained incident vertices of a certain vertex excluding an omit vertex
void Enhanced_constrained_delaunay_triangulation_2::get_constrained_incident_vertices(
	const Vertex_handle& targetVh, vector<CDTPoint>& constrained_incident_vertices_vec, Vertex_handle omitVh) const
{
	if (targetVh == omitVh)
	{
		std::cout << "the omitting vertex handle can not be the same as the target vertex handle " << '\n';
		return;
	}

	// get incident edges of target vertex handle
	Edge_circulator ec = incident_edges(targetVh), done(ec);
	if (ec != nullptr)
	{
		// loop through each edge
		do {
			Edge edge;
			edge.first = ec->first;
			edge.second = ec->second; // manually convert Edge_circulator to Edge

			if (is_constrained(edge)) // process the constraints
			{
				// get the incident face of each incident edge and the opposite vertex index
				Face_handle incident_face = edge.first;
				int opposite_vertex = edge.second;
				int cw = incident_face->cw(opposite_vertex);
				int ccw = incident_face->ccw(opposite_vertex);

				// for each incident edge, there are two vertices
				Vertex_handle incidentEdge_vh_a = incident_face->vertex(ccw);
				Vertex_handle incidentEdge_vh_b = incident_face->vertex(cw);

				// add constrained incident vertices to the vector
				// the vertices of the constrained incident edge must not be the same with targetVh or omitVh
				if (incidentEdge_vh_a != targetVh && incidentEdge_vh_a != omitVh)
					constrained_incident_vertices_vec.push_back(incidentEdge_vh_a->point());
				if (incidentEdge_vh_b != targetVh && incidentEdge_vh_b != omitVh)
					constrained_incident_vertices_vec.push_back(incidentEdge_vh_b->point());
			}

		} while (++ec != done);
	}
}

// if two edges are the same by location
bool Enhanced_constrained_delaunay_triangulation_2::is_same_edge_by_location(const Edge& edgeA, const Edge& edgeB) const
{
	Vertex_handle va = vertices_of_edge(edgeA).first;
	Vertex_handle vb = vertices_of_edge(edgeA).second;
	Vertex_handle vc = vertices_of_edge(edgeB).first;
	Vertex_handle vd = vertices_of_edge(edgeB).second;

	// operator== of Vertex handle only compares the underlying CDTPointers, not coordinates
	// if comparison of coordinates is required, operator== of CDTPoint should be used
	return (
		((va->point() == vc->point()) && (vb->point() == vd->point())) ||
		((va->point() == vd->point()) && (vb->point() == vc->point()))
		);
}

// find common edge
Edge Enhanced_constrained_delaunay_triangulation_2::common_edge(const Face_handle& faceA, const Face_handle& faceB) const
{
	Edge e; // pair<Face_handle, int i>
	e.second = snapoly::constants::NOT_EXIST;

	if (!faceA->has_neighbor(faceB)) {
		cout << "two faces are not adjacent \n";
		cout << "An edge with index -1 will be returned \n";
		return e;
	}

	// if they are neighors, find the common edge

	// find edges
	Edge eA[3];
	for (int i = 0; i < 3; ++i) {
		eA[i] = Edge(faceA, i);
	}

	Edge eB[3];
	for (int i = 0; i < 3; ++i) {
		eB[i] = Edge(faceB, i);
	}

	// find common edge
	for (auto const& ea : eA) {
		for (auto const& eb : eB) {
			if (is_same_edge_by_location(ea, eb)) { // comparing the underlying CDTPoint
				return ea;
			}
		}
	}

	// if no common edge is found
	cout << "no common edge is found between the two faces: \n";
	cout << "An edge with index -1 will be returned \n";
	return e;
}

// remove an edge from the triangulation
void Enhanced_constrained_delaunay_triangulation_2::remove_edge(Edge& edge)
{
	// get two vertices of this edge
	Face_handle incident_face = edge.first;
	int opposite_vertex = edge.second;
	int ccw = incident_face->ccw(opposite_vertex);
	int cw = incident_face->cw(opposite_vertex);
	Vertex_handle v_ccw = incident_face->vertex(ccw);
	Vertex_handle v_cw = incident_face->vertex(cw);

	// remove the two vertices of this edge
	remove(v_ccw);
	remove(v_cw);
}

// if two vertices connected by Edge(f, opposite_vertex) are close within the tolerance
bool Enhanced_constrained_delaunay_triangulation_2::is_degenerate_edge(const Face_handle& f, int opposite_vertex, double squared_tolerance) const
{
	Edge edge;
	edge.first = f;
	edge.second = opposite_vertex;

	int ccw = f->ccw(opposite_vertex);
	int cw = f->cw(opposite_vertex);
	auto squared_length = CGAL::squared_distance(f->vertex(ccw)->point(), f->vertex(cw)->point()); // return type: Kernel::FT

	return ((squared_length + snapoly::constants::EPSILON) < squared_tolerance || std::abs(squared_length - squared_tolerance) < snapoly::constants::EPSILON);
}

bool Enhanced_constrained_delaunay_triangulation_2::is_degenerate_edge(const Edge& edge, double squared_tolerance) const
{
	Face_handle incident_face = edge.first;
	int opposite_vertex = edge.second;
	return is_degenerate_edge(incident_face, opposite_vertex, squared_tolerance);
}


// find the edge with minimum length
std::tuple<Face_handle, int, Kernel::FT> Enhanced_constrained_delaunay_triangulation_2::find_minimum_degenerate_edge(double squared_tolerance) const
{
	Edge e;
	e.first = infinite_face();
	e.second = snapoly::constants::NOT_EXIST; // indicates this edge does not exist

	Kernel::FT min_squared_length = snapoly::constants::DOUBLE_MAX;

	// find the minimum edge
	for (auto& edge : finite_edges()) {
		bool is_found = is_degenerate_edge(edge, squared_tolerance);
		if (is_found) { // if a degenerate constrained edge is found
			auto computed_squared_length = squared_length(edge); // compute the squared length of this edge
			if ((computed_squared_length + snapoly::constants::EPSILON) < min_squared_length) { // if computed squared length is smaller than current minimum squared length
				min_squared_length = computed_squared_length; // update minimum squared length
				e = edge; // update e
			}
			//cout << "length: " << computed_length << '\n';
		}
	}//cout << "the minimum squared length is: " << min_squared_length << '\n';

	return std::make_tuple(e.first, e.second, min_squared_length);
}

// better alignment with concepts
std::tuple<Face_handle, int, Kernel::FT> Enhanced_constrained_delaunay_triangulation_2::find_minimum_vertex_to_vertex(double squared_tolerance) const
{
	return find_minimum_degenerate_edge(squared_tolerance);
}

// if an edge is a sliver base edge (or a constrained sliver base edge)
bool Enhanced_constrained_delaunay_triangulation_2::is_sliver_base(const Face_handle& face, int i, double squared_tolerance, bool constrained_flag) const
{
	// get the height -> based on the edge (face, i)
	auto squared_height_ = squared_height(face, i);

	// build the possibly sliver base edge
	Edge possibly_sliver_base_edge;
	possibly_sliver_base_edge.first = face;
	possibly_sliver_base_edge.second = i;


	// if the base edge's height <= tolerance
	// (height + epsilon) < tolerance indicates height < tolerance
	// std::abs(height - tolerance) < epsilon indicates that height can be considered equal to tolerance (=)
	bool condition_1 = false;

	// if the projection of the opposite vertex is on the constrained base edge
	bool condition_2 = false;

	// if the constrained_flag is specified as true
	if (constrained_flag) {
		condition_1 = ((squared_height_ + snapoly::constants::EPSILON) < squared_tolerance ||
			std::abs(squared_height_ - squared_tolerance) < snapoly::constants::EPSILON) &&
			is_constrained(possibly_sliver_base_edge);
	}
	else {
		condition_1 = ((squared_height_ + snapoly::constants::EPSILON) < squared_tolerance ||
			std::abs(squared_height_ - squared_tolerance) < snapoly::constants::EPSILON);
	}

	// if the projection of the opposite vertex is on the base edge
	// if the projection is within the base edge, then squared_distance will be the same as squared_height
	// otherwise CGAL::squared_distance() calculates the distance between a CDTPoint and the nearest end CDTPoint of a segment
	Segment_2 segment(vertices_of_edge(face, i).first->point(), vertices_of_edge(face, i).second->point());
	Kernel::FT squared_distance_ = CGAL::squared_distance(face->vertex(i)->point(), segment);	
	auto difference = std::abs(squared_height_ - squared_distance_); // type of return value of std::abs(): double
	if (difference < snapoly::constants::EPSILON)
		condition_2 = true;

	return (condition_1 && condition_2);
}

// if a face is a sliver triangle
std::pair<bool, int> Enhanced_constrained_delaunay_triangulation_2::is_sliver_triangle(const Face_handle& face, double squared_tolerance) const
{
	if (is_infinite(face)) // if infinite face
	{
		std::cout << "this is an infinite face, a sliver triangle must be a finite face" << '\n';
		return std::make_pair(false, 0);
	}

	// ONLY one sliver constrained base
	// the other two edges must not be sliver base edges (they may be constrained)

	// Edge (face, 0) is the constrained sliver base
	if (is_sliver_base(face, 0, squared_tolerance, true) &&
		!is_sliver_base(face, 1, squared_tolerance) &&
		!is_sliver_base(face, 2, squared_tolerance)) {
		return std::make_pair(true, 0);
	}

	// Edge (face, 1) is the constrained sliver base
	if (is_sliver_base(face, 1, squared_tolerance, true) &&
		!is_sliver_base(face, 2, squared_tolerance) &&
		!is_sliver_base(face, 0, squared_tolerance)) {
		return std::make_pair(true, 1);
	}

	// Edge (face, 2) is the constrained sliver base
	if (is_sliver_base(face, 2, squared_tolerance, true) &&
		!is_sliver_base(face, 0, squared_tolerance) &&
		!is_sliver_base(face, 1, squared_tolerance)) {
		return std::make_pair(true, 2);
	}

	// if the face is not a sliver triangle
	return std::make_pair(false, 0);
}

// find the sliver triangle with the minimum height (within the tolerance) in the triangulation
std::tuple<Face_handle, int, Kernel::FT> Enhanced_constrained_delaunay_triangulation_2::find_minimum_sliver_triangle(double squared_tolerance) const
{
	Edge e;
	e.first = infinite_face();
	e.second = snapoly::constants::NOT_EXIST; // initialize edge e, DOESNT_EXIST indicates no sliver triangle is found

	Kernel::FT min_squared_height = snapoly::constants::DOUBLE_MAX;
	Face_handle min_triangle = infinite_face();

	for (auto& face : finite_face_handles()) {

		// find the sliver constrained triangle
		std::pair<bool, int> is_sliver = is_sliver_triangle(face, squared_tolerance); // pair<bool, int>

		// if a sliver constrained triangle is found
		if (is_sliver.first) {
			auto computed_squared_height = squared_height(face, is_sliver.second); // computed_height is positive
			if ((computed_squared_height + snapoly::constants::EPSILON) < min_squared_height) { // if min_height gets updated, we update face handle and the opposite vertex at the same time
				min_squared_height = computed_squared_height;
				e.first = face;
				e.second = is_sliver.second;
			}
		}
	} // end for: all finite faces in the triangulation

	//cout << "the squared height of the sliver triangle is: " << min_squared_height << '\n';

	return std::make_tuple(e.first, e.second, min_squared_height);
}

std::tuple<Face_handle, int, Kernel::FT> Enhanced_constrained_delaunay_triangulation_2::find_minimum_vertex_to_boundary(double squared_tolerance) const
{
	return find_minimum_sliver_triangle(squared_tolerance);
}

// mark the edge as constrained
void Enhanced_constrained_delaunay_triangulation_2::mark_constrained(Face_handle f, int i)
{
	Edge edge;
	edge.first = f;
	edge.second = i;
	if (!is_constrained(edge))
		mark_constraint(f, i); // protected member function from CDT class
}
// ====================================================================================================================================================================
