#include "pch.h"
#include "Snap_rounding_2.h"

using namespace snapoly::constants;

// merge the id collections
// the idCollection of rhs will be merged into the idCollection of *this
// the idCollection of rhs will not be changed
// e.g.
// idCollection of *this: [A, B], idCollection of rhs: [B, C, D]
// after the merge operation
// idCollection of *this: [A, B, C, D]
void Constraint::merge_id_collection(const vector<string>& rhs_idCollection)
{
	// unordered set for the ids of *this
	unordered_set<string> idSet;

	// populate the idSet, including the default value "default"
	for (auto& id : idCollection) { // this->idCollection
		idSet.insert(id);
	}

	// merge the ids from rhs.idCollection to idCollection (this->idCollection)
	// if the id from rhs.idCollection is not present in idSet, add it
	for (auto& rhs_id : rhs_idCollection) {
		if (idSet.find(rhs_id) == idSet.end()) { // rhs_id is not in the idCollection
			idCollection.push_back(rhs_id);
		}
	}
}


// is a point inside a CDTPolygon
bool CDTPolygon::is_point_inside_polygon(const CDTPoint& pt, const CDTPolygon& pgn)
{
	// Knernel that will be used in CGAL::bounded_side_2 function
	Kernel kernel{};

	// variables for predicate
	CGAL::Bounded_side outerSide = CGAL::ON_BOUNDED_SIDE; // initialized as inside
	CGAL::Bounded_side innerSide = CGAL::ON_UNBOUNDED_SIDE; // initialized as outside

	// check if point lies inside the exterior (outer boundary)
	outerSide = CGAL::bounded_side_2(pgn.outer_boundary().vertices_begin(), pgn.outer_boundary().vertices_end(), pt, kernel);

	// check if point lies outside each interior (ALL of the interiors)
	for (auto const& hole : pgn.holes()) {
		innerSide = CGAL::bounded_side_2(hole.vertices_begin(), hole.vertices_end(), pt, kernel);
		if (innerSide != CGAL::ON_UNBOUNDED_SIDE)return false;
	}

	// return true only if two conditions are satisfied
	return ((outerSide == CGAL::ON_BOUNDED_SIDE) && (innerSide == CGAL::ON_UNBOUNDED_SIDE));
}


// calculate the area of a CDTPolygon
double CDTPolygon::area() const
{
	const Polygon_2& exterior = m_outerRing;
	double exteriorArea = std::abs(CGAL::polygon_area_2(exterior.vertices_begin(), exterior.vertices_end(), Kernel()));
	double interiorArea = 0;
	if (has_holes()) {
		for (auto const& hole : holes()) {
			interiorArea += std::abs(CGAL::polygon_area_2(hole.vertices_begin(), hole.vertices_end(), Kernel()));
		}
	}
	return std::abs(exteriorArea - interiorArea);
}


// set the tolerance
void Snap_rounding_2::set_tolerance(double tolerance_param)
{
	m_tolerance = tolerance_param;
	m_squared_tolerance = m_tolerance * m_tolerance;
	cout << "the tolerance is set to: " << m_tolerance << '\n';
}


// consider the boundaries of polygons as constraints and insert them to the triangulation
void Snap_rounding_2::insert_polygons_to_triangulation()
{
	cout << "inserting polygons to triangulation ... \n";

	Timer timer; // for logging the run time

	for (auto const& pgn : m_polygons)
	{
		// insert exterior as constraints
		for (auto iter = pgn.outer_boundary().edges_begin(); iter != pgn.outer_boundary().edges_end(); ++iter)
		{
			m_et.insert_constraint(iter->source(), iter->target());
		}

		// if not containing hole(s), continue
		//if (!pgn.has_holes())continue;

		// if containing holes, insert as constraints
		for (auto const& hole : pgn.holes())
		{
			for (auto iter = hole.edges_begin(); iter != hole.edges_end(); ++iter)
			{
				m_et.insert_constraint(iter->source(), iter->target());
			}
		}

	}
	cout << "done \n";
}


// add tag to one polygon
void Snap_rounding_2::add_tag_to_one_polygon(Face_handle& startingFace, const CDTPolygon& refPgn)
{
	if (m_et.is_infinite(startingFace)) {
		std::cerr << "currently the seeding face can not be an infinite face, please check: \n"
			<< "add_tag_to_one_polygon() function in Snap_rounding.hpp file" << endl;
		return;
	}

	queue<Face_handle> faces; // for BFS

	// enqueue the seedingFace
	faces.push(startingFace);

	while (!faces.empty())
	{
		// process the front element
		// when adding elements to the queue, some faces may be added more than once
		// in this sense we need to judge if it has been processed yet before we process it
		Face_handle currentFace = faces.front();

		if (!currentFace->info().processed) { // if not yet been processed, process it

			currentFace->info().faceid_collection.push_back(refPgn.id()); // a face may be tagged more than once

			// === tag the constrained edges of this face (if any) ===
			for (int currentVertex = 0; currentVertex < 3; ++currentVertex) { // find the possible constrained edges of the current face
				Edge e;
				e.first = currentFace;
				e.second = currentVertex;
				if (m_et.is_constrained(e)) { // if it is a constrained edge
					auto vertex_pair = m_et.vertices_of_edge(e); // get two vertices of the constrained edge
					Vertex_handle va = vertex_pair.first;
					Vertex_handle vb = vertex_pair.second;
					Constraint c(va->point(), vb->point());

					// check if the current constraint is present in the list
					// if so we merge the currrent id to it 
					// instead of directly adding - avoid repeatness, e.g. small overlapping area
					auto it = std::find(m_constraintsWithInfo.begin(), m_constraintsWithInfo.end(), c);
					if (it != m_constraintsWithInfo.end()) {
						it->merge_id_collection(currentFace->info().faceid_collection);
						//it->idCollection.push_back(currentFace->info().faceid_collection[1]);
					}
					else { // if the current constraint is not present yet
						// attach the tag info and add it to the constraintsWithID collection
						c.idCollection.push_back(currentFace->info().faceid_collection[1]);
						m_constraintsWithInfo.push_back(c);
					}			
				}
			} // end for: all three edges of the current face
			// === tag the constrained edges of this face (if any) ===

			currentFace->info().processed = true;
		}

		faces.pop(); // if already processed, pop it

		// Add all possible finite neighbours - not crossing constrained edges and not yet processed
		for (int i = 0; i < 3; ++i) {
			Face_handle neighborFace = currentFace->neighbor(i);
			if (!m_et.is_infinite(neighborFace)) { // only add finite neighbors
				Edge commonEdge = m_et.common_edge(currentFace, neighborFace);
				if ((!m_et.is_constrained(commonEdge)) && (!neighborFace->info().processed)) {
					faces.push(neighborFace);
				}
			}
		} // end for: all neighbors (including infinite neighbors)

	} // end while: while the queue is not empty
}


// add tags to the triangulation
void Snap_rounding_2::add_tag_to_triangulation()
{
	cout << "adding tags to triangulation ... \n";

	Timer timer; // for logging the run time
	
	// if empty vector
	if (m_polygons.size() == 0) {
		std::cerr << "no polygons found, please check! \n";
		return;
	}

	// for each polygon, polygon itself should not be changed during the tagging
	// the overlapping case is not handled here
	int count = 0; // Debug
	for (auto const& pgn : m_polygons)
	{
		std::size_t numExterior = static_cast<std::size_t>(pgn.number_of_exterior_points());
		vector<Face_handle> starting_faces;
		starting_faces.reserve(numExterior + 1); // each face handle corresponds to a vertex if possible

		for (auto constIter = pgn.exterior_begin(); constIter != pgn.exterior_end(); ++constIter) {
			CDTPoint p(constIter->x(), constIter->y());
			Face_circulator fc = m_et.incident_faces(m_et.insert(p)), fc_done(fc); // here insert might introduce new point which is removed by remove_overlap
			if (fc != 0) {
				do {
					if (!m_et.is_infinite(fc))
					{
						CDTPoint centroid = CGAL::centroid(fc->vertex(0)->point(), fc->vertex(1)->point(), fc->vertex(2)->point());
						bool centroid_inside = CDTPolygon::is_point_inside_polygon(centroid, pgn); // check inside
						if (centroid_inside) {
							starting_faces.push_back(static_cast<Face_handle>(fc));
							break; // comment this will allows to handle the overlapping area
							// since a starting face can have all 3 edges constrained
							// yet it still depends on the specific arrangements (see overlapping example)
						}
					}

				} while (++fc != fc_done);
			}
		}// cout << "starting faces size: " << starting_faces.size() << '\n';		

		// add tag for each starting face
		if (starting_faces.size()) {
			for (auto& start_face : starting_faces) {
				add_tag_to_one_polygon(start_face, pgn);
			}
		}

		// Debug
		++count;
		if(count % 1000 == 0)cout << count << '\n';
		// Debug

	} // end for: polygons

	cout << "done \n";

}


// remove dangles in the constraintsWithInfo list
//   a   b           centroid
//   /\  /\           / | \ 
//  /  \/  \    ->   /  |  \
//  |      |         |     |
// after snap rounding a and b, there will be redundant constraint dangles in the middle (connecting to the centroid)
// <1> remove dangles from the constraintsWithInfo list
// <2> remove dangles from the triangulation
void Snap_rounding_2::remove_dangles()
{
	cout << "removing dangles ... \n";

	vector<Constraint> danglingConstraintsVec; // store the dangling constraints for removing them from the constraints list
	danglingConstraintsVec.reserve(m_constraintsWithInfo.size());

	vector<CDTPoint> danglingVerticesVec; // store the dangling vertices for removing them from the triangulation

	list<Constraint> constraintsRef = m_constraintsWithInfo; // for references

	int count = 0; // Debug
	for (auto const& c : m_constraintsWithInfo) {

		//Debug
		//++count;
		//cout << count << '\n';
		//Debug


		// if constraint is a dangle:
		// for the ending points of constraint: ea and eb
		// they must be the ending points of other constraints with different ids
		// if not then this constraint(s) is a dangle
		int count_endpoint_p0_present_in_other_constraints = 0;
		int count_endpoint_p1_present_in_other_constraints = 0;

		for (auto const& cRef : constraintsRef) {
			if (c != cRef) { // not itself - just by location, there can be constraints with different ids at the same location
				if (c.p0 == cRef.p0 || c.p0 == cRef.p1)
					++count_endpoint_p0_present_in_other_constraints;
				if (c.p1 == cRef.p0 || c.p1 == cRef.p1)
					++count_endpoint_p1_present_in_other_constraints;
			}
		} // end for: all reference constraints

		// if c is a dangling constraint and p0 is the dangling point
		if (!count_endpoint_p0_present_in_other_constraints) {
			danglingConstraintsVec.push_back(c);
			danglingVerticesVec.push_back(c.p0);
		}

		// if c is a dangling constraint and p1 is the dangling point
		if (!count_endpoint_p1_present_in_other_constraints) {
			danglingConstraintsVec.push_back(c);
			danglingVerticesVec.push_back(c.p1);
		}
	}

	// remove dangles 1 if any: remove dangling constraints from the constraints list
	if (danglingConstraintsVec.size()) {
		for (auto const& dangle : danglingConstraintsVec)
			m_constraintsWithInfo.remove(dangle);
	}
	
	// remove dangles 2 if any: remove dangling constraints from the triangulation
	if (danglingVerticesVec.size()) {
		for (auto const& dangle : danglingVerticesVec) {
			Vertex_handle v = m_et.insert(dangle); // get the vertex handle of the point
			m_et.remove(v); // remove the dangling vertex
			//cout << "dangling vertex: " << v->point() << '\n';
		}
	}

	cout << "done \n";

}


// snap close vertices
void Snap_rounding_2::snap_vertex_to_vertex(Edge& edgeOfVertexToVertex)
{
	// process the geometry - update the constraintsWithID first and then alter the triangulation
	auto vertex_pair = m_et.vertices_of_edge(edgeOfVertexToVertex);
	Vertex_handle va = vertex_pair.first;
	Vertex_handle vb = vertex_pair.second;

	std::vector<CDTPoint> constrained_incident_vertices_va; // incident constrainted points of va
	std::vector<CDTPoint> constrained_incident_vertices_vb; // incident constrainted points of vb

	CDTPoint centroid = m_et.edge_centroid(edgeOfVertexToVertex); // get the centroid of the constrained degenerate edge

	// check if va - vb is constrained, if constrained we need to remove it from the constrintsWithInfo
	if (m_et.is_constrained(edgeOfVertexToVertex)) {
		Constraint c(va->point(), vb->point());
		m_constraintsWithInfo.remove(c);
	}

	// ATTENTION: Will re-introduction of constraints intersect with other constraints?

	// for va and constraints of va, update the corresponding constraint in constraintsWithInfo 
	// if va-vb is constrained, vb is not one of the incident points of constraints of va
	m_et.get_constrained_incident_vertices(va, constrained_incident_vertices_va, vb);
	for (auto& incidentPoint : constrained_incident_vertices_va) {
		Constraint c(va->point(), incidentPoint);
		for (auto& refC : m_constraintsWithInfo) { // find the constraint in the vector
			if (c == refC) { // change the point to the centroid
				if (refC.p0 == va->point())
					refC.p0 = centroid;
				else 
					refC.p1 = centroid;
			}
		}
	} // end for: all incident points of va 

	//for vb and constraints of vb, update the corresponding constraint in constraintsWithInfo 
	// if va-vb is constrained, va is not one of the incident points of constraints of vb
	m_et.get_constrained_incident_vertices(vb, constrained_incident_vertices_vb, va);
	for (auto& incidentPoint : constrained_incident_vertices_vb) {
		Constraint c(vb->point(), incidentPoint);
		for (auto& refC : m_constraintsWithInfo) { // find the constraint in the vector
			if (c == refC) { // change the point to the centroid
				if (refC.p0 == vb->point())
					refC.p0 = centroid;
				else 
					refC.p1 = centroid;
			}
		}
	} // end for: all incident points of vb

	// modification of the original triangulation is a must
	// otherwise the close vertices will always be found

	//test
	/*cout << "va: " << va->point() << '\n';
	cout << "incident points of constraints: " << '\n';
	for (auto const& p : incident_points_of_constraints_va)
		cout << p << '\n';
	cout << '\n';
	cout << "vb: " << vb->point() << '\n';
	cout << "incident points of constraints: " << '\n';
	for (auto const& p : incident_points_of_constraints_vb)
		cout << p << '\n';*/
	//test

	// after updating the constraintsWithInfo, alter the triangulation
	m_et.remove_edge(edgeOfVertexToVertex); //remove the edge - remove the vertices of the constrained degenerate edge

	// re-introduce the constraints with the centroid and other incident vertices of constraints

	// re-introduce the constraints incident to va
	if (!constrained_incident_vertices_va.empty())
	{
		for (auto const& p : constrained_incident_vertices_va)
			m_et.insert_constraint(centroid, p);
	}

	// re-introduce the constraints incident to vb
	if (!constrained_incident_vertices_vb.empty()) {
		for (auto const& p : constrained_incident_vertices_vb)
			m_et.insert_constraint(centroid, p);
	}

	// if no constrained points found, insert the centroid
	if (constrained_incident_vertices_va.empty() || constrained_incident_vertices_vb.empty())
		m_et.insert(centroid);

	// check redundacy and dangles:
	//   a   b           centroid
	//   /\  /\           / | \ 
	//  /  \/  \    ->   /  |  \
	//  |      |         |     |
	// after snap rounding a and b, there will be redundant constraint dangles in the middle (connecting to the centroid)
	// and there will also be dangling vertex

	// remove dangles - the dangles / redundant dangles will be removed
	// remove operation includes remove the dangling constraints from the constraints list
	// and remove the dangling vertices from the triangulation - yet this is a very expensive function especially to the large datasets
	//remove_dangles();

	//cout << "processing geometry done\n"; cout << '\n';

}


// snap close vertex to boundary
void Snap_rounding_2::snap_vertex_to_boundary(Face_handle& sliverFace, int capturingVertexIndex)
{
	// we should also avoided altering the triangulation while traversing it

	// get the capturingVertex
	// (sliverFace, oppositeVertexIndex) represents for the constrained base edge
	Vertex_handle capturingVertex = sliverFace->vertex(capturingVertexIndex);

	auto vertexPairOfConstrainedBase = m_et.vertices_of_edge(sliverFace, capturingVertexIndex);
	Vertex_handle va = vertexPairOfConstrainedBase.first;
	Vertex_handle vb = vertexPairOfConstrainedBase.second;

	// update constraintsWithInfo first

	// constraint of the constrained base
	Constraint c(va->point(), vb->point());

	// constraints of the other two edges, note they can be constrained and already have id(s)
	Constraint ca(va->point(), capturingVertex->point());
	Constraint cb(vb->point(), capturingVertex->point());

	// first we find the constraint "c" in the constraints with info list
	auto it = std::find(m_constraintsWithInfo.begin(), m_constraintsWithInfo.end(), c);
	if (it != m_constraintsWithInfo.end()) { // "c" is found and returned the pointer pointing to it
		c.idCollection = it->idCollection; // populate the idCollection of c
	}
	

	//Debug
	//cout << "id: " << c.idCollection[1] << '\n';
	//Debug

	// then we use the same way to find the constraint "ca" and "cb"
	// "ca" and "cb" may / may not be present in the constraints with info list
	// if they are present, we merge the idCollection of c into their idCollections


	// find "ca"
	bool find_ca = false;
	auto ita = std::find(m_constraintsWithInfo.begin(), m_constraintsWithInfo.end(), ca);
	if (ita != m_constraintsWithInfo.end()) { // if "ca" is found, merge the idCollection of c
		ita->merge_id_collection(c.idCollection);
		find_ca = true;
	} // if not found, the idCollection of ca only contains "unknown"

	// find "cb"
	bool find_cb = false;
	auto itb = std::find(m_constraintsWithInfo.begin(), m_constraintsWithInfo.end(), cb);
	if (itb != m_constraintsWithInfo.end()) { // if "cb" is found, merge the idCollection of c
		itb->merge_id_collection(c.idCollection);
		find_cb = true;
	} // if not found, the idCollection of cb only contains "unknown"


	// if ca / cb is not found - that means ca / cb is not constrained
	// and it doesn't contain any id information
	// we assign the idCollection of c, and we add it to the constraints with info list


	if (!find_ca) {
		ca.idCollection = c.idCollection;
		m_constraintsWithInfo.push_back(ca);
	}

	if (!find_cb) {
		cb.idCollection = c.idCollection;
		m_constraintsWithInfo.push_back(cb);
	}

	// remove the constraint "c" from the constraintsWithInfo
	// the remove function of std::list is to use the argument's value to find and remove the element in the list
	m_constraintsWithInfo.remove(c);


	// after update of the constraiants with info list
	// we modify the geometry of the triangulation


	//cout << "processing geometry ... \n";
	int ccw = sliverFace->ccw(capturingVertexIndex);
	int cw = sliverFace->cw(capturingVertexIndex);
	m_et.mark_constrained(sliverFace, ccw);
	m_et.mark_constrained(sliverFace, cw); // ensure the other edges constrained (if not make constrained)
	// remove the constrained of the constrained base, just the constraint status is changed, the local structure does not change
	m_et.remove_constrained_edge(sliverFace, capturingVertexIndex);
	//cout << "processing geometry done\n";

	// at last step we remove possible dangles from the constraints with info list and the triangulation
	// this step is important, and to avoid the possible cascading effects
	// during the snap rounding - yet this is a very expensive function especially to the large datasets
	//remove_dangles();
}


// snap rounding
void Snap_rounding_2::snap_rounding()
{
	cout << '\n';
	cout << "snap rounding... \n";

	Timer timer; // for logging the run time

	//Debug
	int count = 0;
	//Debug

	//-------------------------------------------------------------------------------------------------------------------------
	while (true) {

		// find minimum vertex to vertex
		//cout << "finding minimum vertex to vertex ... \n";
		std::tuple<Face_handle, int, Kernel::FT>
			findMinimumVertexToVertex = m_et.find_minimum_vertex_to_vertex(m_squared_tolerance);
		Face_handle incidentFaceVertexToVertex = std::get<0>(findMinimumVertexToVertex);
		int indexVertexToVertex = std::get<1>(findMinimumVertexToVertex);
		Kernel::FT minSquaredDistVertexToVertex = std::get<2>(findMinimumVertexToVertex);
		Edge edgeOfVertexToVertex(incidentFaceVertexToVertex, indexVertexToVertex); // edge connecting two close vertices
		//cout << "done \n";

		// find minimum vertex to boundary
		//cout << "finding minimum vertex to boundary ... \n";
		std::tuple<Face_handle, int, Kernel::FT>
			findMinimumVertexToBoundary = m_et.find_minimum_vertex_to_boundary(m_squared_tolerance);
		Face_handle sliverFace = std::get<0>(findMinimumVertexToBoundary);
		int capturingVertexIndex = std::get<1>(findMinimumVertexToBoundary);
		Kernel::FT minSquaredDistVertexToBoundary = std::get<2>(findMinimumVertexToBoundary);
		//cout << "done \n";

		//cout << indexVertexToVertex << " " << capturingVertexIndex << '\n';
		//cout << "minimum distance of Vertex to Vertex: " << std::sqrt(minSquaredDistVertexToVertex) << '\n';
		//cout << "minimum distance of Vertex to Boundary: " << std::sqrt(minSquaredDistVertexToBoundary) << '\n';

		// exit condition
		// if close vertex to vertex and vertex to boundary are not found
		if (indexVertexToVertex == snapoly::constants::NOT_EXIST && capturingVertexIndex == snapoly::constants::NOT_EXIST) {
			cout << "no more snap rounding cases found \n";
			break;
		}		
		else if (indexVertexToVertex == snapoly::constants::NOT_EXIST && capturingVertexIndex != snapoly::constants::NOT_EXIST) {
			// call snap vertex to boundary 
			cout << "no more close vertex - vertex is found under given tolerance, snap vertex to boundary\n";
			snap_vertex_to_boundary(sliverFace, capturingVertexIndex); cout << '\n';
		}
		else if (indexVertexToVertex != snapoly::constants::NOT_EXIST && capturingVertexIndex == snapoly::constants::NOT_EXIST) {
			// call snap vertex to vertex
			cout << "no more close vertex - boundary is found under given tolerance, snap vertex to vertex\n";
			snap_vertex_to_vertex(edgeOfVertexToVertex); cout << '\n';
		}
		else {

			// compare the squared distance
			// if vertex to vertex is closer than vertex to boundary, snap vertex to vertex first
			// if not (> or ==), snap vertex to boundary (when ==, snap vertex to boundary is a more robust way)
			if ((minSquaredDistVertexToVertex + snapoly::constants::EPSILON) < minSquaredDistVertexToBoundary) {
				cout << "snap vertex to vertex first \n"; cout << '\n';
				snap_vertex_to_vertex(edgeOfVertexToVertex);	
			}
			else {
				cout << "snap vertex to boundary first \n"; cout << '\n';
				snap_vertex_to_boundary(sliverFace, capturingVertexIndex);	
			}

			//cout << "done, next loop \n";

		}

		//Debug
		//++count;
		//cout << "sr: " << count << '\n';
		//if (count == 1)break;
		//Debug

	} // end while: until no cases are found under given tolerance

	//-------------------------------------------------------------------------------------------------------------------------

	cout << "done \n";
	cout << '\n';
}


// measure the distortions by area
double Snap_rounding_2::measure_distortions() const
{
	// first make a polygon map
	unordered_map<string, CDTPolygon> originalPolygons;

	// populate the map
	for (auto const& poly : m_polygons)
		originalPolygons[poly.id()] = poly;

	// calculate the area difference
	double area_diff_sum = 0;
	for (auto const& resPoly : m_result_polygons) {
		const string& id = resPoly.id();
		area_diff_sum += std::abs(originalPolygons[id].area() - resPoly.area());
	}

	// get the area sum of the input polygons
	double area_original_sum = 0;
	for (auto const& originalPoly : m_polygons) {
		area_original_sum += originalPoly.area();
	}

	// get the preportion
	double area_diff_percentage =
		area_original_sum > snapoly::constants::EPSILON ?
		(area_diff_sum / area_original_sum) : 0;

	cout << "area_diff: " << area_diff_percentage << '\n';	

	return area_diff_percentage;
}


// check the minimum distance between <1> point to point <2> point to boundary
double Snap_rounding_2::minimum_distance() const
{
	// minimum distance between point to point
	auto minimum_vertex_to_vertex = m_et.find_minimum_vertex_to_vertex(m_squared_tolerance); // the squared_distance is returned
	double minimum_squared_dist_vertex_to_vertex = std::get<2>(minimum_vertex_to_vertex);

	// minimum distance between point to boundary
	auto minimum_vertex_to_boundary = m_et.find_minimum_vertex_to_boundary(m_squared_tolerance); // the squared_distance is returned
	double minimum_squared_dist_vertex_to_boundary = std::get<2>(minimum_vertex_to_boundary);

	// check if no minimum is found
	if (std::abs(minimum_squared_dist_vertex_to_vertex - DOUBLE_MAX) < EPSILON ||
		std::abs(minimum_squared_dist_vertex_to_boundary - DOUBLE_MAX) < EPSILON) {
		cout << "no minimum distance is found under the given tolerance: " << m_tolerance << '\n';
		cout << "0 will be returned \n";
		cout << "squared tolerance: " << m_squared_tolerance << '\n';
		return 0;
	}

	// get the minimum squared distance
	double minimum_squared_dist = std::min(minimum_squared_dist_vertex_to_vertex, minimum_squared_dist_vertex_to_boundary);

	return std::sqrt(minimum_squared_dist);
}


// find the distances between point to point / boundary and store them in a priority queue
void Snap_rounding_2::find_tolerance(std::priority_queue<double>& lengthQueue)
{
	// vertex and vertex are connected by edges
	for (auto& edge : m_et.finite_edges()) {
		auto computed_squared_length = m_et.squared_length(edge);
		if(computed_squared_length + EPSILON < m_squared_tolerance)
			lengthQueue.push(std::sqrt(computed_squared_length));
	}

	// vertex to boundary - can be a vertex to another polygon's boundary or to itself's boundary
	// the distance of vertex to boundary is defined as the distance between a capturing vertex and 
	// the constrained sliver base in a sliver triangle
	// a sliver triangle is a triangle ONLY has one constrained sliver base
	// the other two edges must not be sliver edges (whether they are constrained or not)
	for (auto& face : m_et.finite_face_handles()) {

		// find the sliver constrained triangle
		std::pair<bool, int> is_sliver = m_et.is_sliver_triangle(face, m_squared_tolerance); // pair<bool, int>

		// if a sliver constrained triangle is found
		if (is_sliver.first) {
			auto computed_squared_height = m_et.squared_height(face, is_sliver.second); // computed_height is positive
			if(computed_squared_height + EPSILON < m_squared_tolerance)
				lengthQueue.push(std::sqrt(computed_squared_height));
		}
	} // end for: all finite faces in the triangulation
}


// print a Polygon_2
void snapoly::printer::print(const Polygon_2& polygon_2)
{
	cout << "=== CGAL::Polygon_2 (Ring) ===\n";
	for (auto iter = polygon_2.vertices_begin(); iter != polygon_2.vertices_end(); ++iter)
	{
		cout << "(" << iter->x() << ", " << iter->y() << ")" << '\n';
	}
	cout << "======\n";
}


