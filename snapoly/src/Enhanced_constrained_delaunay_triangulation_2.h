#ifndef ENHANCED_CONSTRAINED_DELAUNAY_TRINGULATION_2_H
#define ENHANCED_CONSTRAINED_DELAUNAY_TRINGULATION_2_H

//#include "pch.h"

template <class T>
class Enhanced_constrained_triangulation_2;

// for convenient access
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::queue;
using std::map;
using std::list;

// constants
namespace snapoly {
	namespace constants {
		const double DOUBLE_MAX = std::numeric_limits<double>::max(); // maximum value of type double
		const int NOT_EXIST = -1; // inicates an Edge object does not exist, for example edge.second = ERROR_INDEX indicating the edge does not exist
		const double EPSILON = 1e-15; // epsilon: tolerance 1e-8 by default? or use a more accurate value such as 1e-15?
	}
}


/*
* customized VertexInfo and FaceInfo - for using constrained delaunay triangulation
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
struct VertexInfo2
{};
struct FaceInfo2
{
	vector<string> faceid_collection; // can we use const char*? will there be any memory issues?
	bool processed; // used for tagging process
	//bool local_constraints_tagged; // used for adding tags to local constraints of polygons

	FaceInfo2()
		: faceid_collection(), processed(false)
	{
		faceid_collection.emplace_back("unknown"); // at least the faceid collection contains one id
	} // reserve? how many ids can a face have at most?

};

// [[maybe used]] Edge with info
struct Edge_with_info
{
	double x1;
	double y1;
	double x2;
	double y2;
	bool is_constrained;
	
	Edge_with_info() : x1(0), y1(0), x2(0), y2(0), is_constrained(false) {}
	Edge_with_info(double param_x1, double param_y1, double param_x2, double param_y2) :
		x1(param_x1), y1(param_y1), x2(param_x2), y2(param_y2), is_constrained(false) {}

	/*
	* compare if two edges are the same
	*
	* @param rhs: an rhs object which will be compared
	* @return bool: return true if same by location, otherwise false
	*/
	bool operator==(const Edge_with_info& rhs) const;	
};
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */



/*
* CGAL definitions
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Exact_predicates_tag Tag; // prerequisites for CDT
typedef CGAL::Triangulation_vertex_base_with_info_2<VertexInfo2, Kernel> Vertex_base_with_info;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> Face_base; // must use Constrained_triangulation face base class in CDT
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, Kernel, Face_base> Face_base_with_info;
typedef CGAL::Triangulation_data_structure_2<Vertex_base_with_info, Face_base_with_info> TDS; // TDS
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, Tag> CDT; // CDT
typedef CGAL::Polygon_2<Kernel> Polygon_2; // polygon_2
typedef CGAL::Polygon_with_holes_2<Kernel> Polygon_with_holes_2; // polygon with holes
typedef CGAL::Segment_2<Kernel> Segment_2; // for line segments [p,q] connecting two points p,q
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */



/*
* helpers
*/
typedef CDT::Point Point; // attention: Point stands for Kernel::Point_2, different from geom::Point
typedef CDT::Vertex Vertex;
typedef CDT::Edge Edge; // Edge: pair<Face_handle, int>
typedef CDT::Edge_circulator Edge_circulator;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Vertex_circulator Vertex_circulator;
typedef CDT::Line_face_circulator Line_face_circulator;
typedef CDT::Face_handle Face_handle;
typedef CDT::Face_circulator Face_circulator;
typedef CDT::Locate_type Locate_type; // Locate_type: enum
typedef CDT::List_faces List_faces; // list of facets
typedef CDT::List_edges List_edges; // list of edges



/*
* Enhanced constrained Delaunay triangulation class, inherited from CDT class
* for management of the constrained Delaunay triangulation
* and provide necessary functions for snap rounding.
*/
class Enhanced_constrained_delaunay_triangulation_2 : public CDT {

public:

	/*
	* calculate and return the length of an edge(f, i)
	* Principally the Face_handle can be an infinite face
	* 
	* @param f: a CDT::Face_handle
	* @param i: the opposite vertex index, Edge(f, i) represnets 
	* an edge opposite to vertex f->vertex(i) in the face f.
	* 
	* @return Kernel::FT the internal Field Type of CGAL.
	*/
	Kernel::FT squared_length(const Face_handle& f, int i) const;


	/*
	* calculate and return the length of an edge.
	* 
	* @return Kernel::FT the internal Field Type of CGAL.
	*/
	Kernel::FT squared_length(const Edge& edge) const;


	/*
	* calculate and return the unsigned area of a triangle formed by three points: p, q, r
	* 
	* @return Kernel::FT the internal Field Type of CGAL.
	*/
	Kernel::FT area(const Point& p, const Point& q, const Point& r) const;


	/*
	* calculate and return the unsigned area of a face
	* in Constrained Delaunay triangulation, a face represents for a triangle.
	* 
	* @note the face must be checked if it is an infinite face, if
	* it is an infinite face, snapoly::constant::DOUBLE_MAX will be returned
	* 
	* @return Kernel::FT the internal Field Type of CGAL.
	*/
	Kernel::FT area(const Face_handle& face) const;


	/*
	* return the squared height of a face based on a certain base edge
	* in Constrained Delaunay triangulation, a face means a triangle
	* the base edge is Edge(face, i), according to which the height is calculated.
	* 
	* @return Kernel::FT the internal Field Type of CGAL.
	*/
	Kernel::FT squared_height(const Face_handle& face, int i) const;


	/*
	* get the centroid of a face
	* pa, pb, pc are three vertices of a face
	* yields a new point
	* 
	* @return: Point
	*/
	Point face_centroid(const Point& pa, const Point& pb, const Point& pc) const;


	/*
	* get the centroid of a face
	* va, vb, vc are three vertex handles of a face
	* yields a new point
	* 
	* @return: Point
	*/
	Point face_centroid(const Vertex_handle& va, const Vertex_handle& vb, const Vertex_handle& vc) const;


	/*
	* get the centroid of a face
	* yields a new point
	* 
	* @return: Point
	*/
	Point face_centroid(const Face_handle& face) const;


	/*
	* get the centroid(midpoint) of an edge
	* yields a new point
	* 
	* @return: Point
	*/
	Point edge_centroid(const Edge& edge) const;


	/*
	* get the two vertex handles of an Edge
	* Edge(incidentFace, oppositeVertex) represents for an edge
	*/
	std::pair<Vertex_handle, Vertex_handle> vertices_of_edge(const Face_handle& incidentFace, int oppositeVertex) const;


	/*
	* get the two vertex handles of an Edge
	*/
	std::pair<Vertex_handle, Vertex_handle> vertices_of_edge(const Edge& edge) const;


	/*
	* get incident constrained vertices of a certain vertex
	* the edges between the certain vertex and its incident vertices are constrained
	* and we call these incident vertices as "constrained incident vertices"
	* 
	* this function also allows to omit a certain incident constrained vertex.
	* 
	* @param targetVh a Vertec_handle indicating a certain vertex
	* @param omitVh a Vertex_handle if it is provided, then among all constrained vertices of targetVertex, this vertex would be omitted
	* if there is an constrained edge between targetVertex and omitVertex
	* @param constrained_incident_vertices_vec where all constrained incident vertices are stored
	* 
	* @note providing an omitVertex parameter is useful when snap rounding two close vertices
	* for example, in snap rounding, we want to keep the constrained incident vertices of two close vertices respectively
	* suppose targetVertex and omitVertex are two close vertices, when we find constrained incident vertices of targetVertex
	* if the edge between targetVertex and omitVertex is constrained, then omitVertex will also be added into the vec
	* later when we re-introduce the constraints, there will be a new constraint between the midpoint and the omitVertex
	* then they become close vertices again
	* Therefore we need to omit that vertex when performing snap rounding on close vertices.
	* 
	* @return void the constrained incident vertices will be added into the constrained_incident_vertices_vec
	*/
	void get_constrained_incident_vertices(
		const Vertex_handle& targetVh, 
		vector<Point>& constrained_incident_vertices_vec, 
		Vertex_handle omitVh = Vertex_handle()) const;


	/*
	* judge if two edges are in the same location
	* va, vb are the two vertices of edge A
	* vc, vd are the two vertices of edge B
	* if only:
	* va->point() == vc->point() && vb->point() == vd->point()
	* or
	* va->point() == vd->point() && vb->point() == vc->point()
	* the two edges are the same
	* 
	* @return bool return true if two edges are in the same location otherwise false
	*/
	bool is_same_edge_by_location(const Edge& edgeA, const Edge& edgeB) const;


	/*
	* get the common edge of two faces
	* 
	* @return Edge if common edge is found return the common edge, if not return an edge with edge.second = -1 
	* indicating there is no common edge between faceA and faceB
	*/
	Edge common_edge(const Face_handle& faceA, const Face_handle& faceB) const;


	/*
	* remove an edge
	* the removing operation is done by removing two vertices of an edge from the triangulation
	*/
	void remove_edge(Edge& edge);


	/*
	* if the length of an edge(f, opposite_vertex) in the triangulation is smaller than the given tolerance
	* this means the two vertices are considered close within the tolerance
	* 
	* @param f Face_handle
	* @param opposite_vertex int
	* @param squared_tolerance double note that tolerance = sqrt(squared_tolerance)
	* 
	* @return bool return true if an edge's length is smaller than a certain tolerance
	*/
	bool is_degenerate_edge(const Face_handle& f, int opposite_vertex, double squared_tolerance) const;


	/*
	* return true if an edge's length is smaller than a certain tolerance
	* same as above
	*/
	bool is_degenerate_edge(const Edge& edge, double squared_tolerance) const;


	/*
	* find and return the degenerate edge with minimum length (<= tolerance) in each traversal
	* if not found, return an Edge object with edge.second = NOT_EXIST
	* 
	* @param squared_tolerance the squared tolerance
	* 
	* @return std::tuple<Face_handle, int, Kernel::FT>
	* Edge(Face_handle, int) and the minimum squared length will be returned
	*/
	std::tuple<Face_handle, int, Kernel::FT> find_minimum_degenerate_edge(double squared_tolerance) const;


	/*
	* for better alignment with concepts
	*/
	std::tuple<Face_handle, int, Kernel::FT> find_minimum_vertex_to_vertex(double squared_tolerance) const;


	/*
	* if an edge (face, i) is a sliver base edge / constrained sliver base edge
	* 
	* in a finite face there are 3 edges: Edge(face, 0), Edge(face, 1), Edge(face, 2)
	* the corresponding heights are: height0, height1, height2
	* for example, if height0 <= tolerance, that indicates the face is sliver (or skinny)
	* then Edge(face, 0) is called a "sliver base" - the same applies on height1 and hright2
	* if Edge(face, 0) is also constrained, then it is called a "constrained sliver base"
	* 
	* @param face a Face_handle
	* @param i opposite vertex
	* (face, i) forms an edge: the edge which is opposite to the vertex i of the face
	* @param squared_tolerance squared tolerance
	* @param constrained_flag indicates whether an edge is a constrained sliver base edge / sliver base edge
	* e.g. 
	* is_sliver_base(face, i, squared_tolerance) = true -> Edge(face, i) is a sliver base but not constrained
	* is_sliver_base(face, i, squared_tolerance, true) = true -> Edge(face, i) is a constrained sliver base
	* @return: bool return true if it is a sliver constrained base edge return false otherwise
	*
	* ATTENTION
	* Note that the projection of the opposite vertex onto the constrained base might not be within it
	* thus this type of triangle should not be snapped -
	* - currently by comparing the height and CGAL::squared_distance(Point_2, Segment_2).
	* if the projection of the opposite vertex is within the segnment
	* then std::abs(height - CGAL::squared_distance(Point_2, Segment_2)) < epsilon
	* if not, the CGAL::squared_distance(Point_2, Segment_2) would be the squared distance
	* from the Point_2 and the closest point of Segment_2, which should be one of its two ends
	* and this value will be larger than squared height
	*/
	bool is_sliver_base(const Face_handle& face, int i, double squared_tolerance, bool constrained_flag = false) const;


	/*
	* A sliver triangle is a triangle which satisifes the following conditions:
	* <1> ONLY one constrained sliver base
	*     the height of the constrained base < tolerance -> constrained sliver base
	* <2> the other two edges must not be sliver base edges
	*
	* @param face a Face_handle
	* @param squared_tolerance squared tolerance
	* 
	* @return:
	* std::pair<bool, int>
	* first indicates if it is a sliver triangle - return true if it is otherwise false
	* second indicates the opposite index of the constrained sliver base
	* i.e. if it is a sliver triangle then Edge(face, index) indicates the constrained sliver base
	*/
	std::pair<bool, int> is_sliver_triangle(const Face_handle& face, double squared_tolerance) const;


	/*
	* return the sliver triangle with the minimum height (within the tolerance) in the triangulation
	* 
	* @param squared_tolerance squared tolerance
	* 
	* @return Face_handle the minimum sliver triangle
	* @return int the opposite vertex of the constrained sliver base
	* i.e. Edge(Face_handle, int) represents the constrained sliver base
	* if such triangle is not found, this index will be -1 (snapoly::constants::NOT_EXIST)
	* @return Kernel::FT the minimum height
	*/
	std::tuple<Face_handle, int, Kernel::FT> find_minimum_sliver_triangle(double squared_tolerance) const;


	/*
	* for better alignment with concepts
	*/
	std::tuple<Face_handle, int, Kernel::FT> find_minimum_vertex_to_boundary(double squared_tolerance) const;


	/*
	* make the edge (f, i) constrained
	* mark_constraint() is protected, so it can only be accessed within the derived class
	*/
	void mark_constrained(Face_handle f, int i);


};


typedef Enhanced_constrained_delaunay_triangulation_2 Enhanced_triangulation;


namespace snapoly {
	namespace printer {

		/*
		* print out a point in the form of: (x, y)
		*/
		void print(const Vertex_handle& v);


		/*
		* print out an Edge(Face_handle, int)
		*/
		void print(const Edge& edge);


		/*
		* print out an Face_handle
		*/
		void print(const Face_handle& face);
	}
}


#endif // !ENHANCED_CONSTRAINED_DELAUNAY_TRINGULATION_2_H
