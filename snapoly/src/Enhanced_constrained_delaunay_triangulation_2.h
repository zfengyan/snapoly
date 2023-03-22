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
	namespace constant {
		const double DOUBLE_MAX = std::numeric_limits<double>::max(); // maximum value of type double
		const int DOESNT_EXIST = -1; // inicates an Edge object does not exist, for example edge.second = ERROR_INDEX indicating the edge does not exist
		const double epsilon = 1e-15; // epsilon: tolerance 1e-8 by default? or use a more accurate value such as 1e-15?
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
// typedef CDT::Point Point; // Kernel::Point_2 would be conflict with geom::Point, thus we explicitly use Kernel::Point_2
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



};


#endif // !ENHANCED_CONSTRAINED_DELAUNAY_TRINGULATION_2_H
