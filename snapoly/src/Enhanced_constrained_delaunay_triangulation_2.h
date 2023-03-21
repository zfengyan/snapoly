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
const double epsilon = 1e-15; // epsilon: tolerance 1e-8 by default? or use a more accurate value such as 1e-15?
const double DOUBLE_MAX = std::numeric_limits<double>::max(); // maximum value of type double
const int DOESNT_EXIST = -1; // inicates an Edge object does not exist, for example edge.second = ERROR_INDEX indicating the edge does not exist


/*
* customized VertexInfo and FaceInfo - for using constrained delaunay triangulation
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
struct VertexInfo2
{
	VertexInfo2() { someAttri = 1; }
	int someAttri;
	void print() const {
		std::cout << "someAttri: " << someAttri << std::endl;
	}
};
struct FaceInfo2
{
	int nesting_level;
	vector<string> faceid_collection; // can we use const char*? will there be any memory issues?
	bool processed; // used for tagging process
	//bool local_constraints_tagged; // used for adding tags to local constraints of polygons

	FaceInfo2()
		: nesting_level(-1), faceid_collection(), processed(false)
	{
		faceid_collection.emplace_back("unknown"); // at least the faceid collection contains one id
	} // reserve? how many ids can a face have at most?

	bool in_domain() const {
		return nesting_level % 2 == 1; // if nesting_level is odd numbers
	}
	void print() const {
		std::cout << "nesting_level " << nesting_level << std::endl;
	}

};



// [[maybe used]] Edge with info
struct Edge_with_info
{
	Edge_with_info() : x1(0), y1(0), x2(0), y2(0), is_constrained(false) {}
	Edge_with_info(double param_x1, double param_y1, double param_x2, double param_y2) :
		x1(param_x1), y1(param_y1), x2(param_x2), y2(param_y2), is_constrained(false) {}

	bool operator==(const Edge_with_info& rhs) const {
		return (
			((std::abs(this->x1 - rhs.x1) < epsilon) &&
				(std::abs(this->y1 - rhs.y1) < epsilon) &&
				(std::abs(this->x2 - rhs.x2) < epsilon) &&
				(std::abs(this->y2 - rhs.y2) < epsilon)) ||
			((std::abs(this->x1 - rhs.x2) < epsilon) &&
				(std::abs(this->y1 - rhs.y2) < epsilon) &&
				(std::abs(this->x2 - rhs.x1) < epsilon) &&
				(std::abs(this->y2 - rhs.y1) < epsilon))
			);
	} // overloading == operator, for find() function with std::vector<Edge_with_info>

	double x1;
	double y1;
	double x2;
	double y2;
	bool is_constrained;
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

#endif // !ENHANCED_CONSTRAINED_DELAUNAY_TRINGULATION_2_H
