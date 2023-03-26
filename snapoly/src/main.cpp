#include "pch.h"
#include "io.h"


/*
* if there are more than one polygons
* the output would be two polygons
* when we write to gpkg, we need to check the polygons' relationship
* i.e. which is the biggest polygon -> the biggest polygon should contain the other polygons
*/




int main()
{
	//Edge e;
	Enhanced_triangulation et;
	et.insert(Kernel::Point_2(2, 1));
	et.insert(Kernel::Point_2(5, 1));
	et.insert(Kernel::Point_2(4.0, 1.2));
	et.insert_constraint(Kernel::Point_2(2, 1), Kernel::Point_2(5, 1));

	//cout << Kernel::Point_2(2, 1) << '\n';
	

	double tolerance = 0.2;
	double squared_tolerance = tolerance * tolerance;

	for (auto const& face : et.finite_face_handles()) {
		for (int i = 0; i < 3; ++i) {
			if (et.is_sliver_base(face, i, squared_tolerance, true))
				cout << face->vertex(i)->point() << '\n';
		}
		cout << et.is_sliver_triangle(face, squared_tolerance).first << '\n';
		snapoly::printer::print(face);
	}

	auto find = et.find_minimum_sliver_triangle(squared_tolerance);
	cout << std::get<2>(find) << '\n';


	// -----------------------------------------------------------------------------------------

	Snap_rounding_2 sr;
	const char* input_file = R"(D:\snapoly\data\ab.gpkg)"; // Andorra_buildings_1
	snapoly::io::add_polygons_from_input_file(input_file, sr.polygons());

	sr.insert_polygons_to_triangulation();

	sr.add_tag_to_triangulation();

	sr.set_tolerance(1.50659e-05); // 9.50158e-06

	sr.snap_rounding();
	
	const char* output_file = R"(D:\snapoly\data\ab_res.gpkg)"; // Andorra_buildings_1
	snapoly::io::export_to_gpkg(output_file, sr.constraintsWithInfo());

	cout << "file saved at: " << output_file << '\n';
	
	// -----------------------------------------------------------------------------------------
	
	const char* res_file = R"(D:\snapoly\data\abres.gpkg)"; // Andorra_buildings_1

	vector<CDTPolygon> resPolygonsVec;
	snapoly::io::build_polygons_from_constraints(sr.constraintsWithInfo(), resPolygonsVec);

	snapoly::io::export_to_gpkg(res_file, resPolygonsVec);


	Kernel::Point_2 p0(1, 1), p1(1, 3);
	cout << CGAL::midpoint(p0, p1);

	//Point p;

	Polygon_2 pgn;
	pgn.push_back(CDTPoint(0, 0));
	pgn.push_back(CDTPoint(2, 0));
	pgn.push_back(CDTPoint(1, 1));
	snapoly::printer::print(pgn);
	


	return 0;
}