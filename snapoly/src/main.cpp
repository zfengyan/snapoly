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

	Snap_rounding_2 sr;
	const char* input_file = R"(D:\snapoly\data\ab.gpkg)"; // Andorra_buildings_1
	snapoly::io::add_polygons_from_input_file(input_file, sr.polygons());

	sr.insert_polygons_to_triangulation();

	sr.add_tag_to_triangulation();

	sr.set_tolerance(1.50659e-05); // 9.50158e-06

	sr.snap_rounding();

	const char* output_file = R"(D:\snapoly\data\abcon.gpkg)"; // Andorra_buildings_1
	snapoly::io::export_to_gpkg(output_file, sr.constraintsWithInfo());

	cout << "file saved at: " << output_file << '\n';

	// -----------------------------------------------------------------------------------------

	const char* res_file = R"(D:\snapoly\data\abres.gpkg)"; // Andorra_buildings_1

	snapoly::io::build_polygons_from_constraints(sr.constraintsWithInfo(), sr.result_polygons());

	snapoly::io::export_to_gpkg(res_file, sr.result_polygons());

	const char* tri_file = R"(D:\snapoly\data\abtri.gpkg)"; // Andorra_buildings_1
	snapoly::io::export_to_gpkg(tri_file, sr.triangulation());

	// -----------------------------------------------------------------------------------------

	cout << sr.measure_distortions();

	return 0;
}