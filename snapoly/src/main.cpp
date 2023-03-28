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
	// files
	const char* input_file = R"(D:\snapoly\data\netherlands\Denhaag.gpkg)";
	const char* tri_file = R"(D:\snapoly\data\netherlands\Denhaag_tri.gpkg)";
	const char* output_boundaries_file = R"(D:\snapoly\data\netherlands\Denhaag_boundaries.gpkg)";
	const char* res_file = R"(D:\snapoly\data\netherlands\Denhaag_res.gpkg)";

	// Snap rounding
	Snap_rounding_2 sr;

	io::add_polygons_from_input_file(input_file, sr.polygons());

	sr.insert_polygons_to_triangulation();

	sr.add_tag_to_triangulation();

	//sr.set_tolerance(1.50659e-05); // 9.50158e-06

	//sr.snap_rounding();

	//TODO:
	//1. coordinates shifting
	//2. when exporting to the GPKG file, attach the CRS information

	io::export_to_gpkg(tri_file, sr.triangulation());

	io::export_to_gpkg(output_boundaries_file, sr.constraintsWithInfo());

	//cout << "file saved at: " << output_file << '\n';

	// -----------------------------------------------------------------------------------------

	io::build_polygons_from_constraints(sr.constraintsWithInfo(), sr.result_polygons());

	io::export_to_gpkg(res_file, sr.result_polygons());




	// -----------------------------------------------------------------------------------------

	//cout << sr.measure_distortions();

	return 0;
}