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
	
	const char* input_file = R"(D:\snapoly\data\examples\triangle.gpkg)";
	const char* tri_file = R"(D:\snapoly\data\examples\triangle_tri.gpkg)";
	const char* output_boundaries_file = R"(D:\snapoly\data\examples\triangle_boundaries.gpkg)";
	const char* res_file = R"(D:\snapoly\data\examples\triangle_res.gpkg)";

	// Timer
	Timer timer;

	// Snap rounding
	Snap_rounding_2 sr;
	sr.set_tolerance(0.01); // for Delft: 0.01m

	io::add_polygons_from_input_file(input_file, sr.polygons());

	sr.insert_polygons_to_triangulation();

	sr.add_tag_to_triangulation();

	// find the minimum tolerance ---------------------------------------------
	std::priority_queue<double> lengthQueue;
	sr.find_tolerance(lengthQueue);
	cout << "snap rounding cases uner the current tolerance: " << lengthQueue.size() << '\n';
	cout << "distances under the current tolerance: " << '\n';
	while (!lengthQueue.empty()) {
		//cout << lengthQueue.top() << '\n';
		lengthQueue.pop();
	}
	//cout << "the squared height of the sliver triangle is: " << min_squared_height << '\n';
	
	// find the minimum tolerance ---------------------------------------------

	//sr.snap_rounding(); // not using remove dangles()



	//double minimum_distance = sr.minimum_distance();
	//cout << "tolerance: " << sr.tolerance() << '\n';
	//cout << " minimum distance under the current tolerance: " << minimum_distance << '\n';

	io::export_to_gpkg(tri_file, sr.triangulation());

	io::export_to_gpkg(output_boundaries_file, sr.constraintsWithInfo());

	//cout << "file saved at: " << output_file << '\n';

	// -----------------------------------------------------------------------------------------

	io::build_polygons_from_constraints(sr.constraintsWithInfo(), sr.result_polygons());

	io::export_to_gpkg(res_file, sr.result_polygons());

	sr.measure_distortions(); // this function must be called after the io::build_polygons_from_constraints() function


	// -----------------------------------------------------------------------------------------



	return 0;


}