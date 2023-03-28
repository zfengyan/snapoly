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
	
	const char* input_file = R"(D:\snapoly\data\netherlands\dencase.gpkg)";
	const char* tri_file = R"(D:\snapoly\data\netherlands\dencase_tri.gpkg)";
	const char* output_boundaries_file = R"(D:\snapoly\data\netherlands\dencase_boundaries.gpkg)";
	const char* res_file = R"(D:\snapoly\data\netherlands\dencase_res.gpkg)";

	// Snap rounding
	Snap_rounding_2 sr;
	sr.set_tolerance(0.0);

	io::add_polygons_from_input_file(input_file, sr.polygons());

	sr.insert_polygons_to_triangulation();

	sr.add_tag_to_triangulation();

	// find the minimum tolerance ---------------------------------------------
	std::priority_queue<double> lengthQueue;
	sr.find_tolerance(lengthQueue);
	cout << "distances under the current tolerance: " << '\n';
	while (!lengthQueue.empty()) {
		cout << lengthQueue.top() << '\n';
		lengthQueue.pop();
	}
	// find the minimum tolerance ---------------------------------------------

	sr.snap_rounding(); 

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

	std::set<Constraint> constraintsSet;

	Constraint c1(CDTPoint(0.0, 0.0), CDTPoint(1.0, 0.0));
	Constraint c2(CDTPoint(2.0, 0.0), CDTPoint(0.0, 0.0));

	constraintsSet.insert(c1);
	constraintsSet.insert(c2);

	cout << constraintsSet.size() << '\n';

	auto it = constraintsSet.find(c1);
	if (it != constraintsSet.end())
		cout << "found!\n";


	return 0;


}