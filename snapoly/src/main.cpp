#include "pch.h"
#include "io.h"



/*
* if there are more than one polygons
* the output would be two polygons
* when we write to gpkg, we need to check the polygons' relationship
* i.e. which is the biggest polygon -> the biggest polygon should contain the other polygons
*/



int main(int argc, char** argv)
{
	//cmdline::parser p;
	//p.set_program_name("snapoly(v1.0)"); // set the program name in the console
	//p.add<std::string>("dataset", 'd', "input dataset (GPKG)", true, ""); // dataset file

	/* run parser -----------------------------------------------------------------------------------------------------------*/
	/*bool ok = p.parse(argc, argv);

	if (argc == 1 || p.exist("help")) {
		std::cerr << p.usage();
		return 0;
	}

	if (!ok) {
		std::cerr << p.error() << std::endl << p.usage();
		return 0;
	}*/
	/* ----------------------------------------------------------------------------------------------------------------------*/


	// Snap rounding
	Snap_rounding_2 sr;
	sr.set_tolerance(0.01); // for Delft: 0.01m


	/* ----------------------------------------------------------------------------------------------------------------------*/
	// usage:
	// snapoly 100 0.01 -- process 100 polygons with the tolerance 0.01
	// snapoly 1000 0.02 -- process 1000 polygons with the tolerance 0.02
	
	if (argc < 2) {
		std::cout << "default tolerance is set to 0.01m" << std::endl;
	}
	else {
		try {

			// how many polygons will be processed
			// long long control_numbers = std::stod(argv[1]);

			// the value of the tolerance, default: 0.01
			double tol = std::stod(argv[2]);
			sr.set_tolerance(tol);
		}
		catch (const std::exception& e) {
			std::cout << "Error: Invalid double value provided." << e.what() << std::endl;
			return 1;
		}
	}
	/* ----------------------------------------------------------------------------------------------------------------------*/

	const char* input_file = R"(D:\snapoly\data\Faroe_islands\faroe.gpkg)";
	const char* tri_file = R"(D:\snapoly\data\Faroe_islands\benchmarking\faroe_tri.gpkg)";
	const char* output_boundaries_file = R"(D:\snapoly\data\Faroe_islands\benchmarking\faroe_boundaries.gpkg)";
	const char* res_file = R"(D:\snapoly\data\Faroe_islands\benchmarking\faroe_res.gpkg)";

	// Timer
	Timer timer;

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

	sr.snap_rounding(); // not using remove dangles()


	//double minimum_distance = sr.minimum_distance();
	//cout << "tolerance: " << sr.tolerance() << '\n';
	//cout << " minimum distance under the current tolerance: " << minimum_distance << '\n';

	//io::export_to_gpkg(tri_file, sr.triangulation());

	//io::export_to_gpkg(output_boundaries_file, sr.constraintsWithInfo());

	//cout << "file saved at: " << output_file << '\n';

	// -----------------------------------------------------------------------------------------

	io::build_polygons_from_constraints(sr.constraintsWithInfo(), sr.result_polygons());

	io::export_to_gpkg(res_file, sr.result_polygons());

	sr.measure_distortions(); // this function must be called after the io::build_polygons_from_constraints() function

	// numbers

	// number of faces
	long long count_faces = 0;
	for (auto const& face : sr.triangulation().finite_face_handles())
		++count_faces;
	cout << "number of faces: " << count_faces << '\n';

	// number of vertices
	long long count_vertices = 0;
	for (auto const& vh : sr.triangulation().finite_vertex_handles())
		++count_vertices;
	cout << "number of vertices: " << count_vertices << '\n';

	// number of edges
	long long count_edges = 0;
	for (auto const& vh : sr.triangulation().finite_edges())
		++count_edges;
	cout << "number of edges: " << count_edges << '\n';
	// -----------------------------------------------------------------------------------------



	return 0;


}