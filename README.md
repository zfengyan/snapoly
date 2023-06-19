Example usage:

`snapoly path-to-gpkg-file tolerance`

`snapoly D:\snapoly\data\examples\very_close_polygons.gpkg 0.01`

>>>
default snap rounding tolerance is: 0.01
the tolerance is set to: 0.01
reading polygons ...
        Path: D:\snapoly\data\examples\very_close_polygons.gpkg
        Type: GPKG
        Num of Layers: 1
        number of polygons: 2
        number of fields: 1
extent:
min X: 85732.3  max X: 85745.9
min Y: 446956   max Y: 446963
done
Time: 0.00186978min
inserting polygons to triangulation ...
done
Time: 2.79667e-06min
adding tags to triangulation ...
done
Time: 2.28833e-06min
snap rounding cases uner the current tolerance: 1
distances under the current tolerance:

snap rounding...
no more close vertex - vertex is found under given tolerance, snap vertex to boundary

no more snap rounding cases found
total snap rounded: 1
done

Time: 7.985e-05min
building polygons from constraints ...
done
Time: 0.0001428min
file snaprounded.gpkg has already existed, overwriting ...
-- output gpkg, write polygons
file saved at: snaprounded.gpkg
Time: 0.00099999min
area_diff: 2.82503e-05
number of faces: 18
number of vertices: 13
number of edges: 30
Time: 0.00326494min