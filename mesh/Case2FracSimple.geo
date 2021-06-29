// Gmsh project created on Mon Jun 28 11:23:31 2021
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.0};
//+
Point(2) = {2, 0, 0, 0.0};
//+
Point(3) = {2, 2, 0, 0.0};
//+
Point(4) = {0, 2, 0, 0.0};
//+
Point(5) = {1, 0, 0, 0.0};
//+
Point(6) = {1, 2, 0, 0.0};
//+
Line(1) = {1, 5};
//+
Line(2) = {5, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 6};
//+
Line(5) = {6, 5};
//+
Line(6) = {6, 4};
//+
Line(7) = {4, 1};
//+
Curve Loop(1) = {7, 1, -5, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, 3, 4, 5};
//+
Plane Surface(2) = {2};
//+
Physical Curve("inlet", 8) = {7};
//+
Physical Curve("outlet", 9) = {3};
//+
Physical Curve("noflux", 10) = {6, 4, 2, 1};
//+
Physical Curve("intersection", 11) = {5};
//+
Physical Surface("frac", 12) = {1, 2};

//+
Point(7) = {1, 0, -1, 0.0};
//+
Point(8) = {1, 2, -1, 0.0};
//+
Point(9) = {1, 0, 1, 0.0};
//+
Point(10) = {1, 2, 1, 0.0};
//+
Line(8) = {6, 10};
//+
Line(9) = {10, 9};
//+
Line(10) = {9, 5};
//+
Line(11) = {5, 7};
//+
Line(12) = {7, 8};
//+
Line(13) = {8, 6};
//+
Curve Loop(3) = {9, 10, -5, 8};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {5, 11, 12, 13};
//+
Plane Surface(4) = {4};
//+
Physical Curve("noflux", 10) += {8, 9, 10, 11, 12, 13};
//+
Physical Surface("frac", 12) += {3, 4};

Coherence Mesh;
Transfinite Curve{:} = 2;
Transfinite Surface{:};
Transfinite Volume{:};
Recombine Surface{:};
Recombine Volume{:};//+