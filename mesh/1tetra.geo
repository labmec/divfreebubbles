//+
a=21;
Point(1) = {0, 0, 0, a};
//+
Point(2) = {1, 0, 0, a};
//+
Point(3) = {0, 1, 0, a};
//+
Point(4) = {0, 0, 1, a};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {3, 1};
//+
Line(6) = {1, 2};
//+
Curve Loop(1) = {4, 2, 3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, 5, 6};
//+
Plane Surface(2) = {-2};
//+
Curve Loop(3) = {4, -1, -5};
//+
Plane Surface(3) = {-3};
//+
Curve Loop(4) = {2, -6, 1};
//+
Plane Surface(4) = {-4};
//+
Surface Loop(1) = {1, 3, 4, 2};
//+
Volume(1) = {1};
//+
Physical Volume("Domain") = {1};
//+
Physical Surface("Surfaces") = {1, 2, 3, 4};


