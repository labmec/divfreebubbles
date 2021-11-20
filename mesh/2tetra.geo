//+
a=10;
Point(1) = {0, 0, 0, a};
//+
Point(2) = {1, 0, 0, a};
//+
Point(3) = {0, 1, 0, a};
//+
Point(4) = {0, 0, 1, a};
Point(5) = {0.5, 0, 0.5, a};
//+
Line(1) = {1, 4};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {3, 1};
//+
Line(5) = {1, 2};


//+
Line(6) = {4, 5};
//+
Line(7) = {5, 2};
//+
Line(8) = {5, 1};
//+
Line(9) = {3, 5};
//+
Curve Loop(1) = {5, 2, 4};
//+
Plane Surface(1) = {-1};
//+
Curve Loop(2) = {2, 9, 7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {9, -6, -3};
//+
Plane Surface(3) = {-3};
//+
Curve Loop(4) = {3, -1, -4};
//+
Plane Surface(4) = {-4};
//+
Curve Loop(5) = {6, 8, 1};
//+
Plane Surface(5) = {-5};
//+
Curve Loop(6) = {7, -5, -8};
//+
Plane Surface(6) = {-6};
//+
Surface Loop(1) = {1, 6, 2, 3, 5, 4};
//+
Volume(1) = {1};
//+
Physical Volume("Domain") = {1};
//+
Physical Surface("S1") = {1,2,3,4,5,6};
//Physical Surface("S2") = {2,3}
//Physical Surface("S3") = {4};
//Physical Surface("S4") = {5,6};
