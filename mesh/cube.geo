a = 0.2;
layers = 1;

Point(1) = {0, 0, 0, a};
Point(2) = {1, 0, 0, a};
Point(3) = {0, 1, 0, a};
Point(4) = {0, 0, 1, a};
Point(5) = {1, 1, 0, a};
Point(6) = {1, 0, 1, a};
Point(7) = {0, 1, 1, a};
Point(8) = {1, 1, 1, a};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 6};
//+
Line(3) = {6, 8};
//+
Line(4) = {8, 5};
//+
Line(5) = {5, 2};
//+
Line(6) = {2, 6};
//+
Line(7) = {7, 4};
//+
Line(8) = {3, 1};
//+
Line(9) = {3, 7};
//+
Line(10) = {7, 8};
//+
Line(11) = {1, 2};
//+
Line(12) = {3, 5};

//Transfinite Curve {1,2,3,4,5,6,7,8,9,10,11,12} = layers Using Progression 1;
//+
Curve Loop(1) = {12, 5, -11, -8};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {4, 5, 6, 3};
//+
Plane Surface(2) = {-2};
//+
Curve Loop(3) = {3, -10, 7, 2};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {9, 7, -1, -8};
//+
Plane Surface(4) = {-4};
//+
Curve Loop(5) = {11, 6, -2, -1};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {10, 4, -12, 9};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {6, 3, 2, 1, 5, 4};

//Transfinite Surface {1,2,3,4,5,6};
//+
Volume(1) = {1};

//Transfinite Volume {1};
//+
Physical Volume("Domain") = {1};
Physical Surface("Surfaces") = {1,2,3,4,5,6};
