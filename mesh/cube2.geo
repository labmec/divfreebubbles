a = 0.5;
layers = 3;
layers2 = 3;
layers3 = 3;
p1 =0;
p2 = 1;
Point(1) = {p1,p1,p1, a};
Point(2) = {p2,p1,p1, a};
Point(3) = {p1,p2,p1, a};
Point(4) = {p1,p1,p2, a};
Point(5) = {p2,p2, p1, a};
Point(6) = {p2, p1,p2, a};
Point(7) = {p1,p2,p2, a};
Point(8) = {p2,p2,p2, a};
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

Transfinite Curve {3,5,7,8} = layers Using Progression 1;
Transfinite Curve {2,10,11,12} = layers2 Using Progression 1;
Transfinite Curve {1,4,6,9} = layers3 Using Progression 1;

//+
Curve Loop(1) = {-7, -2, -3, 10};
//+
Plane Surface(1) = {-1};
//+
Curve Loop(2) = {6, 3, 4, 5};
//+
Plane Surface(2) = {-2};
//+
Curve Loop(3) = {5, -11, -8, 12};
//+
Plane Surface(3) = {-3};
//+
Curve Loop(4) = {8, 1, -7, -9};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {9, 10, 4, -12};
//+
Plane Surface(5) = {-5};
//+
Curve Loop(6) = {2, -6, -11, 1};
//+
Plane Surface(6) = {-6};
//+
Surface Loop(1) = {4, 3, 2, 6, 1, 5};
//+
Transfinite Surface {1,2,3,4,5,6};
//Recombine Surface {1,2,3,4,5,6};
Volume(1) = {1};
Transfinite Volume {1};
//Recombine Volume {1};
//+
Physical Volume("Domain") = {1};
//+
Physical Surface("Surfaces") = {3, 6, 2, 5, 1, 4};
