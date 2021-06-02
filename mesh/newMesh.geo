// Gmsh project created on Wed Jun  2 14:01:13 2021
a = 0.02;
radius = 0.05;
Point(1) = {0, 0, 0, a};
Point(2) = {0, 1, 0, a};
Point(3) = {1, 1, 0, a};
Point(4) = {1, 0, 0, a};
Point(5) = {0+radius, 0, 0, a};
Point(6) = {0, 0+radius, 0, a};
Point(7) = {1-radius, 1, 0, a};
Point(8) = {1, 1-radius, 0, a};
//+
Circle(1) = {6,1,5};
Circle(2) = {8,3,7};
Line(3) = {5, 4};
Line(4) = {4, 8};
Line(5) = {7, 2};
Line(6) = {2, 6};
Line(7) = {2, 4};

nh = 11; ph = 1; 
nv = 11; pv = 1.; pa = 1;
Transfinite Curve {1, 2, 7} = nh Using Progression ph;
Transfinite Curve {3,4,5,6} = nh Using Progression 1/ph;
//+
Curve Loop(1) = {6, 1, 3, -7};
Plane Surface(1) = {-1};
Curve Loop(2) = {4, 2, 5, 7};
Plane Surface(2) = {-2};
Transfinite Surface {1, 2} ;
Recombine Surface{1,2};
//+
Physical Curve("InjectionWell") = {1};
//+
Physical Curve("ProductionWell") = {2};
//+
Physical Curve("Bottom") = {3};
//+
Physical Curve("Top") = {5};
//+
Physical Curve("Left") = {6};
//+
Physical Curve("Right") = {4};
//+
Physical Surface("Surface") = {1, 2};
//+
Physical Point("Point") = {6};
