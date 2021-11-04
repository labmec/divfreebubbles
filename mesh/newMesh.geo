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
nv = 11; pv = 1.15; pa = 1;
Transfinite Curve {1, 2, 7} = nh Using Progression ph;
Transfinite Curve {4,6} = nv Using Progression 1/pv;
Transfinite Curve {3,5} = nv Using Progression pv;
//+
Curve Loop(1) = {6, 1, 3, -7};
Plane Surface(1) = {1};
Curve Loop(2) = {4, 2, 5, 7};
Plane Surface(2) = {2};
Transfinite Surface {1, 2} ;
Recombine Surface{1,2};
//+
Physical Surface("Surface") = {1, 2};//1

Physical Curve("InjectionWell") = {1};//2
//+
Physical Curve("ProductionWell") = {2};//3
//+
Physical Curve("BottomLine") = {3};//4
//+
Physical Curve("TopLine") = {5};//5
//+
Physical Curve("LeftLine") = {6};//6
//+
Physical Curve("RightLine") = {4};//7


//+
//Physical Point("Point") = {5};
