// Gmsh project created on Mon May 24 18:40:06 2021
//+
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Curve {1,3} = 2 Using Progression 1;
Transfinite Curve {2,4} = 2 Using Progression 1;
Transfinite Surface {1};
Recombine Surface{1};


//+
Physical Surface("1") = {1};
Physical Curve("2") = {1};//Bottom Line
Physical Curve("3") = {2};//Right
Physical Curve("4") = {3};//Top
Physical Curve("5") = {4};//Left
Physical Point("6") = {1};
