//Points
a = 0.2;
Point(1) = {0, 0, 0, a};
Point(2) = {0, 1, 0, a};
Point(3) = {1, 1, 0, a};
Point(4) = {1, 0, 0, a};

//Lines
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};

//Transfinite Lines
nh = 11; ph = 1.;
nv = 11; pv = 1.;
//Horizontal Lines
Transfinite Curve {1, 3} = nh Using Progression ph;
//Vertical Lines
Transfinite Curve {4, 2} = nv Using Progression pv;

//Surface
Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};
Transfinite Surface {1};
//Uncomment for quadrilateral elements
//Recombine Surface {1};

//Physical groups
Physical Point("BottomLeftPoint") = {1}; //1
Physical Point("BottomRightPoint") = {4};//2
Physical Point("TopRightPoint") = {3};   //3
Physical Point("TopLeftPoint") = {2};	 //4
Physical Curve("BottomLine") = {1}; 	 //5
Physical Curve("TopLine") = {3};	 //6
Physical Curve("RightLine") = {2};	 //7
Physical Curve("LeftLine") = {4};	 //8
Physical Surface("Surface") = {1};	 //9

