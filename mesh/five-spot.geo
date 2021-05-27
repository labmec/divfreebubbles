//Points
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

Point(9) = {0.5, 0.5, 0, a};
Point(10) = {radius*Sqrt(2)/2, radius*Sqrt(2)/2, 0, a};
Point(11) = {1-radius*Sqrt(2)/2, 1-radius*Sqrt(2)/2, 0, a};

//Lines
Line(1) = {2, 6};
Line(2) = {5, 4};
Line(3) = {4, 8};
Line(4) = {7, 2};
Line(5) = {2, 9};
Line(6) = {9, 4};
Circle(7) = {5, 1, 10};
Circle(8) = {10, 1, 6};
Circle(9) = {7, 3, 11};
Circle(10) = {11, 3, 8};
Line(11) = {10, 9};
Line(12) = {9, 11};

//Transfinite Lines
nh = 31; ph = 1.1; 
nv = 11; pv = 1.; pa = 1.2;
//External Lines
Transfinite Curve {2, 4, 11} = nh Using Progression ph;
Transfinite Curve {1, 3, 12} = nh Using Progression 1/ph;
Transfinite Curve {5} = nv Using Progression 1/pa;
Transfinite Curve {6} = nv Using Progression pa;
//Curves
Transfinite Curve {7, 8, 9, 10} = nv Using Progression pv;

//Surfaces
Curve Loop(1) = {-1, 8, -11, 5};
Plane Surface(1) = {1};
Curve Loop(2) = {11, 6, -2, 7};
Plane Surface(2) = {2};
Curve Loop(3) = {-6, -3, 10, 12};
Plane Surface(3) = {3};
Curve Loop(4) = {9, -12, -5, -4};
Plane Surface(4) = {4};
Transfinite Surface {1, 2} Left;
Transfinite Surface {3, 4} Right;
Recombine Surface{1,2,3,4};

//Physical Groups
Physical Curve("InjectionWell") = {8, 7};   //1
Physical Curve("ProductionWell") = {9, 10}; //2
Physical Curve("BottomLine") = {2};         //3
Physical Curve("TopLine") = {4};            //4
Physical Curve("LeftLine") = {1};           //5
Physical Curve("RightLine") = {3};          //6
Physical Surface("Surface") = {1, 2, 3, 4}; //7
Physical Point("Point") = {10};		    //8
