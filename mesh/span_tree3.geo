a = 0.1;
layers1=3;
layers2=5;
layers3=7;
r1 = 0.25;
c2 = 1.25;

Point(1) = {0,0.5,0, a};
Point(2) = {0.5,0.5,0, a};
Point(3) = {0.5,0,0, a};
Point(4) = {0.5,1,0, a};
Point(5) = {0.5,0.5+r1,0, a};
Point(6) = {0.5-r1,0.5,0, a};
Point(7) = {0.5,0.5-r1,0, a};
Point(8) = {0.5+r1,0.5,0, a};
Point(9) = {c2,0,0, a};
Point(10) = {c2,1,0, a};
Point(11) = {c2+0.5,.5,0, a};
Point(12) = {c2,0.5,0, a};
Point(13) = {c2,0.5-r1,0, a};
Point(14) = {c2,0.5+r1,0, a};
Point(15) = {c2-r1,0.5,0, a};
Point(16) = {c2+r1,0.5,0, a};


Circle(1) = {4, 2, 1};
//+
Circle(2) = {1, 2, 3};
//+
Circle(3) = {5, 2, 6};
//+
Circle(4) = {6, 2, 7};
//+
Circle(5) = {7, 2, 8};
//+
Circle(6) = {8, 2, 5};
//+
Circle(7) = {9, 12, 11};
//+
Circle(8) = {11, 12, 10};
//+
Line(9) = {4, 10};
//+
Line(10) = {9, 3};
//+
Circle(11) = {13, 12, 16};
//+
Circle(12) = {16, 12, 14};
//+
Circle(13) = {14, 12, 15};
//+
Circle(14) = {15, 12, 13};
//+
Curve Loop(1) = {2, -10, 7, 8, -9, 1};
//+
Curve Loop(2) = {5, 6, 3, 4};
//+
Curve Loop(3) = {14, 11, 12, 13};
//+
Plane Surface(1) = {1, 2, 3};
Recombine Surface{1};
//+
Extrude {0, 0, 0.05} {
  Surface{1}; 
}

Recombine Volume {1};
//+
Physical Surface("Surfaces") = {77, 73, 85, 81, 69, 57, 61, 65, 86, 1, 49, 53, 45, 41, 37, 33};
//+
Physical Volume("Domain") = {1};
