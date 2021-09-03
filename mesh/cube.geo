layers =1;

Point(1) = {0, 0, 0, 1.0};
Extrude {1, 0, 0} {
  Point{1}; Layers{layers};//Recombine;
}
Extrude {0, 1, 0} {
  Curve{1}; Layers{layers};// Recombine;
}
Extrude {0, 0, 1} {
  Surface{5}; Layers{layers};// Recombine;
}

//+
Physical Volume("Domain") = {1};
//+
Physical Surface("Surfaces") = {22, 5, 26, 18, 14, 27};
//+
//Physical Curve("Lines") = {21, 2, 3, 17, 9, 8, 10, 12, 7, 1, 13, 4};
