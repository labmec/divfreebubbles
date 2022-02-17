a = 0.5;
layers1=3;
layers2=5;
layers3=7;

Point(1) = {0,0,0, a};

//+
Extrude {.1, 0, 0} {
  Point{1}; Layers{layers1};Recombine;
}
Extrude {.3, 0, 0} {
  Point{2}; Layers{layers3};Recombine;
}
Extrude {.1, 0, 0} {
  Point{3}; Layers{layers1};Recombine;
}
Extrude {.3, 0, 0} {
  Point{4}; Layers{layers3};Recombine;
}
Extrude {.1, 0, 0} {
  Point{5}; Layers{layers1};Recombine;
}

Extrude {0,0.1, 0} {
  Line{1,2,3,4,5}; Layers{layers1};Recombine;
}
Extrude {0,0.2, 0} {
  Line{6,10,18,14,22}; Layers{layers2};Recombine;
}
Extrude {0,0.1, 0} {
  Line{26,30,38,34,42}; Layers{layers1};Recombine;
}
Extrude {0,0, 0.1} {
  Surface{9,13,17,21,25,29,41,45,49,53,57,61,65}; Layers{layers1};Recombine;
}
//+
Physical Volume("Surfaces") += {9, 10, 11, 7, 6, 1, 2, 3, 4, 5, 8, 13, 12};
