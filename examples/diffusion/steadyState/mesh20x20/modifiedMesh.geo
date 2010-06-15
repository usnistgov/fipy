cellSize = 0.5;
Point(2) = {0, 0, 0, cellSize};
Point(3) = {20, 0, 0, cellSize};
Point(4) = {20, 20, 0, cellSize};
Point(5) = {0, 20, 0, cellSize};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 5};
Line(9) = {5, 2};
Line Loop(10) = {6, 7, 8, 9};
Plane Surface(11) = {10};

