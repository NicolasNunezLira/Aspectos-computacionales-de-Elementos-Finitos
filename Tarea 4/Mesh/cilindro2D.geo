//+
lc = 0.025;
lc2= 0.005;
L  = 2.2;
H  = 0.41;

Point(1) = {0, 0, 0, lc};
//+
Point(2) = {L, 0, 0, lc};
//+
Point(3) = {L, H, 0, lc};
//+
Point(4) = {0, H, 0, lc};
//
//+
Point(5) = {0.2, 0.2, 0, lc2};
//+
Point(6) = {0.25, 0.2, 0, lc2};
//+
Point(7) = {0.2, 0.25, 0, lc2};
//+
Point(8) = {0.15, 0.2, 0, lc2};
//+
Point(9) = {0.2, 0.15, 0, lc2};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Circle(5) = {6, 5, 7};
//+
Circle(6) = {7, 5, 8};
//+
Circle(7) = {8, 5, 9};
//+
Circle(8) = {9, 5, 6};
//+
Line Loop(1) = {3, 4, 1, 2};
//+
Line Loop(2) = {8, 5, 6, 7};
//+
Plane Surface(1) = {1, 2};
//+
Physical Line("inlet", 1) = {4};
//+
Physical Line("wall", 2) = {3, 1, 8, 7, 6, 5};
//+
Physical Line("outlet", 3) = {2};
//+
Physical Surface("domain", 1) = {1};
