//+
lc = 0.025;
H = 1.0;

Point(1) = {H, -H, 0, lc};
//+
Point(2) = {H, H, 0, lc};
//+
Point(3) = {-H, H, 0, lc};
//+
Point(4) = {-H, -H, 0, lc};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Line("Fijo",1) = {4,2,3,1};
//+
Physical Surface("domain", 1) = {1};

