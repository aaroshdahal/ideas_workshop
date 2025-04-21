// This code was created by pygmsh v0.7.5.
SetFactory("OpenCASCADE");
// Mesh.CharacteristicLengthMin = 1.0;
// Mesh.CharacteristicLengthMax = 1.0;
Point(1) = {0.0, 5.0,  0.0, 2.0};
Point(2) = {30.0, 5.0,  0, 2.0};
Point(3) = {30.0, -5.0, 0, 2.0};
Point(4) = {0.0, -5.0, 0, 2.0};
Point(5) = {0.0, 0.00, 0, 2.0};
Point(6) = {5.0, 0.00, 0, 2.0};
//+
Line(7) = {1, 2};
Line(8) = {2, 3};
Line(9) = {3, 4};
Line(10) = {4, 5};
Line(11) = {5, 1};
Line Loop(12) = {7, 8, 9, 10, 11};
Plane Surface(13) = {12};
Physical Surface(14) = {13};
//+
Line(15) = {5, 6};
Physical Curve(16) = {15};
Physical Point(17) = {5};
Curve{15} In Surface{13};
//+
Field[1] = Box;
Field[1].VIn = 0.05;
Field[1].VOut = 0.5;
Field[1].XMin = 0.0;
Field[1].XMax = 30.0;
Field[1].YMin = -1.0;
Field[1].YMax = 1.0;
Field[1].ZMin = 0.0;
Field[1].ZMax = 0.0;
Field[1].Thickness = 0.25;
//+
Field[2] = Box;
Field[2].VIn = 0.5;
Field[2].VOut = 0.5;
Field[2].XMin = 0.0;
Field[2].XMax = 30.0;
Field[2].YMin = -5.0;
Field[2].YMax = 5.0;
Field[2].ZMin = 0.0;
Field[2].ZMax = 0.0;
Field[2].Thickness = 0.25;
//+
Field[6] = Min;
Field[6].FieldsList = {1,2};
Background Field = 6;
//+
Mesh.ElementOrder = 1;
Mesh 2;
Mesh.MshFileVersion = 2.0;
//+
Plugin(Crack).Dimension = 1;
Plugin(Crack).PhysicalGroup = 16;
Plugin(Crack).OpenBoundaryPhysicalGroup = 17;
Plugin(Crack).Run;
//+
Save "surfing_mesh.msh";

