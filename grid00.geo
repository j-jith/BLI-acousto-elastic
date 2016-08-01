lc=1;

nr_1 = 21;
n_ex= 10;
nh_1 = 2; //fluid
nh_2 = 2; //solid
nd = 3;

Point(1) = {0.05, 0, 0, lc};
Point(2) = {0.05, 0, 0.006, lc};
Point(3) = {0.05, 0, 0.014, lc};
Point(4) = {0.05, 0, 0.02, lc};
Point(5) = {0.177, 0, 0.02, lc};
Point(6) = {0.177, 0, 0.014, lc};
Point(7) = {0.227, 0, 0.014, lc};
Point(8) = {0.227, 0, 0.006, lc};
Point(9) = {0.177, 0, 0.006, lc};
Point(10) = {0.177, 0, 0, lc};

Line(1) = {1,2};
Transfinite Line{1} = nh_1 Using Progression 1;
Line(2) = {2,3};
Transfinite Line{2} = nh_2 Using Progression 1;
Line(3) = {3,4};
Transfinite Line{3} = nh_1 Using Progression 1;
Line(4) = {4,5};
Transfinite Line{4} = nr_1 Using Progression 1;
Line(5) = {5,6};
Transfinite Line{5} = nh_1 Using Progression 1;
Line(6) = {6,7};
Transfinite Line{6} = nd Using Progression 1;
Line(7) = {7,8};
Transfinite Line{7} = nh_2 Using Progression 1;
Line(8) = {8,9};
Transfinite Line{8} = nd Using Progression 1;
Line(9) = {9,10};
Transfinite Line{9} = nh_1 Using Progression 1;
Line(10) = {10,1};
Transfinite Line{10} = nr_1 Using Progression 1;
Line(11) = {6,9};
Transfinite Line{11} = nh_2 Using Progression 1;

Line Loop(12) = {1,2,3,4,5,11,9,10};
Plane Surface(13) = {12};
Transfinite Surface{13} = {1,4,5,10};

Line Loop(14) = {6,7,8,-11};
Plane Surface(15) = {14};
Transfinite Surface{15} = {6,7,8,9};

Recombine Surface{13,15};

extrude_1[] = Extrude{{0,0,1}, {0,0,0}, Pi/2}{Surface{13}; Layers{n_ex}; Recombine;};
extrude_2[] = Extrude{{0,0,1}, {0,0,0}, Pi/2}{Surface{extrude_1[0]}; Layers{n_ex}; Recombine;};
extrude_3[] = Extrude{{0,0,1}, {0,0,0}, Pi/2}{Surface{extrude_2[0]}; Layers{n_ex}; Recombine;};
extrude_4[] = Extrude{{0,0,1}, {0,0,0}, Pi/2}{Surface{extrude_3[0]}; Layers{n_ex}; Recombine;};

extrude_5[] = Extrude{{0,0,1}, {0,0,0}, Pi/2}{Surface{15}; Layers{n_ex}; Recombine;};
extrude_6[] = Extrude{{0,0,1}, {0,0,0}, Pi/2}{Surface{extrude_5[0]}; Layers{n_ex}; Recombine;};
extrude_7[] = Extrude{{0,0,1}, {0,0,0}, Pi/2}{Surface{extrude_6[0]}; Layers{n_ex}; Recombine;};
extrude_8[] = Extrude{{0,0,1}, {0,0,0}, Pi/2}{Surface{extrude_7[0]}; Layers{n_ex}; Recombine;};

Physical Surface(1) = {74, 32, 116, 158}; //solid fixed end
Physical Volume(9) = {1, 2, 3, 4, 5, 6, 7, 8}; //full volume
