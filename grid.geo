// 3D structured mesh template for basic configuration(BC)
// Author: Jithin Jith <j.jith@outlook.com>

lc=1; //characterisitc length 

t_cav_hub = 0.004; //depth of hub-side cavity
t_disk = 0.008; //depth of impeller disk
t_cav_shroud = 0.008; //depth of shroud-side cavity
t_diff = 0.005; //depth of diffuser cavity

h_cav_hub = t_cav_hub;
h_disk = h_cav_hub + t_disk;
h_cav_shroud = h_disk + t_cav_shroud;

h_diff_h = h_cav_hub + (t_disk - t_diff)/2;
h_diff_s = h_disk - (t_disk - t_diff)/2;

r_offset = 0.053; //axis offset from origin for revolution
r_disk = r_offset+0.122; //radius of disk from origin
r_gap = r_disk + 0.004; //radius after radial gap 
r_diff = 0.203; //radius at the end of diffuser cavity

nr_i = 13; //number of mesh nodes along the length of the impeller disk
nh_i = 3; //number of mesh nodes along the impeller depth

nr_g = 2; //number of lengthwise mesh nodes in the radial gap

nh_c_h = 2; //number of depthwise mesh nodes - cavity hub-side
nh_c_s = 3; //number of depthwise mesh nodes - cavity shroud-side

nr_d = 3; //number of lengthwise mesh nodes - diffuser
nh_d = 2; //number of depthwise mesh nodes - diffuser

n_ex = 5; //number of circumferential mesh divisions in a pi/2 circ. extrusion


//Points

Point(1) = {r_offset, 0, 0, lc};
Point(2) = {r_offset, 0, h_cav_hub, lc};
Point(3) = {r_offset, 0, h_disk, lc};
Point(4) = {r_offset, 0, h_cav_shroud, lc};
Point(5) = {r_gap, 0, h_cav_shroud, lc};
Point(6) = {r_gap, 0, h_diff_s, lc};
Point(7) = {r_diff, 0, h_diff_s, lc};
Point(8) = {r_diff, 0, h_diff_h, lc};
Point(9) = {r_gap, 0, h_diff_h, lc};
Point(10) = {r_gap, 0, 0, lc};

Point(11) = {r_disk, 0, 0, lc};
Point(12) = {r_disk, 0, h_cav_hub, lc};
Point(13) = {r_disk, 0, h_disk, lc};
Point(14) = {r_disk, 0, h_cav_shroud, lc};

//Lines

Line(101) = {1,2};
Line(102) = {2,3};
Line(103) = {3,4};
Line(104) = {4,14};
Line(105) = {14,5};
Line(106) = {5,6};
Line(107) = {6,7};
Line(108) = {7,8};
Line(109) = {8,9};
Line(110) = {9,10};
Line(111) = {10,11};
Line(112) = {11,1};
Line(113) = {3,13};
Line(114) = {13,12};
Line(115) = {12,2};
Line(116) = {12,11};
Line(117) = {14,13};
Line(118) = {9,6};

Transfinite Line{101} = nh_c_h Using Progression 1;
Transfinite Line{102} = nh_i Using Progression 1;
Transfinite Line{103} = nh_c_s Using Progression 1;
Transfinite Line{104} = nr_i Using Progression 1;
Transfinite Line{105} = nr_g Using Progression 1;
//Transfinite Line{106} = nh_c_s Using Progression 1;
Transfinite Line{107} = nr_d Using Progression 1;
Transfinite Line{108} = nh_d Using Progression 1;
Transfinite Line{109} = nr_d Using Progression 1;
//Transfinite Line{110} = nh_c_h Using Progression 1;
Transfinite Line{111} = nr_g Using Progression 1;
Transfinite Line{112} = nr_i Using Progression 1;
Transfinite Line{113} = nr_i Using Progression 1;
Transfinite Line{114} = nh_i Using Progression 1;
Transfinite Line{115} = nr_i Using Progression 1;
Transfinite Line{116} = nh_c_h Using Progression 1;
Transfinite Line{117} = nh_c_s Using Progression 1;
Transfinite Line{118} = nh_d Using Progression 1;

//Surfaces

//hub-side cavity
Line Loop(201) = {101,-115,116,112};
Plane Surface(301) = {201};
Transfinite Surface{301} = {1,2,12,11};

//disk
Line Loop(202) = {102,113,114,115};
Plane Surface(302) = {202};
Transfinite Surface{302} = {2,3,13,12};

//shroud-side cavity
Line Loop(203) = {103,104,117,-113};
Plane Surface(303) = {203};
Transfinite Surface{303} = {3,4,14,13};

//radial gap
Line Loop(204) = {-116,-114,-117,105,106,-118,110,111};
Plane Surface(304) = {204};
Transfinite Surface{304} = {11,14,5,10};

//diffuser
Line Loop(205) = {118,107,108,109};
Plane Surface(305) = {205};
Transfinite Surface{305} = {9,6,7,8};

Recombine Surface{301,302,303,304,305};


//Extrusion to volume

extrude_q1[] = Extrude{{0,0,1}, {0,0,0}, Pi/2}{
                      Surface{301,302,303,304,305};
                      Layers{n_ex};
                      Recombine; };

extrude_q2[] = Extrude{{0,0,1}, {0,0,0}, Pi/2}{
                      Surface{extrude_q1[0],extrude_q1[6],extrude_q1[12],extrude_q1[18],extrude_q1[28]};
                      Layers{n_ex};
                      Recombine; };

extrude_q3[] = Extrude{{0,0,1}, {0,0,0}, Pi/2}{
                      Surface{extrude_q2[0],extrude_q2[6],extrude_q2[12],extrude_q2[18],extrude_q2[28]};
                      Layers{n_ex};
                      Recombine; };

extrude_q4[] = Extrude{{0,0,1}, {0,0,0}, Pi/2}{
                      Surface{extrude_q3[0],extrude_q3[6],extrude_q3[12],extrude_q3[18],extrude_q3[28]};
                      Layers{n_ex};
                      Recombine; };

//Label important physical surfaces/volumes

Physical Volume(101) = {2,7,12,17}; //disk
Physical Volume(102) = {1,3,4,5,6,8,9,10,11,13,14,15,16,18,19,20}; //cavity
Physical Surface(201) = {336, 466, 596, 725}; //disk inner fixed face
