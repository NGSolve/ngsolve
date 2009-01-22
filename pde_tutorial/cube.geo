#
## A cube
#
algebraic3d

# cube consisting of 6 planes:

solid p1 = plane (0, 0, 0; -1, 0, 0) -bc=1;
solid p2 = plane (0, 0, 0; 0, -1, 0) -bc=1;
solid p3 = plane (0, 0, 0; 0, 0, -1) -bc=1;
solid p4 = plane (1, 1, 1; 1, 0, 0) -bc=1;
solid p5 = plane (1, 1, 1; 0, 1, 0) -bc=1;
solid p6 = plane (1, 1, 1; 0, 0, 1) -bc=2;

solid cube = p1 and p2 and p3 and p4 and p5 and p6;

tlo cube;


