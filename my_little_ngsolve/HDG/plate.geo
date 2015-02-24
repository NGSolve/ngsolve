algebraic3d

solid bot = plane (0, 0, -0.01; 0, 0, -1);
solid top = plane (0, 0,  0.01; 0, 0,  1);

solid p1 = plane (0, 0, 0; -1,  0, 0) -bc=1;
solid p2 = plane (0, 0, 0;  0, -1, 0);
solid p3 = plane (1, 1, 0;  1,  0, 0) -bc=3;
solid p4 = plane (1, 1, 0;  0,  1, 0);

solid brick = bot and top and p1 and p2 and p3 and p4 -bc=2;

tlo brick -maxh=0.3;

identify closesurfaces bot top;
