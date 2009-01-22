algebraic3d

solid p1 =  plane (0, 0, 0; -1, 0, 0);

solid cube = plane (0, 0, 0; 0, 0, -1) 
         and plane (0, 0, 0; 0, -1, 0)
         and p1
         and plane (10, 2, 1; 0, 0, 1)
         and plane (10, 2, 1; 0, 1, 0)
         and plane (10, 2, 1; 1, 0, 0) -bc=2 -maxh=2;

tlo cube;

#singular edge p1 cube;
