algebraic3d

solid ocyl = cylinder (0, 0, 0; 1, 0, 0; 0.3);
solid icyl = cylinder (0, 0, 0; 1, 0, 0; 0.29);

solid p1 = plane (0, 0, 0; -1,  0, 0) -bc=1;
solid p2 = plane (1, 0, 0;  1,  0, 0) -bc=3;

solid shell = ocyl and not icyl and p1 and p2 -bc=2;

tlo shell -maxh=0.1;

identify closesurfaces ocyl icyl;
