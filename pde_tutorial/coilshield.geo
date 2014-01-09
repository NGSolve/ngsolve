algebraic3d

solid box = orthobrick (-2, -2, -2; 3, 2, 2) -bc=1;

solid pl1 = plane (1.1, 0, 0; 1, 0, 0);
solid pl2 = plane (1, 0, 0; -1, 0, 0);

solid shield = orthobrick (0, -1.5, -1.5; 2, 1.5, 1.5)
	and pl1 and pl2;

solid coil = cylinder (0, 0, -1; 0, 0, 1; 0.4)
	and not cylinder  (0, 0, -1; 0, 0, 1; 0.2)
	and plane (0, 0, 0.4; 0, 0, 1)
	and plane (0, 0, -0.4; 0, 0, -1);

solid air =  box and not coil and not shield -bc=2;

tlo shield -col=[1, 0, 0] -material=iron;
tlo coil -col=[0, 1, 0]   -material=coil; 
tlo air  -col=[0, 0, 1] -transparent -material=air; 

identify closesurfaces pl1 pl2;


