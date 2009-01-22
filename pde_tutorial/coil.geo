algebraic3d
solid coil = cylinder (0, 0, -1; 0, 0, 1; 0.4)
	and not cylinder  (0, 0, -1; 0, 0, 1; 0.2)
	and plane (0, 0, 0.4; 0, 0, 1)
	and plane (0, 0, -0.4; 0, 0, -1) -bc=2;

solid box = orthobrick (-2, -2, -2; 3, 2, 2) -bc=1;

solid air =  box and not coil -bc=3;

tlo coil -col=[0, 1, 0] -material=matcoil; 
tlo air  -col=[0, 0, 1] -transparent -material=air;
