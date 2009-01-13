algebraic3d
#
# a transformer
#
solid core = 
	    plane (-8,  0,  0; -1,  0,  0)
	and plane ( 8,  0,  0;  1,  0,  0)
	and plane ( 0, -6,  0;  0, -1,  0)
	and plane ( 0,  6,  0;  0,  1,  0)
	and plane ( 0,  0, -1;  0,  0, -1)
	and plane ( 0,  0,  1;  0,  0,  1)
	and not 
	(    plane (-6,  0,  0; -1,  0,  0)
	 and plane (-1,  0,  0;  1,  0,  0)
	 and plane ( 0, -4,  0;  0, -1,  0)
	 and plane ( 0,  4,  0;  0,  1,  0) )
	and not 
	(    plane ( 6,  0,  0;  1,  0,  0)
	 and plane ( 1,  0,  0; -1,  0,  0)
	 and plane ( 0, -4,  0;  0, -1,  0)
	 and plane ( 0,  4,  0;  0,  1,  0) );

solid coil1 =
	        cylinder (-7, -3, 0;-7, 3, 0; 3)
	and not cylinder (-7, -3, 0;-7, 3, 0; 2)
	and plane (0, -3, 0; 0, -1, 0) 
	and plane (0,  3, 0; 0,  1, 0);
solid coil2 =
	        cylinder ( 0, -3, 0; 0, 3, 0; 3)
	and not cylinder ( 0, -3, 0; 0, 3, 0; 2)
	and plane (0, -3, 0; 0, -1, 0) 
	and plane (0,  3, 0; 0,  1, 0);
solid coil3 =
	        cylinder ( 7, -3, 0; 7, 3, 0; 3)
	and not cylinder ( 7, -3, 0; 7, 3, 0; 2)
	and plane (0, -3, 0; 0, -1, 0) 
	and plane (0,  3, 0; 0,  1, 0);
	
solid box =
	    plane (-12, 0,  0; -1,  0,  0)
	and plane ( 12, 0,  0;  1,  0,  0)
	and plane (  0, 8,  0;  0,  1,  0)
	and plane (  0,-8,  0;  0, -1,  0)
	and plane (  0, 0,  5;  0,  0,  1)
	and plane (  0, 0, -5;  0,  0, -1);

solid air = box and not core and not coil1 and not coil2 and not coil3;

tlo coil1 -col=[0,1,0];
tlo coil2 -col=[0,1,0];
tlo coil3 -col=[0,1,0];
tlo air -col=[0,0,1] -transparent;
tlo core -col=[1,1,0];




