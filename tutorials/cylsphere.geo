#
## Cylinder and Spehre
#
algebraic3d
solid cyl = cylinder ( 3, 0, 0; -3, 0, 0; 0.5 )
	and plane (-2, 0, 0; -1, 0, 0)
	and plane (2, 0, 0; 1, 0, 0);
solid sp = sphere (0, 0, 0; 1);

solid main = sp or cyl;

tlo main;
