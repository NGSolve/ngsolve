#
## An elliptic cylinder
#
algebraic3d

solid cutcone = ellipticcylinder ( 0, 0, 0; 1, 0, 0; 0, 0.5, 0)
	and plane (0, 0, 0; 0, 0, -1)
	and plane (0, 0, 1; 0, 0, 1);

tlo cutcone;
