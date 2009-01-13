##
## Example with periodic boundary conditions
##    by Joachim Schoeberl
##
##

algebraic3d

solid p1 = plane (0, 0, 0; 0, 0, -1);
solid p2 = plane (1, 1, 1; 0, 0, 1);

solid p3 = plane (0, 0, 0; 0, -1, 0);
solid p4 = plane (1, 1, 1; 0, 1, 0);

solid p5 = plane (0, 0, 0; -1, 0, 0);
solid p6 = plane (1, 1, 1; 1, 0, 0);

 
solid cube = p1 and p2 and p3 and p4 and p5 and p6;

solid cyls = 
	cylinder (0.5, 0.5, 0; 0.5, 0.5, 1; 0.3)
	or cylinder (0, 0.5, 0.2; 1, 0.5, 0.2; 0.1);
solid matrix = cube and not cyls;
solid inner = cube and cyls;

tlo matrix -transparent;
tlo inner -col=[1,0,0];

identify periodic p1 p2;
identify periodic p3 p4;
identify periodic p5 p6;

