algebraic3d

#
# two intersecting cylinderes
#

solid cyl1 = cylinder ( 1, 0, 0; -1, 0, 0; 0.5 )
	and plane (-1, 0, 0; -1, 0, 0)
	and plane (1, 0, 0; 1, 0, 0);
solid cyl2 = cylinder ( 0, 1, 0.3; 0, -1, 0.3; 0.5 )
	and plane (0, -1, 0; 0, -1, 0)
	and plane (0, 1, 0; 0, 1, 0);
solid main = cyl1 or cyl2;


tlo main;
