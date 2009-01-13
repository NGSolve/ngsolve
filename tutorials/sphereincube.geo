algebraic3d
#
# Example with two sub-domains: 
#
solid cube = plane (0, 0, 0; 0, 0, -1)
         and plane (0, 0, 0; 0, -1, 0)
         and plane (0, 0, 0; -1, 0, 0)
         and plane (1, 1, 1; 0, 0, 1)
         and plane (1, 1, 1; 0, 1, 0)
         and plane (1, 1, 1; 1, 0, 0);
solid sph = sphere (0.5, 0.5, 0.5; 0.3);

solid rest = cube and not sph;

tlo rest -transparent -col=[0,0,1];
tlo sph -col=[1,0,0];

