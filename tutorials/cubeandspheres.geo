#
## Cube and Spheres
#
algebraic3d

# a cube
solid cube = plane (0, 0, 0; 0, 0, -1)
         and plane (0, 0, 0; 0, -1, 0)
         and plane (0, 0, 0; -1, 0, 0)
         and plane (1, 1, 1; 0, 0, 1)
         and plane (1, 1, 1; 0, 1, 0)
         and plane (1, 1, 1; 1, 0, 0);

# two shperes
solid sph1 = sphere (0.5, 0.5, 0.5; 0.58);
solid sph2 = sphere (0.5, 0.5, 0.5; 0.75);

# cut cube with inner and outer sphere
solid main = cube and sph2 and not sph1;

tlo main;
