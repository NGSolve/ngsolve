#
## Cube minus Cylinder
#
algebraic3d

solid cube = plane (0, 0, 0; 0, 0, -1)
         and plane (0, 0, 0; 0, -1, 0)
         and plane (0, 0, 0; -1, 0, 0)
         and plane (1, 1, 1; 0, 0, 1)
         and plane (1, 1, 1; 0, 1, 0)
         and plane (1, 1, 1; 1, 0, 0);
solid sp = sphere (0.5, 0.5, 0; 0.001);

# cut off small sphere:

solid main = cube and not sp;

tlo main;

