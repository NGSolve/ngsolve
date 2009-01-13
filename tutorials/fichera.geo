#
## Fichera Cube
#
algebraic3d
solid c1 = plane (0, 0, 0; 0, 0, -1)
         and plane (0, 0, 0; 0, -1, 0)
         and plane (0, 0, 0; -1, 0, 0)
         and plane (1, 1, 1; 0, 0, 1)
         and plane (1, 1, 1; 0, 1, 0)
         and plane (1, 1, 1; 1, 0, 0) -bc=2;

solid c2 = plane (-0.5, -0.5, -0.5; 0, 0, -1)
         and plane (-0.5, -0.5, -0.5; 0, -1, 0)
         and plane (-0.5, -0.5, -0.5; -1, 0, 0)
         and plane (0.5, 0.5, 0.5; 0, 0, 1)
         and plane (0.5, 0.5, 0.5; 0, 1, 0)
         and plane (0.5, 0.5, 0.5; 1, 0, 0);

# cut off small cube

solid main = c1 and not c2 -bc=1;

tlo main;

