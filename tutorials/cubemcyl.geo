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
solid cyl = cylinder (0.5, 0.5, 0; 0.5, 0.5, 1; 0.03);

# cut off small cylinder from cube:

solid main = cube and not cyl;

tlo main;
