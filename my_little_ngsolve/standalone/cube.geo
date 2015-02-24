#
## A cube
#
algebraic3d

# cube consisting of 6 planes:

solid cube = plane (0, 0, 0; 0, 0, -1) 
         and plane (0, 0, 0; 0, -1, 0)
         and plane (0, 0, 0; -1, 0, 0)
         and plane (1, 1, 1; 0, 0, 1)
         and plane (1, 1, 1; 0, 1, 0)
         and plane (1, 1, 1; 1, 0, 0);

tlo cube;

