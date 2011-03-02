algebraic3d

solid p1 = plane (0.5, 0, 0; 1, 0, 0);

# since surfaces of both bricks are identic they get the same bc id:
solid brick1 = orthobrick (0,0,0; 1,1,1) and p1 -bc=1;
solid brick2 = orthobrick (0,0,-1; 1,1,0) and p1 -bc=2;


tlo brick1;
tlo brick2;

# override bc number: 
# all faces of solid p1 belonging to the boundary of tlo brick1 get bc=3

boundarycondition p1 brick1 3;
