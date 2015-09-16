algebraic3d

# the first plane gets boundary condition nr 1:

solid p1 =  plane (0, 0, 0; -1, 0, 0) -bc=1;

# and the other faces b.c. number = 2:
solid beam = plane (0, 0, 0; 0, 0, -1) 
         and plane (0, 0, 0; 0, -1, 0)
         and p1
         and plane (10, 2, 1; 0, 0, 1)
         and plane (10, 2, 1; 0, 1, 0)
         and plane (10, 2, 1; 1, 0, 0) -bc=2 -maxh=2;

tlo beam;


