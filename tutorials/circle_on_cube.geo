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


solid top = plane (1,1,1; 0, 0, 1);
solid cyl = top 
	and plane (0,0,0; 0, 0, -1)
	and cylinder (0.5, 0.5, 0; 0.5, 0.5, 1; 0.2);
	

tlo cube;

# take just surface 'top' of solid 'cyl'
tlo cyl top -col=[1,0,0];


