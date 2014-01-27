algebraic3d

solid arperture = 
        plane (0.1,0,0; 1, 0, 0) 
        and plane (-1.1, 0, 0; -1, 0, 0)
        and cylinder (-1.1, 0, 0; 0.1, 0, 0; 0.5);

solid pin = plane (-3, 0, 0; -1, 0, 0) -bc=1;
solid cyl = cylinder (-3, 0, 0; 1, 0, 0; 2)
        and pin
        and plane (-1, 0, 0; 1, 0, 0);


solid interface = sphere (0, 0, 0; 2) -bc=3;

solid free = interface and plane (0, 0, 0; -1, 0, 0) -bc = 2;
solid pml =  sphere (0, 0, 0; 3) and not interface and plane (0, 0, 0; -1, 0, 0) -bc=2;

solid vol = cyl or arperture or free -bc=2;


tlo vol -maxh=0.4;
tlo pml -col=[0,1,0] -transparent -maxh=0.4;
