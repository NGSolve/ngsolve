algebraic3d

curve2d procurve1=(8;
	-1,0;
	-0.7,0.7;
	0,1;
	0.7,0.7;
	1,0;
	0.7,-0.7;
	0,-1;
	-0.7,-0.7;
	4;
	3,1,2,3;
	3,3,4,5;
	3,5,6,7;
	3,7,8,1);

curve2d procurve2=(4;
	1,1;
	1,-1;
	-1,-1;
	-1,1;
	4;
	2,1,2;
	2,2,3;
	2,3,4;
	2,4,1);



curve3d pathcurve1=(9;
	0,0,0;
	10,0,5;
	10,10,10;
	10,20,15;
	0,20,20;
	-10,20,25;
	-10,10,30;
	-10,0,35;
	0,0,40;
	4;
	3,1,2,3;
	3,3,4,5;
	3,5,6,7;
	3,7,8,9);

curve3d pathcurve2=(2;
	0,0,0;
	0,10,0;
	1;
	2,1,2);


curve3d pathcurve3=(3;
	0,0,0;
	10,0,5;
	10,10,10;
	1;
	3,1,2,3);

curve3d pathcurve4=(9;
	0,0,0;
	10,0,0;
	10,10,0;
	10,20,0;
	0,20,0;
	-10,20,0;
	-10,10,0;
	-10,0,0;
	0,0,0;
	4;
	3,1,2,3;
	3,3,4,5;
	3,5,6,7;
	3,7,8,9);


solid p1 = plane(1,0,0;-1,0,0);
solid p2 = plane(10,9,10;0,1,0);
solid p3 = plane(0,1,0;0,-1,0);
solid p4 = plane(0,9,0;0,1,0);

solid ob1 = orthobrick(-1,-5,-5;1,5,45);

solid ext = extrusion(pathcurve1;procurve2;0,0,1) and not ob1;

#solid ext = extrusion(pathcurve4;procurve2;0,0,1);

#solid ext = extrusion(pathcurve3;procurve1;0,0,1) and p1 and p2;

#solid ext = extrusion(pathcurve2;procurve2;0,0,1) and p3 and p4;

solid sp = sphere(0,0,0;4);

solid comb = sp or ext;

#tlo ext;

tlo comb;