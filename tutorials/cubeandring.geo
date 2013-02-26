algebraic3d

curve2d testcurve=(8;
	-0.5,1;
	-0.55,1.5;
#	-0.275,1.775;
	-0.5,2;
	0,2.05;
	0.5,2;
	0.55,1.5;
	0.5,1;
	0,0.95;
	4;
	3,8,1,2;
	3,2,3,4;	
	3,4,5,6;
	3,6,7,8);

#curve2d testcurve=(8;
#	-0.5,1;
#	-0.55,1.5;
#	-0.5,2;
#	0,2.05;
#	0.5,2;
#	0.55,1.5;
#	0.5,1;
#	0,0.95;
#	4;
#	3,8,1,2;
#	3,2,3,4;	
#	3,4,5,6;
#	3,6,7,8);

curve2d testcurve1=(4;
	-0.55,1.5;
	0,2.05;
	0.55,1.5;
	0,0.95;
	4;
	2,1,2;
	2,2,3;
	2,3,4;
	2,4,1);


solid mytorus = revolution(0,0,0.5;1,0,0.5;testcurve);
#solid mytorus = revolution(0,0,0.5;1,0,0.5;testcurve1);

solid bbb = orthobrick(-4,-4,-4;4,4,0.1);

solid brickandring = mytorus or bbb;


tlo brickandring;

