algebraic3d

curve2d testcurve=(7;
	1,0;
	1,1;
	1.5,1.5;
	2,3;
	2.5,1.5;
	3,1;
	3,0;
	3;
	3,1,2,3;
	3,3,4,5;
	3,5,6,7);

solid something = revolution(0,0,0;1,0,0;testcurve);

tlo something;