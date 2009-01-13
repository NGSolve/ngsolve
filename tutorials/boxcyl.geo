#
## Two cylinders on a box
#
algebraic3d

#define box:

solid box = plane (0, 0, 0.5; -1, 0, 0)
	and plane (0, 0, 0.5; 0, -1, 0)
	and plane (0, 0, 0.5; 0, 0, -1)
	and plane (2, 1.5, 1; 1, 0, 0)
	and plane (2, 1.5, 1; 0, 1, 0)
	and plane (2, 1.5, 1; 0, 0, 1);

#define cylinders:

solid cyls = (cylinder (0.5, 0.75, 0; 0.5, 0.75, 2; 0.3) 
	or cylinder (1.5, 0.75, 0; 1.5, 0.75, 2; 0.3) )
	and plane (0, 0, 0.7; 0, 0, -1)
	and plane (0, 0, 1.5; 0, 0, 1);

#combine both:

solid main = box or cyls;


#define sub-domains:
tlo main;

singular edge box cyls;

	
