algebraic3d
#
## CSG feature copy
#


# define a axis parallel brick:

solid br = orthobrick (0, 0, 0; 20, 20, 1);

# define reference cylinder:

solid cyl1 = cylinder (0.5, 0.5, -1; 0.5, 0.5, 3; 0.2);


# make copies:
solid cylx = multitranslate (1, 0, 0; 19; cyl1);
solid cyls = multitranslate (0, 1, 0; 19; cylx);

solid main = br and not cyls;

tlo main;

# provide bounding-box for fastening bisection alg:

boundingbox (-1, -1, -1; 21, 21, 2);
