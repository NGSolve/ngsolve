#
## A cone
#
algebraic3d

# Cone given by bottom circle and top circle
# and cut by planes:

solid cutcone = cone ( 0, 0, 0; 1; 3, 0, 0; 0.1 )
	and plane (0, 0, 0; -1, 0, 0)
	and plane (3, 0, 0; 1, 0, 0);

tlo cutcone;
