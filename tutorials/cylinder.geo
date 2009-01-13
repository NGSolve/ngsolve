#
## a cylinder
#
algebraic3d

# cut cylinder by planes:

solid fincyl = cylinder ( 3, 0, 0; -1, 0, 0; 0.5 )
	and plane (0, 0, 0; -1, 0, 0)
	and plane (2, 0, 0; 1, 0, 0);

tlo fincyl;
