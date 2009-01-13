#
# a matrix with holes
#
algebraic3d

solid holes = sphere (0.3, 0.4, 0.4; 0.1)
           or sphere (0.7, 0.2, 0.8; 0.15)
           or sphere (0.8, 0.5, 0.4; 0.11)
           or sphere (0.6, 0.2, 0.8; 0.13)
           or sphere (0.4, 0.3, 0.6; 0.14)
           or sphere (0.6, 0.3, 0.4; 0.16)
           or sphere (0.2, 0.8, 0.6; 0.17)
           or sphere (0.4, 0.6, 0.5; 0.2);

solid cube = plane (0, 0, 0; 0, 0, -1)
         and plane (0, 0, 0; 0, -1, 0)
         and plane (0, 0, 0; -1, 0, 0)
         and plane (1, 1, 1; 0, 0, 1)
         and plane (1, 1, 1; 0, 1, 0)
         and plane (1, 1, 1; 1, 0, 0);

solid rest = cube and not holes;

# two sub-domains

tlo holes -col=[1,0,0];
tlo rest -col=[0,0,1] -transparent;
