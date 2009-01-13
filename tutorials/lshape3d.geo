algebraic3d

#
## 3D Lshape - domain
#

solid c1 = plane (-1, -1, 0; 0, 0, -1)
         and plane (-1, -1, 0; 0, -1, 0)
         and plane (-1, -1, 0; -1, 0, 0)
         and plane (1, 1, 1; 0, 0, 1)
         and plane (1, 1, 1; 0, 1, 0)
         and plane (1, 1, 1; 1, 0, 0);

solid f1 = plane (0, 0, 0; -1, 0, 0);
solid f2 = plane (0, 0, 0; 0, 1, 0);

solid main = c1 and not (f1 and f2);

tlo main;







