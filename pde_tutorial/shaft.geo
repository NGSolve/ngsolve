#
## Crankshaft
#
algebraic3d
solid p1 = plane (0, 0, 0; -1, 0, 0)
         and plane (10, 0, 0; 1, 0, 0)
         and
         (     plane (35, 0, 28; 0, -1, 3)
           and plane (35, 0, -28; 0, -1, -3)
           and plane (35, 0, 0; 0, 1, 0)
           and plane (35, -30, 0; 0, -1, 0)
           or cylinder (-10, 0, 0; 20, 0, 0; 30)
           or cylinder (-10, -30, 0; 20, -30, 0; 20)
         );
solid p2 = plane (35, 0, 0; -1, 0, 0)
         and plane (45, 0, 0; 1, 0, 0)
         and
         (     plane (35, 0, 28; 0, -1, 3)
           and plane (35, 0, -28; 0, -1, -3)
           and plane (35, 0, 0; 0, 1, 0)
           and plane (35, -30, 0; 0, -1, 0)
           or cylinder (30, 0, 0; 50, 0, 0; 30)
           or cylinder (30, -30, 0; 50, -30, 0; 20)
         );
solid p3 = plane (80, 0, 0; -1, 0, 0)
         and plane (90, 0, 0; 1, 0, 0)
         and
         (     plane (0, 0, 28; 0, 1, 3)
           and plane (0, 0, -28; 0, 1, -3)
           and plane (0, 0, 0; 0, -1, 0)
           and plane (0, 30, 0; 0, 1, 0)
           or cylinder (70, 0, 0; 100, 0, 0; 30)
           or cylinder (70, 30, 0; 100, 30, 0; 20)
         );
solid p4 = plane (115, 0, 0; -1, 0, 0)
         and plane (125, 0, 0; 1, 0, 0)
         and
         (     plane (35, 0, 28; 0, 1, 3)
           and plane (35, 0, -28; 0, 1, -3)
           and plane (35, 0, 0; 0, -1, 0)
           and plane (35, 30, 0; 0, 1, 0)
           or cylinder (110, 0, 0; 130, 0, 0; 30)
           or cylinder (110, 30, 0;130, 30, 0; 20)
         );

# bc for main bearing:
solid cyl1a = cylinder (-50, 0, 0; 200, 0, 0; 15) -bc=1;
solid cyl1b = cylinder (-50, 0, 0; 200, 0, 0; 15.0001) -bc=2;
solid cyl1c = cylinder (-50, 0, 0; 200, 0, 0; 15.0002) -bc=3;

solid sh1 = cyl1a
            and plane (-40, 0, 0; -1, 0, 0)
            and plane (5, 0, 0; 1, 0, 0) -bc=5;
solid sh2 = cyl1b
            and plane (40, 0, 0; -1, 0, 0)
            and plane (85, 0, 0; 1, 0, 0);
solid sh3 = cyl1c
            and plane (120, 0, 0; -1, 0, 0)
            and plane (165, 0, 0; 1, 0, 0);

# bc for load:
solid pl1 = cylinder (0, -30, 0; 50, -30, 0; 10)
            and plane (5, 0, 0; -1, 0, 0)
            and plane (40, 0, 0; 1, 0, 0);

solid pl2 = cylinder (80, 30, 0; 130, 30, 0; 10)
            and plane (85, 0, 0; -1, 0, 0)
            and plane (120, 0, 0; 1, 0, 0) -bc=4;
#
#
solid main = p1 or p2 or p3 or p4 or sh1 or sh2 or sh3 or pl1 or pl2 -bc=6;

tlo main;
