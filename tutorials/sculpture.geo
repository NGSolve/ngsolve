algebraic3d
#
# intersection of sphere and cylinders
# motivated by a sculpture found in St. Gallen
#

solid cyls = cylinder ( -100, 0, 0; 200, 0, 0; 40 )
          or cylinder ( 100, -100, 100; 100, 200, 100; 40)
          or cylinder ( 0, 100, -100; 0, 100, 200; 40);
solid sculpture = sphere (50, 50, 50; 80) and not cyls
        and not sphere (50, 50, 50; 50);

tlo sculpture -col=[0.5, 0.5, 0.5]; 
