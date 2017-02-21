Constructive Solid Geometry CSG
===============================

The constructive solid geometry format allows to define geometric primitives such as spheres and cylinders and perform to boolean operations on them. Such objects are of type Solid.

A cube an be defined and meshed as follows:

.. code-block:: python

   from netgen.csg import *

   left  = Plane (Pnt(0,0,0), Vec(-1,0,0) )
   right = Plane (Pnt(1,1,1), Vec( 1,0,0) )
   front = Plane (Pnt(0,0,0), Vec(0,-1,0) )
   back  = Plane (Pnt(1,1,1), Vec(0, 1,0) )
   bot   = Plane (Pnt(0,0,0), Vec(0,0,-1) )
   top   = Plane (Pnt(1,1,1), Vec(0,0, 1) )
   
   cube = left * right * front * back * bot * top
   geo = CSGeometry()
   geo.Add (cube)
   
   mesh = geo.GenerateMesh(maxh=0.1)
   mesh.Save("cube.vol")

The Plane primitive defines the half-space behind the plane given by an arbitrary point in the plane, and the outgoing normal vector.

Boolean operations are defined as follows:

+-----------+-----------------------------+
| operator  |   set operation             |
+===========+=============================+
|     \*    |   intersection              |
+-----------+-----------------------------+
|     \+    |         union               |
+-----------+-----------------------------+
|     \-    | intersection with complement|
+-----------+-----------------------------+

Available primitives are

==============    =============================== ========================================================================
 primitive         csg syntax 	                   meaning
==============    =============================== ========================================================================
half-space        Plane(Pnt p, Vec n)             point p in plane, normal vector
sphere 	          Sphere(Pnt c,float r)           sphere with center c and radius r
cylinder          Cylinder(Pnt a, Pnt b, float r) points a and b define the axes of a infinite cylinder of radius r
brick 	          OrthoBrick ( Pnt a, Pnt b ) 	  axes parallel brick with minimal coordinates a and maximal coordinates b
==============    =============================== ========================================================================


Using the orthobrick primitive, the cube above can be defined by one statement. Drilling a hole through the cube is defined as follows:

.. code-block:: python

   from netgen.csg import *

   cube = OrthoBrick( Pnt(0,0,0), Pnt(1,1,1) )
   hole = Cylinder ( Pnt(0.5, 0.5, 0), Pnt(0.5, 0.5, 1), 0.2)
   
   geo = CSGeometry()
   geo.Add (cube-hole)
   mesh = geo.GenerateMesh(maxh=0.1)
   mesh.Save("cube_hole.vol")

If the whole domain consists of several regions, we give several solid to the geometry object. Now, the first domain is the cube with the hole cut out, the second region is the hole. Don't forget that the cylinder is an infinite cylinder, and must be cut to finite length:

.. code-block:: python

   geo = CSGeometry()
   geo.Add (cube-hole)
   geo.Add (cube*hole)

