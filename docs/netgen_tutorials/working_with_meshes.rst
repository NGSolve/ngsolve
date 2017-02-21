Working with meshes
===================


In this example we create two geometries (a cube and a sphere fitting inside), and mesh them. Then, we manually merge the surface meshes, and create a unified volume mesh, where the sphere and its complement are two different sub-domains.

.. code-block:: python
   
   from netgen.meshing import *
   from netgen.csg import *
   
   # generate brick and mesh it
   geo1 = CSGeometry()
   geo1.Add (OrthoBrick( Pnt(0,0,0), Pnt(1,1,1) ))
   m1 = geo1.GenerateMesh (maxh=0.1)
   m1.Refine()
   
   # generate sphere and mesh it
   geo2 = CSGeometry()
   geo2.Add (Sphere (Pnt(0.5,0.5,0.5), 0.1))
   m2 = geo2.GenerateMesh (maxh=0.05)
   m2.Refine()
   m2.Refine()
   
   print ("***************************")
   print ("** merging suface meshes **")
   print ("***************************")
   
   # create an empty mesh
   mesh = Mesh()
   
   # a face-descriptor stores properties associated with a set of surface elements
   # bc .. boundary condition marker,
   # domin/domout .. domain-number in front/back of surface elements (0 = void)
   
   fd_outside = mesh.Add (FaceDescriptor(bc=1,domin=1))
   fd_inside = mesh.Add (FaceDescriptor(bc=2,domin=2,domout=1))
   
   # copy all boundary points from first mesh to new mesh.
   # pmap1 maps point-numbers from old to new mesh
   
   pmap1 = { }
   for e in m1.Elements2D():
   for v in e.vertices:
   if (v not in pmap1):
   pmap1[v] = mesh.Add (m1[v])
   
   
   # copy surface elements from first mesh to new mesh
   # we have to map point-numbers:
   
   for e in m1.Elements2D():
   mesh.Add (Element2D (fd_outside, [pmap1[v] for v in e.vertices]))
   
                
   
   # same for the second mesh:
   
   pmap2 = { }
   for e in m2.Elements2D():
   for v in e.vertices:
   if (v not in pmap2):
   pmap2[v] = mesh.Add (m2[v])
   
   for e in m2.Elements2D():
   mesh.Add (Element2D (fd_inside, [pmap2[v] for v in e.vertices]))
   
   
   print ("******************")
   print ("** merging done **")
   print ("******************")
   
   
   mesh.GenerateVolumeMesh()
   mesh.Save ("newmesh.vol")


   
