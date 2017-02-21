.. _manual-mesh-generation:

Manual mesh generation
======================

In this example we create a mesh by hand, i.e. by prescribing all vertices, edges, faces and elements on our own. This example creates a structured grid for the unit square [0,1]x[0,1] using triangles or quadrilateral (choice via switch):

We first include what we need from netgen:

.. code-block:: python

   from netgen.geom2d import unit_square, MakeCircle, SplineGeometry
   from netgen.meshing import Element0D, Element1D, Element2D, MeshPoint, FaceDescriptor, Mesh
   from netgen.csg import Pnt

Next, we decide on the parameters for the mesh:

.. code-block:: python

   quads = True
   N=5

We create an empty mesh and initialize the geometry and the dimension:

.. code-block:: python

   mesh = Mesh()
   mesh.SetGeometry(unit_square)
   mesh.dim = 2

Then, we add all mesh points that we will need for the final mesh. Note that these MeshPoints are added to the mesh with 'mesh.Add(..)' which return the point index. This index is then stored in the array 'pnums'. In our case we have the simple structure that we will have (N+1)*(N+1) points in total.

.. code-block:: python

   pnums = []
   for i in range(N + 1):
       for j in range(N + 1):
           pnums.append(mesh.Add(MeshPoint(Pnt(i / N, j / N, 0))))

Next, we add the area elements. Between four neighboring points we either span a quadrilateral (if quads==True) or divide the area into two triangle. These are then simply added to the mesh:

.. code-block:: python

   mesh.SetMaterial(1, "mat")
   for j in range(N):
       for i in range(N):
           if quads:
               mesh.Add(Element2D(1, [pnums[i + j * (N + 1)], pnums[i + (j + 1) * (N + 1)], pnums[i + 1 + (j + 1) * (N + 1)], pnums[i + 1 + j * (N + 1)]]))
           else:
               mesh.Add(Element2D(1, [pnums[i + j * (N + 1)], pnums[i + (j + 1) * (N + 1)], pnums[i + 1 + j * (N + 1)]]))
               mesh.Add(Element2D(1, [pnums[i + (j + 1) * (N + 1)], pnums[i + 1 + (j + 1) * (N + 1)], pnums[i + 1 + j * (N + 1)]]))

Now, we have to add boundary elements and boundary conditions. Therefore we add a FaceDescriptor:

.. code-block:: python

   mesh.Add (FaceDescriptor(surfnr=1,domin=1,bc=1))

followed by the horizontal boundary elements

.. code-block:: python

   for i in range(N):
      mesh.Add(Element1D([pnums[N + i * (N + 1)], pnums[N + (i + 1) * (N + 1)]], index=1))
      mesh.Add(Element1D([pnums[0 + i * (N + 1)], pnums[0 + (i + 1) * (N + 1)]], index=1))

and the vertical boundary elements:

.. code-block:: python

   for i in range(N):
      mesh.Add(Element1D([pnums[i], pnums[i + 1]], index=1))
      mesh.Add(Element1D([pnums[i + N * (N + 1)], pnums[i + 1 + N * (N + 1)]], index=1))

This, together results in a valid mesh. Note that we have chosen boundary condition bc=1 on all boundaries.
