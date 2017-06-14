# # Working with meshes
# ## Two-dimensional meshes
# As example we mesh a unit square [0,1]x[0,1] using quadrilaterals.
# 

# In[7]:

from netgen.meshing import *
from netgen.geom2d import unit_square


# We create an empty mesh and inititalize the geometry and the dimension

# In[8]:


ngmesh = Mesh()
ngmesh.SetGeometry(unit_square)
ngmesh.dim = 2


# and add all the `MeshPoint`'s we will need for the final mesh. Similar to the one-dimensional mesh we store the `PointId`'s in the `pnums` array.

# In[9]:

N = 5
pnums = []
for i in range(N + 1):
    for j in range(N + 1):
        pnums.append(ngmesh.Add(MeshPoint(Pnt(i / N, j / N, 0))))


# Next, we add the quadrilaterals to the mesh. Before that we have so add a `FaceDescriptor` to the mesh.

# In[10]:


ngmesh.Add (FaceDescriptor(surfnr=1,domin=1,bc=1))
ngmesh.SetMaterial(1, "mat")
for j in range(N):
    for i in range(N):
        ngmesh.Add(Element2D(1, [pnums[i + j * (N + 1)],
                                 pnums[i + (j + 1) * (N + 1)],
                                 pnums[i + 1 + (j + 1) * (N + 1)],
                                 pnums[i + 1 + j * (N + 1)]]))


# Finally we have to add boundary elements and set boundary conditions.

# In[11]:


# horizontal boundaries
for i in range(N):
   ngmesh.Add(Element1D([pnums[N + i * (N + 1)], pnums[N + (i + 1) * (N + 1)]], index=1))
   ngmesh.Add(Element1D([pnums[0 + i * (N + 1)], pnums[0 + (i + 1) * (N + 1)]], index=1))

# vertical boundaries
for i in range(N):
   ngmesh.Add(Element1D([pnums[i], pnums[i + 1]], index=2))
   ngmesh.Add(Element1D([pnums[i + N * (N + 1)], pnums[i + 1 + N * (N + 1)]], index=2))


# In[12]:


from ngsolve.solve import Draw
import ngsolve
mesh = ngsolve.Mesh(ngmesh)
Draw(mesh)
