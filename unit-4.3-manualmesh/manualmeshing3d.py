# # Working with meshes
# ## Merge three-dimensional meshes
# 
# In the following example we will merge two surface meshes and generate a unified volume mesh.

# In[20]:


from netgen.meshing import *
from netgen.csg import *

from ngsolve import ngsglobals
ngsglobals.msg_level = 2


# As starting point we create two geometries and mesh them.

# In[21]:


# generate brick and mesh it
geo1 = CSGeometry()
geo1.Add (OrthoBrick( Pnt(0,0,0), Pnt(1,1,1) ))
m1 = geo1.GenerateMesh (maxh=0.1)
# m1.Refine()

# generate sphere and mesh it
geo2 = CSGeometry()
geo2.Add (Sphere (Pnt(0.5,0.5,0.5), 0.1))
m2 = geo2.GenerateMesh (maxh=0.05)
m2.Refine()
# m2.Refine()


# Now we start the merging process. Therefore we create an empty mesh and add a `FaceDescriptor` for each of the surfaces.

# In[22]:


# create an empty mesh
ngmesh = Mesh()

# a face-descriptor stores properties associated with a set of surface elements
# bc .. boundary condition marker,
# domin/domout .. domain-number in front/back of surface elements (0 = void),
# surfnr .. number of the surface described by the face-descriptor

fd_outside = ngmesh.Add (FaceDescriptor(bc=1,domin=1,surfnr=1))
fd_inside = ngmesh.Add (FaceDescriptor(bc=2,domin=2,domout=1,surfnr=2))


# Since the surface elements stay the same in the merged mesh, we copy the points on the surface an the surface elments to the new mesh.

# In[23]:


# copy all boundary points from first mesh to new mesh.
# pmap1 maps point-numbers from old to new mesh
pmap1 = { }
for e in m1.Elements2D():
    for v in e.vertices:
        if (v not in pmap1):
            pmap1[v] = ngmesh.Add (m1[v])


# copy surface elements from first mesh to new mesh
# we have to map point-numbers:

for e in m1.Elements2D():
    ngmesh.Add (Element2D (fd_outside, [pmap1[v] for v in e.vertices]))

# same for the second mesh:
pmap2 = { }
for e in m2.Elements2D():
    for v in e.vertices:
        if (v not in pmap2):
            pmap2[v] = ngmesh.Add (m2[v])

for e in m2.Elements2D():
    ngmesh.Add (Element2D (fd_inside, [pmap2[v] for v in e.vertices]))


# Finally we have to generate the new volume mesh. 

# In[24]:


ngmesh.GenerateVolumeMesh()
import ngsolve
mesh = ngsolve.Mesh(ngmesh)
Draw(mesh)

