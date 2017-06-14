# # Working with meshes
# ## One-dimensional meshes
# Meshes in one-dimension can be constructed using the `netgen.meshing` module.
# We just have to add segments (`Element1D`) and the boundaries (`Element0D`).

# In[1]:


from netgen.meshing import *


# First we create a new `Mesh` and set the spatial dimension to 1.

# In[2]:


m = Mesh()
m.dim = 1


# Then we define and add the `MeshPoint`'s to the mesh. The function `m.Add` returns `PoindId`'s which we store in an array to be able to construct the segments in the next step.

# In[3]:


N = 10
pnums = []
for i in range(0, N+1):
    pnums.append (m.Add (MeshPoint (Pnt(i/N, 0, 0))))

type(pnums[0])


# Now we can loop over the array with the `PointId`'s and add one-dimensional elements to the mesh. Further we can set the material for our domain.

# In[4]:


for i in range(0,N):
    m.Add (Element1D ([pnums[i],pnums[i+1]], index=1))

m.SetMaterial(1,'material')


# Finally we have to add the boundary elements and set boundary conditions.

# In[5]:


m.Add (Element0D (pnums[0], index=1))
m.Add (Element0D (pnums[N], index=2))


# To be able to visualize one-dimensional meshes and solution activate `Show edges` in the menu `View > Viewing options > Mesh`.

# In[6]:


import ngsolve
mesh = ngsolve.Mesh(m)
