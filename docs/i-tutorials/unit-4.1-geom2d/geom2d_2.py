# # Geometries in 2D
# 
# We have to import the `SplineGeometry` class from the `geom2d` module to be able to generate two-dimensional geomentries. After importing the module we can create a new `SplineGeometry`.

# In[1]:

from netgen.geom2d import SplineGeometry


# ## Using lines and splines
# We define a new geometry and and write a list of points we want to use for out geometry and add them to geo geometry.

# In[6]:


geo = SplineGeometry()

pnts =[(0,0),
       #(0,0,0.05), # define a local mesh refinement for one point
       (1,0),
       (1,0.5),
       (1,1),
       (0.5,1),
       (0,1)]

p1,p2,p3,p4,p5,p6 = [geo.AppendPoint(*pnt) for pnt in pnts]


# Then we define the curves which define our geometry and add them to the geometry using `Append`.

# In[7]:


curves = [[["line",p1,p2],"bottom"],
          [["line",p2,p3],"right"],
          [["spline3",p3,p4,p5],"curve"],
          [["line",p5,p6],"top"],
          [["line",p6,p1],"left"]]

[geo.Append(c,bc=bc) for c,bc in curves]


# In[8]:


ngmesh = geo.GenerateMesh(maxh=0.2)


# In[9]:


# from ngsolve import *
# mesh = Mesh(ngmesh)
# mesh.Curve(3)


# Additonally to the boundary condition one can set a maximal mesh size for a whole curve with an optional argument `maxh`.

# In[10]:


geo = SplineGeometry()

p1,p2,p3,p4,p5,p6 = [geo.AppendPoint(*pnt) for pnt in pnts]

geo.Append(["line",p1,p2],maxh=0.02)
geo.Append(["line",p2,p4])
geo.Append(["line",p4,p6])
geo.Append(["line",p6,p1])

ngmesh = geo.GenerateMesh(maxh=0.2)#,quad_dominated=True)

from ngsolve import *
mesh = Mesh(ngmesh)
Draw(mesh)
