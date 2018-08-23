# # Geometries in 2D
# 
# We have to import the `SplineGeometry` class from the `geom2d` module to be able to generate two-dimensional geomentries. After importing the module we can create a new `SplineGeometry`.

# In[1]:


from netgen.geom2d import SplineGeometry
geo = SplineGeometry()


# Now we can use one the prededined objects (Rectangle,Circle) or generate our own geometry with lines or rational splines of 2nd order.
# 
# ## Predefined geometries
# First we use the predefined ones and add a rectangle and a circle to our geomentry with the boundary conditions `rectangle` and `circle`.

# In[2]:


geo.AddRectangle((-1,-1),(1,1),bc="rectangle")
geo.AddCircle((0,0),0.5,bc="circle")


# In[3]:


ngmesh = geo.GenerateMesh()


# To get a proper geometry we have to set domain numbers for the domain on the left side of the curve and for the one in its right side. In this case the curves are parametrized in a mathematical positve sense.
# Additionally we can use `SetMaterial` to identify the domains with names.

# In[4]:


geo = SplineGeometry()
geo.AddRectangle(p1=(-1,-1),
                 p2=(1,1),
                 bc="rectangle",
                 leftdomain=1,
                 rightdomain=0)
geo.AddCircle(c=(0,0),
              r=0.5,
              bc="circle",
              leftdomain=2,
              rightdomain=1)
geo.SetMaterial (1, "outer")
geo.SetMaterial (2, "inner")


# In[5]:


geo.SetDomainMaxH(2, 0.02)
ngmesh = geo.GenerateMesh(maxh=0.1)

from ngsolve import *
mesh = Mesh(ngmesh)
Draw(mesh)
