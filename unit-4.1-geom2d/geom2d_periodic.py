# # Geometries in 2D
# 
# We have to import the `SplineGeometry` class from the `geom2d` module to be able to generate two-dimensional geomentries. After importing the module we can create a new `SplineGeometry`.

# ## Periodic geometries
# The following example showa how construct a geometry for a periodic $L_2$ finite element space.
# 
# Again we start with adding the points to the geometry. In this case the points of a hexgon . For the first three segments we save the return value (the line number) of `geo.Append`. Now we can use those line numbers to identify  each of the last three segments with the opposite, already added one. This identification is done with the optional argument `copy`. The meshing algorithm then just copies the boundary mesh to the opposite segment. Thus the segments have to have the same orientation.

# In[11]:


from math import pi, cos, sin
geo = SplineGeometry()
pnums = [ geo.AddPoint(cos(phi),sin(phi)) for phi in [x*pi/3 for x in range(6)] ]
l1 = geo.Append(["line", 0, 1], leftdomain=1, rightdomain=0, bc="upperRight")
l2 = geo.Append(["line", 1, 2], leftdomain=1, rightdomain=0, bc="upperCenter")
l3 = geo.Append(["line", 2, 3], leftdomain=1, rightdomain=0, bc="upperLeft")
geo.Append(["line", 0, 5], leftdomain=0, rightdomain=1, bc="lowerRight", copy = l3)
geo.Append(["line", 5, 4], leftdomain=0, rightdomain=1, bc="lowerCenter", copy = l2)
geo.Append(["line", 4, 3], leftdomain=0, rightdomain=1, bc="lowerLeft", copy = l1)
ngmesh = geo.GenerateMesh(maxh=0.1)

