from netgen.geom2d import unit_square
from ngsolve import *
from math import pi
mesh = Mesh (unit_square.GenerateMesh(maxh=0.05))

# from netgen.meshing import *
# from netgen.csg import *


# # generate a 1D mesh

# m = Mesh()
# m.dim = 1

# nel = 10

# pnums = []
# for i in range(0, nel+1):
    # pnums.append (m.Add (MeshPoint (Pnt(i/nel, 0, 0))))

# for i in range(0,nel):
    # m.Add (Element1D ([pnums[i],pnums[i+1]], index=1))

# m.Add (Element0D (pnums[0], index=1))
# m.Add (Element0D (pnums[nel], index=2))
# # m.Save("test.vol")




# from ngsolve import *
# mesh = Mesh(m)




order=4
fes = L2(mesh, order=order, dgjumps = True)
fes2 = H1(mesh, order=order)

u = fes.TrialFunction()
v = fes.TestFunction()

n = specialcf.normal(mesh.dim)
h = specialcf.mesh_size
ugf = GridFunction(fes)
uex = GridFunction(fes)
uex.Set(cos(2*pi*x)*sin(2*pi*y))


a = BilinearForm(fes)
cf1 = -0.5 * InnerProduct(grad(u    ), n)*(v-v.Other(bnd=0))
cf2 = -0.5 * InnerProduct(grad(v), n)*(u-u.Other(bnd=0))
cf3 =    2 * ( (order+1)**2)/h * (u-u.Other(bnd=0)) * v

a += grad(u) * grad(v) * dx
a += (cf1+cf2+cf3)*dx(element_boundary=True)
a += (-cf1-cf2-cf3)*ds(skeleton=True)
#a += (-cf1-cf2-cf3)*ds("top|bottom", skeleton=True)
#a += h*u*v*dx

l = LinearForm(fes)
l += 8.0*pi*pi*uex*v*dx                   # volume force
#l += 1e-1*uex*v*ds("left",skeleton=True) # dirichlet data
l += n*grad(uex)*v*ds(skeleton=True)      # neumann data

diff = GridFunction(fes)
ddiff = GridFunction(fes)
l.Assemble()
a.Assemble()
ugf.vec[:]= 0
ugf.vec[fes.GetDofNrs(ElementId(50))[0]] = 1
s=fes.FreeDofs()
s.Clear(fes.GetDofNrs(ElementId(50))[0])
invmat = a.mat.Inverse(s)
w = ugf.vec.CreateVector()
a.Apply(ugf.vec,w)
helper = ugf.vec.CreateVector()
helper.data = l.vec - w
ugf.vec.data = invmat*helper

diff.vec.data=ugf.vec-uex.vec
ddiff.Set((grad(ugf)-grad(uex))*(grad(ugf)-grad(uex)))
err = Integrate(diff*diff,mesh)
err = sqrt(err)
derr = Integrate(ddiff,mesh)
derr = sqrt(derr)
#print("Error = ",err)
#print("D(Error) = ",derr)

Draw (ugf, mesh, "ugf")
Draw(uex,mesh,"uex")
Draw (grad(ugf), mesh, "dugf")
Draw(grad(uex),mesh,"duex")
derr = Integrate(ddiff,mesh)
derr = sqrt(derr)
Draw(diff,mesh,"err")
Draw(ddiff,mesh,"derr")
