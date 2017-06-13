from netgen.geom2d import SplineGeometry
from ngsolve import Mesh, H1, BilinearForm, SymbolicBFI, GridFunction
from ngsolve import x, y, exp, grad, LinearForm, SymbolicLFI, Draw

# Geometry
geo = SplineGeometry()
geo.AddCircle((0.5, 0.5), 0.8,  bc="outer")
geo.AddRectangle((0.7, 0.3), (0.75, 0.7),
                 leftdomain=0, rightdomain=1, bc="scat")
mesh = Mesh(geo.GenerateMesh(maxh=0.05))

# Spaces
fes = H1(mesh, order=4, complex=True)
u = fes.TrialFunction()
v = fes.TestFunction()

# Wavenumber & source
omega = 100
pulse = 1e3*exp(-(100**2)*((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5)))

# Forms
a = BilinearForm(fes)
a += SymbolicBFI(grad(u)*grad(v)-omega**2*u*v)
a += SymbolicBFI(-omega*1j*u*v, definedon=mesh.Boundaries("outer"))
a.Assemble()

f = LinearForm(fes)
f += SymbolicLFI(pulse * v)
f.Assemble()

# Solve
gfu = GridFunction(fes, name="u")
gfu.vec.data = a.mat.Inverse() * f.vec
Draw(gfu)
