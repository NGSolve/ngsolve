from ngsolve import *
from netgen.geom2d import unit_square
from math import pi
from math import log

def test_detect_symmetric():
  deg = 3
  mesh = Mesh(unit_square.GenerateMesh(maxh=0.4))
  def dxx(u): return hess(u)[0,0]
  def dxy(u): return hess(u)[0,1]
  def dyx(u): return hess(u)[1,0]
  def dyy(u): return hess(u)[1,1]
  def dxxavg(u): return 0.5*(dxx(u)+dxx(u.Other()))
  def dxyavg(u): return 0.5*(dxy(u)+dxy(u.Other()))
  def dyxavg(u): return 0.5*(dyx(u)+dyx(u.Other()))
  def dyyavg(u): return 0.5*(dyy(u)+dyy(u.Other()))
  def jumpdn(u): return (grad(u)-grad(u).Other())*n
  def hess(u): return u.Operator("hesse")
  def Lap(u): return dxx(u)+dyy(u)
  def signfct(a): return a/sqrt(a*a)
  d2udxx = -pi*pi*sin(pi*x)*sin(pi*y)
  d2udxy = pi*pi*cos(pi*x)*cos(pi*y)
  d2udyy = -pi*pi*sin(pi*x)*sin(pi*y)
  A00 = 2.0
  A01 = signfct((x-0.5)*(y-0.5))
  A11 = 2.0
  fsource = A00*d2udxx+2*A01*d2udxy+A11*d2udyy
  fes = H1(mesh, order=deg, dirichlet=[1,2,3,4],autoupdate=True,dgjumps=True)
  n = specialcf.normal(2)
  h = specialcf.mesh_size
  mu = 100/h
  def J_h_inner(u,v): return mu*(jumpdn(u)*jumpdn(v)) * dx(bonus_intorder = 10,skeleton=True)
  def B_star_1(u,v):
    Bstar1 = (dxx(u)*dxx(v)+dyx(u)*dyx(v)+dxy(u)*dxy(v)+dyy(u)*dyy(v))*dx(bonus_intorder = 10)
    return Bstar1
  def B_star_2(u,v):
    Bstar2 = (dxxavg(u)+dyyavg(u)\
    -dxxavg(u)*n[0]*n[0]\
    -dxyavg(u)*n[0]*n[1]\
    -dyxavg(u)*n[1]*n[0]\
    -dyyavg(u)*n[1]*n[1])\
    * \
    jumpdn(v)\
    * dx(skeleton=True)\
    +jumpdn(u)*\
    (dxxavg(v)+dyyavg(v)\
    -dxxavg(v)*n[0]*n[0]\
    -dxyavg(v)*n[0]*n[1]\
    -dyxavg(v)*n[1]*n[0]\
    -dyyavg(v)*n[1]*n[1])\
    * dx(skeleton=True)
    return Bstar2
  def B_star(u,v):
    Bstar = B_star_1(u,v)+B_star_2(u,v)
    return Bstar
  def Deltainner(u,v):
      delt = Lap(u)*Lap(v)*dx(bonus_intorder = 10)
      return delt
  def B(u,v):
    B = 0.5*B_star(u,v)+0.5*Deltainner(u,v)
    return B
  #renormalisation parameter
  gamma = (A00+A11)/(A00**2+2*A01**2+A11**2)

  u = fes.TrialFunction()
  v = fes.TestFunction()

  # the right hand side
  f = LinearForm(fes)
  f += fsource * gamma * Lap(v) * dx
  # the bilinear-forms with orders changed
  a1 = BilinearForm(fes, symmetric=False)
  a1 +=   B(u,v) - Deltainner(u,v) + gamma*(A00*dxx(u)+2*A01*dxy(u)+A11*dyy(u))*Lap(v)*dx + J_h_inner(u,v)
  a2 = BilinearForm(fes, symmetric=False)
  a2 +=  gamma*(A00*dxx(u)+2*A01*dxy(u)+A11*dyy(u))*Lap(v)*dx + B(u,v) - Deltainner(u,v) + J_h_inner(u,v)
  a3 = BilinearForm(fes, symmetric=False)
  a1.Assemble()
  a2.Assemble()
  gfu1 = GridFunction(fes,autoupdate=True)
  gfu2 = GridFunction(fes,autoupdate=True)
  f.Assemble()
  gfu1.vec.data = a1.mat.Inverse(fes.FreeDofs())*f.vec
  gfu2.vec.data = a2.mat.Inverse(fes.FreeDofs())*f.vec
  #### value should be zero but it is not
  L2err = sqrt(Integrate ( (gfu1-gfu2)*(gfu1-gfu2), mesh))
  print(L2err)
  assert L2err < 1e-14
