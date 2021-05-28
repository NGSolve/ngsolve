#!/usr/bin/env python
# coding: utf-8

import numpy as np
from ngsolve import *
from ngsolve.comp import IntegrationRuleSpace
from ngsolve.fem import MinimizationCF, NewtonCF
from netgen.csg import *
import netgen


# 'minimal' mesh
m = netgen.meshing.Mesh(dim=1)

N = 1
pnums = []
for i in range(0, N+1):
    pnums.append (m.Add (netgen.meshing.MeshPoint (Pnt(2*i/N, 0, 0))))
    
idx = m.AddRegion("material", dim=1)
for i in range(0,N):
    m.Add (netgen.meshing.Element1D ([pnums[i],pnums[i+1]], index=idx))
    
idx_left = m.AddRegion("left", dim=0)
idx_right = m.AddRegion("right", dim=0)

m.Add (netgen.meshing.Element0D (pnums[0], index=idx_left))
m.Add (netgen.meshing.Element0D (pnums[N], index=idx_right))

mesh = Mesh(m)

d = mesh.dim

int_order = 0

fes_ir = IntegrationRuleSpace(mesh, order=int_order)
fes_ir_vec = fes_ir**2
fes_compound_a = fes_ir * fes_ir
fes_M = MatrixValued(fes_ir, dim=3, symmetric=True)
fes_compound_b = fes_ir * fes_M

irs = fes_ir.GetIntegrationRules()
irs_dx = dx(intrules=irs)


# In[ ]:


irs # missing for line elements


# # Simplest case: a equation in a single scalar variable

# In[ ]:


u = GridFunction(fes_ir)
du = fes_ir.TrialFunction()


# ## Linear equation

# In[ ]:


pot = du**2 - du/2
eq = 2*du - 1/2

expected = np.array([1/4])


# In[ ]:


u.Interpolate(CoefficientFunction(3))


# In[ ]:


ncf = NewtonCF(eq, u)


# In[ ]:


u.Interpolate(ncf)


# In[ ]:


u.vec.FV().NumPy()


# In[ ]:


assert np.allclose(u.vec.FV().NumPy(), expected, atol=1e-14, rtol=0)


# In[ ]:





# In[ ]:


u.Interpolate(CoefficientFunction(3))


# In[ ]:


mcf = MinimizationCF(pot, u)


# In[ ]:


u.Interpolate(mcf)


# In[ ]:


u.vec.FV().NumPy()


# In[ ]:


assert np.allclose(u.vec.FV().NumPy(), expected, atol=1e-14, rtol=0)


# ## Nonlinear equation

# In[ ]:


pot = du**3 - du**2/2
eq = 3*du**2 - du

expected = np.array([1/3])


# In[ ]:


u.Interpolate(CoefficientFunction(3))


# In[ ]:


ncf = NewtonCF(eq, u)


# In[ ]:


u.Interpolate(ncf)


# In[ ]:


u.vec.FV().NumPy()


# In[ ]:


assert np.allclose(u.vec.FV().NumPy(), expected, atol=1e-14, rtol=0)


# In[ ]:





# In[ ]:


u.Interpolate(CoefficientFunction(3))


# In[ ]:


mcf = MinimizationCF(pot, u)


# In[ ]:


u.Interpolate(mcf)


# In[ ]:


u.vec.FV().NumPy()


# In[ ]:


assert np.allclose(u.vec.FV().NumPy(), expected, atol=1e-14, rtol=0)


# # Simplest compound case: two linear equations in two scalar variables

# In[ ]:


u = GridFunction(fes_compound_a)
du1, du2 = fes_compound_a.TrialFunction()

uvec = GridFunction(fes_ir_vec)
duvec = fes_ir_vec.TrialFunction()


# In[ ]:


u.Interpolate(CoefficientFunction((3, 3)))
uvec.Interpolate(CoefficientFunction((3, 3)))


# In[ ]:


u.vec.FV().NumPy(), uvec.vec.FV().NumPy()


# ## Linear equation with SPD Jacobian

# In[ ]:


du = CoefficientFunction((du1, du2))
a = CoefficientFunction((2/3, 1))
M = CoefficientFunction((2, 1/2, 1/2, 4), dims=(2, 2))

def pot_func(u):
    return 1/2 * InnerProduct(M, OuterProduct(u, u)) + 4 * InnerProduct(u, a)

def res_func(u):
    return M * u + 4 * a

eq = res_func(du)
eq_vec = res_func(duvec)
check = res_func(uvec)
pot = pot_func(duvec)

expected = np.array([-1.11827957, -0.86021505])


# ### Compound space version

# In[ ]:


u.Interpolate(CoefficientFunction((3, 3)))
uvec.Interpolate(CoefficientFunction((3, 3)))
ncf = NewtonCF(eq, u, maxiter=1)
uvec.Interpolate(ncf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# In[ ]:


u.Interpolate(CoefficientFunction((3, 3)))
uvec.Interpolate(CoefficientFunction((3, 3)))
ncf = NewtonCF(eq, u.components, maxiter=1)
uvec.Interpolate(ncf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# ### Vector space version

# In[ ]:


uvec.Interpolate(CoefficientFunction((3, 3)))
ncf = NewtonCF(eq_vec, uvec, maxiter=1)
uvec.Interpolate(ncf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# ### Minimization version

# In[ ]:


uvec.Interpolate(CoefficientFunction((3, 3)))
mcf = MinimizationCF(pot, uvec)
uvec.Interpolate(mcf, maxiter=1)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# ### Check solution

# In[ ]:


uvec2 = GridFunction(uvec.space)
uvec2.Interpolate(check)
assert np.allclose(uvec2.vec.FV().NumPy(), 0)


# ## Linear equation with non-symmetric Jacobian

# In[ ]:


du = CoefficientFunction((du1, du2))
a = CoefficientFunction((2/3, 1))
M = CoefficientFunction((2, 1/2, 1, 4), dims=(2, 2))

def res_func(u):
    return M * u + 4 * a

eq = res_func(du)
eq_vec = res_func(duvec)
check = res_func(uvec)

expected = np.array([-1.15555556, -0.71111111])


# ### Compound space version

# In[ ]:


u.Interpolate(CoefficientFunction((3, 3)))
uvec.Interpolate(CoefficientFunction((3, 3)))
ncf = NewtonCF(eq, u, maxiter=1)
uvec.Interpolate(ncf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# In[ ]:


u.Interpolate(CoefficientFunction((3, 3)))
uvec.Interpolate(CoefficientFunction((3, 3)))
ncf = NewtonCF(eq, u.components, maxiter=1)
uvec.Interpolate(ncf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# ### Vector space version

# In[ ]:


uvec.Interpolate(CoefficientFunction((3, 3)))
ncf = NewtonCF(eq_vec, uvec, maxiter=1)
uvec.Interpolate(ncf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# ### Check solution

# In[ ]:


uvec2 = GridFunction(uvec.space)
uvec2.Interpolate(check)
assert np.allclose(uvec2.vec.FV().NumPy(), 0)


# ## Nonlinear equation with nonsymmetric Jacobian

# In[ ]:


du = CoefficientFunction((du1, du2))
a = CoefficientFunction((2/3, 1))
M = CoefficientFunction((2, 1/2, 1, 4), dims=(2, 2))

def res_func(u):
    return InnerProduct(M, OuterProduct(u, u)) * (u + a) + 4 * a

eq = res_func(du)
eq_vec = res_func(duvec)
check = res_func(uvec)

expected = np.array([-0.90980601, -1.36470902])


# ### Compound space version

# In[ ]:


u.Interpolate(CoefficientFunction((-1, -1)))
ncf = NewtonCF(eq, u)
uvec.Interpolate(ncf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# ### Vector space version

# In[ ]:


uvec.Interpolate(CoefficientFunction((-1, -1)))
ncf = NewtonCF(eq_vec, uvec)
uvec.Interpolate(ncf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# ### Check solution

# In[ ]:


uvec2 = GridFunction(uvec.space)
uvec2.Interpolate(check)
assert np.allclose(uvec2.vec.FV().NumPy(), 0)


# # Compound case: one scalar and one vector-valued variable

# In[ ]:


fes_compound = fes_ir * (fes_ir**2)
u = GridFunction(fes_compound)
u1, u2 = u.components
du1, du2 = fes_compound.TrialFunction()

# for comparison
fes_ir_vec = fes_ir ** 3
uvec = GridFunction(fes_ir_vec)
duvec = fes_ir_vec.TrialFunction()


# In[ ]:


# u.Interpolate(CoefficientFunction((3, 3, 3))) # fails
cf1 = CoefficientFunction(1)
cf2 = CoefficientFunction((2, 3))
cf = CoefficientFunction((1, 2, 3))
u1.Interpolate(cf1)
u2.Interpolate(cf2)
uvec.Interpolate(cf)


# In[ ]:


u.vec.FV().NumPy(), uvec.vec.FV().NumPy()


# ## Linear equation with symmetric Jacobian

# In[ ]:


du = CoefficientFunction((du1, du2))
a = CoefficientFunction((2/3, 1, -1))
M = CoefficientFunction((2,    1/2, 1/3, 
                         1/2,    4,   1/6,
                         1/3,  1/6,   3), dims=(3, 3))

def pot_func(u):
    return 1/2 * InnerProduct(u, M * u) + 4 * InnerProduct(a, u)

def res_func(u):
    return M * u + 4 * a

pot = pot_func(du)
eq = res_func(du)
eq_vec = res_func(duvec)
check = res_func(uvec)

expected = np.array([-1.36581405, -0.89321965,  1.53471376])


# ### Compound space version

# In[ ]:


u1.Interpolate(cf1)
u2.Interpolate(cf2)
uvec.Interpolate(cf)
ncf = NewtonCF(eq, u, maxiter=1)
uvec.Interpolate(ncf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# In[ ]:


u1.Interpolate(cf1)
u2.Interpolate(cf2)
uvec.Interpolate(cf)
mcf = MinimizationCF(pot, u, maxiter=1)
uvec.Interpolate(mcf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# ### Vector space version

# In[ ]:


u1.Interpolate(cf1)
u2.Interpolate(cf2)
uvec.Interpolate(cf)
ncf = NewtonCF(eq_vec, uvec, maxiter=1)
uvec.Interpolate(ncf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# ### Check solution

# In[ ]:


uvec2 = GridFunction(uvec.space)
uvec2.Interpolate(check)
assert np.allclose(uvec2.vec.FV().NumPy(), 0)


# ## Linear equation with nonsymmetric Jacobian

# In[ ]:


du = CoefficientFunction((du1, du2))
a = CoefficientFunction((2/3, 1, -1))
M = CoefficientFunction((2,    1/2, 1/3, 
                         1/6,    4,   0,
                         1/4,  0,   2), dims=(3, 3))

def res_func(u):
    return M * u + 4 * a

eq = res_func(du)
eq_vec = res_func(duvec)
check = res_func(uvec)

expected = np.array([-1.46236559, -0.9390681 ,  2.1827957])


# ### Compound space version

# In[ ]:


u1.Interpolate(cf1)
u2.Interpolate(cf2)
uvec.Interpolate(cf)
ncf = NewtonCF(eq, u, maxiter=1)
uvec.Interpolate(ncf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# ### Vector space version

# In[ ]:


u1.Interpolate(cf1)
u2.Interpolate(cf2)
uvec.Interpolate(cf)
ncf = NewtonCF(eq_vec, uvec, maxiter=1)
uvec.Interpolate(ncf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# ### Check solution

# In[ ]:


uvec2 = GridFunction(uvec.space)
uvec2.Interpolate(check)
assert np.allclose(uvec2.vec.FV().NumPy(), 0)


# In[ ]:





# # Compound case: one scalar, one symmetric matrix-valued and one vector-valued variable

# In[ ]:


fes_M = MatrixValued(fes_ir, dim=3, symmetric=True)

fes_compound = fes_ir * fes_M * (fes_ir**2)

u = GridFunction(fes_compound)
u1, u2, u3 = u.components
du1, du2, du3 = fes_compound.TrialFunction()

fes_vec = fes_ir**12
uvec = GridFunction(fes_vec)


# In[ ]:


cf1 = CoefficientFunction(1/2)
cf2 = CoefficientFunction(tuple(range(1, 10)), dims=(3, 3))
cf3 = CoefficientFunction((2/3, 3/4))
u1.Interpolate(cf1)
u2.Interpolate(1/2 * (cf2 + cf2.trans))
u3.Interpolate(cf3)


# In[ ]:


u.vec.FV().NumPy()


# ## Linear equation with nonsymmetric Jacobian

# In[ ]:


def sym(M):
    return 1/2 * (M + M.trans)


# In[ ]:


d = CoefficientFunction(1/2)
a = CoefficientFunction((2/3, 4/9))
M22 = CoefficientFunction((2,    1/2,
                           1/6,    4), dims=(2, 2))
M33 = CoefficientFunction((2,    1/2, 1/3, 
                           1/2,    4,   0,
                           1/3,    0,   2), dims=(3, 3))

# scalar eq component
def res_func_1(u1, u2, u3):
    return 4 * u1 + 2 * u2[0, 1] + 5 * u2[1,2] + (M22 * u3) * a 

# matrix eq component
def res_func_2(u1, u2, u3):
    return (M33 * u2 + M33 
        + CoefficientFunction((M33[0], 0, 0, 0, 0, 0, 0, 0, 0), dims=(3, 3))
        + CoefficientFunction((u1, u3[0], u3[1]/3, 
                               0, 0, 0, 
                               0, 0, 0), dims=(3, 3))
        )

# vector eq component
def res_func_3(u1, u2, u3):
    return a * u1 + M22 * u3 + 2 * a * u3[0] + M22 * a / 2 + M22 * CoefficientFunction((u2[0, 0], u2[2,2]))

    

def res_func(u1, u2, u3):
    return CoefficientFunction(
        (res_func_1(u1, u2, u3), 
         res_func_2(u1, u2, u3), 
         res_func_3(u1, u2, u3)))


eq = res_func(du1, du2, du3)

uv1 = CoefficientFunction(uvec[0])
uv2 = CoefficientFunction(tuple([uvec[i] for i in range(1, 10)]), dims=(3,3))
uv3 = CoefficientFunction(tuple([uvec[i] for i in range(10, 12)]))
check = res_func(uv1, uv2, uv3)

expected = np.array([-6.32030540e-01, -1.65811403e+00, -1.03161645e-01, -4.81744750e-04,
       -1.03161645e-01, -9.87104794e-01,  5.77134791e-03, -4.81744750e-04,
        5.77134791e-03, -9.99919709e-01,  9.39655500e-01,  6.55157652e-01])


# In[ ]:


u1.Interpolate(cf1)
u2.Interpolate(1/2 * (cf2 + cf2.trans))
u3.Interpolate(cf3)
ncf = NewtonCF(eq, u.components, maxiter=1)
uvec.Interpolate(ncf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# In[ ]:


u1.Interpolate(cf1)
u2.Interpolate(1/2 * (cf2 + cf2.trans))
u3.Interpolate(cf3)
ncf = NewtonCF(eq, u, maxiter=1)
uvec.Interpolate(ncf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# In[ ]:


u1.Interpolate(cf1)
u2.Interpolate(1/2 * (cf2 + cf2.trans))
u3.Interpolate(cf3)
ncf = NewtonCF(eq, CoefficientFunction((u1, u2, u3)), maxiter=1)
uvec.Interpolate(ncf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# In[ ]:


u1.Interpolate(cf1)
u2.Interpolate(1/2 * (cf2 + cf2.trans))
u3.Interpolate(cf3)
ncf = NewtonCF(eq, [u1, u2, u3], maxiter=1)
uvec.Interpolate(ncf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# In[ ]:


u1.Interpolate(cf1)
u2.Interpolate(1/2 * (cf2 + cf2.trans))
u3.Interpolate(cf3)
ncf = NewtonCF(eq, [u1, u2, u3], maxiter=1)
uvec.Interpolate(ncf)
uvec.vec.FV().NumPy()
# This interpolation is currently not supported
#u.Interpolate(ncf)
#u.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# ### Check solution
# 
# Not that the error for the symmetric-matrix-valued component `u2` can have a non-vanishing skew-symmetric component.

# In[ ]:


uvec_check = GridFunction(uvec.space)
uvec_check.Interpolate(check)
nvec = uvec_check.vec.FV().NumPy()
allres = nvec.reshape((12, int(nvec.size/12))).T


# In[ ]:


# Only the symmetric part of the error for u2 is relevant and indeed vanishes
_res = allres[0][1:10].reshape((3,3))
assert np.allclose(1/2 * (_res + _res.T), 0)


# ### Interpolate component-wise to GF via an intermediate GF (`uvec`)

# In[ ]:


u1.Interpolate(uvec[0]), u2.Interpolate(uvec[1:10]), u3.Interpolate(uvec[10:])


# In[ ]:


u2.vec.FV().NumPy(), uvec.vec.FV().NumPy()[1:10]


# ### Interpolation via explicit projection

# NOTE: SolveM and ApplyM are currently not implemented for IRSpace but could be done easily.
# Also, currently Mass() does neighter respect quadrature weights nor VS embeddings -- just seems to be identity (see below)!

# In[ ]:


v1, v2, v3 = fes_compound.TestFunction()


# In[ ]:


MM = fes_compound.Mass()


# In[ ]:


M = BilinearForm(fes_compound, symmetric=True)
M += InnerProduct(CoefficientFunction((du1, du2, du3)), CoefficientFunction((v1, v2, v3))) * irs_dx
M.Assemble()


# In[ ]:


# CacheCF crashes because ProxyUserData->caches is not set and not properly initialized. 
# The problem is the non-default initialization of size!
# Fixed locally by default-inialization of FlatArray members!
L = LinearForm(fes_compound)
ff = CoefficientFunction((1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
L += InnerProduct(CacheCF(ncf), CoefficientFunction((v1, v2, v3))) * irs_dx
L.Assemble()


# In[ ]:


L.vec.FV().NumPy()


# In[ ]:


u.vec.FV().NumPy()[:] = 0
u.vec.data = M.mat.Inverse() * L.vec
u.vec.FV().NumPy()


# In[ ]:





# In[ ]:


# fes.Mass() --> wrong!
u.vec.FV().NumPy()[:] = 0
u.vec.data = MM.Inverse() * L.vec
u.vec.FV().NumPy()


# ### Interpolation with [dual basis](https://ngsolve.org/docu/nightly/i-tutorials/unit-2.10-dualbasis/dualbasis.html)
# 
# Does not work ATM.

# In[ ]:


for v in fes_compound.TestFunction():
    print(type(v))
for v in fes_compound.TestFunction():
    print(v.Operator("dual"))
Du1, Du2, Du3 = [v.Operator("dual") for v in fes_compound.TestFunction()]


# In[ ]:


fes_compound.components[2].TestFunction().Operator("dual")


# ## Linear equation with symmetric Jacobian

# In[ ]:


def sym(M):
    return 1/2 * (M + M.trans)


# In[ ]:


d = CoefficientFunction(3)
a = CoefficientFunction((2/3, 4/9))
M22 = CoefficientFunction((1/8,    1/5,
                           1/5,    3), dims=(2, 2))
M33 = CoefficientFunction((2,    1/2, 1/3, 
                           1/2,    4,   1/2,
                           1/3,    1/2,   2), dims=(3, 3))

def pot_func(u1, u2, u3):
    return d * u1**2 + InnerProduct(u2, M33 * u2) + InnerProduct(u3, M22 * u3)         + 2 * u1 + 4 * InnerProduct(M33, u2) + InnerProduct(u3, a)

def res_func(u1, u2, u3):
    return CoefficientFunction(
        (2 * d *u1 + 2, 2 * M33 * u2 + 4 * M33, 2 * M22 * u3 + a)
    )

pot = pot_func(du1, du2, du3)
res = res_func(du1, du2, du3)

uv1 = CoefficientFunction(uvec[0])
uv2 = CoefficientFunction(tuple([uvec[i] for i in range(1, 10)]), dims=(3,3))
uv3 = CoefficientFunction(tuple([uvec[i] for i in range(10, 12)]))
check_pot = pot_func(uv1, uv2, uv3)
check_res = res_func(uv1, uv2, uv3)

expected = np.array([-3.33333333e-01, -2.00000000e+00,  8.88178420e-16,  0.00000000e+00,
        8.88178420e-16, -2.00000000e+00,  8.88178420e-16,  0.00000000e+00,
        8.88178420e-16, -2.00000000e+00, -2.85240464e+00,  1.16086235e-01])


# ### Newton version

# In[ ]:


u1.Interpolate(cf1)
u2.Interpolate(1/2 * (cf2 + cf2.trans))
u3.Interpolate(cf3)
print(u.vec.FV().NumPy())
ncf = NewtonCF(res, u.components, maxiter=1)
uvec.Interpolate(ncf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# In[ ]:


uvec_check = GridFunction(uvec.space)
uvec_check.Interpolate(check_res)
nvec = uvec_check.vec.FV().NumPy()
allres = nvec.reshape((12, int(nvec.size/12))).T


# In[ ]:


assert np.allclose(allres[0][0], 0)
assert np.allclose(allres[0][10:], 0)
# Only the symmetric part of the error for u2 is relevant and indeed vanishes
_res = allres[0][1:10].reshape((3,3))
assert np.allclose(1/2 * (_res + _res.T), 0)


# ### Minimization version

# In[ ]:


u1.Interpolate(cf1)
u2.Interpolate(1/2 * (cf2 + cf2.trans))
u3.Interpolate(cf3)
print(u.vec.FV().NumPy())
mcf = MinimizationCF(pot, u.components, maxiter=1)
uvec.Interpolate(mcf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# In[ ]:





# ## Nonlinear equation with symmetric Jacobian

# In[ ]:


def sym(M):
    return 1/2 * (M + M.trans)


# In[ ]:


d = CoefficientFunction(3)
a = CoefficientFunction((2/3, 4/9))
M22 = CoefficientFunction((1/8,    1/5,
                           1/5,    3), dims=(2, 2))
M33 = CoefficientFunction((2,    1/2, 1/3, 
                           1/2,    4,   1/2,
                           1/3,    1/2,   2), dims=(3, 3))

wc1 = CoefficientFunction(1)
wc2 = CoefficientFunction(1)
wd = CoefficientFunction(50)

def pot_func(u1, u2, u3):
    return d * (u1**4 + u1**2) + wc1 * u1 * Det(u2) + wd * (sqrt(Det(u2)) - 1)**2 + 10 * InnerProduct(u2, u2)         + wc2 * InnerProduct(u3, u3) * Trace(u2) + InnerProduct(u3, M22 * u3)         + 2 * u1 + 4 * InnerProduct(M33, u2) + InnerProduct(u3, a)

def res_func(u1, u2, u3):
    return CoefficientFunction((
        4 * d *u1**3 + 2 * d * u1 + wc1 * Det(u2) + 2, 
        (wc1 * u1 * Det(u2) + 2 * wd * (sqrt(Det(u2)) - 1) * (1/2 * sqrt(Det(u2)))) * Inv(u2) 
             + 20 * u2 + wc2 * InnerProduct(u3, u3) * Id(3) + 4 * M33, 
        2 * wc2 * Trace(u2) * u3 + 2 * M22 * u3 + a
    ))

pot = pot_func(du1, du2, du3)
res = res_func(du1, du2, du3)

uv1 = CoefficientFunction(uvec[0])
uv2 = CoefficientFunction(tuple([uvec[i] for i in range(1, 10)]), dims=(3,3))
uv3 = CoefficientFunction(tuple([uvec[i] for i in range(10, 12)]))
check_pot = pot_func(uv1, uv2, uv3)
check_res = res_func(uv1, uv2, uv3)

expected = np.array([-0.30564104,  0.60862982, -0.03132745, -0.02373416, -0.03132745,
        0.48047082, -0.03132745, -0.02373416, -0.03132745,  0.60862982,
       -0.17851929, -0.03970393])


# ### Minimization version

# In[ ]:


u1.Interpolate(CoefficientFunction(0))
u2.Interpolate(Id(3))
u3.Interpolate(CoefficientFunction((0, 0)))
print(u.vec.FV().NumPy())
mcf = MinimizationCF(pot, u.components)
uvec.Interpolate(mcf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# In[ ]:


uvec_check = GridFunction(uvec.space)
uvec_check.Interpolate(check_res)
nvec = uvec_check.vec.FV().NumPy()
allres = nvec.reshape((12, int(nvec.size/12))).T


# In[ ]:


assert np.allclose(allres[0][0], 0)
assert np.allclose(allres[0][10:], 0)
# Only the symmetric part of the error for u2 is relevant and indeed vanishes
_res = allres[0][1:10].reshape((3,3))
assert np.allclose(1/2 * (_res + _res.T), 0)


# ### Newton version

# In[ ]:


u1.Interpolate(CoefficientFunction(0))
u2.Interpolate(Id(3))
u3.Interpolate(CoefficientFunction((0, 0)))
print(u1.vec.FV().NumPy())
print(u.vec.FV().NumPy())
ncf = NewtonCF(res, u.components)
uvec.Interpolate(ncf)
uvec.vec.FV().NumPy()


# In[ ]:


assert np.allclose(uvec.vec.FV().NumPy(), expected, atol=1e-8, rtol=0)


# In[ ]:


uvec_check = GridFunction(uvec.space)
uvec_check.Interpolate(check_res)
nvec = uvec_check.vec.FV().NumPy()
allres = nvec.reshape((12, int(nvec.size/12))).T


# In[ ]:


assert np.allclose(allres[0][0], 0)
assert np.allclose(allres[0][10:], 0)
# Only the symmetric part of the error for u2 is relevant and indeed vanishes
_res = allres[0][1:10].reshape((3,3))
assert np.allclose(1/2 * (_res + _res.T), 0)


# In[ ]:




