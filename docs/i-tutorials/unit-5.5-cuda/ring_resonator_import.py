import math
from ngsolve import *
import netgen.geom2d as gm

SetHeapSize(50*1000*1000)
# geometrical and material properties
geo = gm.SplineGeometry()

xneg  =-0.43
xpos  = 0.43
yneg  =-0.48
ypos  = 0.48
wslab = 0.04
cringx = 0.0
cringy = 0.0
rring = 0.4
gap   = 0.005

pntx = [xneg,xpos]
pnty = [yneg,-rring-gap-wslab,-rring-gap,rring+gap,rring+gap+wslab,ypos]


pts = []
for yi in pnty:
    for xi in pntx:
        pts.append (geo.AddPoint(xi,yi))

#### parameters for source position


### geometry #######################################################
#inner rects
geo.Append (["line", pts[0], pts[1]], leftdomain=1, rightdomain=0)
geo.Append (["line", pts[1], pts[3]], leftdomain=1, rightdomain=0)
geo.Append (["line", pts[3], pts[2]], leftdomain=1, rightdomain=2)
geo.Append (["line", pts[2], pts[0]], leftdomain=1, rightdomain=0)

geo.Append (["line", pts[3], pts[5]], leftdomain=2, rightdomain=0,bc="normal_wg_rightbottom")
geo.Append (["line", pts[5], pts[4]], leftdomain=2, rightdomain=3)
geo.Append (["line", pts[4], pts[2]], leftdomain=2, rightdomain=0,bc="normal_wg_leftbottom")

geo.Append (["line", pts[5], pts[7]], leftdomain=3, rightdomain=0)
geo.Append (["line", pts[7], pts[6]], leftdomain=3, rightdomain=4)
geo.Append (["line", pts[6], pts[4]], leftdomain=3, rightdomain=0)

geo.Append (["line", pts[7], pts[9]], leftdomain=4, rightdomain=0,bc="normal_wg_righttop")
geo.Append (["line", pts[9], pts[8]], leftdomain=4, rightdomain=5)
geo.Append (["line", pts[8], pts[6]], leftdomain=4, rightdomain=0,bc="normal_wg_lefttop")

geo.Append (["line", pts[9], pts[11]], leftdomain=5, rightdomain=0)
geo.Append (["line", pts[11], pts[10]], leftdomain=5, rightdomain=0)
geo.Append (["line", pts[10], pts[8]], leftdomain=5, rightdomain=0)

geo.AddCircle(c=(cringx,cringy), r=rring, leftdomain=6, rightdomain=3)
geo.AddCircle(c=(cringx,cringy), r=rring-wslab, leftdomain=7, rightdomain=6)

for i in (1,3,5,7):
    geo.SetMaterial(i, "air")
for i in (2,4,6):
    geo.SetMaterial(i, "eps_nine")
data = geo.CreatePML(0.05)
normals = data["normals"]


mesh = Mesh(geo.GenerateMesh(maxh=0.05)).Curve(3)

eps_r = {"air" : 1, "eps_nine" : 3**3}

order = 3

for mat in mesh.GetMaterials():
    if mat.startswith("pml_normal_wg"):
        eps_r[mat] = eps_r["eps_nine"]

for mat in mesh.GetMaterials():
    if mat not in eps_r:
        eps_r[mat] = eps_r["air"]



### Parameters for Source field ##########################################################################
wavelength = 1.542
fcen       = 5/wavelength
df         = 0.1
tpeak      = 1

fes_facet = FacetFESpace(mesh, order=order+1)
gfsource = GridFunction(fes_facet)

source_cf =  sin( (pi/wslab)*(y-pnty[3]) ) 
gfsource.Set(source_cf,definedon=mesh.Boundaries("normal_wg_lefttop"))


##########################################################################################################


fes_u = VectorL2(mesh, order=order, piola=True, order_policy=ORDER_POLICY.CONSTANT)
fes_p = L2(mesh, order=order+1, all_dofs_together=True, order_policy=ORDER_POLICY.CONSTANT)
fes_tr = FacetFESpace(mesh, order=order+1)
fes_hdiv = HDiv(mesh, order=order+1, orderinner=1)

n = specialcf.normal(2) 
# h = specialcf.mesh_size

p,q = fes_p.TnT()
u,v = fes_u.TnT()
ptr,qtr = fes_tr.TnT()
uhdiv,vhdiv = fes_hdiv.TnT()

Bel = BilinearForm(grad(p)*v*dx - p*(v*n)*dx(element_boundary=True), geom_free=True).Assemble()
Btr = BilinearForm(0.5*ptr*(n*v)*dx(element_boundary=True), geom_free=True).Assemble()
Bstab = BilinearForm(p*(vhdiv*n)*dx(element_boundary=True),  geom_free=True).Assemble()
Mstab = BilinearForm(uhdiv*vhdiv*dx, diagonal=True).Assemble()

Mstabinv = Mstab.mat.Inverse()


fes = fes_p*fes_p*fes_u*fes_u
q,qhat, v,vhat = fes.TestFunction()
Lsrc  = LinearForm(gfsource*q*dx(element_boundary=True)).Assemble()


nvec = { mat : ((normals[mat][0], normals[mat][1]) if mat in normals else (0,0)) for mat in mesh.GetMaterials() }

cfn = CF( [CF(nvec[mat]) for mat in mesh.GetMaterials()])
cft = CF( ( cfn[1], -cfn[0] ) )

pml1d = mesh.Materials("pml_default.*|pml_normal.*")

eps = CF([eps_r[mat] for mat in mesh.GetMaterials()])


# damping matrices
sigma = 10   # pml damping parameter
dampp1 = fes_p.Mass (eps, definedon=pml1d)
dampp2 = fes_p.Mass (eps, definedon=mesh.Materials("pml_corner"))
dampu1 = fes_u.Mass (OuterProduct(cfn,cfn), definedon=pml1d)
dampu2 = fes_u.Mass (OuterProduct(cft,cft), definedon=pml1d)

gfu = GridFunction(fes)
gfstab = GridFunction(fes_hdiv)
gfu.vec[:] = 0

# For the time stepping
tau = 2e-4
tend = 100
t = 0


w = gfu.vec.CreateVector()

emb_p, emb_phat, emb_u, emb_uhat = fes.embeddings

traceop = fes_p.TraceOperator(fes_tr, False)
fullB = emb_u @ (Bel.mat + Btr.mat @ traceop) @ emb_p.T

dampingu = emb_u @ dampu1 @ emb_u.T + (-emb_u + emb_uhat) @ dampu2 @ (emb_u.T + emb_uhat.T)
dampingp = emb_p @ dampp1 @ emb_p.T + emb_p @ dampp2 @ (2*emb_p.T-emb_phat.T) + emb_phat @ dampp2 @ emb_p.T

invmassp = fes_p.Mass(eps).Inverse()
invmassu = fes_u.Mass(Id(mesh.dim)).Inverse()
invp = emb_p @ invmassp @ emb_p.T + emb_phat @ invmassp @ emb_phat.T
invu = emb_u @ invmassu @ emb_u.T + emb_uhat @ invmassu @ emb_uhat.T

# time-envelope:
def Envelope(t):
    if abs((t-tpeak)/tpeak) < 1:
        return (2*exp(1)/sqrt(math.pi))*sin(2*math.pi*fcen*t)*exp (-1/(1-((t-tpeak)/tpeak)**2))
    else:
        return 0
