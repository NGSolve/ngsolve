# Convergence study of the TDNNS note: writes conv.pdf and prints the table of rates.
#   python3 figures.py [cache.json]      (about 6 minutes without cache)
import sys, os, json, tempfile
from math import pi, log
from ngsolve import *
from ngsolve.meshes import MakeStructured2DMesh, MakeStructured3DMesh
from netgen.meshing import Mesh as NM, Element3D, Element2D, Element1D, MeshPoint, Pnt, FaceDescriptor

SetHeapSize(int(2e9))

def trigquad(n, alltrig=False):
    # unit square, every third cell split into two trigs
    m = NM(dim=2)
    m.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    p = {(i, j): m.Add(MeshPoint(Pnt(i/n, j/n, 0))) for i in range(n+1) for j in range(n+1)}
    for i in range(n):
        for j in range(n):
            a, b, c, d = p[i, j], p[i+1, j], p[i+1, j+1], p[i, j+1]
            if alltrig or (i+j) % 3 == 0:
                m.Add(Element2D(1, [a, b, c])); m.Add(Element2D(1, [a, c, d]))
            else:
                m.Add(Element2D(1, [a, b, c, d]))
    for i in range(n):
        m.Add(Element1D([p[i, 0], p[i+1, 0]], index=1)); m.Add(Element1D([p[n, i], p[n, i+1]], index=1))
        m.Add(Element1D([p[i+1, n], p[i, n]], index=1)); m.Add(Element1D([p[0, i+1], p[0, i]], index=1))
    return Mesh(m)

def prismhex(n):
    # unit cube, columns of hexes alternating with columns of prism pairs
    m = NM(dim=3)
    m.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    P = {(i, j, k): m.Add(MeshPoint(Pnt(i/n, j/n, k/n))) for i in range(n+1) for j in range(n+1) for k in range(n+1)}
    bnd = lambda vs: m.Add(Element2D(1, vs))
    for i in range(n):
        for j in range(n):
            for k in range(n):
                a, b, c, d = P[i, j, k], P[i+1, j, k], P[i+1, j+1, k], P[i, j+1, k]
                e, f, g, h = P[i, j, k+1], P[i+1, j, k+1], P[i+1, j+1, k+1], P[i, j+1, k+1]
                hexcell = (i+j) % 2 == 0
                if hexcell:
                    m.Add(Element3D(1, [a, b, c, d, e, f, g, h]))
                else:
                    m.Add(Element3D(1, [a, b, c, e, f, g])); m.Add(Element3D(1, [a, c, d, e, g, h]))
                if k == 0:
                    if hexcell: bnd([a, d, c, b])
                    else: bnd([a, c, b]); bnd([a, d, c])
                if k == n-1:
                    if hexcell: bnd([e, f, g, h])
                    else: bnd([e, f, g]); bnd([e, g, h])
                if j == 0: bnd([a, b, f, e])
                if j == n-1: bnd([c, d, h, g])
                if i == 0: bnd([d, a, e, h])
                if i == n-1: bnd([b, c, g, f])
    return Mesh(m)

def tetprism(n):
    # unstructured tets with one prism layer at the boundary
    from netgen.occ import Box, OCCGeometry
    from netgen.meshing import BoundaryLayerParameters
    bl = BoundaryLayerParameters(boundary=".*", thickness=[1/(3*n)])
    return Mesh(OCCGeometry(Box((0, 0, 0), (1, 1, 1))).GenerateMesh(maxh=1/n, boundary_layers=[bl]))

MESHES = {
    "trig":       lambda n: trigquad(n, alltrig=True),
    "quad":       lambda n: MakeStructured2DMesh(quads=True, nx=n, ny=n),
    "trig+quad":  trigquad,
    "tet":        lambda n: MakeStructured3DMesh(hexes=False, nx=n, ny=n, nz=n),
    "prism":      lambda n: MakeStructured3DMesh(hexes=False, prism=True, nx=n, ny=n, nz=n),
    "hex":        lambda n: MakeStructured3DMesh(hexes=True, nx=n, ny=n, nz=n),
    "prism+hex":  prismhex,
    "tet+prism":  tetprism,
}
LEVELS2D = {1: [4, 8, 16, 32, 64], 2: [4, 8, 16, 32, 64], 3: [4, 8, 16, 32, 64]}
LEVELS3D = {1: [2, 4, 8, 16], 2: [2, 4, 8, 12], 3: [2, 4, 8]}

def solve(mesh, m, type1, nu=0.3):
    # hybrid TDNNS, Dirichlet data of a smooth solution, relative L2 errors
    dim = mesh.dim
    mu, lam = 1/2/(1+nu), nu/((1+nu)*(1-2*nu))
    A = lambda s: 1/(2*mu)*(s - lam/(2*mu+dim*lam)*Trace(s)*Id(dim))
    cs = (x, y, z)[:dim]
    b = sin(pi*x)*sin(pi*y)*(sin(pi*z) if dim == 3 else 1)
    uex = CF((b, 2*x*b) if dim == 2 else (b, 2*x*b, (1+y*y)*b))
    gu = CF(tuple(uex[i].Diff(c) for i in range(dim) for c in cs), dims=(dim, dim))
    sex = mu*(gu+gu.trans) + lam*Trace(gu)*Id(dim)
    f = -CF(tuple(sum(sex[i, j].Diff(c) for j, c in enumerate(cs)) for i in range(dim)))

    Q = HCurl(mesh, order=m, type1=type1)
    S = Discontinuous(Q.CreateStressSpace())
    M = S.CreateMultiplierSpace()
    X = S*Q*M
    (sigma, u, un), (tau, v, vn) = X.TnT()
    n = specialcf.normal(dim)
    tang = lambda w: w - (w*n)*n
    nn = lambda s: (s*n)*n
    a = -(InnerProduct(A(sigma), tau) + div(sigma)*v + div(tau)*u)*dx \
        + ((sigma*n)*tang(v) + (tau*n)*tang(u) + nn(sigma)*(vn*n) + nn(tau)*(un*n))*dx(element_boundary=True)
    with TaskManager():
        gf = Solve(a == f*v*dx, u[BND(".*")] == uex, un[BND(".*")] == uex)
    sh, uh, _ = gf.components
    nrm = lambda c: sqrt(Integrate(InnerProduct(c, c), mesh, order=2*m+6))
    return dict(ne=mesh.ne, ncoupling=X.FreeDofs(True).NumSet(), err_s=nrm(sh-sex)/nrm(sex), err_u=nrm(uh-uex)/nrm(uex))

def compute():
    res = []
    for kind, mk in MESHES.items():
        levels = LEVELS2D if kind in ("trig", "quad", "trig+quad") else LEVELS3D
        for type1 in (False, True):
            for m in (1, 2, 3):
                for n in levels[m]:
                    if kind in ("prism+hex", "tet+prism") and n > 8: continue
                    r = solve(mk(n), m, type1)
                    r.update(kind=kind, type1=type1, m=m, n=n)
                    res.append(r)
                    print(r, flush=True)
    return res

def rates(res, kind, type1, m):
    rs = [r for r in res if r["kind"] == kind and r["type1"] == type1 and r["m"] == m]
    d = 2 if kind in ("trig", "quad", "trig+quad") else 3
    a, b = rs[-2], rs[-1]
    hr = (b["ne"]/a["ne"])**(1/d)
    return log(a["err_s"]/b["err_s"])/log(hr), log(a["err_u"]/b["err_u"])/log(hr)

def table(res):
    print("rates sigma/u between the two finest meshes")
    print(f"{'mesh':10} " + " ".join(f"{'default m='+str(m):>14}" for m in (1, 2, 3)) + " " + " ".join(f"{'type1 m='+str(m):>14}" for m in (1, 2, 3)))
    for kind in MESHES:
        print(f"{kind:10} " + " ".join("{:6.2f}/{:5.2f}  ".format(*rates(res, kind, t1, m)) for t1 in (False, True) for m in (1, 2, 3)))

def plot(res, fname):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 7, "axes.linewidth": 0.5, "lines.linewidth": 0.9, "lines.markersize": 3.2,
                         "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.minor.width": 0.3, "ytick.minor.width": 0.3})
    col = {1: "#2a78d6", 2: "#eb6834", 3: "#1baf7a"}
    mark = {"trig": "^", "quad": "s", "trig+quad": "o", "tet": "^", "prism": "D", "hex": "s", "prism+hex": "o", "tet+prism": "v"}
    fig, axs = plt.subplots(1, 4, figsize=(7.0, 1.8), sharey=False)
    panels = [(("trig", "quad", "trig+quad"), "err_s", r"2D: $\|\sigma-\sigma_h\|/\|\sigma\|$"),
              (("trig", "quad", "trig+quad"), "err_u", r"2D: $\|u-u_h\|/\|u\|$"),
              (("tet", "prism", "hex", "prism+hex", "tet+prism"), "err_s", r"3D: $\|\sigma-\sigma_h\|/\|\sigma\|$"),
              (("tet", "prism", "hex", "prism+hex", "tet+prism"), "err_u", r"3D: $\|u-u_h\|/\|u\|$")]
    for ax, (kinds, key, title) in zip(axs, panels):
        d = 2 if "trig" in kinds else 3
        for m in (1, 2, 3):
            for kind in kinds:
                rs = [r for r in res if r["kind"] == kind and not r["type1"] and r["m"] == m]
                ax.loglog([r["ncoupling"] for r in rs], [r[key] for r in rs], marker=mark[kind], color=col[m],
                          markerfacecolor="white", markeredgewidth=0.7)
            # reference slope h^p, p = m (sigma), m+1 (u)
            p = m if key == "err_s" else m+1
            rs = [r for r in res if r["kind"] == kinds[1] and not r["type1"] and r["m"] == m]
            n0, n1, e1 = rs[1]["ncoupling"], rs[-1]["ncoupling"], rs[-1][key]
            ax.loglog([n0, n1], [0.4*e1*(n1/n0)**(p/d), 0.4*e1], color="0.45", lw=0.5, ls="--")
            ax.annotate(f"$h^{p}$", (n1, 0.4*e1), textcoords="offset points", xytext=(4, -6), color="0.3", fontsize=6)
            first = min((r for k in kinds for r in res if r["kind"] == k and not r["type1"] and r["m"] == m),
                        key=lambda r: r["ncoupling"])
            ax.annotate(f"$m={m}$", (first["ncoupling"], first[key]), textcoords="offset points", xytext=(-5, -2),
                        ha="right", color="0.15", fontsize=6)
        ax.set_title(title, fontsize=7, pad=3)
        ax.set_xlabel("coupling dofs", labelpad=1)
        ax.grid(True, which="major", lw=0.3, color="0.88")
        ax.tick_params(length=2, pad=1.5)
        ax.margins(x=0.24, y=0.07)
        hs = [plt.Line2D([], [], marker=mark[k], color="0.25", lw=0, markerfacecolor="white", markeredgewidth=0.7, label=k) for k in kinds]
        if key == "err_s":
            ax.legend(handles=hs, loc="lower left", frameon=False, fontsize=5.6, handletextpad=0.1, borderpad=0.1, labelspacing=0.2)
    fig.subplots_adjust(left=0.045, right=0.985, bottom=0.185, top=0.895, wspace=0.27)
    fig.savefig(fname)
    fig.savefig(os.path.splitext(fname)[0] + ".png", dpi=250)  # for the html note

if __name__ == "__main__":
    cache = sys.argv[1] if len(sys.argv) > 1 else os.path.join(tempfile.gettempdir(), "tdnns_conv.json")
    if os.path.exists(cache):
        res = json.load(open(cache))
    else:
        res = compute()
        json.dump(res, open(cache, "w"))
    table(res)
    plot(res, os.path.join(os.path.dirname(os.path.abspath(__file__)), "conv.pdf"))
