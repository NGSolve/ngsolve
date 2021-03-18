
from netgen.csg import *
from ngsolve import *
import numpy as np

def test_eval_curved_surface():
    geo = CSGeometry()
    cyl = Cylinder((0,0,-1), (0,0,2), 0.3)
    box = OrthoBrick((-1,-1,0), (1,1,1))
    geo.Add(cyl*box)
    mesh = Mesh(geo.GenerateMesh(maxh=0.1))
    mesh.Curve(5)

    fes = H1(mesh, order=3)
    u = GridFunction(fes)
    u.Set(sin(atan2(y,x)), dual=True)

    phivals = np.linspace(0, 2*np.pi, 1000)
    xvals = 0.3 * np.cos(phivals)
    yvals = 0.3 * np.sin(phivals)
    zvals = 0.5 * np.ones(5)
    pts = mesh(xvals, yvals, 0.5, VOL_or_BND=BND)
    uvals = u(pts).flatten()
    sinvals = np.sin(phivals)
    assert max(uvals-sinvals) < 1e-5

