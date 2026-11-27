"""
ngsolve.convert: the old release evaluates a pickled GridFunction at points,
the current one projects the values onto its space.

Without network the current interpreter plays the old release, which tests
everything but the change of basis. With NGSOLVE_CONVERT_OLD_PYTHON set to
the interpreter of an older release, the pickles are written by that one.
"""

import os, pickle, subprocess, sys

import pytest
from ngsolve import *
from ngsolve.convert import ConvertGridFunction, PickleVersions
import ngsolve

oldpython = os.environ.get("NGSOLVE_CONVERT_OLD_PYTHON", sys.executable)

# written to a file and run by the old interpreter: old API only
writer = r'''
import pickle, sys
import numpy as np
from ngsolve import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve.meshes import MakeStructured3DMesh, MakeStructured2DMesh
case, out = sys.argv[1:3]
trigs = lambda: Mesh(unit_square.GenerateMesh(maxh=0.3))
tets = lambda: Mesh(unit_cube.GenerateMesh(maxh=0.5))
hexes = lambda: MakeStructured3DMesh(hexes=True, nx=2, ny=2, nz=2, mapping=lambda x, y, z: (x + 0.1*y*z, y + 0.1*x*z, z))
quads = lambda: MakeStructured2DMesh(quads=True, nx=3, ny=3, mapping=lambda x, y: (x + 0.1*x*y, y))
spaces = {
    "h1_trig": lambda: H1(trigs(), order=3, dirichlet="left|bottom"),
    "h1_complex": lambda: H1(trigs(), order=2, complex=True),
    "l2_quad": lambda: L2(quads(), order=2),
    "vectorh1_trig": lambda: VectorH1(trigs(), order=2),
    "hcurl_tet": lambda: HCurl(tets(), order=2),
    "hcurl_hex": lambda: HCurl(hexes(), order=2),
    "hdiv_trig": lambda: HDiv(trigs(), order=2),
    "hdiv_tet": lambda: HDiv(tets(), order=1),
    "hdivdiv_trig": lambda: HDivDiv(trigs(), order=1),
}
fes = spaces[case]()
gf = GridFunction(fes)
rng = np.random.default_rng(3)
vec = gf.vec.FV().NumPy()
vec[:] = rng.standard_normal(len(vec))
if fes.is_complex:
    vec[:] = vec + 1j * rng.standard_normal(len(vec))
with open(out, "wb") as f:
    pickle.dump(gf, f)
'''

# the order of the current space that contains the old one
cases = {"h1_trig": 0, "h1_complex": 0, "l2_quad": 0, "vectorh1_trig": 0, "hcurl_tet": 0,
         "hcurl_hex": 1, "hdiv_trig": 0, "hdiv_tet": 0, "hdivdiv_trig": 0}


@pytest.mark.parametrize("case", sorted(cases))
def test_convert(case, tmp_path):
    script = tmp_path / "writer.py"
    script.write_text(writer)
    old = str(tmp_path / "old.pkl")
    env = None if oldpython == sys.executable else {k: v for k, v in os.environ.items()
                                                   if k not in ("PYTHONPATH", "PYTHONHOME", "DYLD_LIBRARY_PATH")}
    subprocess.run([oldpython, str(script), case, old], check=True, env=env)
    assert "ngsolve" in PickleVersions(old)

    raise_order = cases[case] if oldpython != sys.executable else 0
    new = str(tmp_path / "new.pkl")
    gf = ConvertGridFunction(old, new, python=oldpython, raise_order=raise_order, verbose=False)
    with open(new, "rb") as f:
        gf2 = pickle.load(f)
    assert gf2.space.ndof == gf.space.ndof
    if oldpython == sys.executable:
        # same release, same basis: the coefficients themselves come back
        with open(old, "rb") as f:
            gfold = pickle.load(f)
        diff = gf.vec.CreateVector()
        diff.data = gf.vec - gfold.vec
        assert Norm(diff) < 1e-8 * Norm(gfold.vec)
