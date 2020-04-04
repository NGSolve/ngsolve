
import pytest
from elements import *
from ngsolve import *
import scipy.sparse as sps
import numpy as np

@pytest.mark.parametrize("element", [Trig, Quad, Tet, Prism])
def test_dual_hcurl(element):
    mesh = Mesh(element())
    fes = HCurl(mesh, order=0)
    u,v = fes.TnT()
    vdual = v.Operator("dual")
    a = BilinearForm(fes)
    print("mesh.dim = ", mesh.dim)
    a += u * vdual * dx(element_vb=BBND if mesh.dim == 3 else BND)
    a.Assemble()
    dense = sps.csr_matrix(a.mat.CSR()).todense()
    assert dense == pytest.approx(np.identity(a.mat.height)), "mat = " + str(dense)
