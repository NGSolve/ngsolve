
import pytest

from netgen.occ import *
from ngsolve import *
from ngsolve.bem import *
from ngsolve.bla import Vec3D

def test_transform():
    box = Box((-10,-10,-10), (10,10,10))
    mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=5))

    kappa = 0.01
    order = 20
    mp = SingularMultiPoleCF(order, kappa, Vec3D(0,0,0), rad=1)
    mp.AddCharge(Vec3D(0.3, -0.1,0.5), 1)

    dx,dy,dz = 0, 0.2, 0.3
    mp2 = SingularMultiPoleCF(order, kappa, Vec3D(dx,dy,dz), rad=2)
    mp.Transform(mp2)

    meshpnt = mesh(1,1,1)
    assert mp(meshpnt) == pytest.approx(mp2(meshpnt))

