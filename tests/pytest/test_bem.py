
import pytest

from netgen.occ import *
from ngsolve import *
from ngsolve.bem import *

def test_transform():
    box = Box((-10,-10,-10), (10,10,10))
    mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=5))

    kappa = 20
    order = 80
    mp = SingularMultiPoleCF(order, kappa, (0,0,0), rad=1)
    mp.AddCharge((0.3, -0.1,0.5), 1)

    dx,dy,dz = 0.1, 0.2, 0.3
    mp2 = SingularMultiPoleCF(order, kappa, (dx,dy,dz), rad=2)
    mp.Transform(mp2)

    meshpnt = mesh(1,1,1)
    assert mp(meshpnt) == pytest.approx(mp2(meshpnt), rel=1e-10)

    regmp = RegularMultiPoleCF(order, kappa, (3,3,3), rad=2)
    mp2.Transform(regmp)
    meshpnt = mesh(2,2,2)
    assert mp(meshpnt) == pytest.approx(regmp(meshpnt), rel=1e-10)
    

    

    
def test_singularmltransform():
    box = Box((-10,-10,-10), (10,10,10))
    mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=5))

    kappa = 0.01
    order = 20
    mp = SingularMLMultiPoleCF((0,0,0), r=1, kappa=kappa)
    num = 200
    for i in range(num):
        z = i/num
        mp.mlmp.AddCharge((0.1, 0, z), 1/num)

    val1 = mp(mesh(1,1,4))
    mp.mlmp.Calc()
    val2 = mp(mesh(1,1,4))

    assert val1 == pytest.approx(val2)

