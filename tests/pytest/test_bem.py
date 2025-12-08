
import pytest

from netgen.occ import *
from ngsolve import *
from ngsolve.bem import *


@pytest.mark.parametrize(
    "order_source_1, order_source_2, order_target_1, order_target_2",
    [
        (80, 80, 80, 80),
        (80, 82, 79, 81),
    ],
)
def test_transform(order_source_1, order_source_2, order_target_1, order_target_2):
    box = Box((-10,-10,-10), (10,10,10))
    mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=5))

    kappa = 20
    S = SingularExpansionCF(order_source_1, kappa, (0,0,0), rad=1)
    S.AddCharge((0.3, -0.1,0.5), 1)

    dx,dy,dz = 0.1, 0.2, 0.3
    S2 = SingularExpansionCF(order_source_2, kappa, (dx,dy,dz), rad=2)
    S.Transform(S2)

    meshpnt = mesh(1,1,1)
    assert S(meshpnt) == pytest.approx(S2(meshpnt), rel=1e-10)

    R = RegularExpansionCF(order_target_1, kappa, (3,3,3), rad=2)
    S2.Transform(R)
    meshpnt = mesh(2,2,2)
    assert S(meshpnt) == pytest.approx(R(meshpnt), rel=1e-10)
    
    R2 = RegularExpansionCF(order_target_2, kappa, (3+dx,3+dy,3+dz), rad=1)
    R.Transform(R2)
    meshpnt = mesh(3,3,3)
    assert S(meshpnt) == pytest.approx(R2(meshpnt), rel=1e-10)

    

    
def test_singularmltransform():
    box = Box((-10,-10,-10), (10,10,10))
    mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=5))

    kappa = 0.01
    order = 20
    S = SingularMLExpansionCF((0,0,0), r=1, kappa=kappa)
    num = 200
    for i in range(num):
        z = i/num
        S.expansion.AddCharge((0.1, 0, z), 1/num)

    val1 = S(mesh(1,1,4))
    S.expansion.Calc()
    val2 = S(mesh(1,1,4))

    assert val1 == pytest.approx(val2)
