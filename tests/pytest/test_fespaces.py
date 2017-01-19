import pytest
from ngsolve import *
from netgen.geom2d import unit_square

def test_GetFE():
    # 2d tests
    mesh = Mesh(unit_square.GenerateMesh())
    spaces = {"h1ho" : [VOL,BND,BBND],"facet" : [VOL,BND] ,"HDG" : [VOL,BND],"nodal" : [VOL],"nonconforming" : [VOL],"VectorH1" : [VOL,BND],"hcurlho" : [VOL,BND,BBND],"hdiv" : [VOL,BND],"hdivho" : [VOL,BND],"l2" : [VOL,BND],"l2ho" : [VOL,BND],"l2surf" : [BND],"vectorfacet" : [VOL,BND],"number" : [VOL],"hdivdiv" : [VOL,BND], "hcurl" : [VOL,BND,BBND]}
    
    for spacename in spaces.keys():
        space = FESpace(type=spacename,mesh=mesh,order=3)
        for vb in spaces[spacename]:
            for el in space.Elements(vb):
                assert(list(el.dofs) == list(space.GetDofNrs(el)))
