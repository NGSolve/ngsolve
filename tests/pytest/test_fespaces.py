import pytest
from ngsolve import *
from netgen.geom2d import unit_square

# define spaces and options for testing
spaces = {"h1ho" :        { "vorb" : [VOL,BND,BBND], "order" : [3] }, \
          "facet" :       { "vorb" : [VOL,BND]     , "order" : [3] }, \
          "HDG" :         { "vorb" : [VOL,BND]     , "order" : [3] }, \
          "nodal" :       { "vorb" : [VOL]         , "order" : [1] }, \
          "nonconforming":{ "vorb" : [VOL]         , "order" : [3] }, \
          "VectorH1" :    { "vorb" : [VOL,BND]     , "order" : [3] }, \
          "hcurlho" :     { "vorb" : [VOL,BND,BBND], "order" : [3] }, \
          "hdiv" :        { "vorb" : [VOL,BND]     , "order" : [1] }, \
          "hdivho" :      { "vorb" : [VOL,BND]     , "order" : [3] }, \
          "l2" :          { "vorb" : [VOL,BND]     , "order" : [1] }, \
          "l2ho" :        { "vorb" : [VOL,BND]     , "order" : [3] }, \
          "l2surf" :      { "vorb" : [BND]         , "order" : [3] }, \
          "vectorfacet" : { "vorb" : [VOL,BND]     , "order" : [3] }, \
          "number" :      { "vorb" : [VOL]         , "order" : [3] }, \
          "hdivdiv" :     { "vorb" : [VOL]         , "order" : [3] }, \
          # "hcurl" :       { "vorb" : [VOL,BND]     , "order" : [1] }, \
}


def test_GetFE():
    # 2d tests
    mesh = Mesh(unit_square.GenerateMesh())
    for spacename in spaces.keys():
        for order in spaces[spacename]["order"]:
            space = FESpace(type=spacename,mesh=mesh,order=order)
            for vb in spaces[spacename]["vorb"]:
                for el in space.Elements(vb):
                    assert space.GetFE(el).ndof == len(space.GetDofNrs(el)), [spacename,vb,order]
                    
