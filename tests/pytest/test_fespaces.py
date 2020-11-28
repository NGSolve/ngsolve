import pytest
from ngsolve import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube

# define spaces and options for testing
spaces2d={"h1ho" :        { "vorb" : [VOL,BND,BBND], "order" : [3], "quad" : True }, \
          "facet" :       { "vorb" : [VOL,BND]     , "order" : [3], "quad" : True }, \
          "HDG" :         { "vorb" : [VOL,BND]     , "order" : [3], "quad" : True }, \
          "nodal" :       { "vorb" : [VOL]         , "order" : [1], "quad" : True }, \
          "nonconforming":{ "vorb" : [VOL]         , "order" : [3], "quad" : True }, \
          "VectorH1" :    { "vorb" : [VOL,BND,BBND], "order" : [3], "quad" : True }, \
          "hcurlho" :     { "vorb" : [VOL,BND,BBND], "order" : [3], "quad" : True }, \
          "hdiv" :        { "vorb" : [VOL,BND]     , "order" : [1], "quad" : True }, \
          "hdivho" :      { "vorb" : [VOL,BND]     , "order" : [3], "quad" : True }, \
          "l2" :          { "vorb" : [VOL,BND]     , "order" : [1], "quad" : True }, \
          "l2ho" :        { "vorb" : [VOL,BND]     , "order" : [3], "quad" : True }, \
          "l2surf" :      { "vorb" : [BND]         , "order" : [3], "quad" : True }, \
          "vectorfacet" : { "vorb" : [VOL,BND]     , "order" : [3], "quad" : True }, \
          "number" :      { "vorb" : [VOL]         , "order" : [3], "quad" : True }, \
          "hdivdiv" :     { "vorb" : [VOL]         , "order" : [3], "quad" : True }, \
          # "hcurl" :       { "vorb" : [VOL,BND]     , "order" : [1], "quad" : True }, \
}

spaces3d={"h1ho" :        { "vorb" : [VOL,BND,BBND,BBBND], "order" : [4] }, \
          "facet" :       { "vorb" : [VOL,BND]           , "order" : [3] }, \
          "HDG" :         { "vorb" : [VOL,BND]           , "order" : [3] }, \
          "nodal" :       { "vorb" : [VOL]               , "order" : [1] }, \
          "nonconforming":{ "vorb" : [VOL]               , "order" : [3] }, \
          "VectorH1" :    { "vorb" : [VOL,BND,BBND,BBBND], "order" : [4] }, \
          "hcurlho" :     { "vorb" : [VOL,BND,BBND]      , "order" : [3] }, \
          "hdiv" :        { "vorb" : [VOL,BND]           , "order" : [1] }, \
          "hdivho" :      { "vorb" : [VOL,BND]           , "order" : [3] }, \
          "l2" :          { "vorb" : [VOL,BND]           , "order" : [1] }, \
          "l2ho" :        { "vorb" : [VOL,BND]           , "order" : [3] }, \
          "l2surf" :      { "vorb" : [BND]               , "order" : [3] }, \
          "vectorfacet" : { "vorb" : [VOL,BND]           , "order" : [3] }, \
          "number" :      { "vorb" : [VOL]               , "order" : [3] }, \
          "hdivdiv" :     { "vorb" : [VOL,BND]           , "order" : [3] }, \
          # "hcurl" :       { "vorb" : [VOL,BND]           , "order" : [1] }, \
}

surfacespaces = {"h1ho" :          { "vorb" : [BND,BBND,BBBND], "order" : [3], "quad" : True }, \
                 "facetsurface" :  { "vorb" : [BND,BBND]      , "order" : [3], "quad" : True }, \
                 "VectorH1" :      { "vorb" : [BND,BBND,BBBND], "order" : [3], "quad" : True }, \
                 "hcurlho" :       { "vorb" : [BND,BBND]      , "order" : [3], "quad" : True }, \
                 "hdivhosurface" : { "vorb" : [BND,BBND]      , "order" : [3], "quad" : True }, \
                 "l2surf" :        { "vorb" : [BND,BBND]      , "order" : [3], "quad" : True }, \
                 "number" :        { "vorb" : [BND,BBND]      , "order" : [3], "quad" : True }, \
                 "hdivdivsurf" :   { "vorb" : [BND,BBND]      , "order" : [3], "quad" : True }, \
}


# 2d tests
def test_2DGetFE(quads=False):
    mesh = Mesh(unit_square.GenerateMesh(quad=quads))
    for spacename in spaces2d.keys():
        if not quads or (quads and spaces2d[spacename]["quad"]):
            for order in spaces2d[spacename]["order"]:
                space = FESpace(type=spacename,mesh=mesh,order=order)
                for vb in spaces2d[spacename]["vorb"]:
                    for el in space.Elements(vb):
                        assert space.GetFE(el).ndof == len(space.GetDofNrs(el)), [spacename,vb,order]

    return

# 3d tests
def test_3DGetFE():
    mesh = Mesh(unit_cube.GenerateMesh())
    for spacename in spaces3d.keys():
        for order in spaces3d[spacename]["order"]:
            space = FESpace(type=spacename,mesh=mesh,order=order)
            for vb in spaces3d[spacename]["vorb"]:
                for el in space.Elements(vb):
                    assert space.GetFE(el).ndof == len(space.GetDofNrs(el)), [spacename,vb,order]
    return

def test_SurfaceGetFE(quads=False):
    # 2d surface in 3d tests
    import netgen.meshing as meshing
    import netgen.csg as csg 
    geo = csg.CSGeometry()
    bottom   = csg.Plane (csg.Pnt(0,0,0), csg.Vec(0,0, 1) )
    surface = csg.SplineSurface(bottom)
    pts = [(0,0,0),(0,1,0),(1,1,0),(1,0,0)]
    geopts = [surface.AddPoint(*p) for p in pts]
    for p1,p2,bc in [(0,1,"left"), (1, 2,"top"),(2,3,"right"),(3,0,"bottom")]:
        surface.AddSegment(geopts[p1],geopts[p2],bc)
    geo.AddSplineSurface(surface)
    
    mesh = Mesh(geo.GenerateMesh(perfstepsend=meshing.MeshingStep.MESHSURFACE,quad=quads))
    
    for spacename in surfacespaces.keys():
        if not quads or (quads and spaces2d[spacename]["quad"]):
            for order in surfacespaces[spacename]["order"]:
                space = FESpace(type=spacename,mesh=mesh,order=order)
                for vb in surfacespaces[spacename]["vorb"]:
                    for el in space.Elements(vb):
                        assert space.GetFE(el).ndof == len(space.GetDofNrs(el)), [spacename,vb,order]
    return

if __name__ == "__main__":
    test_2DGetFE(quads=False)
    test_2DGetFE(quads=True)
    test_3DGetFE()
    test_SurfaceGetFE(quads=False)
    test_SurfaceGetFE(quads=True)
