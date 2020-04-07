
import pytest
from netgen.meshing import *

def Trig():
    mesh = Mesh(2)
    mesh.AddRegion("",2)
    pnts = [mesh.Add(MeshPoint(Pnt(0,0,0))), mesh.Add(MeshPoint(Pnt(1,0,0))), mesh.Add(MeshPoint(Pnt(0,1,0)))]
    mesh.Add(Element2D(1, pnts))
    return mesh

def Quad():
    mesh = Mesh(2)
    mesh.AddRegion("",2)
    pnts = [mesh.Add(MeshPoint(Pnt(0,0,0))), mesh.Add(MeshPoint(Pnt(1,0,0))), mesh.Add(MeshPoint(Pnt(1,1,0))), mesh.Add(MeshPoint(Pnt(0,1,0)))]
    mesh.Add(Element2D(1, pnts))
    return mesh

def Tet():
    mesh = Mesh(3)
    mesh.AddRegion("",2)
    mesh.AddRegion("",3)
    pnts = [mesh.Add(MeshPoint(Pnt(0,0,0))), mesh.Add(MeshPoint(Pnt(1,0,0))), mesh.Add(MeshPoint(Pnt(0,0,1))), mesh.Add(MeshPoint(Pnt(0,1,0)))]
    mesh.Add(Element2D(1, pnts[:3]))
    mesh.Add(Element2D(1, list(reversed(pnts[1:4]))))
    mesh.Add(Element2D(1, list(reversed(pnts[:2] + [pnts[3]]))))
    mesh.Add(Element2D(1, [pnts[0]] + pnts[2:4]))
    mesh.Add(Element3D(1, pnts))
    return mesh

def Prism():
    mesh = Mesh(3)
    mesh.AddRegion("",2)
    mesh.AddRegion("",3)
    pnts = [mesh.Add(MeshPoint(Pnt(0,0,0))), mesh.Add(MeshPoint(Pnt(1,0,0))), mesh.Add(MeshPoint(Pnt(0,1,0))),
            mesh.Add(MeshPoint(Pnt(0,0,1))), mesh.Add(MeshPoint(Pnt(1,0,1))), mesh.Add(MeshPoint(Pnt(0,1,1)))]
    mesh.Add(Element2D(1, list(reversed(pnts[:3]))))
    mesh.Add(Element2D(1, pnts[3:]))
    mesh.Add(Element2D(1, pnts[:2] + list(reversed(pnts[3:5]))))
    mesh.Add(Element2D(1, pnts[1:3] + list(reversed(pnts[4:6]))))
    mesh.Add(Element2D(1, [pnts[2], pnts[0], pnts[3], pnts[5]]))
    mesh.Add(Element3D(1, list(reversed(pnts[:3])) + list(reversed(pnts[3:]))))
    return mesh

if __name__ == "__main__":
    mesh = Prism()
    import ngsolve
    ngsolve.Draw(ngsolve.Mesh(mesh))
