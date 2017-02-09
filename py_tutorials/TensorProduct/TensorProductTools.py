from netgen.meshing import *
from netgen.geom2d import unit_square
from netgen.csg import *
from ngsolve import *
import netgen.meshing as ngmeshing

x1 = CoordCF(3)
y1 = CoordCF(4)
z1 = CoordCF(5)

def MakeTensorProductMesh(mesh1,mesh2):
    if mesh1.dim + mesh2.dim >3:
        print('MakeMesh only possible if dim<=3!!')
        return
    if mesh1.dim + mesh2.dim == 3:
        tpmesh = MakeMesh3D(mesh1,mesh2)
    if mesh1.dim + mesh2.dim == 2:
        tpmesh = MakeMesh2D(mesh1,mesh2)
    AddSurfElements(tpmesh,mesh1,mesh2)
    return tpmesh

def AddSurfElements(tpmesh,mesh1,mesh2):
    if mesh1.dim + mesh2.dim == 3:
        tpmesh = AddSurfElements2D(tpmesh,mesh1,mesh2)
    if mesh1.dim + mesh2.dim == 2:
        tpmesh = AddSurfElements1D(tpmesh,mesh1,mesh2)
    return tpmesh

    
def MakeMesh2D(mesh1,mesh2):
    tpmesh = ngmeshing.Mesh(dim = mesh1.dim + mesh2.dim)
    ngm1 = mesh1.ngmesh
    ngm2 = mesh2.ngmesh

    els1 = ngm1.Elements1D()
    els2 = ngm2.Elements1D()
    vert1 = ngm1.Points()
    vert2 = ngm2.Points()
    pids = []
    for i in range(len(vert1)):
        for j in range(len(vert2)):
            pids.append(tpmesh.Add(MeshPoint(Pnt(vert1[ngmeshing.PointId(i+1)].p[0],vert2[ngmeshing.PointId(j+1)].p[0],0     ))))
    tpmesh.Add (FaceDescriptor(surfnr=1,domin=1,bc=1))
    for elx in els1:
        for ely in els2:
            pnum = [(elx.vertices[1].nr-1) * len(vert2) + ely.vertices[1].nr-1,
                    (elx.vertices[0].nr-1) * len(vert2) + ely.vertices[1].nr-1,
                    (elx.vertices[0].nr-1) * len(vert2) + ely.vertices[0].nr-1,
                    (elx.vertices[1].nr-1) * len(vert2) + ely.vertices[0].nr-1]
                    
            elpids = [pids[p] for p in pnum]
            tpmesh.Add(Element2D(1,elpids))
    return tpmesh

def MakeMesh3D(mesh1,mesh2):
    tpmesh = ngmeshing.Mesh(dim = mesh1.dim + mesh2.dim)
    ngm1 = mesh1.ngmesh
    ngm2 = mesh2.ngmesh
    if mesh1.dim==2:
        els1 = ngm1.Elements2D()
        els2 = ngm2.Elements1D()
        vert1 = ngm1.Points()
        vert2 = ngm2.Points()
        pids = []
        for i in range(len(vert1)):
            for j in range(len(vert2)):
                pids.append(tpmesh.Add(MeshPoint(Pnt(vert1[ngmeshing.PointId(i+1)].p[0],vert1[ngmeshing.PointId(i+1)].p[1],vert2[ngmeshing.PointId(j+1)].p[0]))))
        for elx in els1:
            for ely in els2:
                pnum = []
                for j in reversed(ely.vertices):
                    for i in elx.vertices:
                        pnum.append((i.nr-1) * len(vert2) + j.nr-1)
                elpids = [pids[p] for p in pnum]
                tpmesh.Add( Element3D(1,elpids) )
        return tpmesh
    else:
        els1 = ngm1.Elements1D()
        els2 = ngm2.Elements2D()
        vert1 = ngm1.Points()
        vert2 = ngm2.Points()
        pids = []
        for i in range(len(vert1)):
            for j in range(len(vert2)):
                pids.append(tpmesh.Add(MeshPoint(Pnt(vert1[ngmeshing.PointId(i+1)].p[0],vert2[ngmeshing.PointId(j+1)].p[0],vert2[ngmeshing.PointId(j+1)].p[1]))))
        for elx in els1:
            for ely in els2:
                pnum = []
                for i in reversed(elx.vertices):
                    for j in (ely.vertices):
                        pnum.append((i.nr-1) * len(vert2) + j.nr-1)
                elpids = [pids[p] for p in pnum]
                tpmesh.Add( Element3D(1,elpids) )
        return tpmesh

def AddSurfElements1D(tpmesh,mesh1,mesh2):
    ngm1 = mesh1.ngmesh;
    ngm2 = mesh2.ngmesh;
    els1 = ngm1.Elements1D()
    els2 = ngm2.Elements1D()
    for ely in els2:
        elpids = ely.vertices
        elpids1=[]
        for i in elpids:
            elpids1.append(PointId((i.nr-1)+(len(ngm1.Points()) -1 )*(len(ngm2.Points())) + 1 ))
        tpmesh.Add(Element1D(elpids))
        tpmesh.Add(Element1D(elpids1))
    for elx in els1:
        elpids = elx.vertices
        elpids1=[]
        elpids2=[]
        for i in elpids:
            elpids1.append(PointId( (i.nr-1)*(len(ngm2.Points()))+1) )
            elpids2.append(PointId( (i.nr-1)*(len(ngm2.Points()))+(len(ngm2.Points()))) )
        tpmesh.Add(Element1D(elpids1))
        tpmesh.Add(Element1D(elpids2))       

def AddSurfElements2D(tpmesh,mesh1,mesh2):
    ngm1 = mesh1.ngmesh;
    ngm2 = mesh2.ngmesh;
    if mesh1.dim==2:
        els1 = ngm1.Elements2D()
        els2 = ngm2.Elements1D()
        tpmesh.Add (FaceDescriptor(surfnr=1,domin=1,bc=1))
        for elx in els1:
            vert_loc = elx.vertices
            vert_glob = []
            for vx in vert_loc:
                vert_glob.append(PointId((vx.nr-1)*len(ngm2.Points())+len(ngm2.Points())))
            tpmesh.Add(Element2D(1,vert_glob))
        for elx in els1:
            vert_loc = elx.vertices
            vert_glob = []
            for vx in vert_loc:
                vert_glob.append(PointId((vx.nr-1)*len(ngm2.Points())+1))
            tpmesh.Add(Element2D(1,vert_glob))
        els1 = ngm1.Elements1D()
        for elx in els1:
            for ely in els2:
                vert_glob=[]
#            for vy in ely.vertices:
#                for vx in elx.vertices:
                vx = elx.vertices
                vy = ely.vertices
                vert_glob = [PointId((vx[1].nr-1)*len(ngm2.Points())+vy[0].nr),
                            PointId((vx[1].nr-1)*len(ngm2.Points())+vy[1].nr),
                            PointId((vx[0].nr-1)*len(ngm2.Points())+vy[1].nr),
                            PointId((vx[0].nr-1)*len(ngm2.Points())+vy[0].nr)]
                tpmesh.Add(Element2D(1,vert_glob))
    else:
        els1 = ngm1.Elements1D()
        els2 = ngm2.Elements2D()
        tpmesh.Add (FaceDescriptor(surfnr=1,domin=1,bc=1))
        for ely in els2:
            vert_loc = ely.vertices
            vert_glob = []
            for vy in vert_loc:
                vert_glob.append(PointId((vy.nr)+(len(ngm1.Points())-1)*(len(ngm2.Points()))))
            tpmesh.Add(Element2D(1,vert_glob))
        for ely in els2:
            vert_loc = ely.vertices
            vert_glob = []
            for vy in reversed(vert_loc):
                vert_glob.append(PointId( vy.nr))
            tpmesh.Add(Element2D(1,vert_glob))
        els2 = ngm2.Elements1D()
        for elx in els1:
            for ely in els2:
                vert_glob=[]
                vx = elx.vertices
                vy = ely.vertices
                vert_glob = [PointId((vx[0].nr-1)*len(ngm2.Points())+vy[0].nr),
                             PointId((vx[0].nr-1)*len(ngm2.Points())+vy[1].nr),
                             PointId((vx[1].nr-1)*len(ngm2.Points())+vy[1].nr),
                             PointId((vx[1].nr-1)*len(ngm2.Points())+vy[0].nr)]
                tpmesh.Add(Element2D(1,vert_glob))
    return tpmesh