from netgen.meshing import *
from netgen.csg import *
import ngsolve

def Make1DMesh(n, mapping = None, periodic=False):
    """
    Generate an equidistant 1D mesh with N cells

    Parameters
    ----------
    n : int
      Number of cells.

    mapping: lamda
      Mapping to transform the generated points. If None, the identity mapping is used.

    periodic: bool
      If True, the endpoints are identified to generate a periodic mesh.

    Returns
    -------
    (ngsolve.mesh)
      Returns generated 1D NGSolve mesh

    """
    mesh = Mesh(dim=1)
    pids = []
    for i in range(n+1):
        x = i/n
        if mapping:
            x = mapping(x)
        pids.append (mesh.Add (MeshPoint(Pnt(x, 0, 0))))

    idx_inner = mesh.AddRegion("dom", dim=1)
    idx_left = mesh.AddRegion("left", dim=0)
    idx_right = mesh.AddRegion("right", dim=0)
        
    for i in range(n):
        mesh.Add(Element1D([pids[i+1],pids[i]],index=idx_inner))
    mesh.Add (Element0D( pids[0], index=idx_left))
    mesh.Add (Element0D( pids[n], index=idx_right))
    if periodic:
        mesh.AddPointIdentification(pids[0],pids[n],1,2)
    ngsmesh = ngsolve.Mesh(mesh)
    return ngsmesh

def MakeStructured2DMesh(quads=True, nx=10, ny=10, secondorder=False, periodic_x=False, periodic_y=False, mapping = None, bbpts=None, bbnames=None, flip_triangles=False, boundarylayer=None, hppnts=None):
    """
    Generate a structured 2D mesh

    Parameters
    ----------
    quads : bool
      If True, a quadrilateral mesh is generated. If False, the quads are split to triangles.

    nx : int
      Number of cells in x-direction.

    ny : int
      Number of cells in y-direction.

    secondorder : bool
      If True, second order curved elements are used.

    periodic_x: bool
      If True, the left and right boundaries are identified to generate a periodic mesh in x-direction.

    periodic_y: bool
      If True, the top and bottom boundaries are identified to generate a periodic mesh in y-direction.

    mapping: lamda
      Mapping to transform the generated points. If None, the identity mapping is used.
    
    bbpts : list
      List of points which should be handled as BBND and are named with bbnames. The mesh (nx, ny and mapping) must be constructed in such a way that the bbpts coincide with generated points. Otherwise an Exception is thrown.

    bbnames : list
      List of bbnd names as strings. Size must coincide with size of bbpts. Otherwise an Exception is thrown.

    flip_triangles : bool
      If set to True together with quads=False the quads are cut the other way round

    boundarylayer : dict
      If not None it expects a dictionary of the form { "boundaryname" : [t1,...,tn] } where ti denote the thickness of layer i. The number of layers are included in nx/ny. After the layers are placed the remaining number of cells are used to divide the remaining grid uniformly.

    hppnts : list
      If not None it expects a list of the form [ (px1,py1, hpref1), (px2,py2, hpref2), ... ] where px,py are the point coordinates which have to be resolved in the mesh and hpref the refinement factor

    Returns
    -------
    (ngsolve.mesh)
      Returns generated 2D NGSolve mesh

    """
    mesh = Mesh()
    mesh.dim=2

    if (bbpts and bbnames) and len(bbpts) != len(bbnames):
        raise Exception("Lenght of bbnames does not coincide with length of bbpts!")

    found = []
    indbbpts = []
    if bbpts:
        for i in range(len(bbpts)):
            found.append(False)
            indbbpts.append(None)
    foundhp = [ False for i in hppnts] if hppnts else []
        

    pids = []
    if periodic_y:
        minioni = []
        masteri = []
    if periodic_x:        
        minionj = []
        masterj = []
        
    numlayerleft  = len(boundarylayer["left"]) if (boundarylayer and boundarylayer.get("left")) else 0
    numlayerright = len(boundarylayer["right"]) if (boundarylayer and boundarylayer.get("right")) else 0
    numlayertop   = len(boundarylayer["top"]) if (boundarylayer and boundarylayer.get("top")) else 0
    numlayerbot   = len(boundarylayer["bottom"]) if (boundarylayer and boundarylayer.get("bottom")) else 0

    thicknessleft  = [0]
    thicknessright = [0]
    thicknesstop   = [0]
    thicknessbot   = [0]
    for i in range(numlayerleft):
        thicknessleft.append(thicknessleft[-1]+boundarylayer["left"][i])
    for i in range(numlayerright):
        thicknessright.append(thicknessright[-1]+boundarylayer["right"][i])
    for i in range(numlayertop):
        thicknesstop.append(thicknesstop[-1]+boundarylayer["top"][i])
    for i in range(numlayerbot):
        thicknessbot.append(thicknessbot[-1]+boundarylayer["bottom"][i])

        
    for i in range(ny+1):
        for j in range(nx+1):
            x = thicknessleft[j] if j < numlayerleft else ((thicknessleft[-1]+(j-numlayerleft)/(nx-numlayerleft-numlayerright)*(1-thicknessleft[-1]-thicknessright[-1])) if j < nx-numlayerright else 1-thicknessright[nx-j])
            y = thicknessbot[i] if i < numlayerbot else ((thicknessbot[-1]+(i-numlayerbot)/(ny-numlayerbot-numlayertop)*(1-thicknessbot[-1]-thicknesstop[-1])) if i< ny-numlayertop else 1-thicknesstop[ny-i])
            pids.append(mesh.Add (MeshPoint(Pnt(x,y,0))))
            if periodic_y:
                if i == 0:
                    minioni.append(pids[-1])
                if i == ny:
                    masteri.append(pids[-1])  
            if periodic_x:                       
                if j == 0:
                    minionj.append(pids[-1])
                if j == nx:
                    masterj.append(pids[-1])        
    if periodic_y:
        for i in range(len(minioni)):   
            mesh.AddPointIdentification(masteri[i],minioni[i],identnr=1,type=2)
    if periodic_x:            
        for j in range(len(minionj)):        
            mesh.AddPointIdentification(masterj[j],minionj[j],identnr=2,type=2)                                       

    # mesh.Add(FaceDescriptor(surfnr=1,domin=1,bc=1))
    idx_dom = mesh.AddRegion("dom", dim=2)
    idx_bottom = mesh.AddRegion("bottom", dim=1)
    idx_right  = mesh.AddRegion("right", dim=1)
    idx_top    = mesh.AddRegion("top", dim=1)
    idx_left   = mesh.AddRegion("left", dim=1)
    
    for i in range(ny):
        for j in range(nx):
            base = i * (nx+1) + j
            if quads:
                pnum = [base,base+1,base+nx+2,base+nx+1]
                elpids = [pids[p] for p in pnum]
                el = Element2D(idx_dom,elpids)
                if not mapping:
                    el.curved=False
                mesh.Add(el)
            else:
                if flip_triangles:
                    pnum1 = [base,base+1,base+nx+2]
                    pnum2 = [base,base+nx+2,base+nx+1]
                else:
                    pnum1 = [base,base+1,base+nx+1]
                    pnum2 = [base+1,base+nx+2,base+nx+1]
                elpids1 = [pids[p] for p in pnum1]
                elpids2 = [pids[p] for p in pnum2]
                mesh.Add(Element2D(idx_dom,elpids1)) 
                mesh.Add(Element2D(idx_dom,elpids2))                          

    for i in range(nx):
        mesh.Add(Element1D([pids[i], pids[i+1]], index=idx_bottom))
    for i in range(ny):
        mesh.Add(Element1D([pids[i*(nx+1)+nx], pids[(i+1)*(nx+1)+nx]], index=idx_right))
    for i in range(nx):
        mesh.Add(Element1D([pids[ny*(nx+1)+i+1], pids[ny*(nx+1)+i]], index=idx_top))
    for i in range(ny):
        mesh.Add(Element1D([pids[(i+1)*(nx+1)], pids[i*(nx+1)]], index=idx_left))

    # mesh.SetBCName(0, "bottom")        
    # mesh.SetBCName(1, "right")        
    # mesh.SetBCName(2, "top")        
    # mesh.SetBCName(3, "left")  

    mesh.Compress()       
    
    if secondorder:
        mesh.SecondOrder()
    
    if mapping:
        for p in mesh.Points():
            x,y,z = p.p
            x,y = mapping(x,y)
            p[0] = x
            p[1] = y

    for k in range(len(found)):
        i = 0
        for p in mesh.Points():
            if abs(p.p[0]-bbpts[k][0])+abs(p.p[1]-bbpts[k][1]) < 1e-6:
                indbbpts[k] = pids[i]
                found[k] = True
            i += 1
    for k in range(len(found)):
        if found[k] == False:
            raise Exception("bbpnt[",k,"] not in structured mesh!")
    for i in range(len(indbbpts)):
        mesh.Add(Element0D(indbbpts[i], index=i+1))
        mesh.SetCD2Name(i+1, bbnames[i])

    for k in range(len(foundhp)):
        i = 0
        for p in mesh.Points():
            if abs(p.p[0]-hppnts[k][0])+abs(p.p[1]-hppnts[k][1]) < 1e-6:
                mesh.AddSingularity(pids[i],hppnts[k][-1])
                foundhp[k] = True
            i += 1
    for k in range(len(foundhp)):
        if foundhp[k] == False:
            raise Exception("hppnts[",k,"] not in structured mesh!")
    
            
    ngsmesh = ngsolve.Mesh(mesh)
    return ngsmesh

def MakeQuadMesh(nx=10, ny=10, periodic_x=False, periodic_y=False, mapping = None):
    """
    Generate a structured quadrilateral 2D mesh

    Parameters
    ----------
    nx : int
      Number of cells in x-direction.

    ny : int
      Number of cells in y-direction.

    periodic_x: bool
      If True, the left and right boundaries are identified to generate a periodic mesh in x-direction.

    periodic_y: bool
      If True, the top and bottom boundaries are identified to generate a periodic mesh in y-direction.

    mapping: lambda
      Mapping to transform the generated points. If None, the identity mapping is used.
    

    Returns
    -------
    (ngsolve.mesh)
      Returns generated 2D NGSolve mesh

    """
    return MakeStructured2DMesh(quads=True, nx=nx, ny=ny, periodic_x=periodic_x, periodic_y=periodic_y, mapping=mapping)    

def MakeStructured3DMesh(hexes=True, nx=10, ny=None, nz=None, secondorder=False, periodic_x=False, periodic_y=False, periodic_z=False, mapping = None, cuboid_mapping=False, prism=False):
    """
    Generate a structured quadrilateral 2D mesh

    Parameters
    ----------
    hexes: bool
      If True, a mesh consisting of hexahedra is generated.

    nx : int
      Number of cells in x-direction.

    ny : int
      Number of cells in y-direction.

    nz : int
      Number of cells in z-direction.

    secondorder : bool
      If True, second order curved elements are used.
 
    periodic_x: bool
      If True, the left and right boundaries are identified to generate a periodic mesh in x-direction.

    periodic_y: bool
      If True, the top and bottom boundaries are identified to generate a periodic mesh in y-direction.

    periodic_z: bool
      If True, the top and bottom boundaries are identified to generate a periodic mesh in z-direction.

    mapping: lambda
      Mapping to transform the generated points. If None, the identity mapping is used.
    
    cuboid_mapping: bool
      If True, a straight geometry is assumed.

    prism : bool
      If True, a mesh consisting of prism is generated. If hexes and prism is set to True, also a mesh consisting of prism is generated.


    Returns
    -------
    (ngsolve.mesh)
      Returns generated 3D NGSolve mesh

    """
    if nz == None:
        if ny == None:
            nz = nx
        else:
            raise Exception("MakeStructured3DMesh: No default value for nz if nx and ny are provided")
    if ny == None:
        ny = nx
        
    netmesh = Mesh()
    netmesh.dim = 3

    if cuboid_mapping:
        P1 = (0,0,0)
        P2 = (1,1,1)
        if mapping:
            P1 = mapping(*P1)
            P2 = mapping(*P2)
        cube = OrthoBrick(Pnt(P1[0], P1[1], P1[2]), Pnt(P2[0], P2[1], P2[2])).bc(1)
        geom = CSGeometry()
        geom.Add(cube)
        netmesh.SetGeometry(geom)

    pids = []
    if periodic_x:
        minioni = []
        masteri = []
    if periodic_y:        
        minionj = []
        masterj = []
    if periodic_z:        
        minionk = []
        masterk = []        
    for i in range(nx+1):
        for j in range(ny+1):
            for k in range(nz+1):
                # x,y,z = mapping(i / nx, j / ny, k / nz)
                x,y,z = i / nx, j / ny, k / nz
                # if mapping:
                #   x,y,z = mapping(x,y,z)
                pids.append(netmesh.Add(MeshPoint(Pnt( x,y,z ))))
                if periodic_x:
                    if i == 0:
                        minioni.append(pids[-1])
                    if i == nx:
                        masteri.append(pids[-1])  
                if periodic_y:           
                    if j == 0:
                        minionj.append(pids[-1])
                    if j == ny:
                        masterj.append(pids[-1]) 
                if periodic_z:                    
                    if k == 0:
                        minionk.append(pids[-1])
                    if k == nz:
                        masterk.append(pids[-1])
    if periodic_x:
        for i in range(len(minioni)):   
            netmesh.AddPointIdentification(masteri[i],minioni[i],identnr=1,type=2)     
    if periodic_y:        
        for j in range(len(minionj)):            
            netmesh.AddPointIdentification(masterj[j],minionj[j],identnr=2,type=2) 
    if periodic_z:        
        for k in range(len(minionk)):            
            netmesh.AddPointIdentification(masterk[k],minionk[k],identnr=3,type=2)                                                      

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                base = i * (ny+1)*(nz+1) + j*(nz+1) + k
                baseup = base+(ny+1)*(nz+1)
                pnum = [base, base+1, base+(nz+1)+1, base+(nz+1),
                        baseup, baseup+1, baseup+(nz+1)+1, baseup+(nz+1)]
                if prism:
                    for qarr in [[0,3,4,1,2,5],[3,7,4,2,6,5]]:
                        elpids = [pids[p] for p in [pnum[q] for q in qarr]]
                        netmesh.Add(Element3D(1, elpids))
                elif hexes:
                    elpids = [pids[p] for p in pnum]
                    el = Element3D(1, elpids)
                    if not mapping:
                        el.curved = False
                    netmesh.Add(el)
                else:
                    #  a poor mans kuhn triangulation of a cube
                    for qarr in [[0, 4, 5, 6],
                                 [0, 6, 7, 4],
                                 [0, 3, 7, 6],
                                 [0, 1, 6, 5],
                                 [0, 1, 2, 6],
                                 [0, 3, 6, 2]]:
                        elpids = [pids[p] for p in [pnum[q] for q in qarr]]
                        netmesh.Add(Element3D(1, elpids))

    def AddSurfEls(p1, dxi, nxi, deta, neta, facenr):
        def add_seg(i, j, os):
            base = p1 + i*dxi + j*deta
            pnum = [base, base+os]
            elpids = [pids[p] for p in pnum]
            netmesh.Add(Element1D(elpids, index=facenr))
        for i in range(nxi):
            for j in [0,neta]:
                add_seg(i,j,dxi)
        for i in [0,nxi]:
            for j in range(neta):
                add_seg(i,j,deta)
        for i in range(nxi):
            for j in range(neta):
                base = p1 + i*dxi+j*deta
                pnum = [base, base+dxi, base+dxi+deta, base+deta]
                if prism:
                    if facenr <= 4:
                        qarr = [1,2,3,0]
                        elpids = [pids[pnum[p]] for p in qarr]
                        netmesh.Add(Element2D(facenr, elpids))
                    else:
                        qarrs = [[3, 1, 2], [3, 0, 1]]
                        for qarr in qarrs:
                            elpids = [pids[p] for p in [pnum[q] for q in qarr]]
                            netmesh.Add(Element2D(facenr, elpids))
                elif hexes:
                    elpids = [pids[p] for p in pnum]
                    netmesh.Add(Element2D(facenr, elpids))
                else:
                    qarrs = [[0, 1, 2], [0, 2, 3]]
                    for qarr in qarrs:
                        elpids = [pids[p] for p in [pnum[q] for q in qarr]]
                        netmesh.Add(Element2D(facenr, elpids))

    #order is important!
    netmesh.Add(FaceDescriptor(surfnr=4, domin=1, bc=1))
    netmesh.Add(FaceDescriptor(surfnr=2, domin=1, bc=2))
    netmesh.Add(FaceDescriptor(surfnr=5, domin=1, bc=3))
    netmesh.Add(FaceDescriptor(surfnr=3, domin=1, bc=4))
    netmesh.Add(FaceDescriptor(surfnr=0, domin=1, bc=5))
    netmesh.Add(FaceDescriptor(surfnr=1, domin=1, bc=6))
        
    # y-z-plane, smallest x-coord: ("back")
    AddSurfEls(0, 1, nz,  nz+1, ny, facenr=1) # y-z-plane
    # x-z-plane, smallest y-coord: ("left")
    AddSurfEls(0, (ny+1)*(nz+1), nx, 1, nz,facenr=2)
    # y-z-plane, largest x-coord: ("front")
    AddSurfEls((nx+1)*(ny+1)*(nz+1)-1, -(nz+1), ny, -1, nz, facenr=3) 
    # x-z-plane, largest y-coord: ("right")
    AddSurfEls((nx+1)*(ny+1)*(nz+1)-1, -1, nz, -(ny+1)*(nz+1), nx, facenr=4)
    # x-y-plane, smallest z-coord: ("bottom")
    AddSurfEls(0, nz+1, ny, (ny+1)*(nz+1), nx,facenr=5) 
    # x-y-plane, largest z-coord: ("top")
    AddSurfEls((nx+1)*(ny+1)*(nz+1)-1, -(ny+1)*(nz+1), nx, -(nz+1), ny, facenr=6) 

    netmesh.SetBCName(0,"back")
    netmesh.SetBCName(1,"left")
    netmesh.SetBCName(2,"front")
    netmesh.SetBCName(3,"right")
    netmesh.SetBCName(4,"bottom")
    netmesh.SetBCName(5,"top")
    
    netmesh.Compress()

    if secondorder:
        netmesh.SecondOrder()
    
    if mapping:
        for p in netmesh.Points():
            x,y,z = p.p
            x,y,z = mapping(x,y,z)
            p[0] = x
            p[1] = y
            p[2] = z
            
    ngsmesh = ngsolve.Mesh(netmesh)
    # ngsmesh.ngmesh.Save("tmp.vol.gz")
    # ngsmesh = ngsolve.Mesh("tmp.vol.gz")
    return ngsmesh

def MakeHexMesh(nx=10, ny=10, nz=10, secondorder=False, periodic_x=False, periodic_y=False, periodic_z=False, mapping = None, cuboid_mapping=False):
    """
    Generate a structured hexahedra 3D mesh

    Parameters
    ----------
    nx : int
      Number of cells in x-direction.

    ny : int
      Number of cells in y-direction.

    nz : int
      Number of cells in z-direction.

    periodic_x: bool
      If True, the left and right boundaries are identified to generate a periodic mesh in x-direction.

    periodic_y: bool
      If True, the top and bottom boundaries are identified to generate a periodic mesh in y-direction.

    periodic_z: bool
      If True, the top and bottom boundaries are identified to generate a periodic mesh in z-direction.

    mapping: lambda
      Mapping to transform the generated points. If None, the identity mapping is used.
    
    cuboid_mapping: bool
      If True, a straight geometry is assumed.


    Returns
    -------
    (ngsolve.mesh)
      Returns generated 3D NGSolve mesh consisting of only hexahedra

    """
    return MakeStructured3DMesh(hexes=True, nx=nx, ny=ny, nz=nz, secondorder=secondorder, periodic_x=periodic_x, periodic_y=periodic_y, periodic_z=periodic_z, mapping=mapping, cuboid_mapping=cuboid_mapping)

def MakePrismMesh(nx=10, ny=None, nz=None, secondorder=False, periodic_x=False, periodic_y=False, periodic_z=False, mapping = None, cuboid_mapping=False):
    """
    Generate a structured prism 3D mesh

    Parameters
    ----------
    nx : int
      Number of cells in x-direction.

    ny : int
      Number of cells in y-direction.

    nz : int
      Number of cells in z-direction.

    periodic_x: bool
      If True, the left and right boundaries are identified to generate a periodic mesh in x-direction.

    periodic_y: bool
      If True, the top and bottom boundaries are identified to generate a periodic mesh in y-direction.

    periodic_z: bool
      If True, the top and bottom boundaries are identified to generate a periodic mesh in z-direction.

    mapping: lambda
      Mapping to transform the generated points. If None, the identity mapping is used.
    
    cuboid_mapping: bool
      If True, a straight geometry is assumed.


    Returns
    -------
    (ngsolve.mesh)
      Returns generated 3D NGSolve mesh consisting of only prism

    """
    return MakeStructured3DMesh(hexes=False, nx=nx, ny=ny, nz=nz, secondorder=secondorder, periodic_x=periodic_x, periodic_y=periodic_y, periodic_z=periodic_z, mapping=mapping, cuboid_mapping=cuboid_mapping, prism=True)

from math import pi
from ngsolve import Draw, sin, cos
if __name__ == "__main__":

    mesh = Make1DMesh(n=4)
    print("simple 1D mesh -- no visualization -- ")

    mesh = MakeQuadMesh(nx=4, ny=4)
    Draw(mesh)
    input("simple quad mesh -- press any key to continue -- ")

    mesh = MakeStructured2DMesh(quads=False, nx=4, ny=4)
    Draw(mesh)
    input("simple trig mesh -- press any key to continue -- ")    

    mesh = MakeQuadMesh(nx=4, ny=4, periodic_x=True, periodic_y=False)
    Draw(mesh)
    input("x-periodic quad mesh -- press any key to continue -- ") 

    mesh = MakeStructured2DMesh(quads=False, nx=16, ny=4, periodic_x=False, periodic_y=True)
    Draw(mesh)
    input("y-periodic non-uniform trig mesh -- press any key to continue -- ")        
    
    mesh = MakeStructured2DMesh(quads=True, nx=32, ny=16,
                            mapping = lambda x,y : ((1.1+sin(pi*(y-0.5)))*sin(pi*x),
                                                (1.1+sin(pi*(y-0.5)))*cos(pi*x)))
    Draw(mesh)
    input("structured quad half circle with a whole -- press any key to continue -- ")
    
    mesh = MakeHexMesh(nx=8)
    Draw(mesh)
    input("simple cube mesh -- press any key to continue -- ")

    mesh = MakeHexMesh(nx=8, periodic_x=True, periodic_y=True, periodic_z=True)
    Draw(mesh)
    input("periodic cube mesh -- press any key to continue -- ")    
    
    mesh = MakeStructured3DMesh(hexes=False, nx=3, ny=6, nz=10,
                            mapping = lambda x,y,z : (x,0.5*y*(y+x),exp(z)),
                            cuboid_mapping=False )
    Draw(mesh)
    input("mapped, anisotropic linear non-cuboid mesh -- press any key to continue -- ")
    
    mesh = MakeStructured3DMesh(hexes=True, nx=8, ny=16, nz=8, periodic_x=True,
                            mapping = lambda x,y,z : (5*x*x*(0.5-x/3),10*y*y*(0.5-y/3),5*z*z*(0.5-z/3)),
                            cuboid_mapping=True )
    Draw(mesh)
    input("x-periodic, mapped, anisotropic non-linear cuboid mesh -- press any key to continue --")
    
    mesh = MakeStructured3DMesh(hexes=False, nx=8, ny=8, nz=8, periodic_x=True,
                            cuboid_mapping=True )
    mesh.Refine()
    Draw(mesh)
    input("structured tet mesh (periodic in x-dir.), refined -- finished.")





def MakeStructuredSurfaceMesh(quads=True, nx=10, ny=10, mapping = None, secondorder=False, bbbpts=None, bbbnames=None, flip_triangles=False):
    """
    Generate a structured 2D surface mesh in 3D

    Parameters
    ----------
    quads : bool
      If True, a quadrilateral mesh is generated. If False, the quads are split to triangles.

    nx : int
      Number of cells in x-direction.

    ny : int
      Number of cells in y-direction.

    mapping : lamda
      Mapping to transform the generated points. If None, the identity mapping is used.
    
    secondorder : bool
      If True, use quadratic elements, else linear elements are used.

    bbbpts : list
      List of points which should be handled as BBBND and are named with bbbnames. The mesh (nx, ny and mapping) must be constructed in such a way that the bbbpts coincide with generated points. Otherwise an Exception is thrown.

    bbbnames : list
      List of bbbnd names as strings. Size must coincide with size of bbbpts. Otherwise an Exception is thrown.

    flip_triangles : bool
      If set tot True together with quads=False the quads are cut the other way round

    Returns
    -------
    (ngsolve.mesh)
      Returns generated NGSolve mesh

    """
    mesh = Mesh(dim=3)

    if (bbbpts and bbbnames) and len(bbbpts) != len(bbbnames):
        raise Exception("Lenght of bbbnames does not coincide with length of bbbpts!")

    found = []
    indbbbpts = []
    if bbbpts:
        for i in range(len(bbbpts)):
            found.append(False)
            indbbbpts.append(None)

    pids = []
    for i in range(ny+1):
        for j in range(nx+1):
            x,y,z = j/nx, i/ny, 0
            pids.append(mesh.Add (MeshPoint(Pnt(x,y,z))))
            
    mesh.Add(FaceDescriptor(surfnr=1,domin=1,bc=1))
    
    for i in range(ny):
        for j in range(nx):
            base = i * (nx+1) + j
            if quads:
                pnum = [base,base+1,base+nx+2,base+nx+1]
                elpids = [pids[p] for p in pnum]
                el = Element2D(1,elpids)
                if not mapping:
                    el.curved=False
                mesh.Add(el)
            else:
                if flip_triangles:
                    pnum1 = [base,base+1,base+nx+2]
                    pnum2 = [base,base+nx+2,base+nx+1]
                else:
                    pnum1 = [base,base+1,base+nx+1]
                    pnum2 = [base+1,base+nx+2,base+nx+1]
                elpids1 = [pids[p] for p in pnum1]
                elpids2 = [pids[p] for p in pnum2]
                mesh.Add(Element2D(1,elpids1)) 
                mesh.Add(Element2D(1,elpids2))                          

    for i in range(nx):
        mesh.Add(Element1D([pids[i], pids[i+1]], index=1))
    for i in range(ny):
        mesh.Add(Element1D([pids[i*(nx+1)+nx], pids[(i+1)*(nx+1)+nx]], index=2))
    for i in range(nx):
        mesh.Add(Element1D([pids[ny*(nx+1)+i+1], pids[ny*(nx+1)+i]], index=3))
    for i in range(ny):
        mesh.Add(Element1D([pids[(i+1)*(nx+1)], pids[i*(nx+1)]], index=4))

    mesh.SetCD2Name(1, "bottom")        
    mesh.SetCD2Name(2, "right")        
    mesh.SetCD2Name(3, "top")        
    mesh.SetCD2Name(4, "left")

    if secondorder:
        mesh.SecondOrder()
    
    if mapping:
        for p in mesh.Points():
            x,y,z = p.p
            x,y,z = mapping(x,y,z)
            p[0],p[1],p[2] = x,y,z
           
    for k in range(len(found)):
        i = 0
        for p in mesh.Points():
            if abs(p.p[0]-bbbpts[k][0])+abs(p.p[1]-bbbpts[k][1])+abs(p.p[2]-bbbpts[k][2]) < 1e-6:
                indbbbpts[k] = pids[i]
                found[k] = True
            i += 1
    for k in range(len(found)):
        if found[k] == False:
            raise Exception("bbbpnt[",k,"] not in structured mesh!")

    for i in range(len(indbbbpts)):
        mesh.Add(Element0D(indbbbpts[i], index=i+1))
        mesh.SetCD3Name(i+1, bbbnames[i])

    mesh.Compress()       
    ngsmesh = ngsolve.Mesh(mesh)
    return ngsmesh
