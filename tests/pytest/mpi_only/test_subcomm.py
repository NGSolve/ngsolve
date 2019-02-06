from ngsolve import *

def find_group(comm, groups):
    mg = []
    for g in groups:
        if comm.rank in g:
            mg = g
            break
    if not len(mg):
        mg = [p for p in range(comm.size) if not p in set(p for g in groups for p in g)]
    if not len(mg) or len(mg)==2:
        mg = [comm.rank]
    return mg

#load ngsolve-mesh
def test_load_ngs():
    comm = MPI_Init()
    mesh = Mesh('square.vol.gz')
    comm.Barrier()
    
#load ngsolve-mesh into sub-communicator
def test_load_ngs_sub1():
    comm = MPI_Init()
    assert comm.size>=5
    groups = [[1,2,4]]
    mygrp = find_group(comm, groups)
    sub_comm = comm.SubComm(mygrp)
    assert sub_comm.size == len(mygrp)
    assert sub_comm.rank == mygrp.index(comm.rank)
    mesh = Mesh('square.vol.gz', sub_comm)
    assert mesh.comm.size == sub_comm.size
    assert mesh.comm.rank == sub_comm.rank
    comm.Barrier()

#set global comm to sub-comm, then load ngsolve-mesh
def test_load_ngs_sub2():
    comm = MPI_Init()
    assert comm.size>=5
    groups = [[1,2,4]]
    sub_comm = comm.SubComm(find_group(comm, groups))
    from ngsolve.comp import SetNGSComm
    SetNGSComm(sub_comm)
    newcomm = MPI_Init()
    assert newcomm.size == sub_comm.size
    assert newcomm.rank == sub_comm.rank
    mesh = Mesh('square.vol.gz')
    assert mesh.comm.size == sub_comm.size
    assert mesh.comm.rank == sub_comm.rank
    SetNGSComm(comm)
    newnewcomm = MPI_Init()
    assert newnewcomm.size == comm.size
    assert newnewcomm.rank == comm.rank
    comm.Barrier()

    
#load netgen-mesh, then make ngsolve-mesh
def test_load_ng():
    comm = MPI_Init()
    import netgen
    ngmesh = netgen.meshing.Mesh(dim=2)
    ngmesh.Load('square.vol.gz')
    mesh = Mesh(ngmesh)
    comm.Barrier()

#load netgen-mesh into sub-comm 
def test_load_ng_sub1():
    comm = MPI_Init()
    assert comm.size>=5
    groups = [[0,2,3,4]]
    sub_comm = comm.SubComm(find_group(comm, groups))
    import netgen
    ngmesh = netgen.meshing.Mesh(dim=2, comm=sub_comm)
    ngmesh.Load('square.vol.gz')
    # communicator now wrapped on netgen-side!
    assert ngmesh.comm.size == sub_comm.size
    assert ngmesh.comm.rank == sub_comm.rank
    mesh = Mesh(ngmesh)
    assert mesh.comm.size == sub_comm.size
    assert mesh.comm.rank == sub_comm.rank
    comm.Barrier()

#set global comm to sub-comm, then load netgen-mesh
def test_load_ng_sub2():
    comm = MPI_Init()
    assert comm.size>=5
    groups = [[0,2,3,4]]
    sub_comm = comm.SubComm(find_group(comm, groups))
    from ngsolve.comp import SetNGSComm
    SetNGSComm(sub_comm)
    import netgen
    ngmesh = netgen.meshing.Mesh(dim=2)
    ngmesh.Load('square.vol.gz')
    # communicator now wrapped on netgen-side!
    assert ngmesh.comm.size == sub_comm.size
    assert ngmesh.comm.rank == sub_comm.rank
    mesh = Mesh(ngmesh)
    assert mesh.comm.size == sub_comm.size
    assert mesh.comm.rank == sub_comm.rank
    SetNGSComm(comm)
    comm.Barrier()

def test_load_dist():
    comm = MPI_Init()
    mecomm = comm.SubComm([comm.rank])
    import netgen
    ngmesh = netgen.meshing.Mesh(dim=2, comm=mecomm)
    if comm.rank == 0:
        ngmesh.Load('square.vol.gz')
        assert ngmesh.comm.rank==0
        assert ngmesh.comm.size==1
    ngmesh.Distribute(comm)
    assert ngmesh.comm.rank==comm.rank
    assert ngmesh.comm.size==comm.size
    comm.Barrier()

def test_load_dist_sub():
    comm = MPI_Init()
    groups = [[2,3,4]]
    sub_comm = comm.SubComm(find_group(comm, groups))
    mecomm = comm.SubComm([comm.rank])
    import netgen
    ngmesh = netgen.meshing.Mesh(dim=2, comm=mecomm)
    if sub_comm.rank == 0:
        ngmesh.Load('square.vol.gz')
        assert ngmesh.comm.rank==0
        assert ngmesh.comm.size==1
    ngmesh.Distribute(sub_comm)
    assert ngmesh.comm.rank==sub_comm.rank
    assert ngmesh.comm.size==sub_comm.size
    comm.Barrier()

    
#mesh on master and then distribute
def test_mesh_dist():
    comm = MPI_Init()
    import netgen
    ngmesh2d = netgen.meshing.Mesh(dim=2)
    if comm.rank==0:
        from netgen.geom2d import unit_square
        ngmesh2d = unit_square.GenerateMesh(maxh=0.1)
    ngmesh2d.Distribute()
    mesh2d = Mesh(ngmesh2d)
    ngmesh3d = netgen.meshing.Mesh(dim=3)
    if comm.rank==0:
        from netgen.csg import unit_cube
        ngmesh3d = unit_cube.GenerateMesh(maxh=0.2)
    ngmesh3d.Distribute()
    mesh3d = Mesh(ngmesh3d)
    comm.Barrier()
    
#mesh on master and then distribute into subcomm
def test_mesh_dist_sub1():
    comm = MPI_Init()
    groups = [[0,1,3,4]]
    assert comm.size>=5
    sub_comm = comm.SubComm(find_group(comm, groups))
    import netgen
    ngmesh2d = netgen.meshing.Mesh(dim=2)
    if sub_comm.rank==0:
        from netgen.geom2d import unit_square
        ngmesh2d = unit_square.GenerateMesh(maxh=0.1)
    ngmesh2d.Distribute(sub_comm)
    assert ngmesh2d.comm.size == sub_comm.size
    assert ngmesh2d.comm.rank == sub_comm.rank
    mesh2d = Mesh(ngmesh2d)
    assert mesh2d.comm.size == sub_comm.size
    assert mesh2d.comm.rank == sub_comm.rank
    ngmesh3d = netgen.meshing.Mesh(dim=3)
    if sub_comm.rank==0:
        from netgen.csg import unit_cube
        ngmesh3d = unit_cube.GenerateMesh(maxh=0.2)
    ngmesh3d.Distribute(sub_comm)
    assert ngmesh3d.comm.size == sub_comm.size
    assert ngmesh3d.comm.rank == sub_comm.rank
    mesh3d = Mesh(ngmesh3d)
    assert mesh3d.comm.size == sub_comm.size
    assert mesh3d.comm.rank == sub_comm.rank
    comm.Barrier()

#set global comm, then mesh on master and then distribute
def test_mesh_dist_sub2():
    comm = MPI_Init()
    groups = [[1,2,3]]
    assert comm.size>=5
    sub_comm = comm.SubComm(find_group(comm, groups))
    from ngsolve.comp import SetNGSComm
    SetNGSComm(sub_comm)
    import netgen
    ngmesh2d = netgen.meshing.Mesh(dim=2)
    if sub_comm.rank==0:
        from netgen.geom2d import unit_square
        ngmesh2d = unit_square.GenerateMesh(maxh=0.1)
    ngmesh2d.Distribute()
    assert ngmesh2d.comm.size == sub_comm.size
    assert ngmesh2d.comm.rank == sub_comm.rank
    mesh2d = Mesh(ngmesh2d)
    assert mesh2d.comm.size == sub_comm.size
    assert mesh2d.comm.rank == sub_comm.rank
    ngmesh3d = netgen.meshing.Mesh(dim=3)
    if sub_comm.rank==0:
        from netgen.csg import unit_cube
        ngmesh3d = unit_cube.GenerateMesh(maxh=0.2)
    ngmesh3d.Distribute()
    assert ngmesh3d.comm.size == sub_comm.size
    assert ngmesh3d.comm.rank == sub_comm.rank
    mesh3d = Mesh(ngmesh3d)
    assert mesh3d.comm.size == sub_comm.size
    assert mesh3d.comm.rank == sub_comm.rank
    SetNGSComm(comm)
    comm.Barrier()

    

if __name__ == "__main__":
    comm = MPI_Init()
    #need at least NP=5 for this test to make sense

    test_load_ngs()
    test_load_ngs_sub1()
    test_load_ngs_sub2()

    test_load_ng()
    test_load_ng_sub1()
    test_load_ng_sub2()

    test_load_dist()
    test_load_dist_sub()

    test_mesh_dist()
    test_mesh_dist_sub1()
    test_mesh_dist_sub2()
