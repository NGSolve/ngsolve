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

#test loading standard mesh
def test_simple():
    comm = MPI_Init()
    mesh = Mesh('square.vol.gz')

#test loading into sub-communicator
def test_sub1():
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
    
#test setting ngs_comm to sub-comm,
#then standard loading
def test_sub2():
    comm = MPI_Init()
    assert comm.size>=5
    groups = [[1,2,4]]
    sub_comm = comm.SubComm(find_group(comm, groups))
    from ngsolve.ngstd import SetNGSComm
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

#test loading standard mesh
def test_ngmesh_simple():
    comm = MPI_Init()
    import netgen
    ngmesh = netgen.meshing.Mesh(dim=2)
    ngmesh.Load('square.vol.gz')
    mesh = Mesh(ngmesh)

#test setting ngs_comm to sub-comm,
#then standard loading
def test_ngmesh_sub1():
    comm = MPI_Init()
    assert comm.size>=5
    groups = [[0,2,3,4]]
    sub_comm = comm.SubComm(find_group(comm, groups))
    from ngsolve.ngstd import SetNGSComm
    SetNGSComm(sub_comm)
    import netgen
    ngmesh = netgen.meshing.Mesh(dim=2)
    ngmesh.Load('square.vol.gz')
    # communicator now wrapped on netgen-side!
    # assert ngmesh.comm.size == sub_comm.size
    # assert ngmesh.comm.rank == sub_comm.rank
    mesh = Mesh(ngmesh)
    assert mesh.comm.size == sub_comm.size
    assert mesh.comm.rank == sub_comm.rank
    SetNGSComm(comm)


if __name__ == "__main__":
    comm = MPI_Init()
    #need at least NP=5 for this test to make sense
    print('\n------ simple ------- ')
    test_simple()
    print('\n------ sub1 ------- ')
    test_sub1()
    print('\n------ sub2 ------- ')
    test_sub2()
    print('\n------ ng_simple ------- ')
    test_ngmesh_simple()
    print('\n------ ng_sub ------- ')
    test_ngmesh_sub1()
