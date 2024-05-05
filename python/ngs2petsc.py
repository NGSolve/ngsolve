import ngsolve as ngs
import netgen.meshing as ngm
import petsc4py.PETSc as psc
from mpi4py import MPI
import numpy as np

#PETSc Matrix

def CreatePETScMatrix (ngs_mat, freedofs=None):
    locmat = ngs_mat.local_mat
    eh, ew = locmat.entrysizes
    if eh != ew: raise Exception ("only square entries are allowed")
    
    val,col,ind = locmat.CSR()
    ind = np.array(ind).astype(psc.IntType)
    col = np.array(col).astype(psc.IntType)    
    apsc_loc = psc.Mat().createBAIJ(size=(eh*locmat.height, eh*locmat.width), bsize=eh, csr=(ind,col,val), comm=MPI.COMM_SELF)

    if freedofs is not None:
        locfree = np.flatnonzero(freedofs).astype(psc.IntType)
        isfree_loc = psc.IS().createBlock(indices=locfree, bsize=eh)
        apsc_loc = apsc_loc.createSubMatrices(isfree_loc)[0]

    comm = MPI.COMM_WORLD
    if comm.Get_size() > 1: 
        pardofs = ngs_mat.row_pardofs
        comm = pardofs.comm.mpi4py
        globnums, nglob = pardofs.EnumerateGlobally(freedofs)
        if freedofs is not None:
            globnums = np.array(globnums, dtype=psc.IntType)[freedofs]
            lgmap = psc.LGMap().create(indices=globnums, bsize=eh, comm=comm)
            mat = psc.Mat().create(comm=comm)
            mat.setSizes(size=nglob*eh, bsize=eh)
            mat.setType(psc.Mat.Type.IS)
            mat.setLGMap(lgmap)
            mat.setISLocalMat(apsc_loc)
            mat.assemble()
            mat.convert("mpiaij")
            return mat
    else:
        if freedofs is not None:
            mat = apsc_loc
            mat.assemble()
            mat.convert("seqaij")
            return mat


#PETSc Vector

class VectorMapping:
    def __init__ (self, pardofs, freedofs=None):
        self.pardofs = pardofs
        self.freedofs = freedofs
        comm = pardofs.comm.mpi4py
        globnums, self.nglob = pardofs.EnumerateGlobally(freedofs)
        self.es = self.pardofs.entrysize
        if self.freedofs is not None:
            globnums = np.array(globnums, dtype=psc.IntType)[freedofs]
            self.locfree = np.flatnonzero(freedofs).astype(psc.IntType)            
            self.isetlocfree = psc.IS().createBlock (indices=self.locfree, bsize=self.es, comm=comm)            
        else:
            self.isetlocfree = None
        self.iset = psc.IS().createBlock (indices=globnums, bsize=self.es, comm=comm)

    def CreatePETScVector (self):
        return psc.Vec().createMPI(self.pardofs.ndofglobal, comm=MPI.COMM_WORLD)        

    def CreateNGSolveVector (self):
        return ngs.la.CreateParallelVector(self.pardofs)
    
    def N2P (self, ngs_vector, psc_vector=None):
        if psc_vector is None:
            psc_vector = psc.Vec().createMPI(self.nglob*self.es, bsize=self.es, comm=MPI.COMM_WORLD)
        ngs_vector.Distribute()
        locvec = psc.Vec().createWithArray(ngs_vector.FV().NumPy(), comm=MPI.COMM_SELF)
        if "n2p_scat" not in self.__dict__:
            self.n2p_scat = psc.Scatter().create(locvec, self.isetlocfree, psc_vector, self.iset)
        psc_vector.set(0)
        self.n2p_scat.scatter (locvec, psc_vector, addv=psc.InsertMode.ADD)  # 1 max, 2 sum+keep
        return psc_vector

    def P2N (self, psc_vector, ngs_vector=None):
        if ngs_vector is None:
            ngs_vector = ngs.la.CreateParallelVector(self.pardofs)

        ngs_vector.SetParallelStatus(ngs.la.PARALLEL_STATUS.CUMULATED)
        ngs_vector[:] = 0.0
        locvec = psc.Vec().createWithArray(ngs_vector.FV().NumPy(), comm=MPI.COMM_SELF)
        if "p2n_scat" not in self.__dict__:
            self.p2n_scat = psc.Scatter().create(psc_vector, self.iset, locvec, self.isetlocfree)
        self.p2n_scat.scatter (psc_vector, locvec, addv=psc.InsertMode.INSERT)
        return ngs_vector

    
#PETSc Preconditioner

class PETScPreconditioner(ngs.BaseMatrix):
    def __init__(self,mat,freedofs=None, solverParameters=None):
        ngs.BaseMatrix.__init__(self)
        self.ngsmat = mat
        self.vecmap = VectorMapping (mat.row_pardofs, freedofs)
        self.pscmat = CreatePETScMatrix (mat, freedofs)
        self.precond = psc.PC().create(comm=self.pscmat.getComm())

        self.precond.setOperators(self.pscmat)
        options_object = psc.Options()
        if solverParameters is not None:
            for optName, optValue in solverParameters.items():
                options_object[optName] = optValue

        self.precond.setOptionsPrefix(None) # optionsPrefix)
        self.precond.setFromOptions()
        self.precond.setUp()
        self.pscx, self.pscy = self.pscmat.createVecs()

    def Shape(self):
        return self.ngsmat.shape
    def CreateVector(self,col):
        return self.ngsmat.CreateVector(not col)
    
    def Mult(self,x,y):
        self.vecmap.N2P(x,self.pscx)
        self.precond.apply(self.pscx, self.pscy)
        self.vecmap.P2N(self.pscy, y)
        
    def MultTrans(self,x,y):
        self.vecmap.N2P(x,self.pscx)
        self.precond.applyTranspose(self.pscx, self.pscy)
        self.vecmap.P2N(self.pscy, y)
        

from ngsolve.comp import RegisterPreconditioner

# for backward compatibility
def MakeGAMGPreconditioner(mat, freedofs, flags):
    flags.Set("pc_type", "gamg")
    return PETScPreconditioner(mat, freedofs, flags)

RegisterPreconditioner ("gamg", MakeGAMGPreconditioner)


def MakePETScPreconditioner(mat, freedofs, flags):
    return PETScPreconditioner(mat, freedofs, flags)

RegisterPreconditioner ("petsc", MakePETScPreconditioner, docflags = { \
                "pc_type" : "type of PETSc preconditioner",
                "levels" : "AMG levels" })


#PETSc DMPlex

FACE_SETS_LABEL = "Face Sets"
CELL_SETS_LABEL = "Cell Sets"
EDGE_SETS_LABEL = "Edge Sets"

class DMPlexMapping:
    def __init__(self,mesh=None,name="Default"):
        self.name = name
        print(type(mesh))
        if type(mesh) is ngs.comp.Mesh or type(mesh) is ngm.Mesh:
           self.createPETScDMPlex(mesh) 
        elif type(mesh) is psc.DMPlex:
           self.createNGSMesh(mesh)
        else:
            raise ValueError("Mesh format not recognised.")
    def createNGSMesh(self,plex):
        ngmesh = ngm.Mesh(dim=plex.getCoordinateDim())
        self.ngmesh = ngmesh
        if plex.getDimension() == 2:
            coordinates = plex.getCoordinates().getArray().reshape([-1,2])
            self.ngmesh.AddPoints(coordinates)  
            cstart,cend = plex.getHeightStratum(0)
            vstart, vend = plex.getHeightStratum(2)
            cells = []
            for i in range(cstart,cend):
                sIndex  = plex.getCone(i)
                s1 = plex.getCone(sIndex[0])-vstart
                s2 = plex.getCone(sIndex[1])-vstart
                if np.linalg.det(np.array([coordinates[s1[1]]-coordinates[s1[0]],coordinates[s2[1]]-coordinates[s1[1]]])) < 0.:
                    cells = cells+[[s1[1],s1[0],s2[1]]]
                else:
                    cells = cells+[[s1[0],s1[1],s2[1]]]

        fd = ngmesh.Add(ngm.FaceDescriptor(bc=1))
        self.ngmesh.AddElements(dim=plex.getDimension(), index=1, data=np.asarray(cells,dtype=np.int32), base=0)
        
    def createPETScDMPlex(self,mesh):
        if type(mesh) is ngs.comp.Mesh:
            self.ngmesh = mesh.ngmesh
        else:
            self.ngmesh = mesh
        comm = mesh.comm
        if self.ngmesh.dim == 3:
            if comm.rank == 0:
                V = self.ngmesh.Coordinates()
                T = self.ngmesh.Elements3D().NumPy()["nodes"]
                T = np.array([list(np.trim_zeros(a, 'b')) for a in list(T)])-1
                surfMesh, dim = False, 3
                if len(T) == 0:
                    surfMesh, dim = True, 2
                    T = self.ngmesh.Elements2D().NumPy()["nodes"]
                    T = np.array([list(np.trim_zeros(a, 'b')) for a in list(T)])-1
                plex = psc.DMPlex().createFromCellList(dim, T, V)
                plex.setName(self.name)
                vStart, vEnd = plex.getDepthStratum(0)
                if surfMesh:
                    for e in self.ngmesh.Elements1D():
                        join = plex.getJoin([vStart+v.nr-1 for v in e.vertices])
                        plex.setLabelValue(FACE_SETS_LABEL, join[0], int(e.surfaces[1]))
                else:
                    for e in self.ngmesh.Elements2D():
                        join = plex.getFullJoin([vStart+v.nr-1 for v in e.vertices])
                        plex.setLabelValue(FACE_SETS_LABEL, join[0], int(e.index))
                    for e in self.ngmesh.Elements1D():
                        join = plex.getJoin([vStart+v.nr-1 for v in e.vertices])
                        plex.setLabelValue(EDGE_SETS_LABEL, join[0], int(e.index))
                self.plex = plex
            else:
                plex = psc.DMPlex().createFromCellList(3,
                                                        np.zeros((0, 4), dtype=np.int32),
                                                        np.zeros((0, 3), dtype=np.double))
                self.plex = plex
        elif self.ngmesh.dim == 2:
            if comm.rank == 0:
                V = self.ngmesh.Coordinates()
                T = self.ngmesh.Elements2D().NumPy()["nodes"]
                T = np.array([list(np.trim_zeros(a, 'b')) for a in list(T)])-1
                plex = psc.DMPlex().createFromCellList(2, T, V)
                plex.setName(self.name)
                vStart, vEnd = plex.getDepthStratum(0)   # vertices
                for e in self.ngmesh.Elements1D():
                    join = plex.getJoin([vStart+v.nr-1 for v in e.vertices])
                    plex.setLabelValue(FACE_SETS_LABEL, join[0], int(e.index))
                if not ((1 == self.ngmesh.Elements2D().NumPy()["index"]).all()):
                    for e in self.ngmesh.Elements2D():
                        join = plex.getFullJoin([vStart+v.nr-1 for v in e.vertices])
                        plex.setLabelValue(CELL_SETS_LABEL, join[0], int(e.index))

                self.plex = plex
            else:
                plex = psc.DMPlex().createFromCellList(2,
                                                        np.zeros((0, 3), dtype=np.int32),
                                                        np.zeros((0, 2), dtype=np.double))
                self.plex = plex


#Krylov Solver

counter = 0

class KrylovSolver():
    """
    Inspired by Firedrake solver class.
    """    
    global counter
    counter += 1
    def __init__(self, a, fes, p=None, solver_parameters=None,options_prefix=None):
        self.fes = fes
        a.Assemble()
        Amat = a.mat
        if p is not None:
            p.Assemble()
            Pmat = p.mat
        else:
            Pmat = None
        if not isinstance(Amat, ngs.la.SparseMatrixd):
            raise TypeError("Provided operator is a '%s', not an la.SparseMatrixd" % type(Amat).__name__)
        if Pmat is not None and not isinstance(Pmat, ngs.la.SparseMatrixd):
            raise TypeError("Provided preconditioner is a '%s', not an la.SparseMatrixd" % type(Pmat).__name__)

        self.solver_parameters = solver_parameters
        self.options_prefix = options_prefix
        options_object = psc.Options() 
        for optName, optValue in self.solver_parameters.items():
            options_object[optName] = optValue
	
	#Creating the PETSc Matrix
        Asc = CreatePETScMatrix(Amat, fes.FreeDofs())
        self.A = Asc
        self.comm = MPI.COMM_WORLD
        #Setting up the preconditioner
        if Pmat is not None:
            Psc = CreatePETScMatrix(Pmat, fes.FreeDofs())
            self.P = Psc
        else:
            self.P = Asc
        #Setting options prefix
        self.A.setOptionsPrefix(self.options_prefix)
        self.P.setOptionsPrefix(self.options_prefix)
        self.A.setFromOptions()
        self.P.setFromOptions()

        self.ksp = psc.KSP().create(comm=self.comm)

        # Operator setting must come after null space has been
        # applied
        self.ksp.setOperators(A=self.A, P=self.P)
        # Set from options now (we're not allowed to change parameters
        # anyway).
        self.ksp.setOptionsPrefix(self.options_prefix)
        self.ksp.setFromOptions()
    def solve(self, f):
        f.Assemble()
        u = ngs.GridFunction(self.fes)
        self.vmap = VectorMapping(self.fes.ParallelDofs(), self.fes.FreeDofs())
        upsc, fpsc = self.A.createVecs()
        self.vmap.N2P(f.vec, fpsc)
        self.ksp.solve(fpsc, upsc)
        self.vmap.P2N(upsc, u.vec);
        return u
