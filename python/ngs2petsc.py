import ngsolve as ngs
import petsc4py.PETSc as psc
import slepc4py.SLEPc as spc
from mpi4py import MPI
import numpy as np

def CreatePETScMatrix (ngs_mat, freedofs=None):
    pardofs = ngs_mat.row_pardofs
    comm = pardofs.comm.mpi4py

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
    


    
class PETScPreconditioner(ngs.BaseMatrix):
    def __init__(self,mat,freedofs=None):
        ngs.BaseMatrix.__init__(self)
        self.ngsmat = mat
        self.vecmap = VectorMapping (mat.row_pardofs, freedofs)
        self.mat = CreatePETScMatrix (mat, freedofs)

        self.precond = psc.PC().create()
        self.precond.setType("gamg")
        self.precond.setOperators(self.mat)
        self.precond.setUp()
        self.pscx, self.pscy = self.mat.createVecs()
        
    def Height(self):
        return self.ngsmat.height
    def Width(self):
        return self.ngsmat.width
    def CreateRowVector(self):
        return self.ngsmat.CreateColVector()
    def CreateColVector(self):
        return self.ngsmat.CreateRowVector()
    
    def Mult(self, x,y):
        self.vecmap.N2P(x,self.pscx)
        self.precond.apply(self.pscx, self.pscy)
        self.vecmap.P2N(self.pscy, y)
        

def MakePreconditioner(mat, freedofs):
    return PETScPreconditioner(mat, freedofs)

class SLEPcEigenProblem():
    def __init__(self,PType,SType):
        self.E = spc.EPS().create()
        if PType == "GHEP":
            self.E.setProblemType(spc.EPS.ProblemType.GHEP) 
        elif PType == "HEP":
            self.E.setProblemType(spc.EPS.ProblemType.HEP)
        elif PType == "NHEP":
            self.E.setProblemType(spc.EPS.ProblemType.NHEP)
        elif PType == "GNHEP":
            self.E.setProblemType(spc.EPS.ProblemType.GNHEP)
        self.E.setType(SType)
    def setOperators(self,mats,Dofs):
        if len(mats)==1:
            self.A = CreatePETScMatrix(mats[0],Dofs);
            self.E.setOperators(self.A)
        elif len(mats)==2:
            self.A = CreatePETScMatrix(mats[0],Dofs);
            self.M = CreatePETScMatrix(mats[1],Dofs);
            self.E.setOperators(self.A,self.M)
    def SpectralTransformation(self,opt):
        self.ST = self.E.getST();
        self.ST.setType(opt);
        self.KSP = self.ST.getKSP();
    def setWhich(self,n):
        self.E.setDimensions(n,psc.DECIDE)
    def getPair(self,fes,s):
        xr, xi = self.A.createVecs()
        lam = self.E.getEigenpair(s, xr, xi)
        vecmap = VectorMapping(fes.ParallelDofs(), fes.FreeDofs())
        vr =  ngs.GridFunction(fes)
        vecmap.P2N(xr, vr.vec)
        vi =  ngs.GridFunction(fes)
        vecmap.P2N(xi, vi.vec)
        return lam,vr,vi
    def Solve(self):
        self.E.setST(self.ST);
        self.E.solve();
        self.Iterations = self.E.getIterationNumber()
        self.Converged = self.E.getConverged()
        sol_type = self.E.getType()
        nev, ncv, mpd = self.E.getDimensions()
        self.tolerance, Mit = self.E.getTolerances()
from ngsolve.comp import RegisterPreconditioner
RegisterPreconditioner ("gamg", MakePreconditioner)


    
