import ngsolve as ngs
import petsc4py.PETSc as psc
from mpi4py import MPI
import numpy as np

def CreatePETScMatrix (ngs_mat):
    pardofs = ngs_mat.row_pardofs
    # comm = MPI.COMM_WORLD  
    comm = pardofs.comm.mpi4py
    globnums, nglob = pardofs.EnumerateGlobally()
    iset = psc.IS().createGeneral (indices=globnums, comm=comm)
    lgmap = psc.LGMap().createIS(iset)

    
    locmat = ngs_mat.local_mat
    val,col,ind = locmat.CSR()
    ind = np.array(ind, dtype='int32')
    apsc_loc = psc.Mat().createAIJ(size=(locmat.height, locmat.width), csr=(ind,col,val), comm=MPI.COMM_SELF)
    
    mat = psc.Mat().createPython(size=nglob, comm=comm)
    mat.setType(psc.Mat.Type.IS)
    mat.setLGMap(lgmap)
    mat.setISLocalMat(apsc_loc)
    mat.assemble()
    mat.convert("mpiaij")
    return mat



class VectorMapping:
    def __init__ (self, pardofs):
        self.pardofs = pardofs
        comm = pardofs.comm.mpi4py        
        globnums, nglob = pardofs.EnumerateGlobally()
        self.iset = psc.IS().createGeneral (indices=globnums, comm=comm)

    def CreatePETScVector (self):
        return psc.Vec().createMPI(self.pardofs.ndofglobal, comm=MPI.COMM_WORLD)        

    def CreateNGSolveVector (self):
        return ngs.la.CreateParallelVector(self.pardofs)
    
    def N2P (self, ngs_vector, psc_vector=None):
        if psc_vector is None:
            psc_vector = psc.Vec().createMPI(self.pardofs.ndofglobal, comm=MPI.COMM_WORLD)
        ngs_vector.Cumulate()
        psc_loc = psc_vector.getSubVector(self.iset)
        psc_loc.getArray()[:] = ngs_vector.FV()
        psc_vector.restoreSubVector(self.iset, psc_loc)
        return psc_vector

    def P2N (self, psc_vector, ngs_vector=None):
        if ngs_vector is None:
            ngs_vector = ngs.la.CreateParallelVector(self.pardofs)
        psc_loc = psc_vector.getSubVector(self.iset)
        ngs_vector.FV()[:] = ngs.Vector(psc_loc.getArray())
        return ngs_vector
    

    
