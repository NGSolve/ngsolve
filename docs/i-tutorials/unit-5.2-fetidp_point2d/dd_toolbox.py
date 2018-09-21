from ngsolve import *
from ngsolve.la import ParallelMatrix, ParallelDofs, SparseMatrixd
from ngsolve.la import DISTRIBUTED, CUMULATED, CreateParallelVector
from functools import reduce

def BlockOp2SPM(A, h, w, take_dofs):
    cv = A.CreateRowVector()
    cv2 = A.CreateRowVector()
    nblocks = cv.nblocks
    import numpy as np
    firstis = np.cumsum([0,*[len(cv[l]) for l in range(nblocks)]])
    rcv_list = []
    for C in range(nblocks):
        N = len(cv[C])
        cv[:] = 0
        evec = cv[C]
        for c in range(N):
            if take_dofs[firstis[C]+c]:
                cv[:] = 0
                cv[C][c] = 1
                cv2.data = A * cv
                rcv_list += [ (firstis[R]+j, firstis[C]+c, cv2[R][j]) \
                              for R in range(nblocks) \
                              for j in range(len(cv2[R])) \
                              if take_dofs[firstis[R]+j] ]
    [rows, cols, vals] = zip(*rcv_list) if len(rcv_list) else [[],[],[]]
    return SparseMatrixd.CreateFromCOO(list(rows), list(cols), list(vals), h, w)

def Op2SPM(A, take_rows=None, take_cols=None):
    evec = A.CreateRowVector()
    W = len(evec)
    w = W if take_rows==None else take_rows.NumSet()
    tv = A.CreateColVector()
    H = len(tv)
    h = H if take_cols==None else take_cols.NumSet()
    vals = [0 for k in range(w*h)]
    rows = [0 for k in range(w*h)]
    cols = [0 for k in range(w*h)]
    c = 0
    rdofs = [k for k in range(H)] if take_rows==None \
            else [v[0] for v in enumerate(take_rows) if v[1]]
    cdofs = [k for k in range(W)] if take_cols==None \
            else [v[0] for v in enumerate(take_cols) if v[1] ]
    for k in range(h):
        for j in range(w):
            rows[c] = rdofs[k]
            cols[c] = cdofs[j]
            c = c+1
    c = 0
    for k in range(h):
        evec[:] = 0
        evec[rdofs[k]] = 1.0
        tv.data = A * evec
        for j in range(w):
            vals[c] = tv[cdofs[j]]
            c = c+1
    return SparseMatrixd.CreateFromCOO(rows, cols, vals, H, W)

class ZeroMat(BaseMatrix):
    def __init__(self, h):
        super(ZeroMat, self).__init__()
        self.h = h
        self.w = h
    def Mult(self, x, y):
        y[:] = 0
    def MultTransAdd(self, s, x, y):
        y[:] = 0
    def CreateRowVector(self):
        return CreateVVector(self.w)
    def CreateColVector(self):
        return CreateVVector(self.h)
    def VHeight(self):
        return self.h
    def VWidth(self):
        return self.w

class SPInv(BaseMatrix):
    def __init__(self, A, Ainv, B, invtype_loc='sparsecholesky'):
        super(SPInv, self).__init__()
        self.A = A
        self.Ainv = Ainv
        self.B = B
        S = B @ Ainv @ B.T
        self.S_sparse = Op2SPM(S)
        self.Sinv = self.S_sparse.Inverse(inverse=invtype_loc)
        Sinv_sparse = Op2SPM(self.Sinv)
        self.tlam = self.B.CreateColVector()
        self.tu = self.B.CreateRowVector()        
    def CreateVector(self):
        return self.CreateRowVector()
    def CreateRowVector(self):
        return BlockVector([self.A.CreateRowVector(), self.B.CreateColVector()])
    def CreateColVector(self):
        return BlockVector([self.A.CreateColVector(), self.B.CreateColVector()])
    def VHeight(self):
        return self.A.height+self.B.height
    def VWidth(self):
        return self.A.width+self.B.height
    def Mult(self, x, y):
        self.tu.data = self.Ainv * x[0]
        self.tlam.data = self.B * self.tu
        self.tlam.data -= x[1]
        y[1].data = self.Sinv * self.tlam
        self.tu.data = x[0]
        self.tu.data -= self.B.T * y[1]
        y[0].data = self.Ainv * self.tu
    def MultTransAdd(self, s, x, y):
        self.MultAdd(s, x, y)

class DPSpace_Inverse(BaseMatrix):
    def CreateRowVector(self):
        return self.A_dp.CreateRowVector()
    def CreateColVector(self):
        return self.A_dp.CreateColVector()
    def __init__ (self, mat, freedofs, c_points, c_mat, c_pardofs, \
                  invtype_loc='sparsecholesky', invtype_glob='masterinverse'):
        super(DPSpace_Inverse, self).__init__()
        self.A = mat.local_mat
        pardofs = mat.row_pardofs
        self.c_mat = c_mat
        comm = MPI_Init()
        ess_pardofs = pardofs.SubSet(c_points)
        self.A_dp = ParallelMatrix(self.A, ess_pardofs)
        self.nu = pardofs.ndoflocal
        self.nup = self.c_mat.height
        self.zm = ZeroMat(self.c_mat.height)
        self.Mgg = BlockMatrix([[self.A, None], \
                           [None, self.zm]])
        self.Mgl = BlockMatrix([[self.A, None], \
                           [None, self.zm-IdentityMatrix()]])
        self.Mll = BlockMatrix([[self.A, c_mat.T], \
                           [c_mat, None]])
        self.Add_inv = self.A.Inverse(BitArray(~c_points & freedofs), inverse=invtype_loc)
        self.Mll_inv = SPInv(self.A, self.Add_inv, self.c_mat, invtype_loc)
        self.Sgg = self.Mgg - self.Mgl @ self.Mll_inv @ self.Mgl.T
        s_distprocs = [list(pardofs.Dof2Proc(k)) for k in range(pardofs.ndoflocal)] \
                      + \
                      [list(c_pardofs.Dof2Proc(k)) for k in range(c_pardofs.ndoflocal)]
        s_freedofs = BitArray([b for b in c_points]+[True for k in range(self.nup)])
        Spardofs = ParallelDofs(s_distprocs, comm)
        Sgg_sparse = BlockOp2SPM(self.Sgg, Spardofs.ndoflocal, Spardofs.ndoflocal, s_freedofs)
        self.S = ParallelMatrix(Sgg_sparse, Spardofs)
        self.Sinv = self.S.Inverse(s_freedofs, inverse=invtype_glob)
        self.RhsOp = IdentityMatrix() - (self.Mgl.T @ self.Mll_inv)
        self.bx = BlockVector([self.A.CreateRowVector(), \
                               self.c_mat.CreateColVector()])
        self.brhs = BlockVector([self.A.CreateRowVector(), \
                                 self.c_mat.CreateColVector()])
        self.rhs = self.S.CreateRowVector()
        self.sol = self.S.CreateRowVector()
    def Mult(self, x, y):
        x.Distribute()
        self.bx[0].data = x
        self.bx[1][:] = 0
        self.brhs.data = self.RhsOp * self.bx
        self.rhs.local_vec[0:self.nu] = self.brhs[0]
        self.rhs.local_vec[self.nu:self.nu+self.nup] = self.brhs[1]
        self.rhs.SetParallelStatus(DISTRIBUTED)
        self.sol.data = self.Sinv * self.rhs
        self.sol.Cumulate()
        self.brhs[0].data = self.sol.local_vec[0:self.nu]
        self.brhs[1].data = self.sol.local_vec[self.nu:self.nu+self.nup]
        self.bx.data -= self.Mgl * self.brhs
        self.brhs.data = self.Mll_inv * self.bx
        y.local_vec.data = self.sol.local_vec[0:self.nu] + self.brhs[0]
        y.SetParallelStatus(CUMULATED)

class ScaledMat(BaseMatrix):
    def __init__ (self, A, vals):
        super(ScaledMat, self).__init__()
        self.A = A
        self.scale = [k for k in enumerate(vals) if k[1]!=1.0]
    def CreateRowVector(self):
        return self.A.CreateRowVector()
    def CreateColVector(self):
        return self.A.CreateColVector()
    def Mult(self, x, y):
        y.data = self.A * x
        for k in self.scale:
            y[k[0]] = y[k[0]] * k[1]

class LocGlobInverse(BaseMatrix):
    def __init__ (self, Aglob, freedofs, invtype_loc='sparsecholesky',\
                  invtype_glob='masterinverse'):
        super(LocGlobInverse, self).__init__()
        self.Aglob = Aglob  
        self.A = Aglob.local_mat  
        pardofs = Aglob.col_pardofs 
        local_dofs = BitArray([len(pardofs.Dof2Proc(k))==0 for k in range(pardofs.ndoflocal)]) & freedofs
        global_dofs = BitArray(~local_dofs & freedofs)
        self.All_inv = self.A.Inverse(local_dofs, inverse=invtype_loc) 
        sp_loc_asmult = self.A @ (IdentityMatrix() - self.All_inv @ self.A)
        sp_loc = Op2SPM(sp_loc_asmult, global_dofs, global_dofs)
        sp = ParallelMatrix(sp_loc, pardofs)
        self.Sp_inv = sp.Inverse(global_dofs, inverse=invtype_glob)
        self.tv = self.Aglob.CreateRowVector()
        self.btilde = self.Aglob.CreateRowVector()
    def CreateRowVector(self):
        return self.Aglob.CreateRowVector()
    def CreateColVector(self):
        return self.CreateRowVector()
    def Mult(self, b, u):
        self.tv.data = self.All_inv * b
        b.Distribute()
        self.btilde.data = b - self.A * self.tv
        u.data = self.Sp_inv * self.btilde
        self.tv.data =  b - self.A * u
        u.data += self.All_inv * self.tv

def FindFEV(dim, NV, pardofs, freedofs):
    if freedofs==None:
        freedofs = BitArray([True for k in range(pardofs.ndoflocal)])
    if dim==3:
        faces = [set(d for d in pardofs.Proc2Dof(p) if freedofs[d] and d<NV ) for p in pardofs.ExchangeProcs()]
        edges = sorted([tuple(sorted(e)) for e in set(tuple(f1.intersection(f2)) for f1 in faces for f2 in faces if f1 is not f2) if len(e)>1])
    else:
        edges = [sorted(tuple(d for d in pardofs.Proc2Dof(p) if freedofs[d] and d<NV )) for p in pardofs.ExchangeProcs()]
        faces = edges

    vertices = sorted(set([ v for e1 in edges for e2 in edges if e1 is not e2 for v in set(e1).intersection(set(e2)) ]))

    vec = CreateParallelVector(pardofs)
    vec.local_vec[:] = 0.0
    for v in vertices:
        vec.local_vec[v] = 1
    vec.SetParallelStatus(DISTRIBUTED)
    vec.Cumulate()
    vertices = [ v for v in range(NV) if vec.local_vec[v]!=0 ]

    if dim==3:
        all_e = reduce(lambda x,y: set(x).union(set(y)), faces) if len(faces) else {}
        faces2 = [[v for v in f if not v in all_e] for f in faces]
        faces = [f for f in faces2 if len(f)]
    
    edges2 = [[v for v in e if not v in vertices] for e in edges] 
    edges = [e for e in edges2 if len(e)]
    return [faces, edges, vertices]



# Okay, this is a bit annoying:
# This class is a wrapper around a block-operator,
# that simply linearizes the blockvectors and itself
# works on standard ParallelVectors.
# Hard coded for 2x2 blocks, and does currently not
# work on multidim vectors.
#
# This is needed for the (C++-side) GMRes
#
# TODO: implement something like this on C++ side ASAP
#       that way we can also get rid of unnecessary vector
#       copying
#
#    (at the very least put a ParallelFlatVector(loc_vec, pardofs)
#     into the python interface)
#
#        -- rant over ---
#
class LinMat(BaseMatrix):
    def __init__ (self, A, pd_array):
        super(LinMat, self).__init__()
        self.A = A
        comm = MPI_Init()
        dist_procs = reduce(lambda x,y:x+y, [[list(p.Dof2Proc(k)) for k in range(p.ndoflocal)] for p in pd_array])
        self.long_pardofs = ParallelDofs(dist_procs, comm)
        self.n1 = pd_array[0].ndoflocal
        self.n2 = self.n1+pd_array[1].ndoflocal

        self.x = self.A.CreateRowVector()
        self.y = self.A.CreateRowVector()

    def IsComplex(self):
        return False
    def Height(self):
        return self.n2
    def Width(self):
        return self.n2
    def CreateRowVector(self):
        return CreateParallelVector(self.long_pardofs)
    def CreateColVector(self):
        return CreateParallelVector(self.long_pardofs)
    def Mult(self, x, y):
        self.x[0].local_vec.data = x.local_vec[0:self.n1]
        self.x[1].local_vec.data = x.local_vec[self.n1:self.n2]
        self.x[0].SetParallelStatus(x.GetParallelStatus())
        self.x[1].SetParallelStatus(x.GetParallelStatus())
        self.y.data = self.A * self.x
        self.y[0].Cumulate()
        self.y[1].Cumulate()
        y.local_vec[0:self.n1] = self.y[0].local_vec
        y.local_vec[self.n1:self.n2] = self.y[1].local_vec
        y.SetParallelStatus(CUMULATED)
