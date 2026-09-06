from ngsolve import *
from ngsolve.la import *
from pyngcore import MPI_Comm
import mpi4py.MPI as mpi

# formats of parallel vectors and of agnostic operators next to a ParallelMatrix
def test_parallel_formats():
    comm = MPI_Comm(mpi.COMM_WORLD)
    mesh = Mesh('square.vol.gz', comm)
    fes = H1(mesh, order=2, dirichlet=".*")
    u, v = fes.TnT()
    a = BilinearForm(u*v*dx + grad(u)*grad(v)*dx).Assemble()
    f = LinearForm(v*dx).Assemble()

    fmt = a.mat.RowFormat()
    assert fmt.is_parallel == (comm.size > 1)
    assert fmt.scalar == "double"
    assert f.vec.GetFormat().is_parallel == (comm.size > 1)

    I0 = IdentityMatrix()
    for op in [I0 @ a.mat, a.mat @ I0, a.mat + I0, 2.0*a.mat, a.mat.T, I0 @ (a.mat @ I0)]:
        vec = op.CreateColVector()
        assert vec.GetFormat().is_parallel == (comm.size > 1), str(op.ColFormat())
        assert vec.size == f.vec.size

    # identity applied to a parallel vector keeps distribution and values
    y = (I0 * f.vec).Evaluate()
    assert y.GetFormat().is_parallel == (comm.size > 1)
    assert abs(Norm(y - f.vec)) < 1e-14 * max(1.0, Norm(f.vec))

    # a product with an agnostic factor multiplies like the matrix
    y1 = (a.mat * f.vec).Evaluate()
    y2 = ((I0 @ a.mat) * f.vec).Evaluate()
    assert abs(Norm(y1 - y2)) < 1e-12 * max(1.0, Norm(y1))
