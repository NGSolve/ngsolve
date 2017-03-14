============================
 Interfacing to numpy/scipy
============================

In some occasions or for some users it might be interesting to access NGSolve data from python in a fashion which is compatible with numpy and/or scipy. 
We give a few examples of possible use cases.

Working with small vectors and dense matrices:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
see [Vectors and matrices](ngspy-howto-linalg)

Working with large vectors
^^^^^^^^^^^^^^^^^^^^^^^^^^

You can get a "view" on an NGSolve-BaseVector using `.FV()` (which
will give you a FlatVector) combined with `.NumPy()` which will give a
numpy array which operates on the NGSolve-Vector-Data. For example the
following works, assuming b to be an NGSolve-Vector:

.. code-block:: python

                b.FV().NumPy()[:] = abs(b.FV().NumPy()) - 1.0

which will give you the component-wise operation (absolute value minus one) applied on the vector b. During this operation data does not need to be copied.

Working with sparse matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can access the sparse matrix information of a BaseMatrix using

.. code-block:: python

                rows,cols,vals = a.mat.COO()

Note that a bilinear form with flag `symmetric==True` will only give you one half of the matrix.
These information can be put into a standard scipy-matrix, e.g. with 

.. code-block:: python

                import scipy.sparse as sp
                A = sp.csr_matrix((vals,(rows,cols)))

You can use this, for instance, to examine the sparsity pattern of your matrix:

.. code-block:: python

                import matplotlib.pylab as plt
                plt.spy(A)
                plt.show()

or to compute the condition number (note that we export to a dense matrix here):

.. code-block:: python

                import numpy as np
                np.linalg.cond(A.todense())


Using iterative solvers from scipy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To use iterative solvers from scipy we have to wrap a `LinearOperator` around the NGSolve-matrix. The crucial component is the application of matrix vector product. Here is a very simple example where the scipy-cg-solver is used to solve the linear system (no preconditioner, no Dirichlet-dofs):

.. code-block:: python

                import scipy
                import scipy.sparse.linalg

                tmp1 = f.vec.CreateVector()
                tmp2 = f.vec.CreateVector()
                def matvec(v):
                    tmp1.FV().NumPy()[:] = v
                    tmp2.data = a.mat * tmp1
                    return tmp2.FV().NumPy()

                A = scipy.sparse.linalg.LinearOperator( (a.mat.height,a.mat.width), matvec)

                u.vec.FV().NumPy()[:], succ = scipy.sparse.linalg.cg(A, f.vec.FV().NumPy())

                
---------------------------------------------------------------------------

You can also use a sparse matrix format from python to run to previous example, see above. However, for preconditioning actions a sparse matrix is not necessarily set up such that the `LinearOperator` is often more useful.
