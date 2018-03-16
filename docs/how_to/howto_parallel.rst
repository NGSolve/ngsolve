Parallel computing with NGS-Py
==============================

There are several options to run NGS-Py in parallel, either in a shared-memory, or distributed memory paradigm.

Shared memory parallelisation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

NGSolve shared memory parallelisation is based on a the task-stealing
paradigm. On entering a parallel execution block, worker threads are
created. The master thread executes the algorithm, and whenever a
parallelized function is executed, it creates tasks. The waiting
workers pick up and process these tasks. Since the threads stay alive
for a longer time, these paradigm allows to parallelize also very
small functions, practically down to the range of 10 micro seconds.

The task parallelization is also available in NGS-Py. By the *with
Taskmanager* statement one creates the threads to be used in the
following code-block. At the end of the block, the threads are stopped.


.. code-block:: python

                with Taskmanager():
                    a = BilinearForm(fespace)
                    a += SymbolicBFI(u*v)
                    a.Assemble()

                    
Here, the assembling operates in parallel. The finite element space
provides a coloring such that elements of the same color can be
processed simultaneously. Also helper functions such as sparse matrix
graph creation uses parallel loops.

The default number of threads is the (logical) number of cores.
It can be overwritten by the environment variable NGS_NUM_THREADS,
or by calling the python function

.. code-block:: python

                SetNumThreads ( num_threads )


Another typical example for parallel execution are equation
solvers. Here is a piece of code of the conjugate gradient solver from
NGS-Py:

.. code-block:: python

                with Taskmanager():

                  ...
                  for it in range(maxsteps):
                      w.data = mat * s
                      wd = wdn
                      as_s = InnerProduct (s, w)
                      alpha = wd / as_s
                      u.data += alpha * s
                      d.data += (-alpha) * w


The master thread executes the algorithm. In matrix - vector product
function calls, and also in vector updates and innner products tasks
are created and picked up by workers. 
        


Distributed memory
^^^^^^^^^^^^^^^^^^


The distributed memory paradigm requires to build Netgen as well as NGSolve with MPI - support, which must be enabled during the cmake configuration step. 

Many ngsolve features can be used in the MPI-parallel version, some
features are work in progress, some others may take for longer. The
following list shows what is available:

.. csv-table:: Distributed parallel features
               :header: feature , y/n , comment
               :widths:  20,10,50

               TaskManager, yes, allows for hybrid parallelization
               mesh generation, no,
               mesh distribution, yes, using metis
               uniform refinement, yes,  refine the netgen-mesh
               adaptive refinement, no,
               finite element spaces, yes, all spaces should work
               "Forms, Gridfunction", yes
               periodic boundaries, no
               discontinuous Galerkin, yes, apply operator
               diagonal preconditioner, yes
               multigrid preconditioner, no
               BDDC preconditioner, yes
               direct solvers, yes, "MUMPs, Masterinverse"
               pickling, in progress
               

