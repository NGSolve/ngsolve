Define and update preconditioners
=================================

A preconditioner is defined for a bilinear-form, and aims at providing a cheap, approximative inverse of the matrix. The matrix is restricted to the non-dirichlet (free) degrees of freedom, provided by the underlying FESpace.

The canonical way is to define the preconditioner after the bilinear-form, but before calling Assemble:

.. code-block:: python

   a = BilinearForm(fes)
   a += SymbolicBFI(grad(u)*grad(v))
   c = Preconditioner(a, "local")
   a.Assemble()
   

The preconditioner registers itself with the bilinear-form. Whenever the form is updated, the preconditioner is updated as well.

You can define the preconditioner after assembling, but then you have to call manually c.Update()

The ratio if this ordering is that some preconditioners (e.g. bddc, amg, ...) require access to the element-matrices, which are only available during assembling.

The preconditioners included in NGSolve are the following. Additional user-defined preconditioners can be implemented in plug-ins. An example is given in MyLittleNGSolve



+---------------+---------------------------------------------+
| name          | preconditioner                              |
+===============+=============================================+
|local          |Jacobi / block-Jacobi                        |
+---------------+---------------------------------------------+
|direct         |a sparse direct factorization                |
+---------------+---------------------------------------------+
|multigrid      |h-version and high-order/low-order multigrid |
+---------------+---------------------------------------------+
|bddc           |p-version domain decomposition               |
+---------------+---------------------------------------------+


Block smoothers
---------------

A block-Jacobi or block-Gauss-Seidel preconditioner inverts the sub-matrices belonging to small, overlapping sets of degrees of freedom. The sets are chosen by the flag ``blocktype``:

.. code-block:: python

   c = Preconditioner(a, "local", blocktype="vertexpatch")
   c = Preconditioner(a, "local", blocktype=["vertexedge", "face", "inner"], GS=True)
   c = Preconditioner(a, "multigrid", blocktype="vertexpatch")

The first line gives an additive, the second one a symmetric multiplicative (Gauss-Seidel) preconditioner. The multigrid preconditioner uses the blocks for its smoother on the fine levels.

The named block types are defined by the mesh topology and work for every finite element space:

+---------------+------------------------------------------------------------------------+
| blocktype     | one block per ..., containing                                          |
+===============+========================================================================+
|vertex         |vertex: its dofs                                                        |
+---------------+------------------------------------------------------------------------+
|edge           |edge: its dofs                                                          |
+---------------+------------------------------------------------------------------------+
|face           |face: its dofs                                                          |
+---------------+------------------------------------------------------------------------+
|cell           |cell (3D): its dofs                                                     |
+---------------+------------------------------------------------------------------------+
|facet          |facet, codimension 1 (edge in 2D, face in 3D): its dofs                 |
+---------------+------------------------------------------------------------------------+
|inner          |element, codimension 0: its inner dofs                                  |
+---------------+------------------------------------------------------------------------+
|element        |element: all its dofs, including the ones on its boundary               |
+---------------+------------------------------------------------------------------------+
|vertexedge     |vertex: its dofs and the dofs of the edges containing it                |
+---------------+------------------------------------------------------------------------+
|vertexpatch    |vertex: its dofs and the dofs of all edges, faces, cells containing it  |
+---------------+------------------------------------------------------------------------+
|edgepatch      |edge: its dofs and the dofs of all faces, cells containing it           |
+---------------+------------------------------------------------------------------------+
|facepatch      |face: its dofs and the dofs of the cells containing it                  |
+---------------+------------------------------------------------------------------------+

A list of names combines several types. The blocks follow each other in the order of the list, which is the order of a Gauss-Seidel sweep. The patches are open: the block of a vertex patch does not contain the dofs on the boundary of the patch.

Two options can be appended to a name, separated by a colon: ``"vertexpatch:periodic"`` merges the blocks of periodic vertices, and ``"edge:hdivlo"`` restricts the blocks to a set of dofs known to the space (here the lowest order dofs of an HDiv space).

Dirichlet dofs never enter the blocks. If the bilinear-form uses static condensation, the eliminated dofs are left out as well.

Without a named type, ``block=True`` selects blocks chosen by the space, for example the dofs of a vertex and its edges, the faces, and the element interiors for H1. Some spaces provide further types selected by an integer ``blocktype``.

The blocks are also available directly, for programming your own preconditioners:

.. code-block:: python

   blocks = fes.CreateSmoothingBlocks("vertexpatch")
   smoother = a.mat.CreateBlockSmoother(blocks)

``smoother`` can be applied as a block-Jacobi preconditioner, and provides the methods ``Smooth`` and ``SmoothBack`` for block-Gauss-Seidel sweeps.
