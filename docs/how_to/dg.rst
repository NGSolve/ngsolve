Discontinuous Galerkin methods
==============================

Discontinuous Galerkin (DG) methods have certain advantages: One can apply upwinding for convection dominated problems, and explicit time-stepping methods are cheap due to block-diagonal or even diagonal mass matrices.

The bilinear-form of a DG method involves integrals over functions defined on neighbouring elements. In NGS-Py we can access the neighbouring element via the *.Other()* method. The following example gives the boundary integral of an upwind scheme for the convection equation. The *CoefficientFunction b* is a given vector-field, the wind. With specialcf.normal(2) we get the outer element normal vector. We have to provide the space dimension of the mesh. Depending on the sign of *<b,n>* we choose the trial function from the current element, or from the neighbour using the *IfPos* function. If the edge is on the domain-boundary, a given boundary value may be specified by the *bnd* argument:

.. code-block:: python

                b = CoefficientFunction( (y-0.5,0.5-x) ) 
                bn = b*specialcf.normal(2)
                a += SymbolicBFI (bn*IfPos(bn, u, u.Other(bnd=ubnd)) * v, element_boundary=True)

One can also use the test-function from the neighbour element. The coefficient functions are always evaluated on the current element.

When one works with assembled matrices, there is a drawback of DG methdos: The matrix stencil becomes larger. We have to tell NGSolve to reserver more entries in the matrix using

.. code-block:: python

                FESpace( ... , flags = { "dgjumps" : True })

If we don't assemble the matrix, but work with operator application on the fly, we don't have to specify it.


The above expression leads to a loop over elements, and the boundary integrals are evaluated for the whole element-boundary, which consists of internal and boundary edges (or faces). But sometimes we need different terms for internal and boundary facets. Here we can use the *skeleton* flag. This leads to separate loops for internal and boundary facets. The following definition is mathematically equivalent to the method above. The VOL and BND specifier tell whether we want to loop over internal or boundary edges. Since here every edge is processed only once, the skeleton formulation is slightly more efficient:

.. code-block:: python

                a += SymbolicBFI ( bn*IfPos(bn, u, u.Other()) * (v-v.Other()), VOL, skeleton=True)
                a += SymbolicBFI ( bn*IfPos(bn, u, ubnd) * v, BND, skeleton=True)


To solve with the block-diagonal (or even diagonal) mass matrix of an *L2* -finite element space, we can use the *SolveM* method of the FESpace. The *rho* argument allows to specify a density coefficientfunction for the mass-matrix. The operation is performed inplace for the given vector.
                
 .. code-block:: python
                 
                density = CoefficientFunction(1)
                fes.SolveM (rho=density, vec=u)


Several examples of DG methods are given in the DG directory of the py_tutorials.
