Working with CoefficientFunctions
=================================

A CoefficientFunction is a function which can be evaluated on a mesh, and may be used to provide the coefficient or a right-hand-side to the variational formulation. Because typical finite element procedures iterate over elements, and map integration points from a reference element to a physical element, the evaluation of a CoefficientFunction requires a mapped integration point. The mapped integration point contains the coordinate, as well as the Jacobian, but also element number and region index. This allows the efficient evaluation of a wide class of functions needed for finite element computations.

.. autoclass:: ngsolve.CoefficientFunction

Basic Coefficient Functions
----------------------------

Python objects implicitly convertible to *float* or *Complex* are implicitly converted to *CoefficientFunction* if needed. The Cartesian coordinates *x*, *y*, and *z* are pre-defined coordinate coefficient functions:

.. code-block:: python

   f = LinearForm(V)
   f += SymbolicLFI( x * V.TestFunction())

A *CoefficientFunction* initialized with a *list* stores one value per region, depending where the *CoefficientFunction* is used this is a domain or boundary region. Here one usually uses generator expressions in combination with ``mesh.GetMaterials()`` and ``mesh.GetBoundaries()``:

.. code-block:: python

   alpha_d = {"air" : 1, "box" : 100}
   alpha = CoefficientFunction([alpha_d[mat] for mat in mesh.GetMaterials()])

A *CoefficientFunction* initialized with a *tuple* gives a vector coefficient function. Components of the coefficient function can be accessed by the bracket operator:

.. code-block:: python

   vec_cf = CoefficientFunction( (fx,fy) )  # Creates a vector CF
   fx_again = vec_cf[0]                     # Creates a scalar CF from the first component

A matrix valued CF can be created with the additional *dims* argument:

.. code-block:: python

   mat_cf = CoefficientFunction((f11,f12,f21,f22),dims=(2,2))

ProxyFunctions
---------------

The finite element spaces provide proxy placeholder functions for :any:`Symbolic Integrators<symbolic-integrators>` 
when assembling the system matrices the proxy functions are replaced with the FEM basis functions.

.. automethod:: ngsolve.FESpace.TrialFunction

.. automethod:: ngsolve.FESpace.TestFunction

Operations on CoefficientFunctions
-----------------------------------

CoefficientFunctions can be combined with algebraic operations (+,-,*,/,**).
Math functions are available as well for CoefficientFunctions: *sin*, *cos*,
*tan*, *exp*, *log*, *atan*, *sqrt*, *Conj*.

Special CoefficientFunctions
-------------------------------

NGSolve provides special coefficient functions needed for finite element computations.
Special coefficient functions are mesh dependent.

You can get the normal vector on an interface with

.. autoattribute:: ngsolve.specialcf.normal

this can be used in :any:`discontinuous Galerkin methods<discontinuous-galerkin>`.

The tangential vector can be obtained in the same way:

.. autoattribute:: ngsolve.specialcf.tangential

The local mesh size can be obtained by a coefficient function as well:

.. autoattribute:: ngsolve.specialcf.mesh_size

This has application i.e. in :download:`hybrid DG methods</../py_tutorials/hybrid_dg.py>`.


Additional CoefficientFunctions
--------------------------------

IfPos
^^^^^^

The function *IfPos* provides CoefficientFunctions depending on some condition:

.. autofunction:: ngsolve.IfPos

Parameter CoefficientFunction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want a CoefficientFunction with a variable parameter, instead of using a 
new *BilinearForm* every time you can use a *Parameter* CF. The parameter can be 
modified with the *Set* method and when assembling the BilinearForm again, the 
updated parameter is used.

.. autoclass:: ngsolve.Parameter
   :members:

BSpline CoefficientFunction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can create a *BSpline* as a CoefficientFunction as well. BSplines are
differentiable and integrable:

.. autoclass:: ngsolve.BSpline
   :members:

Compiling a CoefficientFunctions
----------------------------------

.. automethod:: ngsolve.CoefficientFunction.Compile

Evaluating CoefficientFunctions
--------------------------------

Sometimes it can be useful to evaluate a *CoefficientFunction*. Since some CoefficientFunctions are only defined on the mesh (like a *GridFunction*) this can only be done with information from the mesh. For this we request a mapped integration point from the mesh and plug it into the *CoefficientFunction*:


.. code-block:: python

  >>> mip = mesh(0.2,0.4)
  >>> cf = x*x*y
  >>> cf(mip)
  0.016000000000000004
