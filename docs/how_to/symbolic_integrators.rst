.. _symbolic-integrators:

Symbolic Integrators
=====================

Symbolic integrators can be used to define arbitrary (multi-)physical problems. The variational formulation of the (non-)linear problem can be implemented in a very natural way. Examples of the usage of these integrators can i.e. be found in :any:`Navier Stokes Equation <how-to-navierstokes>`, :any:`Magnetic Field <whetting-cmagnet>` or :any:`Nonlinear elasticity <how-to-elasticity>`.

The finite element space provides placeholder coefficient functions by the methods ``TestFunction`` and ``TrialFunction``. They have implemented canonical derivatives and traces of the finite element space and insert the basis functions of the FESpace during assembling of the system matrices.

Linear Problems
----------------

For linear problems we use the function ``Assemble`` of the ``BilinearForm`` to assemble the matrix and vector. For example for the :math:`H_\text{curl}` linear problem 

.. math::

   \int_\Omega \mu^{-1} \nabla \times u  \cdot \nabla \times v + 10^{-6} u \cdot v \, dx = \int_C \begin{pmatrix} y \\ -x \\ 0 \end{pmatrix} \cdot v \, dx 

from example :any:`Magnetic Fields <whetting-cmagnet>` we have to define the space

.. literalinclude:: /../py_tutorials/cmagnet.py
   :start-after: ngsglobals.msg_level
   :end-before: # u and v refer to

and the BilinearForm

.. literalinclude:: /../py_tutorials/cmagnet.py
   :start-after: nu =
   :end-before: c =
   :append: a.Assemble()

as well as the LinearForm

.. literalinclude:: /../py_tutorials/cmagnet.py
   :start-after: a.Assemble()
   :end-before: u =

The argument of the symbolic integrator must be a coefficient function depending linearly on the test and trial function.

.. automethod:: ngsolve.BilinearForm.Assemble


Nonlinear Problems
-------------------

If your left hand side of the variational formulation is nonlinear there are multiple ways to get a discretisation, depending on what you want.

Apply
^^^^^

The function ``Apply`` applies the formulation to the given ``BaseVector``. You can get a ``BaseVector`` form your ``GridFunction`` with ``GridFunction.vec``. The output vector can be created with ``BaseVector.CreateVector``.

.. automethod:: ngsolve.BilinearForm.Apply


AssembleLinearization
^^^^^^^^^^^^^^^^^^^^^^

For a variational formulation

.. math::

   \int_\Omega f(u) v \, dx

the method ``AssembleLinearization`` computes

.. math::

   \int_\Omega f'(u_\text{lin}) u v \, dx

with automatic differentiation of :math:`f(u)` and an input ``BaseVector`` :math:`u_\text{in}`.

.. automethod:: ngsolve.BilinearForm.AssembleLinearization


Assemble
^^^^^^^^^

You can do your own linearization as well using :any:`Assemble` and a ``GridFunction`` as a ``CoefficientFunction`` in your integrator. Let ``gfu_old`` be this gridfunction then

.. code-block:: python

   a = BilinearForm(fes)
   a += SymbolicBFI(gfu_old * u * v)

will be a linearization for

.. math::

   \int_\Omega u^2 v \, dx

Every time you call ``Assemble`` the bilinearform is updated with the new values of the GridFunction.


Symbolic Energy
----------------

``SymbolicEnergy`` can be used to solve a minimization problem. In :download:`this </../py_tutorials/symbolic_energy.py>` tutorial we show how to solve the nonlinear problem 

.. math::

   \min_{u \in V} 0.05 \nabla u + u^4 - 100u

For this we use ``SymbolicEnergy``:

.. literalinclude:: /../py_tutorials/symbolic_energy.py
   :start-after: u = V.TrialFunction()
   :end-before: u = GridFunction

from the ``GridFunction`` we can create new ``BaseVector``:

.. literalinclude:: /../py_tutorials/symbolic_energy.py
   :start-after: a += Symbolic
   :end-before: Draw

With this we can use :any:`AssembleLinearization` to do a Newton iteration to solve the problem:

.. literalinclude:: /../py_tutorials/symbolic_energy.py
   :start-after: Draw
