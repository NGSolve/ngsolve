.. _how-to-elasticity:

Nonlinear elasticity
====================

We solve the geometric nonlinear elasticity equation using a
hyper-elastic energy density. We solve the stationary equation using
the incremental load method.

The example teaches how to

* Define a non-linear variational formulation using SymbolicEnergy
* Solve it by Newton's method
* Use a parameter for increasing the load

Download: :download:`elasticity.py</../py_tutorials/elasticity.py>`

.. literalinclude:: /../py_tutorials/elasticity.py
