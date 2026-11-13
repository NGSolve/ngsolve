:orphan:

.. _how-to-elasticity:

Nonlinear elasticity
====================

.. note::

   This page shows the classic Python-script version of this example.
   The current, interactive version is the i-tutorial :doc:`Solving nonlinear elasticity </i-tutorials/wta/elasticity>`.

We solve the geometric nonlinear elasticity equation using a
hyper-elastic energy density. We solve the stationary equation using
the incremental load method.

The example teaches how to

* Define a non-linear variational formulation by an energy functional
* Solve it by Newton's method
* Use a parameter for increasing the load

Download: :download:`elasticity.py</../py_tutorials/elasticity.py>`

.. literalinclude:: /../py_tutorials/elasticity.py
