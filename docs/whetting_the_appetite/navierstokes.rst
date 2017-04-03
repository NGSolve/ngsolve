.. _how-to-navierstokes:

Navier Stokes Equation
======================

We solve the time-dependent incompressible Navier Stokes Equation. For
that

* we use the P3/P2 Taylor-Hood mixed finite element pairing
* and perform operator splitting time-integration with the non-linear term explicit, but time-dependent Stokes fully implicit.

The example is from the Schäfer-Turek `benchmark <http://www.mathematik.tu-dortmund.de/lsiii/cms/papers/SchaeferTurek1996.pdf>`
a two-dimensional cylinder, at Reynolds number 100

Download :download:`navierstokes.py</../py_tutorials/navierstokes.py>`

.. literalinclude:: /../py_tutorials/navierstokes.py

The absolute value of velocity:

.. image:: res_navierstokes/navierstokes_solution.jpg
