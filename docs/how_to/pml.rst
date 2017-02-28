Perfectly Matched Layers
========================

Perfectly Matched Layers are now implemented in NGSolve via a complex mesh 
deformation. All pre-implemented PMLs are of the form

.. math::

  \hat x(x)&:=x+i\alpha d(x)

where :math:`d(x)` is some distance function and :math:`\alpha` is the scaling 
parameter.

Pre-Implemented Scalings
------------------------

.. autofunction:: ngsolve.comp.pml.Radial


