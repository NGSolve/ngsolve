Perfectly Matched Layers
========================

Perfectly Matched Layers are now implemented in NGSolve via a complex mesh 
deformation. All pre-implemented PMLs are of the form

.. math::

  \hat x(x)&:=x+i\alpha d(x)

where :math:`d(x)` is some distance function and :math:`\alpha` is the scaling 
parameter.




Creating a PML Transformation
----------------------------------

A PML transformation is a mesh independent object, which can be created using its
generator function defined in the python module ``comp.pml``.

.. autoclass:: ngsolve.comp.pml.PML
  :members:




Pre-implemented scalings
^^^^^^^^^^^^^^^^^^^^^^^^^^

The following scalings are available by default:

.. autofunction:: ngsolve.comp.pml.Radial

.. autofunction:: ngsolve.comp.pml.Cartesian

.. autofunction:: ngsolve.comp.pml.BrickRadial

.. autofunction:: ngsolve.comp.pml.HalfSpace


Creating your own scaling
^^^^^^^^^^^^^^^^^^^^^^^^^

Aside from the pre-implemented scalings, one can also create scalings using
coefficient functions or combine available PMLs.

.. autofunction:: ngsolve.comp.pml.Custom

.. autofunction:: ngsolve.comp.pml.Compound

PML transformations can be added using the ``+`` operator.


Applying the Transformation to a Mesh
-------------------------------------------

A PML object ``pmlobj`` can be applied to a mesh ``m`` on domain ``domain`` 
using

.. code-block:: python

  m.SetPML(pmlobj,'domain')

After this is done all (symbolic) integrators will respect the additional complex
transformation.
Note that right now PMLs are tested on H1 spaces only.

.. caution::

  Evaluating coordinate coefficient functions on complex mapped 
  integration points will result only in evaluation of their real part.
  Thus, using non constant coefficient functions in the PML domain should be
  handled with caution.




