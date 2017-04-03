The Trace() operator
====================

Mathematically, the trace operator restricts a domain-function to the boundary of the domain. For function spaces like H(curl) or H(div), the trace operator delivers only tangential, and normal components, respectively.

The trace operator is used in NGS-Py as follows:

.. code-block:: python

                a.Assemble(u.Trace()*v.Trace(), BND)
                
The evaluation of boundary values involves only degrees of freedom geometrically located on the boundary.

Traces of derivatives can are formed by either one of the following.
For example, the derivative of the trace of an H1-function gives the tangential derivative:


.. code-block:: python

                u.Trace().Deriv()
                u.Deriv().Trace()


As an popular exception, in NGS-Py, the Trace()-operator for H1-functions is optional.

For element-boundary integrals, the Trace()-operator must not be used. Here, the function is evaluated at a point in the volume-element.








