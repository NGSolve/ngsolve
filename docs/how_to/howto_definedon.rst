Spaces and forms on sub-domains
======================

We can define finite element spaces and forms on sub-domains. This is in particular useful for multi-physics.

For that, we define a geometry and assign domain lables to the different regions. In the example below we have an "outer" and an "inner" sub-domain.


We define the first H1-finite element space fes1 only on the inner sub-domain by writing 

.. code-block:: python
                
                fes1 = H1(mesh, definedon="inner")


For the definedon argument we may give one (!) string which is interpreted as a regex filter. The space is defined on all sub-domains matching the regex. Alternatively, we may give a (1-based) list of sub-domain numbers.
                

We can generate a region via

.. code-block:: python
                
                mesh.Materials("inner")
                mesh.Boundaries("b|r|t")

A region has a link to the mesh, and an array of flags defining if a sub-domain (respectively boundary part) belongs to the region. A region may be a volume region, or are boundary region. The boundary region in the statement above is the union of the boundaries labeled "b", "r", and "t".

If an integrator shall be defined only on a sub-domain, we can give the region to it:

.. code-block:: python
                
                f += SymbolicLFI (u1*v, definedon=mesh.Materials("inner"))


Why cannot we give just say:  definedon="inner" ? An integrator internally creates an array of flags specifying whether it is defined on a region, or not. The expression SymbolicFLI(.., definedon="inner") would have no information of the mesh, and could not generated this array of flags. Creating the region from the mesh provides the connection to the mesh object.
                

You find a complete example using sub-domains is here (
:download:`example</how_to/res_definedon/definedon.py>`)
                
.. literalinclude:: res_definedon/definedon.py




