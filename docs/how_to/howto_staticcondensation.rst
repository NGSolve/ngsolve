Static condensation of internal bubbles
=======================================

Element-internal unknowns of higher order finite elements can be eliminated from the global linear system. This corresponds to the Schur-complement system for the coupling unknowns on the element boundaries.

To enable this option, set the `eliminate_internal` flag for the bilinear-form. This assembles the global Schur-complement system, and stores harmonic extensions, and internal inverses.

The user is responsible to transform the right-hand side and the solution vector. An example is as follows:

.. code-block:: python

                from netgen.geom2d import unit_square
                from ngsolve import *

                mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))

                fes = H1(mesh, order=10, dirichlet=[1,2])
                u = fes.TestFunction()
                v = fes.TrialFunction()
                
                a = BilinearForm(fes, flags = { "eliminate_internal" : True } )
                a += SymbolicBFI (grad(u) * grad(v))
                a.Assemble()
                
                f = LinearForm(fes)
                f += SymbolicLFI (1 * v)
                f.Assemble()
                
                u = GridFunction(fes)

modify right-hand side

.. code-block:: python

                f.vec.data += a.harmonic_extension_trans * f.vec

solve for external unknowns 

.. code-block:: python

                u.vec.data = a.mat.Inverse(fes.FreeDofs(True)) * f.vec

and find element-internal solution

.. code-block:: python

                u.vec.data += a.harmonic_extension * u.vec
                u.vec.data += a.inner_solve * f.vec

                Draw (u)








