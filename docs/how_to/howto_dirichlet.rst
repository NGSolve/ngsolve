Setting inhomogeneous Dirichlet boundary conditions
=========================================================

A mesh stores boundary elements, which know the *bc* name given in the geometry.
The Dirichlet boundaries are given as a regular expression of these names to the finite element space:

.. code-block:: python

                V = FESpace(mesh,order=3,dirichlet="top|back")
                u = GridFunction(V)


The BitArray of free (i.e. unconstrained) dof numbers can be obtained via

.. code-block:: python

                freedofs = V.FreeDofs()
                print (freedofs)


Inhomogeneous Dirichlet values are set via

.. code-block:: python

                u.Set(x*y, definedon=mesh.Boundaries("top|back"))

This function performs an element-wise L2 projection combined with arithmetic averaging of coupling dofs.

As usual we define biform a and liform f. Here, the full Neumann matrix and non-modified right hand sides are stored.

Boundary constraints are treated by the preconditioner. For example, a Jacobi preconditioner created via

.. code-block:: python

                c = Preconditioner(a, "local")

inverts only the unconstrained diagonal elements, and fills the remaining values with zeros.

The boundary value solver keeps the Dirichlet-values unchanged, and solves only for the free values

.. code-block:: python

                BVP(bf=a,lf=f,gf=u, pre=c).Do()


A do-it-yoursolve version of homogenization is:

.. code-block:: python

                res = f.vec.CreateVector()
                res.data = f.vec - a.mat * u.vec
                u.vec.data += a.mat.Inverse(v.FreeDofs()) * res


