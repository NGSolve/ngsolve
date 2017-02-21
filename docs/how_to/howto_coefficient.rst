Working with CoefficientFunctions
=================================


A CoefficientFunction is a function which can be evaluated on a mesh, and may be used to provide the coefficient or a right-hind-side to the variational formulation. As typical finite element procedures iterate over elements, and map integration points from a reference element to a physical element, the evaluation of a CoefficientFunction requires a mapped integration point. The mapped integration point contains the coordinate, as well as the Jacobian, but also element number and region index. This allows the efficient evaluation of a wide class of functions needed for finite element computations.

Some examples of using coefficient functions are:

.. code-block:: python
  
  a += Laplace (CoefficientFunction([1,17,5]))
  f += Source (x*(1-x))
  
A *CoefficientFunction* initialized with an array stores one value per region, which may be a domain or a boundary region. Note that domain and boundary indices have 1-based counting. The Cartesian coordinates *x*, *y*, and *z* are pre-defined coordinate coefficient functions. Algebraic operations as - or * combine existing CoefficientFunctions to new ones.

A CoefficientFunction may be real or complex valued, and may be scalar or vector-valued. 
(Matrix valued CoefficientFunctions will be available in future versions.) Scalar CFs can be combined to vectorial CF by creating a new object from a tuple, and components of a vectorial CF can be accessed by the bracket operator:

.. code-block:: python
                
                fx = CoefficientFunction (0)
                fy = CoefficientFunction ([1,2,3])
                vec_cf = CoefficientFunction( (fx,fy) )
                fx_again = vec_cf[0]


Since CFs require a mapped integration-point as argument, we first
have to generate one. Calling the mesh-object with x, y, and
optionally z coordinates generates one for us. Note, this requires a
global search for the element containing the point:

.. code-block:: python
                
  >>> mip = mesh(0.2,0.4)
  >>> (x(mip), y(mip))
  (0.2, 0.4)


Also a GridFunction is a CoefficientFunction by inheritance, and can
be used whenever a CF is required. If the finite element space
provides differentiation, then

.. code-block:: python

                u.Deriv()

creates a new CoefficientFunction. Depending on the space, *Deriv* means the gradient, the curl, or the divergence. Direct evaluation of higher order derivatives is not supported, but can be approximated via interpolation:


.. code-block:: python

                V1 = H1(mesh,order=5)
                V2 = H1(mesh,order=4,dim=2)
                u = GridFunction(V1)
                du = GridFunction(V2)
                u.Set(y*sin(3.1415*x))
                du.Set(u.Deriv())
                uxx,uxy,uyx,uyy = du.Deriv()
                Draw (uxx, mesh, "uxx")
                Draw (uxy-uyx, mesh, "not_perfect")


The symbolic definition of bilinear- and linear-forms uses also
CoefficientFunctions, as in this example (*b* is a regular, vectorial
CF):

.. code-block:: python

                u = V.TrialFunction()
                v = V.TestFunction()
                a += SymbolicBFI (b * u.Deriv() * v)

Here, *u* and *v* are so called ProxyFunctions, which do not evaluate to values, but are used as arguments in the definition of the forms. Also *u.Deriv()* is a ProxyFunction. When the form is evaluated in an integration point, the ProxyFunctions are set to unit-vectors, and the CoefficientFunction evaluates to an actual number.

