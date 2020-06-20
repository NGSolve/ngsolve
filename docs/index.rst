NGS-Py Finite Element Tool
==================================

Netgen/NGSolve 6 contains a rich Python interface. Program flow as
well as geometry description and equation setup can be controlled from
Python. You should be familiar with weak formulations of partial
differential equations and the finite element method (NGSolve-oriented
lecture notes are here: `Scientific Computing <https://www.asc.tuwien.ac.at/~schoeberl/wiki/lva/notes/scicomp.pdf>`_)
and the `Python programming language <https://www.python.org>`_. The
interface to Python is inspired by the `FEniCS project <https://fenicsproject.org>`_.

We acknowledge support from the `TU Wien <https://www.tuwien.at>`_ and by the Austrian Science Foundation `FWF <https://www.fwf.ac.at/>`_ within project grant `SFB65 Taming Complexity in Partial Differential Systems <https://www.univie.ac.at/sfb65/>`_.
             
.. toctree::
  :caption: Whetting the appetite

  whetting_the_appetite/poisson.rst
  
  whetting_the_appetite/adaptive.rst
  
  whetting_the_appetite/cmagnet.rst
  
  whetting_the_appetite/navierstokes.rst
  
  whetting_the_appetite/elasticity.rst

.. * `Installation instructions using binaries <https://ngsolve.org/downloads>`_ for Windows/Mac/Linux

`Tutorial on Using NGSpy <http://web.pdx.edu/%7Egjay/teaching/mth610_2015.2/TutorialNGSpy.html>`_ by Jay Gopalakrishnan  

.. toctree::
   :maxdepth: 1
   :caption: Installation

   Download installer <https://ngsolve.org/downloads>
   install/install_sources.rst
   install/gettingstarted.rst
   install/usejupyter.rst

.. toctree::
   :maxdepth: 2
   :caption: i-Tutorials

   i-tutorials/index.rst

.. toctree::
  :maxdepth: 1
  :caption: Central Concepts

  central/meshes.rst

.. toctree::
  :maxdepth: 1
  :caption: How-to

  how_to/coefficient.rst

  how_to/howto_dirichlet.rst
  
  how_to/howto_preconditioners.rst

  how_to/howto_traceoperator.rst

  how_to/howto_linalg.rst

  how_to/howto_staticcondensation.rst

  how_to/dg.rst 
  
  how_to/howto_parallel.rst  

  how_to/howto_numpy.rst

  how_to/periodic.rst

  how_to/symbolic_integrators.rst

  how_to/howto_definedon.rst

  how_to/pml.rst

..  * [Iterating over elements]

.. include:: mylittlengs/README.rst

.. toctree::
   :maxdepth: 1
   :caption: C++ Tutorials

   mylittlengs/1_Basic/README
   mylittlengs/2_Advanced/README

.. toctree::
   :maxdepth: 1
   :caption: Netgen Tutorials

   netgen_tutorials/define_2d_geometries.rst

   netgen_tutorials/define_3d_geometries.rst

   netgen_tutorials/working_with_meshes.rst

   netgen_tutorials/manual_mesh_generation.rst

   netgen_tutorials/meshsize.rst
   
