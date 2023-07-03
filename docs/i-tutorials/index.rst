
Interactive NGSolve Tutorial
============================

What are the i-tutorials ?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The i-tutorials are interactive tutorials to NGS-Py, the Python
front-end to NGSolve. The i-tutorials are Jupyter notebooks which allow
you to explore the features of NGS-Py.

The i-tutorials have been setup for the 2017 NGSolve user meeting. The
authors of the sections are Jay Gopalakrishnan (Getting Started),
Joachim Schöberl (Advanced Topics), Christoph Lehrenfeld (Time-dependent
and non-linear problems), Christoph Wintersteiger (Geometric modeling
and mesh generation). Big thanks to Matthias Hochsteger for integrating
the Netgen-GUI into Jupyter.

Copyright: The i-tutorials are a part of NGSolve and are covered by the
LGPL open source license. You may extend and modify them for your use,
but you have to refer to the original source.

We acknowledge support from the `TU Wien <https://www.tuwien.at>`_ and by the Austrian Science Foundation `FWF <https://www.fwf.ac.at/>`_ within project grant `SFB65 Taming Complexity in Partial Differential Systems <https://www.univie.ac.at/sfb65/>`_.


Installation
~~~~~~~~~~~~

To work through the i-tutorials, you first have to install

- Netgen/NGSolve, see http://www.ngsolve.org/downloads 
- Jupyter from http://www.jupyter.org. You can use the pip package manager: "pip3 install jupyter" 
- download and unpack the i-tutorials from `here <../i-tutorials.zip>`__.

To use the webgui visualization within jupyter you need

.. code-block:: bash

    pip3 install webgui_jupyter_widgets
    jupyter nbextension install --user --py widgetsnbextension
    jupyter nbextension enable --user --py widgetsnbextension
    jupyter nbextension install --user --py webgui_jupyter_widgets
    jupyter nbextension enable --user --py webgui_jupyter_widgets

  
Some of the tutorials require packages from scipy and matplotlib, so it
is a good idea to install them as well:

.. code-block:: bash

    pip3 install scipy matplotlib matplotlib

i-tutorials on Youtube
----------------------

We have recorded the tutorial sessions from the 3rd NGSolve Usermeeting in which most of these tutorial files were presented. You can watch the presentations on the `NGSolve Youtube channel <https://www.youtube.com/playlist?list=PL_5FauasEdy01ftt2BH9XCsDjgW_moaBi>`__.

Starting
~~~~~~~~

You start with the interactive tutorial by opening a terminal, go to the
main folder containing the i-tutorials, and start

.. code-block:: bash

   jupyter-notebook index.ipynb

Whetting the Appetite
=====================

.. toctree::
   :maxdepth: 1

   wta/poisson.ipynb
   wta/adaptivity.ipynb
   wta/maxwell.ipynb 
   wta/coil.ipynb 
   wta/navierstokes.ipynb
   wta/elasticity.ipynb 
   wta/elasticity3D.ipynb 

   
1. Getting started
===================

.. toctree::
   :maxdepth: 1
             
   unit-1.1-poisson/poisson.ipynb
   unit-1.2-coefficient/coefficientfunction.ipynb
   unit-1.3-dirichlet/dirichlet.ipynb
   unit-1.4-staticcond/staticcond.ipynb
   unit-1.5-subdomains/subdomains.ipynb
   unit-1.6-adaptivity/adaptivity.ipynb
   unit-1.7-helmholtz/helmholtz.ipynb
   unit-1.7-helmholtz/pml.ipynb
   unit-1.8-meshtopology/meshtopology.ipynb
   
2. Advanced Topics
==================

.. toctree::
   :maxdepth: 1

   unit-2.1.1-preconditioners/preconditioner.ipynb
   unit-2.1.2-blockjacobi/blockjacobi.ipynb
   unit-2.1.3-bddc/bddc.ipynb
   unit-2.2-eigenvalues/pinvit.ipynb
   unit-2.3-hcurlhdiv/hcurlhdiv.ipynb
   unit-2.4-Maxwell/Maxwell.ipynb
   unit-2.4-Maxwell/Maxwellevp.ipynb
   unit-2.5-mixed/mixed.ipynb
   unit-2.6-stokes/stokes.ipynb
   unit-2.7-hybrid/hybrid.ipynb
   unit-2.8-DG/DG.ipynb
   unit-2.9-fourthorder/fourthorder.ipynb
   unit-2.10-dualbasis/dualbasis.ipynb
   unit-2.11-matrixfree/matrixfree.ipynb

3. Time-dependent and non-linear problems
==========================================

.. toctree::
   :maxdepth: 1

   unit-3.1-parabolic/parabolic.ipynb
   unit-3.2-navierstokes/navierstokes.ipynb
   unit-3.3-scalardg/scalardg.ipynb
   unit-3.3.1-wavedg/wavedg.ipynb
   unit-3.4-simplehyp/shallow2D.ipynb
   unit-3.5-surfacehdg/surfacehdg.ipynb
   unit-3.6-opsplit/opsplit.ipynb
   unit-3.7-nonlinear/nonlinear.ipynb
   unit-3.8-nonlmin/nonlmin.ipynb

4. Geometric modeling and mesh generation
==========================================

.. toctree::
   :maxdepth: 1

   unit-4.1.1-geom2d/geom2d.ipynb
   unit-4.1.2-csg2d/csg2d.ipynb
   unit-4.2-csg/csg.ipynb
   unit-4.3-manualmesh/manualmeshing.ipynb 
   unit-4.4-occ/occ.ipynb
   unit-4.4-occ/bottle.ipynb
   unit-4.4-occ/workplane.ipynb               


5. MPI - Parallelization and CUDA Support
=========================================

.. toctree::
   :maxdepth: 1

   unit-5a.1-mpi/poisson_mpi.ipynb 
   unit-5a.2-pardofs/pardofs.ipynb 
   unit-5a.3-petsc/petsc.ipynb 
   unit-5a.3-petsc/PETSc_interface.ipynb

CUDA Device support:

.. toctree::
   :maxdepth: 1

   unit-5.5-cuda/poisson_cuda.ipynb 
   unit-5.5-cuda/wave_cuda.ipynb 
   unit-5.5-cuda/EulerEquations.ipynb 

some more MPI examples:

.. toctree::
   :maxdepth: 1

   historic/unit-5.0-mpi_basics/MPI-Parallelization_in_NGSolve.ipynb
   historic/unit-5.1-mpi_ngsolve/mpi_basics.ipynb
   historic/unit-5.2-fetidp_point2d/feti-dp-i.ipynb
   historic/unit-5.3-fetidp_point3d/feti-dp-ii.ipynb
   historic/unit-5.4-fetidp_edge/feti-dp-iii.ipynb
   historic/ unit-5.5-fetidp_inexact/feti-dp-iv.ipynb
   
   
6. Various Topics
=========================

.. toctree::
    :maxdepth: 1
 
    unit-6.1.1-surfacemeshes/surface_meshes.ipynb
    unit-6.1.2-surfacepde/surface_pdes.ipynb
    unit-6.1.3-rmplate/Reissner_Mindlin_plate.ipynb
    unit-6.1.4-shells/shell.ipynb
    
    unit-6.2-contact/contact.ipynb 
    unit-6.3-plasticity/plasticity.ipynb


7.  Shape- and Topology Optimization
=====================================
Peter Gangl and Kevin Sturm

.. toctree::
   :maxdepth: 1

   unit-7-optimization/01_Shape_Derivative_Levelset.ipynb
   unit-7-optimization/02_Shape_Derivative_Laplace.ipynb
   unit-7-optimization/03_Shape_Derivative_Laplace_SemiAuto.ipynb
   unit-7-optimization/03b_Shape_Derivative_Laplace_FullyAuto.ipynb

   unit-7-optimization/04_Topological_Derivative_Levelset.ipynb
   unit-7-optimization/05_Topological_Derivative_Transmission.ipynb


8. Unfitted Finite Elements
===========================
C. Lehrenfeld and the `ngsxfem <https://github.com/ngsxfem/ngsxfem/>`_ authors

These units require the Add-on `ngsxfem <https://github.com/ngsxfem/ngsxfem/>`_ to be installed.
There are further ngsxfem-tutorials `here <https://github.com/ngsxfem/ngsxfem-jupyter/>`_.

.. toctree::
   :maxdepth: 1

   unit-8.1-basics/basics.ipynb
   unit-8.2-intlset/intlset.ipynb
   unit-8.3-cutfem/cutfem.ipynb
   unit-8.4-spacetime_fitted/spacetime_fitted.ipynb
   unit-8.5-spacetime_unfitted/spacetime_unfitted.ipynb
   unit-8.6-mlset_basic/mlset_basic.ipynb
   unit-8.7-mlset_pde/mlset_pde.ipynb
   unit-8.8-aggregation/aggregation.ipynb
   unit-8.9-unfmixed/unfmixed.ipynb

9. Extending by C++ programming
===============================

.. toctree::
   :maxdepth: 1

   unit-9.1-C++FE/CppExtension.ipynb
   unit-9.2-C++Assemble/cppassembling.ipynb
   unit-9.3-highorder/highorder.ipynb
   

10. NGSolve and ...
===================

.. toctree::
   :maxdepth: 1
              
   unit-10.1-ngspice/NGSpiceNGSolve.ipynb
   unit-10.2-tensorflow/TensorFlowNGSolve.ipynb


Appendix
========

.. toctree::
  :maxdepth: 1
  
  appendix-taskmanager/taskmanager.ipynb 
  appendix-webgui/webgui-internal.ipynb
  
       
