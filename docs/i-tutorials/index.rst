
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

- pip3 install ipywidgets
- jupyter nbextension install --py widgetsnbextension
- jupyter nbextension enable --py widgetsnbextension
- jupyter nbextension install --user --py ngsolve
- jupyter nbextension enable --user --py ngsolve

  
Some of the tutorials require packages from scipy and matplotlib, so it
is a good idea to install them as well:

-  pip3 install scipy
-  pip3 install matplotlib

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
   wta/navierstokes.ipynb
   wta/elasticity.ipynb

   
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

   unit-4.1-geom2d/geom2d.ipynb
   unit-4.2-csg/csg.ipynb
   unit-4.3-manualmesh/manualmeshing.ipynb


5. MPI - Parallelization
=========================

.. toctree::
   :maxdepth: 1

   unit-5.0-mpi_basics/MPI-Parallelization in NGSolve.ipynb
   unit-5.1-mpi_ngsolve/mpi_basics.ipynb
   unit-5.2-fetidp_point2d/feti-dp-i.ipynb
   unit-5.3-fetidp_point3d/feti-dp-ii.ipynb
   unit-5.4-fetidp_edge/feti-dp-iii.ipynb
   unit-5.5-fetidp_inexact/feti-dp-iv.ipynb
   unit-5a.1-mpi/poisson_mpi.ipynb 
   unit-5a.2-pardofs/pardofs.ipynb 
   unit-5a.3-petsc/petsc.ipynb 
   
   
6. Various Topics
=========================

.. toctree::
    :maxdepth: 1
 
    unit-6.1.1-surfacemeshes/surface_meshes.ipynb
    unit-6.1.2-surfacepde/surface_pdes.ipynb
    unit-6.1.3-rmplate/Reissner_Mindlin_plate.ipynb
    unit-6.1.4-shells/shell.ipynb
