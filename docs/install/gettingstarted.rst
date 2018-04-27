Getting started with Netgen/NGSolve
###################################

Ubuntu
=======
After the installation you can open a command prompt and type

.. code:: bash

   netgen

to open the GUI. To test the installation go to
   
- `Load a geometry and generate a mesh`_
- `Load and execute python files`_

Windows
=======
The netgen GUI can be started by clicking on the "netgen.exe" in your installation folder. When using the msi-installer the executatble is added to the start menu. Another option to start the GUI is to open a command prompt and type

.. code:: bash

   netgen

To test the installation go to

- `Load a geometry and generate a mesh`_
- `Load and execute python files`_

Mac OS X
========

After installing Netgen/NGSolve you find the "Netgen.app" in your Applications folder.
By setting some environment variables according to the following lines you can execute netgen from the command prompt (you should also add them to your bash start-up file).

.. literalinclude:: setenviron_mac.sh

To test the installation go to

- `Load a geometry and generate a mesh`_
- `Load and execute python files`_

Load a geometry and generate a mesh
===================================

First we need to locate the geometry files (e.g. \*.geo, \*.stl, \*.in2d). Usually they can be found in:

- Linux: ``/usr/share/netgen``
- Windows: ``$INSTALL_DIR/share/netgen``
- Mac OS X: ``/Applications/Netgen.app/Contents/Resources/share/netgen``

Now you can load a geometry with "File/Load Geometry" from the menu in the GUI and generate a mesh by clicking the "Generate Mesh" button.

If you are using the command prompt, you can start netgen and load a geometry simultaneously.

.. code:: bash

   netgen sculpture.geo

Load and execute python files
===================================

Again we need to locate the python tutorial. Usually they can be found in:

- Linux: ``/usr/share/ngsolve/py_tutorials/intro``
- Windows: ``$INSTALL_DIR/share/ngsolve/py_tutorials/intro``
- Mac OS X: ``/Applications/Netgen.app/Contents/Resources/share/ngsolve/py_tutorials/intro``
    
These python-files can be executed with "Solve/Load Python" from the menu.

Another option is to start netgen and execute a python file right from the command prompt by

.. code:: bash

   netgen poisson.py
