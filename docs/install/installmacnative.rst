Build on Mac OS X.
##################

This page describes how to build Netgen and NGsolve on Apple computers
running OS X 10.9 or greater. These steps will give you a native OSX
style *framework* and a *Netgen.app*. These steps have been tested on
MacOS version 10.12.

Here we describe how to build the latest version |version|-dev.

Prerequisites
*************

* You need to have Xcode installed. We will build using the default
  ``clang`` compiler. If it is not in the default location/path, please
  set

  .. code:: bash
     
     export CC=<MyPathtoClang>
     export CXX=<MyPathtoClang++>

* You need to have ``make`` installed for command line use. Xcode
  Command Line Tools package comes with make. You can install Xcode
  Command Line Tools by ``xcode-select --install``.

* You need to have ``cmake`` installed. You can either download and
  install cmake from source or install the CMake App from
  http://www.cmake.org. If you install the CMake App, make sure to
  install the command line tools. Open CMake, click on "How to Install For Command Line Use"
  in the "Tools" menu and execute one of the suggested options.

* Install `Python |python_version| <https://www.python.org/downloads/mac-osx/>`_


Getting the source
******************

Create a folder where your Netgen/NGsolve sources and builds will reside.
  
.. code:: bash
   
   export NGROOT=<PathToTheFolderYouJustMade>

Move to the folder, git clone the NGSolve repository.

.. code:: bash
   
   cd $NGROOT
   git clone git://git.code.sf.net/p/ngsolve/git ngsolve-src
	    
To also fetch the dependencies (Netgen) we must tell git to load the submodules

.. code:: bash
   
   cd $NGROOT/ngsolve-src
   git submodule update --init --recursive
	    
Building from the source
************************

Create a folder for the builds.

.. code:: bash
   
   mkdir $NGROOT/ngsolve-build


Configuring
===========

Change into the directory for builds and call cmake with a link to the source directory

.. code:: bash

   cd $NGROOT/ngsolve-build
   cmake XYZ $NGROOT/ngsolve-src

and cmake parameters XYZ, which can be used to set a lot of options.

Building
========

Now, call

.. code:: bash

    make 

You may want to add "-jx" with x the number of threads you want to use
for the compilation. If everything goes smooth you can install the
resulting build calling

.. code:: bash

   make install

Finishing the installation
==========================

Add the following line to your ``.bashrc`` file in your home directory

.. code:: bash

   export PYTHONPATH=$PYTHONPATH:/Applications/Netgen.app/Contents/Resources/lib/python|python_version|/site-packages:.
   export NETGENDIR=/Applications/Netgen.app/Contents/MacOS
   export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$NETGENDIR
   export DYLD_FRAMEWORK_PATH=$DYLD_FRAMEWORK_PATH:$NETGENDIR/../Frameworks
   export PATH=$NETGENDIR:$PATH
	  
and execute the file with

.. code:: bash

   source .bashrc

to set all environment variable needed to start Netgen/NGSolve from the command line.

Test the installation
=====================

Navigate to your "Applications" folder and start "Netgen".
Now you can:

* load geometries ("File/Load Geometry") from ``Resources/share/netgen`` and generate a mesh
* load PDE-files ("Solve/Load PDE") from ``Resources/share/ngsolve`` and solve
* execute python-files ("Solve/Load Python") ``Resources/share/ngsolve/py_tutorials/intro``
    

Test the installation from the command line
===========================================

Netgen
------
Now the installation should be finished. Test it with calling netgen

.. code:: bash

    netgen

in ``/Applications/Netgen.app/Contents/Resources/share/netgen`` you can find several geometry and mesh
files which you can use to try if netgen does what it should do.

NGSolve
-------
Test NGSolve with calling netgen

.. code:: bash

    netgen

and see if you get a message saying that the module NGSolve-|version|-dev has
been loaded. In ``/Applications/Netgen.app/Contents/Resources/share/ngsolve`` you can find example PDE
problems which you can use to try if Netgen/NGSolve does what it should do.
