.. |Python| replace:: Python 3.8
.. _Python: https://www.python.org/downloads/mac-osx/
		    
Build on Mac OS X.
##################

This page describes how to build Netgen and NGsolve on Apple computers
running OS X 10.9 or greater. These steps will give you a native OSX
style *framework* and a *Netgen.app*. These steps have been tested on
MacOS version 10.12.

Here we describe how to build the latest version |version|.

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
  install `CMake <http://www.cmake.org>`_. Make sure to
  install the command line tools: Open CMake, click on "How to Install For Command Line Use"
  in the "Tools" menu and execute one of the suggested options.

* Install |Python|_


Getting the source
******************

Create a folder where your Netgen/NGsolve sources and builds will reside.
  
.. code:: bash
   
   export NGROOT=<PathToTheFolderYouJustMade>

Move to the folder, git clone the NGSolve repository.

.. code:: bash
   
   cd $NGROOT
   git clone https://github.com/NGSolve/ngsolve.git ngsolve-src
	    
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
   cmake $NGROOT/ngsolve-src




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

.. literalinclude:: setenviron_mac.sh

and execute the file with

.. code:: bash

   source .bashrc

to set all environment variable needed to start Netgen/NGSolve from the command line.
