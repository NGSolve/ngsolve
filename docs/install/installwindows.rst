Build on Windows
################

In case you want to change stuff in the source code or want to achieve a
local installation, you should proceed as follows. The recommended way
is to checkout the source code of NGSolve. During the process Netgen
will be automatically downloaded and built.

Prerequisites
*************

Make sure that you have the following packages installed

- You need a recent compiler, we advise Microsoft Visual Studio 2017 or newer
- Python 3: You can get it `here <https://www.python.org/downloads/windows/>`__ if
  you don't already have it. Make sure to download one of the "x86-64"
  installers. Also make sure to select the option "Add Python 3 to PATH"
  during the installation.
- `cmake <http://www.cmake.org>`__ to generate
  Visual Studio Projects Files (Win32 Installer works for both
  architectures) we recommend the install option 'Add CMake to the system
  PATH for all users'
- `git <https://git-scm.com/downloads>`__ for downloading Netgen/NGSolve sources

Directory structure
*******************

For ease of presentation we use a directory structure with three main
directories:

- "src" for the source code NGSolve
- "build" for the output of CMake and the compiler
- "install" for the actual installation

In the following we assume that you have chosen a base directory where
three (empty) directories "src", "build" and "install" have been
created. This base directory will be denoted as ``BASEDIR`` in the
following.

Getting the source
******************

Start "Git Bash", navigate to ``BASEDIR`` and execute

.. code:: bash

   git clone --recurse-submodules https://github.com/NGSolve/ngsolve.git src

Go to the source directory and initialize the git submodules to download
the source files for Netgen and Pybind11.

.. code:: bash

       cd src
       git submodule update --init --recursive
       cd ..

After this step you should have the source code of Netgen in
src/external\_dependencies/netgen.

Configuring
===========

Next, change to the "build" directory and use cmake to configure from
the sources. The standard configuring command looks like this:

.. code:: bash

       mkdir build
       cd build
       cmake "../src" -DCMAKE_INSTALL_PREFIX="BASEDIR/install" -G "Visual Studio 15 Win64"

There are many options for cmake, which you can find using

.. code:: bash

    cmake --help

Building
========

To build Netgen/NGSolve, exectue

.. code:: bash

       cmake --build . --config Release --target install

Alternatively you can open ``BASEDIR\build\SUPERBUILD.sln`` in
Visual Studio and build the "INSTALL" project. Make sure to select the
desired solution configuration (e.g. Release).

Finishing the installation
==========================

Finally you have to set the environment variables PATH, NETGENDIR and
PYTHONPATH to the appropriate values. This is done by building the
"set\_environment\_variables" target:

.. code:: bash

       cmake --build . --config Release --target set_environment_variables
