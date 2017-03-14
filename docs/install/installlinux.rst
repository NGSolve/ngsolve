Build on Linux
##############

In case you want to change stuff in the source code or want to achieve a
local installation, you should proceed as follows. This concerns an
installation of Netgen and NGSolve.

Prerequisites
*************

Make sure that you have the following packages installed:

- You need a recent compiler, we advise gcc in version 4.9 or higher. 
- We advise to have Python installed, in version 3.4 or higher (you can compile
  netgen/ngsolve also without python support). Make sure to install
  according packages in their "dev"-version to have the suitable header
  files installed. 
- You will need tcl / tk in version >=8.5 Make sure to
  install according packages in their "dev"-version to have the suitable
  header files installed.
- git (to get the sources)
- cmake (>=2.8.12) for the build system
- libxmu-dev (you might also need the xorg-dev package)
- libglu (again, the "dev"-version)
- liblapacke-dev

The following line (executed as root) should update/install all prerequisites:

.. code:: bash
	    
   apt-get update && apt-get -y install python3 libpython3-dev libxmu-dev tk-dev tcl-dev cmake git g++ libglu1-mesa-dev liblapacke-dev

Getting the source
******************

Choose a directory where you want to install Netgen/NGSolve, which will be denoted as ``${BASEDIR}`` in the following installation guide. Change into your base directory and clone the git repository.

.. code:: bash

   export BASEDIR=~/ngsuite
   mkdir -p $BASEDIR
   cd $BASEDIR
   git clone git://git.code.sf.net/p/ngsolve/git ngsolve-src

To also fetch the dependencies (Netgen) we must tell git to load the
submodules

.. code:: bash

   cd $BASEDIR/ngsolve-src
   git submodule update --init --recursive

Building from the source
************************

Now we create folders for builds and installations.

.. code:: bash

   mkdir $BASEDIR/ngsolve-build
   mkdir $BASEDIR/ngsolve-install
	  
Configuring
===========

Change into the directory for builds and call cmake with a link to the source directory

.. code:: bash

   cd $BASEDIR/ngsolve-build
   cmake XYZ ${BASEDIR}/ngsolve-src

and parameters XYZ, which can be used to set a lot of options. You can also use
ccmake or cmake-gui (if installed) to select the options there.
Important options are 

- The installation directory, here ``-DINSTALL_DIR=${BASEDIR}/ngsolve-install``
- Release type, here a reasonable choice is (that's default)
  ``-DCMAKE_BUILDTYPE=RelWithDebInfo``
- If you have intel mkl installed you have to activate it with
  ``-DUSE_MKL=ON -DMKL_ROOT=/opt/intel/mkl/`` where you should replace the
  ``MKL_ROOT`` with your correct path

So your configuring command could look like this:

.. code:: bash

    cd $BASEDIR/ngsolve-build && cmake -DINSTALL_DIR=${BASEDIR}/ngsolve-install ${BASEDIR}/ngsolve-src
          -DUSE_MKL=ON -DMKL_ROOT=/opt/intel/mkl/

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

Finally you have to set the environment variable ``NETGENDIR`` to the
location of the executable, eg. by

.. code:: bash

    export NETGENDIR="${BASEDIR}/ngsolve-install/bin"

or

.. code:: bash

    setenv NETGENDIR "${BASEDIR}/ngsolve-install/bin"

(depends on your linux distribution). You may want to add the
corresponding line to your .bashrc, s.t. it is automatically set
whenever you log in into your bash-shell. To be able to start netgen
from the command line, you have to add ``NETGENDIR`` to the ``PATH``, eg. by

.. code:: bash

    export PATH=$NETGENDIR:$PATH

or

.. code:: bash

    setenv PATH "$NETGENDIR:$PATH"

When you want to run Python scripts, you also have to set ``PYTHONPATH`` to
the appropriate directory:

.. code:: bash

    export PYTHONPATH=$NETGENDIR/../`python3 -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(1,0,''))"`

Test the installation
=====================

Netgen
------

Now the installation should be finished. Test it with calling netgen

.. code:: bash

    netgen

in ``${BASEDIR}/ngsolve-install/share/netgen`` you can find several geometry and mesh
files which you can use to try if netgen does what it should do.

NGSolve
-------

Test NGSolve with calling netgen

.. code:: bash

    netgen

and see if you get a message saying that the module NGSolve-|version|-dev has
been loaded. In ``${BASEDIR}/ngsolve-install/share/ngsolve`` you can find example PDE
problems which you can use to try if Netgen/NGSolve does what it should do.

Everything together
*******************

For brave ones here is a complete copy/paste script to download,
configure, build and install Netgen/NGSolve when all dependencies are
already installed. Note that bash is assumed as environment here. There
is absolutely no warranty that this will work, so use it at your own
risk! In case this script is failing, please follow the procedure above
before asking for help.

.. code:: bash

   export BASEDIR=~/ngsuite
   mkdir -p $BASEDIR
   cd $BASEDIR && git clone git://git.code.sf.net/p/ngsolve/git ngsolve-src
   cd $BASEDIR/ngsolve-src && git submodule update --init --recursive
   mkdir $BASEDIR/ngsolve-build
   mkdir $BASEDIR/ngsolve-install
   cd $BASEDIR/ngsolve-build
   cmake -DINSTALL_DIR=${BASEDIR}/ngsolve-install ${BASEDIR}/ngsolve-src
   make -j4
   make install
   echo "export NETGENDIR=${BASEDIR}/ngsolve-install/bin" >> ~/.bashrc
   echo "export PATH=\$NETGENDIR:\$PATH" >> ~/.bashrc
   export PYTHONPATH_TMP=`python3 -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(1,0,''))"`
   echo "export PYTHONPATH=\$NETGENDIR/../${PYTHONPATH_TMP}:\$PATH" >> ~/.bashrc
   source ~/.bashrc
   cd ${BASEDIR}/ngsolve-install/share/ngsolve/py_tutorials/intro
   netgen navierstokes.py

Keep it up to date
******************

To update Netgen/NGSolve go to the source directory (``${BASEDIR}/ngsolve-src``) and fetch the latest
sources. The second line is needed to update all dependencies provided
as git submodules (such as Netgen).

.. code:: bash

   git pull
   git submodule update --recursive --init

After that, go to the build directory (``${BASEDIR}/ngsolve-build``) and build/install again

.. code:: bash

   make
   make install
