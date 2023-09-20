Build on Linux
##############

In case you want to change stuff in the source code or want to achieve a
local installation, you should proceed as follows. This concerns an
installation of Netgen and NGSolve.

Prerequisites
*************

Make sure that you have the following packages installed:

- You need a recent compiler (GCC 7 or later)
- We advise to have Python installed, in version 3.4 or higher (you can compile
  Netgen/NGSolve also without python support). Make sure to install
  according packages in their "dev"-version to have the suitable header
  files installed. 
- You will need tcl / tk in version >=8.5 Make sure to
  install according packages in their "dev"-version to have the suitable
  header files installed.
- git (to get the sources)
- cmake (>=3.3) for the build system
- libxmu-dev (you might also need the xorg-dev package)
- libglu (again, the "dev"-version)
- liblapacke-dev

The following line should update/install all prerequisites on Ubuntu (you need root privileges):

.. code:: bash
	    
   sudo apt-get update && sudo apt-get -y install python3 python3-distutils python3-tk libpython3-dev libxmu-dev tk-dev tcl-dev cmake git g++ libglu1-mesa-dev liblapacke-dev libocct-data-exchange-dev libocct-draw-dev occt-misc libtbb-dev libxi-dev

Getting the source
******************

Choose a directory where you want to install Netgen/NGSolve, which will be denoted as ``${BASEDIR}`` in the following installation guide. Change into your base directory and clone the git repository.

.. code:: bash

   export BASEDIR=~/ngsuite
   mkdir -p $BASEDIR
   cd $BASEDIR
   git clone --recurse-submodules https://github.com/NGSolve/ngsolve.git ngsolve-src

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
   cmake -DCMAKE_INSTALL_PREFIX=${BASEDIR}/ngsolve-install ${BASEDIR}/ngsolve-src

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


Everything together
*******************

For brave ones here is a complete copy/paste script to download,
configure, build and install Netgen/NGSolve when all dependencies are
already installed. Note that bash is assumed as environment here. There
is absolutely no warranty that this will work, so use it at your own
risk! In case this script is failing, please follow the procedure above
before asking for help.

.. literalinclude:: compilengsolve.sh

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
