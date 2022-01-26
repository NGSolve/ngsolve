.. role:: scrollable
	  
CMake options for building Netgen/NGSolve
*****************************************

.. table::
   :class: rows

   +-----------------------------+--------+-------------------+-----------------------------------------------------+
   |        CMake flag           | Default| Description                                                             |
   +=============================+========+===================+=====================================================+
   | CMAKE_INSTALL_PREFIX        |        | Directs the installation to this folder. The default values are         |
   +                             +        +-------------------+-----------------------------------------------------+
   |                             |        | Linux             | /opt/netgen                                         |
   +                             +        +-------------------+-----------------------------------------------------+
   |                             |        | Windows           | C:\\netgen                                          |
   +                             +        +-------------------+-----------------------------------------------------+
   |                             |        | MacOS             | /Applications/Netgen.app                            |
   +-----------------------------+--------+-------------------+-----------------------------------------------------+
   |  USE_MKL                    | OFF    | Some Intel MKL function are used (e.g. Lapack, Pardiso,...)             |
   +                             +        +-----------------------+-------------------------------------------------+
   |                             |        | MKL_ROOT              | Path to MKL installation (e.g. /opt/intel/mkl)  |
   +                             +        +-----------------------+-------------------------------------------------+
   |                             |        | MKL_STATIC            | Link statically if set to ON (OFF by default)   |
   +                             +        +-----------------------+-------------------------------------------------+
   |                             |        | MKL_SDL               | Link with single dynamic library                |
   |                             |        |                       | (only if MKL_STATIC is OFF)                     |
   +                             +        +-----------------------+-------------------------------------------------+
   |                             |        | MKL_INTERFACE_LAYER   | Link with single dynamic library                |
   |                             |        |                       | (only if MKL_STATIC is OFF)                     |
   +-----------------------------+--------+-----------------------+-------------------------------------------------+
   |  USE_OCC                    | OFF    | Compiles with OpenCascade support to enable                             |
   |                             |        | various geometry types (e.g. \*.iges, \*.step).                         |
   |                             |        | Currently not tested on MacOS.                                          |
   +                             +        +-------------------+-----------------------------------------------------+
   |                             |        | Ubunutu           | The required library can be installed with          |
   |                             |        |                   |                                                     |
   |                             |        |                   | .. code:: bash                                      |
   |                             |        |                   |                                                     |
   |                             |        |                   |   sudo apt-add-repository universe                  |
   |                             |        |                   |   sudo apt-get update                               |
   |                             |        |                   |   sudo apt-get install libocct-data-exchange-dev    |
   |                             |        |                   |                        libocct-draw-dev occt-misc   |
   +                             +        +-------------------+-----------------------------------------------------+
   |                             |        | Windows           | A precompiled version is downloaded and used.       |
   +-----------------------------+--------+-------------------+-----------------------------------------------------+
   |  USE_GUI                    | ON     | Build Netgen with GUI                                                   |
   +-----------------------------+--------+-------------------+-----------------------------------------------------+
   |  USE_PYTHON                 | ON     | Enable Python bindings                                                  |
   +-----------------------------+--------+-------------------+-----------------------------------------------------+
   |  USE_MPI                    | OFF    | Enable MPI parallelization support                                      |
   +-----------------------------+--------+-------------------+-----------------------------------------------------+
   |  USE_JPEG                   | OFF    | Build with JPEG support to make screenshots from the GUI                |
   +-----------------------------+--------+-------------------+-----------------------------------------------------+
   |  USE_MPEG                   | OFF    | Enable video recording using libavcodec                                 |
   +-----------------------------+--------+-------------------+-----------------------------------------------------+
   |  USE_LAPACK                 | ON     | Link with LAPACK/BLAS libraries                                         |
   +-----------------------------+--------+-------------------+-----------------------------------------------------+
   |  USE_MUMPS                  | OFF    | Enable sparse direct solver MUMPS                                       |
   +-----------------------------+--------+-------------------+-----------------------------------------------------+
   |  USE_PARDISO                | OFF    | Enable sparse direct solver PARDISO                                     |
   +-----------------------------+--------+-------------------+-----------------------------------------------------+
   |  USE_UMFPACK                | OFF    | Enable sparse direct solver SuiteSparse/UMFPACK                         |
   +-----------------------------+--------+-------------------+-----------------------------------------------------+
   |  USE_VTUNE                  | OFF    | Enable Intel VTune pause/resume numproc                                 |
   +-----------------------------+--------+-------------------+-----------------------------------------------------+
   |  USE_NUMA                   | OFF    | Compile with NUMA-aware code on multisocket machines (requires libnuma) |
   +-----------------------------+--------+-------------------+-----------------------------------------------------+
