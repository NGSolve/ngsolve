CMake options for building Netgen/NGSolve
*****************************************

+----------------+------------+-------------------+-----------------------------------------------------+
| CMake flag     | Default    |                   | Description                                         |
+================+============+===================+=====================================================+
| -DINSTALL_DIR  | unused     | Directs the installation to this folder                                 |
|                |            |                                                                         |
+----------------+------------+-------------------+-----------------------------------------------------+
| -DUSE_MKL      | OFF        |                   | Some Intel MKL function are used                    |
|                |            |                   | (e.g. Lapack, Pardiso,...)                          |
+----------------+------------+-------------------+-----------------------------------------------------+
| -DMKL_ROOT     | unused     |                   | defines the path to your	                        |
|                |            |                   | MKL installation                                    |
+----------------+------------+-------------------+-----------------------------------------------------+
| -DUSE_OCC      | OFF        | Requires the installation of OpenCascade and enables the functionality  |
|                |            | to use various geometry types (e.g. \*.iges, \*.step).                  |
+                +            +-------------------+-----------------------------------------------------+
|                |            | Ubunutu/Linux     | The required library can be installed with          |
|                |            |                   |                                                     |
|                |            |                   | .. code:: bash                                      |
|                |            |                   |                                                     |
|                |            |                   |   sudo apt-add-repository universe                  |
|                |            |                   |   sudo apt-get update                               |
|                |            |                   |   sudo apt-get install liboce-ocaf-dev              |
+                +            +-------------------+-----------------------------------------------------+
|                |            | Windows           | A precompiled version is downloaded and used.       |
+                +            +-------------------+-----------------------------------------------------+
|                |            | Mac OS X          |                                                     |
+----------------+------------+-------------------+-----------------------------------------------------+
| -DUSE_CACHE    | OFF        |                   |                                                     |
|                |            |                   |                                                     |
+----------------+------------+-------------------+-----------------------------------------------------+
