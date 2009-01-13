AC_INIT([netgen], [4.9.2], [], [])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_ARG_ENABLE([occ], 
        [  --enable-occ            Compile with OpenCascade geometry kernel],
        [occon=true, occfl=-DOCC],
        [occon=false] )


AM_CONFIG_HEADER(config.h)

AC_LANG_CPLUSPLUS
AC_PROG_CXX
AM_PROG_LIBTOOL

AM_CONDITIONAL([OCC], [test x$occon = xtrue])
AC_SUBST([OCCFLAGS], [$occfl])

AC_OUTPUT(Makefile libsrc/Makefile libsrc/csg/Makefile libsrc/general/Makefile \
	libsrc/geom2d/Makefile libsrc/gprim/Makefile libsrc/include/Makefile libsrc/interface/Makefile \
	libsrc/linalg/Makefile libsrc/meshing/Makefile libsrc/occ/Makefile libsrc/opti/Makefile \
	libsrc/parallel/Makefile libsrc/stlgeom/Makefile libsrc/visualization/Makefile ng/Makefile \
	tutorials/Makefile)




