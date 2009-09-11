AC_INIT([ngsolve],[4.9.11],[])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AC_PROG_CXX

AC_PREFIX_DEFAULT(["/opt/netgen"])
AC_DISABLE_STATIC
AC_PROG_LIBTOOL

# ACX_LAPACK


# Tcl/Tk configuration:

TEA_INIT([3.6])

TEA_PATH_TCLCONFIG
TEA_LOAD_TCLCONFIG 

# TEA_PATH_TKCONFIG
# TEA_LOAD_TKCONFIG

TEA_PUBLIC_TCL_HEADERS




AC_ARG_WITH([netgen], 
        [  --with-netgen=dir       use Netgen installed in directory dir],
        [netgendir=$withval],
        [netgendir="/opt/netgen"]
        )


AC_SUBST([LAPACK_INCLUDES], ["-DLAPACK"])
AC_ARG_WITH([lapack],
        [  --with-lapack=lapacklib   linker flags for linking dense linear algebra library LAPACK, default is '-llapack -lblas'],          
        [AC_SUBST([LAPACK_LIBS], [$withval])],
        [AC_SUBST([LAPACK_LIBS], "-llapack -lblas")],
        )

AC_ARG_ENABLE([pardiso],
        [  --enable-pardiso           enable sparse direct solver pardiso],          
        [AC_SUBST([PARDISO_INCLUDES], ["-DUSE_PARDISO"])]
        [AC_SUBST([PARDISO_LIBS], ["-lpardiso_GNU_IA32 -llapack"])]
#        [AC_SUBST([PARDISO_LIBS], ["-lpardiso_GNU42_EM64T_INT64_P -lgfortran"])]
        )


AC_ARG_WITH([pardiso],
        [  --with-pardiso=libpardiso.so  enable sparse direct solver pardiso],          
        [AC_SUBST([PARDISO_INCLUDES], ["-DUSE_PARDISO"])]
        [AC_SUBST([PARDISO_LIBS], ["$withval -llapack"])]
#        [AC_SUBST([PARDISO_LIBS], ["-lpardiso_GNU42_EM64T_INT64_P -lgfortran"])]
        )

AC_ARG_ENABLE([MKLpardiso],
        [  --enable-MKLpardiso        enable sparse direct solver pardiso from MKL],          
        [AC_SUBST([PARDISO_INCLUDES], ["-DUSE_PARDISO -DUSE_MKL"])]
        )



AC_ARG_ENABLE([mumps],
        [  --enable-mumps             enable sparse direct solver MUMPS],          
        [AC_SUBST([MUMPS_INCLUDES], ["-DUSE_MUMPS"])]
        [AC_SUBST([MUMPS_LIBS], ["-ldmumps -lzmumps -lmumps_common -lpord -lmpiseq -lgfortran"])]
        )



# AC_ARG_WITH([pardiso],
#       [  --with-pardiso=libpardiso.so  enable sparse direct solver pardiso],          
#       [AC_SUBST([PARDISO_INCLUDES], ["-DUSE_PARDISO"])]
#       [AC_SUBST([PARDISO_LIBS], [$withval])]
#       )

#	[AC_SUBST([LAPACK_LIB], ["-L/opt/intel/mkl/10.0.2.018/lib/em64t -lmkl_core -lmkl_def -lmkl_sequential -lmkl_gf_lp64"])]
#        [AC_SUBST([LAPACK_LIB], ["-llapack"])]
#	[AC_SUBST([LAPACK_LIB], ["-L/opt/intel/mkl/10.0.2.018/lib/em64t -lmkl_sequential"])]


AC_SUBST([NETGEN_INCLUDES], ["-DNETGEN_ELTRANS -I$netgendir/include"])




AC_CONFIG_HEADERS(config.h) AC_CONFIG_FILES(Makefile include/Makefile
ngstd/Makefile basiclinalg/Makefile linalg/Makefile fem/Makefile
comp/Makefile solve/Makefile multigrid/Makefile parallel/Makefile 
programming_demos/Makefile pde_tutorial/Makefile windows/Makefile 
my_little_ngsolve/Makefile)

AC_OUTPUT
