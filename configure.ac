AC_INIT([ngsolve],[4.9.14-dev],[])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_CONFIG_MACRO_DIR([m4])

# AC_PROG_CC
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
        [AC_HELP_STRING([--with-lapack=lapacklib][linker flags for linking dense linear algebra library LAPACK, default is '-llapack -lblas'])],          
        [AC_SUBST([LAPACK_LIBS], [$withval])],
        [AC_SUBST([LAPACK_LIBS], "-llapack -lblas")],
        )

AC_ARG_ENABLE([pardiso],
        [AC_HELP_STRING([--enable-pardiso],[enable sparse direct solver pardiso])], 
        [AC_SUBST([PARDISO_INCLUDES], ["-DUSE_PARDISO"])]
        [AC_SUBST([PARDISO_LIBS], ["-lpardiso_GNU_IA32 -llapack"])]
#        [AC_SUBST([PARDISO_LIBS], ["-lpardiso_GNU42_EM64T_INT64_P -lgfortran"])]
        )

AC_ARG_ENABLE([pardiso4],
        [AC_HELP_STRING([--enable-pardiso4],[enable sparse direct solver pardiso])],
        [AC_SUBST([PARDISO_INCLUDES], ["$PARDISO_INCLUDES -DUSE_PARDISO400"])]
        )

AC_ARG_WITH([pardiso],
        [AC_HELP_STRING([--with-pardiso=libpardiso.so],[enable sparse direct solver pardiso])],
        [AC_SUBST([PARDISO_INCLUDES], ["$PARDISO_INCLUDES -DUSE_PARDISO"])]
        [AC_SUBST([PARDISO_LIBS], ["$withval -llapack -lblas -lgomp -gfortran"])]
#        [AC_SUBST([PARDISO_LIBS], ["-lpardiso_GNU42_EM64T_INT64_P -lgfortran"])]
        )

AC_ARG_ENABLE([MKLpardiso],
        [  --enable-MKLpardiso        enable sparse direct solver pardiso from MKL],          
        [AC_SUBST([PARDISO_INCLUDES], ["-DUSE_PARDISO -DUSE_MKL"])]
        )



AC_ARG_ENABLE([mumps],
        [AC_HELP_STRING([--enable-mumps],[enable sparse direct solver MUMPS])],
        [AC_SUBST([MUMPS_INCLUDES], ["-DUSE_MUMPS"])]
        [AC_SUBST([MUMPS_LIBS], ["-L/opt/mumps/lib -ldmumps -lzmumps -lmumps_common -lpord -lgfortran -lmpi_f77 -lscalapack"])]
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

AC_ARG_ENABLE([parallel],
        [AC_HELP_STRING([--enable-parallel],[enable mpi parallelization])],
        [AC_SUBST([MPI_INCLUDES], "-DPARALLEL -I/usr/local/include/vampirtrace -DnoVTRACE")
        ]
        )       

#        [AC_SUBST([MPI_INCLUDES], "-I/usr/lib64/mpi/gcc/openmpi/include -DPARALLEL")
#        [AC_SUBST([MPI_INCLUDES], "-I/opt/mpich/include -I/usr/include/mpich2 -I/usr/include/metis -DPARALLEL")
#        AC_SUBST([MPI_LIBS], "-lmetis -L/opt/mpich/ch-p4/lib -lmpich -lvt")

 


AC_CONFIG_HEADERS(config.h) AC_CONFIG_FILES(Makefile include/Makefile
ngstd/Makefile basiclinalg/Makefile linalg/Makefile fem/Makefile
comp/Makefile solve/Makefile multigrid/Makefile parallel/Makefile 
programming_demos/Makefile pde_tutorial/Makefile windows/Makefile)

AC_OUTPUT
