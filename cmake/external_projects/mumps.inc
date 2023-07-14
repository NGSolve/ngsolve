ORDERINGSF  = -Dmetis -Dparmetis
ORDERINGSC  = $(ORDERINGSF)

LORDERINGS = -L@PARMETIS_DIR@/lib -lparmetis
IORDERINGSF = 
IORDERINGSC = -I@PARMETIS_DIR@/include

PLAT    =
LIBEXT  = .a
OUTC    = -o 
OUTF    = -o
RM = /bin/rm -f
CC = @MPI_C_COMPILER@
FC = @MPI_Fortran_COMPILER@
FL = @MPI_Fortran_COMPILER@
AR = ar vr 
RANLIB = ranlib
SCALAP  = @SCALAPACK_LIBRARY@ @BLACS_LIBRARY@ @LAPACK_LIBRARY@

INCPAR = -I../include -I/usr/include

LIBPAR = $(SCALAP) @MPI_C_LIBRARIES@ @MPI_Fortran_LIBRARIES@

# INCSEQ = -I$(topdir)/libseq
# LIBSEQ  =  -L$(topdir)/libseq -lmpiseq

LIBBLAS = -lblas
LIBOTHERS = -lpthread

#Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
CDEFS   = -DAdd_

#Begin Optimized options
# uncomment -fopenmp in lines below to benefit from OpenMP
OPTF    = -fallow-argument-mismatch -O  -fPIC -DALLOW_NON_INIT # -fopenmp
OPTL    = -O -fPIC # -fopenmp
OPTC    = -O -fPIC # -fopenmp
#End Optimized options
INCS = $(INCPAR)
LIBS = $(LIBPAR)
LIBSEQNEEDED =
