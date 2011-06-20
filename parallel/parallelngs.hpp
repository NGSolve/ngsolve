#ifndef FILE_PARALLELNGS
#define FILE_PARALLELNGS



#ifdef VTRACE
#include "vt_user.h"
#else
#define VT_USER_START(n)
#define VT_USER_END(n)
#define VT_TRACER(n)
#endif


#ifdef PARALLEL

#include <comp.hpp>

/*
#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
*/

#include <mpi.h>

namespace netgen {
  extern MPI_Group MPI_HIGHORDER_WORLD;
  extern MPI_Comm MPI_HIGHORDER_COMM;
}
using netgen::MPI_HIGHORDER_COMM;

#include "mpiwrapper.hpp"
#include "paralleldofs.hpp"

#include "parallelvector.hpp"
#include "parallel_matrices.hpp"

extern void Parallel_Exit ();

#endif

#endif




