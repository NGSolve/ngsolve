#ifndef FILE_PARALLELNGS
#define FILE_PARALLELNGS


#include <ngstd.hpp>
#include <la.hpp>

namespace ngparallel
{
  using namespace ngstd;
  using namespace ngla;
}

#ifdef PARALLEL

/*
#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
*/

#define OMPI_SKIP_MPICXX
#include <mpi.h>

#endif




#include "mpiwrapper.hpp"
#include "paralleldofs.hpp"

#include "parallelvector.hpp"
#include "parallel_matrices.hpp"


// extern void Parallel_Exit ();


#endif




