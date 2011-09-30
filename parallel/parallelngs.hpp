#ifndef FILE_PARALLELNGS
#define FILE_PARALLELNGS




#ifdef VTRACE
#include "vt_user.h"
#else
#define VT_USER_START(n)
#define VT_USER_END(n)
#define VT_TRACER(n)
#endif

#include <ngstd.hpp>
#include <la.hpp>

using namespace ngstd;
using namespace ngla;
// #include <comp.hpp>


#ifdef PARALLEL

/*
#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
*/

#include <mpi.h>

#endif



#include "mpiwrapper.hpp"
#include "paralleldofs.hpp"

#include "parallelvector.hpp"
#include "parallel_matrices.hpp"

extern void Parallel_Exit ();





#endif




