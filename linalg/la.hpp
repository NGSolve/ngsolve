#ifndef FILE_NGLA
#define FILE_NGLA

#include <bla.hpp>


/*
#ifdef LAPACK
#include "../lapack/LapackGEP.hpp"
#endif
*/


#ifdef PARALLEL
#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
#include <mpi.h>
#endif





namespace ngparallel
{
  class ParallelDofs;
}

namespace ngcomp
{
  class Preconditioner;
  class LocalPreconditioner;
}
/** namespace for linear algebra.
 */

namespace ngla
{
  using namespace std;
  using namespace ngstd;
  using namespace ngbla;

  using ngstd::BitArray;
  using ngstd::FlatArray;
  using ngstd::Array;
  using ngstd::Table;
  using ngstd::INT;
}

#include "basevector.hpp"
#include "vvector.hpp"
#include "basematrix.hpp"
#include "sparsematrix.hpp"
#include "order.hpp"
#include "sparsecholesky.hpp"
#include "pardisoinverse.hpp"
#include "superluinverse.hpp"
#include "mumpsinverse.hpp"
#include "jacobi.hpp"
#include "blockjacobi.hpp"
#include "commutingAMG.hpp"
#include "special_matrix.hpp"
#include "cg.hpp"
#include "chebyshev.hpp"
#include "eigen.hpp"




#endif
