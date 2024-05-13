#ifndef FILE_NGLA
#define FILE_NGLA

#include <bla.hpp>
#include <core/mpi_wrapper.hpp>

/*
namespace ngcomp
{
  class Preconditioner;
  class LocalPreconditioner;
}
*/

/** namespace for linear algebra.
 */

namespace ngla
{
  using namespace std;
  using namespace ngstd;
  using namespace ngbla;
}


#include "paralleldofs.hpp"
#include "basevector.hpp"
#include "vvector.hpp"
#include "multivector.hpp"
#include "basematrix.hpp"
#include "sparsematrix.hpp"
#include "sparsematrix_dyn.hpp"
#include "order.hpp"
#include "sparsecholesky.hpp"
#include "pardisoinverse.hpp"
// include these only from c++-files
// #include "umfpackinverse.hpp"
// #include "superluinverse.hpp"
// #include "mumpsinverse.hpp"
#include "jacobi.hpp"
#include "blockjacobi.hpp"
#include "commutingAMG.hpp"
#include "special_matrix.hpp"
#include "elementbyelement.hpp"
#include "cg.hpp"
#include "chebyshev.hpp"
#include "eigen.hpp"
#include "arnoldi.hpp"

#endif
