#ifndef FILE_NGLA
#define FILE_NGLA

#include <bla.hpp>

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

  using ngstd::BitArray;
  using ngstd::FlatArray;
  using ngstd::Array;
  using ngstd::Table;
  using ngstd::INT;
}


#include "paralleldofs.hpp"
#include "basevector.hpp"
#include "vvector.hpp"
#include "basematrix.hpp"
#include "sparsematrix.hpp"
#include "order.hpp"
#include "sparsecholesky.hpp"
// include these only from c++-files
// #include "pardisoinverse.hpp"
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

#include "cuda_linalg.hpp"


#endif
