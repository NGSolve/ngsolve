#ifndef FILE_NGLA
#define FILE_NGLA

#include <bla.hpp>

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
#include "jacobi.hpp"
#include "blockjacobi.hpp"
#include "commutingAMG.hpp"
#include "diagonalmatrix.hpp"
#include "special_matrix.hpp"
#include "elementbyelement.hpp"
#include "cg.hpp"
#include "chebyshev.hpp"
#include "eigen.hpp"
#include "arnoldi.hpp"

#endif
