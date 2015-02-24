/*

Discretization of Stokes equation

-nu Laplace u + grad p = f
div u = 0

 */


#include <solve.hpp>

using namespace ngsolve;




/*
  Build one finite element space for the three fields  u_x, u_y, p.
  Taylor-Hood element: P2 for the velocities u, and P1 for the pressure
*/

class FESpaceStokes : public CompoundFESpace
{
public:

  FESpaceStokes (const MeshAccess & ama, const Flags & flags)
    : CompoundFESpace (ama, flags)
  {
    int order = int (flags.GetNumFlag ("order", 2));
    // if (order < 2)
    // throw Exception ("Taylor-Hood elements need order 2 or higher");

    Flags uflags, pflags;
    uflags.SetFlag ("order", order+1);
    uflags.SetFlag ("orderinner", order+1);
    AddSpace (new H1HighOrderFESpace (ma, uflags));
    AddSpace (new H1HighOrderFESpace (ma, uflags));
    
    pflags.SetFlag ("order", order);
    AddSpace (new H1HighOrderFESpace (ma, pflags));

    // pflags.SetFlag ("order", order-1);
    // AddSpace (new L2HighOrderFESpace (ma, pflags));
  }
  
  virtual string GetClassName () const { return "Demo-StokesFESpace"; }
};



/*
  The weak form is written as

  \int (Bv)^T D Bu dx, 

  where B is a differential operator, and D is a matrix
*/

class DiffOpStokes : public DiffOp<DiffOpStokes>
{
  // 5 components:
  // du1/dx, du1/dy, du2/dx, du2/dy, p

public:
  enum { DIM = 1 };          // just one copy of the finite element spaces
  enum { DIM_SPACE = 2 };    // 2D domain 
  enum { DIM_ELEMENT = 2 };  // 2D elements (in contrast to 1D boundary elements)
  enum { DIM_DMAT = 5 };     // D-matrix is 5x5
  enum { DIFFORDER = 0 };    // minimal differential order (to determine integration order)
  

  /*
    Computes the 5 x N matrix consisting of B_j(phi_i).
    The phi_i, with j = 0,...N-1 are the shape functions,
    and j=0...4 are the 5 components of the diff-op:

 
    The arguments are generic.
    bfel  must be a proper finite element, i.e., a compound element consiting of 3 scalar parts.
    mip  is the mapped integration point, containing the Jacobi matrix of the mapping
    mat  is a matrix. Can be a general dense matrix, or a matrix of fixed height 5 (thus generic)
  */

  template <typename FEL, typename MIP, typename MAT>
  static void GenerateMatrix (const FEL & bfel, const MIP & mip,
			      MAT & mat, LocalHeap & lh)
  {
    // must get the right elements, otherwise an exception is thrown.

    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (bfel);

    // a scalar H1 element
    const ScalarFiniteElement<2> & fel_u = 
      dynamic_cast<const ScalarFiniteElement<2>&> (cfel[0]);
    const ScalarFiniteElement<2> & fel_p = 
      dynamic_cast<const ScalarFiniteElement<2>&> (cfel[2]);
    
    int nd_u = fel_u.GetNDof();
    int nd_p = fel_p.GetNDof();
    
    // transformation of derivatives from reference element to general element:
    FlatMatrixFixWidth<2> gradu(nd_u, lh);
    fel_u.CalcMappedDShape (mip, gradu);

    // the shape functions of the pressure
    FlatVector<> vecp(nd_p, lh);
    fel_p.CalcShape (mip.IP(), vecp);

    mat = 0;

    // the first nd_u shape functions belong to u_x, the next nd_u belong to u_y:
    mat.Rows(0,2).Cols(cfel.GetRange(0)) = Trans (gradu);
    mat.Rows(2,4).Cols(cfel.GetRange(1)) = Trans (gradu);

    // ... and finally nd_p shape functions for the pressure:
    mat.Row(4).Range(cfel.GetRange(2)) = vecp;
  }
};



/*
  The 5x5 coefficient matrix:

  nu 0  0  0  1
  0  nu 0  0  0
  0  0  nu 0  0
  0  0  0  nu 1
  1  0  0  1  0

*/

class StokesDMat : public DMatOp<StokesDMat,5>
{
  CoefficientFunction * nu;
public:

  enum { DIM_DMAT = 5 };

  StokesDMat (CoefficientFunction * anu) : nu(anu) { ; }
  
  template <typename FEL, typename MIP, typename MAT>
  void GenerateMatrix (const FEL & fel, const MIP & mip,
		       MAT & mat, LocalHeap & lh) const
  {
    mat = 0;

    // (nu  grud u, grad v)
    double val = nu -> Evaluate (mip);
    for (int i = 0; i < 4; i++)
      mat(i, i) = val;

    // (du1/dx+du2/dy, p)
    mat(4,0) = mat(0,4) = 1;
    mat(4,3) = mat(3,4) = 1;
  }  
};


/*
  Combine differential operator and D-matrix to an integrator
*/
class StokesIntegrator 
  : public T_BDBIntegrator<DiffOpStokes, StokesDMat, FiniteElement>
{
public:

  StokesIntegrator (Array<CoefficientFunction*> & coeffs)
    :  T_BDBIntegrator<DiffOpStokes, StokesDMat, FiniteElement>
  (StokesDMat (coeffs[0]))
  { ; }

  virtual string Name () const { return "Stokes"; }
};






/*
  Evaluate (u_x, u_y) 
 */
class DiffOpIdU : public DiffOp<DiffOpIdU>
{
  // 2 components:
  // u1 u2

public:
  enum { DIM = 1 };
  enum { DIM_SPACE = 2 };
  enum { DIM_ELEMENT = 2 };
  enum { DIM_DMAT = 2 };
  enum { DIFFORDER = 0 };
  
  template <typename FEL, typename MIP, typename MAT>
  static void GenerateMatrix (const FEL & bfel, const MIP & mip,
			      MAT & mat, LocalHeap & lh)
  {
    const CompoundFiniteElement & cfel = 
      dynamic_cast<const CompoundFiniteElement&> (bfel);
    const ScalarFiniteElement<2> & fel_u = 
      dynamic_cast<const ScalarFiniteElement<2>&> (cfel[0]);
    
    int nd_u = fel_u.GetNDof();

    FlatVector<> vecu(nd_u, lh);
    fel_u.CalcShape (mip.IP(), vecu);

    mat = 0;
    mat.Row(0).Range(cfel.GetRange(0)) = vecu;
    mat.Row(1).Range(cfel.GetRange(1)) = vecu;
  }
};


class StokesUIntegrator 
  : public T_BDBIntegrator<DiffOpIdU, DiagDMat<2>, FiniteElement>
{
public:
  ///
  StokesUIntegrator (Array<CoefficientFunction*> & /* coeffs */)
    :  T_BDBIntegrator<DiffOpIdU, DiagDMat<2>, FiniteElement>
  (DiagDMat<2> (new ConstantCoefficientFunction(1)))
  { ; }

  ///
  virtual string Name () const { return "Stokes IdU"; }
};



static RegisterFESpace<FESpaceStokes> initfes ("stokes");
static RegisterBilinearFormIntegrator<StokesIntegrator> initstokes ("stokes", 2, 1);
static RegisterBilinearFormIntegrator<StokesUIntegrator> initstokesu ("stokesu", 2, 0);
