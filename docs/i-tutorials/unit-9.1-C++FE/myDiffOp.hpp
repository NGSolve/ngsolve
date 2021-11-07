#ifndef FILE_MYDIFFOP_HPP
#define FILE_MYDIFFOP_HPP

#include <fem.hpp>

#include "myElement.hpp"

namespace ngfem
{
  /* 
     DiffOps provide the link between function evaluation, and finite elements.
     Different DiffOps evaluate either shape functions, or derivatives.
     Typically, a FiniteElementSpace creates a DiffOp to for function evaluation, 
     and for the canonical derivative, as well as DiffOps for evaluation at the
     boundary. These DiffOps are used when evaluating GridFunctions, and setting 
     up element-matrices from trial- and test-functions.
     DiffOps use static polymorphism, aka Curiously Recurring Template Pattern (CRTP).
   */

  // Our implementation of the identity operator, in two space dimensions.
  class MyDiffOpId : public DiffOp<MyDiffOpId>
  {
  public:
    // some constants for the diffop:
    
    static constexpr int DIM = 1;       // dimension of the input
    static constexpr int DIM_SPACE = 2; // dimension of the space
    static constexpr int DIM_ELEMENT = 2; // spatial dimension of the element
    static constexpr int DIM_DMAT = 1;  // dimension of the output
    static constexpr int DIFFORDER = 0; // order of differentiation

    static bool SupportsVB (VorB checkvb) { return true; } // can do VOL and BND terms

    // fill evaluation matrix of dimension 1 times fel.ndof 
    // the input mip is a mapped integration point, which has also
    // access to the integration point on the reference element
    template<typename MIP, typename MAT>
    static void GenerateMatrix(const FiniteElement & fel, const MIP & mip,
                               MAT & mat, LocalHeap & lh)
    {
      dynamic_cast<const MyBaseElement&> (fel).CalcShape(mip.IP(), mat.Row(0));
    }

    // can overload more functionality for performance optimization,
    // like evaluation in the whole integration rule
  };
    
   

    
  // Gradient DiffOp: implements the chain rule for mapping gradients
  // from the reference element to the pysical element
  class MyDiffOpGradient : public DiffOp<MyDiffOpGradient>
  {
  public:    
    static constexpr int DIM = 1;       // dimension of the input
    static constexpr int DIM_SPACE = 2; // dimension of the space
    static constexpr int DIM_ELEMENT = 2; // dimension of the element
    static constexpr int DIM_DMAT = 2;  // dimension of the output
    static constexpr int DIFFORDER = 1; // order of differentiation
    
    // so that you can call grad(u)
    static string Name() { return "grad"; }

    static const MyBaseElement& Cast(const FiniteElement & fel)
    { return static_cast<const MyBaseElement&> (fel); }

    // computes the 2 times ndof evaluation matrix
    template<typename MIP, typename MAT>
    static void GenerateMatrix(const FiniteElement & fel, const MIP & mip,
                               MAT & mat, LocalHeap & lh)
    {
      // the gradient on the physical element follows by the chain-rule
      // from the gradient on the reference element
      // F is the Jacobian matrix of the mapping
      // (phi_i)' = (phiref_i)' F^{-1}

      // The LocalHeap provides cheap dynamic memory allocation
      // for temporary matrices
      HeapReset hr(lh);
      FlatMatrix<double> dshape(fel.GetNDof(), 2, lh);
      
      // gradients of basis functions on the reference element
      Cast(fel).CalcDShape(mip.IP(), dshape);
      mat = Trans(dshape * mip.GetJacobianInverse());
    }
  };
}
#endif // FILE_MYDIFFOP_HPP
