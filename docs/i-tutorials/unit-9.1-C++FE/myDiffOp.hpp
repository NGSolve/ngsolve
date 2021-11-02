#ifndef FILE_MYDIFFOP_HPP
#define FILE_MYDIFFOP_HPP

#include <fem.hpp>

#include "myElement.hpp"

namespace ngfem
{
  /* 
     DiffOps are the static polymorphism classes implementing
     Operators for GridFunctions and Test- and TrialFunctions.
     They only have static member functions and are never created
     but only used as a template argument for DifferentialOperators.
     They use CRTP for "static virtual" methods.
     DifferentialOperator is the runtime polymorphism class wrapping
     them. Use T_DifferentialOperator<DiffOp> as a coupling
     mechanism.
   */

  // Our implementation of the identity
  class MyDiffOpId : public DiffOp<MyDiffOpId>
  {
  public:
    // ******************** Necessary stuff *************************
    
    static constexpr int DIM = 1; // dimension of the input
    static constexpr int DIM_SPACE = 2; // dimension of the space
    static constexpr int DIM_ELEMENT = 2; // dimension of the element
    static constexpr int DIM_DMAT = 1; // dimension of the output
    static constexpr int DIFFORDER = 0; // order of differentiation

    template<typename MIP, typename MAT>
    static void GenerateMatrix(const FiniteElement& fel, const MIP& mip,
                               MAT& mat, LocalHeap& lh)
    {
      HeapReset hr(lh);
      Cast(fel).CalcShape(mip.IP(), mat.Row(0));
    }

    // ******************** Helper function ***************************
    static const MyBaseElement& Cast(const FiniteElement& fel)
    { return static_cast<const MyBaseElement&> (fel); }


    // ******************** Performance improvements *******************
    // ...
  };

    // Our implementation of the gradient of a function
  class MyDiffOpGradient : public DiffOp<MyDiffOpGradient>
  {
  public:
    // ******************** Necessary stuff *************************
    
    static constexpr int DIM = 1; // dimension of the input
    static constexpr int DIM_SPACE = 2; // dimension of the space
    static constexpr int DIM_ELEMENT = 2; // dimension of the element
    static constexpr int DIM_DMAT = 2; // dimension of the output
    static constexpr int DIFFORDER = 1; // order of differentiation
    
    // so that you can call grad(u)
    static string Name() { return "grad"; }

    template<typename MIP, typename MAT>
    static void GenerateMatrix(const FiniteElement& fel, const MIP& mip,
                               MAT& mat, LocalHeap& lh)
    {
      // We need to compute the gradient on the reference element and
      // then transform it to the using the jacobian of the mapped
      // integration point
      HeapReset hr(lh);
      // create a temporary matrix on the local heap
      FlatMatrixFixWidth<2> dshape(fel.GetNDof(), lh);
      Cast(fel).CalcDShape(mip.IP(), dshape);
      mat = Trans(dshape * mip.GetJacobianInverse());
    }

    // ******************** Helper function ***************************
    static const MyBaseElement& Cast(const FiniteElement& fel)
    { return static_cast<const MyBaseElement&> (fel); }


    // ******************** Performance improvements *******************
    // ...
  };
}
#endif // FILE_MYDIFFOP_HPP
