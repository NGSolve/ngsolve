#include <comp.hpp>

namespace ngcomp
{

  class NumberFiniteElement : public FiniteElement
  {
  public:
    NumberFiniteElement ()
      : FiniteElement(1, 0) { ; }
    HD virtual ELEMENT_TYPE ElementType() const { return ET_POINT; }
  };


  class NumberDiffOp : public DiffOp<NumberDiffOp>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 0 };
    enum { DIM_ELEMENT = 0 };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 0 };
  
    template <typename MIP, typename MAT>
    static void GenerateMatrix (const FiniteElement & fel, const MIP & mip,
				MAT && mat, LocalHeap & lh)
    {
      mat(0,0) = 1;
    }
  };





  class NumberFESpace : public FESpace
  {
  public:
    NumberFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool checkflags=false)
      : FESpace (ama, flags)
    { 
      evaluator = make_shared<T_DifferentialOperator<NumberDiffOp>>();
    }

    virtual int GetNDof() const { return 1; }

    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const
    {
      return *new (lh) NumberFiniteElement();
    }

    virtual void GetDofNrs (int elnr, Array<int> & dnums) const
    {
      dnums.SetSize(1);
      dnums[0] = 0;
    }

    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const
    {
      dnums.SetSize(1);
      dnums[0] = 0;
    }

  };


  static RegisterFESpace<NumberFESpace> init ("number");

}
