#ifndef FILE_FMMOPERATOR
#define FILE_FMMOPERATOR
#include "fmminterface.hpp"
#include <ngstd.hpp>
#include <bla.hpp>
#include <variant>
#include "../linalg/basematrix.hpp"
#include "../linalg/vvector.hpp"

namespace ngsbem
{
  using namespace ngla;

  struct FMMOperatorInfo
  {
    string kernel_name;                          // Kernel name of the inspected FMM operator.
    size_t source_size = 0;                      // Number of source quadrature points used by the FMM operator.
    size_t target_size = 0;                      // Number of target quadrature points used by the FMM operator.
    size_t source_dofs = 0;                      // Number of trial-space degrees of freedom.
    size_t target_dofs = 0;                      // Number of test-space degrees of freedom.
    std::variant<double, Complex> kappa = double(0); // Wave number used by the kernel.
    size_t total_num_nodes = 0;                  // Combined number of source and target tree nodes.
    size_t total_num_leaves = 0;                 // Combined number of source and target tree leaves.
    size_t total_multipole_coefficients = 0;     // Combined number of stored source and target expansion coefficients.
    size_t total_memory_bytes = 0;               // Combined memory used by stored expansion coefficients.
    size_t num_s2r = 0;                          // Number of recorded source-to-regular box translations.
    size_t num_direct_evaluations = 0;           // Estimated quadrature-point interactions handled by direct fallback.
    double direct_fallback_fraction = 0;         // Fraction of dense quadrature-point pairs handled by direct fallback.
    size_t nearfield_nze = 0;                    // Nonzero entries in the sparse Sauter-Schwab nearfield correction.
    double nearfield_fraction = 0;               // nearfield_nze divided by target_dofs*source_dofs.
    double source_bbox_radius = 0;               // Radius of the source root box.
    double target_bbox_radius = 0;               // Radius of the target root box.
    FMMTreeStats source_tree;                    // Source tree statistics.
    FMMTreeStats target_tree;                    // Target tree statistics.
    FMM_Parameters parameters;                   // FMM parameters used to build the trees.
  };

  inline std::tuple<Vec<3>, double> GetCenterAndRadius(const Array<Vec<3>>& xpts)
  {
    Vec<3> xmax(-1e99, -1e99, -1e99);
    Vec<3> xmin(1e99, 1e99, 1e99);
    for (auto xi : xpts)
      {
        for (int j = 0; j < 3; j++)
          {
            xmin(j) = min(xmin(j), xi(j));
            xmax(j) = max(xmax(j), xi(j));
          }
      }

    Vec<3> cx = 0.5*(xmin+xmax);
    double rx = MaxNorm(xmax-xmin);

    return { cx, rx };
  }
  
  
  template <typename TSCAL>
  class Base_FMM_Operator : public BaseMatrix
  {
  protected:
    Array<Vec<3>> xpts, ypts, xnv, ynv;
    Vec<3> cx, cy;
    double rx, ry;
    IVec<2> kernelshape;
    FMM_Parameters fmm_params;
  public:
    Base_FMM_Operator(Array<Vec<3>> _xpts, Array<Vec<3>> _ypts,
                      Array<Vec<3>> _xnv, Array<Vec<3>> _ynv, IVec<2> _kernelshape,
                      const FMM_Parameters & _params)
      : xpts(std::move(_xpts)), ypts(std::move(_ypts)),
        xnv(std::move(_xnv)), ynv(std::move(_ynv)), kernelshape(_kernelshape), fmm_params(_params)
    {
      std::tie(cx, rx) = GetCenterAndRadius(xpts);      
      std::tie(cy, ry) = GetCenterAndRadius(ypts);
    }

    int VHeight() const override { return  ypts.Size()*kernelshape[0]; }
    int VWidth() const override { return  xpts.Size()*kernelshape[1]; }
      
    VecFormat RowFormat () const override { return VVectorFormat<TSCAL> (xpts.Size() * kernelshape[1]); }
    VecFormat ColFormat () const override { return VVectorFormat<TSCAL> (ypts.Size() * kernelshape[0]); }

    virtual FMMOperatorInfo GetFMMInfo () const = 0;
  };


  
  template <typename TSCAL>
  class NGS_DLL_HEADER FMM_Operator : public Base_FMM_Operator<TSCAL>
  {
    string kernel_name;
    shared_ptr<const BaseFMMInterface> source, target;
    Timer<> tapply, teval, ttrans;
    typedef Base_FMM_Operator<TSCAL> BASE;
    using BASE::xpts, BASE::ypts, BASE::xnv, BASE::ynv, BASE::cx, BASE::cy, BASE::rx, BASE::ry;
    using BASE::fmm_params, BASE::kernelshape;
  public:
    FMM_Operator(string name, IVec<2> shape,
                 shared_ptr<const BaseFMMInterface> asource, shared_ptr<const BaseFMMInterface> atarget,
                 Array<Vec<3>> xpts, Array<Vec<3>> ypts, Array<Vec<3>> xnv, Array<Vec<3>> ynv,
                 const FMM_Parameters & fmm_params);

    void Mult(const BaseVector & x, BaseVector & y) const override;
    void MultTrans(const BaseVector & x, BaseVector & y) const override;
    BaseMatrix::OperatorInfo GetOperatorInfo () const override;
    FMMOperatorInfo GetFMMInfo () const override;
  };

  extern template class FMM_Operator<double>;
  extern template class FMM_Operator<Complex>;

  // Walk a BaseMatrix tree and return the first Base_FMM_Operator found, or nullptr.
  template <typename TSCAL>
  inline const Base_FMM_Operator<TSCAL> * FindFMMOperator (const BaseMatrix * mat)
  {
    if (!mat) return nullptr;
    if (auto fmm = dynamic_cast<const Base_FMM_Operator<TSCAL>*>(mat))
      return fmm;
    if (auto prod = dynamic_cast<const ProductMatrix*>(mat))
      {
        if (auto found = FindFMMOperator<TSCAL>(prod->SPtrA().get())) return found;
        if (auto found = FindFMMOperator<TSCAL>(prod->SPtrB().get())) return found;
      }
    if (auto sum = dynamic_cast<const SumMatrix*>(mat))
      {
        if (auto found = FindFMMOperator<TSCAL>(sum->SPtrA().get())) return found;
        if (auto found = FindFMMOperator<TSCAL>(sum->SPtrB().get())) return found;
      }
    return nullptr;
  }

}


#endif
