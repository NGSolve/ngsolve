#ifndef FILE_FMMOPERATOR
#define FILE_FMMOPERATOR
#include <mptools.hpp>
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
      
    AutoVector CreateRowVector () const override
    {
      return make_unique<VVector<TSCAL>>(xpts.Size() * kernelshape[1]);
    }
    AutoVector CreateColVector () const override
    {
      return make_unique<VVector<TSCAL>>(ypts.Size() * kernelshape[0]);
    }

    virtual FMMOperatorInfo GetFMMInfo () const = 0;
  };


  
  template <typename KERNEL>
  class FMM_Operator : public Base_FMM_Operator<typename KERNEL::value_type>
  {
    KERNEL kernel;
    typedef Base_FMM_Operator<typename KERNEL::value_type> BASE;
    using BASE::xpts, BASE::ypts, BASE::xnv, BASE::ynv, BASE::cx, BASE::cy, BASE::rx, BASE::ry;
    using BASE::fmm_params;
  public:
    FMM_Operator(KERNEL _kernel, Array<Vec<3>> _xpts, Array<Vec<3>> _ypts,
                 Array<Vec<3>> _xnv, Array<Vec<3>> _ynv, const FMM_Parameters & fmm_params)
      : BASE(std::move(_xpts), std::move( _ypts), std::move(_xnv), std::move(_ynv), KERNEL::Shape(), fmm_params),
      kernel(_kernel)
    {

      /*
      // build matrix block for testing
      int nump = 3;
      auto shape = KERNEL::Shape();
      Vector<typename KERNEL::value_type> x(nump*shape[1]);
      Vector<typename KERNEL::value_type> y(nump*shape[0]);
      auto matx = x.AsMatrix(nump, shape[1]);
      auto maty = y.AsMatrix(nump, shape[0]);      

      Matrix<typename KERNEL::value_type> mat(y.Size(),x.Size());
      for (int i = 0; i < x.Size(); i++)
        {
          x = 0.0;
          x(i) = 1;
          
          maty = 0;
          auto singmp = kernel.CreateMultipoleExpansion (cx, rx);
          ParallelFor (nump, [&](int i){
            kernel.AddSource(*singmp, xpts[i], xnv[i], matx.Row(i));
          });
          singmp->CalcMP();
          
          auto regmp = kernel.CreateLocalExpansion (cy, ry);
          ParallelFor (nump, [&](int i){
            regmp->AddTarget(xpts[i]);
          });
          regmp->CalcMP(singmp);
          ParallelFor (nump, [&](int i) {
            kernel.EvaluateMP(*regmp, xpts[i], ynv[i], maty.Row(i)); 
          });

          mat.Col(i) = y;
        }
      *testout << "fmm-matrix = " << endl << mat << endl;
      */
    }

    void Mult(const BaseVector & x, BaseVector & y) const override
    {
      static Timer tall("ngbem fmm apply "+KERNEL::Name()); RegionTimer reg(tall);
      
      auto shape = KERNEL::Shape();
      
      auto matx = x.FV<typename KERNEL::value_type>().AsMatrix(xpts.Size(), shape[1]);
      auto maty = y.FV<typename KERNEL::value_type>().AsMatrix(ypts.Size(), shape[0]);      
      
      maty = 0;
      auto singmp = kernel.source.CreateMultipoleExpansion (cx, rx, fmm_params);
      ParallelFor (xpts.Size(), [&](int i){
        kernel.source.AddSource(*singmp, xpts[i], xnv[i], matx.Row(i));
      });
      singmp->CalcMP();
      auto regmp = kernel.target.CreateLocalExpansion (cy, ry, fmm_params);
      ParallelFor (ypts.Size(), [&](int i){
        regmp->AddTarget(ypts[i]);
      });
      regmp->CalcMP(singmp);

      static Timer teval("ngbem fmm apply "+KERNEL::Name() + " eval");
      teval.Start();
      ParallelFor (ypts.Size(), [&](int i) {
        kernel.target.EvaluateMP(*regmp, ypts[i], ynv[i], maty.Row(i));
      });
      teval.Stop();
    }

    void MultTrans(const BaseVector & x, BaseVector & y) const override
    {
        static Timer tall("ngbem fmm apply Trans "+KERNEL::Name()); RegionTimer reg(tall);

        auto shape = KERNEL::Shape();
        auto matx = x.FV<typename KERNEL::value_type>().AsMatrix(xpts.Size(), shape[0]);
        auto maty = y.FV<typename KERNEL::value_type>().AsMatrix(ypts.Size(), shape[1]);

        maty = 0;
        auto singmp = kernel.target.CreateMultipoleExpansion (cy, ry, fmm_params);
        ParallelFor (ypts.Size(), [&](int i){
          kernel.target.AddSource(*singmp, ypts[i], ynv[i], matx.Row(i));
        });
        singmp->CalcMP();
        auto regmp = kernel.source.CreateLocalExpansion (cx, rx, fmm_params);
        ParallelFor (xpts.Size(), [&](int i){
          regmp->AddTarget(xpts[i]);
        });
        regmp->CalcMP(singmp);
        ParallelFor (xpts.Size(), [&](int i) {
          kernel.source.EvaluateMP(*regmp, xpts[i], xnv[i], maty.Row(i));
        });
    }

    BaseMatrix::OperatorInfo GetOperatorInfo () const override
    {
      return { string("FMM_Operator ")+KERNEL::Name(), this->Height(), this->Width() };
    }

    FMMOperatorInfo GetFMMInfo () const override
    {
      static Timer t("ngbem fmm build info trees"); RegionTimer reg(t);

      auto singmp = kernel.source.CreateMultipoleExpansion (cx, rx, fmm_params);
      auto shape = KERNEL::Shape();
      Vector<typename KERNEL::value_type> zero_vals(shape[1]);
      zero_vals = typename KERNEL::value_type(0);
      for (size_t i = 0; i < xpts.Size(); i++)
        kernel.source.AddSource(*singmp, xpts[i], xnv[i], zero_vals);

      auto regmp = kernel.target.CreateLocalExpansion (cy, ry, fmm_params);
      for (size_t i = 0; i < ypts.Size(); i++)
        regmp->AddTarget(ypts[i]);

      // S2R walk: allocates target multipoles and counts translations + direct-evaluation pairs.
      auto counts = regmp->CollectM2LStatistics(singmp);

      FMMOperatorInfo info;
      info.kernel_name = KERNEL::Name();
      info.source_size = xpts.Size();
      info.target_size = ypts.Size();
      auto k = singmp->Kappa();
      if constexpr (std::is_same_v<decltype(k), Complex>)
        info.kappa = k;
      else
        info.kappa = double(k);
      info.parameters = fmm_params;

      singmp->CollectStatistics(info.source_tree);
      regmp->CollectStatistics(info.target_tree);

      info.total_num_nodes  = info.source_tree.num_nodes  + info.target_tree.num_nodes;
      info.total_num_leaves = info.source_tree.num_leaves + info.target_tree.num_leaves;
      info.total_multipole_coefficients = info.source_tree.total_coefficients + info.target_tree.total_coefficients;
      info.total_memory_bytes = info.source_tree.multipole_bytes + info.target_tree.multipole_bytes;

      info.num_s2r = counts.num_s2r;
      info.num_direct_evaluations = counts.num_direct_evaluations;

      info.source_bbox_radius = rx;
      info.target_bbox_radius = ry;

      double dense = double(info.source_size) * double(info.target_size);
      info.direct_fallback_fraction = (dense > 0) ? double(info.num_direct_evaluations) / dense : 0.0;

      return info;
    }
  };

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
