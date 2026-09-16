#include "fmmoperator.hpp"
#include <map>
#include <mutex>

namespace ngsbem
{
  namespace
  {
    Timer<> GetFMMTimer(const string & name)
    {
      // Register each runtime kernel name once, not once per matrix instance.
      static std::mutex mutex;
      static std::map<string,Timer<>> timers;
      std::lock_guard<std::mutex> guard(mutex);
      return timers.try_emplace(name, name).first->second;
    }
  }

  template <typename TSCAL>
  FMM_Operator<TSCAL> ::
  FMM_Operator(string name, IVec<2> shape,
               shared_ptr<const BaseFMMInterface> asource, shared_ptr<const BaseFMMInterface> atarget,
               Array<Vec<3>> xpts, Array<Vec<3>> ypts, Array<Vec<3>> xnv, Array<Vec<3>> ynv,
               const FMM_Parameters & fmm_params)
    : BASE(std::move(xpts), std::move(ypts), std::move(xnv), std::move(ynv), shape, fmm_params),
      kernel_name(std::move(name)), source(std::move(asource)), target(std::move(atarget)),
      tapply(GetFMMTimer("ngbem fmm apply "+kernel_name)), teval(GetFMMTimer("ngbem fmm apply "+kernel_name+" eval")),
      ttrans(GetFMMTimer("ngbem fmm apply Trans "+kernel_name))
  { }


  template <typename TSCAL>
  void FMM_Operator<TSCAL> :: Mult(const BaseVector & x, BaseVector & y) const
  {
    RegionTimer reg(tapply);

    auto shape = kernelshape;
    auto matx = x.FV<TSCAL>().AsMatrix(xpts.Size(), shape[1]);
    auto maty = y.FV<TSCAL>().AsMatrix(ypts.Size(), shape[0]);

    maty = 0;
    auto singmp = source->CreateMultipoleExpansion(cx, rx, fmm_params);
    ParallelFor (xpts.Size(), [&](int i){
      source->AddSource(*singmp, xpts[i], xnv[i], matx.Row(i));
    });
    singmp->CalcMP();
    auto regmp = target->CreateLocalExpansion(cy, ry, fmm_params);
    ParallelFor (ypts.Size(), [&](int i){
      regmp->AddTarget(ypts[i]);
    });
    regmp->CalcMP(singmp);

    teval.Start();
    ParallelFor (ypts.Size(), [&](int i) {
      target->EvaluateMP(*regmp, ypts[i], ynv[i], maty.Row(i));
    }, TasksPerThread(10));
    teval.Stop();
  }


  template <typename TSCAL>
  void FMM_Operator<TSCAL> :: MultTrans(const BaseVector & x, BaseVector & y) const
  {
    RegionTimer reg(ttrans);

    auto shape = kernelshape;
    auto matx = x.FV<TSCAL>().AsMatrix(ypts.Size(), shape[0]);
    auto maty = y.FV<TSCAL>().AsMatrix(xpts.Size(), shape[1]);

    maty = 0;
    auto singmp = target->CreateMultipoleExpansion(cy, ry, fmm_params);
    ParallelFor (ypts.Size(), [&](int i){
      target->AddSource(*singmp, ypts[i], ynv[i], matx.Row(i));
    });
    singmp->CalcMP();
    auto regmp = source->CreateLocalExpansion(cx, rx, fmm_params);
    ParallelFor (xpts.Size(), [&](int i){
      regmp->AddTarget(xpts[i]);
    });
    regmp->CalcMP(singmp);
    ParallelFor (xpts.Size(), [&](int i) {
      source->EvaluateMP(*regmp, xpts[i], xnv[i], maty.Row(i));
    }, TasksPerThread(10));
  }


  template <typename TSCAL>
  BaseMatrix::OperatorInfo FMM_Operator<TSCAL> :: GetOperatorInfo () const
  {
    return { string("FMM_Operator ")+kernel_name, this->Height(), this->Width() };
  }


  template <typename TSCAL>
  FMMOperatorInfo FMM_Operator<TSCAL> :: GetFMMInfo () const
  {
    static Timer t("ngbem fmm build info trees"); RegionTimer reg(t);

    auto singmp = source->CreateMultipoleExpansion(cx, rx, fmm_params);
    auto shape = kernelshape;
    Vector<TSCAL> zero_vals(shape[1]);
    zero_vals = TSCAL(0);
    for (size_t i = 0; i < xpts.Size(); i++)
      source->AddSource(*singmp, xpts[i], xnv[i], zero_vals);

    auto regmp = target->CreateLocalExpansion(cy, ry, fmm_params);
    for (size_t i = 0; i < ypts.Size(); i++)
      regmp->AddTarget(ypts[i]);

    // S2R walk: allocates target multipoles and counts translations + direct-evaluation pairs.
    auto counts = regmp->CollectM2LStatistics(singmp);

    FMMOperatorInfo info;
    info.kernel_name = kernel_name;
    info.source_size = xpts.Size();
    info.target_size = ypts.Size();
    info.kappa = singmp->GetKappa();
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

  template class FMM_Operator<double>;
  template class FMM_Operator<Complex>;
}
