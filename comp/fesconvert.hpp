#ifndef FILE_FESCONVERT
#define FILE_FESCONVERT

namespace ngcomp
{

  shared_ptr<BaseMatrix> ConvertOperator (shared_ptr<FESpace> spacea, shared_ptr<FESpace> spaceb, VorB vb, LocalHeap & lh,
					  shared_ptr<DifferentialOperator> diffop = nullptr, shared_ptr<CoefficientFunction> trial_cf = nullptr,
					  const Region * reg = NULL, shared_ptr<BitArray> range_dofs = nullptr, bool localop = false, bool parmat = true,
					  bool use_simd = true, int bonus_intorder_ab = 0, int bonus_intorder_bb = 0, bool geom_free = false);

} // namespace ngcomp

#endif
