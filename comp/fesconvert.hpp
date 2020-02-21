#ifndef FILE_FESCONVERT
#define FILE_FESCONVERT

namespace ngcomp
{

  shared_ptr<BaseMatrix> ConvertOperator (shared_ptr<FESpace> spacea, shared_ptr<FESpace> spaceb, VorB vb, LocalHeap & lh,
					  shared_ptr<DifferentialOperator> diffop = nullptr, const Region * reg = NULL,
					  bool localop = false, bool parmat = true, bool use_simd = true);

} // namespace ngcomp

#endif
