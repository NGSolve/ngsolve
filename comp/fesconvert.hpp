#ifndef FILE_FESCONVERT
#define FILE_FESCONVERT

namespace ngcomp
{

  shared_ptr<BaseMatrix> ConvertOperator (shared_ptr<FESpace> spacea, shared_ptr<FESpace> spaceb,
					  DifferentialOperator * diffop,
					  VorB vb, const Region * reg, LocalHeap & lh,
					  bool localop = false, bool parmat = true, bool use_simd = true);

} // namespace ngcomp

#endif
