#ifndef NGS_CUDSS_HPP
#define NGS_CUDSS_HPP

// partial override of overloaded function (MultAdd)
#pragma nv_diag_suppress 611
#pragma nv_diag_suppress 20013

#include <la.hpp>

#include <cudss.h>

#include "cuda_ngstd.hpp"
using namespace ngs_cuda;

namespace ngla {
cudssHandle_t Get_Cudss_Handle();
}

#include "unifiedvector.hpp"

namespace ngla {
using namespace ngs_cuda;

class CudssMatrix : public BaseMatrix {
public:
  CudssMatrix() = delete;
  CudssMatrix(const CudssMatrix &) = delete;
  CudssMatrix(shared_ptr<BaseMatrix> m);

  AutoVector CreateRowVector() const override {
    return make_unique<UnifiedVector>(Width());
  }
  AutoVector CreateColVector() const override {
    return make_unique<UnifiedVector>(Height());
  }

  void Analyze();
  void Factor();
  void MultAdd(double s, const BaseVector &x, BaseVector &y) const override;

private:
  void checkCall(cudssStatus_t status) const;
  cudssHandle_t cudss_handle;
  cudssConfig_t dev_config;
  cudssData_t dev_data;
  cudssMatrix_t dev_matrix, dev_solution, dev_rhs;

  mutable Vector<Dev<double>> solution, rhs, values;
  Array<Dev<int>> col_indices, row_start;
};

} // namespace ngla

#endif // NGS_CUDSS_HPP
