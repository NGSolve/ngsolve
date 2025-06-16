#include <cudss.h>
#include <la.hpp>
#include <memory>

#include "cuda_linalg.hpp"
#include "ngs_cudss.hpp"

namespace ngla {
cudssHandle_t Get_Cudss_Handle() {
  static Timer tblashandle("CUDA create cublas handle");
  RegionTimer reg(tblashandle);

  static cudssHandle_t handle;
  static bool first_call = true;

  if (first_call) {
    first_call = false;
    cudssCreate(&handle);
  }
  return handle;
}

CudssMatrix::CudssMatrix(shared_ptr<BaseMatrix> m)
    : solution(m->Width()), rhs(m->Height()), values(m->NZE()) {
  cudss_handle = Get_Cudss_Handle();

  auto p_sparse_mat = dynamic_cast<const SparseMatrix<double> *>(m.get());

  if (p_sparse_mat == nullptr) {
    throw Exception(
        "CudssMatrix: Input matrix not derived from ScparseMatrix<double>");
  }

  auto mview_type = CUDSS_MVIEW_FULL;

  if (dynamic_cast<const SparseMatrixSymmetric<double> *>(m.get())) {
    cout << "cudss: have symmetric storage" << endl;
    mview_type = CUDSS_MVIEW_LOWER;
  }

  auto &sparse_mat = *p_sparse_mat;

  auto height = sparse_mat.Height();
  auto width = sparse_mat.Width();

  col_indices = sparse_mat.GetColIndices();
  Array<int> h_row_start(sparse_mat.Height() + 1);
  for (auto i : sparse_mat.GetFirstArray().Range())
    h_row_start[i] = sparse_mat.GetFirstArray()[i];
  row_start = h_row_start;
  values.H2D(sparse_mat.GetValues());

  auto mat_type = CUDSS_MTYPE_GENERAL;
  if (sparse_mat.IsSPD()) {
    cout << "cudss: is spd" << endl;
    mat_type = CUDSS_MTYPE_SPD;
  } else if (sparse_mat.IsSymmetric().IsTrue()) {
    cout << "cudss: is symmetric" << endl;
    mat_type = CUDSS_MTYPE_SYMMETRIC;
  } else {
    cout << "cudss: is general" << endl;
    mat_type = CUDSS_MTYPE_SPD;
  }
  checkCall(cudssMatrixCreateCsr(&dev_matrix, height, width, sparse_mat.NZE(),
                                 row_start.Data(), nullptr, col_indices.Data(),
                                 values.Data(), CUDA_R_32I, CUDA_R_64F,
                                 mat_type, mview_type, CUDSS_BASE_ZERO));

  checkCall(cudssMatrixCreateDn(&dev_rhs, width, 1, width, rhs.Data(),
                                CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR));

  checkCall(cudssMatrixCreateDn(&dev_solution, width, 1, width, solution.Data(),
                                CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR));

  checkCall(cudssConfigCreate(&dev_config));
  checkCall(cudssDataCreate(cudss_handle, &dev_data));

  Analyze();
  Factor();
}

void CudssMatrix::Analyze() {
  static Timer t("CudssMatrix::Analyze");
  RegionTimer reg(t);
  checkCall(cudssExecute(cudss_handle, CUDSS_PHASE_ANALYSIS, dev_config,
                         dev_data, dev_matrix, dev_solution, dev_rhs));
}

void CudssMatrix::Factor() {
  static Timer t("CudssMatrix::Factor");
  RegionTimer reg(t);
  checkCall(cudssExecute(cudss_handle, CUDSS_PHASE_FACTORIZATION, dev_config,
                         dev_data, dev_matrix, dev_solution, dev_rhs));
}

void CudssMatrix::MultAdd(double s, const BaseVector &x, BaseVector &y) const {
  static Timer t("CudssMatrix::MultAdd"),
      tsolve("CudssMatrix::MultAdd - solve");
  RegionTimer reg(t);

  rhs.H2D(x.FV<double>());

  {
    RegionTimer reg(tsolve);
    checkCall(cudssExecute(cudss_handle, CUDSS_PHASE_SOLVE, dev_config,
                           dev_data, dev_matrix, dev_solution, dev_rhs));
  }

  y.FV<double>() += s * D2H(solution);
}

void CudssMatrix::checkCall(cudssStatus_t status) const {
  if (status != CUDSS_STATUS_SUCCESS) {
    throw Exception("CudssMatrix: cudss call failed with status " +
                    ToString(status));
  }
}
} // namespace ngla
