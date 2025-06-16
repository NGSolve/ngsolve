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

CudssMatrix::CudssMatrix(shared_ptr<BaseMatrix> m,
                         shared_ptr<BitArray> freedofs,
                         shared_ptr<const Array<int>> cluster)
    : map_inner_dofs(freedofs, cluster) {
  cudss_handle = Get_Cudss_Handle();

  width = m->Width();
  height = m->Height();

  is_complex = m->IsComplex();
  auto mview_type = CUDSS_MVIEW_FULL;
  auto mat_type = CUDSS_MTYPE_GENERAL;
  auto data_type = CUDA_R_64F;

  shared_ptr<MatrixGraph> p_mat = nullptr;

  if (is_complex) {
    data_type = CUDA_C_64F;
    shared_ptr<SparseMatrixTM<Complex>> p_sparse_mat =
        dynamic_pointer_cast<SparseMatrixTM<Complex>>(m);

    if (p_sparse_mat == nullptr) {
      throw Exception("CudssMatrix: Input matrix not derived from "
                      "ScparseMatrixTM<Complex>");
    }

    if (dynamic_cast<const SparseMatrixSymmetric<Complex> *>(m.get())) {
      cout << "cudss: have complex symmetric storage" << endl;
      mview_type = CUDSS_MVIEW_LOWER;
    }

    if (map_inner_dofs)
      p_sparse_mat = map_inner_dofs.ProjectMatrix(p_sparse_mat);

    c_values = make_unique<Vector<Dev<Complex>>>(p_sparse_mat->NZE());
    c_values->H2D(p_sparse_mat->GetValues());
    c_rhs = make_unique<Vector<Dev<Complex>>>(p_sparse_mat->Width());
    c_solution = make_unique<Vector<Dev<Complex>>>(p_sparse_mat->Width());
    nze = c_values->Size();
    p_mat = dynamic_pointer_cast<MatrixGraph>(p_sparse_mat);
  } else {
    shared_ptr<SparseMatrixTM<double>> p_sparse_mat =
        dynamic_pointer_cast<SparseMatrixTM<double>>(m);

    if (p_sparse_mat == nullptr) {
      throw Exception("CudssMatrix: Input matrix not derived from "
                      "ScparseMatrixTM<double>");
    }

    if (dynamic_cast<const SparseMatrixSymmetric<double> *>(m.get())) {
      cout << "cudss: have symmetric storage" << endl;
      mview_type = CUDSS_MVIEW_LOWER;
    }

    if (map_inner_dofs)
      p_sparse_mat = map_inner_dofs.ProjectMatrix(p_sparse_mat);

    values = make_unique<Vector<Dev<double>>>(p_sparse_mat->NZE());
    values->H2D(p_sparse_mat->GetValues());
    nze = values->Size();
    rhs = make_unique<Vector<Dev<double>>>(p_sparse_mat->Width());
    solution = make_unique<Vector<Dev<double>>>(p_sparse_mat->Width());

    auto &mat = *p_sparse_mat;

    if (mat.IsSPD()) {
      cout << "cudss: is spd" << endl;
      mat_type = CUDSS_MTYPE_SPD;
    } else if (mat.IsSymmetric().IsTrue()) {
      cout << "cudss: is symmetric" << endl;
      mat_type = CUDSS_MTYPE_SYMMETRIC;
    } else {
      cout << "cudss: is general" << endl;
      mat_type = CUDSS_MTYPE_SPD;
    }
    cout << "p_sparse mat " << p_sparse_mat.get() << endl;
    p_mat = dynamic_pointer_cast<MatrixGraph>(p_sparse_mat);
  }

  auto &mat = *p_mat;

  inner_height = mat.Height();
  inner_width = mat.Width();

  col_indices = mat.GetColIndices();
  Array<int> h_row_start(mat.Height() + 1);
  for (auto i : mat.GetFirstArray().Range())
    h_row_start[i] = mat.GetFirstArray()[i];
  row_start = h_row_start;

  void *val_data, *rhs_data, *solution_data;
  if (is_complex) {
    val_data = c_values->Data();
    rhs_data = c_rhs->Data();
    solution_data = c_solution->Data();
  } else {
    val_data = values->Data();
    rhs_data = rhs->Data();
    solution_data = solution->Data();
  }

  checkCall(cudssMatrixCreateCsr(&dev_matrix, inner_height, inner_width, nze,
                                 row_start.Data(), nullptr, col_indices.Data(),
                                 val_data, CUDA_R_32I, data_type, mat_type,
                                 mview_type, CUDSS_BASE_ZERO));

  checkCall(cudssMatrixCreateDn(&dev_rhs, inner_width, 1, inner_width, rhs_data,
                                data_type, CUDSS_LAYOUT_COL_MAJOR));

  checkCall(cudssMatrixCreateDn(&dev_solution, inner_width, 1, inner_width,
                                solution_data, data_type,
                                CUDSS_LAYOUT_COL_MAJOR));

  checkCall(cudssConfigCreate(&dev_config));
  checkCall(cudssDataCreate(cudss_handle, &dev_data));

  cout << "cudss start analyze, nze = " << nze << ", height = " << inner_height
       << ", width = " << inner_width << endl;
  auto t0 = WallTime();
  Analyze();
  auto t1 = WallTime();
  cout << "Time to analyze: " << t1 - t0 << " seconds" << endl;
  Factor();
  auto t2 = WallTime();
  cout << "Time to factor:  " << t2 - t1 << " seconds" << endl;
  size_t lu_nnz = 0;
  size_t bytes_written = 0;
  checkCall(cudssDataGet(cudss_handle, dev_data, CUDSS_DATA_LU_NNZ, &lu_nnz,
                         sizeof(lu_nnz), &bytes_written));
  cout << "nze in factorization " << lu_nnz / (1024. * 1024.) << "M" << endl;
  cout << "memory usage (in GiB): "
       << (1.0 * lu_nnz * (is_complex ? sizeof(Complex) : sizeof(double)) /
           (1024.0 * 1024.0 * 1024.0))
       << endl;
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

  if (is_complex)
    return MultAdd(Complex(s), x, y);

  if (map_inner_dofs) {
    Vector<double> h_rhs(map_inner_dofs.size);
    map_inner_dofs.Project(h_rhs, x.FV<double>());
    rhs->H2D(h_rhs);
  } else {
    rhs->H2D(x.FV<double>());
  }

  {
    RegionTimer reg(tsolve);
    checkCall(cudssExecute(cudss_handle, CUDSS_PHASE_SOLVE, dev_config,
                           dev_data, dev_matrix, dev_solution, dev_rhs));
  }

  if (map_inner_dofs) {
    map_inner_dofs.EmbedAdd(y.FV<double>(), D2H(*solution), s);
  } else {
    y.FV<double>() += s * D2H(*solution);
  }
}

void CudssMatrix::MultAdd(Complex s, const BaseVector &x, BaseVector &y) const {
  static Timer t("CudssMatrix::MultAdd"),
      tsolve("CudssMatrix::MultAdd - solve");
  RegionTimer reg(t);

  if (map_inner_dofs) {
    Vector<Complex> h_rhs(map_inner_dofs.size);
    map_inner_dofs.Project(h_rhs, x.FV<Complex>());
    c_rhs->H2D(h_rhs);
  } else {
    c_rhs->H2D(x.FV<Complex>());
  }

  {
    RegionTimer reg(tsolve);
    checkCall(cudssExecute(cudss_handle, CUDSS_PHASE_SOLVE, dev_config,
                           dev_data, dev_matrix, dev_solution, dev_rhs));
  }

  if (map_inner_dofs) {
    auto sol_h = D2H(*c_solution);
    map_inner_dofs.EmbedAdd(y.FV<Complex>(), sol_h, s);
  } else {
    y.FV<Complex>() += s * D2H(*c_solution);
  }
}

void CudssMatrix::checkCall(cudssStatus_t status) const {
  if (status != CUDSS_STATUS_SUCCESS) {
    throw Exception("CudssMatrix: cudss call failed with status " +
                    ToString(status));
  }
}

} // namespace ngla
