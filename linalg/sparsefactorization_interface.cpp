#include "sparsefactorization_interface.hpp"
#include <la.hpp>
#include <pybind11/pytypes.h>

namespace ngla {

SparseFactorizationInterface::SparseFactorizationInterface(
    shared_ptr<const BaseMatrix> m, shared_ptr<BitArray> inner_dofs,
    shared_ptr<const Array<int>> cluster)
    : SparseFactorization(dynamic_pointer_cast<const BaseSparseMatrix>(m),
                          inner_dofs, cluster),
      map_inner_dofs(inner_dofs, cluster) {
  width = m->Width();
  height = m->Height();

  is_complex = m->IsComplex();
}

void SparseFactorizationInterface::Update() {
  auto m = matrix.lock();

  shared_ptr<const MatrixGraph> p_mat = nullptr;

  if (is_complex) {
    shared_ptr<const SparseMatrixTM<Complex>> p_sparse_mat =
        dynamic_pointer_cast<const SparseMatrixTM<Complex>>(m);

    if (p_sparse_mat == nullptr) {
      throw Exception(
          "SparseFactorizationInterface: Input matrix not derived from "
          "SparseMatrixTM<Complex>");
    }

    if (dynamic_cast<const SparseMatrixSymmetric<Complex> *>(m.get())) {
      is_symmetric_storage = true;
      is_symmetric = true;
    }

    if (map_inner_dofs)
      p_sparse_mat = map_inner_dofs.ProjectMatrix(p_sparse_mat);

    p_mat = dynamic_pointer_cast<const MatrixGraph>(p_sparse_mat);
  } else {
    shared_ptr<const SparseMatrixTM<double>> p_sparse_mat =
        dynamic_pointer_cast<const SparseMatrixTM<double>>(m);

    if (p_sparse_mat == nullptr) {
      throw Exception(
          "SparseFactorizationInterface: Input matrix not derived from "
          "ScparseMatrixTM<double>");
    }

    if (dynamic_cast<const SparseMatrixSymmetric<double> *>(m.get())) {
      is_symmetric_storage = true;
      is_symmetric = true;
    }

    if (map_inner_dofs)
      p_sparse_mat = map_inner_dofs.ProjectMatrix(p_sparse_mat);

    inner_mat = p_sparse_mat;
    inner_rhs = p_sparse_mat->CreateColVector();
    inner_solution = p_sparse_mat->CreateColVector();

    auto &mat = *p_sparse_mat;

    p_mat = dynamic_pointer_cast<const MatrixGraph>(p_sparse_mat);
  }

  auto &mat = *p_mat;

  inner_height = mat.Height();
  inner_width = mat.Width();

  if (!is_analyzed) {
    Analyze();
    is_analyzed = true;
  }
  Factor();
}

void SparseFactorizationInterface::MultAdd(double s, const BaseVector &x,
                                           BaseVector &y) const {
  if (is_complex)
    return MultAdd(Complex(s), x, y);

  if (map_inner_dofs)
    map_inner_dofs.Project(inner_rhs->FV<double>(), x.FV<double>());
  else
    inner_rhs->FV<double>() = x.FV<double>();

  Solve(*inner_rhs, *inner_solution);

  if (map_inner_dofs) {
    map_inner_dofs.EmbedAdd(y.FV<double>(), inner_solution->FV<double>(), s);
  } else {
    y.FV<double>() += s * inner_solution->FV<double>();
  }
}

void SparseFactorizationInterface::MultAdd(Complex s, const BaseVector &x,
                                           BaseVector &y) const {
  if (map_inner_dofs)
    map_inner_dofs.Project(inner_rhs->FV<Complex>(), x.FV<Complex>());
  else
    inner_rhs->FV<Complex>() = x.FV<Complex>();

  Solve(*inner_rhs, *inner_solution);

  if (map_inner_dofs)
    map_inner_dofs.EmbedAdd(y.FV<Complex>(), inner_solution->FV<Complex>(), s);
  else
    y.FV<Complex>() += s * inner_solution->FV<Complex>();
}

} // namespace ngla
