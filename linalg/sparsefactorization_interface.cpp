#include "sparsefactorization_interface.hpp"
#include <la.hpp>
#include <pybind11/pytypes.h>

namespace ngla {


bool IsMatrixSymmetric(shared_ptr<const BaseSparseMatrix> mat, double tol)
{
    static Timer timer("IsMatrixSymmetric");
    RegionTimer rt(timer);
    if(mat->Height() != mat->Width())
        return false;

    auto p_mat = dynamic_pointer_cast<const MatrixGraph>(mat);
    auto p_real = dynamic_pointer_cast<const SparseMatrixTM<double>>(mat);
    auto p_complex = dynamic_pointer_cast<const SparseMatrixTM<Complex>>(mat);

    if(!p_real && !p_complex)
        return false;

    auto is_close = [tol](auto a, auto b) {
        auto diff = std::abs(a-b);
        return diff <= tol * std::max(1.0, std::max(std::abs(a), std::abs(b)));
    };

    auto is_complex = p_complex != nullptr;
    // size_t nze = mat->NZE();

    auto first = mat->GetFirstArray();
    auto cols = mat->GetColIndices();

    for(auto irow : Range(mat->Height()))
    {
        for(auto i : Range(first[irow], first[irow+1]))
        {
            auto icol = cols[i];
            if(icol == irow)
                continue;
            auto pos_symm = mat->GetPosition(icol, irow);
            if(pos_symm < 0)
                return false;

            if(is_complex && !is_close(p_complex->GetValues()[i], p_complex->GetValues()[pos_symm]))
                return false;
            if(!is_complex && !is_close(p_real->GetValues()[i], p_real->GetValues()[pos_symm]))
                return false;
        }
    }
    return true;
}


shared_ptr<BaseSparseMatrix> ExtractTri(shared_ptr<const BaseSparseMatrix> mat, bool lower)
{
    static Timer timer("ExtractTri");
    RegionTimer rt(timer);
    auto p_mat = dynamic_pointer_cast<const MatrixGraph>(mat);
    auto p_real = dynamic_pointer_cast<const SparseMatrixTM<double>>(mat);
    auto p_complex = dynamic_pointer_cast<const SparseMatrixTM<Complex>>(mat);

    if(!p_real && !p_complex)
        throw Exception("ExtractTri: matrix not derived from SparseMatrixTM");

    Array<int> rowi, coli;
    Array<double> vals_real;
    Array<Complex> vals_complex;

    auto first = mat->GetFirstArray();
    auto cols = mat->GetColIndices();

    for(auto irow : Range(mat->Height()))
    {
        for(auto i : Range(first[irow], first[irow+1]))
        {
            auto icol = cols[i];
            if(lower ? icol > irow : icol < irow)
                continue;

            rowi.Append(irow);
            coli.Append(icol);
            if(p_complex)
                vals_complex.Append(p_complex->GetValues()[i]);
            else
                vals_real.Append(p_real->GetValues()[i]);
        }
    }

    if(p_complex)
        return SparseMatrixTM<Complex>::CreateFromCOO(rowi, coli, vals_complex, mat->Height(), mat->Width());
    else
        return SparseMatrixTM<double>::CreateFromCOO(rowi, coli, vals_real, mat->Height(), mat->Width());
}

SparseFactorizationInterface::SparseFactorizationInterface(
    shared_ptr<const BaseMatrix> m, shared_ptr<BitArray> inner_dofs,
    shared_ptr<const Array<int>> cluster)
    : SparseFactorization(dynamic_pointer_cast<const BaseSparseMatrix>(m),
                          inner_dofs, cluster)
{
  width = m->Width();
  height = m->Height();

  is_complex = m->IsComplex();
}

void SparseFactorizationInterface::SetSubset(shared_ptr<BitArray> inner, shared_ptr<const Array<int>> cluster) {
    SparseFactorization::SetSubset(inner, cluster);
    map_inner_dofs.Init(inner, cluster);
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

    inner_mat = p_sparse_mat;
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

    if(is_symmetric.IsMaybeFalse())
        is_symmetric = IsMatrixSymmetric(inner_mat);

    p_mat = dynamic_pointer_cast<const MatrixGraph>(p_sparse_mat);
  }

  inner_rhs = inner_mat->CreateColVector();
  inner_solution = inner_mat->CreateColVector();

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
