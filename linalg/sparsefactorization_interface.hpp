#ifndef SPARSEFACTORIZATION_INTERFACE_HPP
#define SPARSEFACTORIZATION_INTERFACE_HPP

#include "basevector.hpp"
#include "sparsecholesky.hpp"
#include "sparsematrix.hpp"
#include "../ngstd/python_ngstd.hpp"

namespace ngla {

  void ExportSparseFactorizationInterface(py::module &m);

struct MapInnerDofs {
  shared_ptr<BitArray> inner;
  shared_ptr<const Array<int>> cluster;
  Array<int> project;
  Array<int> embed;
  size_t size = 0;
  Array<size_t> value_map;
  bool have_value_map = false;
  int bs = 1;                   // entry size of a blocked matrix

  MapInnerDofs() {}

  void Init(shared_ptr<BitArray> ainner,
               shared_ptr<const Array<int>> acluster = nullptr,
               int abs = 1, size_t nblocks = 0)
  {
    inner = ainner;
    cluster = acluster;
    bs = abs;
    have_value_map = false;
    if (!inner && !cluster) {
      if (bs == 1) {
        size = 0;
        return;
      }
      size = nblocks;
      project.SetSize(size);
      embed.SetSize(size);
      for (size_t i : Range(size)) {
        project[i] = i;
        embed[i] = i;
      }
      return;
    }
    if (inner) {
      size = inner->NumSet();
      project.SetSize(size);
      embed.SetSize(inner->Size());
      int j = 0;
      for (int i = 0; i < inner->Size(); i++) {
        if ((*inner)[i]) {
          project[j] = i;
          embed[i] = j++;
        } else
          embed[i] = -1;
      }
      return;
    }

    int j = 0;
    for (int i = 0; i < cluster->Size(); i++) {
      if ((*cluster)[i]) {
        project.Append(i);
        embed.Append(j++);
      } else {
        embed.Append(-1);
      }
    }
    size = project.Size();
  }

  operator bool() const { return inner || cluster || bs > 1; }

  template <typename T>
  void Project(FlatVector<T> dst, FlatVector<T> src) const {
    for (size_t i = 0; i < project.Size(); i++)
      for (int k = 0; k < bs; k++)
        dst[i*bs+k] = src[project[i]*bs+k];
  }

  template <typename T> void Embed(T &dst, const T &src) const {
    for (size_t i : Range(embed))
      for (int k = 0; k < bs; k++) {
        if (embed[i] >= 0)
          dst[i*bs+k] = src[embed[i]*bs+k];
        else
          dst[i*bs+k] = 0.0;
      }
  }

  template <typename T>
  void EmbedAdd(FlatVector<T> dst, FlatVector<T> src, T scale) const {
    for (size_t i : Range(embed))
      if (embed[i] >= 0)
        for (int k = 0; k < bs; k++)
          dst[i*bs+k] += scale * src[embed[i]*bs+k];
  }

  template <typename T>
  shared_ptr<SparseMatrixTM<T>>
  ProjectMatrix(shared_ptr<const SparseMatrixTM<T>> m,
                shared_ptr<SparseMatrixTM<T>> cached = nullptr) {
    auto vals_ori = m->GetValues();

    if (cached && have_value_map && cached->GetValues().Size() == value_map.Size()) {
      auto dst = cached->GetValues();
      for (size_t k : Range(value_map))
        dst[k] = vals_ori[value_map[k]];
      cached->SetSPD(m->IsSPD());
      return cached;
    }

    Array<int> rowi, coli;
    Array<T> vals;
    // auto &dofs = *inner;

    auto is_used = [this](int i, int j) {
      if (inner)
        return (*inner)[i] && (*inner)[j];
      return (*cluster)[i] == (*cluster)[j];
    };

    for (auto i : project)
      for (auto j : m->GetRowIndices(i))
        if (is_used(i, j)) {
          rowi.Append(embed[i]);
          coli.Append(embed[j]);
          vals.Append(vals_ori[m->GetPosition(i, j)]);
        }

    auto res = SparseMatrixTM<T>::CreateFromCOO(rowi, coli, vals,
                                                project.Size(), project.Size());
    res->SetSPD(m->IsSPD());

    shared_ptr<SparseMatrixTM<T>> result;
    if(dynamic_cast<const SparseMatrixSymmetric<T>*>(m.get()))
        result = make_shared<SparseMatrixSymmetric<T>>(*res);
    else
        result = res;

    value_map.SetSize(result->GetValues().Size());
    for (size_t r : Range(size))
      for (auto c : result->GetRowIndices(r))
        value_map[result->GetPosition(r, c)] = m->GetPosition(project[r], project[c]);
    have_value_map = true;

    return result;
  }

  template <int N, typename T>
  shared_ptr<SparseMatrixTM<T>>
  ProjectMatrixBlocked(shared_ptr<const SparseMatrix<Mat<N, N, T>>> m,
                       shared_ptr<SparseMatrixTM<T>> cached = nullptr) {
    auto block_vals = m->GetValues();
    FlatArray<T> vals_ori(block_vals.Size()*N*N, (T*)(void*)block_vals.Data());

    if (cached && have_value_map && cached->GetValues().Size() == value_map.Size()) {
      auto dst = cached->GetValues();
      for (size_t k : Range(value_map))
        dst[k] = vals_ori[value_map[k]];
      cached->SetSPD(m->IsSPD());
      return cached;
    }

    // symmetric storage keeps the lower triangle only, expand both of them
    bool sym_storage =
        dynamic_cast<const SparseMatrixSymmetric<Mat<N, N, T>>*>(m.get()) != nullptr;

    auto is_used = [this](int i, int j) {
      if (inner)
        return bool((*inner)[i] && (*inner)[j]);
      if (cluster)
        return bool((*cluster)[i] == (*cluster)[j]);
      return true;
    };

    Array<int> rowi, coli;
    Array<T> vals;
    auto append_block = [&](int bi, int bj, size_t pos, bool transposed) {
      for (int k = 0; k < N; k++)
        for (int l = 0; l < N; l++) {
          rowi.Append(embed[bi]*N + k);
          coli.Append(embed[bj]*N + l);
          vals.Append(vals_ori[pos*N*N + (transposed ? l*N+k : k*N+l)]);
        }
    };

    for (auto i : project)
      for (auto j : m->GetRowIndices(i)) {
        if (!is_used(i, j)) continue;
        auto pos = m->GetPosition(i, j);
        append_block(i, j, pos, false);
        if (sym_storage && i != j)
          append_block(j, i, pos, true);
      }

    auto result = SparseMatrixTM<T>::CreateFromCOO(rowi, coli, vals,
                                                   size*N, size*N);
    result->SetSPD(m->IsSPD());

    value_map.SetSize(result->GetValues().Size());
    for (size_t r : Range(size*N))
      for (auto c : result->GetRowIndices(r)) {
        int bi = project[r/N], k = r%N;
        int bj = project[c/N], l = c%N;
        bool transposed = sym_storage && bj > bi;
        auto pos = transposed ? m->GetPosition(bj, bi) : m->GetPosition(bi, bj);
        value_map[result->GetPosition(r, c)] =
            pos*N*N + (transposed ? l*N+k : k*N+l);
      }
    have_value_map = true;

    return result;
  }
};

bool IsMatrixSymmetric(shared_ptr<const BaseSparseMatrix> mat, double tol = 0);
shared_ptr<BaseSparseMatrix> ExtractTri(shared_ptr<const BaseSparseMatrix> mat, bool lower=true);

class SparseFactorizationInterface : public SparseFactorization {
protected:
  shared_ptr<const BaseSparseMatrix> inner_mat;
  shared_ptr<BaseVector> inner_rhs, inner_solution;
  MapInnerDofs map_inner_dofs;
  bool is_complex = false;
  xbool is_symmetric = maybe;
  bool is_symmetric_storage = false;
  bool is_analyzed = false;
  int width, height, inner_width, inner_height;

public:
  SparseFactorizationInterface() = delete;
  SparseFactorizationInterface(shared_ptr<const BaseMatrix> m,
                               shared_ptr<BitArray> ainner = nullptr,
                               shared_ptr<const Array<int>> acluster = nullptr);

  virtual ~SparseFactorizationInterface() {}

  void SetSubset(shared_ptr<BitArray> inner, shared_ptr<const Array<int>> cluster) override;

  AutoVector CreateRowVector() const override {
    return matrix.lock()->CreateColVector();
  }

  AutoVector CreateColVector() const override {
    return matrix.lock()->CreateRowVector();
  }

  shared_ptr<const BaseSparseMatrix> GetInnerMatrix() const {
    return inner_mat;
  }

  void MultAdd(double s, const BaseVector &x, BaseVector &y) const override;
  void MultAdd(Complex s, const BaseVector &x, BaseVector &y) const override;

  void MultTransAdd(double s, const BaseVector &x, BaseVector &y) const override;
  void MultTransAdd(Complex s, const BaseVector &x, BaseVector &y) const override;
  void MultConjTransAdd(Complex s, const BaseVector &x,
                        BaseVector &y) const override;

  virtual void Update() override;

  virtual void Analyze() {}
  virtual void Factor() {}
  virtual void Solve(const BaseVector &rhs, BaseVector &solution) const = 0;
  virtual void SolveTrans(const BaseVector &rhs, BaseVector &solution) const;
  virtual void SolveConjTrans(const BaseVector &rhs,
                              BaseVector &solution) const;

  bool IsSymmetricStorage() const { return is_symmetric_storage; }
  xbool IsSymmetric() const override { return is_symmetric; }
  bool IsSPD() const { return inner_mat && inner_mat->IsSPD(); }
};

} // namespace ngla

#endif // SPARSEFACTORIZATION_INTERFACE_HPP
