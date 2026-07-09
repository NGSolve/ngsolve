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

  MapInnerDofs() {}

  void Init(shared_ptr<BitArray> ainner,
               shared_ptr<const Array<int>> acluster = nullptr)
  {
    inner = ainner;
    cluster = acluster;
    have_value_map = false;
    if (!inner && !cluster) {
      size = 0;
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

  operator bool() const { return inner || cluster; }

  template <typename T>
  void Project(FlatVector<T> dst, FlatVector<T> src) const {
    for (size_t i = 0; i < project.Size(); i++)
      dst[i] = src[project[i]];
  }

  template <typename T> void Embed(T &dst, const T &src) const {
    for (size_t i : Range(embed)) {
      if (embed[i] >= 0)
        dst[i] = src[embed[i]];
      else
        dst[i] = 0.0;
    }
  }

  template <typename T>
  void EmbedAdd(FlatVector<T> dst, FlatVector<T> src, T scale) const {
    for (size_t i : Range(embed))
      if (embed[i] >= 0)
        dst[i] += scale * src[embed[i]];
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
    return make_unique<VVector<double>>(Width());
  }

  AutoVector CreateColVector() const override {
    return make_unique<VVector<double>>(Height());
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
