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

  MapInnerDofs() {}

  void Init(shared_ptr<BitArray> ainner,
               shared_ptr<const Array<int>> acluster = nullptr)
  {
    inner = ainner;
    cluster = acluster;
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
  ProjectMatrix(shared_ptr<const SparseMatrixTM<T>> m) {
    Array<int> rowi, coli;
    Array<T> vals;
    // auto &dofs = *inner;

    auto vals_ori = m->GetValues();

    auto &cluster_array = *cluster;
    auto &inner_bitarray = *inner;
    auto is_used = [this, &inner_bitarray, &cluster_array](int i, int j) {
      if (inner)
        return inner_bitarray[i] && inner_bitarray[j];
      return cluster_array[i] == cluster_array[j];
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
    return res;
  }
};

class SparseFactorizationInterface : public SparseFactorization {
protected:
  shared_ptr<const BaseSparseMatrix> inner_mat;
  shared_ptr<BaseVector> inner_rhs, inner_solution;
  MapInnerDofs map_inner_dofs;
  bool is_complex = false;
  bool is_symmetric = false;
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

  virtual void Update() override;

  virtual void Analyze() {}
  virtual void Factor() {}
  virtual void Solve(const BaseVector &rhs, BaseVector &solution) const = 0;
};

} // namespace ngla

#endif // SPARSEFACTORIZATION_INTERFACE_HPP
