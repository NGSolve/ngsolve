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

struct MapInnerDofs {
  shared_ptr<BitArray> inner;
  shared_ptr<const Array<int>> cluster;
  Array<int> project;
  Array<int> embed;
  size_t size = 0;

  MapInnerDofs(shared_ptr<BitArray> ainner,
               shared_ptr<const Array<int>> acluster = nullptr)
      : inner(ainner), cluster(acluster) {
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
  shared_ptr<SparseMatrixTM<T>> ProjectMatrix(shared_ptr<SparseMatrixTM<T>> m) {
    Array<int> rowi, coli;
    Array<T> vals;
    auto &dofs = *inner;

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

class CudssMatrix : public BaseMatrix {
public:
  CudssMatrix() = delete;
  CudssMatrix(const CudssMatrix &) = delete;
  CudssMatrix(shared_ptr<BaseMatrix> m, shared_ptr<BitArray> freedofs = nullptr,
              shared_ptr<const Array<int>> cluster = nullptr);

  AutoVector CreateRowVector() const override {
    return make_unique<UnifiedVector>(Width());
  }
  AutoVector CreateColVector() const override {
    return make_unique<UnifiedVector>(Height());
  }

  void Analyze();
  void Factor();
  void MultAdd(double s, const BaseVector &x, BaseVector &y) const override;
  void MultAdd(Complex s, const BaseVector &x, BaseVector &y) const override;
  int VWidth() const override { return width; }
  int VHeight() const override { return height; }

private:
  void checkCall(cudssStatus_t status) const;
  MapInnerDofs map_inner_dofs;
  cudssHandle_t cudss_handle;
  cudssConfig_t dev_config;
  cudssData_t dev_data;
  cudssMatrix_t dev_matrix, dev_solution, dev_rhs;
  bool is_complex = false;
  int width, height, inner_width, inner_height;
  size_t nze;

  unique_ptr<Vector<Dev<double>>> solution, rhs, values;
  unique_ptr<Vector<Dev<Complex>>> c_solution, c_rhs, c_values;
  Array<Dev<int>> col_indices, row_start;
};

} // namespace ngla

#endif // NGS_CUDSS_HPP
