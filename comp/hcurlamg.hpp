#ifndef HCURLAMG_HPP
#define HCURLAMG_HPP

#include <comp.hpp>

namespace ngcomp
{
  template<typename SCAL>
  class HCurlAMG_Matrix : public BaseMatrix
  {
  protected:
    size_t size;
    shared_ptr<SparseMatrixTM<SCAL>> mat;
    shared_ptr<BaseJacobiPrecond> smoother;
    shared_ptr<BaseBlockJacobiPrecond> blocksmoother;
    shared_ptr<SparseMatrixTM<double>> prolongation, restriction;
    shared_ptr<SparseMatrixTM<double>> vert_prolongation, vert_restriction;
    shared_ptr<SparseMatrixTM<double>> gradient, trans_gradient;
    shared_ptr<BaseMatrix> node_h1;
    shared_ptr<BaseMatrix> coarse_precond;
    int smoothing_steps = 1;
    bool node_on_each_level;
    bool use_smoothed_prolongation;
    int coarsenings_per_level;
    bool blockjacobi_smoother;
    bool need_vertex_prolongation = false;
    int nv = 0;

  public:
    HCurlAMG_Matrix (bool _node_on_each_level = false,
                     bool _use_smoothed_prolongation = true,
                     int _coarsenings_per_level = 3,
                     bool _blockjacobi_smoother = true)
      : node_on_each_level(_node_on_each_level),
        use_smoothed_prolongation(_use_smoothed_prolongation),
        coarsenings_per_level(_coarsenings_per_level),
        blockjacobi_smoother(_blockjacobi_smoother)
    {}

    HCurlAMG_Matrix(shared_ptr<SparseMatrixTM<SCAL>> _mat,
                    shared_ptr<BitArray> freedofs,
                    FlatArray<IVec<3>> f2e,
                    FlatArray<IVec<2>> e2v,
                    FlatArray<double> edge_weights,
                    FlatArray<double> face_weights,
                    size_t level)
      : HCurlAMG_Matrix()
    {
      for(const auto& verts : e2v)
        nv = max3(nv, verts[0], verts[1]);
      nv++;
      Init(_mat, freedofs, f2e, e2v, edge_weights, face_weights, level);
    }

    void Init(shared_ptr<SparseMatrixTM<SCAL>> _mat,
              shared_ptr<BitArray> freedofs,
              FlatArray<IVec<3>> f2e,
              FlatArray<IVec<2>> e2v,
              FlatArray<double> edge_weights,
              FlatArray<double> face_weights,
              size_t level);

    int VHeight() const override { return size; }
    int VWidth() const override { return size; }
    bool IsComplex() const override { return is_same<SCAL, Complex>(); }
    AutoVector CreateRowVector() const override
    { return mat->CreateColVector(); }
    AutoVector CreateColVector() const override
    { return mat->CreateRowVector(); }

    void Mult(const BaseVector& f, BaseVector& u) const override;
  protected:

    Array<double> CalcEdgeCollapse(FlatArray<IVec<3>> f2e,
                                   FlatArray<IVec<2>> e2v,
                                   FlatArray<double> edge_weights,
                                   FlatArray<double> face_weights,
                                   FlatTable<int> e2f) const;
    struct AMGInfo
    {
      Array<IVec<3>> f2e;
      Array<IVec<2>> e2v;
      Array<double> edge_weights;
      Array<double> face_weights;
      shared_ptr<BitArray> freedofs;
      shared_ptr<SparseMatrixTM<double>> prolongation;
      shared_ptr<SparseMatrixTM<double>> vert_prolongation;
    };
    virtual void BuildCoarseMat(const AMGInfo& cinfo, int level);

    virtual shared_ptr<BitArray>
    GetHCurlFreeDofs(shared_ptr<BitArray> freedofs) const
    {
      return freedofs;
    }

    virtual shared_ptr<BitArray>
    CreateCoarseFreedofs(shared_ptr<BitArray> freedofs,
                         int nce, int ncv,
                         FlatArray<int> edge_map,
                         FlatArray<int> vert_map,
                         FlatArray<IVec<2>> e2v,
                         FlatArray<IVec<2>> ce2v) const;

    virtual shared_ptr<BitArray> GetH1FreeDofs(FlatArray<IVec<2>> e2v,
                                               FlatTable<int> e2f,
                                               shared_ptr<BitArray> freedofs) const;

    AMGInfo CalcCoarsening(FlatArray<double> coll_weights,
                           shared_ptr<BitArray> freedofs,
                           FlatArray<IVec<3>> f2e,
                           FlatArray<IVec<2>> e2v,
                           FlatTable<int> e2f,
                           FlatArray<double> edge_weights,
                           FlatArray<double> face_weights,
                           int nv,
                           int level,
                           int coarsening) const;
  };

  class HCurlAMG : public Preconditioner
  {
  protected:
    shared_ptr<BitArray> freedofs;
    shared_ptr<BaseMatrix> mat;
    ParallelHashTable<IVec<1>, double> edge_weights_ht;
    ParallelHashTable<IVec<3>, double> face_weights_ht;
    shared_ptr<FESpace> fes;
    bool node_on_each_level;
  public:
    HCurlAMG(shared_ptr<BilinearForm> _bfa, const Flags& _flags,
             const string& name = "HCurl_AMG");

    void InitLevel(shared_ptr<BitArray> _freedofs) override
    { freedofs = _freedofs; }

    void AddElementMatrix(FlatArray<int> dnums, const FlatMatrix<double> & elmat,
                          ElementId id, LocalHeap & lh) override
    {
      AddElementMatrixCommon(dnums, elmat, id, lh);
    }
    void AddElementMatrix(FlatArray<int> dnums, const FlatMatrix<Complex> & elmat,
                          ElementId id, LocalHeap & lh) override
    {
      AddElementMatrixCommon(dnums, elmat, id, lh);
    }

    void FinalizeLevel(const BaseMatrix* matrix) override;
    void Update() override {}

    const BaseMatrix& GetMatrix() const override
    { return *mat; }

  private:
    template<typename SCAL>
    void AddElementMatrixCommon(FlatArray<int> dnums,
                                const FlatMatrix<SCAL> & elmat,
                                ElementId id, LocalHeap & lh);
  };

  template<typename SCAL>
  class APhiMatrix : public HCurlAMG_Matrix<SCAL>
  {
    using BASE = HCurlAMG_Matrix<SCAL>;
    size_t hcurlsize;
    size_t h1size;
  public:
    APhiMatrix(shared_ptr<SparseMatrixTM<SCAL>> _mat,
               shared_ptr<BitArray> freedofs,
               FlatArray<IVec<3>> f2e,
               FlatArray<IVec<2>> e2v,
               FlatArray<double> edge_weights,
               FlatArray<double> face_weights,
               size_t level);
  protected:
    void BuildCoarseMat(const typename BASE::AMGInfo& cinfo, int level) override;
    shared_ptr<BitArray>
    GetHCurlFreeDofs(shared_ptr<BitArray> freedofs) const override
    {
      auto fd = make_shared<BitArray>(freedofs->Size());
      fd->Clear();
      for(auto i : Range(hcurlsize))
        if(freedofs->Test(i))
          fd->SetBit(i);
      return fd;
    }
    shared_ptr<BitArray> GetH1FreeDofs(FlatArray<IVec<2>> e2v,
                                       FlatTable<int> e2f,
                                       shared_ptr<BitArray> freedofs) const override
    {
      auto fd = make_shared<BitArray>(h1size);
      fd->Clear();
      for(auto i : Range(h1size))
        if(freedofs->Test(hcurlsize + i))
          fd->SetBit(i);
      return fd;
    }
  shared_ptr<BitArray> CreateCoarseFreedofs(shared_ptr<BitArray> freedofs,
                                            int nce, int ncv,
                                            FlatArray<int> edge_map,
                                            FlatArray<int> vert_map,
                                            FlatArray<IVec<2>> e2v,
                                            FlatArray<IVec<2>> ce2v) const override;
  };

  class APhiHCurlAMG : public HCurlAMG
  {
    shared_ptr<BilinearForm> bfa;
  public:
    APhiHCurlAMG(shared_ptr<BilinearForm> _bfa,
                 const Flags& _flags,
                 const string& name = "A-Phi_HcurlAMG");

    void InitLevel(shared_ptr<BitArray> _freedofs) override
    { freedofs = _freedofs; }

    void FinalizeLevel(const BaseMatrix* matrix) override;
    void Update() override {}

    const BaseMatrix& GetMatrix() const override
    { return *mat; }

    void AddElementMatrix(FlatArray<int> dnums,
                          const FlatMatrix<double> & elmat,
                          ElementId id, LocalHeap & lh) override;
    void AddElementMatrix(FlatArray<int> dnums,
                          const FlatMatrix<Complex> & elmat,
                          ElementId id, LocalHeap & lh) override;
    template<typename SCAL>
    void AddElementMatrixCommon(FlatArray<int> dnums,
                                const FlatMatrix<SCAL> & belmat,
                                ElementId id, LocalHeap & lh);
  };
} // namespace ngcomp

#endif // HCURLAMG_HPP
