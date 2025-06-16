#ifndef HCURLAMG_HPP
#define HCURLAMG_HPP

// #include <comp.hpp>
#include <basematrix.hpp>
#include <sparsematrix.hpp>
#include "preconditioner.hpp"

namespace ngcomp
{

  struct HCurlAMG_Parameters
  {
    int verbose = 0;
    int smoothing_steps = 3;
    int coarsenings_per_level = 3; // number of binary coarsenings per level
    bool block_smoother = true;    // block or point smoother ?
    bool use_smoothed_prolongation = true;  
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
    HCurlAMG_Parameters param;

  public:
    HCurlAMG(shared_ptr<BilinearForm> _bfa, const Flags& _flags,
             const string& name = "HCurl_AMG");

    static DocInfo GetDocu ();    
    
    void InitLevel(shared_ptr<BitArray> _freedofs) override
    { freedofs = _freedofs; }

    void AddElementMatrix(FlatArray<int> dnums, FlatMatrix<double> elmat,
                          ElementId id, LocalHeap & lh) override
    {
      AddElementMatrixCommon(dnums, elmat, id, lh);
    }
    void AddElementMatrix(FlatArray<int> dnums, FlatMatrix<Complex> elmat,
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
                                FlatMatrix<SCAL> elmat,
                                ElementId id, LocalHeap & lh);
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
                          FlatMatrix<double> elmat,
                          ElementId id, LocalHeap & lh) override;
    void AddElementMatrix(FlatArray<int> dnums,
                          FlatMatrix<Complex> elmat,
                          ElementId id, LocalHeap & lh) override;
    template<typename SCAL>
    void AddElementMatrixCommon(FlatArray<int> dnums,
                                FlatMatrix<SCAL> belmat,
                                ElementId id, LocalHeap & lh);
  };
} // namespace ngcomp

#endif // HCURLAMG_HPP
