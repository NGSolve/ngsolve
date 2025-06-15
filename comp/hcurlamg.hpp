#ifndef HCURLAMG_HPP
#define HCURLAMG_HPP

// #include <comp.hpp>
#include <basematrix.hpp>
#include <sparsematrix.hpp>
#include "preconditioner.hpp"

namespace ngcomp
{

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
