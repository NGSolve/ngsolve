#ifndef FILE_SPECIALELEMENTGROUP
#define FILE_SPECIALELEMENTGROUP

#include <ngstd.hpp>
#include <bla.hpp>
#include <specialelement.hpp>

namespace ngla { class BaseVector; }

namespace ngcomp
{
  class FESpace;

  class SpecialElementGroup
  {
  protected:
    // set in Update, the form rebuilds its graph when it changes
    size_t timestamp = 0;

  public:
    template <typename SCAL>
    using TAddElMat = std::function<void(FlatArray<DofId>,FlatArray<DofId>,FlatMatrix<SCAL>,ElementId,LocalHeap&)>;
    typedef std::function<void(FlatArray<DofId>,FlatVector<double>,LocalHeap&)> TAddElVec;
    typedef std::function<void(int,FlatArray<DofId>)> TElDofs;

    virtual ~SpecialElementGroup() { }

    // (re)create the elements
    virtual void Update() { }
    size_t GetTimeStamp() const { return timestamp; }

    // if set, they must be the spaces of the form
    virtual shared_ptr<FESpace> GetTrialSpace() const { return nullptr; }
    virtual shared_ptr<FESpace> GetTestSpace() const { return nullptr; }

    virtual int GetNElements() = 0;
    // trial-space dofs
    virtual void GetDofNrs(std::function<void(int,FlatArray<DofId>)> eldofs) = 0;
    // test-space dofs
    virtual void GetDofNrs2(std::function<void(int,FlatArray<DofId>)> eldofs)
    {
      GetDofNrs(eldofs);
    }

    // addelmat (test dofs, trial dofs, elmat, ..), may be called in parallel
    virtual void Assemble(std::function<void(FlatArray<DofId>,FlatArray<DofId>,FlatMatrix<double>,ElementId,LocalHeap&)> addelmat, LocalHeap& lh) = 0;

    virtual void Assemble(std::function<void(FlatArray<DofId>,FlatArray<DofId>,FlatMatrix<float>,ElementId,LocalHeap&)> addelmat, LocalHeap& lh)
    {
      throw Exception("SpecialElementGroup::Assemble(float) not implemented");
    }

    virtual void Assemble(std::function<void(FlatArray<DofId>,FlatArray<DofId>,FlatMatrix<Complex>,ElementId,LocalHeap&)> addelmat, LocalHeap& lh)
    {
      Assemble([&](FlatArray<DofId> dnumsr, FlatArray<DofId> dnumsc, FlatMatrix<double> elmat, ElementId ei, LocalHeap &lh) {
        HeapReset hr(lh);
        FlatMatrix<Complex> cm(elmat.Height(), elmat.Width(), lh);
        cm = elmat;
        addelmat(dnumsr, dnumsc, cm, ei, lh);
      }, lh);
    }

    virtual void AssembleLinearization(const ngla::BaseVector & lin, TAddElMat<double> addelmat, LocalHeap& lh)
    {
      throw Exception("SpecialElementGroup::AssembleLinearization not implemented");
    }

    // addelvec (test dofs, elvec, ..), may be called in parallel
    virtual void Apply(const ngla::BaseVector & x, TAddElVec addelvec, LocalHeap& lh)
    {
      throw Exception("SpecialElementGroup::Apply not implemented");
    }

    virtual double Energy(const ngla::BaseVector & x, LocalHeap& lh)
    {
      throw Exception("SpecialElementGroup::Energy not implemented");
    }
  };


  // a group given by a list of special elements
  class NGS_DLL_HEADER SpecialElementList : public SpecialElementGroup
  {
  protected:
    Array<unique_ptr<SpecialElement>> elements;

  public:
    int GetNElements() override { return elements.Size(); }
    void GetDofNrs(std::function<void(int,FlatArray<DofId>)> eldofs) override;
    void GetDofNrs2(std::function<void(int,FlatArray<DofId>)> eldofs) override;

    using SpecialElementGroup::Assemble;
    void Assemble(std::function<void(FlatArray<DofId>,FlatArray<DofId>,FlatMatrix<double>,ElementId,LocalHeap&)> addelmat, LocalHeap& lh) override;
    void AssembleLinearization(const ngla::BaseVector & lin, TAddElMat<double> addelmat, LocalHeap& lh) override;
    void Apply(const ngla::BaseVector & x, TAddElVec addelvec, LocalHeap& lh) override;
    double Energy(const ngla::BaseVector & x, LocalHeap& lh) override;
  };
}


#endif
