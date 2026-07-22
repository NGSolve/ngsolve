#ifndef FILE_SPECIALELEMENTGROUP
#define FILE_SPECIALELEMENTGROUP

#include <bla.hpp>
using namespace ngbla;



namespace ngcomp
{
  class SpecialElementGroup
  {
  public:
    virtual ~SpecialElementGroup() { }

    virtual void Update() { } 

    virtual int GetNElements() = 0;
    virtual void GetDofNrs(std::function<void(int,FlatArray<DofId>)> eldofs) = 0;
    
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
  };
}


#endif
