#ifndef FILE_NUMBERFESPACE
#define FILE_NUMBERFESPACE

/*********************************************************************/
/* File:   numberfespace.hpp                                           */
/* Author: Start                                                     */
/* Date:   30.FJan. 2018                                              */
/*********************************************************************/

#include <fespace.hpp>

namespace ngcomp
{

  class NumberFESpace : public FESpace
  {
  public:
    NumberFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool checkflags=false);
    void Update() override;

    // virtual size_t GetNDof() const { return 1; }

    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;
    
    virtual void GetDofNrs (ElementId ei, Array<int> & dnums) const override;
    
    virtual void GetGlobalDofNrs (int gnr, Array<int> & dnums) const override;
    
  };

}
#endif
