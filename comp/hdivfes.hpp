#ifndef FILE_HDIVFES
#define FILE_HDIVFES

/*********************************************************************/
/* File:   hdivfes.hpp                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   12. Jan. 2002                                             */
/*********************************************************************/

#include "fespace.hpp"

namespace ngcomp
{

  /** 
      Finite Element Space for H(div)
  */

  class RaviartThomasFESpace : public FESpace
  {
    ///
    Array<int> ndlevel;
  public:
    ///
    /*
      RaviartThomasFESpace (shared_ptr<MeshAccess> ama,
      int adim, bool acomplex);
    */
    ///
    RaviartThomasFESpace (shared_ptr<MeshAccess> ama, const Flags& flags, bool parseflags=false);
    ///
    virtual ~RaviartThomasFESpace ();
  
    ///
    virtual const char * GetType() 
    { return "RaviartThomas"; }

    static shared_ptr<FESpace> Create (shared_ptr<MeshAccess> ma, const Flags & flags);

    ///
    void Update() override;

    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;
    ///
    virtual size_t GetNDof () const throw() override;
    ///
    virtual size_t GetNDofLevel (int level) const override;

    ///
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    ///
    // virtual Table<int> * CreateSmoothingBlocks (int type = 0) const;


    virtual void VTransformMR (ElementId ei, 
			       const SliceMatrix<double> mat, TRANSFORM_TYPE tt) const override;
    virtual void VTransformMC (ElementId ei, 
			       const SliceMatrix<Complex> mat, TRANSFORM_TYPE tt) const override { ; }

    virtual void VTransformVR (ElementId ei, 
			       const SliceVector<double> vec, TRANSFORM_TYPE tt) const override;
    virtual void VTransformVC (ElementId ei, 
			       const SliceVector<Complex> vec, TRANSFORM_TYPE tt) const override { ; }

    void GetTransformationFactors (ElementId ei, FlatVector<> fac) const;

    virtual string GetClassName () const override
    {
      return "RaviartThomasFESpace";
    }
  };

}

#endif
