#ifndef FILE_HDIVFES
#define FILE_HDIVFES

/*********************************************************************/
/* File:   hdivfes.hpp                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   12. Jan. 2002                                             */
/*********************************************************************/

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
    virtual void Update(LocalHeap & lh);

    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;
    ///
    virtual size_t GetNDof () const throw();
    ///
    virtual size_t GetNDofLevel (int level) const;

    ///
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const;
    ///
    // virtual Table<int> * CreateSmoothingBlocks (int type = 0) const;


    virtual void VTransformMR (ElementId ei, 
			       const SliceMatrix<double> mat, TRANSFORM_TYPE tt) const;
    virtual void VTransformMC (ElementId ei, 
			       const SliceMatrix<Complex> mat, TRANSFORM_TYPE tt) const { ; }

    virtual void VTransformVR (ElementId ei, 
			       const SliceVector<double> vec, TRANSFORM_TYPE tt) const;
    virtual void VTransformVC (ElementId ei, 
			       const SliceVector<Complex> vec, TRANSFORM_TYPE tt) const { ; }

    void GetTransformationFactors (int elnr, FlatVector<> & fac) const;

    virtual string GetClassName () const
    {
      return "RaviartThomasFESpace";
    }
  };

}

#endif
