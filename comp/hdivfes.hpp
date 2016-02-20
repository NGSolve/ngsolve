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

    ///
    virtual int GetNDof () const throw();
    ///
    virtual int GetNDofLevel (int level) const;

    ///
    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;
    ///
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;

    ///
    // virtual Table<int> * CreateSmoothingBlocks (int type = 0) const;


    virtual void VTransformMR (int elnr, bool boundary,
			       const SliceMatrix<double> & mat, TRANSFORM_TYPE tt) const;
    virtual void VTransformMC (int elnr, bool boundary,
			       const SliceMatrix<Complex> & mat, TRANSFORM_TYPE tt) const { ; }

    virtual void VTransformVR (int elnr, bool boundary,
			       const FlatVector<double> & vec, TRANSFORM_TYPE tt) const;
    virtual void VTransformVC (int elnr, bool boundary,
			       const FlatVector<Complex> & vec, TRANSFORM_TYPE tt) const { ; }

    void GetTransformationFactors (int elnr, FlatVector<> & fac) const;

    virtual string GetClassName () const
    {
      return "RaviartThomasFESpace";
    }
  };

}

#endif
