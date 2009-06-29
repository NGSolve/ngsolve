#ifndef FILE_H1LOFE
#define FILE_H1LOFE

/*********************************************************************/
/* File:   h1lofe.hpp                                                */
/* Author: Start                                                     */
/* Date:   29. Jun. 2009                                              */
/*********************************************************************/

namespace ngfem
{




  template <int DIM>
  class DShapeElement
  {
    double * data;
  public:
    DShapeElement (double * adata) : data(adata) { ; }
    void operator= (AutoDiff<DIM> ad) 
    { 
      for (int i = 0; i < DIM; i++) 
        data[i] = ad.DValue(i); 
    }
  };

  template <int DIM>
  class DShapeAssign
  {
    double * dshape;
  public:
    DShapeAssign (FlatMatrixFixWidth<DIM> mat)
    { dshape = &mat(0,0); }

    DShapeAssign (double * adshape)
    { dshape = adshape; }

    DShapeElement<DIM> operator[] (int i) const
    { return DShapeElement<DIM> (dshape + i*DIM); }

    const DShapeAssign Addr (int i) const
    { return DShapeAssign (dshape+i*DIM); } 
  };






  /**
     Base-element for template polymorphism.
     Barton and Nackman Trick
  */

  template <class FEL, int SDIM, int NDOF>
  class T_ScalarFiniteElement2 : public ScalarFiniteElement<SDIM>
  {

  public:
    
  protected:
    enum { DIM = SDIM };

    T_ScalarFiniteElement2 ()
      : ScalarFiniteElement<SDIM> (ELEMENT_TYPE(FEL::ELTYPE), NDOF, int(FEL::ORDER))
    {
      try
	{
	  CalcIPData (ELEMENT_TYPE(FEL::ELTYPE), Spec().ipdata);
	  /*
	    if (!Spec().ipdata.Size())
	    CalcIPData ();
	  */
	}
      catch (Exception & e)
	{
	  e.Append ("In Constructor of finite element ");
	  e.Append (typeid(FEL).name());
	  throw e;
	}
    }

    virtual ~T_ScalarFiniteElement2() { ; }

  public:
    virtual const FlatVector<> GetShapeV (const IntegrationPoint & ip) const
    {
      return FlatVector<> (Spec().ipdata[ip.IPNr()] . shape);
    }

    virtual const FlatMatrix<> GetDShapeV (const IntegrationPoint & ip) const
    {
      return FlatMatrix<> (Spec().ipdata[ip.IPNr()] . dshape);
    }


    const Vec<NDOF> & GetShape (const IntegrationPoint & ip,
				LocalHeap & lh) const
    {
      if (ip.IPNr() != -1)
	// return Spec().ipdata[ip.IPNr()] . shape;
	return  reinterpret_cast<const Vec<NDOF> & >
	  ( Spec().ipdata[ip.IPNr()] . shape(0) );
      else
	{
	  throw Exception ("GetDShape, ipnr == -1");
	}
    }

    const Mat<NDOF,SDIM> & GetDShape (const IntegrationPoint & ip,
				      LocalHeap & lh) const
    {
      if (ip.IPNr() != -1)
	{
	  // return Spec().ipdata[ip.IPNr()] . dshape;
	  return  reinterpret_cast<const Mat<NDOF,SDIM> & >
	    ( Spec().ipdata[ip.IPNr()] . dshape(0,0) );
	}
      else
	{
	  throw Exception ("GetDShape, ipnr == -1");
	}
    }





    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const
    {
      double pt[DIM];
      for (int i = 0; i < DIM; i++) pt[i] = ip(i);
      Spec().T_CalcShape (pt, shape); 
    }

    static void CalcShapeStat (const IntegrationPoint & ip, 
                               FlatVector<> shape) 
    {
      double pt[DIM];
      for (int i = 0; i < DIM; i++) pt[i] = ip(i);
      FEL::T_CalcShape (pt, shape); 
    }
    
    
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<SDIM> dshape) const
    {
      AutoDiff<DIM> adp[DIM];
      for (int i = 0; i < DIM; i++)
        adp[i] = AutoDiff<DIM> (ip(i), i);
      
      DShapeAssign<DIM> ds(dshape); 
      static_cast<const FEL*> (this) -> T_CalcShape (adp, ds);
      
      // FEL::CalcDShapeStat (ip, dshape);
    }


    virtual void 
    CalcMappedDShape (const SpecificIntegrationPoint<DIM,DIM> & sip, 
                      FlatMatrixFixWidth<DIM> dshape) const
    {
      AutoDiff<DIM> adp[DIM];
      
      for (int i = 0; i < DIM; i++)
        adp[i].Value() = sip.IP()(i);
      
      for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
          adp[i].DValue(j) = sip.GetJacobianInverse()(i,j);
      
      DShapeAssign<DIM> ds(dshape); 
      Spec().T_CalcShape (adp, ds);
    }




    //  static  Array<IPDataFix> ipdata;
  private:

    FEL & Spec() { return static_cast<FEL&> (*this); }
    const FEL & Spec() const { return static_cast<const FEL&> (*this); }

    /*
      void CalcIPData () 
      {
      const Array<IntegrationPoint*> & ipts = 
      GetIntegrationRules().GetIntegrationPoints (ELEMENT_TYPE(FEL::ELTYPE));
    
      (*testout) << "New: calc IP Data for element type  " << FEL::ELTYPE 
      << ", ndof = " << GetNDof() << ": " << ipts.Size() << endl;
    
      Spec().ipdata.SetSize (ipts.Size());
      for (int i = 0; i < ipts.Size(); i++)
      {
      FEL::CalcShapeStat (*ipts[i], Spec().ipdata[i] . shape);
      FEL::CalcDShapeStat (*ipts[i], Spec().ipdata[i] . dshape);
      }
      }
    */
  };




  ///
  class FE_Trig0 : public T_ScalarFiniteElement2<FE_Trig0,2,1>
  {
  public:
    enum { SDIM = 2 };
    enum { NDOF = 1 };
    enum { ORDER = 0 };
    enum { ELTYPE = ET_TRIG };

    static IPDataArray ipdata;

    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx x[2], TFA & shape) 
    {
      shape[0] = 1;
    }

    virtual const IntegrationRule & NodalIntegrationRule() const;
  }; 

  ///
  class FE_Trig1 : public T_ScalarFiniteElement2<FE_Trig1,2,3>
  {
  public:
    enum { SDIM = 2 };
    enum { NDOF = 3 };
    enum { ORDER = 1 };
    enum { ELTYPE = ET_TRIG };

    static IPDataArray ipdata;

    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx x[2], TFA & shape) 
    {
      shape[0] = x[0];
      shape[1] = x[1];      
      shape[2] = 1-x[0]-x[1];
    }

    ///
    virtual const IntegrationRule & NodalIntegrationRule() const;
  }; 


  class FE_Trig2 : public T_ScalarFiniteElement2<FE_Trig2,2,6>
  {
  public:
    enum { SDIM = 2 };
    enum { NDOF = 6 };
    enum { ORDER = 2 };
    enum { ELTYPE = ET_TRIG };

    static IPDataArray ipdata;

    ///
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[2], TFA & shape) 
    {
      Tx x = hx[0];
      Tx y = hx[1];
      Tx lam3 = 1-x-y;
    
      shape[0] = x * (2*x-1);
      shape[1] = y * (2*y-1);
      shape[2] = lam3 * (2*lam3-1);
      shape[3] = 4 * y * lam3;
      shape[4] = 4 * x * lam3;
      shape[5] = 4 * x * y;
    }

    virtual const IntegrationRule & NodalIntegrationRule() const;
  }; 







  /* ***************************** Tet *************************************** */



  ///
  class FE_Tet0 : public T_ScalarFiniteElement2<FE_Tet0,3,1>
  {
  public:
    enum { SDIM = 3 };
    enum { NDOF = 1 };
    enum { ORDER = 0 };
    enum { ELTYPE = ET_TET };

    static IPDataArray ipdata;

    ///
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[3], TFA & shape) 
    {
      shape[0] = 1;
    }
  
    ///
    virtual const IntegrationRule & NodalIntegrationRule() const;
  };


  ///
  class FE_Tet1 : public T_ScalarFiniteElement2<FE_Tet1,3,4>
  {
  public:
    enum { SDIM = 3 };
    enum { NDOF = 4 };
    enum { ORDER = 1 };
    enum { ELTYPE = ET_TET };

    static IPDataArray ipdata;

    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[3], TFA & shape) 
    {
      Tx x = hx[0];
      Tx y = hx[1];
      Tx z = hx[2];

      shape[0] = x;
      shape[1] = y;
      shape[2] = z;
      shape[3] = 1-x-y-z;
    }

    ///
    virtual const IntegrationRule & NodalIntegrationRule() const;
    virtual void GetDofs (Array<Dof> & dofs) const;
  };


  ///
  class FE_Tet2 : public T_ScalarFiniteElement2<FE_Tet2,3,10>
  {
  public:
    enum { SDIM = 3 };
    enum { NDOF = 10 };
    enum { ORDER = 2 };
    enum { ELTYPE = ET_TET };

    static IPDataArray ipdata;

    ///
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[3], TFA & shape) 
    {
      Tx x = hx[0];
      Tx y = hx[1];
      Tx z = hx[2];
      Tx lam4 = 1 - x - y - z;
    
      shape[0] = 2 * x * x - x;  
      shape[1] = 2 * y * y - y;
      shape[2] = 2 * z * z - z;
      shape[3] = 2 * lam4 * lam4 - lam4;

      shape[4] = 4 * x * y;
      shape[5] = 4 * x * z;
      shape[6] = 4 * x * lam4;
      shape[7] = 4 * y * z;
      shape[8] = 4 * y * lam4;
      shape[9] = 4 * z * lam4;
    }
  };





  ///
  class FE_Tet2HB : public T_ScalarFiniteElement2<FE_Tet2HB,3,10>
  {
  public:
    enum { SDIM = 3 };
    enum { NDOF = 10 };
    enum { ORDER = 2 };
    enum { ELTYPE = ET_TET };

    static IPDataArray ipdata;

    ///
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[3], TFA & shape) 
    {
      Tx x = hx[0];
      Tx y = hx[1];
      Tx z = hx[2];
      Tx lam4 = 1 - x - y - z;

      shape[0] = x;
      shape[1] = y;
      shape[2] = z;
      shape[3] = lam4;

      shape[4] = 4 * x * y;
      shape[5] = 4 * x * z;
      shape[6] = 4 * x * lam4;
      shape[7] = 4 * y * z;
      shape[8] = 4 * y * lam4;
      shape[9] = 4 * z * lam4;
    }
  };







}




#endif
