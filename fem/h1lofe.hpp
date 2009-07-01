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

  template <class FEL, ELEMENT_TYPE ET, int NDOF, int ORDER>
  class T_ScalarFiniteElement2 : public ScalarFiniteElement<ET_trait<ET>::DIM>
  {

  public:
    
  protected:
    enum { DIM = ET_trait<ET>::DIM };

    T_ScalarFiniteElement2 ()
    // : ScalarFiniteElement<SDIM> (ELEMENT_TYPE(FEL::ELTYPE), NDOF, int(FEL::ORDER))
      : ScalarFiniteElement<DIM> (ET, NDOF, ORDER)
    {
      try
	{
	  // CalcIPData (ELEMENT_TYPE(FEL::ELTYPE), Spec().ipdata);
	  CalcIPData (ET, Spec().ipdata);
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

    const Mat<NDOF,DIM> & GetDShape (const IntegrationPoint & ip,
				     LocalHeap & lh) const
    {
      if (ip.IPNr() != -1)
	{
	  // return Spec().ipdata[ip.IPNr()] . dshape;
	  return  reinterpret_cast<const Mat<NDOF,DIM> & >
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
			     FlatMatrixFixWidth<DIM> dshape) const
    {
      AutoDiff<DIM> adp[DIM];
      for (int i = 0; i < DIM; i++)
        adp[i] = AutoDiff<DIM> (ip(i), i);
      
      DShapeAssign<DIM> ds(dshape); 
      static_cast<const FEL*> (this) -> T_CalcShape (adp, ds);
      
      // FEL::CalcDShapeStat (ip, dshape);
    }

    static void CalcDShapeStat (const IntegrationPoint & ip, 
				FlatMatrixFixWidth<DIM> dshape) 
    {
      AutoDiff<DIM> adp[DIM];
      for (int i = 0; i < DIM; i++)
        adp[i] = AutoDiff<DIM> (ip(i), i);
      
      DShapeAssign<DIM> ds(dshape); 
      FEL::T_CalcShape (adp, ds);
      
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
  class FE_Segm0 : public T_ScalarFiniteElement2<FE_Segm0,ET_SEGM,1,0>
  {
  public:
    static IPDataArray ipdata;

    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx x[1], TFA & shape) 
    {
      shape[0] = 1;
    }
  }; 

  /*
  ///
  class FE_SegmDummy : public T_ScalarFiniteElement2<FE_SegmDummy,ET_SEGM,0,0>
  {
  public:
    static IPDataArray ipdata;

    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx x[1], TFA & shape) 
    {
      ;
    }

  }; 
  */

  ///
  class FE_Segm1 : public T_ScalarFiniteElement2<FE_Segm1,ET_SEGM,2,1>
  {
  public:
    static IPDataArray ipdata;

    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx x[1], TFA & shape) 
    {
      shape[0] = x[0];
      shape[1] = 1-x[0];
    }
  }; 


  ///
  class FE_Segm1L2 : public T_ScalarFiniteElement2<FE_Segm1L2,ET_SEGM,2,1>
  {
  public:
    static IPDataArray ipdata;

    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx x[1], TFA & shape) 
    {
      shape[0] = 1;
      shape[1] = 2*x[0]-1;
    }
  }; 


  ///
  class FE_Segm2 : public T_ScalarFiniteElement2<FE_Segm2,ET_SEGM,3,2>
  {
  public:

    static IPDataArray ipdata;

    ///
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[1], TFA & shape) 
    {
      Tx x = hx[0];

      shape[0] = 2*x*x - x;
      shape[1] = 2*x*x - 3*x + 1;  
      shape[2] = 4 * x * (1-x);
    }
  }; 


  ///
  class FE_Segm2HB : public T_ScalarFiniteElement2<FE_Segm2HB,ET_SEGM,3,2>
  {
  public:

    static IPDataArray ipdata;

    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[1], TFA & shape) 
    {
      Tx x = hx[0];

      shape[0] = x;
      shape[1] = 1-x;
      shape[2] = 4 * x * (1-x);
    }
  }; 



  ///
  class FE_Segm2L2 :public T_ScalarFiniteElement2<FE_Segm2L2,ET_SEGM,3,2>
  {
  public:

    static IPDataArray ipdata;


    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[1], TFA & shape) 
    {
      Tx x = hx[0];

      shape[0] = 1;
      shape[1] = 2*x-1;
      shape[2] = (2*x-1)*(2*x-1)-1.0/3.0;
    }
  }; 











  /* ***************************** Trig Elements *************************************** */

  ///
  class FE_Trig0 : public T_ScalarFiniteElement2<FE_Trig0,ET_TRIG,1,0>
  {
  public:
    static IPDataArray ipdata;

    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx x[2], TFA & shape) 
    {
      shape[0] = 1;
    }

    virtual const IntegrationRule & NodalIntegrationRule() const;
  }; 

  ///
  class FE_Trig1 : public T_ScalarFiniteElement2<FE_Trig1,ET_TRIG,3,1>
  {
  public:
    static IPDataArray ipdata;

    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx x[2], TFA & shape) 
    {
      shape[0] = x[0];
      shape[1] = x[1];      
      shape[2] = 1-x[0]-x[1];
    }

    virtual const IntegrationRule & NodalIntegrationRule() const;
  }; 


  ///
  class FE_Trig2 : public T_ScalarFiniteElement2<FE_Trig2,ET_TRIG,6,2>
  {
  public:
    static IPDataArray ipdata;

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


  /* ***************************** Quad *************************************** */


  class FE_Quad0 : public T_ScalarFiniteElement2<FE_Quad0,ET_QUAD,1,0>
  {
  public:
    static IPDataArray ipdata;

    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[2], TFA & shape) 
    {
      shape[0] = 1.0;
    }
			   
    virtual const IntegrationRule & NodalIntegrationRule() const;
  };


  /// quad of order 1
  class FE_Quad1 : public T_ScalarFiniteElement2<FE_Quad1,ET_QUAD,4,1>
  {
  public:
    static IPDataArray ipdata;
		
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[2], TFA & shape) 
    {
      Tx x = hx[0];
      Tx y = hx[1];

      shape[0] = (1-x) * (1-y);
      shape[1] =    x  * (1-y);
      shape[2] =    x  *  y;
      shape[3] = (1-x) *  y;
    }
	  
    virtual const IntegrationRule & NodalIntegrationRule() const;
  }; 


  /// quad or order 2
  class FE_Quad2 : public T_ScalarFiniteElement2<FE_Quad2,ET_QUAD,9,2>
  {
  public:
    static IPDataArray ipdata;

    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[2], TFA & shape) 
    {
      Tx x = hx[0];
      Tx y = hx[1];

      Vec<3,Tx> px, py;
      px(0) = (1-x) * (1-2*x);
      px(1) = 4 * x * (1-x);
      px(2) = x * (2*x-1);
      py(0) = (1-y) * (1-2*y);
      py(1) = 4 * y * (1-y);
      py(2) = y * (2*y-1);
      
      int ii = 0;
      for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
	  shape[ii++] = px(i) * py(j);
    }
			  
    virtual const IntegrationRule & NodalIntegrationRule() const;
  }; 


  class FE_Quad2aniso :  public T_ScalarFiniteElement2<FE_Quad2aniso,ET_QUAD,6,2>
  {
  public:
    static IPDataArray ipdata;

    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[2], TFA & shape) 
    {
      Tx x = hx[0];
      Tx y = hx[1];

      shape[0] = (1-x)*(1-2*x) * (1-y);
      shape[1] = x*(2*x-1) * (1-y);
      shape[2] = x*(2*x-1) * y;
      shape[3] = (1-x)*(1-2*x) * y;
      shape[4] = 4*x*(1-x) * (1-y);
      shape[5] = 4*x*(1-x) * y;
    }
  }; 
  


  /* ***************************** Tet *************************************** */



  ///
  class FE_Tet0 : public T_ScalarFiniteElement2<FE_Tet0,ET_TET,1,0>
  {
  public:
    static IPDataArray ipdata;

    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[3], TFA & shape) 
    {
      shape[0] = 1;
    }
  
    virtual const IntegrationRule & NodalIntegrationRule() const;
  };


  ///
  class FE_Tet1 : public T_ScalarFiniteElement2<FE_Tet1,ET_TET,4,1>
  {
  public:
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

    virtual const IntegrationRule & NodalIntegrationRule() const;
    virtual void GetDofs (Array<Dof> & dofs) const;
  };


  ///
  class FE_Tet2 : public T_ScalarFiniteElement2<FE_Tet2,ET_TET,10,2>
  {
  public:
    static IPDataArray ipdata;

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
  class FE_Tet2HB : public T_ScalarFiniteElement2<FE_Tet2HB,ET_TET,10,2>
  {
  public:
    static IPDataArray ipdata;

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
