#ifndef FILE_SCALARFE
#define FILE_SCALARFE

/*********************************************************************/
/* File:   scalarfe.hpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

namespace ngfem
{


  /**
     Scalar finite element.
     Provides shape functions and derivaties.
     Values of shape functions and derivatives in integration points 
     might be stored as static data (IPData).
  */
  template <int D>
  class ScalarFiniteElement : public FiniteElement
  {
  
  protected:

    /// stored information in integration points
    class IPData
    {
    public:
      FlatVector<> shape;
      FlatMatrix<> dshape;
    };

    class IPDataArray
    {
    public:
      Array<IPData> data;
      int Size() const { return data.Size(); }
      const IPData & operator[] (int i) { return data[i]; }

      DynamicMem<double> block;
    };


    IPData * p_ipdata;

  public:
    virtual string ClassName(void) const {return "ScalarFiniteElement";}
    ///
    ScalarFiniteElement () { dimspace = D; p_ipdata = 0; }
    ///
    ScalarFiniteElement (ELEMENT_TYPE aeltype, 
			 int andof = 0, int aorder = 0)
      : FiniteElement (D, aeltype, andof, aorder) { p_ipdata = 0; /* block = 0; */ }
    ///
    virtual ~ScalarFiniteElement () { ; }
    ///
    virtual const IntegrationRule & NodalIntegrationRule() const;

    /**
       returns shape functions in point ip.
       returns stored values for valid ip.IPNr(), else computes values
    */
    const FlatVector<> GetShape (const IntegrationPoint & ip, 
				 LocalHeap & lh) const
    {
      if (p_ipdata && ip.IPNr() > -1) 
	return p_ipdata[ip.IPNr()].shape;
      else
	{
	  FlatVector<> shape(ndof, lh);
	  CalcShape (ip, shape);
	  return shape;
	}
    }

    /**
       returns derivatives in point ip.
       returns stored values for valid ip.IPNr(), else computes values
    */
    const FlatMatrix<> GetDShape (const IntegrationPoint & ip, LocalHeap & lh) const
    {
      if (p_ipdata && ip.IPNr() > -1)
	{
	  return p_ipdata[ip.IPNr()].dshape;
	}
      else
	{
	  FlatMatrix<> dshape(ndof, dimspace, lh);
	  CalcDShape (ip, dshape);
	  return dshape;
	}
    }


    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const = 0;
  
    /// compute dshape, matrix: ndof x spacedim
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<D> dshape) const;

    /// compute dshape, matrix: ndof x spacedim
    virtual void CalcMappedDShape (const SpecificIntegrationPoint<D,D> & sip, 
				   FlatMatrixFixWidth<D> dshape) const;

    virtual double
    Evaluate (const IntegrationPoint & ip, 
	      FlatVector<double> x, LocalHeap & lh) const
    {
      return InnerProduct (GetShape(ip, lh), x);
    }  


    /**
       returns second derivatives in point ip.
       returns stored values for valid ip.IPNr(), else computes values
    */
    const FlatMatrix<> GetDDShape (const IntegrationPoint & ip, LocalHeap & lh) const
    {
      FlatMatrix<> ddshape(ndof, dimspace*dimspace, lh);
      CalcDDShape (ip, ddshape);
      return ddshape;
    }

    /// compute dshape, matrix: ndof x (spacedim spacedim)
    virtual void CalcDDShape (const IntegrationPoint & ip, 
			      FlatMatrix<> ddshape) const;



    virtual void EvaluateShapeGrid (const IntegrationRuleTP<D> & ir,
				    const FlatVector<double> coefs,
				    FlatVector<double> gridvalues,
				    LocalHeap & lh) const;
				  
    virtual void EvaluateShapeGridTrans (const IntegrationRuleTP<D> & ir,
					 const FlatVector<double> gridvalues,
					 FlatVector<double> coefs,
					 LocalHeap & lh) const;
				  
    virtual void EvaluateDShapeGrid (const IntegrationRuleTP<D> & ir,
				     const FlatVector<double> coefs,
				     FlatMatrixFixWidth<D> gridvalues,
				     LocalHeap & lh) const;
				  
    virtual void EvaluateDShapeGridTrans (const IntegrationRuleTP<D> & ir,
					  const FlatMatrixFixWidth<D> gridvalues,
					  FlatVector<double> coefs,
					  LocalHeap & lh) const;
				  

  protected:
    ///
    void CalcIPData (ELEMENT_TYPE et,
		     IPDataArray & ipdata);
  };















  /**
     Base-element for template polymorphism.
     Barton and Nackman Trick
  */

  template <class FEL, int SDIM, int NDOF>
  class T_ScalarFiniteElement : public ScalarFiniteElement<SDIM>
  {

  public:

  protected:
    T_ScalarFiniteElement ()
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

    virtual ~T_ScalarFiniteElement() { ; }

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
      FEL::CalcShapeStat (ip, shape);
    }


    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<SDIM> dshape) const
    {
      FEL::CalcDShapeStat (ip, dshape);
    }



    template <class TP, class T>
    static void CalcDShapeStat (const TP & p, T & dshape)
    {

      int i;
      double eps = 2e-5;
      Vec<NDOF> shape1, shape2, shape3, shape4;

      IntegrationPoint ip;
      for (i = 0; i < SDIM; i++)
	ip(i) = p(i);

      for (i = 0; i < SDIM; i++)
	{
	  IntegrationPoint ip1 = ip;
	  IntegrationPoint ip2 = ip;
	  ip1(i) -= eps;
	  ip2(i) += eps;
	  FEL::CalcShapeStat (ip1, shape1);
	  FEL::CalcShapeStat (ip2, shape2);

	  ip1(i) -= eps;
	  ip2(i) += eps;
	  FEL::CalcShapeStat (ip1, shape3);
	  FEL::CalcShapeStat (ip2, shape4);

	  for (int j = 0; j < NDOF; j++)
	    dshape(j, i) = 
	      2/(3*eps) * (shape2(j) - shape1(j)) 
	      -1/(12*eps) * (shape4(j) - shape3(j));
	}    
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















  /* ********************************** Segm ********************************* */

  /*
 ///
 class FE_Segm0 : public T_ScalarFiniteElement<FE_Segm0,1,1>
 {
 public:
 enum { SDIM = 1 };
 enum { NDOF = 1 };
 enum { ORDER = 0 };
 enum { ELTYPE = ET_SEGM };

 static IPDataArray ipdata;

 template <class TP, class T>
 static void CalcShapeStat (TP & p, T & shape)
 {
 shape(0) = 1;
 }

 template <class TP, class T>
 static void CalcDShapeStat (TP & p, T & dshape)
 {
 dshape(0,0) = 0;
 }
 }; 

 ///
 class FE_Segm1 : public T_ScalarFiniteElement<FE_Segm1,1,2>
 {
 public:
 enum { SDIM = 1 };
 enum { NDOF = 2 };
 enum { ORDER = 1 };
 enum { ELTYPE = ET_SEGM };

 static IPDataArray ipdata;

 template <class TP, class T>
 static void CalcShapeStat (TP & p, T & shape)
 {
 double x = p(0);

 shape(0) = x;
 shape(1) = 1-x;
 }

 template <class TP, class T>
 static void CalcDShapeStat (TP & p, T & dshape)
 {
 dshape(0,0) = 1;
 dshape(1,0) = -1;
 }
 }; 







 ///
 class FE_Segm1L2 : public T_ScalarFiniteElement<FE_Segm1L2,1,2>
 {
 public:
 enum { SDIM = 1 };
 enum { NDOF = 2 };
 enum { ORDER = 1 };
 enum { ELTYPE = ET_SEGM };

 static IPDataArray ipdata;


 FE_Segm1L2 ();

 template <class TP, class T>
 static void CalcShapeStat (TP & p, T & shape)
 {
 double x = p(0);

 shape(0) = 1;
 shape(1) = 2*x-1;
 }

 template <class TP, class T>
 static void CalcDShapeStat (TP & p, T & dshape)
 {
 dshape(0,0) = 0;
 dshape(1,0) = 2;
 }
 }; 



 ///
 class FE_Segm2 : public T_ScalarFiniteElement<FE_Segm2,1,3>
 {

 public:
 enum { SDIM = 1 };
 enum { NDOF = 3 };
 enum { ORDER = 2 };
 enum { ELTYPE = ET_SEGM };

 static IPDataArray ipdata;

 ///
 template <class TP, class T>
 static void CalcShapeStat (TP & p, T & shape)
 {
 double x = p(0);

 shape(0) = 2*x*x - x;
 shape(1) = 2*x*x - 3*x + 1;  
 shape(2) = 4 * x * (1-x);
 }

  
 ///
 template <class TP, class T>
 static void CalcDShapeStat (TP & p, T & dshape)
 {
 double x = p(0);

 dshape(0, 0) = 4*x - 1;
 dshape(1, 0) = 4*x - 3;
 dshape(2, 0) = 4 - 8 * x;
 }


 }; 


 ///
 class FE_Segm2HB : public T_ScalarFiniteElement<FE_Segm2HB,1,3>
 {
 public:
 enum { SDIM = 1 };
 enum { NDOF = 3 };
 enum { ORDER = 2 };
 enum { ELTYPE = ET_SEGM };

 static IPDataArray ipdata;

 ///
 template <class TP, class T>
 static void CalcShapeStat (TP & p, T & shape)
 {
 double x = p(0);

 shape(0) = x;
 shape(1) = 1-x;
 shape(2) = 4 * x * (1-x);
 }
 }; 



 ///
 class FE_Segm2L2 : public T_ScalarFiniteElement<FE_Segm2L2,1,3>
 {
 public:
 enum { SDIM = 1 };
 enum { NDOF = 3 };
 enum { ORDER = 2 };
 enum { ELTYPE = ET_SEGM };

 static IPDataArray ipdata;

 FE_Segm2L2();
  
 ///
 template <class TP, class T>
 static void CalcShapeStat (TP & p, T & shape)
 {
 double x = p(0);

 shape(0) = 1;
 shape(1) = 2*x-1;
 shape(2) = (2*x-1)*(2*x-1)-1.0/3.0;
 }
 }; 
  */


  /// segment of fixed order
  template <int ORDER>
  class FE_TSegmL2 : public ScalarFiniteElement<1>
  {
    static IPDataArray ipdata;
  public:
    ///
    FE_TSegmL2();
    ///
    virtual ~FE_TSegmL2();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
    ///
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<1> shape) const;
  };



  /*
  /// Non-conforming finite elements
  class FE_NcSegm1 : public T_ScalarFiniteElement<FE_NcSegm1,1,1>
  {
  public:
    enum { SDIM = 1 };
    enum { NDOF = 1 };
    enum { ORDER = 1 };
    enum { ELTYPE = ET_SEGM };

    static IPDataArray ipdata;

    ///
    template <class TP, class T>
    static void CalcShapeStat (TP & p, T & shape)
    {
      shape(0) = 1;
    }
  
    ///
    template <class TP, class T>
    static void CalcDShapeStat (TP & p, T & dshape)
    {
      dshape(0,0) = 0;
    }

  }; 




  /// potential space for Nedelec IIb
  class FE_Segm3Pot : public ScalarFiniteElement<1>
  {
    ///
    static IPDataArray ipdata;
  public:
    ///
    FE_Segm3Pot();
    ///
    virtual ~FE_Segm3Pot();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
  }; 
  */









  /* ********************************* Trigs ******************************* */

  /*
 ///
 class FE_Trig0 : public ScalarFiniteElement<2>
 {
 ///
 static IPDataArray ipdata;
 public:
 ///
 FE_Trig0();
 ///
 virtual ~FE_Trig0();

 ///
 virtual void CalcShape (const IntegrationPoint & ip, 
 FlatVector<> shape) const;
 ///
 virtual void CalcDShape (const IntegrationPoint & ip, 
 FlatMatrixFixWidth<2> dshape) const;
 ///
 virtual const IntegrationRule & NodalIntegrationRule() const;
 };

 ///
 class FE_Trig1 : public T_ScalarFiniteElement<FE_Trig1,2,3>
 {
 public:
 enum { SDIM = 2 };
 enum { NDOF = 3 };
 enum { ORDER = 1 };
 enum { ELTYPE = ET_TRIG };

 static IPDataArray ipdata;

 template<typename Tx, typename TFA>  
 void T_CalcShape (Tx x[2], TFA & shape) const
 {
 shape[0] = x[0];
 shape[1] = x[1];      
 shape[2] = 1-x[0]-x[1];
 }


 ///
 template <class TP, class T>
 static void CalcShapeStat (TP & p, T & shape)
 {
 double x = p(0);
 double y = p(1);

 shape(0) = x;
 shape(1) = y;
 shape(2) = 1-x-y;
 }

 ///
 template <class TP, class T>
 static void CalcDShapeStat (TP & p, T & dshape)
 {
 dshape = 0;
 dshape(0,0) = 1;
 dshape(1,1) = 1;
 dshape(2,0) = -1;
 dshape(2,1) = -1;
 }

 ///
 virtual const IntegrationRule & NodalIntegrationRule() const;
 }; 


 ///
 class FE_Trig2 : public T_ScalarFiniteElement<FE_Trig2,2,6>
 {
 public:
 enum { SDIM = 2 };
 enum { NDOF = 6 };
 enum { ORDER = 2 };
 enum { ELTYPE = ET_TRIG };

 static IPDataArray ipdata;

 ///
 template <class TP, class T>
 static void CalcShapeStat (TP & p, T & shape)
 {
 double x = p(0);
 double y = p(1);
 double lam3 = 1-x-y;
    
 shape(0) = x * (2*x-1);
 shape(1) = y * (2*y-1);
 shape(2) = lam3 * (2*lam3-1);
 shape(3) = 4 * y * lam3;
 shape(4) = 4 * x * lam3;
 shape(5) = 4 * x * y;
 }

 virtual const IntegrationRule & NodalIntegrationRule() const;
 }; 
  */


  ///
  class FE_Trig2HB : public ScalarFiniteElement<2>
  {
    ///
    static IPDataArray ipdata;

  public:
    ///
    FE_Trig2HB();
    ///
    virtual ~FE_Trig2HB();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
  }; 



  /// potential space for Nedelec IIb
  class FE_Trig3Pot : public ScalarFiniteElement<2>
  {
    ///
    static IPDataArray ipdata;
  public:
    ///
    FE_Trig3Pot();
    ///
    virtual ~FE_Trig3Pot();

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
			  
  }; 



  class FE_NC_Trig1 : public T_ScalarFiniteElement<FE_NC_Trig1,2,3>
  {
  public:
    enum { SDIM = 2 };
    enum { NDOF = 3 };
    enum { ORDER = 1 };
    enum { ELTYPE = ET_TRIG };

    static IPDataArray ipdata;

    ///
    template <class TP, class T>
    static void CalcShapeStat (TP & p, T & shape)
    {
      double x = p(0);
      double y = p(1);

      shape(0) = 1-2*y;
      shape(1) = 1-2*x;
      shape(2) = 2*(x+y)-1;
    }

    ///
    template <class TP, class T>
    static void CalcDShapeStat (TP & p, T & dshape)
    {
      // double x = p(0);
      // double y = p(1);

      dshape = 0;
      dshape(0,1) = -2;
      dshape(1,0) = -2;
      dshape(2,0) = 1;
      dshape(2,1) = 1;
    }
			   
    ///
    virtual const IntegrationRule & NodalIntegrationRule() const;
  }; 




  /* ***************************** Tet *************************************** */

  /*
 ///
 class FE_Tet0 : public T_ScalarFiniteElement<FE_Tet0,3,1>
 {
 public:
 enum { SDIM = 3 };
 enum { NDOF = 1 };
 enum { ORDER = 0 };
 enum { ELTYPE = ET_TET };

 static IPDataArray ipdata;

 ///
 template <class TP, class T>
 static void CalcShapeStat (TP & p, T & shape)
 {
 shape(0) = 1;
 }
  
 /// 
 template <class TP, class T>
 static void CalcDShapeStat (TP & p, T & dshape)
 {
 dshape = 0;
 }
  
 ///
 virtual const IntegrationRule & NodalIntegrationRule() const;
 };


 ///
 class FE_Tet1 : public T_ScalarFiniteElement<FE_Tet1,3,4>
 {
 public:
 enum { SDIM = 3 };
 enum { NDOF = 4 };
 enum { ORDER = 1 };
 enum { ELTYPE = ET_TET };

 static IPDataArray ipdata;

 ///
 template <class TP, class T>
 static void CalcShapeStat (TP & p, T & shape)
 {
 double x = p(0);
 double y = p(1);
 double z = p(2);

 shape(0) = x;
 shape(1) = y;
 shape(2) = z;
 shape(3) = 1-x-y-z;
 }
  
 ///
 template <class TP, class MAT>
 static void CalcDShapeStat (TP & p, MAT & dshape)
 {
 dshape = 0;
 dshape(0,0) = 1;
 dshape(1,1) = 1;
 dshape(2,2) = 1;
 dshape(3,0) = -1;
 dshape(3,1) = -1;
 dshape(3,2) = -1;
 }


 ///
 virtual const IntegrationRule & NodalIntegrationRule() const;
 virtual void GetDofs (Array<Dof> & dofs) const;
 };



 ///
 class FE_Tet2 : public T_ScalarFiniteElement<FE_Tet2,3,10>
 {
 public:
 enum { SDIM = 3 };
 enum { NDOF = 10 };
 enum { ORDER = 2 };
 enum { ELTYPE = ET_TET };

 static IPDataArray ipdata;

 ///
 template <class TP, class T>
 static void CalcShapeStat (TP & p, T & shape)
 {
 double x = p(0);
 double y = p(1);
 double z = p(2);
 double lam4 = 1 - x - y - z;
    
 shape(0) = 2 * x * x - x;  
 shape(1) = 2 * y * y - y;
 shape(2) = 2 * z * z - z;
 shape(3) = 2 * lam4 * lam4 - lam4;

 shape(4) = 4 * x * y;
 shape(5) = 4 * x * z;
 shape(6) = 4 * x * lam4;
 shape(7) = 4 * y * z;
 shape(8) = 4 * y * lam4;
 shape(9) = 4 * z * lam4;
 }
 };


 ///
 class FE_Tet2HB : public ScalarFiniteElement<3>
 {
 ///
 static IPDataArray ipdata;
 public:
 ///
 FE_Tet2HB();
 ///
 virtual ~FE_Tet2HB();
 ///
 virtual void CalcShape (const IntegrationPoint & ip, 
 FlatVector<> shape) const;
 };
  */





  /// tet of order 3 (obsolete)
  class FE_Tet3Pot : public ScalarFiniteElement<3>
  {
    ///
    static IPDataArray ipdata;
  public:
    ///
    FE_Tet3Pot();
    ///
    virtual ~FE_Tet3Pot();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
  }; 





  /* ***************************** Quads ********************************* */

  /*
 /// quad of order 0
 class FE_Quad0 : public ScalarFiniteElement<2>
 {
 static IPDataArray ipdata;
 public:
 FE_Quad0();
 virtual ~FE_Quad0();

 virtual void CalcShape (const IntegrationPoint & ip, 
 FlatVector<> shape) const;
			  
 virtual void CalcDShape (const IntegrationPoint & ip, 
 FlatMatrixFixWidth<2> dshape) const;
			   
 virtual const IntegrationRule & NodalIntegrationRule() const;
 };


 /// quad of order 1
 class FE_Quad1 : public ScalarFiniteElement<2>
 {
 static IPDataArray ipdata;

 public:
 FE_Quad1();
 virtual ~FE_Quad1();
 virtual void CalcShape (const IntegrationPoint & ip, 
 FlatVector<> shape) const;
 virtual void CalcDShape (const IntegrationPoint & ip, 
 FlatMatrixFixWidth<2> dshape) const;
			  
 virtual const IntegrationRule & NodalIntegrationRule() const;
 }; 


 /// quad or order 2
 class FE_Quad2 : public ScalarFiniteElement<2>
 {
 static IPDataArray ipdata;

 public:
 FE_Quad2();
 virtual ~FE_Quad2();
 virtual void CalcShape (const IntegrationPoint & ip, 
 FlatVector<> shape) const;
 virtual void CalcDShape (const IntegrationPoint & ip, 
 FlatMatrixFixWidth<2> dshape) const;
			  
 virtual const IntegrationRule & NodalIntegrationRule() const;
 }; 
  */

  /// quad of order 3
  class FE_Quad3 : public ScalarFiniteElement<2>
  {
    static IPDataArray ipdata;

  public:
    FE_Quad3();
    virtual ~FE_Quad3();
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<2> dshape) const;
  }; 


  /*
 /// second order x, first order y
 class FE_Quad2aniso : public ScalarFiniteElement<2>
 {
 static IPDataArray ipdata;

 public:
 FE_Quad2aniso();
 virtual ~FE_Quad2aniso();
 virtual void CalcShape (const IntegrationPoint & ip, 
 FlatVector<> shape) const;
 }; w
  */


  /* **************************** Pyramid Elements *********************** */


  /// pyramid of order 0
  class FE_Pyramid0 : public ScalarFiniteElement<3>
  {
    ///
    static IPDataArray ipdata;

  public:

    ///
    FE_Pyramid0();
    ///
    virtual ~FE_Pyramid0();

    ///
    virtual int SpatialDim () const { return 3; }
    ///
    virtual int GetNDof () const { return 1; }
    ///
    virtual int Order () const { return 0; }
    ///
    virtual ELEMENT_TYPE ElementType() const { return ET_PYRAMID; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
			  
    ///
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> dshape) const;

    ///
    virtual const IntegrationRule & NodalIntegrationRule() const;

  };




  /// pyramid of order 1
  class FE_Pyramid1 : public ScalarFiniteElement<3>
  {
    ///
    static IPDataArray ipdata;

  public:

    ///
    FE_Pyramid1();
    ///
    virtual ~FE_Pyramid1();

    ///
    virtual int SpatialDim () const { return 3; }
    ///
    virtual int GetNDof () const { return 5; }
    ///
    virtual int Order () const { return 3; }
    ///
    virtual ELEMENT_TYPE ElementType() const { return ET_PYRAMID; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
			  
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> dshape) const;
			  

    ///
    virtual const IntegrationRule & NodalIntegrationRule() const;
  };



  /// pyramid of order 2
  class FE_Pyramid2 : public ScalarFiniteElement<3>
  {
    ///
    static IPDataArray ipdata;

  public:

    ///
    FE_Pyramid2();
    ///
    virtual ~FE_Pyramid2();

    ///
    virtual int SpatialDim () const { return 3; }
    ///
    virtual int GetNDof () const { return 13; }
    ///
    virtual int Order () const { return 4; }
    ///
    virtual ELEMENT_TYPE ElementType() const { return ET_PYRAMID; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
			  
  };





  /* ******************************* Prism Elements ********************* */



  /// prism of order 0
  class FE_Prism0 : public ScalarFiniteElement<3>
  {
    ///
    static IPDataArray ipdata;

  public:

    ///
    FE_Prism0();
    ///
    virtual ~FE_Prism0();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;

    ///
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> dshape) const;

    ///
    virtual const IntegrationRule & NodalIntegrationRule() const;
  };


  /// prism of order 1
  class FE_Prism1 : public ScalarFiniteElement<3>
  {
    ///
    static IPDataArray ipdata;

  public:

    ///
    FE_Prism1();
    ///
    virtual ~FE_Prism1();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
  
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> dshape) const;
			  
    virtual const IntegrationRule & NodalIntegrationRule() const;
  };




  /// prism of order 2
  class FE_Prism2 : public ScalarFiniteElement<3>
  {
    ///
    static IPDataArray ipdata;

  public:

    ///
    FE_Prism2();
    ///
    virtual ~FE_Prism2();
    /*
   ///
   virtual int SpatialDim () const { return 3; }
   ///
   virtual int GetNDof () const { return 18; }
   ///
   virtual int Order () const { return 3; }
   ///
   virtual ELEMENT_TYPE ElementType() const { return ET_PRISM; }
    */
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
			  
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> dshape) const;
  };





  /// in plane second order
  class FE_Prism2aniso : public ScalarFiniteElement<3>
  {
    ///
    static IPDataArray ipdata;

  public:

    ///
    FE_Prism2aniso();
    ///
    virtual ~FE_Prism2aniso();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
			  
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> dshape) const;
  };



  /// in plane second order
  class FE_Prism2HBaniso : public ScalarFiniteElement<3>
  {
    ///
    static IPDataArray ipdata;

  public:

    ///
    FE_Prism2HBaniso();
    ///
    virtual ~FE_Prism2HBaniso();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
			  

  };






  /// in plane third order
  class FE_Prism3aniso : public ScalarFiniteElement<3>
  {
    ///
    static IPDataArray ipdata;

  public:

    ///
    FE_Prism3aniso();
    ///
    virtual ~FE_Prism3aniso();

    ///
    virtual int SpatialDim () const { return 3; }
    ///
    virtual int GetNDof () const { return 20; }
    ///
    virtual int Order () const { return 3; }
    ///
    virtual ELEMENT_TYPE ElementType() const { return ET_PRISM; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
			  
  };





  /* **************************** Hex elements ************************* */



  ///
  class FE_Hex0 : public ScalarFiniteElement<3>
  {
    ///
    static IPDataArray ipdata;
  public:
    ///
    FE_Hex0();
    ///
    virtual ~FE_Hex0();

    ///
    virtual int SpatialDim () const { return 3; }
    ///
    virtual int GetNDof () const { return 1; }
    ///
    virtual int Order () const { return 1; }
    ///
    virtual ELEMENT_TYPE ElementType() const { return ET_HEX; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
  };


  ///
  class FE_Hex1 : public ScalarFiniteElement<3>
  {
    ///
    static IPDataArray ipdata;
  public:
    ///
    FE_Hex1();
    ///
    virtual ~FE_Hex1();

    ///
    virtual int SpatialDim () const { return 3; }
    ///
    virtual int GetNDof () const { return 8; }
    ///
    virtual int Order () const { return 1; }
    ///
    virtual ELEMENT_TYPE ElementType() const { return ET_HEX; }
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
			  
  };










  ///
  class FE_NcTrig1 : public ScalarFiniteElement<2>
  {
    ///
    static IPDataArray ipdata;

  public:
    ///
    FE_NcTrig1();
    ///
    virtual ~FE_NcTrig1();

    ///
    virtual int SpatialDim () const { return 2; }
    ///
    virtual int GetNDof () const { return 3; }
    ///
    virtual int Order () const { return 1; }
    ///
    virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
			  
  }; 





  ///
  class FE_NcTet1 : public ScalarFiniteElement<3>
  {
    ///
    static IPDataArray ipdata;

  public:
    ///
    FE_NcTet1();
    ///
    virtual ~FE_NcTet1();

    ///
    virtual int SpatialDim () const { return 3; }
    ///
    virtual int GetNDof () const { return 4; }
    ///
    virtual int Order () const { return 1; }
    ///
    virtual ELEMENT_TYPE ElementType() const { return ET_TET; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
			  
  }; 


  /* ************************** */
  /*    Dummy Elements          */
  /* ************************** */

  class FE_SegmDummy : public T_ScalarFiniteElement<FE_SegmDummy,1,0>
  {
  public:
    enum { SDIM = 1 };
    enum { NDOF = 0 };
    enum { ORDER = 0 };
    enum { ELTYPE = ET_SEGM };

    static IPDataArray ipdata;

    template <class TP, class T>
    static void CalcShapeStat (TP & p, T & shape)
    {
      ;
    }

    template <class TP, class T>
    static void CalcDShapeStat (TP & p, T & dshape)
    {
      ;
    }
  }; 

  class FE_TrigDummy : public T_ScalarFiniteElement<FE_TrigDummy,2,0>
  {
  public:
    enum { SDIM = 2 };
    enum { NDOF = 0 };
    enum { ORDER = 0 };
    enum { ELTYPE = ET_TRIG };

    static IPDataArray ipdata;

    template <class TP, class T>
    static void CalcShapeStat (TP & p, T & shape)
    {
      ;
    }

    template <class TP, class T>
    static void CalcDShapeStat (TP & p, T & dshape)
    {
      ;
    }
  }; 


  class FE_QuadDummy : public T_ScalarFiniteElement<FE_QuadDummy,2,0>
  {
  public:
    enum { SDIM = 2 };
    enum { NDOF = 0 };
    enum { ORDER = 0 };
    enum { ELTYPE = ET_QUAD };

    static IPDataArray ipdata;

    template <class TP, class T>
    static void CalcShapeStat (TP & p, T & shape)
    {
      ;
    }

    template <class TP, class T>
    static void CalcDShapeStat (TP & p, T & dshape)
    {
      ;
    }
  }; 


  class FE_TetDummy : public T_ScalarFiniteElement<FE_TetDummy,3,0>
  {
  public:
    enum { SDIM = 3 };
    enum { NDOF = 0 };
    enum { ORDER = 0 };
    enum { ELTYPE = ET_TET };

    static IPDataArray ipdata;

    template <class TP, class T>
    static void CalcShapeStat (TP & p, T & shape)
    {
      ;
    }

    template <class TP, class T>
    static void CalcDShapeStat (TP & p, T & dshape)
    {
      ;
    }
  }; 


  class FE_HexDummy : public T_ScalarFiniteElement<FE_HexDummy,3,0>
  {
  public:
    enum { SDIM = 3 };
    enum { NDOF = 0 };
    enum { ORDER = 0 };
    enum { ELTYPE = ET_HEX };

    static IPDataArray ipdata;

    template <class TP, class T>
    static void CalcShapeStat (TP & p, T & shape)
    {
      ;
    }

    template <class TP, class T>
    static void CalcDShapeStat (TP & p, T & dshape)
    {
      ;
    }
  }; 

  class FE_PrismDummy : public T_ScalarFiniteElement<FE_PrismDummy,3,0>
  {
  public:
    enum { SDIM = 3 };
    enum { NDOF = 0 };
    enum { ORDER = 0 };
    enum { ELTYPE = ET_PRISM };

    static IPDataArray ipdata;

    template <class TP, class T>
    static void CalcShapeStat (TP & p, T & shape)
    {
      ;
    }

    template <class TP, class T>
    static void CalcDShapeStat (TP & p, T & dshape)
    {
      ;
    }
  }; 


  class FE_PyramidDummy : public T_ScalarFiniteElement<FE_PyramidDummy,3,0>
  {
  public:
    enum { SDIM = 3 };
    enum { NDOF = 0 };
    enum { ORDER = 0 };
    enum { ELTYPE = ET_PYRAMID };

    static IPDataArray ipdata;

    template <class TP, class T>
    static void CalcShapeStat (TP & p, T & shape)
    {
      ;
    }

    template <class TP, class T>
    static void CalcDShapeStat (TP & p, T & dshape)
    {
      ;
    }
  }; 


}

#endif
