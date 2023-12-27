#ifndef FILE_H1LOFE
#define FILE_H1LOFE

/*********************************************************************/
/* File:   h1lofe.hpp                                                */
/* Author: Start                                                     */
/* Date:   29. Jun. 2009                                              */
/*********************************************************************/

#include "tscalarfe.hpp"
#include "coefficient.hpp"

#ifdef FILE_H1LOFE_CPP
#define H1LOFE_EXTERN
#include <tscalarfe_impl.hpp>
#else
#define H1LOFE_EXTERN extern
#endif


namespace ngfem
{


  template <ELEMENT_TYPE ET, int ORDER>
  class ScalarFE : public T_ScalarFiniteElement<ScalarFE<ET,ORDER>,ET>
  {
    typedef T_ScalarFiniteElement<ScalarFE<ET,ORDER>,ET> BASE;
  public:
    INLINE ScalarFE ()
    {
      this->ndof = ET_trait<ET>::PolDimension(ORDER);
      this->order = ORDER; 
    }

    static constexpr int DIM = ngfem::Dim(ET);

    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<DIM,Tx> ip, TFA & shape);

    template<typename Tx, typename TFA>  
    INLINE void T_CalcDualShape (const TIP<DIM,Tx> ip, TFA & shape) const
    {
      if (ORDER == 0)
        {
          if (ip.vb == VOL)
            shape[0] = 1;
          else
            shape[0] = 0;
        }
      else if (ORDER == 1)
        {
          if (int(ip.vb) == ET_trait<ET>::DIM)
            shape[ip.facetnr] = 1;
        }
      else
        throw Exception (string("CalcDualShape not overloaded for element ") + typeid(*this).name());
    }

    virtual bool DualityMassDiagonal () const override { return true; }
    virtual bool GetDiagDualityMassInverse (FlatVector<> diag) const override
    {
      diag = 1;
      return true;
    }

    
    virtual void Interpolate (const ElementTransformation & trafo, 
                              const class CoefficientFunction & func, SliceMatrix<> coefs,
                              LocalHeap & lh) const override
    {
      if (auto ipts = GetNodalPoints(); ipts.Size())
        {
          HeapReset hr(lh);
          IntegrationRule ir(ipts.Size(), ipts.Data());
          auto & mir = trafo(ir, lh);
          func.Evaluate (mir, coefs);
        }
      else
        BASE::Interpolate (trafo, func, coefs, lh);
    }

    FlatArray<IntegrationPoint> GetNodalPoints() const
    {
      return { 0, nullptr };
    }
    
    
    void CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
    {
      if (ORDER == 0)
        {
          if (mip.IP().VB() == VOL)
            shape[0] = 1;
          else
            shape[0] = 0;
        }
      else if (ORDER == 1)
        {
          shape = 0.0;
          if (int(mip.IP().VB()) == ET_trait<ET>::DIM)
            shape[mip.IP().FacetNr()] = 1;
        }
      else
        throw Exception (string("CalcDualShape not overloaded for element ") + typeid(*this).name());
    }
  };
  
  /*
  class FE_Point : public T_ScalarFiniteElementFO<FE_Point,ET_POINT,1,0>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (Tx x[], TFA & shape) 
    {
      shape[0] = 1;
    }
  }; 
  */

  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_POINT,0> :: T_CalcShape (TIP<0,Tx> ip, TFA & shape) 
  {
    shape[0] = Tx(1.0);
  }
  using FE_Point = ScalarFE<ET_POINT,0>;

  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_POINT,1> :: T_CalcShape (TIP<0,Tx> ip, TFA & shape) 
  {
    shape[0] = Tx(1.0);
  }

  /*
  class FE_Segm0 : public T_ScalarFiniteElementFO<FE_Segm0,ET_SEGM,1,0>
  {
  public:   
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (Tx x[1], TFA & shape) 
    {
      shape[0] = 1;
    }
  }; 
  */
  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_SEGM,0> :: T_CalcShape (TIP<1,Tx> ip, TFA & shape) 
  {
    shape[0] = Tx(1.0);
  }
  using FE_Segm0 = ScalarFE<ET_SEGM,0>;


  /*
  class FE_Segm1 : public T_ScalarFiniteElementFO<FE_Segm1,ET_SEGM,2,1>
  {
  public:
        NGS_DLL_HEADER FE_Segm1() { ; }
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (Tx x[1], TFA & shape) 
    {
      shape[0] = x[0];
      shape[1] = 1-x[0];
    }
  }; 
  */

  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_SEGM,1> :: T_CalcShape (TIP<1,Tx> ip, TFA & shape) 
  {
    shape[0] = ip.x;
    shape[1] = 1-ip.x;
  }
  using FE_Segm1 = ScalarFE<ET_SEGM,1>;


  /*
  class FE_Segm2 : public T_ScalarFiniteElementFO<FE_Segm2,ET_SEGM,3,2>
  {
  public:
        NGS_DLL_HEADER FE_Segm2() { ; }
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (Tx hx[1], TFA & shape) 
    {
      Tx x = hx[0];

      shape[0] = 2*x*x - x;
      shape[1] = 2*x*x - 3*x + 1;  
      shape[2] = 4 * x * (1-x);
    }
  }; 
  */
  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_SEGM,2> :: T_CalcShape (TIP<1,Tx> ip, TFA & shape) 
  {
    Tx x = ip.x;
    shape[0] = 2*x*x - x;
    shape[1] = 2*x*x - 3*x + 1;  
    shape[2] = 4 * x * (1-x);
  }
  using FE_Segm2 = ScalarFE<ET_SEGM,2>;


  ///
  class FE_Segm2HB : public T_ScalarFiniteElementFO<FE_Segm2HB,ET_SEGM,3,2>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<1,Tx> ip, TFA & shape) 
    {
      Tx x = ip.x;

      shape[0] = x;
      shape[1] = 1-x;
      shape[2] = 4 * x * (1-x);
    }
  }; 


  ///
  class FE_Segm1L2 : public T_ScalarFiniteElementFO<FE_Segm1L2,ET_SEGM,2,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<1,Tx> ip, TFA & shape) 
    {
      shape[0] = Tx(1.0);
      shape[1] = 2*ip.x-1;
    }
  }; 

  ///
  class FE_Segm2L2 :public T_ScalarFiniteElementFO<FE_Segm2L2,ET_SEGM,3,2>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<1,Tx> ip, TFA & shape) 
    {
      Tx x = ip.x;

      shape[0] = Tx(1.0);
      shape[1] = 2*x-1;
      shape[2] = (2*x-1)*(2*x-1)-1.0/3.0;
    }
  }; 

  ///
  class FE_NcSegm1 :public T_ScalarFiniteElementFO<FE_NcSegm1,ET_SEGM,1,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<1,Tx> ip, TFA & shape) 
    {
      shape[0] = Tx(1.0);
    }
  }; 

  /// potential space for Nedelec IIb
  class FE_Segm3Pot :public T_ScalarFiniteElementFO<FE_Segm3Pot,ET_SEGM,4,3>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<1,Tx> ip, TFA & shape) 
    {
      Tx x = ip.x;
      Tx lam2 = 1-x;

      shape[0] = x;
      shape[1] = lam2;
      shape[2] = 3 * x * lam2 * (lam2+x);
      shape[3] = 7.5 * x * lam2 * (x-lam2);
    }
  }; 


  /// segment of fixed order
  template <int ORDER>
  class FE_TSegmL2 : public T_ScalarFiniteElementFO<FE_TSegmL2<ORDER>, ET_SEGM, ORDER+1, ORDER>
  {
  public:
    template<typename Tx, typename TFA>  
      static INLINE void T_CalcShape (TIP<1,Tx> ip, TFA & shape) 
    {
      Tx x = ip.x;

      if (ORDER >= 4)
        LegendrePolynomial(ORDER, 2*x-1, shape);
      else
        {
          if (ORDER >= 0) shape[0] = Tx(1.0);
          if (ORDER >= 1) shape[1] = 2*x-1;
          if (ORDER >= 2) shape[2] = (2*x-1)*(2*x-1)-1.0/3.0;
          if (ORDER >= 3) shape[3] = (2*x-1)*(2*x-1)*(2*x-1);
        }
    }
  };






  /* ***************************** Trig Elements *************************************** */

  /*
  ///
  class FE_Trig0 : public T_ScalarFiniteElement<FE_Trig0,ET_TRIG>
  {
  public:
    FE_Trig0() { ndof = 1; order = 0; }
    template<typename Tx, typename TFA>  
    INLINE void T_CalcShape (Tx x[2], TFA & shape) const
    {
      shape[0] = 1;
    }
  }; 
  */

  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_TRIG,0> :: T_CalcShape (TIP<2,Tx> ip, TFA & shape) 
  {
    shape[0] = Tx(1.0);
  }
  using FE_Trig0 = ScalarFE<ET_TRIG,0>;

  /*
  ///
  class FE_Trig1 : public T_ScalarFiniteElement<FE_Trig1,ET_TRIG>
  {
  public:
    FE_Trig1() { ndof = 3; order = 1;}
    template<typename Tx, typename TFA>  
    INLINE void T_CalcShape (Tx x[2], TFA & shape) const
    {
      shape[0] = x[0];
      shape[1] = x[1];      
      shape[2] = 1-x[0]-x[1];
    }
  }; 
  */
  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_TRIG,1> :: T_CalcShape (TIP<2,Tx> ip, TFA & shape) 
  {
    Tx x = ip.x;
    Tx y = ip.y;

    shape[0] = x;
    shape[1] = y;      
    shape[2] = 1-x-y;
  }
  using FE_Trig1 = ScalarFE<ET_TRIG,1>;

  /*
  class FE_Trig2 : public T_ScalarFiniteElement<FE_Trig2,ET_TRIG>
  {
  public:
    FE_Trig2() { ndof = 6; order = 2;}
    template<typename Tx, typename TFA>  
    INLINE void T_CalcShape (Tx hx[2], TFA & shape) const
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
  }; 
  */

  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_TRIG,2> :: T_CalcShape (TIP<2,Tx> ip, TFA & shape) 
  {
      Tx x = ip.x;
      Tx y = ip.y;
      Tx lam3 = 1-x-y;
    
      shape[0] = x * (2*x-1);
      shape[1] = y * (2*y-1);
      shape[2] = lam3 * (2*lam3-1);
      shape[3] = 4 * y * lam3;
      shape[4] = 4 * x * lam3;
      shape[5] = 4 * x * y;
  }
  using FE_Trig2 = ScalarFE<ET_TRIG,2>;


  class FE_Trig2HB : public T_ScalarFiniteElement<FE_Trig2HB,ET_TRIG>
  {
  public:
    FE_Trig2HB() { ndof = 6; order = 2;}
    template<typename Tx, typename TFA>  
    INLINE void T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
    {
      Tx x = ip.x;
      Tx y = ip.y;
      Tx lam3 = 1-x-y;
    
      shape[0] = x;
      shape[1] = y;
      shape[2] = lam3;
      shape[3] = 4 * y * lam3;
      shape[4] = 4 * x * lam3;
      shape[5] = 4 * x * y;
    }
  }; 



  class FE_NcTrig1 : public T_ScalarFiniteElement<FE_NcTrig1,ET_TRIG>
  {
  public:
    FE_NcTrig1() { ndof = 3; order = 1;}
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<2,Tx> ip, TFA & shape) 
    {
      Tx x = ip.x;
      Tx y = ip.y;
      Tx lam3 = 1-x-y;
    
      shape[0] = 1-2*y;
      shape[1] = 1-2*x;
      shape[2] = 1-2*lam3;
    }
  };




  /* ***************************** Quad *************************************** */

  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_QUAD,0> :: T_CalcShape (TIP<2,Tx> ip, TFA & shape) 
  {
    shape[0] = Tx(1.0);
  }

  /*
  class FE_Quad0 : public T_ScalarFiniteElementFO<FE_Quad0,ET_QUAD,1,0>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (Tx hx[2], TFA & shape) 
    {
      shape[0] = 1.0;
    }
  };
  */
  /*
  /// quad of order 1
  class FE_Quad1 : public T_ScalarFiniteElementFO<FE_Quad1,ET_QUAD,4,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (Tx hx[2], TFA & shape) 
    {
      Tx x = hx[0];
      Tx y = hx[1];

      shape[0] = (1-x) * (1-y);
      shape[1] =    x  * (1-y);
      shape[2] =    x  *  y;
      shape[3] = (1-x) *  y;
    }
  }; 
  */

  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_QUAD,1> :: T_CalcShape (TIP<2,Tx> ip, TFA & shape) 
  {
    Tx x = ip.x;
    Tx y = ip.y;
    
    shape[0] = (1-x) * (1-y);
    shape[1] =    x  * (1-y);
    shape[2] =    x  *  y;
    shape[3] = (1-x) *  y;
  }


  /// quad or order 2
  class FE_Quad2 : public T_ScalarFiniteElementFO<FE_Quad2,ET_QUAD,9,2>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<2,Tx> ip, TFA & shape) 
    {
      Tx x = ip.x;
      Tx y = ip.y;

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
  };
  
  /// quad or order 2
  class FE_Quad2Serendipity : public T_ScalarFiniteElementFO<FE_Quad2Serendipity,ET_QUAD,8,2>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<2,Tx> ip, TFA & shape) 
    {
      Tx x = ip.x;
      Tx y = ip.y;
      
      shape[0] = (1-x)*(1-y) - 2 * x*(1-x) * (1-y) - 2 * y*(1-y)*(1-x);
      shape[1] =     x*(1-y) - 2 * x*(1-x) * (1-y) - 2 * y*(1-y)*x;
      shape[2] =     x*    y - 2 * x*(1-x) * y     - 2 * y*(1-y)*x;
      shape[3] = (1-x)*    y - 2 * x*(1-x) * y     - 2 * y*(1-y)*(1-x);

      shape[4] = 4 * x*(1-x) * (1-y);
      shape[5] = 4 * x*(1-x) * y;
      shape[6] = 4 * y*(1-y) * (1-x);
      shape[7] = 4 * y*(1-y) * x;
    }
  }; 


  class FE_Quad2aniso :  public T_ScalarFiniteElementFO<FE_Quad2aniso,ET_QUAD,6,2>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<2,Tx> ip, TFA & shape) 
    {
      Tx x = ip.x;
      Tx y = ip.y;

      shape[0] = (1-x)*(1-2*x) * (1-y);
      shape[1] = x*(2*x-1) * (1-y);
      shape[2] = x*(2*x-1) * y;
      shape[3] = (1-x)*(1-2*x) * y;
      shape[4] = 4*x*(1-x) * (1-y);
      shape[5] = 4*x*(1-x) * y;
    }
  }; 
  


  /* ***************************** Tet *************************************** */


  /*
  ///
  class FE_Tet0 : public T_ScalarFiniteElementFO<FE_Tet0,ET_TET,1,0>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (Tx hx[3], TFA & shape)
    {
      shape[0] = 1;
    }
  };
  */
  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_TET,0> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
  {
    shape[0] = Tx(1.0);
  }

  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_TET,1> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
  {
    Tx x = ip.x;
    Tx y = ip.y;
    Tx z = ip.z;
    
    shape[0] = x;
    shape[1] = y;
    shape[2] = z;
    shape[3] = 1-x-y-z;
  }

  /*
  class FE_Tet1 : public T_ScalarFiniteElementFO<FE_Tet1,ET_TET,4,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (Tx hx[3], TFA & shape) 
    {
      Tx x = hx[0];
      Tx y = hx[1];
      Tx z = hx[2];

      shape[0] = x;
      shape[1] = y;
      shape[2] = z;
      shape[3] = 1-x-y-z;
    }
  };
  */

  ///
  // class FE_Tet2 : public T_ScalarFiniteElementFO<FE_Tet2,ET_TET,10,2>
  // public:
  //  template<typename Tx, typename TFA>  
  //  static INLINE void T_CalcShape (TIP<3,Tx> ip, TFA & shape)
  
  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_TET,2> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
  {
      Tx x = ip.x;
      Tx y = ip.y;
      Tx z = ip.z;
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
  };



  ///
  class FE_Tet2HB : public T_ScalarFiniteElementFO<FE_Tet2HB,ET_TET,10,2>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
    {
      Tx x = ip.x;
      Tx y = ip.y;
      Tx z = ip.z;
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



  class FE_NcTet1 : public T_ScalarFiniteElementFO<FE_NcTet1,ET_TET,4,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
    {
      Tx x = ip.x;
      Tx y = ip.y;
      Tx z = ip.z;
      Tx lam4 = 1-x-y-z;
    
      shape[0] = 1-3*x;
      shape[1] = 1-3*y;
      shape[2] = 1-3*z;
      shape[3] = 1-3*lam4;
    }
  };


  /* ***************************** Prism *********************************** */

  ///
  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_PRISM,0> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
  {
    shape[0] = Tx(1);
  }
  using FE_Prism0 = ScalarFE<ET_PRISM,0>;
  
  ///
  /*
  class FE_Prism1 : public T_ScalarFiniteElementFO<FE_Prism1,ET_PRISM,6,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
    {
      Tx x = ip.x;
      Tx y = ip.y;
      Tx z = ip.z;
      
      shape[0] = x * (1-z);
      shape[1] = y * (1-z);
      shape[2] = (1-x-y) * (1-z);
      shape[3] = x * z;
      shape[4] = y * z;
      shape[5] = (1-x-y) * z;
    }
  };
  */

  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_PRISM,1> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
  {
    Tx x = ip.x;
    Tx y = ip.y;
    Tx z = ip.z;
    
    shape[0] = x * (1-z);
    shape[1] = y * (1-z);
    shape[2] = (1-x-y) * (1-z);
    shape[3] = x * z;
    shape[4] = y * z;
    shape[5] = (1-x-y) * z;
  }
  using FE_Prism1 = ScalarFE<ET_PRISM,1>;

  ///
  class FE_Prism2 : public T_ScalarFiniteElementFO<FE_Prism2,ET_PRISM,18,2>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
    {
      Tx x = ip.x;
      Tx y = ip.y;
      Tx z = ip.z;
      

      shape[0] = x * (1-z);
      shape[1] = y * (1-z);
      shape[2] = (1-x-y) * (1-z);
      shape[3] = x * z;
      shape[4] = y * z;
      shape[5] = (1-x-y) * z;
      
      shape[6] = 4 * x * (1-x-y) * (1-z);
      shape[7] = 4 * x * y       * (1-z);
      shape[8] = 4 * y * (1-x-y) * (1-z);
      shape[9] = 4 * x * (1-x-y) * z;
      shape[10] = 4 * x * y       * z;
      shape[11] = 4 * y * (1-x-y) * z;
      
      shape[12] = x * (1-z) * z;
      shape[13] = y * (1-z) * z;
      shape[14] = (1-x-y) * (1-z) * z;
      shape[15] = 4 * x * (1-x-y) * (1-z) * z;
      shape[16] = 4 * x * y       * (1-z) * z;
      shape[17] = 4 * y * (1-x-y) * (1-z) * z;
    }
  };



  /// prism, order P2 in plane and P1 in z direction
  class FE_Prism2aniso : public T_ScalarFiniteElementFO<FE_Prism2aniso,ET_PRISM,12,2>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
    {
      Tx x = ip.x;
      Tx y = ip.y;
      Tx z = ip.z;
      
      Tx lam3 = 1-x-y;
      
      shape[0] = x * (2*x-1) * (1-z);
      shape[1] = y * (2*y-1) * (1-z);
      shape[2] = lam3 * (2*lam3-1) * (1-z);
      shape[3] = x * (2*x-1) * z;
      shape[4] = y * (2*y-1) * z;
      shape[5] = lam3 * (2*lam3-1) * z;
      
      shape[6] = 4 * x * lam3 * (1-z);
      shape[7] = 4 * x * y       * (1-z);
      shape[8] = 4 * y * lam3 * (1-z);
      shape[9] = 4 * x * lam3 * z;
      shape[10] = 4 * x * y       * z;
      shape[11] = 4 * y * lam3 * z;
    }
  };



  class FE_Prism2HBaniso : public T_ScalarFiniteElementFO<FE_Prism2HBaniso,ET_PRISM,12,2>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
    {
      Tx x = ip.x;
      Tx y = ip.y;
      Tx z = ip.z;
      
      shape[0] = x * (1-z);
      shape[1] = y * (1-z);
      shape[2] = (1-x-y) * (1-z);
      shape[3] = x * z;
      shape[4] = y * z;
      shape[5] = (1-x-y) * z;
      
      shape[6] = 4 * x * (1-x-y) * (1-z);
      shape[7] = 4 * x * y       * (1-z);
      shape[8] = 4 * y * (1-x-y) * (1-z);
      shape[9] = 4 * x * (1-x-y) * z;
      shape[10] = 4 * x * y       * z;
      shape[11] = 4 * y * (1-x-y) * z;
    }
  };








  /* ***************************** Hex *********************************** */

  /// P0 hex element
  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_HEX,0> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
  {
    shape[0] = Tx(1);
  }
  using FE_Hex0 = ScalarFE<ET_HEX,0>;
  
  /*
  /// trilinear hex element
  class FE_Hex1 : public T_ScalarFiniteElementFO<FE_Hex1,ET_HEX,8,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
    {
      Tx x = ip.x;
      Tx y = ip.y;
      Tx z = ip.z;

      shape[0] = (1-x) * (1-y) * (1-z);
      shape[1] =    x  * (1-y) * (1-z);
      shape[2] =    x  *    y  * (1-z);
      shape[3] = (1-x) *    y  * (1-z);
      shape[4] = (1-x) * (1-y) *    z ;
      shape[5] =    x  * (1-y) *    z ;
      shape[6] =    x  *    y  *    z ;
      shape[7] = (1-x) *    y  *    z ;
    }
  };
  */
      
  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_HEX,1> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
  {
    Tx x = ip.x;
    Tx y = ip.y;
    Tx z = ip.z;
    
    shape[0] = (1-x) * (1-y) * (1-z);
    shape[1] =    x  * (1-y) * (1-z);
    shape[2] =    x  *    y  * (1-z);
    shape[3] = (1-x) *    y  * (1-z);
    shape[4] = (1-x) * (1-y) *    z ;
    shape[5] =    x  * (1-y) *    z ;
    shape[6] =    x  *    y  *    z ;
    shape[7] = (1-x) *    y  *    z ;
  };
  using FE_Hex1 = ScalarFE<ET_HEX,1>;


  class FE_Hex20 : public T_ScalarFiniteElementFO<FE_Hex20,ET_HEX,20,2>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
    {
      Tx x = ip.x;
      Tx y = ip.y;
      Tx z = ip.z;

      Tx lam[8] = {(1-x)*(1-y)*(1-z),x*(1-y)*(1-z),x*y*(1-z),(1-x)*y*(1-z),
                   (1-x)*(1-y)*z,x*(1-y)*z,x*y*z,(1-x)*y*z}; 
      Tx sigma[8]={(1-x)+(1-y)+(1-z),x+(1-y)+(1-z),x+y+(1-z),(1-x)+y+(1-z),
                   (1-x)+(1-y)+z,x+(1-y)+z,x+y+z,(1-x)+y+z}; 
      
      Tx shapee[12];

      static const int e[12][2] =
        {
          { 0, 1 }, { 2, 3 }, { 3, 0 }, { 1, 2 },
          { 4, 5 }, { 6, 7 }, { 7, 4 }, { 5, 6 },
          { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 },
        };
      
      for (int i = 0; i < 12; i++)
        {
          auto lame = lam[e[i][0]]+lam[e[i][1]];
          auto xi = sigma[e[i][1]]-sigma[e[i][0]];
          shapee[i] = (1-xi*xi)*lame;
        }
      for (int i = 0; i < 12; i++)
        {
          lam[e[i][0]] -= 0.5 * shapee[i];
          lam[e[i][1]] -= 0.5 * shapee[i];
        }
      for (int i = 0; i < 8; i++)
        shape[i] = lam[i];
      for (int i = 0; i < 12; i++)
        shape[i+8] = shapee[i];
    }
  }; 

  
  /* ***************************** Pyramid *********************************** */

  ///
  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_PYRAMID,0> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
  {
    shape[0] = Tx(1.0);
  }
  using FE_Pyramid0 = ScalarFE<ET_PYRAMID,0>;
  
  ///
  /*
  class FE_Pyramid1 : public T_ScalarFiniteElementFO<FE_Pyramid1,ET_PYRAMID,5,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
    {
      Tx x = ip.x;
      Tx y = ip.y;
      Tx z = ip.z;

      // if (z == 1) z -= 1e-10;
      z -= 1e-10;

      shape[0] = (1-z-x)*(1-z-y) / (1-z);
      shape[1] = x*(1-z-y) / (1-z);
      shape[2] = x*y / (1-z);
      shape[3] = (1-z-x)*y / (1-z);
      shape[4] = z;
    }
  };
  */

  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_PYRAMID,1> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
  {
    Tx x = ip.x;
    Tx y = ip.y;
    Tx z = ip.z;
    
    // if (z == 1) z -= 1e-10;
    z -= 1e-10;
    
    shape[0] = (1-z-x)*(1-z-y) / (1-z);
    shape[1] = x*(1-z-y) / (1-z);
    shape[2] = x*y / (1-z);
    shape[3] = (1-z-x)*y / (1-z);
    shape[4] = z;
  }
  using FE_Pyramid1 = ScalarFE<ET_PYRAMID,1>;




  
  /* ***************************** Hexamid *********************************** */

  ///
  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_HEXAMID,0> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
  {
    shape[0] = Tx(1.0);
  }

  template<> template<typename Tx, typename TFA>  
  void ScalarFE<ET_HEXAMID,1> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
  {
    Tx y = ip.y;
    Tx z = ip.z;
    Tx den = (1-y)*(1-z);
    den += Tx(1e-12);
    Tx x = ip.x / den;    
          
    shape[0] = (1-x)*(1-y)*(1-z);
    shape[1] = (  x)*(1-y)*(1-z);
    shape[2] = (  x)*(  y)*(1-z);
    shape[3] = (1-x)*(  y)*(1-z);
    shape[4] = (1-x)*(1-y)*(  z);
    shape[5] = (  x)*(1-y)*(  z);
    shape[6] =       (  y)*(  z);
    
          /*
    // if (z == 1) z -= 1e-10;
    z -= 1e-10;
    
    shape[0] = (1-z-x)*(1-z-y) / (1-z);
    shape[1] = x*(1-z-y) / (1-z);
    shape[2] = x*y / (1-z);
    shape[3] = (1-z-x)*y / (1-z);
    shape[4] = z;
          */
  }


  

#ifdef FILE_H1LOFE_CPP

  template<>
  FlatArray<IntegrationPoint> ScalarFE<ET_TRIG,1> :: GetNodalPoints() const
  {
    static IntegrationPoint ipts[] = {
      { 1., 0. }, { 0., 1. }, { 0., 0. }
    };

    return { 3, &ipts[0] };
  }

  template<>
  FlatArray<IntegrationPoint> ScalarFE<ET_TET,1> :: GetNodalPoints() const
  {
    static IntegrationPoint ipts[] = {
      { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 0, 0, 0 }
    };

    return { 4, &ipts[0] };
  }

  template<>
  FlatArray<IntegrationPoint> ScalarFE<ET_TET,2> :: GetNodalPoints() const
  {
    static IntegrationPoint ipts[] = {
      { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 0, 0, 0 },
      { 0.5, 0.5, 0 }, { 0.5, 0, 0.5 }, { 0.5, 0, 0 },
      { 0, 0.5, 0.5 }, { 0, 0.5, 0 }, { 0, 0, 0.5 } 
    };

    return { 10, &ipts[0] };
  }

  /*
  template<> 
  void ScalarFE<ET_TET,1> :: Interpolate (const ElementTransformation & trafo, 
                                          const class CoefficientFunction & func, SliceMatrix<> coefs,
                                          LocalHeap & lh) const
  {
    HeapReset hr(lh);
    IntegrationPoint ipts[] = {
      { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 0, 0, 0 }
    };
    IntegrationRule ir(4, ipts);
    auto & mir = trafo(ir, lh);
    func.Evaluate (mir, coefs);
  }
  */
#endif
  




  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Point,ET_POINT>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Segm0,ET_SEGM>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Segm1,ET_SEGM>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Segm2,ET_SEGM>;


  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Segm2HB,ET_SEGM>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Segm1L2,ET_SEGM>;

  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Segm2L2,ET_SEGM>;

  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_NcSegm1,ET_SEGM>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Segm3Pot,ET_SEGM>;


  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_TSegmL2<0>,ET_SEGM>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_TSegmL2<1>,ET_SEGM>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_TSegmL2<2>,ET_SEGM>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_TSegmL2<3>,ET_SEGM>;

  H1LOFE_EXTERN template class T_ScalarFiniteElement<ScalarFE<ET_TRIG,0>,ET_TRIG>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<ScalarFE<ET_TRIG,1>,ET_TRIG>;
  // H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Trig2,ET_TRIG>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Trig2HB,ET_TRIG>;

  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_NcTrig1,ET_TRIG>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<ScalarFE<ET_QUAD,0>,ET_QUAD>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<ScalarFE<ET_QUAD,1>,ET_QUAD>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Quad2,ET_QUAD>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Quad2aniso,ET_QUAD>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Quad2Serendipity,ET_QUAD>;

  H1LOFE_EXTERN template class T_ScalarFiniteElement<ScalarFE<ET_TET,0>,ET_TET>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<ScalarFE<ET_TET,1>,ET_TET>;
  // H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Tet1,ET_TET>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<ScalarFE<ET_TET,2>,ET_TET>;  
  // H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Tet2,ET_TET>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Tet2HB,ET_TET>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_NcTet1,ET_TET>;


  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Prism0,ET_PRISM>;
  // H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Prism1,ET_PRISM>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<ScalarFE<ET_PRISM,1>,ET_PRISM>;  
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Prism2,ET_PRISM>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Prism2aniso,ET_PRISM>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Prism2HBaniso,ET_PRISM>;

  H1LOFE_EXTERN template class T_ScalarFiniteElement<ScalarFE<ET_HEXAMID,0>,ET_HEXAMID>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<ScalarFE<ET_HEXAMID,1>,ET_HEXAMID>;

  
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Hex0,ET_HEX>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Hex1,ET_HEX>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Hex20,ET_HEX>;
  
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Pyramid0,ET_PYRAMID>;
  // H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Pyramid1,ET_PYRAMID>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<ScalarFE<ET_PYRAMID,1>,ET_PYRAMID>;    










  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Point,ET_POINT,1,0>;

  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Segm0,ET_SEGM,1,0>;
  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Segm1,ET_SEGM,2,1>;
  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Segm2,ET_SEGM,3,2>;

  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Segm2HB,ET_SEGM,3,2>;
  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Segm1L2,ET_SEGM,2,1>;

  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Segm2L2,ET_SEGM,3,2>;

  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_NcSegm1,ET_SEGM,1,1>;
  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Segm3Pot,ET_SEGM,4,3>;


  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_TSegmL2<0>,ET_SEGM,1,0>;
  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_TSegmL2<1>,ET_SEGM,2,1>;
  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_TSegmL2<2>,ET_SEGM,3,2>;
  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_TSegmL2<3>,ET_SEGM,4,3>;

  /*
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Trig0,ET_TRIG>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Trig1,ET_TRIG>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Trig2,ET_TRIG>;
  H1LOFE_EXTERN template class T_ScalarFiniteElement<FE_Trig2HB,ET_TRIG>;
  */
  // H1LOFE_EXTERN template class T_ScalarFiniteElementFO<ScalarFE<ET_TRIG,0>,ET_TRIG,1,0>;
  // H1LOFE_EXTERN template class T_ScalarFiniteElementFO<ScalarFE<ET_TRIG,1>,ET_TRIG,3,1>;
  // H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_NcTrig1,ET_TRIG,3,1>;
  // H1LOFE_EXTERN template class T_ScalarFiniteElementFO<ScalarFE<ET_QUAD,0>,ET_QUAD,1,0>;
  // H1LOFE_EXTERN template class T_ScalarFiniteElementFO<ScalarFE<ET_QUAD,1>,ET_QUAD,4,1>;
  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Quad2,ET_QUAD,9,2>;
  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Quad2aniso,ET_QUAD,6,2>;
  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Quad2Serendipity,ET_QUAD,8,2>;

  // H1LOFE_EXTERN template class T_ScalarFiniteElementFO<ScalarFE<ET_TET,0>,ET_TET,1,0>;
  // H1LOFE_EXTERN template class T_ScalarFiniteElementFO<ScalarFE<ET_TET,1>,ET_TET,4,1>;
  // H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Tet1,ET_TET,4,1>;
  // H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Tet2,ET_TET,10,2>;
  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Tet2HB,ET_TET,10,2>;
  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_NcTet1,ET_TET,4,1>;


  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Prism0,ET_PRISM,1,0>;
  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Prism1,ET_PRISM,6,1>;
  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Prism2,ET_PRISM,18,2>;
  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Prism2aniso,ET_PRISM,12,2>;
  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Prism2HBaniso,ET_PRISM,12,2>;


  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Hex0,ET_HEX,1,0>;
  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Hex1,ET_HEX,8,1>;
  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Hex20,ET_HEX,20,2>;
  
  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Pyramid0,ET_PYRAMID,1,0>;
  H1LOFE_EXTERN template class T_ScalarFiniteElementFO<FE_Pyramid1,ET_PYRAMID,5,1>;

  H1LOFE_EXTERN template class ScalarFE<ET_POINT,0>;
  H1LOFE_EXTERN template class ScalarFE<ET_POINT,1>;

  H1LOFE_EXTERN template class ScalarFE<ET_SEGM,0>;
  H1LOFE_EXTERN template class ScalarFE<ET_SEGM,1>;

  H1LOFE_EXTERN template class ScalarFE<ET_TRIG,0>;
  H1LOFE_EXTERN template class ScalarFE<ET_TRIG,1>;
  H1LOFE_EXTERN template class ScalarFE<ET_TRIG,2>;

  H1LOFE_EXTERN template class ScalarFE<ET_QUAD,0>;
  H1LOFE_EXTERN template class ScalarFE<ET_QUAD,1>;

  H1LOFE_EXTERN template class ScalarFE<ET_TET,0>;
  H1LOFE_EXTERN template class ScalarFE<ET_TET,1>;
  H1LOFE_EXTERN template class ScalarFE<ET_TET,2>;

  H1LOFE_EXTERN template class ScalarFE<ET_PRISM,0>;
  H1LOFE_EXTERN template class ScalarFE<ET_PRISM,1>;

  H1LOFE_EXTERN template class ScalarFE<ET_PYRAMID,0>;
  H1LOFE_EXTERN template class ScalarFE<ET_PYRAMID,1>;

  H1LOFE_EXTERN template class ScalarFE<ET_HEXAMID,0>;
  H1LOFE_EXTERN template class ScalarFE<ET_HEXAMID,1>;
  
  H1LOFE_EXTERN template class ScalarFE<ET_HEX,0>;
  H1LOFE_EXTERN template class ScalarFE<ET_HEX,1>;
}




#endif
