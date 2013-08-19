#ifndef FILE_H1LOFE
#define FILE_H1LOFE

/*********************************************************************/
/* File:   h1lofe.hpp                                                */
/* Author: Start                                                     */
/* Date:   29. Jun. 2009                                              */
/*********************************************************************/



namespace ngfem
{

  template <ELEMENT_TYPE ET>
  class ScalarDummyFE : public T_ScalarFiniteElementFO<ScalarDummyFE<ET>,ET,0,0>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx x[1], TFA & shape) 
    { ; }
  };


  class FE_Point : public T_ScalarFiniteElementFO<FE_Point,ET_POINT,1,0>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx x[], TFA & shape) 
    {
      shape[0] = 1;
    }
  }; 


  ///
  class FE_Segm0 : public T_ScalarFiniteElementFO<FE_Segm0,ET_SEGM,1,0>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx x[1], TFA & shape) 
    {
      shape[0] = 1;
    }
  }; 

  ///
  class FE_Segm1 : public T_ScalarFiniteElementFO<FE_Segm1,ET_SEGM,2,1>
  {
  public:
    NGS_DLL_HEADER FE_Segm1() { ; };
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx x[1], TFA & shape) 
    {
      shape[0] = x[0];
      shape[1] = 1-x[0];
    }
  }; 

  ///
  class FE_Segm2 : public T_ScalarFiniteElementFO<FE_Segm2,ET_SEGM,3,2>
  {
  public:
        NGS_DLL_HEADER  FE_Segm2() { ; };
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
  class FE_Segm2HB : public T_ScalarFiniteElementFO<FE_Segm2HB,ET_SEGM,3,2>
  {
  public:
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
  class FE_Segm1L2 : public T_ScalarFiniteElementFO<FE_Segm1L2,ET_SEGM,2,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx x[1], TFA & shape) 
    {
      shape[0] = 1;
      shape[1] = 2*x[0]-1;
    }
  }; 

  ///
  class FE_Segm2L2 :public T_ScalarFiniteElementFO<FE_Segm2L2,ET_SEGM,3,2>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[1], TFA & shape) 
    {
      Tx x = hx[0];

      shape[0] = 1;
      shape[1] = 2*x-1;
      shape[2] = (2*x-1)*(2*x-1)-1.0/3.0;
    }
  }; 

  ///
  class FE_NcSegm1 :public T_ScalarFiniteElementFO<FE_NcSegm1,ET_SEGM,1,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[1], TFA & shape) 
    {
      shape[0] = 1;
    }
  }; 

  /// potential space for Nedelec IIb
  class FE_Segm3Pot :public T_ScalarFiniteElementFO<FE_Segm3Pot,ET_SEGM,4,3>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[1], TFA & shape) 
    {
      Tx x = hx[0];
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
    // static IPDataArray ipdata;
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[1], TFA & shape) 
    {
      Tx x = hx[0];
      // Tx lam2 = 1-x;

      // very primitive ...
      // shape = 0;
      
      if (ORDER >= 0) shape[0] = 1;
      if (ORDER >= 1) shape[1] = 2*x-1;
      if (ORDER >= 2) shape[2] = (2*x-1)*(2*x-1)-1.0/3.0;
      if (ORDER >= 3) shape[3] = (2*x-1)*(2*x-1)*(2*x-1);
      if (ORDER >= 4)
        {
          throw Exception ("TSegmL2: Legendre polynomials not implemented");
        }
    }
  };






  /* ***************************** Trig Elements *************************************** */

  ///
  class FE_Trig0 : public T_ScalarFiniteElement<FE_Trig0,ET_TRIG>
  {
  public:
    FE_Trig0() { ndof = 1; order = 0; }
    template<typename Tx, typename TFA>  
    inline void T_CalcShape (Tx x[2], TFA & shape) const
    {
      shape[0] = 1;
    }
  }; 

  ///
  class FE_Trig1 : public T_ScalarFiniteElement<FE_Trig1,ET_TRIG>
  {
  public:
    FE_Trig1() { ndof = 3; order = 1;}
    template<typename Tx, typename TFA>  
    inline void T_CalcShape (Tx x[2], TFA & shape) const
    {
      shape[0] = x[0];
      shape[1] = x[1];      
      shape[2] = 1-x[0]-x[1];
    }
  }; 


  ///
  class FE_Trig2 : public T_ScalarFiniteElement<FE_Trig2,ET_TRIG>
  {
  public:
    FE_Trig2() { ndof = 6; order = 2;}
    template<typename Tx, typename TFA>  
    inline void T_CalcShape (Tx hx[2], TFA & shape) const
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

  class FE_Trig2HB : public T_ScalarFiniteElement<FE_Trig2HB,ET_TRIG>
  {
  public:
    FE_Trig2HB() { ndof = 6; order = 2;}
    template<typename Tx, typename TFA>  
    inline void T_CalcShape (Tx hx[2], TFA & shape) const
    {
      Tx x = hx[0];
      Tx y = hx[1];
      Tx lam3 = 1-x-y;
    
      shape[0] = x;
      shape[1] = y;
      shape[2] = lam3;
      shape[3] = 4 * y * lam3;
      shape[4] = 4 * x * lam3;
      shape[5] = 4 * x * y;
    }
  }; 



  class FE_NcTrig1 : public T_ScalarFiniteElementFO<FE_NcTrig1,ET_TRIG,3,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[2], TFA & shape) 
    {
      Tx x = hx[0];
      Tx y = hx[1];
      Tx lam3 = 1-x-y;
    
      shape[0] = 1-2*y;
      shape[1] = 1-2*x;
      shape[2] = 1-2*lam3;
    }
  };




  /* ***************************** Quad *************************************** */


  class FE_Quad0 : public T_ScalarFiniteElementFO<FE_Quad0,ET_QUAD,1,0>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[2], TFA & shape) 
    {
      shape[0] = 1.0;
    }
  };


  /// quad of order 1
  class FE_Quad1 : public T_ScalarFiniteElementFO<FE_Quad1,ET_QUAD,4,1>
  {
  public:
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
  }; 


  /// quad or order 2
  class FE_Quad2 : public T_ScalarFiniteElementFO<FE_Quad2,ET_QUAD,9,2>
  {
  public:
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
  }; 


  class FE_Quad2aniso :  public T_ScalarFiniteElementFO<FE_Quad2aniso,ET_QUAD,6,2>
  {
  public:
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
  class FE_Tet0 : public T_ScalarFiniteElementFO<FE_Tet0,ET_TET,1,0>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[3], TFA & shape) 
    {
      shape[0] = 1;
    }
  };


  ///
  class FE_Tet1 : public T_ScalarFiniteElementFO<FE_Tet1,ET_TET,4,1>
  {
  public:
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
  };


  ///
  class FE_Tet2 : public T_ScalarFiniteElementFO<FE_Tet2,ET_TET,10,2>
  {
  public:
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
  class FE_Tet2HB : public T_ScalarFiniteElementFO<FE_Tet2HB,ET_TET,10,2>
  {
  public:
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



  class FE_NcTet1 : public T_ScalarFiniteElementFO<FE_NcTet1,ET_TET,4,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[2], TFA & shape) 
    {
      Tx x = hx[0];
      Tx y = hx[1];
      Tx z = hx[2];
      Tx lam4 = 1-x-y-z;
    
      shape[0] = 1-2*x;
      shape[1] = 1-2*y;
      shape[2] = 1-2*z;
      shape[3] = 1-2*lam4;
    }
  };


  /* ***************************** Prism *********************************** */

  ///
  class FE_Prism0 : public T_ScalarFiniteElementFO<FE_Prism0,ET_PRISM,1,0>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[3], TFA & shape) 
    {
      shape[0] = 1;
    }
  };

  ///
  class FE_Prism1 : public T_ScalarFiniteElementFO<FE_Prism1,ET_PRISM,6,1>
  {
  public:
    NGS_DLL_HEADER FE_Prism1() { ; }

    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[3], TFA & shape) 
    {
      Tx x = hx[0];
      Tx y = hx[1];
      Tx z = hx[2];
      
      shape[0] = x * (1-z);
      shape[1] = y * (1-z);
      shape[2] = (1-x-y) * (1-z);
      shape[3] = x * z;
      shape[4] = y * z;
      shape[5] = (1-x-y) * z;
    }
  };


  ///
  class FE_Prism2 : public T_ScalarFiniteElementFO<FE_Prism2,ET_PRISM,18,2>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[3], TFA & shape) 
    {
      Tx x = hx[0];
      Tx y = hx[1];
      Tx z = hx[2];
      

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
    static void T_CalcShape (Tx hx[3], TFA & shape) 
    {
      Tx x = hx[0];
      Tx y = hx[1];
      Tx z = hx[2];
      
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
    static void T_CalcShape (Tx hx[3], TFA & shape) 
    {
      Tx x = hx[0];
      Tx y = hx[1];
      Tx z = hx[2];
      
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
  class FE_Hex0 : public T_ScalarFiniteElementFO<FE_Hex0,ET_HEX,1,0>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[3], TFA & shape) 
    {
      shape[0] = 1;
    }
  };

  /// trilinear hex element
  class FE_Hex1 : public T_ScalarFiniteElementFO<FE_Hex1,ET_HEX,8,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[3], TFA & shape) 
    {
      Tx x = hx[0];
      Tx y = hx[1];
      Tx z = hx[2];

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




  /* ***************************** Pyramid *********************************** */

  ///
  class FE_Pyramid0 : public T_ScalarFiniteElementFO<FE_Pyramid0,ET_PYRAMID,1,0>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[3], TFA & shape) 
    {
      shape[0] = 1;
    }
  };


  ///
  class FE_Pyramid1 : public T_ScalarFiniteElementFO<FE_Pyramid1,ET_PYRAMID,5,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx hx[3], TFA & shape) 
    {
      Tx x = hx[0];
      Tx y = hx[1];
      Tx z = hx[2];

      // if (z == 1) z -= 1e-10;
      z -= 1e-10;

      shape[0] = (1-z-x)*(1-z-y) / (1-z);
      shape[1] = x*(1-z-y) / (1-z);
      shape[2] = x*y / (1-z);
      shape[3] = (1-z-x)*y / (1-z);
      shape[4] = z;
    }
  };




}




#endif
