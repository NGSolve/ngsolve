#ifndef FILE_H1HOFEFO_IMPL
#define FILE_H1HOFEFO_IMPL

/*********************************************************************/
/* File:   h1hofefo_impl.hpp                                         */
/* Author: Joachim Schoeberl                                         */
/* Date:   2009                                                      */
/*********************************************************************/

#include "recursive_pol_tet.hpp"



namespace ngfem
{

  template <int ORDER, int I = ORDER>
  class TrigProduct
  {
  public:
    template <class PX, class PY, class TRes>
    static void Do (const PX & polx, const PY & poly, TRes && res)
    {
      TrigProduct<ORDER, I-1>::Do (polx,poly, res);

      int ii = (ORDER+1)*(ORDER+2)/2 - (ORDER-I+1)*(ORDER-I+2)/2;

      for (int j = 0; j <= ORDER-I; j++)
        res[ii++] = polx[I] * poly[j];
    }
  };

  template <int ORDER>
  class TrigProduct<ORDER,-1>
  {
  public:
    template <class PX, class PY, class TRes>
    static void Do (const PX & polx, const PY & poly, TRes && res) { ; }
  };

  


  template <int ORDER>   template<typename Tx, typename TFA>  
  void H1HighOrderFEFO<ET_TRIG, ORDER> :: T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
  {
    Tx lam[3] = { ip.x, ip.y, 1-ip.x-ip.y };

    for (int i = 0; i < 3; i++) shape[i] = lam[i];

    int ii = 3;

    for (int i = 0; i < 3; i++)
      { 
        IVec<2> e = GetEdgeSort (i, vnums);
	LegendrePolynomial::EvalScaledMult (ORDER-2, 
					    lam[e[1]]-lam[e[0]], lam[e[0]]+lam[e[1]], 
					    lam[e[0]]*lam[e[1]], shape.Addr(ii));
	ii += ORDER-1;
      }

    // inner dofs
    if (ORDER >= 3)
      {
        IVec<4> f = GetFaceSort (0, vnums);
	DubinerBasis::EvalMult (ORDER-3, 
				 lam[f[0]], lam[f[1]], 
				 lam[f[0]]*lam[f[1]]*lam[f[2]], shape+ii);
      }
  }

  template <>   template<typename Tx, typename TFA>  
  void H1HighOrderFEFO<ET_TRIG, 2> :: T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
  {
    Tx lam[3] = { ip.x, ip.y, 1-ip.x-ip.y };
    for (int i = 0; i < 3; i++) shape[i] = lam[i];
    
    int ii = 3;
    for (int i = 0; i < 3; i++)
      { 
        IVec<2> e = GetEdge (i);
        shape[ii] = lam[e[0]] * lam[e[1]];
        ii++;
      }
  }

  template <> template<typename Tx, typename TFA>  
  void H1HighOrderFEFO<ET_TRIG, 1> :: T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
  {
    Tx lam[3] = { ip.x, ip.y, 1-ip.x-ip.y };
    for (int i = 0; i < 3; i++) shape[i] = lam[i];
  }







  
  template <int ORDER>   template<typename Tx, typename TFA>  
  void H1HighOrderFEFO<ET_TET, ORDER> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    Tx lam[4] = { ip.x, ip.y, ip.z, 1-ip.x-ip.y-ip.z };

    for (int i = 0; i < 4; i++) shape[i] = lam[i];

    int ii = 4;
    for (int i = 0; i < 6; i++)
      { 
        IVec<2> e = GetEdgeSort (i, vnums);
        LegendrePolynomial::
          EvalScaledMultFO<ORDER-2> (lam[e[1]]-lam[e[0]], lam[e[0]]+lam[e[1]], 
                                     lam[e[0]]*lam[e[1]], shape+ii);
        ii += ORDER-1;
      }

    // face dofs
    for (int i = 0; i < 4; i++)
      if (ORDER >= 3)
	{
          IVec<4> f = GetFaceSort (i, vnums);
	  int vop = 6 - f[0] - f[1] - f[2];  	

	  DubinerBasis::EvalScaledMult (ORDER-3, lam[f[0]], lam[f[1]], 1-lam[vop], 
                                        lam[f[0]]*lam[f[1]]*lam[f[2]], shape+ii);
	  ii += (ORDER-2)*(ORDER-1)/2;
	}

    if (ORDER >= 4)
      ii +=  TetShapesInnerLegendre::
        Calc (ORDER, lam[0]-lam[3], lam[1], lam[2], 
              shape.Addr(ii) );
  }



  template <> template<typename Tx, typename TFA>  
  void H1HighOrderFEFO<ET_TET, 2> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    // Tx lam[4] = { hx[0], hx[1], hx[2], 1-hx[0]-hx[1]-hx[2] };
    Tx lam[4] = { ip.x, ip.y, ip.z, 1-ip.x-ip.y-ip.z };    

    for (int i = 0; i < 4; i++) shape[i] = lam[i];

    int ii = 4;

    for (int i = 0; i < 6; i++)
      { 
        IVec<2> e = GetEdge (i);
        shape[ii] = lam[e[0]] * lam[e[1]];
        ii++;
      }
  }


  template <>   template<typename Tx, typename TFA>  
  void H1HighOrderFEFO<ET_TET, 3> :: T_CalcShape (TIP<3,Tx> ip, TFA & shape) const
  {
    Tx lam[4] = { ip.x, ip.y, ip.z, 1-ip.x-ip.y-ip.z };    
    // Tx lam[4] = { hx[0], hx[1], hx[2], 1-hx[0]-hx[1]-hx[2] };    

    for (int i = 0; i < 4; i++) shape[i] = lam[i];

    int ii = 4;
    for (int i = 0; i < 6; i++)
      { 
        IVec<2> e = GetEdgeSort (i, vnums);
        Tx bub = lam[e[0]]*lam[e[1]];
        shape[ii] = bub;
        shape[ii+1] = bub * (lam[e[1]]-lam[e[0]]);
        ii += 2;
      }

    for (int i = 0; i < 4; i++)
      {
        IVec<4> f = GetFace (i);
        shape[ii++] = lam[f[0]]*lam[f[1]]*lam[f[2]];
      }
  }

}
 
#endif

