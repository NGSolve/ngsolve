/*********************************************************************/
/* File:   h1hofefo.cpp                                              */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/


#define FILE_H1HOFEFO_CPP
 
#include <fem.hpp>
#include <h1hofefo.hpp>

#include <tscalarfe_impl.hpp>


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
  void H1HighOrderFEFO<ET_TRIG, ORDER> :: T_CalcShape (Tx hx[2], TFA & shape) const
  {
    Tx lam[3] = { hx[0], hx[1], 1-hx[0]-hx[1] };

    for (int i = 0; i < 3; i++)
      shape[i] = lam[i];

    int ii = 3;

    // edge dofs
    const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
    for (int i = 0; i < 3; i++)
      { 
        int es = edges[i][0], ee = edges[i][1];
        if (vnums[es] > vnums[ee]) swap (es, ee);
	INT<2> e (es, ee);

	// INT<2> e = GetEdgeSort (i, vnums);
	
	LegendrePolynomial::EvalScaledMult (ORDER-2, 
					    lam[e[1]]-lam[e[0]], lam[e[0]]+lam[e[1]], 
					    lam[e[0]]*lam[e[1]], shape.Addr(ii));
	ii += ORDER-1;
	/*
        int es = edges[i][0], ee = edges[i][1];
        if (vnums[es] > vnums[ee]) swap (es, ee);
        
        ii += T_ORTHOPOL::CalcScaled<ORDER> 
          (lam[ee]-lam[es], lam[es]+lam[ee], shape.Addr(ii));
	  */
      }

    // inner dofs

    if (ORDER >= 3)
      {
        INT<4> f = GetFaceSort (0, vnums);

	DubinerBasis3::EvalMult (ORDER-3, 
				 lam[f[0]], lam[f[1]], 
				 lam[f[0]]*lam[f[1]]*lam[f[2]], shape+ii);
      }
    /*
    int fav[3] = { 0, 1, 2 }; 
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
    if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	
    
    Tx polx[ORDER-2], poly[ORDER-2];
    
    T_TRIGSHAPES::CalcSplitted<ORDER> (lam[fav[2]]-lam[fav[1]],
                                       lam[fav[0]], polx, poly);
    
    TrigProduct<ORDER-3>::Do (polx, poly, shape.Addr(ii));
    */
  }

  template <>   template<typename Tx, typename TFA>  
  void H1HighOrderFEFO<ET_TRIG, 2> :: T_CalcShape (Tx hx[2], TFA & shape) const
  {
    Tx lami[3] = { hx[0], hx[1], 1-hx[0]-hx[1] };
    
    for (int i = 0; i < 3; i++)
      shape[i] = lami[i];
    
    int ii = 3;

    // edge dofs
    const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
    for (int i = 0; i < 3; i++)
      { 
        int es = edges[i][0], ee = edges[i][1];
        if (vnums[es] > vnums[ee]) swap (es, ee);
	INT<2> e (es, ee);

	// INT<2> e = GetEdgeSort (i, vnums);
	
	LegendrePolynomial::EvalScaledMult (0, 
					    lami[e[1]]-lami[e[0]], lami[e[0]]+lami[e[1]], 
					    lami[e[0]]*lami[e[1]], shape+ii);
	ii += 1;

	/*
        int es = edges[i][0], ee = edges[i][1];
        if (vnums[es] > vnums[ee]) swap (es, ee);
        
        ii += T_ORTHOPOL::CalcScaled<2> 
          (lami[ee]-lami[es], lami[es]+lami[ee], shape.Addr(ii));
	*/
      }
  }

  template <> template<typename Tx, typename TFA>  
  void H1HighOrderFEFO<ET_TRIG, 1> :: T_CalcShape (Tx hx[2], TFA & shape) const
  {
    Tx lami[3] = { hx[0], hx[1], 1-hx[0]-hx[1] };

    for (int i = 0; i < 3; i++)
      shape[i] = lami[i];
  }











  static Timer ttet("tetfo");
  
  template <int ORDER>   template<typename Tx, typename TFA>  
  void H1HighOrderFEFO<ET_TET, ORDER> :: T_CalcShape (Tx hx[3], TFA & shape) const
  {
    Tx lam[4] = { hx[0], hx[1], hx[2], 1-hx[0]-hx[1]-hx[2] };

    for (int i = 0; i < 4; i++)
      shape[i] = lam[i];

    int ii = 4;

    // edge dofs
    // const EDGE * edges = ElementTopology::GetEdges (ET_TET);
    for (int i = 0; i < 6; i++)
      { 
        INT<2> e = GetEdgeSort (i, vnums);
        LegendrePolynomial::
          EvalScaledMultFO<ORDER-2> (lam[e[1]]-lam[e[0]], lam[e[0]]+lam[e[1]], 
                                     lam[e[0]]*lam[e[1]], shape+ii);
        ii += ORDER-1;

        // ii += T_ORTHOPOL::CalcScaled<ORDER> 
        // (lam[e[1]]-lam[e[0]], lam[e[1]]+lam[e[0]], shape.Addr(ii));
      }

    // face dofs
    for (int i = 0; i < 4; i++)
      if (ORDER >= 3)
	{
          INT<4> f = GetFaceSort (i, vnums);
	  int vop = 6 - f[0] - f[1] - f[2];  	

	  DubinerBasis3::EvalScaledMult (ORDER-3, lam[f[0]], lam[f[1]], 1-lam[vop], 
					 lam[f[0]]*lam[f[1]]*lam[f[2]], shape+ii);
	  ii += (ORDER-2)*(ORDER-1)/2;
          /*
	  ii += T_FACESHAPES::Calc (ORDER, 
				    lam[f[2]]-lam[f[1]],
				    lam[f[0]], lam[vop],  shape.Addr(ii));
          */
	}

    if (ORDER >= 4)
      ii += T_INNERSHAPES::Calc (ORDER,
                                 lam[0]-lam[3], lam[1], lam[2], 
                                 shape.Addr(ii) );
  }











  int link_it_h1hofefo;
}
 
