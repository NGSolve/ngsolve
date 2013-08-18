/*********************************************************************/
/* File:   h1hofe.cpp                                                */
/* Author: Start                                                      */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/

 
#include <fem.hpp>
#include <h1hofefo.hpp>

#include "tscalarfe_impl.hpp"


namespace ngfem
{

  /*
  template <ELEMENT_TYPE ET, int ORDER>
  T_H1HighOrderFiniteElementFO<ET,ORDER> ::
  ~T_H1HighOrderFiniteElementFO ()
  {
    ;
  }

  template <int ORDER>
  H1HighOrderFEFO<ET_TRIG, ORDER> :: ~H1HighOrderFEFO()
  { ; }

  template <int ORDER>
  H1HighOrderFEFO<ET_TET, ORDER> :: ~H1HighOrderFEFO()
  { ; }
  */




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
	
	LegendrePolynomial::EvalScaledMult (ORDER-2, 
					    lami[e[1]]-lami[e[0]], lami[e[0]]+lami[e[1]], 
					    lami[e[0]]*lami[e[1]], shape.Addr(ii));
	ii += ORDER-1;
	/*
        int es = edges[i][0], ee = edges[i][1];
        if (vnums[es] > vnums[ee]) swap (es, ee);
        
        ii += T_ORTHOPOL::CalcScaled<ORDER> 
          (lami[ee]-lami[es], lami[es]+lami[ee], shape.Addr(ii));
	  */
      }

    // inner dofs
    int fav[3] = { 0, 1, 2 }; 
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
    if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
    if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	
    
    Tx polx[ORDER-2], poly[ORDER-2];
    
    T_TRIGSHAPES::CalcSplitted<ORDER> (lami[fav[2]]-lami[fav[1]],
                                       lami[fav[0]], polx, poly);
    
    TrigProduct<ORDER-3>::Do (polx, poly, shape.Addr(ii));
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
					    lami[e[0]]*lami[e[1]], shape.Addr(ii));
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
    Tx lami[4] = { hx[0], hx[1], hx[2], 1-hx[0]-hx[1]-hx[2] };

    for (int i = 0; i < 4; i++)
      shape[i] = lami[i];

    int ii = 4;

    // edge dofs
    // const EDGE * edges = ElementTopology::GetEdges (ET_TET);
    for (int i = 0; i < 6; i++)
      { 
        // int es = edges[i][0], ee = edges[i][1];
        // if (vnums[es] > vnums[ee]) swap (es, ee);
        INT<2> e = GetEdgeSort (i, vnums);
        ii += T_ORTHOPOL::CalcScaled<ORDER> 
          (lami[e[1]]-lami[e[0]], lami[e[1]]+lami[e[0]], shape.Addr(ii));
      }

    // face dofs
    for (int i = 0; i < 4; i++)
      if (ORDER >= 3)
	{
          INT<4> f = GetFaceSort (i, vnums);
	  int vop = 6 - f[0] - f[1] - f[2];  	
          
	  ii += T_FACESHAPES::Calc (ORDER, 
				    lami[f[2]]-lami[f[1]],
				    lami[f[0]], lami[vop],  shape.Addr(ii));
	}

    if (ORDER >= 4)
      ii += T_INNERSHAPES::Calc (ORDER,
                                 lami[0]-lami[3], lami[1], lami[2], 
                                 shape.Addr(ii) );
  }











  template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TRIG,1>, ET_TRIG>;
  template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TRIG,2>, ET_TRIG>;
  template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TRIG,3>, ET_TRIG>;
  template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TRIG,4>, ET_TRIG>;
  template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TRIG,5>, ET_TRIG>;
  template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TRIG,6>, ET_TRIG>;

  template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TET,1>, ET_TET>;
  template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TET,2>, ET_TET>;
  template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TET,3>, ET_TET>;
  template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TET,4>, ET_TET>;
  template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TET,5>, ET_TET>;
  template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TET,6>, ET_TET>;



  /*
  template class T_H1HighOrderFiniteElementFO<ET_TRIG>;
  template class T_H1HighOrderFiniteElementFO<ET_TET>;
  */
  template class H1HighOrderFEFO<ET_TRIG,1>;
  template class H1HighOrderFEFO<ET_TRIG,2>;
  template class H1HighOrderFEFO<ET_TRIG,3>;
  template class H1HighOrderFEFO<ET_TRIG,4>;
  template class H1HighOrderFEFO<ET_TRIG,5>;
  template class H1HighOrderFEFO<ET_TRIG,6>;

  template class H1HighOrderFEFO<ET_TET,1>;
  template class H1HighOrderFEFO<ET_TET,2>;
  template class H1HighOrderFEFO<ET_TET,3>;
  template class H1HighOrderFEFO<ET_TET,4>;
  template class H1HighOrderFEFO<ET_TET,5>;
  template class H1HighOrderFEFO<ET_TET,6>;


  int link_it_h1hofefo;
}
 
