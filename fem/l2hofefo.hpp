#ifndef FILE_L2HOFEFO
#define FILE_L2HOFEFO

/*********************************************************************/
/* File:   l2hofefo.hpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   Apr. 2009                                                 */
/*********************************************************************/

#include "tscalarfe.hpp"


namespace ngfem
{

  /**
     High order finite elements for L2 of fixed order
  */

  template<int DIM>
  class L2HighOrderFiniteElementFO : virtual public ScalarFiniteElement<DIM>
  {
  protected:
    int vnums[8];

  public:
    void SetVertexNumber (int nr, int vnum) { vnums[nr] = vnum; }
  };


  template <ELEMENT_TYPE ET, int ORDER> class L2HighOrderFEFO;

  template <ELEMENT_TYPE ET, int ORDER>
  class T_L2HighOrderFiniteElementFO : 
    public L2HighOrderFiniteElementFO<ET_trait<ET>::DIM>,
    public T_ScalarFiniteElement2< L2HighOrderFEFO<ET,ORDER>, ET >
  {
  protected:
    enum { DIM = ET_trait<ET>::DIM };

    using ScalarFiniteElement<DIM>::ndof;
    using ScalarFiniteElement<DIM>::order;
    using ScalarFiniteElement<DIM>::eltype;
    using ScalarFiniteElement<DIM>::dimspace;

    using L2HighOrderFiniteElementFO<DIM>::vnums;

  public:

    T_L2HighOrderFiniteElementFO () 
    {
      for (int i = 0; i < ET_trait<ET>::N_VERTEX; i++)
	vnums[i] = i;
      dimspace = DIM;
      eltype = ET;
      order = ORDER;
    }
  };




  /**
     High order triangular finite element
  */
  template <int ORDER>
  class L2HighOrderFEFO<ET_TRIG, ORDER> : public T_L2HighOrderFiniteElementFO<ET_TRIG, ORDER>
  {
    using T_L2HighOrderFiniteElementFO<ET_TRIG, ORDER>::ndof;
    using T_L2HighOrderFiniteElementFO<ET_TRIG, ORDER>::vnums; 

    typedef IntegratedLegendreMonomialExt T_ORTHOPOL;
    typedef TrigShapesInnerLegendre T_TRIGSHAPES;

  public:
    enum { NDOF = (ORDER+1)*(ORDER+2)/2 };
    L2HighOrderFEFO () { ndof = NDOF; }


    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[2], TFA & shape) const
    {
      Tx lami[3] = { hx[0], hx[1], 1-hx[0]-hx[1] };
      
      int fav[3] = { 0, 1, 2};
      if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
      if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
      if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	
      
      Tx x = lami[fav[0]]; 
      Tx y = lami[fav[1]];
      Tx l3 = lami[fav[2]];
      
      Vec<ORDER+1, Tx> polx, poly;
      
      LegendrePolynomialFO<ORDER>::EvalScaled (x-l3, 1-y, polx);
      LegendrePolynomialFO<ORDER>::Eval (2*y-1, poly);
      
      for (int i = 0, ii = 0; i <= ORDER; i++)
	{
	  // JacobiPolynomial (n-i, 2*y-1, 2*i+1, 0, poly);
	  for (int j = 0; j <= ORDER-i; j++)
	    shape[ii++] = polx[i] * poly[j];
	}
    }

  };

}


#endif
