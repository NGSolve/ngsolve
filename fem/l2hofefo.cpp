/*********************************************************************/
/* File:   l2hofe.cpp                                                */
/* Author: Start                                                      */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/

 
#include <fem.hpp>
#include <l2hofefo.hpp>


namespace ngfem
{

  


  /*

  template <ELEMENT_TYPE ET, int ORDER>
  void T_L2HighOrderFiniteElementFO<ET, ORDER> :: 
  CalcShape (const IntegrationPoint & ip, 
             FlatVector<> shape) const
  {
    double pt[DIM];
    for (int i = 0; i < DIM; i++) pt[i] = ip(i);
    static_cast<const L2HighOrderFEFO<ET,ORDER>*> (this) -> T_CalcShape (pt, shape); 
  }

  template <ELEMENT_TYPE ET, int ORDER>
  void T_L2HighOrderFiniteElementFO<ET, ORDER> :: 
  CalcDShape (const IntegrationPoint & ip, 
              FlatMatrixFixWidth<DIM> dshape) const
  {
    AutoDiff<DIM> adp[DIM];
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);

    DShapeAssign<DIM> ds(dshape); 
    static_cast<const L2HighOrderFEFO<ET,ORDER>*> (this) -> T_CalcShape (adp, ds);
  }


  template <ELEMENT_TYPE ET, int ORDER>
  void T_L2HighOrderFiniteElementFO<ET, ORDER> :: 
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
    static_cast<const L2HighOrderFEFO<ET, ORDER>*> (this) -> T_CalcShape (adp, ds);
  }


  template <ELEMENT_TYPE ET, int ORDER>
  double T_L2HighOrderFiniteElementFO<ET, ORDER> :: 
  Evaluate (const IntegrationPoint & ip, 
            FlatVector<double> x) const
  {
    double sum = 0.0;
    double pt[DIM];
    for (int i = 0; i < DIM; i++) pt[i] = ip(i);

    EvaluateShape ds(x, &sum);
    static_cast<const L2HighOrderFEFO<ET,ORDER>*> (this) -> T_CalcShape (pt, ds);
    return sum;
  }
  */


  /*



  template <int ORDER, int I = ORDER>
  class TrigProduct
  {
  public:
    template <class PX, class PY, class TRes>
    static void Do (const PX & polx, const PY & poly, TRes & res)
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
    static void Do (const PX & polx, const PY & poly, TRes & res) { ; }
  };

  */  

  /*
  template <int ORDER>   template<typename Tx, typename TFA>  
  void L2HighOrderFEFO<ET_TRIG, ORDER> :: T_CalcShape (Tx hx[2], TFA & shape) const
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
  */
  



  /*
  
  template class L2HighOrderFiniteElementFO<1>;
  template class L2HighOrderFiniteElementFO<2>;
  template class L2HighOrderFiniteElementFO<3>;

  template class T_L2HighOrderFiniteElementFO<ET_TRIG,0>;
  template class T_L2HighOrderFiniteElementFO<ET_TRIG,1>;
  template class T_L2HighOrderFiniteElementFO<ET_TRIG,2>;
  template class T_L2HighOrderFiniteElementFO<ET_TRIG,3>;
  template class T_L2HighOrderFiniteElementFO<ET_TRIG,4>;
  template class T_L2HighOrderFiniteElementFO<ET_TRIG,5>;
  template class T_L2HighOrderFiniteElementFO<ET_TRIG,6>;


  template class L2HighOrderFEFO<ET_TRIG,0>;
  template class L2HighOrderFEFO<ET_TRIG,1>;
  template class L2HighOrderFEFO<ET_TRIG,2>;
  template class L2HighOrderFEFO<ET_TRIG,3>;
  template class L2HighOrderFEFO<ET_TRIG,4>;
  template class L2HighOrderFEFO<ET_TRIG,5>; 
  template class L2HighOrderFEFO<ET_TRIG,6>;
  */


}
 
