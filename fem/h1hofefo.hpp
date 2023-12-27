#ifndef FILE_H1HOFEFO
#define FILE_H1HOFEFO

/*********************************************************************/
/* File:   h1hofefo.hpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   Apr. 2009                                                 */
/*********************************************************************/


#include <recursive_pol.hpp>

namespace ngfem
{

  /**
     High order finite elements for H^1 of fixed order
  */

  template <ELEMENT_TYPE ET, int ORDER> class H1HighOrderFEFO;

  /**
     High order triangular finite element
  */
  template <int ORDER>
  class H1HighOrderFEFO<ET_TRIG, ORDER> : 
    public T_ScalarFiniteElement< H1HighOrderFEFO<ET_TRIG,ORDER>, ET_TRIG >,
    public ET_trait<ET_TRIG> 
  {
    using ScalarFiniteElement<2>::ndof;
    using ScalarFiniteElement<2>::order;    

    using ET_trait<ET_TRIG> :: N_VERTEX;

    int vnums[N_VERTEX];

  public:
    enum { NDOF = (ORDER+1)*(ORDER+2)/2 };

    INLINE H1HighOrderFEFO () 
    {
      order = ORDER;
      ndof = NDOF; 
      for (int i = 0; i < N_VERTEX; i++) vnums[i] = i;
    }

    template <typename TA> 
    INLINE H1HighOrderFEFO<ET_TRIG, ORDER> * SetVertexNumbers (const TA & avnums)
    { 
      for (int i = 0; i < N_VERTEX; i++) vnums[i] = avnums[i]; 
      return this;
    }

    template<typename Tx, typename TFA>  
    INLINE void T_CalcShape (TIP<2,Tx> ip, TFA & shape) const; 
  };


  /**
     High order tet finite element
  */
  template <int ORDER>
  class H1HighOrderFEFO<ET_TET, ORDER> : 
    public T_ScalarFiniteElement< H1HighOrderFEFO<ET_TET,ORDER>, ET_TET>,
    public ET_trait<ET_TET> 
  {
    using ET_trait<ET_TET> :: N_VERTEX;

    int vnums[N_VERTEX];
    using ScalarFiniteElement<3>::ndof;
    using ScalarFiniteElement<3>::order;

    typedef IntegratedLegendreMonomialExt T_ORTHOPOL;

  public:
    enum { NDOF = (ORDER+1)*(ORDER+2)*(ORDER+3)/6 };
    NGS_DLL_HEADER H1HighOrderFEFO () 
    {
      order = ORDER; 
      ndof = NDOF; 
      for (int i = 0; i < N_VERTEX; i++) vnums[i] = i;
    }

    template <typename TA> 
    INLINE H1HighOrderFEFO<ET_TET, ORDER> * SetVertexNumbers (const TA & avnums)
    { 
      for (int i = 0; i < N_VERTEX; i++) vnums[i] = avnums[i]; 
      return this;
    }

    template<typename Tx, typename TFA>  
    INLINE void T_CalcShape (TIP<3,Tx> ip, TFA & shape) const; 
  };



}




#ifdef FILE_H1HOFEFO_CPP

#define H1HOFEFO_EXTERN
#include <h1hofefo_impl.hpp>
#include <tscalarfe_impl.hpp>

#else
#define H1HOFEFO_EXTERN extern
#endif


namespace ngfem
{

  H1HOFEFO_EXTERN template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TRIG,1>, ET_TRIG>;
  H1HOFEFO_EXTERN template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TRIG,2>, ET_TRIG>;
  H1HOFEFO_EXTERN template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TRIG,3>, ET_TRIG>;
  H1HOFEFO_EXTERN template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TRIG,4>, ET_TRIG>;
  H1HOFEFO_EXTERN template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TRIG,5>, ET_TRIG>;
  H1HOFEFO_EXTERN template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TRIG,6>, ET_TRIG>;

  H1HOFEFO_EXTERN template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TET,1>, ET_TET>;
  H1HOFEFO_EXTERN template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TET,2>, ET_TET>;
  H1HOFEFO_EXTERN template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TET,3>, ET_TET>;
  H1HOFEFO_EXTERN template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TET,4>, ET_TET>;
  H1HOFEFO_EXTERN template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TET,5>, ET_TET>;
  H1HOFEFO_EXTERN template class T_ScalarFiniteElement<H1HighOrderFEFO<ET_TET,6>, ET_TET>;

  H1HOFEFO_EXTERN template class H1HighOrderFEFO<ET_TRIG,1>;
  H1HOFEFO_EXTERN template class H1HighOrderFEFO<ET_TRIG,2>;
  H1HOFEFO_EXTERN template class H1HighOrderFEFO<ET_TRIG,3>;
  H1HOFEFO_EXTERN template class H1HighOrderFEFO<ET_TRIG,4>;
  H1HOFEFO_EXTERN template class H1HighOrderFEFO<ET_TRIG,5>;
  H1HOFEFO_EXTERN template class H1HighOrderFEFO<ET_TRIG,6>;

  H1HOFEFO_EXTERN template class H1HighOrderFEFO<ET_TET,1>;
  H1HOFEFO_EXTERN template class H1HighOrderFEFO<ET_TET,2>;
  H1HOFEFO_EXTERN template class H1HighOrderFEFO<ET_TET,3>;
  H1HOFEFO_EXTERN template class H1HighOrderFEFO<ET_TET,4>;
  H1HOFEFO_EXTERN template class H1HighOrderFEFO<ET_TET,5>;
  H1HOFEFO_EXTERN template class H1HighOrderFEFO<ET_TET,6>;


}


#endif
