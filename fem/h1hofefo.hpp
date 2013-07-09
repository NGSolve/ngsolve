#ifndef FILE_H1HOFEFO
#define FILE_H1HOFEFO

/*********************************************************************/
/* File:   h1hofefo.hpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   Apr. 2009                                                 */
/*********************************************************************/


namespace ngfem
{

  /**
     High order finite elements for H^1 of fixed order
  */

  template<int DIM>
  class H1HighOrderFiniteElementFO : virtual public ScalarFiniteElement<DIM>
  {
  public:
    int vnums[8];

  public:
    template <typename TA> 
    void SetVertexNumbers (const TA & avnums)
    { for (int i = 0; i < avnums.Size(); i++) vnums[i] = avnums[i]; }

    void SetVertexNumber (int nr, int vnum) { vnums[nr] = vnum; }
  };


  template <ELEMENT_TYPE ET, int ORDER> class H1HighOrderFEFO;

  /*
  template <ELEMENT_TYPE ET>
  class T_H1HighOrderFiniteElementFO : 
    public H1HighOrderFiniteElementFO<ET_trait<ET>::DIM>
  {
  protected:
    enum { DIM = ET_trait<ET>::DIM };

    using ScalarFiniteElement<DIM>::eltype;
    using H1HighOrderFiniteElementFO<DIM>::vnums;

  public:

    T_H1HighOrderFiniteElementFO () 
    {
    }

  };
  */


  /**
     High order triangular finite element
  */
  template <int ORDER>
  class H1HighOrderFEFO<ET_TRIG, ORDER> : 
    public H1HighOrderFiniteElementFO<2>,
    public T_ScalarFiniteElement< H1HighOrderFEFO<ET_TRIG,ORDER>, ET_TRIG >,
    public ET_trait<ET_TET> 
  {
    using ScalarFiniteElement<2>::ndof;
    using H1HighOrderFiniteElementFO<2>::vnums; 

    typedef IntegratedLegendreMonomialExt T_ORTHOPOL;
    typedef TrigShapesInnerLegendre T_TRIGSHAPES;

  public:
    enum { NDOF = (ORDER+1)*(ORDER+2)/2 };

    H1HighOrderFEFO () 
    {
      order = ORDER;
      ndof = NDOF; 
      eltype = ET_TRIG;
      for (int i = 0; i < ET_trait<ET_TRIG>::N_VERTEX; i++) 
        vnums[i] = i;
    }

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[2], TFA & shape) const; 
  };


  /**
     High order tet finite element
  */
  template <int ORDER>
  class H1HighOrderFEFO<ET_TET, ORDER> : 
    public H1HighOrderFiniteElementFO<3>,
    public T_ScalarFiniteElement< H1HighOrderFEFO<ET_TET,ORDER>, ET_TET>,
    public ET_trait<ET_TET> 
  {
    using H1HighOrderFiniteElementFO<3>::vnums; 
    using ScalarFiniteElement<3>::ndof;

    typedef IntegratedLegendreMonomialExt T_ORTHOPOL;

    typedef TetShapesInnerLegendre T_INNERSHAPES;
    typedef TetShapesFaceLegendre T_FACESHAPES;


  public:
    enum { NDOF = (ORDER+1)*(ORDER+2)*(ORDER+3)/6 };
    H1HighOrderFEFO () 
    {
      order = ORDER; 
      ndof = NDOF; 
      eltype = ET_TET;
      for (int i = 0; i < ET_trait<ET_TET>::N_VERTEX; i++)
	vnums[i] = i;
    }

    // ~H1HighOrderFEFO ();

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[3], TFA & shape) const; 
  };

}


#endif
