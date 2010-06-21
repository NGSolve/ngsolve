#ifndef FILE_H1HOFEFO
#define FILE_H1HOFEFO

/*********************************************************************/
/* File:   h1hofefo.hpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   Apr. 2009                                                 */
/*********************************************************************/

#include "tscalarfe.hpp"



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
    void SetVertexNumbers (FlatArray<int> & avnums)
    {
      for (int i = 0; i < avnums.Size(); i++)
	vnums[i] = avnums[i];
    }

    void SetVertexNumber (int nr, int vnum) { vnums[nr] = vnum; }
  };


  template <ELEMENT_TYPE ET, int ORDER> class H1HighOrderFEFO;

  template <ELEMENT_TYPE ET, int ORDER>
  class T_H1HighOrderFiniteElementFO : public H1HighOrderFiniteElementFO<ET_trait<ET>::DIM>
  {
  protected:
    enum { DIM = ET_trait<ET>::DIM };

    using ScalarFiniteElement<DIM>::ndof;
    using ScalarFiniteElement<DIM>::order;
    using ScalarFiniteElement<DIM>::eltype;
    // using ScalarFiniteElement<DIM>::dimspace;

    using H1HighOrderFiniteElementFO<DIM>::vnums;

    typedef IntegratedLegendreMonomialExt T_ORTHOPOL;
    typedef TrigShapesInnerLegendre T_TRIGSHAPES;
    // typedef TrigShapesInnerJacobi T_TRIGSHAPES;

  public:

    T_H1HighOrderFiniteElementFO () 
    {
      for (int i = 0; i < ET_trait<ET>::N_VERTEX; i++)
	vnums[i] = i;
      // dimspace = DIM;
      eltype = ET;
      order = ORDER;
    }

    // virtual void ComputeNDof();
    // virtual void GetInternalDofs (Array<int> & idofs) const;

    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;

    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<DIM> dshape) const;

    virtual void CalcMappedDShape (const SpecificIntegrationPoint<DIM,DIM> & sip, 
				   FlatMatrixFixWidth<DIM> dshape) const;

    virtual double Evaluate (const IntegrationPoint & ip, 
			     FlatVector<double> x) const;
  };



  /**
     High order segment finite element
  */
  /*
    template <int ORDER> 
    class H1HighOrderFE<ET_SEGM> : public T_H1HighOrderFiniteElement<ET_SEGM>
    {
    public:
    H1HighOrderFE () { ; }
    H1HighOrderFE (int aorder);

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[1], TFA & shape) const; 
    };
  */

  /**
     High order triangular finite element
  */
  template <int ORDER>
  class H1HighOrderFEFO<ET_TRIG, ORDER> : public T_H1HighOrderFiniteElementFO<ET_TRIG, ORDER>
  {
    using T_H1HighOrderFiniteElementFO<ET_TRIG, ORDER>::ndof;
    using T_H1HighOrderFiniteElementFO<ET_TRIG, ORDER>::vnums; 

    typedef IntegratedLegendreMonomialExt T_ORTHOPOL;
    typedef TrigShapesInnerLegendre T_TRIGSHAPES;

  public:
    enum { NDOF = (ORDER+1)*(ORDER+2)/2 };
    H1HighOrderFEFO () { ndof = NDOF; }

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[2], TFA & shape) const; 
  };


  /**
     High order triangular finite element
  */
  template <int ORDER>
  class H1HighOrderFEFO<ET_TET, ORDER> : public T_H1HighOrderFiniteElementFO<ET_TET, ORDER>
  {
    using T_H1HighOrderFiniteElementFO<ET_TET, ORDER>::ndof;
    using T_H1HighOrderFiniteElementFO<ET_TET, ORDER>::vnums; 

    typedef IntegratedLegendreMonomialExt T_ORTHOPOL;

    typedef TetShapesInnerLegendre T_INNERSHAPES;
    typedef TetShapesFaceLegendre T_FACESHAPES;


  public:
    enum { NDOF = (ORDER+1)*(ORDER+2)*(ORDER+3)/6 };
    H1HighOrderFEFO () { ndof = NDOF; }

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[3], TFA & shape) const; 
  };



}


#endif
