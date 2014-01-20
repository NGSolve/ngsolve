#ifndef FILE_L2HOFEFO
#define FILE_L2HOFEFO

/*********************************************************************/
/* File:   l2hofefo.hpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   Apr. 2009                                                 */
/*********************************************************************/


// #include "tscalarfe_impl.hpp"
// #include "l2hofe_impl.hpp"


namespace ngfem
{

  /**
     High order finite elements for L2 of fixed order
  */


  template <ELEMENT_TYPE ET, int ORDER> class L2HighOrderFEFO_Shapes;



  template <ELEMENT_TYPE ET, int ORDER,
            typename BASE = L2HighOrderFE<ET, L2HighOrderFEFO_Shapes<ET,ORDER> > >

  class L2HighOrderFEFO : public BASE
  {
  protected:
    // using typename BASE::T_IMPL;
    using typename BASE::T_SHAPES;
    typedef L2HighOrderFEFO_Shapes<ET,ORDER> SHAPES;


    enum { DIM = ET_trait<ET>::DIM };

    using ScalarFiniteElement<DIM>::ndof;
    using ScalarFiniteElement<DIM>::order;
    using DGFiniteElement<DIM>::vnums;

  public:

    L2HighOrderFEFO () 
    {
      for (int i = 0; i < ET_trait<ET>::N_VERTEX; i++)
	vnums[i] = i;
      order = ORDER;
      ndof = SHAPES::NDOF;
    }

    /*
    virtual void PrecomputeShapes (const IntegrationRule & ir) 
    {
      int classnr =  ET_trait<ET>::GetClassNr (vnums);

      PrecomputedScalShapes<DIM> * pre = new  PrecomputedScalShapes<DIM> (ir.GetNIP(), ndof);

      Mat<SHAPES::NDOF, DIM> dshapes;
      for (int i = 0; i < ir.GetNIP(); i++)
	{
	  this->CalcShape (ir[i], pre->shapes.Row(i));
	  this->CalcDShape (ir[i], FlatMatrix<> (dshapes));
	  pre->dshapes.Rows (DIM*i, DIM*(i+1)) = Trans (dshapes);
	}

      SHAPES::precomp.Add (classnr, order, ir.GetNIP(), pre);
    }
    */

    void Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, FlatVector<double> vals) const
    {
      int classnr =  ET_trait<ET>::GetClassNr (vnums);

      PrecomputedScalShapes<DIM> * pre = SHAPES::precomp.Get (classnr, order, ir.GetNIP());
      if (pre)
	vals = FlatMatrixFixWidth<SHAPES::NDOF> (pre->shapes) * coefs;
      else
        this -> BASE::T_IMPL::Evaluate (ir, coefs, vals);
    }

    virtual void EvaluateGradTrans (const IntegrationRule & ir, FlatMatrixFixWidth<DIM> values, FlatVector<> coefs) const
    {
      int classnr =  ET_trait<ET>::GetClassNr (vnums);

      PrecomputedScalShapes<DIM> * pre = SHAPES::precomp.Get (classnr, order, ir.GetNIP());
      if (pre)
	coefs = Trans (FlatMatrixFixWidth<SHAPES::NDOF> (pre->dshapes)) * FlatVector<> (DIM*SHAPES::NDOF, &values(0,0));
      else
        BASE::T_IMPL:: EvaluateGradTrans (ir, values, coefs);
    }

    /*
    NGS_DLL_HEADER virtual void GetTrace (int facet, FlatVector<> coefs, FlatVector<> fcoefs) const 
    { 
      cout << "L2HighOrderFEFO::GetTrace not available" << endl;
    }

    NGS_DLL_HEADER virtual void GetTraceTrans (int facet, FlatVector<> fcoefs, FlatVector<> coefs) const
    { 
      cout << "L2HighOrderFEFO::GetTraceTrans not available" << endl;
    }
    */

    NGS_DLL_HEADER virtual void GetDiagMassMatrix (FlatVector<> mass) const
    {
      if (ET == ET_TRIG)
	{
	  for (int ix = 0, ii = 0; ix <= ORDER; ix++)
	    for (int iy = 0; iy <= ORDER - ix; iy++, ii++)
	      mass(ii) = 1.0 / ((2 * iy + 1) * (2 * ix + 2 * iy + 2));
	}
      else
	{
	  cerr << "L2HighOrderFEFO::getdiagmass not implemented" << endl;
	}
    }

  };




  /**
     High order triangular finite element
  */
  template <int ORDER>
  class L2HighOrderFEFO_Shapes<ET_TRIG, ORDER> : public L2HighOrderFEFO<ET_TRIG, ORDER>
  {
    using L2HighOrderFEFO<ET_TRIG, ORDER>::ndof;
    using L2HighOrderFEFO<ET_TRIG, ORDER>::vnums; 

  public:

    enum { NDOF = (ORDER+1)*(ORDER+2)/2 };

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[2], TFA & shape) const
    {
      Tx lam[3] = { x[0], x[1], 1-x[0]-x[1] };
      INT<4> f = this -> GetFaceSort (0, vnums);
      DubinerBasis::Eval (ORDER, lam[f[0]], lam[f[1]], shape);
      
      /*
      Tx lami[3] = { hx[0], hx[1], 1-hx[0]-hx[1] };
      INT<4> f = this->GetFaceSort (0, vnums);
      Tx x = lami[f[0]],  y = lami[f[1]],  l3 = lami[f[2]];

      Vec<ORDER+1, Tx> polx;
      Mat<ORDER+1, ORDER+1, Tx> polsy;

      LegendrePolynomial leg;
      leg.EvalScaledFO<ORDER> (x-l3, 1-y, polx);

      DubinerJacobiPolynomialsFO<ORDER, ORDER, 1,0>::Eval (2*y-1, polsy);
      
      for (int i = 0, ii = 0; i <= ORDER; i++)
	for (int j = 0; j <= ORDER-i; j++)
	  shape[ii++] = polx[i] * polsy(i,j);
      */
    }

  };




#ifdef FILE_L2HOFEFO_CPP
#define L2HOFEFO_EXTERN

#else
#define L2HOFEFO_EXTERN extern
#endif

  
  L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,0>>;
  L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,1>>;

  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,0>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,1>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,2>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,3>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,4>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,5>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,6>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,7>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,8>;

}


#endif
