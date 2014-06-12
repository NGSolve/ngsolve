#ifndef FILE_L2HOFEFO
#define FILE_L2HOFEFO

/*********************************************************************/
/* File:   l2hofefo.hpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   Apr. 2009                                                 */
/*********************************************************************/


// #include "tscalarfe_impl.hpp"
// #include "l2hofe_impl.hpp"




template <int NUM>
class Cl_Iterate
{
public:
  template <typename FUNC>
  static INLINE void Do (FUNC f)
  {
    Cl_Iterate<NUM-1>::Do(f);
    f(NUM);
  }
};

template <>
class Cl_Iterate<0>
{
public:
  template <typename FUNC>
  static INLINE void Do (FUNC f)  { f(0); }
};

template <int NUM, typename FUNC>
INLINE void Iterate (FUNC f)
{
  Cl_Iterate<NUM>::Do(f);
}








template <int O, int IX, int IY>
class Cl_IterateTrig
{
public:
  template <typename FUNC>
  static INLINE void Do (FUNC f)
  {
    Cl_IterateTrig<O,IX,IY-1>::Do(f);
    f(IX,IY);
  }
};

template <int O, int IX>
class Cl_IterateTrig<O,IX,0>
{
public:
  template <typename FUNC>
  static INLINE void Do (FUNC f)  
  {
    Cl_IterateTrig<O,IX-1,O-IX+1>::Do(f);    
    f(IX,0); 
  }
};

template <int O>
class Cl_IterateTrig<O,0,0>
{
public:
  template <typename FUNC>
  static INLINE void Do (FUNC f)  
  {
    f(0,0); 
  }
};









template <int O, int IX, int IY>
class Cl_IterateTrig3b
{
public:
  template <typename FUNC>
  static INLINE void Do (FUNC && f)
  {
    Cl_IterateTrig3b<O,IX,IY-1>::Do(f);
    f(IX,IY);
  }
};
template <int O, int IX>
class Cl_IterateTrig3b<O,IX,0>
{
public:
  template <typename FUNC>
  static INLINE void Do (FUNC && f)
  {
    f(IX,0);
  }
};


template <int O, int IX>
class Cl_IterateTrig3
{
public:
  template <typename FUNC>
  static INLINE void Do (FUNC && f)
  {
    Cl_IterateTrig3<O,IX-1>::Do(f);
    Cl_IterateTrig3b<O,IX,O-IX>::Do (f);
  }
};
template <int O>
class Cl_IterateTrig3<O,0>
{
public:
  template <typename FUNC>
  static INLINE void Do (FUNC && f)
  {
    Cl_IterateTrig3b<O,0,O>::Do (f);
  }
};












template <int NUM, typename FUNC>
INLINE void IterateTrig (FUNC && f)
{
  // Cl_IterateTrig<NUM,NUM,0>::Do(f);
  Cl_IterateTrig3<NUM,NUM>::Do(f);
}




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
      for (int i = 0; i < ET_trait<ET>::N_VERTEX; i++) vnums[i] = i;
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

    virtual void Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, FlatVector<double> vals) const
    {
      // static Timer t("evaluate");
      // RegionTimer r(t);
      // t.AddFlops (ir.GetNIP()*coefs.Size());

      int classnr =  ET_trait<ET>::GetClassNr (vnums);
      PrecomputedScalShapes<DIM> * pre = SHAPES::precomp.Get (classnr, order, ir.GetNIP());
      if (pre)
        vals = FlatMatrixFixWidth<SHAPES::NDOF> (pre->shapes) * coefs;
      else
        this -> BASE::T_IMPL::Evaluate (ir, coefs, vals);
    }

    virtual void EvaluateGradTrans (const IntegrationRule & ir, FlatMatrixFixWidth<DIM> values, FlatVector<> coefs) const
    {
      /*
        static Timer t("evaluate grad trans");
        RegionTimer r(t);
        t.AddFlops (DIM*ir.GetNIP()*coefs.Size());
      */

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
      if (ET == ET_SEGM)
	{
	  Iterate<ORDER> ([&] (int i) { mass[i] = 1.0 / (2*i+1); });
	}
      else if (ET == ET_TRIG)
	{
          int ii = 0;
          IterateTrig<ORDER> ([&] (int ix, int iy)
            {
              mass[ii] = 1.0 / ((2 * ix + 1) * (2 * ix + 2 * iy + 2));
              ii++;
            });

          /*
	  for (int ix = 0, ii = 0; ix <= ORDER; ix++)
	    for (int iy = 0; iy <= ORDER - ix; iy++, ii++)
	      mass(ii) = 1.0 / ((2 * iy + 1) * (2 * ix + 2 * iy + 2));
          */
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
  class L2HighOrderFEFO_Shapes<ET_SEGM, ORDER> : public L2HighOrderFEFO<ET_SEGM, ORDER>
  {
    using L2HighOrderFEFO<ET_SEGM, ORDER>::ndof;
    using L2HighOrderFEFO<ET_SEGM, ORDER>::vnums; 

  public:

    enum { NDOF = (ORDER+1) };

    template<typename Tx, typename TFA>  
    INLINE void T_CalcShape (Tx hx[2], TFA & shape) const
    {
      Tx lam[2] = { hx[0], 1-hx[0] };
      INT<2> e = this -> GetEdgeSort (0, vnums);
      LegendrePolynomial::EvalFO<ORDER> (lam[e[1]]-lam[e[0]], shape);
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
    INLINE void T_CalcShape (Tx hx[2], TFA & shape) const
    {
      Tx lam[3] = { hx[0], hx[1], 1-hx[0]-hx[1] };
      INT<4> f = this -> GetFaceSort (0, vnums);
      // DubinerBasis3::Eval (ORDER, lam[f[0]], lam[f[1]], shape);


      Tx x = lam[f[0]];
      Tx y = lam[f[1]];
      LegendrePolynomial_Old leg;
      Tx p1, p2, p3, p4;

      int ii = 0;
      IterateTrig<ORDER> ([&] (int ix, int iy) LAMBDA_INLINE
        {
          if (iy == 0)
            leg.EvalScaledNext (ix, y-(1-x-y), 1-x, p1, p2);

          JacobiPolynomialNew jac(1+2*ix, 0);
          shape[ii] = jac.EvalNextMult (iy, 2*x-1, p1, p3, p4);
          ii++;
        });
    }

  };


#ifdef NONE

#ifdef FILE_L2HOFEFO_CPP
#define L2HOFEFO_EXTERN

#else
#define L2HOFEFO_EXTERN extern
#endif

  L2HOFEFO_EXTERN template class T_ScalarFiniteElement< L2HighOrderFEFO_Shapes<ET_TRIG,0>, ET_TRIG, DGFiniteElement<2>>;
  L2HOFEFO_EXTERN template class T_ScalarFiniteElement< L2HighOrderFEFO_Shapes<ET_TRIG,1>, ET_TRIG, DGFiniteElement<2>>;
  L2HOFEFO_EXTERN template class T_ScalarFiniteElement< L2HighOrderFEFO_Shapes<ET_TRIG,2>, ET_TRIG, DGFiniteElement<2>>;
  L2HOFEFO_EXTERN template class T_ScalarFiniteElement< L2HighOrderFEFO_Shapes<ET_TRIG,3>, ET_TRIG, DGFiniteElement<2>>;
  L2HOFEFO_EXTERN template class T_ScalarFiniteElement< L2HighOrderFEFO_Shapes<ET_TRIG,4>, ET_TRIG, DGFiniteElement<2>>;
  L2HOFEFO_EXTERN template class T_ScalarFiniteElement< L2HighOrderFEFO_Shapes<ET_TRIG,5>, ET_TRIG, DGFiniteElement<2>>;
  L2HOFEFO_EXTERN template class T_ScalarFiniteElement< L2HighOrderFEFO_Shapes<ET_TRIG,6>, ET_TRIG, DGFiniteElement<2>>;
  /*
  L2HOFEFO_EXTERN template class T_ScalarFiniteElement< L2HighOrderFEFO_Shapes<ET_TRIG,7>, ET_TRIG, DGFiniteElement<2>>;
  L2HOFEFO_EXTERN template class T_ScalarFiniteElement< L2HighOrderFEFO_Shapes<ET_TRIG,8>, ET_TRIG, DGFiniteElement<2>>;
  */
  

  L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,0>>;
  L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,1>>;
  L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,2>>;
  L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,3>>;
  L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,4>>;
  L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,5>>;
  L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,6>>;
  /*
  L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,7>>;
  L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,8>>;
  // L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,9>>;
  // L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,10>>;
  */

  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,0>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,1>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,2>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,3>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,4>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,5>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,6>;
  /*
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,7>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,8>;
  */
  // L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,9>;
  // L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,10>;

#endif

}


#endif

