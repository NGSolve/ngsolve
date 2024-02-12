#ifndef FILE_L2HOFEFO
#define FILE_L2HOFEFO

/*********************************************************************/
/* File:   l2hofefo.hpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   Apr. 2009                                                 */
/*********************************************************************/


// #include "tscalarfe_impl.hpp"
// #include "l2hofe_impl.hpp"







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

  class GenericOrientation;
  template <int V1, int V2, int V3, int V4=-1> class FixedOrientation;
  
  template <ELEMENT_TYPE ET, int ORDER, typename ORIENTATION = GenericOrientation>
  class L2HighOrderFEFO_Shapes;



  template <ELEMENT_TYPE ET, int ORDER,
            typename ORIENTATION = GenericOrientation,
            typename BASE = L2HighOrderFE<ET, L2HighOrderFEFO_Shapes<ET,ORDER,ORIENTATION> > >

  class L2HighOrderFEFO : public BASE
  {
  protected:
    // using typename BASE::T_IMPL;
    using typename BASE::T_SHAPES;
    typedef L2HighOrderFEFO_Shapes<ET,ORDER,ORIENTATION> SHAPES;


    enum { DIM = ET_trait<ET>::DIM };

    using ScalarFiniteElement<DIM>::ndof;
    using ScalarFiniteElement<DIM>::order;
    using DGFiniteElement<ET>::vnums;

  public:

    INLINE L2HighOrderFEFO () 
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

    HD virtual void Evaluate (const IntegrationRule & ir, BareSliceVector<double> coefs, FlatVector<double> vals) const
    {
      // static Timer t("evaluate");
      // RegionTimer r(t);
      // t.AddFlops (ir.GetNIP()*coefs.Size());

      /*
#ifndef __CUDA_ARCH__
      int classnr =  ET_trait<ET>::GetClassNr (vnums);
      PrecomputedScalShapes<DIM> * pre = SHAPES::precomp.Get (classnr, order, ir.GetNIP());
      if (pre)
        vals = FlatMatrixFixWidth<SHAPES::NDOF> (pre->shapes) * coefs;
      else
#endif
      */
        this -> BASE::T_IMPL::Evaluate (ir, coefs, vals);
    }

    HD virtual void EvaluateGradTrans (const IntegrationRule & ir, BareSliceMatrix<> values, BareSliceVector<> coefs) const
    {
      /*
        static Timer t("evaluate grad trans");
        RegionTimer r(t);
        t.AddFlops (DIM*ir.GetNIP()*coefs.Size());
      */
#ifndef __CUDA_ARCH__
      int classnr =  ET_trait<ET>::GetClassNr (vnums);

      PrecomputedScalShapes<DIM> * pre = SHAPES::precomp.Get (classnr, order, ir.GetNIP());
      if (pre)
	// coefs.Range(0,SHAPES::NDOF) = Trans (FlatMatrixFixWidth<SHAPES::NDOF> (pre->dshapes)) * FlatVector<> (DIM*SHAPES::NDOF, &values(0,0));
        coefs.Range(0,SHAPES::NDOF) = Trans (pre->dshapes) * FlatVector<> (DIM*SHAPES::NDOF, &values(0,0));        
      else
#endif
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

    HD NGS_DLL_HEADER virtual void GetDiagMassMatrix (FlatVector<> mass) const
    {
      if (ET == ET_SEGM)
	{
	  Iterate<ORDER+1> ([&] (int i) { mass[i] = 1.0 / (2*i+1); });
	}
      else if (ET == ET_TRIG)
	{
          int ii = 0;

          IterateTrig<ORDER> ([&] (int ix, int iy) LAMBDA_INLINE
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
      else if (ET == ET_TET)
        {
          /*
          int order = ORDER;
          for (int ix = 0, ii = 0; ix <= order; ix++)
            for (int iy = 0; iy <= order - ix; iy++)
              for (int iz = 0; iz <= order - ix-iy; iz++, ii++)
                mass(ii) = 1.0 / ((2 * ix + 1) * (2 * ix + 2 * iy + 2) * (2 * ix + 2 * iy + 2 * iz + 3));
          */
          int ii = 0;
          Iterate<ORDER+1>
            ([mass, &ii] (auto ix)
             {
               Iterate<ORDER+1-ix.value>
                 ([mass, ix, &ii] (auto iy)
                  {
                    Iterate<ORDER+1-ix.value-iy.value>
                      ([mass, ix, iy, &ii] (auto iz)
                       {
                         mass(ii) = 1.0 / ((2 * ix + 1) * (2 * ix + 2 * iy + 2) * (2 * ix + 2 * iy + 2 * iz + 3));
                         ii++;
                       });
                  });
             });
          
        }
#ifndef __CUDA_ARCH__
      else
	{
	  cerr << "L2HighOrderFEFO::getdiagmass not implemented" << endl;
	}
#endif
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
    using L2HighOrderFEFO<ET_SEGM, ORDER>::DIM; 

  public:

    enum { NDOF = (ORDER+1) };

    template<typename Tx, typename TFA>  
    INLINE void T_CalcShape (const TIP<DIM,Tx> & ip, TFA & shape) const
    {
      Tx lam[2] = { ip.x, 1-ip.x };
      IVec<2> e = this -> GetEdgeSort (0, vnums);
      // LegendrePolynomial_Old::EvalFO<ORDER> (lam[e[1]]-lam[e[0]], shape);
      LegendrePolynomial::EvalFO<ORDER> (lam[e[1]]-lam[e[0]], shape);
    }
  };




  /**
     High order triangular finite element
  */
  template <int ORDER>
  class L2HighOrderFEFO_Shapes<ET_TRIG, ORDER, GenericOrientation>
    : public L2HighOrderFEFO<ET_TRIG, ORDER, GenericOrientation>
  {
    using L2HighOrderFEFO<ET_TRIG, ORDER>::ndof;
    using L2HighOrderFEFO<ET_TRIG, ORDER>::vnums; 

  public:
    enum { NDOF = (ORDER+1)*(ORDER+2)/2 };

    template<typename Tx, typename TFA>  
    INLINE void T_CalcShape (const TIP<2,Tx> & ip, TFA & shape) const
    {
      Tx lam[3] = { ip.x, ip.y, 1-ip.x-ip.y };
      IVec<4> f = this -> GetFaceSort (0, vnums);

      Tx x = lam[f[0]];
      Tx y = lam[f[1]];

      int ii = 0;
      JacobiPolynomialAlpha jac(1);
      LegendrePolynomial::EvalScaled 
        (IC<ORDER>(), 
         y-(1-x-y), 1-x,
         SBLambda ([&] (auto i, Tx val) LAMBDA_INLINE 
                   {
                     // JacobiPolynomialFix<1+2*i,0> jac;
                     jac.EvalMult (IC<ORDER-i.value>(), 2*x-1, val, 
                                   SBLambda([&](auto j, Tx v2) 
                                            {
                                              shape[ii++] = v2;
                                            }));
                     jac.IncAlpha2();
                   }));
    }
  };


  template <int ORDER, int V1, int V2, int V3>
  class L2HighOrderFEFO_Shapes<ET_TRIG, ORDER, FixedOrientation<V1,V2,V3>>
    : public L2HighOrderFEFO<ET_TRIG, ORDER, FixedOrientation<V1,V2,V3>>
  {
    using L2HighOrderFEFO<ET_TRIG, ORDER, FixedOrientation<V1,V2,V3>>::ndof;
    using L2HighOrderFEFO<ET_TRIG, ORDER, FixedOrientation<V1,V2,V3>>::vnums;
    using L2HighOrderFEFO<ET_TRIG, ORDER, FixedOrientation<V1,V2,V3>>::DIM;     
    
  public:
    enum { NDOF = (ORDER+1)*(ORDER+2)/2 };

    template<typename Tx, typename TFA>  
      INLINE void T_CalcShape (const TIP<DIM,Tx> & ip, TFA & shape) const
    {
      Tx lam[3] = { ip.x, ip.y, 1-ip.x-ip.y };

      IVec<3> hvnums(V1,V2,V3);
      IVec<4> f = this -> GetFaceSort (0, hvnums);

      Tx x = lam[f[0]];
      Tx y = lam[f[1]];

      int ii = 0;
      JacobiPolynomialAlpha jac(1);
      LegendrePolynomial::EvalScaled 
        (IC<ORDER>(), 
         y-(1-x-y), 1-x,
         SBLambda ([&] (auto i, Tx val) LAMBDA_INLINE 
                   {
                     jac.EvalMult (IC<ORDER-i.value>(), 2*x-1, val, shape+ii);
                     ii += IC<ORDER-i.value+1>();
                     jac.IncAlpha2();
                   }));
    }
  };

  template <int ORDER, int V1, int V2, int V3, int V4>
  class L2HighOrderFEFO_Shapes<ET_TET, ORDER, FixedOrientation<V1,V2,V3,V4>>
    : public L2HighOrderFEFO<ET_TET, ORDER, FixedOrientation<V1,V2,V3,V4>>
  {
    using L2HighOrderFEFO<ET_TET, ORDER, FixedOrientation<V1,V2,V3,V4>>::ndof;
    using L2HighOrderFEFO<ET_TET, ORDER, FixedOrientation<V1,V2,V3,V4>>::vnums; 

  public:
    enum { NDOF = (ORDER+1)*(ORDER+2)*(ORDER+3)/6 };

    template<typename Tx, typename TFA>  
    INLINE void T_CalcShape (const TIP<3,Tx> & ip, TFA & shape) const
    {
//       Tx lami[4] = { ip.x, ip.y, ip.z, 1-ip.x-ip.y-ip.z };
// 
//       IVec<4> hvnums(V1,V2,V3,V4);
//       unsigned char sort[4] = { 0, 1, 2, 3 };
//       if (hvnums[sort[0]] > hvnums[sort[1]]) Swap (sort[0], sort[1]);
//       if (hvnums[sort[2]] > hvnums[sort[3]]) Swap (sort[2], sort[3]);
//       if (hvnums[sort[0]] > hvnums[sort[2]]) Swap (sort[0], sort[2]);
//       if (hvnums[sort[1]] > hvnums[sort[3]]) Swap (sort[1], sort[3]);
//       if (hvnums[sort[1]] > hvnums[sort[2]]) Swap (sort[1], sort[2]);
// 
//       Tx lamis[4];
//       for (int i = 0; i < 4; i++)
//         lamis[i] = lami[sort[i]];

      // Hack: only working because there are currently only 2 possibilities for V1-V4:
      // 0,1,2,3
      // 0,1,3,2
      const Tx l4 = 1-ip.x-ip.y-ip.z;
      auto const lamis = make_tuple( ip.x, ip.y, V3<V4 ? ip.z : l4, V3<V4 ? l4 : ip.z );
      const auto & lami0 = get<0>(lamis);
      const auto & lami1 = get<1>(lamis);
      const auto & lami2 = get<2>(lamis);
      const auto & lami3 = get<3>(lamis);

      size_t ii = 0;
      LegendrePolynomial leg;
      JacobiPolynomialAlpha jac1(1);    
      leg.EvalScaled 
        (IC<ORDER>(), lami2-lami3, lami2+lami3,
         SBLambda ([&](auto k, Tx polz) LAMBDA_INLINE
                   {
                     JacobiPolynomialAlpha jac2(2*k+2);
                     jac1.EvalScaledMult 
                       (IC<ORDER-k.value>(), lami1-lami2-lami3, 1-lami0, polz, 
                        SBLambda ([&] (auto j, Tx polsy) LAMBDA_INLINE
                                  {
                                    jac2.EvalMult(IC<ORDER-k.value-j.value>(), 2 * lami0 - 1, polsy, shape+ii);
                                    ii += IC<ORDER-k.value-j.value+1>();
                                    jac2.IncAlpha2();
                                  }));
                     jac1.IncAlpha2();
                   }));
      /*
      size_t ii = 0;
      LegendrePolynomial leg;
      // JacobiPolynomialAlpha jac1(1);    
      leg.EvalScaled 
        (IC<ORDER>(), lami[0], lami[0]+lami[1],
         SBLambda ([&](auto k, Tx polz) LAMBDA_INLINE
                   {
                     JacobiPolynomialFix<1+2*k,0> jac1;  
                     // JacobiPolynomialAlpha jac2(2*k+2);
                     jac1.EvalScaledMult 
                       (IC<ORDER-k>(), lami[1], lami[0]+lami[1]+lami[2], polz, 
                        SBLambda ([&] (auto j, Tx polsy) LAMBDA_INLINE
                                  {
                                    JacobiPolynomialFix<2+2*k+2*j,0> jac2;                                      
                                    jac2.EvalMult(IC<ORDER-k-j>(), lami[2], polsy, shape+ii);
                                    ii += IC<ORDER-k-j+1>();
                                    // jac2.IncAlpha2();
                                  }));
                     // jac1.IncAlpha2();
                   }));
      */
    }
  };
  


#ifdef FILE_L2HOFEFO_CPP
#define L2HOFEFO_EXTERN

#else
#define L2HOFEFO_EXTERN extern
#endif

  /*
  L2HOFEFO_EXTERN template class T_ScalarFiniteElement< L2HighOrderFEFO_Shapes<ET_TRIG,0>, ET_TRIG, DGFiniteElement<2>>;
  L2HOFEFO_EXTERN template class T_ScalarFiniteElement< L2HighOrderFEFO_Shapes<ET_TRIG,1>, ET_TRIG, DGFiniteElement<2>>;
  L2HOFEFO_EXTERN template class T_ScalarFiniteElement< L2HighOrderFEFO_Shapes<ET_TRIG,2>, ET_TRIG, DGFiniteElement<2>>;
  L2HOFEFO_EXTERN template class T_ScalarFiniteElement< L2HighOrderFEFO_Shapes<ET_TRIG,3>, ET_TRIG, DGFiniteElement<2>>;
  L2HOFEFO_EXTERN template class T_ScalarFiniteElement< L2HighOrderFEFO_Shapes<ET_TRIG,4>, ET_TRIG, DGFiniteElement<2>>;
  L2HOFEFO_EXTERN template class T_ScalarFiniteElement< L2HighOrderFEFO_Shapes<ET_TRIG,5>, ET_TRIG, DGFiniteElement<2>>;
  L2HOFEFO_EXTERN template class T_ScalarFiniteElement< L2HighOrderFEFO_Shapes<ET_TRIG,6>, ET_TRIG, DGFiniteElement<2>>;
  L2HOFEFO_EXTERN template class T_ScalarFiniteElement< L2HighOrderFEFO_Shapes<ET_TRIG,7>, ET_TRIG, DGFiniteElement<2>>;
  L2HOFEFO_EXTERN template class T_ScalarFiniteElement< L2HighOrderFEFO_Shapes<ET_TRIG,8>, ET_TRIG, DGFiniteElement<2>>;
  */
  
  /*
  L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,0>>;
  L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,1>>;
  L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,2>>;
  L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,3>>;
  L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,4>>;
  L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,5>>;
  L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,6>>;
  L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,7>>;
  L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,8>>;
  // L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,9>>;
  // L2HOFEFO_EXTERN template class L2HighOrderFE<ET_TRIG, L2HighOrderFEFO_Shapes<ET_TRIG,10>>;
  */

  /*
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,0>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,1>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,2>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,3>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,4>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,5>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,6>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,7>;
  L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,8>;
  */
  // L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,9>;
  // L2HOFEFO_EXTERN template class L2HighOrderFEFO<ET_TRIG,10>;


}


#endif

