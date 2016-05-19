#ifdef __CUDA_ARCH__
#include "tscalarfe_impl_cuda.hpp"
#endif


#ifndef FILE_TSCALARFE_IMPL
#define FILE_TSCALARFE_IMPL

 
namespace ngfem
{

  template <int DIM> 
  INLINE Vec<DIM, AutoDiff<DIM>> Ip2Ad (const IntegrationPoint & ip)
  {
    Vec<DIM, AutoDiff<DIM> > adp;
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);
    return adp;
  }

  




  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  CalcShape (const IntegrationPoint & ip, SliceVector<> shape) const
  {
    Vec<DIM> pt = ip.Point();
    T_CalcShape (&pt(0), shape);
  }

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  CalcDShape (const IntegrationPoint & ip, 
              SliceMatrix<> dshape) const
  {
    Vec<DIM, AutoDiff<DIM> > adp = ip; 
    T_CalcShape (&adp(0), SBLambda ([&] (int i, AutoDiff<DIM> shape)
                                    { shape.StoreGradient (&dshape(i,0)) ; }));
  }

#ifndef FASTCOMPILE

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  CalcShape (const IntegrationRule & ir, SliceMatrix<> shape) const
  {
    for (int i = 0; i < ir.Size(); i++)
      {
	Vec<DIM> pt = ir[i].Point();
	T_CalcShape (&pt(0), shape.Col(i));
      }
  }


  template <class FEL, ELEMENT_TYPE ET, class BASE>
  double T_ScalarFiniteElement<FEL,ET,BASE> :: 
  Evaluate (const IntegrationPoint & ip, SliceVector<double> x) const
  {
    Vec<DIM> pt = ip.Point();

    double sum = 0;
    T_CalcShape (&pt(0), SBLambda ( [&](int i, double val) { sum += x(i)*val; } ));
    return sum;
  }  


  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  Evaluate (const IntegrationRule & ir, SliceVector<double> coefs, FlatVector<double> vals) const
  {
    for (int i = 0; i < ir.GetNIP(); i++)
      {
        Vec<DIM> pt = ir[i].Point();

        double sum = 0;
        T_CalcShape (&pt(0), SBLambda ( [&](int i, double shape) { sum += coefs(i)*shape; } ));
        vals(i) = sum;
      }
  }

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  Evaluate (const SIMD_IntegrationRule & ir, BareSliceVector<> coefs, ABareVector<double> values) const
  {
    // static Timer t("ScalarFE::Evaluate", 2); RegionTimer reg(t);
    // t.AddFlops (ir.GetNIP()*ndof);
    
    FlatArray<SIMD<IntegrationPoint>> hir = ir;
    for (int i = 0; i < hir.Size(); i++)
      {
        Vec<DIM,SIMD<double>> pt = hir[i];
        SIMD<double> sum = 0;
        T_CalcShape (&pt(0), SBLambda ( [&](int j, SIMD<double> shape) { sum += coefs(j)*shape; } ));
        values.Get(i) = sum.Data();
      }
  }
    

  

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  Evaluate (const IntegrationRule & ir, SliceMatrix<> coefs, SliceMatrix<> values) const
  {
    for (int i = 0; i < ir.GetNIP(); i++)
      {
        Vec<DIM> pt = ir[i].Point();

        values.Row(i) = 0.0;
        T_CalcShape (&pt(0), SBLambda ( [&](int j, double shape) 
                                        { 
                                          // sum += coefs(i)*shape; 
                                          values.Row(i) += shape * coefs.Row(j); 
                                        } 
                                        ));
      }
  }

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  EvaluateTrans (const IntegrationRule & ir, FlatVector<> vals, SliceVector<double> coefs) const
  {
    coefs = 0.0;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
        Vec<DIM> pt = ir[i].Point();
        T_CalcShape (&pt(0), SBLambda ( [&](int j, double shape) 
                                        { coefs(j) += vals(i)*shape; } ));
      }
  }

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  AddTrans (const SIMD_IntegrationRule & ir, ABareVector<double> values,
            BareSliceVector<> coefs) const
  {
    for (int i = 0; i < ir.Size(); i++)
      {
        Vec<DIM,SIMD<double>> pt = ir[i];
        /*
        for (int j = 0; j < DIM; j++)
          pt(j) = ir[i](j);
        */
        SIMD<double> val = values.Get(i);
        T_CalcShape (&pt(0), SBLambda ( [&](int j, SIMD<double> shape) { coefs(j) += HSum(val*shape); } ));
      }
  }
    

  
  template <class FEL, ELEMENT_TYPE ET, class BASE>
  auto T_ScalarFiniteElement<FEL,ET,BASE> :: 
  EvaluateGrad (const IntegrationPoint & ip, SliceVector<double> coefs) const -> Vec<DIM>
  {
    Vec<DIM, AutoDiff<DIM>> adp = ip;
    AutoDiff<DIM> sum = 0.0;
    T_CalcShape (&adp(0), SBLambda ( [&](int i, AutoDiff<DIM> val) 
                                     { 
                                       sum += coefs(i) * val;
                                     }));
    return AD2Vec<DIM> (sum);
  }

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  EvaluateGrad (const IntegrationRule & ir, SliceVector<double> coefs, 
                FlatMatrixFixWidth<DIM> vals) const
  {
    for (int i = 0; i < ir.GetNIP(); i++)
      {
        Vec<DIM, AutoDiff<DIM>> adp = ir[i]; 

        Vec<DIM> sum = 0.0;
        T_CalcShape (&adp(0), SBLambda ([&] (int j, AD2Vec<DIM> shape)
                                        { sum += coefs(j) * shape; }));
        vals.Row(i) = sum;
      }
  }


  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  EvaluateGrad (const SIMD_BaseMappedIntegrationRule & bmir,
                BareSliceVector<> coefs,
                ABareMatrix<double> values) const
  {
    auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
    for (int i = 0; i < mir.Size(); i++)
      {
        Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
        Vec<DIM,SIMD<double>> sum(0.0);
        T_CalcShape (&adp(0), SBLambda ([&] (int j, AD2Vec<DIM,SIMD<double>> shape)
                                        { sum += coefs(j) * shape; }));
        for (int k = 0; k < DIM; k++)
          values.Get(k,i) = sum(k).Data();
      }
  }


  
  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  EvaluateGradTrans (const IntegrationRule & ir, 
                     FlatMatrixFixWidth<DIM> vals, SliceVector<double> coefs) const
  {
    coefs = 0.0;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
        Vec<DIM, AutoDiff<DIM>> adp = ir[i];
        T_CalcShape (&adp(0), SBLambda ([&] (int j, AD2Vec<DIM> shape)
                                        { coefs(j) += InnerProduct (vals.Row(i), shape); }));
      }
  }
  
  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  EvaluateGradTrans (const IntegrationRule & ir, SliceMatrix<> values, SliceMatrix<> coefs) const
  {
    int nels = coefs.Width();
    coefs = 0.0;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
        Vec<DIM, AutoDiff<DIM>> adp = ir[i];  
        T_CalcShape (&adp(0), SBLambda ([&] (int j, AD2Vec<DIM> shape)
                                        { 
                                          Vec<DIM> grad = shape;
                                          FlatMatrixFixWidth<DIM> mvals(nels, &values(i,0));
                                          coefs.Row(j) += mvals * grad;
                                        }));
      }
  }

  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  AddGradTrans (const SIMD_BaseMappedIntegrationRule & bmir,
                ABareMatrix<double> values,
                BareSliceVector<> coefs) const
  {
    auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIM>&> (bmir);
    for (int i = 0; i < mir.Size(); i++)
      {
        Vec<DIM, AutoDiff<DIM,SIMD<double>>> adp = mir[i];
        T_CalcShape (&adp(0), SBLambda ([&] (int j, AD2Vec<DIM,SIMD<double>> shape)
                                        {
                                          SIMD<double> sum = 0.0;
                                          for (int k = 0; k < DIM; k++)
                                            sum += shape(k) * values.Get(k,i);
                                          coefs(j) += HSum(sum);
                                        }));
      }
  }
  
  /*
  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  CalcDShape (const IntegrationPoint & ip, 
	      const std::function<void(int,Vec<DIM>)> & callback) const
  {
    Vec<DIM, AutoDiff<DIM> > adp;
    for (int i = 0; i < DIM; i++)
      adp[i] = AutoDiff<DIM> (ip(i), i);
      
    // DShapeAssign<DIM> ds(dshape); 
    // T_CalcShape (&adp(0), ds);


    T_CalcShape (&adp(0), SBLambda ([&] (int i, AutoDiff<DIM> shape)
                                    {
				      Vec<DIM> v;
				      shape.StoreGradient (&v(0));
				      callback (i,v);
				    }));
  }
  */



  /*
  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  CalcMappedDShape (const MappedIntegrationPoint<DIM,DIM> & mip, 
		    FlatMatrixFixWidth<DIM> dshape) const
  {
    Vec<DIM, AutoDiff<DIM> > adp;   
    for (int i = 0; i < DIM; i++)
      adp[i].Value() = mip.IP()(i);
      
    for (int i = 0; i < DIM; i++)
      for (int j = 0; j < DIM; j++)
	adp[i].DValue(j) = mip.GetJacobianInverse()(i,j);

    T_CalcShape (&adp(0), SBLambda ([&] (int i, AutoDiff<DIM> shape)
                                    { shape.StoreGradient (&dshape(i,0)) ; }));
  }
  */


  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  CalcMappedDShape (const MappedIntegrationPoint<DIM,DIM> & mip, 
		    SliceMatrix<> dshape) const
  {
    Vec<DIM, AutoDiff<DIM>> adp = mip;

    T_CalcShape (&adp(0), SBLambda ([&] (int i, AutoDiff<DIM> shape)
                                    { shape.StoreGradient (&dshape(i,0)) ; }));
  }


  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  CalcMappedDShape (const MappedIntegrationRule<DIM,DIM> & mir, 
		    SliceMatrix<> dshape) const
  {
    for (int i = 0; i < mir.Size(); i++)
      T_ScalarFiniteElement::CalcMappedDShape (mir[i], dshape.Cols(i*DIM,(i+1)*DIM));
  }


  /*
    ... not yet working
  template <class FEL, ELEMENT_TYPE ET, class BASE>
  void T_ScalarFiniteElement<FEL,ET,BASE> :: 
  GetPolOrders (FlatArray<PolOrder<DIM> > orders) const
  {
    Vec<DIM,PolOrder<DIM>> po;

    switch (ET)
      {
      case ET_TRIG:
        po[0] = INT<DIM> (1,1); 
        po[1] = INT<DIM> (1,1); 
        break;
      case ET_QUAD:
        po[0] = INT<DIM> (1,0); 
        po[1] = INT<DIM> (0,1); 
        break;
      case ET_TET:
        po[0] = INT<DIM> (1,1,1); 
        po[1] = INT<DIM> (1,1,1); 
        po[2] = INT<DIM> (1,1,1); 
        break;
      case ET_PRISM:
        po[0] = INT<DIM> (1,1,0); 
        po[1] = INT<DIM> (1,1,0); 
        po[2] = INT<DIM> (0,0,1); 
        break;

      default:
        for (int i = 0; i < DIM; i++)
          for (int j = 0; j < DIM; j++)
            po[i](j) = 1;
      }

    T_CalcShape (&po[0], orders);
    // did not work for old tensor productelements: order cancelation for lam_e
  }
  */


#endif

}



#endif
