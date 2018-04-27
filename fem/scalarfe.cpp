/*********************************************************************/
/* File:   scalarfe.cpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


/* 
   Finite Element Definitions
*/

#define FILE_SCALARFE_CPP


#include <fem.hpp>
#include "h1lofe.hpp"
#include "l2hofe.hpp"

namespace ngfem
{

  /*
  template <int D>
  ScalarFiniteElement<D> :: ~ScalarFiniteElement () 
  { ; }
  */

  string BaseScalarFiniteElement :: ClassName() const 
  {
    return "ScalarFiniteElement"; 
  }

  /*
  template <int D>
  void ScalarFiniteElement<D> :: 
  GetPolOrders (FlatArray<PolOrder<D> > orders) const
  {
#ifndef __CUDA_ARCH__
    throw Exception (string ("GetPolOrders not implemented for element") + ClassName());
#endif
  }
  */

 
  /*
    a ( eps - (-eps) ) + b ( 2 eps - (-2eps) )  = 1
    a ( eps^3 - (-eps)^3) + b ( 8 eps^3 - -8 eps^3 ) = 0
    
    2 a + 4 b = 1 / eps
    2 a + 16 b = 0  

    b = -1 / 12 eps
    a = 2 / 3 eps
  */

  /*
  /// compute dshape, matrix: ndof x spacedim
  template<int D>
  void ScalarFiniteElement<D> :: 
  CalcMappedDShape (const MappedIntegrationPoint<D,D> & mip, 
                    FlatMatrixFixWidth<D> dshape) const
  {
    CalcDShape (mip.IP(), dshape);
    for (int i = 0; i < dshape.Height(); i++)
      {
        Vec<D> hv = dshape.Row(i);
        dshape.Row(i) = Trans (mip.GetJacobianInverse ()) * hv;
      }
  }
  */

  void BaseScalarFiniteElement ::
  CalcShape (const IntegrationPoint & ip, 
             BareSliceVector<Complex> shape) const
  {
    CalcShape (ip, SliceVector<double> (ndof, 2*shape.Dist(), reinterpret_cast<double*> (&shape(0))));
    SliceVector<double> imag_part(ndof, 2*shape.Dist(), reinterpret_cast<double*> (&shape(0))+1);
    imag_part = 0.0;
  }


  
  void BaseScalarFiniteElement :: 
  CalcShape (const IntegrationRule & ir, 
	     BareSliceMatrix<> shape) const
  {
    for (int i = 0; i < ir.Size(); i++)
      CalcShape (ir[i], shape.Col(i));
  }

  void BaseScalarFiniteElement :: 
  CalcShape (const SIMD_IntegrationRule & ir, 
             BareSliceMatrix<SIMD<double>> shape) const
  {
    throw ExceptionNOSIMD("SIMD - CalcShape not overloaded");
  }


  template<int D>
  void ScalarFiniteElement<D> :: 
  CalcMappedDShape (const BaseMappedIntegrationPoint & bmip, 
                    BareSliceMatrix<> dshape) const
  {
    auto & mip = static_cast<const MappedIntegrationPoint<D,D> &> (bmip);
    CalcDShape (mip.IP(), dshape);
    for (size_t i = 0; i < dshape.Height(); i++)
      {
        Vec<D> hv = dshape.Row(i);
        FlatVec<D> (&dshape(i,0)) = Trans (mip.GetJacobianInverse ()) * hv;
      }
  }



  template<int D>
  void ScalarFiniteElement<D> :: 
  CalcMappedDShape (const BaseMappedIntegrationRule & bmir, 
                    BareSliceMatrix<> dshapes) const
  {
    auto & mir = static_cast<const MappedIntegrationRule<D,D> &> (bmir);    
    for (size_t i = 0; i < mir.Size(); i++)
      CalcMappedDShape (mir[i], dshapes.Cols(i*D,(i+1)*D));
  }

  void BaseScalarFiniteElement :: 
  CalcMappedDShape (const SIMD_BaseMappedIntegrationRule & mir, 
                    BareSliceMatrix<SIMD<double>> dshapes) const
  {
    throw ExceptionNOSIMD("SIMD - CalcDShape not overloaded");    
  }

 /*
  template<int D>
  void ScalarFiniteElement<D> :: 
  CalcDShape (const IntegrationPoint & ip, 
	      const std::function<void(int,Vec<D>)> & callback) const
  {
    cout << "ScalarFE<D>::CalcMappedDShape(callback) not implemented" << endl;
  }  
  */

  double BaseScalarFiniteElement :: 
  Evaluate (const IntegrationPoint & ip, BareSliceVector<double> x) const
  {
    VectorMem<20, double> shape(ndof);
    CalcShape (ip, shape);
    return InnerProduct (shape, x);
  }  
  
  template<int D>
  Vec<D> ScalarFiniteElement<D> :: 
  EvaluateGrad (const IntegrationPoint & ip, BareSliceVector<double> x) const
  {
    MatrixFixWidth<D> dshape(ndof);
    CalcDShape (ip, dshape);
    Vec<D> grad = Trans (dshape) * x;
    return grad;
  }  


  void BaseScalarFiniteElement :: 
  Evaluate (const IntegrationRule & ir, BareSliceVector<double> coefs, BareSliceVector<double> vals) const
  {
    for (size_t i = 0; i < ir.GetNIP(); i++)
      vals(i) = Evaluate (ir[i], coefs);
  }

  void BaseScalarFiniteElement :: 
  Evaluate (const SIMD_IntegrationRule & ir, BareSliceVector<> coefs, BareVector<SIMD<double>> values) const
  {
    throw ExceptionNOSIMD (string("Evaluate (simd) not implemented for class ")+typeid(*this).name());
  }

  void BaseScalarFiniteElement :: 
  Evaluate (const SIMD_IntegrationRule & ir, SliceMatrix<> coefs, BareSliceMatrix<SIMD<double>> values) const
  {
    for (size_t i = 0; i < coefs.Width(); i++)
      Evaluate (ir, coefs.Col(i), values.Row(i));
  }

  
  void BaseScalarFiniteElement :: 
  Evaluate (const IntegrationRule & ir, SliceMatrix<> coefs, BareSliceMatrix<> values) const
  {
    VectorMem<100> shapes(coefs.Height());
    for (size_t i = 0; i < ir.Size(); i++)
      {
        CalcShape (ir[i], shapes);
        values.Row(i).AddSize(coefs.Width()) = Trans(coefs) * shapes;
      }
  }


  template<int D>
  void ScalarFiniteElement<D> :: 
  EvaluateGrad (const IntegrationRule & ir, BareSliceVector<double> coefs, BareSliceMatrix<> vals) const
  {
    for (size_t i = 0; i < ir.GetNIP(); i++)
      vals.Row(i).AddSize(D) = EvaluateGrad (ir[i], coefs);
  }

  void BaseScalarFiniteElement :: 
  EvaluateGrad (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const
  {
    throw ExceptionNOSIMD (string("EvaluateGrad (simd) not implemented for class ")+typeid(*this).name());
  }

  void BaseScalarFiniteElement :: 
  EvaluateGrad (const SIMD_IntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const
  {
    throw ExceptionNOSIMD (string("EvaluateGrad (simd) not implemented for class ")+typeid(*this).name());
  }

  
  void BaseScalarFiniteElement :: 
  EvaluateTrans (const IntegrationRule & ir, FlatVector<double> vals, BareSliceVector<double> coefs) const
  {
    VectorMem<20, double> shape(ndof);
    coefs.AddSize(ndof) = 0.0;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	CalcShape (ir[i], shape);
	coefs.AddSize(ndof) += vals(i) * shape;
      }
  }

  void BaseScalarFiniteElement :: 
  AddTrans (const SIMD_IntegrationRule & ir, BareVector<SIMD<double>> values, BareSliceVector<> coefs) const
  {
    throw ExceptionNOSIMD (string("AddTrans (simd) not implemented for class ")+typeid(*this).name());    
  }

  void BaseScalarFiniteElement :: 
  AddTrans (const SIMD_IntegrationRule & ir, BareSliceMatrix<SIMD<double>> values, SliceMatrix<> coefs) const
  {
    for (int i = 0; i < coefs.Width(); i++)
      AddTrans (ir, values.Row(i), coefs.Col(i));
  }


  void BaseScalarFiniteElement ::
  CalcDualShape (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
  {
    throw Exception (string("CalcDualShape not overloaded for element ") + typeid(*this).name());
  }

  
  template<int D>
  list<tuple<string,double>> ScalarFiniteElement<D> :: Timing () const
  {
    list<tuple<string,double>>timings;
    IntegrationRule ir(ElementType(), 2*Order());
    SIMD_IntegrationRule simdir(ElementType(), 2*Order());
    Vector<> shape(GetNDof()), coefs(GetNDof());
    Vector<> values(ir.Size());
    Matrix<> dvalues(ir.Size(), D);
    Vector<SIMD<double>> avalues(simdir.Size());
    Matrix<SIMD<double>> advalues(D, simdir.Size());
    FE_ElementTransformation<D,D> trafo(ElementType());
    static LocalHeap lh (10000000, "FE - Timing");
    HeapReset hr(lh);
    // auto & mir = trafo(ir, lh);
    auto & simdmir = trafo(simdir, lh);
    
    coefs = 1;
    
    double maxtime = 0.5;
    double time;

    constexpr size_t steps = 1000;
    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> CalcShape(ir[0], shape);
                     });
    timings.push_back(make_tuple("CalcShape", time/steps*1e9/GetNDof()));

    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> Evaluate(ir, coefs, values);
                     }, maxtime);
    timings.push_back(make_tuple("Evaluate",time/steps*1e9/(GetNDof()*ir.GetNIP())));

    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> Evaluate(simdir, coefs, avalues);
                     }, maxtime);
    timings.push_back(make_tuple("Evaluate(SIMD)", time/steps*1e9/(GetNDof()*ir.GetNIP())));

    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> EvaluateGrad(ir, coefs, dvalues);
                     }, maxtime);
    timings.push_back(make_tuple("Evaluate Grad", time/steps*1e9/(D*GetNDof()*ir.GetNIP())));

    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> EvaluateGrad(simdmir, coefs, advalues);
                     }, maxtime);
    timings.push_back(make_tuple("Evaluate Grad(SIMD)", time/steps*1e9/(D*GetNDof()*ir.GetNIP())));

    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> EvaluateTrans(ir, values, coefs);
                     }, maxtime);
    timings.push_back(make_tuple("Evaluate Trans", time/steps*1e9/(GetNDof()*ir.GetNIP())));

    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> AddTrans(simdir, avalues, coefs);
                     }, maxtime);
    timings.push_back(make_tuple("Evaluate Trans (SIMD)", time/steps*1e9/(GetNDof()*ir.GetNIP())));

    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> EvaluateGradTrans(ir, dvalues, coefs);
                     }, maxtime);
    timings.push_back(make_tuple("Evaluate Trans Grad", time/steps*1e9/(D*GetNDof()*ir.GetNIP())));

    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> AddGradTrans(simdmir, advalues, coefs);
                     }, maxtime);
    timings.push_back(make_tuple("Evaluate Trans Grad(SIMD)", time/steps*1e9/(D*GetNDof()*ir.GetNIP())));

    return timings;
  }


  
  template<int D>
  void ScalarFiniteElement<D> :: 
  EvaluateGradTrans (const IntegrationRule & ir, FlatMatrixFixWidth<D,double> vals, 
                     BareSliceVector<double> coefs) const
  {
    MatrixFixWidth<D> dshape(ndof);
    coefs.AddSize(ndof) = 0.0;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	CalcDShape (ir[i], dshape);
	coefs.AddSize(ndof) += dshape * vals.Row(i);
      }
  }


  template<int D>
  void ScalarFiniteElement<D> :: 
  EvaluateGradTrans (const IntegrationRule & ir, SliceMatrix<> values, SliceMatrix<> coefs) const
  {
#ifndef __CUDA_ARCH__
    cout << "EvalGradTrans not overloaded" << endl;
#endif
  }

  void BaseScalarFiniteElement :: 
  AddGradTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values,
                BareSliceVector<> coefs) const
  {
    throw ExceptionNOSIMD (string("AddGradTrans (simd) not implemented for class ")+typeid(*this).name());
  }


  








  template<int D>
  void ScalarFiniteElement<D> :: CalcDDShape (const IntegrationPoint & ip, 
                                              BareSliceMatrix<> ddshape) const
  {
    int nd = GetNDof();
    double eps = 1e-3;

    double pos[4] = { -2, -1, 1, 2 };
    double weight[4] = { 1.0/12, -2.0/3, 2.0/3, -1.0/12 };

    STACK_ARRAY(double, mem, nd*D);
    FlatMatrixFixWidth<D> dshape(nd, &mem[0]);  

    ddshape.AddSize(nd, D*D) = 0.0;
    
    for (int p = 0; p < 4; p++)
      for (int i = 0; i < D; i++)
        {
          IntegrationPoint ip1 = ip;
          ip1(i) += eps * pos[p];

          CalcDShape (ip1, dshape);

          for (int j = 0; j < nd; j++)
            for (int k = 0; k < D; k++)
              ddshape(j,D*i+k) += weight[p]/eps * dshape(j,k);
        }
  }


  template<int D>
  void ScalarFiniteElement<D> :: CalcMappedDDShape (const BaseMappedIntegrationPoint & bmip, 
                                                    BareSliceMatrix<> hddshape) const
  {
    auto & mip = static_cast<const MappedIntegrationPoint<D,D> &> (bmip);    
    int nd = GetNDof();
    auto ddshape = hddshape.AddSize(nd, D*D);
    double eps = 1e-7;
    MatrixFixWidth<D> dshape1(nd), dshape2(nd);
    const ElementTransformation & eltrans = mip.GetTransformation();

    for (int i = 0; i < D; i++)
      {
	IntegrationPoint ip1 = mip.IP();
	IntegrationPoint ip2 = mip.IP();
        ip1(i) -= eps;
        ip2(i) += eps;
        MappedIntegrationPoint<D,D> mip1(ip1, eltrans);
        MappedIntegrationPoint<D,D> mip2(ip2, eltrans);

	CalcMappedDShape (mip1, dshape1);
	CalcMappedDShape (mip2, dshape2);

        ddshape.Cols(D*i,D*(i+1)) = (0.5/eps) * (dshape2-dshape1);
      }  

    for (int j = 0; j < D; j++)
      {
        for (int k = 0; k < nd; k++)
          for (int l = 0; l < D; l++)
            dshape1(k,l) = ddshape(k, l*D+j);
        
        dshape2 = dshape1 * mip.GetJacobianInverse();
        
        for (int k = 0; k < nd; k++)
          for (int l = 0; l < D; l++)
            ddshape(k, l*D+j) = dshape2(k,l);
      }
    
  }



				  



  void BaseScalarFiniteElement :: GetDiagMassMatrix (FlatVector<> mass) const
  {
    throw Exception (string("mass matrix certainly not diagonal for element ")+typeid(*this).name());
  }






  template <int D>
  void DGFiniteElement<D>:: 
  GetDiagMassMatrix (FlatVector<> mass) const
  {
#ifndef __CUDA_ARCH__
    IntegrationRule ir(this->ElementType(), 2*order);
    VectorMem<50> shape(ndof);
    mass = 0;
    for (int i = 0; i < ir.Size(); i++)
      {
        this -> CalcShape (ir[i], shape);
        for (int j = 0; j < ndof; j++)
          mass(j) += ir[i].Weight() * sqr (shape(j));
      }
#endif
  }


  template <int D>
  void DGFiniteElement<D>:: 
  CalcTraceMatrix (int facet, FlatMatrix<> trace) const
  {
    ELEMENT_TYPE ftype = ElementTopology::GetFacetType (this->ElementType(), facet);
    Facet2ElementTrafo f2el(this->ElementType(), FlatArray<int> (8, const_cast<int*> (vnums)) );
    const IntegrationRule & ir = SelectIntegrationRule (ftype, 2*order);

    ScalarFiniteElement<0> * facetfe0 = NULL;
    ScalarFiniteElement<1> * facetfe1 = NULL;
    ScalarFiniteElement<2> * facetfe2 = NULL;
    switch (ftype)
      {
      case ET_POINT : facetfe0 = new FE_Point; break;
      case ET_SEGM : facetfe1 = new L2HighOrderFE<ET_SEGM> (order); break;
      case ET_TRIG : facetfe2 = new L2HighOrderFE<ET_TRIG> (order); break;
      case ET_QUAD : facetfe2 = new L2HighOrderFE<ET_QUAD> (order); break;
      default:
	;
      }

    int ndof_facet = trace.Height();
    Vector<> shape(ndof);
    Vector<> fshape(ndof_facet);
    Vector<> norms(ndof_facet);

    trace = 0.0;
    norms = 0.0;
    for (int i = 0; i < ir.Size(); i++)
      {
	if (D == 1) 
          facetfe0 -> CalcShape (ir[i], fshape);
	else if (D == 2) 
          facetfe1 -> CalcShape (ir[i], fshape);
	else            
          facetfe2 -> CalcShape (ir[i], fshape);

	this -> CalcShape (f2el (facet, ir[i]), shape);

	trace += ir[i].Weight() * fshape * Trans (shape);
	for (int j = 0; j < norms.Size(); j++)
	  norms(j) += ir[i].Weight() * sqr (fshape(j));
      }

    for (int j = 0; j < fshape.Size(); j++)
      trace.Row(j) /= norms(j);

    delete facetfe0;
    delete facetfe1;
    delete facetfe2;
  }

  template <int D>
  void DGFiniteElement<D>:: 
  CalcGradientMatrix (FlatMatrix<> gmat) const
  {
    IntegrationRule ir (this->ElementType(), 2*order);

    Vector<> shape(ndof);
    MatrixFixWidth<D> dshape(ndof);
    Vector<> norms(ndof);
    
    gmat = 0.0;
    norms = 0.0;
    for (int i = 0; i < ir.Size(); i++)
      {
	this -> CalcShape (ir[i], shape);
	this -> CalcDShape (ir[i], dshape);
        
        for (int j = 0; j < ndof; j++)
          for (int k = 0; k < ndof; k++)
            for (int l = 0; l < D; l++)
              gmat(k*D+l, j) += ir[i].Weight() * dshape(j,l) * shape(k);

	for (int j = 0; j < norms.Size(); j++)
	  norms(j) += ir[i].Weight() * sqr (shape(j));
      }
    for (int j = 0; j < ndof; j++)
      gmat.Rows(D*j, D*(j+1)) /= norms(j);
  }


  template <int D>
  void DGFiniteElement<D>:: 
  GetGradient (FlatVector<> coefs, FlatMatrixFixWidth<D> grad) const
  {
    Matrix<> gmat(D*grad.Height(), coefs.Size());
    CalcGradientMatrix (gmat);
    FlatVector<> vgrad(gmat.Height(), &grad(0,0));
    vgrad = gmat * coefs;
  }
  
  template <int D>
  void DGFiniteElement<D>:: 
  GetGradientTrans (FlatMatrixFixWidth<D> grad, FlatVector<> coefs) const 
  {
    Matrix<> gmat(D*grad.Height(), coefs.Size());
    CalcGradientMatrix (gmat);
    FlatVector<> vgrad(gmat.Height(), &grad(0,0));
    coefs = Trans (gmat) * vgrad;
  }
  
  template <int D>
  void DGFiniteElement<D>:: 
  GetTrace (int facet, FlatVector<> coefs, FlatVector<> fcoefs) const
  {
    Matrix<> trace(fcoefs.Size(), coefs.Size());
    CalcTraceMatrix(facet, trace);
    fcoefs = trace * coefs;
  }
  
  template <int D>
  void DGFiniteElement<D>:: 
  GetTraceTrans (int facet, FlatVector<> fcoefs, FlatVector<> coefs) const
  {
    Matrix<> trace(fcoefs.Size(), coefs.Size());
    CalcTraceMatrix(facet, trace);
    coefs = Trans (trace) * fcoefs;
  }



  template class ScalarFiniteElement<0>;
  template class ScalarFiniteElement<1>;
  template class ScalarFiniteElement<2>;
  template class ScalarFiniteElement<3>;


  template class DGFiniteElement<0>;
  template class DGFiniteElement<1>;
  template class DGFiniteElement<2>;
  template class DGFiniteElement<3>;
  

}

