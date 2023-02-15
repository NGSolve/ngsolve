/*********************************************************************/
/* File:   hcurlfe.cpp                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   16. Apr. 2000                                             */
/*********************************************************************/

/* 
   Nedelec's finite elements for
   Maxwell equations
*/

#define FILE_HCURLFE_CPP

#include <fem.hpp>
#include <thcurlfe_impl.hpp>
#include "hcurlhofe_impl.hpp"
#include "hcurllofe.hpp"
#include "hdivlofe.hpp"

namespace ngfem
{

  template <int D>
  string HCurlFiniteElement<D> :: ClassName(void) const
  { 
    return ToString ("HCurlFiniteElement<") + ToString(D) + ">";
  }


  template <int D>
  void HCurlFiniteElement<D> ::
  CalcCurlShape (const IntegrationPoint & ip, 
		 SliceMatrix<> curlshape) const
  {
    if (DIM == 1) return;
    double eps = 1e-6;  

    ArrayMem<double, 200> hm1(DIM*ndof), hm2(DIM*ndof), 
      hm3(DIM*ndof), hm4(DIM*ndof), hmi(DIM*ndof);

    FlatMatrixFixWidth<DIM> shape1(ndof, &hm1[0]);
    FlatMatrixFixWidth<DIM> shape2(ndof, &hm2[0]);
    FlatMatrixFixWidth<DIM> shape3(ndof, &hm3[0]);
    FlatMatrixFixWidth<DIM> shape4(ndof, &hm4[0]);

    FlatMatrixFixWidth<DIM> dshapei(ndof, &hmi[0]);
    curlshape = 0;

    for (int i = 0; i < DIM; i++)
      {
	IntegrationPoint ip1 = ip;
	IntegrationPoint ip2 = ip;
	ip1(i) -= eps;
	ip2(i) += eps;

	CalcShape (ip1, shape1);
	CalcShape (ip2, shape2);

	ip1(i) -= eps;
	ip2(i) += eps;
	CalcShape (ip1, shape3);
	CalcShape (ip2, shape4);
	
	dshapei = 
	  2/(3*eps) * (shape2 - shape1) 
	  - 1/(12*eps) * (shape4 - shape3);


	if (DIM == 3)
	  switch (i)
	    {
	    case 0:
	      {
		for (int j = 0; j < ndof; j++)
		  {
		    curlshape(j,1) -= dshapei(j,2);
		    curlshape(j,2) += dshapei(j,1);
		  }
		break;
	      }
	    case 1:
	      {
		for (int j = 0; j < ndof; j++)
		  {
		    curlshape(j,2) -= dshapei(j,0);
		    curlshape(j,0) += dshapei(j,2);
		  }
		break;
	      }
	    case 2:
	      {
		for (int j = 0; j < ndof; j++)
		  {
		    curlshape(j,0) -= dshapei(j,1);
		    curlshape(j,1) += dshapei(j,0);
		  }
		break;
	      }
	    }
	else
	  {
	    switch (i)
	      {
	      case 0:
		{
		  for (int j = 0; j < ndof; j++)
		    curlshape(j,0) += dshapei(j,1);
		  break;
		}
	      case 1:
		{
		  for (int j = 0; j < ndof; j++)
		    curlshape(j,0) -= dshapei(j,0);
		  break;
		}
	      }
	    
	  }
	
      }
  }


  /// compute shape
  
  template <int D>
  void HCurlFiniteElement<D> ::
  CalcMappedShape (const BaseMappedIntegrationPoint & bmip,
                   SliceMatrix<> shape) const
  {
    CalcShape (bmip.IP(), shape);

    Switch<4-DIM>
      (bmip.DimSpace()-DIM,
       [this,&bmip,shape](auto CODIM)
       {
         auto & mip = static_cast<const MappedIntegrationPoint<DIM,DIM+CODIM.value>&> (bmip);
         auto trans = Trans (mip.GetJacobianInverse());
         
         for (int i = 0; i < ndof; i++)
           shape.Row(i).Range(DIM+CODIM.value) = trans * Vec<DIM> (shape.Row(i));
       });
  }

  template <int D>
  void HCurlFiniteElement<D> ::
  CalcMappedShape (const BaseMappedIntegrationRule & mir, SliceMatrix<> shapes) const
  {
    for (int i = 0; i < mir.Size(); i++)
      CalcMappedShape (mir[i], shapes.Cols(i*D, (i+1)*D));
  }

  template <int D>
  void HCurlFiniteElement<D> ::
  CalcMappedShape (const SIMD<BaseMappedIntegrationPoint> & bmip, BareSliceMatrix<SIMD<double>> shapes) const
  {
    throw ExceptionNOSIMD(string("SIMD - HCurlFE::CalcMappedShape not overloaded, et = ")
                          + typeid(*this).name());
  }

  
  template <int D>
  void HCurlFiniteElement<D> ::
  CalcMappedShape (const SIMD_BaseMappedIntegrationRule & mir, 
                   BareSliceMatrix<SIMD<double>> dshapes) const
  {
    throw ExceptionNOSIMD(string("SIMD - HCurlFE::CalcShape not overloaded, et = ")
                          + typeid(*this).name());
  }


  template <int D>
  void HCurlFiniteElement<D> ::
  Evaluate (const IntegrationRule & ir, BareSliceVector<> coefs, SliceMatrix<> values) const
  {
    LocalHeapMem<100000> lhdummy("hcurlfe-lh");
    for (int i = 0; i < ir.Size(); i++)
      values.Row(i) = EvaluateShape (ir[i], coefs, lhdummy);
  }
  
  template <int D>
  void HCurlFiniteElement<D> ::
  Evaluate (const MappedIntegrationRule<D,D> & mir, BareSliceVector<> coefs, SliceMatrix<> values) const
  {
    LocalHeapMem<100000> lhdummy("hcurlfe-lh");
    for (int i = 0; i < mir.Size(); i++)
      values.Row(i) = Trans(mir[i].GetJacobianInverse()) * EvaluateShape (mir[i].IP(), coefs, lhdummy);
  }
  

  /// compute curl of shape
  template <int D>
  void HCurlFiniteElement<D> ::
  CalcMappedCurlShape (const BaseMappedIntegrationPoint & bmip,
                       SliceMatrix<> curlshape) const
  {
    auto & mip = static_cast<const MappedIntegrationPoint<D,D>&> (bmip);
    CalcCurlShape (mip.IP(), curlshape);
    if (DIM == 2)
      {
        curlshape /= mip.GetJacobiDet();        
      }
    else
      {
        Mat<DIM> trans = (1.0/mip.GetJacobiDet()) * mip.GetJacobian();
        for (int i = 0; i < ndof; i++)
          {
            Vec<DIM> hs = curlshape.Row(i);
            curlshape.Row(i) = trans * hs;
          }
      }
  }

  template <int D>
  void HCurlFiniteElement<D> ::
  CalcMappedCurlShape (const MappedIntegrationRule<DIM,DIM> & mir, 
                       SliceMatrix<> curlshape) const
  {
    for (int i = 0; i < mir.Size(); i++)
      CalcMappedCurlShape (mir[i], curlshape.Cols(i*DIM_CURL_(D), (i+1)*DIM_CURL_(D)));
  }
  
  template <int D>
  void HCurlFiniteElement<D> ::
  CalcMappedCurlShape (const SIMD_BaseMappedIntegrationRule & mir, 
                       BareSliceMatrix<SIMD<double>> dshapes) const
  {
    throw ExceptionNOSIMD("SIMD - HCurlFE::CalcMappedCurlShape not overloaded");
  }
  
  

  template <int D>
  void HCurlFiniteElement<D> ::
  EvaluateCurl (const IntegrationRule & ir, BareSliceVector<> coefs, FlatMatrixFixWidth<DIM_CURL_(D)> curl) const
  {
    LocalHeapMem<10000> lhdummy("hcurlfe-lh");
    for (int i = 0; i < ir.Size(); i++)
      curl.Row(i) = EvaluateCurlShape (ir[i], coefs, lhdummy);
  }

  template <int D>
  void HCurlFiniteElement<D> ::
  EvaluateMappedCurl (const MappedIntegrationRule<D,D> & mir, 
                      BareSliceVector<> coefs, FlatMatrixFixWidth<DIM_CURL_(D)> curl) const
  {
    /*
    LocalHeapMem<1000> lhdummy("dummy");
    if (D == 2)
      for (int i = 0; i < mir.Size(); i++)
        curl.Row(i) = (1.0/mir[i].GetJacobiDet())*EvaluateCurlShape (mir[i].IP(), coefs, lhdummy);
    else
      {
        for (int i = 0; i < mir.Size(); i++)
          {
            Mat<DIM_CURL_TRAIT<D>::DIM> trans = (1.0/mir[i].GetJacobiDet()) * mir[i].GetJacobian();
            curl.Row(i) = trans*EvaluateCurlShape (mir[i].IP(), coefs, lhdummy);
          }
      }
    */
    EvaluateCurl (mir.IR(), coefs, curl);
    if (D == 2)
      for (int i = 0; i < mir.Size(); i++)
        curl.Row(i) *= 1.0/mir[i].GetJacobiDet();
    else
      {
        for (int i = 0; i < mir.Size(); i++)
          {
            Vec<DIM_CURL_(D)> hv = curl.Row(i);
            Mat<DIM_CURL_(D)> trans = (1.0/mir[i].GetJacobiDet()) * mir[i].GetJacobian();
            curl.Row(i) = trans*hv;
          }
      }
  }



  template<int D>
  list<tuple<string,double>> HCurlFiniteElement<D> :: Timing () const
  {
    list<tuple<string,double>>timings;
    IntegrationRule ir(ElementType(), 2*Order());
    SIMD_IntegrationRule simdir(ElementType(), 2*Order());

    constexpr int DIMC = D*(D-1)/2;
    Matrix<> shape(GetNDof(),D);
    Vector<> coefs(GetNDof());
    Matrix<> values(ir.Size(), D);
    Matrix<> dvalues(ir.Size(), DIMC);
    Matrix<SIMD<double>> avalues(D,simdir.Size());
    Matrix<SIMD<double>> advalues(DIMC, simdir.Size());
    Matrix<SIMD<double>> simd_shapes(DIMC*GetNDof(), simdir.Size());
    FE_ElementTransformation<D,D> trafo(ElementType());

    LocalHeap lh (10000000, "FE - Timing");
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
    timings.push_back(make_tuple("CalcShape", time/D/steps*1e9/GetNDof()));

    try
      {
        time = RunTiming([&]() {
            for (size_t i = 0; i < steps; i++)
              this -> CalcMappedShape(simdmir, simd_shapes);
          });
        timings.push_back(make_tuple("CalcShape (SIMD)", time/D/steps*1e9/(simdir.GetNIP()*GetNDof())));
      }
    catch (const ExceptionNOSIMD& e) { };

    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> Evaluate(ir, coefs, values);
                     }, maxtime);
    timings.push_back(make_tuple("Evaluate",time/D/steps*1e9/(GetNDof()*ir.GetNIP())));
    
    try
      {
        time = RunTiming([&]() {
            for (size_t i = 0; i < steps; i++)
              this -> Evaluate(simdmir, coefs, avalues);
          }, maxtime);
        timings.push_back(make_tuple("Evaluate(SIMD)", time/D/steps*1e9/(GetNDof()*ir.GetNIP())));
      }
    catch (const ExceptionNOSIMD& e) { };

    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> EvaluateCurl(ir, coefs, dvalues);
                     }, maxtime);
    timings.push_back(make_tuple("Evaluate Curl", time/DIMC/steps*1e9/(D*GetNDof()*ir.GetNIP())));

    try
      {
        time = RunTiming([&]() {
            for (size_t i = 0; i < steps; i++)
              this -> EvaluateCurl(simdmir, coefs, advalues);
          }, maxtime);
        timings.push_back(make_tuple("Evaluate Curl(SIMD)", time/DIMC/steps*1e9/(D*GetNDof()*ir.GetNIP())));
      }
    catch (const ExceptionNOSIMD& e) { };

    
    /*
    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> EvaluateTrans(ir, values, coefs);
                     }, maxtime);
    timings.push_back(make_tuple("Evaluate Trans", time/steps*1e9/(GetNDof()*ir.GetNIP())));
    */

    try
      {
        time = RunTiming([&]() {
            for (size_t i = 0; i < steps; i++)
              this -> AddTrans(simdmir, avalues, coefs);
          }, maxtime);
        timings.push_back(make_tuple("Evaluate Trans (SIMD)", time/D/steps*1e9/(GetNDof()*ir.GetNIP())));
      }
    catch (const ExceptionNOSIMD& e) { };
    
    /*
    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> EvaluateCurlTrans(ir, dvalues, coefs);
                     }, maxtime);
    timings.push_back(make_tuple("Evaluate Trans Curl", time/steps*1e9/(D*GetNDof()*ir.GetNIP())));
    */

    try
      {
        time = RunTiming([&]() {
            for (size_t i = 0; i < steps; i++)
              this -> AddCurlTrans(simdmir, advalues, coefs);
          }, maxtime);
        timings.push_back(make_tuple("Evaluate Trans Curl(SIMD)", time/DIMC/steps*1e9/(D*GetNDof()*ir.GetNIP())));
      }
    catch (const ExceptionNOSIMD& e) { };

    return timings;
  }

  

  template <int D>
  void HCurlFiniteElement<D> ::
  ComputeEdgeMoments (int enr, ScalarFiniteElement<1> & testfe,
		      FlatMatrix<> moments, int order, int shapenr) const
  {
    int test_ndof = testfe.GetNDof();
    
    MatrixFixWidth<DIM> shape(ndof);
    Vector<> shapetau(ndof);
    Vector<> testshape(test_ndof);
    Vector<> tau(D), p1(D), p2(D), p(D);
    
    const IntegrationRule & linerule = 
      SelectIntegrationRule (ET_SEGM, order);
    
    const POINT3D * points = ElementTopology::GetVertices (ElementType());
    const EDGE & edge = ElementTopology::GetEdges (ElementType()) [enr];
    
    
    for (int j = 0; j < D; j++)
      {
	p1(j) = points[edge[0]][j];
	p2(j) = points[edge[1]][j];
      }
    
    tau = p2 - p1;
    moments = 0;
    
    for (int j = 0; j < linerule.GetNIP(); j++)
      {
	const IntegrationPoint & ip = linerule[j];
	
	p = p1 + ip.Point()[0] * tau;
	IntegrationPoint ip3d(p, 0);
      
	testfe.CalcShape (ip, testshape);

	if (shapenr == 1)
	  CalcShape1 (ip3d, shape);
	else
	  CalcShape2 (ip3d, shape);

	shapetau = shape * tau;

	// moments += ChangeSize(ip.Weight() * (testshape * Trans (shapetau)),moments.Height(),moments.Width());
        moments.Rows(0, testshape.Height()) += ip.Weight() * (testshape * Trans (shapetau));
      }
  }


  template <int D>
  void HCurlFiniteElement<D> ::
  ComputeFaceMoments (int fnr, HDivFiniteElement<2> & testfe,
		      FlatMatrix<> moments, int order, int shapenr) const
  {
    int j;
    int test_ndof = testfe.GetNDof();

    MatrixFixWidth<DIM> shape(ndof);
    Matrix<> shapetau(ndof, 2);
    MatrixFixWidth<2> testshape(test_ndof);
    Matrix<> tau(D, 2);
    
    const IntegrationRule & facerule = 
      SelectIntegrationRule (testfe.ElementType(), order);
    
    const POINT3D * points = ElementTopology::GetVertices (ElementType());
    const FACE & face = ElementTopology::GetFaces (ElementType()) [fnr];
    
    Vector<> p1(D), p2(D), p3(D), p(D);

    for (j = 0; j < D; j++)
      {
	if (testfe.ElementType() == ET_TRIG)
	  {
	    p1(j) = points[face[0]][j];
	    p2(j) = points[face[1]][j];
	    p3(j) = points[face[2]][j];
	  }
	else
	  {
	    p1(j) = points[face[1]][j];
	    p2(j) = points[face[3]][j];
	    p3(j) = points[face[0]][j];
	  }
	tau(j,0) = p1(j) - p3(j);
	tau(j,1) = p2(j) - p3(j);
      }

    moments = 0;
    
    for (j = 0; j < facerule.GetNIP(); j++)
      {
	const IntegrationPoint & ip = facerule[j];
	
	Vec<2> p2d;

	p2d(0) = ip(0); 
	p2d(1) = ip(1);
	p = p3 + tau * p2d;
	
	IntegrationPoint ip3d(p, 0);

	testfe.CalcShape (ip, testshape);
	
	switch (shapenr)
	  {
	  case 1:
	    {
	      CalcShape1 (ip3d, shape);
	      break;
	    }
	  case 2:
	    {
	      CalcShape2 (ip3d, shape);
	      break;
	    }
	  case 3:
	    {
	      CalcShape3 (ip3d, shape);
	      break;
	    }
	  case 4:
	    {
	      CalcShape4 (ip3d, shape);
	      break;
	    }
	  default:
	    throw Exception ("illegal face shape functions class");
	  }
	
	shapetau = shape * tau;
	// moments += ChangeSize(ip.Weight() * testshape * Trans (shapetau),moments.Height(),moments.Width());
        moments.Rows(0, testshape.Height()) += ip.Weight() * (testshape * Trans (shapetau));

      }
  }


  template <int D>
  void HCurlFiniteElement<D> ::
  ComputeVolMoments (HDivFiniteElement<3> & testfe,
		     FlatMatrix<> moments, int order, int shapenr) const
  {
    int j;
    int test_ndof = testfe.GetNDof();

    MatrixFixWidth<DIM> shape(ndof);
    MatrixFixWidth<3> testshape(test_ndof);
    
    const IntegrationRule & rule = 
      SelectIntegrationRule (ElementType(), order);


    moments = 0;
    for (j = 0; j < rule.GetNIP(); j++)
      {
	const IntegrationPoint & ip = rule[j];

	testfe.CalcShape (ip, testshape);

	switch (shapenr)
	  {
	  case 1:
	    CalcShape1 (ip, shape);
	    break;
	  case 2:
	    CalcShape2 (ip, shape);
	    break;
	  case 3:
	    CalcShape3 (ip, shape);
	    break;
	  case 4:
	    CalcShape4 (ip, shape);
	    break;
	  }
	// moments += ChangeSize(ip.Weight() * testshape * Trans (shape),moments.Height(),moments.Width());
        moments.Rows(0, testshape.Height()) += ip.Weight() * (testshape * Trans (shape));
      }
  }
  
  template <int D>
  void HCurlFiniteElement<D> ::
  CalcDualShape (const BaseMappedIntegrationPoint & bmip, SliceMatrix<> shape) const
  {
    // throw Exception(string("CalcDualShape not implemented for H(curl) element ")+typeid(*this).name());
    static bool first = true;
    if (first)
      cerr << "CalcDualShape not implemented for H(curl) element " << typeid(*this).name() << endl;
    first = false;
  }

  template <int D>
  void HCurlFiniteElement<D> ::
  CalcDualShape (const SIMD_BaseMappedIntegrationRule & bmir, BareSliceMatrix<SIMD<double>> shape) const
  {
    //throw ExceptionNOSIMD (string("CalcDualShape SIMD not implemented for H(curl) element ") +typeid(*this).name());
    static bool firstsimd = true;
    if (firstsimd)
      cerr << "CalcDualShape SIMD not implemented for H(curl) element " << typeid(*this).name() << endl;
    firstsimd = false;
  }

  
  template class HCurlFiniteElement<1>;
  template class HCurlFiniteElement<2>;
  template class HCurlFiniteElement<3>;

  /*
    Input: Matrix a, H <= W
    Compute: Regular WxW matrix inv such that A Inv = (I 0)

    decomposes 
    A = L Q

    L is an lower triangular matrix of size H * H
    Q is a orthogonal matrix as product of Householder factors:
    Q = Q_0 ... Q_{W-2}
    Q_i = I - u_i u_i^T / c_i       u_ij = 0 for j < i
  */
  inline void Householder (FlatMatrix<> & a, FlatMatrix<> inv)
  {

    (*testout) << "a = " << a << endl;

    int h = a.Height();
    int w = a.Width();


    Vector<> c(h);
    Vector<> d(h);
    int i, j, k;
    double scale, sum, sigma, tau;
    for (k = 0; k < h; k++)
      {
	scale = 0;
	for (i = k; i < w; i++)
	  if (fabs(a(k,i)) > scale) scale = fabs(a(k,i));

	if (scale == 0)
	  {
	    (*testout) << "breakdown in householder, a = " << a << endl;
	    throw Exception ("Housholder: Matrix not full rank");
	  }

	for (i = k; i < w; i++)
	  a(k,i) /= scale; 

	sum = 0;
	for (i = k; i < w; i++)
	  sum += a(k,i) * a(k,i);

	sum = sqrt (sum);

	sigma = 0;
	if (a(k,k) > 0)
	  sigma = sum;
	else // if (a(k,k) < 0)
	  sigma = -sum;

	a(k,k) += sigma;
	c(k) = sigma*a(k,k);
	d(k) = -scale * sigma;
	for (j = k+1; j < h; j++)
	  {
	    sum = 0;
	    for (i = k; i < w; i++)
	      sum += a(k,i) * a(j,i);
	    tau = sum / c(k);
	    for (i = k; i < w; i++)
	      a(j,i) -= tau * a(k,i);
	  }
      }
    //    d(h-1) = a(h-1, h-1);

    inv = 0;
    for (k = 0; k < w; k++)
      inv(k,k) = 1;
    for (k = 0; k < w; k++)
      {
	for (i = 0; i < h; i++)
	  {
	    for (j = 0; j < i; j++)
	      inv(i,k) -= a(i,j) * inv(j,k);
	    inv(i,k) /= d(i);
	  }
      }

    for (k = 0; k < w; k++)
      {
	for (i = h-1; i >= 0; i--)
	  {
	    sum = 0;
	    for (j = i; j < w; j++)
	      sum += a(i,j) * inv(j,k);
	    sum /= c(i);
	    for (j = i; j < w; j++)
	      inv(j,k) -= sum * a(i,j);
	  }
      }
  }



  /* ************************* Edge potential elements *************** */
  
  class FE_Trig3EdgeBubble : public ScalarFiniteElement<2>
  {
  public:
    ///
    FE_Trig3EdgeBubble()
      : ScalarFiniteElement<2> (6, 3) { ; } 

    virtual ELEMENT_TYPE ElementType() const override { return ET_TRIG; }

    ///
    using ScalarFiniteElement<2>::CalcShape;
    virtual void CalcShape (const IntegrationPoint & ip, 
			    BareSliceVector<> shape) const override
    {
      double x = ip(0);
      double y = ip(1);
      double l3 = 1-x-y;

      shape.Range(0,ndof) = 0.0; //!

      shape(0) = x * l3;
      shape(1) = x * (x-l3) * l3;
      shape(2) = y * l3;
      shape(3) = y * (y-l3) * l3;
      shape(4) = x * y;
      shape(5) = x * (x-y) * y;
    }

    virtual void CalcDShape (const IntegrationPoint & ip, 
			     BareSliceMatrix<> dshape) const override
    {
      double x = ip(0);
      double y = ip(1);
      double l3 = 1-x-y;

      // shape(0) = x * l3;
      dshape(0,0) = l3-x;
      dshape(0,1) = -x;

      // shape(1) = x * (x-l3) * l3;
      dshape(1,0) = 4*x*l3-x*x-l3*l3;
      dshape(1,1) = -x*x+2*x*l3;

      // shape(2) = y * l3;
      dshape(2,0) = -y;
      dshape(2,1) = l3-y;

      // shape(3) = y * (y-l3) * l3;
      dshape(3,0) = -y*y+2*y*l3;
      dshape(3,1) = 4*y*l3-y*y-l3*l3;

      // shape(4) = x * y;
      dshape(4,0) = y;
      dshape(4,1) = x;

      // shape(5) = x * (x-y) * y;
      dshape(5,0) = 2*x*y-y*y;
      dshape(5,1) = x*x-2*x*y;
    }

    ///
    // virtual const Array<IPData> & GetIPData () const { return ipdata; }
  }; 






  /* ***************** Gradient matrix ********************* */


  template <int D>
  void ComputeGradientMatrix (const ScalarFiniteElement<D> & h1fe,
			      const HCurlFiniteElement<D> & hcurlfe,
			      FlatMatrix<> gradient)
  {
    int ndh1 = h1fe.GetNDof();
    int ndhcurl = hcurlfe.GetNDof();
    
    Matrix<> mathchc(ndhcurl);
    Matrix<> invmathchc(ndhcurl);
    Matrix<> mathch1(ndhcurl, ndh1);
    Matrix<> dshapeh1(ndh1, D);
    MatrixFixWidth<D> shapehcurl(ndhcurl);

    const IntegrationRule & ir = 
      SelectIntegrationRule (h1fe.ElementType(), 2*hcurlfe.Order());
    
    mathchc = 0;
    mathch1 = 0;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	const IntegrationPoint & ip = ir[i];
	h1fe.CalcDShape (ip, dshapeh1);
	hcurlfe.CalcShape (ip, shapehcurl);
	
	mathchc += ip.Weight() * (shapehcurl * Trans (shapehcurl));
	mathch1 += ip.Weight() * (shapehcurl * Trans (dshapeh1));
      }

    CalcInverse (mathchc, invmathchc);
    gradient = invmathchc * mathch1;

    (*testout) << " Compute Gradient Matrix H1-HCurl Low order FEs " << endl
	       << gradient << endl; 
  }
  
  template
  void ComputeGradientMatrix<2> (const ScalarFiniteElement<2> & h1fe,
				 const HCurlFiniteElement<2> & hcurlfe,
				 FlatMatrix<> gradient);

  template
  void ComputeGradientMatrix<3> (const ScalarFiniteElement<3> & h1fe,
				 const HCurlFiniteElement<3> & hcurlfe,
				 FlatMatrix<> gradient);




  /* ******************  Workaround for the destruction of static ipdata-members ************* */
  class ControlDestruction
  {
  private:
    Array<bool*> isdestructed;
  public:
    void NewControl(bool * variable)
    {
      isdestructed.Append(variable);
      *variable = false;
    }
    ~ControlDestruction()
    {
      for(int i=0; i<isdestructed.Size(); i++)
	*isdestructed[i] = true;
    }
  };

  ControlDestruction controldestruction;





  /* ******************** Segm elements ********************* */

 
  // Array<HCurlFiniteElement<1>::IPData> FE_NedelecSegm1::ipdata;

  FE_NedelecSegm1 :: FE_NedelecSegm1()
    : HCurlFiniteElement<1>(NDOF, 0) { ; }
  
  void FE_NedelecSegm1 :: 
  CalcShape (const IntegrationPoint & ip, 
	     SliceMatrix<> shape) const
  {
    shape = 0.0; //!
    shape (0,0) = 1;
  }




  // Array<HCurlFiniteElement<1>::IPData> FE_NedelecSegm2::ipdata;

  FE_NedelecSegm2 :: FE_NedelecSegm2()
    : HCurlFiniteElement<1>(NDOF, 1)  { ; }
  
  void FE_NedelecSegm2 :: 
  CalcShape (const IntegrationPoint & ip, 
	     SliceMatrix<> shape) const
  {
    shape = 0.0; //!
    shape (0,0) = 1;
    shape (1,0) = 2*ip(0)-1;
  }






  // Array<HCurlFiniteElement<1>::IPData> FE_NedelecSegm3::ipdata;

  FE_NedelecSegm3 :: FE_NedelecSegm3()
    : HCurlFiniteElement<1>(3, 2)
  { ; }


  void FE_NedelecSegm3 :: 
  CalcShape (const IntegrationPoint & ip, 
	     SliceMatrix<> shape) const
  {
    shape = 0.0; //!
    shape (0,0) = 1;
    shape (1,0) = 2*ip(0)-1;
    shape (2,0) = ip(0) * (1-ip(0));
  }










  /* ******************** triangular elements *********************** */

  /*
    Array<HCurlFiniteElement<2>::IPData> FE_NedelecTrig1::ipdata;

    FE_NedelecTrig1 :: FE_NedelecTrig1()
    : HCurlFiniteElement<2> (ET_TRIG, 3, 1)
    {
    CalcIPData(ipdata);
    controldestruction.NewControl(&ipdatadestructed);
    }

    FE_NedelecTrig1 :: ~FE_NedelecTrig1()
    {
    if(!ipdatadestructed)
    ipdata.DeleteAll();
    }

    void FE_NedelecTrig1 :: 
    CalcShape (const IntegrationPoint & ip, 
    FlatMatrixFixWidth<2> shape) const
    {
    double x = ip(0);
    double y = ip(1);

    shape = 0.0; //!
    shape (0,0) = 1 - y;
    shape (1,0) = -y;
    shape (2,0) = -y;
    
    shape (0,1) = x;
    shape (1,1) = -(1 - x);
    shape (2,1) = x;
    }




  
    Array<HCurlFiniteElement<2>::IPData> FE_NedelecTrig2::ipdata;
    Mat<FE_NedelecTrig2::NDOF> FE_NedelecTrig2::trans;

    FE_NedelecTrig2 :: FE_NedelecTrig2()
    : HCurlFiniteElement<2>(ET_TRIG, NDOF, 1)
    {
    if (!ipdata.Size())
    Orthogonalize();

    CalcIPData(ipdata);
    controldestruction.NewControl(&ipdatadestructed);
    }

    FE_NedelecTrig2 :: ~FE_NedelecTrig2()
    {
    if(!ipdatadestructed)
    ipdata.DeleteAll();
    }


    void FE_NedelecTrig2 :: CalcShape (const IntegrationPoint & ip, 
    FlatMatrixFixWidth<2> shape) const
    {
    Mat<NDOF,2> shape1;
    CalcShape1 (ip, shape1);

    shape = 0.0; //!
    shape = ChangeSize(Trans (trans) * shape1,shape.Height(),2);
    }

    void FE_NedelecTrig2 :: 
    CalcShape1 (const IntegrationPoint & ip, 
    FlatMatrixFixWidth<2> shape) const
    {
    double x = ip(0);
    double y = ip(1);
    shape = 0.0;

    shape (0,0) = 1;
    shape (1,0) = x;
    shape (2,0) = y;
    shape (3,1) = 1;
    shape (4,1) = x;
    shape (5,1) = y;
    }


    void FE_NedelecTrig2 :: Orthogonalize()
    {
    Mat<NDOF> fiphij;

    Matrix<> edgemoments(2, NDOF);
    FE_Segm1L2 segm1;
    
    for (int i = 0; i < 3; i++)
    {
    ComputeEdgeMoments (i, segm1, edgemoments, 2);

    for (int j = 0; j < NDOF; j++)
    {
    fiphij(i, j)    = edgemoments(0, j);
    fiphij(3+i, j)  = edgemoments(1, j);
    }
    }

    CalcInverse (fiphij, trans);
    }
  */


  
  

  /*
    Array<HCurlFiniteElement<2>::IPData> FE_NedelecTrig3::ipdata;
    Mat<FE_NedelecTrig3::NDOF> FE_NedelecTrig3::trans;
    Mat<FE_NedelecTrig3::NEDGEDOF> FE_NedelecTrig3::trans2;
  
    FE_NedelecTrig3 :: FE_NedelecTrig3()
    : HCurlFiniteElement<2>(ET_TRIG, NDOF, 2)
    {
    if (!ipdata.Size())
    Orthogonalize();

    CalcIPData(ipdata);
    controldestruction.NewControl(&ipdatadestructed);
    }

    FE_NedelecTrig3 :: ~FE_NedelecTrig3()
    {
    if(!ipdatadestructed)
    ipdata.DeleteAll();
    }

    void FE_NedelecTrig3 :: CalcShape (const IntegrationPoint & ip, 
    FlatMatrixFixWidth<2> shape) const
    {
    shape = 0.0; //!
    Mat<NDOF,2> hshape;
    CalcShape1 (ip, hshape);
    shape = Trans (trans) * hshape;

    Mat<NEDGEDOF,2> shape2, hshape2;
    CalcShape2 (ip, hshape2);
    shape2 = Trans (trans2) * hshape2;
    
    int i, j;
    for (i = 0; i < NEDGEDOF; i++)
    for (j = 0; j < 2; j++)
    shape(i+3,j) = shape2(i,j);

    Mat<3,2> shape1;
    trig1.CalcShape(ip, shape1);
    for (i = 0; i < 3; i++)
    for (j = 0; j < 2; j++)
    shape(i,j) = shape1(i,j);
    }
  
    void FE_NedelecTrig3 :: CalcShape1 (const IntegrationPoint & ip, 
    FlatMatrixFixWidth<2> shape) const
    {
    double x = ip(0);
    double y = ip(1);
    shape = 0.0;
    
    for (int i = 0; i < 2; i++)
    {
    int base = 6 * i;
    shape (base  ,i) = 1;
    shape (base+1,i) = x;
    shape (base+2,i) = y;
    shape (base+3,i) = x*x;
    shape (base+4,i) = x*y;
    shape (base+5,i) = y*y;
    }
    }

    // base functions which are gradients of edge shape functions
    void FE_NedelecTrig3 :: CalcShape2 (const IntegrationPoint & ip, 
    FlatMatrixFixWidth<2> shape) const
    {
    FE_Trig3EdgeBubble febub;
    febub.CalcDShape (ip, FlatMatrix<> (shape));
    }


    void FE_NedelecTrig3 :: Orthogonalize()
    {
    // int i, j, k, l;

    // Matrix<> fiphij(NDOF);
    Mat<NDOF> fiphij;

    Matrix<> edgemoments(3, NDOF);
    FE_Segm2L2 segm2;
    
    for (int i = 0; i < 3; i++)
    {
    ComputeEdgeMoments (i, segm2, edgemoments, 4);
    for (int j = 0; j < NDOF; j++)
    {
    fiphij(i  , j) = edgemoments(0, j);
    fiphij(3+i, j) = edgemoments(1, j);
    fiphij(6+i, j) = edgemoments(2, j);
    }
    }

    Matrix<> facemoments(3, NDOF);
    FE_RTTrig0 rttrig0;
    ComputeFaceMoments (0, rttrig0, facemoments, 4);

    for (int j = 0; j < NDOF; j++)
    {
    fiphij(9 , j) = facemoments(1, j);
    fiphij(10, j) = facemoments(0, j);
    fiphij(11, j) = facemoments(2, j);
    }

    CalcInverse (fiphij, trans);
    
    // edge shape functions:

    int nd = NEDGEDOF;
    Mat<NEDGEDOF> fiphij2;
    Matrix<> edgemoments2(3, nd);

    for (int i = 0; i < 3; i++)
    {
    ComputeEdgeMoments (i, segm2, edgemoments, 4, 2);
    for (int j = 0; j < nd; j++)
    {
    fiphij2(i, j) = edgemoments(1, j);
    fiphij2(3+i, j) = edgemoments(2, j);
    }
    }
    CalcInverse (fiphij2, trans2);
    }
  
  */







  /* ******************** quad elements *********************** */

  /*
  // Array<HCurlFiniteElement<2>::IPData> FE_NedelecQuad1::ipdata;

  FE_NedelecQuad1 :: FE_NedelecQuad1()
    : HCurlFiniteElement<2> (ET_QUAD, NDOF, 1)

  {
    // CalcIPData(ipdata);
    // controldestruction.NewControl(&ipdatadestructed);
  }

  FE_NedelecQuad1 :: ~FE_NedelecQuad1()
  {
    // if(!ipdatadestructed)
    // ipdata.DeleteAll();
  }

  void FE_NedelecQuad1 :: 
  CalcShape (const IntegrationPoint & ip, 
	     FlatMatrixFixWidth<2> shape) const
  {
    shape = 0;
    double x = ip(0);
    double y = ip(1);

    shape(0,0) = 1-y;
    shape(1,0) = -y;
    shape(2,1) = -1+x;
    shape(3,1) = x;
  }
  */
  

  
  // Q(ORDER-1, ZORDER-2) x Q(ORDER-2, ZORDER-1)
  template <int ORDER, int ZORDER>
  class FE_TFaceTest : public HDivFiniteElement<2>
  {
  private:
    ///
    // Array<IPData> ipdata;

  public:
    enum { NDOF = (ORDER-1) * ZORDER + ORDER * (ZORDER-1) };
    enum { MAXORDER = (ORDER > ZORDER) ? ORDER : ZORDER };

    ///
    FE_TFaceTest()
      : HDivFiniteElement<2> (NDOF, MAXORDER)
    { ; }

    ///
    virtual ~FE_TFaceTest()
    { 
      // ipdata.DeleteAll(); 
    }

    virtual ELEMENT_TYPE ElementType() const { return ET_QUAD; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    SliceMatrix<> shape) const
    {
      double x = ip(0);
      double y = ip(1);

      Vec<ORDER> polx;
      Vec<ZORDER> poly;

      FE_TSegmL2<ORDER-1> fex;
      FE_TSegmL2<ZORDER-1> fey;
      
      IntegrationPoint ipx(x, 0, 0, 1);
      IntegrationPoint ipy(y, 0, 0, 1);

      fex.CalcShape (ipx, polx);
      fey.CalcShape (ipy, poly);

      shape = 0;
      int i, j, ii = 0;
      for (i = 0; i < ORDER; i++)
	for (j = 0; j < ZORDER-1; j++)
	  shape(ii++, 0) = polx(i) * poly(j);
      for (i = 0; i < ORDER-1; i++)
	for (j = 0; j < ZORDER; j++)
	  shape(ii++, 1) = polx(i) * poly(j);
    }

    ///
    // virtual const Array<IPData> & GetIPData () const 
    // { return ipdata; }
  };





  
  // template <int ORDER, int ZORDER>
  // Array<HCurlFiniteElement<2>::IPData> FE_TNedelecQuad<ORDER,ZORDER>::ipdata;
  /*
    template <int ORDER, int ZORDER>
    Mat<FE_TNedelecQuadTraits<ORDER,ZORDER>::NDOF> FE_TNedelecQuad<ORDER,ZORDER>::trans;

    template <int ORDER, int ZORDER>
    Mat<FE_TNedelecQuadTraits<ORDER,ZORDER>::NEDGEDOF> FE_TNedelecQuad<ORDER,ZORDER>::trans2;
  */
  template <int ORDER, int ZORDER>
  Matrix<> FE_TNedelecQuad<ORDER,ZORDER>::trans;

  template <int ORDER, int ZORDER>
  Matrix<> FE_TNedelecQuad<ORDER,ZORDER>::trans2;

  template <int ORDER, int ZORDER>
  FE_TNedelecQuad<ORDER,ZORDER> :: FE_TNedelecQuad()
    : HCurlFiniteElement<2>(NDOF, MAXORDER)
  {
    // if (!ipdata.Size())
    Orthogonalize();
    // CalcIPData(ipdata);
    // controldestruction.NewControl(&ipdatadestructed);
  }

  template <int ORDER, int ZORDER>
  FE_TNedelecQuad<ORDER,ZORDER> :: ~FE_TNedelecQuad()
  {
    // if(!ipdatadestructed)
    // ipdata.DeleteAll();
  }

  template <int ORDER, int ZORDER>
  void FE_TNedelecQuad<ORDER,ZORDER> :: 
  CalcShape (const IntegrationPoint & ip, 
	     SliceMatrix<> shape) const
  {
    int i, j;
    shape = 0.0; //!

    Mat<NDOF,2> hshape;
    CalcShape1 (ip, hshape);
    shape = Trans (trans) * hshape;

    Mat<NEDGEDOF,2> shape2, hshape2;
    CalcShape2 (ip, hshape2);
    shape2 = Trans (trans2) * hshape2;

    for (i = 0; i < NEDGEDOF; i++)
      for (j = 0; j < 2; j++)
	shape(i+4, j) = shape2(i,j);

    Mat<4,2> loshape;
    quad1.CalcShape (ip, loshape);
    for (i = 0; i < 4; i++)
      for (j = 0; j < 2; j++)
	shape(i, j) = loshape(i,j);
  }



  template <int ORDER, int ZORDER>
  void FE_TNedelecQuad<ORDER, ZORDER> ::
  CalcShape1 (const IntegrationPoint & ip, 
	      FlatMatrixFixWidth<2> shape) const
  {
    double x = ip(0);
    double y = ip(1);

    Vec<ORDER+1> polx;
    Vec<ZORDER+1> poly;
    
    int i, j;
    polx(0) = 1;
    for (i = 0; i < ORDER; i++)
      polx(i+1) = x * polx(i);
    poly(0) = 1;
    for (i = 0; i < ZORDER; i++)
      poly(i+1) = y * poly(i);

    shape = 0;

    int ii = 0;
    for (i = 0; i < ORDER; i++)
      for (j = 0; j < ZORDER+1; j++)
	shape(ii++, 0) = polx(i) * poly(j);
    for (i = 0; i < ORDER+1; i++)
      for (j = 0; j < ZORDER; j++)
	shape(ii++, 1) = polx(i) * poly(j);
  }





  template <int ORDER, int ZORDER>
  void FE_TNedelecQuad<ORDER, ZORDER> ::
  CalcShape2 (const IntegrationPoint & ip, 
	      FlatMatrixFixWidth<2> shape) const
  {
    // \nabla of edge-bubble


    double x = ip(0);
    double y = ip(1);

    Vec<ORDER+1> polx, dpolx;
    Vec<ZORDER+1> poly, dpoly;
    
    int i, j;
    polx(0) = 1;
    for (i = 0; i < ORDER; i++)
      polx(i+1) = x * polx(i);
    poly(0) = 1;
    for (i = 0; i < ZORDER; i++)
      poly(i+1) = y * poly(i);

    dpolx(0) = 0;
    for (i = 1; i < ORDER+1; i++)
      dpolx(i) = i * polx(i-1);
    dpoly(0) = 0;
    for (i = 1; i < ZORDER+1; i++)
      dpoly(i) = i * poly(i-1);
    
    shape = 0;
    int ii = 0;
    for (i = 0; i < ORDER-1; i++)
      {
	// \nabla x * (1-x) * polx * y
	shape(ii, 0) = (dpolx(i) * x*(1-x) + polx(i) * (1-2*x)) * y;
	shape(ii, 1) = polx(i) * x * (1-x);
	ii++;
	// \nabla x * (1-x) * polx * (1-y)
	shape(ii, 0) = (dpolx(i) * x*(1-x) + polx(i) * (1-2*x)) * (1-y);
	shape(ii, 1) = -polx(i) * x * (1-x);
	ii++;
      }
    for (j = 0; j < ZORDER-1; j++)
      {
	// \nabla x * y * (1-y) * poly
	shape(ii, 0) = poly(j) * y * (1-y);
	shape(ii, 1) = x * (dpoly(j) * y*(1-y) + poly(j) * (1-2*y));
	ii++;

	// \nabla (1-x) * y * (1-y) * poly
	shape(ii, 0) = -poly(j) * y * (1-y);
	shape(ii, 1) = (1-x) * (dpoly(j) * y*(1-y) + poly(j) * (1-2*y));
	ii++;
      }
  }


  template <int ORDER, int ZORDER>
  void FE_TNedelecQuad<ORDER, ZORDER> :: Orthogonalize()
  {
    int i, j, k;

    Mat<NDOF> fiphij;

    FE_TSegmL2<MAXORDER-1> segm;
    Mat<MAXORDER+1, NDOF> edgemoments;
    
    // horizontal+vertical edges
    int base = 4;
    for (i = 0; i < 4; i++)
      {
	int nedge = (i < 2) ? ORDER : ZORDER;

	// OJEOJEOJE
	ComputeEdgeMoments (i, segm, edgemoments, 2*MAXORDER);

	for (j = 0; j < NDOF; j++)
	  {
	    fiphij(i,j) = edgemoments(0,j);  // lowest order
	    for (k = 1; k < nedge; k++)
	      fiphij(base+k-1, j) = edgemoments(k,j);
	  }
	base += nedge-1;
      }

    FE_TFaceTest<ORDER,ZORDER> facetest;
    enum { NTEST = FE_TFaceTest<ORDER,ZORDER>::NDOF };
    Mat<NTEST,NDOF> facemoments;

    ComputeFaceMoments (0, facetest, facemoments, 2*MAXORDER);
    
    for (j = 0; j < NDOF; j++)
      for (k = 0; k < NTEST; k++)
	fiphij(base+k, j) = facemoments(k,j);

    trans.SetSize (NDOF, NDOF);
    CalcInverse (fiphij, trans);




    Mat<NEDGEDOF> fiphij2;

    // horizontal+vertical edges
    base = 0;
    for (i = 0; i < 4; i++)
      {
	int nedge = (i < 2) ? ORDER : ZORDER;

	ComputeEdgeMoments (i, segm, edgemoments, 2*MAXORDER, 2);
	nedge--;
	for (j = 0; j < NEDGEDOF; j++)
	  for (k = 0; k < nedge; k++)
	    fiphij2(base+k, j) = edgemoments(k+1,j);

	base += nedge;
      }

    trans2.SetSize (NEDGEDOF, NEDGEDOF);
    CalcInverse (fiphij2, trans2);
  }
  



  template class  FE_TNedelecQuad<1,2>;
  template class  FE_TNedelecQuad<1,3>;
  template class  FE_TNedelecQuad<2,1>;
  template class  FE_TNedelecQuad<2,2>;
  template class  FE_TNedelecQuad<2,3>;
  template class  FE_TNedelecQuad<2,4>;
  template class  FE_TNedelecQuad<3,1>;
  template class  FE_TNedelecQuad<3,2>;
  template class  FE_TNedelecQuad<3,3>;




















  /* ******************** Tet elements *********************** */
  
#ifdef OLD
  // Array<HCurlFiniteElement<3>::IPData> FE_NedelecTet1o::ipdata;
  

  FE_NedelecTet1o :: FE_NedelecTet1o()
    : HCurlFiniteElement<3> (ET_TET, NDOF, 1)
  {
    // CalcIPData(ipdata);
    // controldestruction.NewControl(&ipdatadestructed);
  }
  
  FE_NedelecTet1o :: ~FE_NedelecTet1o()
  {
    // if(!ipdatadestructed)
    // ipdata.DeleteAll();
  }
  
  void FE_NedelecTet1o :: 
  CalcShape (const IntegrationPoint & ip, 
	     FlatMatrixFixWidth<3> shape) const
  {
    shape = 0.0; //!
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);

    shape (0,0) = 1 - z - y;
    shape (1,0) = y;
    shape (2,0) = z;
    shape (3,0) = -y;
    shape (4,0) = -z;
    shape (5,0) = 0;
    
    shape (0,1) = x;
    shape (1,1) = 1 - x - z;
    shape (2,1) = z;
    shape (3,1) = x;
    shape (4,1) =  0;
    shape (5,1) = -z;      
    
    shape (0,2) = x;
    shape (1,2) = y;
    shape (2,2) = 1 - y - x;
    shape (3,2) = 0;
    shape (4,2) = x;
    shape (5,2) = y;      
  }


  void FE_NedelecTet1o :: 
  CalcCurlShape (const IntegrationPoint & ip, 
		 FlatMatrixFixWidth<3> curlshape) const
  {
    curlshape (0,0) = 0;
    curlshape (0,1) = -2;
    curlshape (0,2) = 2;

    curlshape (1,0) = 2;
    curlshape (1,1) = 0;
    curlshape (1,2) = -2;

    curlshape (2,0) = -2;
    curlshape (2,1) = 2;
    curlshape (2,2) = 0;

    curlshape (3,0) = 0;
    curlshape (3,1) = 0;
    curlshape (3,2) = 2;

    curlshape (4,0) = 0;
    curlshape (4,1) = -2;
    curlshape (4,2) = 0;

    curlshape (5,0) = 2;
    curlshape (5,1) = 0;
    curlshape (5,2) = 0;
  }




  // Array<HCurlFiniteElement<3>::IPData> FE_NedelecTet2o::ipdata;
  Mat<FE_NedelecTet2o::NDOF> FE_NedelecTet2o::trans;

  FE_NedelecTet2o :: FE_NedelecTet2o()
    : HCurlFiniteElement<3> (ET_TET, NDOF, 1)
  {
    // if (!ipdata.Size())
      Orthogonalize();
      // CalcIPData(ipdata);
      // controldestruction.NewControl(&ipdatadestructed);
  }

  FE_NedelecTet2o :: ~FE_NedelecTet2o()
  {
    // if(!ipdatadestructed)
    // ipdata.DeleteAll();
  }

  void FE_NedelecTet2o :: CalcShape (const IntegrationPoint & ip, 
                                     FlatMatrixFixWidth<3> shape) const
  {
    shape = 0.0; //!
    Mat<NDOF,3> shape1;
    CalcShape1 (ip, shape1);
    shape = Trans (trans) * shape1;
  }

  void FE_NedelecTet2o :: 
  CalcShape1 (const IntegrationPoint & ip, 
	      FlatMatrixFixWidth<3> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);
    shape = 0.0;
    for (int i = 0; i < 3; i++)
      {
	int base = 4 * i;
	shape (base  , i) = 1;
	shape (base+1, i) = x;
	shape (base+2, i) = y;
	shape (base+3, i) = z;
      }
  }

  void FE_NedelecTet2o :: Orthogonalize()
  {
    int i, j;

    Mat<NDOF> fiphij;

    Mat<2,NDOF> edgemoments;
    FE_Segm1L2 segm1;
    
    for (i = 0; i < 6; i++)
      {
	ComputeEdgeMoments (i, segm1, edgemoments, 2);

	for (j = 0; j < NDOF; j++)
	  {
	    fiphij(i, j)    = edgemoments(0, j);
	    fiphij(6+i, j)  = edgemoments(1, j);
	  }
      }

    CalcInverse (fiphij, trans);
  }
  
#endif






  class FE_Tet3EdgeBubble : public ScalarFiniteElement<3>
  {
  public:
    ///
    FE_Tet3EdgeBubble()
      : ScalarFiniteElement<3> (12, 3) { ; }

    virtual ELEMENT_TYPE ElementType() const override { return ET_TET; }

    using ScalarFiniteElement<3>::CalcShape;
    virtual void CalcShape (const IntegrationPoint & ip, 
			    BareSliceVector<> shape) const override
    {
      shape.Range(0,ndof) = 0.0; //!
      double x = ip(0);
      double y = ip(1);
      double z = ip(2);
      double l4 = 1-x-y-z;

      double c1 = -3.0;
      double c2 = 7.5;

      shape(0) = c1 * x * l4;
      shape(1) = c1 * y * l4;
      shape(2) = c1 * z * l4;
      shape(3) = c1 * x * y;
      shape(4) = c1 * x * z;
      shape(5) = c1 * y * z;

      shape(6) = c2 * l4 * (l4-x) * x;
      shape(7) = c2 * l4 * (l4-y) * y;
      shape(8) = c2 * l4 * (l4-z) * z;
      shape(9) = c2 * y * (x-y) * x;
      shape(10) = c2 * z * (x-z) * x;
      shape(11) = c2 * z * (y-z) * y;
    }

    using ScalarFiniteElement<3>::CalcDShape;
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     BareSliceMatrix<> dshape) const override
    {
      double x = ip(0);
      double y = ip(1);
      double z = ip(2);
      double l4 = 1-x-y-z;

      double c1 = -3.0;
      double c2 = 7.5;

      // shape(0) = c1 * x * l4;
      dshape(0,0) = c1 * (l4 - x);
      dshape(0,1) = c1 * (-x);
      dshape(0,2) = c1 * (-x);

      //      shape(1) = c1 * y * l4;
      dshape(1,0) = c1 * (-y);
      dshape(1,1) = c1 * (l4 - y);
      dshape(1,2) = c1 * (-y);

      // shape(2) = c1 * z * l4;
      dshape(2,0) = c1 * (-z);
      dshape(2,1) = c1 * (-z);
      dshape(2,2) = c1 * (l4-z);

      //      shape(3) = c1 * x * y;
      dshape(3,0) = c1 * y;
      dshape(3,1) = c1 * x;
      dshape(3,2) = 0;

      //      shape(4) = c1 * x * z;
      dshape(4,0) = c1 * z;
      dshape(4,1) = 0;
      dshape(4,2) = c1 * x;

      // shape(5) = c1 * y * z;
      dshape(5,0) = 0;
      dshape(5,1) = c1 * z;
      dshape(5,2) = c1 * y;

      // shape(6) = c2 * l4 * (l4-x) * x;
      dshape(6,0) = c2 * (-2*l4*x+l4*l4-2*x*l4+x*x);
      dshape(6,1) = c2 * (-2*l4*x+x*x);
      dshape(6,2) = c2 * (-2*l4*x+x*x);

      // shape(7) = c2 * l4 * (l4-y) * y;
      dshape(7,0) = c2 * (-2*l4*y+y*y);
      dshape(7,1) = c2 * (-2*l4*y+l4*l4-2*y*l4+y*y);
      dshape(7,2) = c2 * (-2*l4*y+y*y);

      // shape(8) = c2 * l4 * (l4-z) * z;
      dshape(8,0) = c2 * (-2*l4*z+z*z);
      dshape(8,1) = c2 * (-2*l4*z+z*z);
      dshape(8,2) = c2 * (-2*l4*z+l4*l4-2*z*l4+z*z);

      // shape(9) = c2 * y * (x-y) * x;
      dshape(9,0) = c2 * (2*x*y-y*y);
      dshape(9,1) = c2 * (x*x-2*x*y);
      dshape(9,2) = 0;

      // shape(10) = c2 * z * (x-z) * x;
      dshape(10,0) = c2 * (2*x*z-z*z);
      dshape(10,1) = 0;
      dshape(10,2) = c2 * (x*x-2*x*z);

      // shape(11) = c2 * z * (y-z) * y;
      dshape(11,0) = 0;
      dshape(11,1) = c2 * (2*y*z-z*z);
      dshape(11,2) = c2 * (y*y-2*y*z);
    }

    ///
    // virtual const Array<IPData> & GetIPData () const { return ipdata; }
  }; 


  /*
    template<typename Tx, typename TFA>  
    void FE_NedelecTet3 :: T_CalcShape (Tx hx[3], TFA & shape)
    {
    Tx x = hx[0], y = hx[1], z = hx[2];
    Tx lami[4] = { x, y, z, 1-x-y-z };
    
    const EDGE * edges = ElementTopology::GetEdges (ET_TET);
    for (int i = 0; i < 6; i++)
    {
    Tx lam1 = lami[edges[i][0]];
    Tx lam2 = lami[edges[i][1]];
    shape[i] = uDv_minus_vDu<3> (lam1, lam2);
    shape[i+6] = Du<3> (lam1*lam2);
    shape[i+12] = Du<3> (lam1*lam2*(lam1-lam2));
    }

    const FACE * faces = ElementTopology::GetFaces (ET_TET); 
    for (int i = 0; i < 4; i++)
    for (int k = 0; k < 3; k++)
    {
    int k1 = (k+1)%3, k2 = (k+2)%3;
    shape[18+3*i+k] = uDv_minus_vDu<3> (lami[faces[i][k]],
    lami[faces[i][k1]]*lami[faces[i][k2]]);
    }
    }

    template class T_HCurlFiniteElement<FE_NedelecTet3,ET_TET,30,2>;
  */

  /*
    Array<HCurlFiniteElement<3>::IPData> FE_NedelecTet3::ipdata;
    Mat<FE_NedelecTet3::NDOF> FE_NedelecTet3::trans;
    Mat<FE_NedelecTet3::NEDGEDOF> FE_NedelecTet3::trans2;
    Mat<FE_NedelecTet3::NFACEDOF> FE_NedelecTet3::trans3;

    FE_NedelecTet3 :: FE_NedelecTet3()
    : HCurlFiniteElement<3> (ET_TET, NDOF, 2)
    {
    if (!ipdata.Size())
    Orthogonalize();
    CalcIPData(ipdata);
    controldestruction.NewControl(&ipdatadestructed);
    }

    FE_NedelecTet3 :: ~FE_NedelecTet3()
    {
    if(!ipdatadestructed)
    ipdata.DeleteAll();
    }


    void FE_NedelecTet3 :: 
    CalcShape (const IntegrationPoint & ip, 
    FlatMatrixFixWidth<3> shape) const
    {
    shape = 0.0; //!
    FlatMatrixFixWidth<3> tet1shape(6, &shape(0,0));
    tet1.FE_NedelecTet1::CalcShape (ip, tet1shape);

    FlatMatrixFixWidth<3> shape2(NEDGEDOF, &shape(6,0));
    FE_NedelecTet3::CalcShape2 (ip, shape2);

    
    Mat<NFACEDOF,3> hshape3;
    FE_NedelecTet3::CalcShape3 (ip, hshape3);
    for (int i = 0; i < 12; i+=3)
    {
    for (int j = 0; j < 3; j++)
    for (int k = 0; k < 3; k++)
    shape(18+i+j,k) = 
    trans3(i  ,i+j) * hshape3(i  ,k) + 
    trans3(i+1,i+j) * hshape3(i+1,k) + 
    trans3(i+2,i+j) * hshape3(i+2,k);
    }
    }


    void FE_NedelecTet3 :: 
    CalcCurlShape (const IntegrationPoint & ip, 
    FlatMatrixFixWidth<3> curlshape) const
    {
    FlatMatrixFixWidth<3> tet1curlshape(6, &curlshape(0,0));
    tet1.FE_NedelecTet1::CalcCurlShape (ip, tet1curlshape);

    FlatMatrixFixWidth<3> curlshape2(NEDGEDOF, &curlshape(6,0));
    curlshape2 = 0;
    
    Mat<NFACEDOF,3> hcurlshape3;
    FE_NedelecTet3::CalcCurlShape3 (ip, hcurlshape3);
    for (int i = 0; i < 12; i+=3)
    {
    for (int j = 0; j < 3; j++)
    for (int k = 0; k < 3; k++)
    curlshape(18+i+j,k) = 
    trans3(i  ,i+j) * hcurlshape3(i  ,k) + 
    trans3(i+1,i+j) * hcurlshape3(i+1,k) + 
    trans3(i+2,i+j) * hcurlshape3(i+2,k);
    }
    }



    void FE_NedelecTet3 :: CalcShape1 (const IntegrationPoint & ip, 
    FlatMatrixFixWidth<3> shape) const
    {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);
    shape = 0.0;
    for (int i = 0; i < 3; i++)
    {
    int base = 10 * i;
    shape (base  , i) = 1;
    shape (base+1, i) = x;
    shape (base+2, i) = y;
    shape (base+3, i) = z;
    shape (base+4, i) = x*x;
    shape (base+5, i) = x*y;
    shape (base+6, i) = x*z;
    shape (base+7, i) = y*y;
    shape (base+8, i) = y*z;
    shape (base+9, i) = z*z;
    }
    }

    // base functions which are gradients of edge shape functions
    void FE_NedelecTet3 :: 
    CalcShape2 (const IntegrationPoint & ip, 
    FlatMatrixFixWidth<3> shape) const
    {
    FE_Tet3EdgeBubble febub;
    febub.CalcDShape (ip, FlatMatrix<> (shape));
    }


    // face shape functions
    void FE_NedelecTet3 :: 
    CalcShape3 (const IntegrationPoint & ip, 
    FlatMatrixFixWidth<3> shape) const
    {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);
    double l4 = 1 - x - y - z;
    shape = 0.0;

    shape(0,1) = z * l4;   // face1
    shape(1,2) = y * l4;   // face1
    shape(2,0) = y * z;    // face1
    shape(2,1) = y * z;
    shape(2,2) = y * z;

    shape(3,0) = z * l4;   // face2
    shape(4,2) = x * l4;   // face2
    shape(5,0) = x * z;    // face2
    shape(5,1) = x * z;
    shape(5,2) = x * z;

    shape(6,0) = y * l4;   // face3 
    shape(7,1) = x * l4;   // face3
    shape(8,0) = x * y;    // face3
    shape(8,1) = x * y;    
    shape(8,2) = x * y;

    shape(9 ,0) = y * z;    // face4
    shape(10,1) = x * z;    // face4
    shape(11,2) = x * y;    // face4
    }




    // face shape functions
    void FE_NedelecTet3 :: 
    CalcCurlShape3 (const IntegrationPoint & ip, 
    FlatMatrixFixWidth<3> curlshape) const
    {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);
    double l4 = 1 - x - y - z;

    // shape(0,1) = z * l4;   // face1
    curlshape(0,0) = -l4 + z;
    curlshape(0,1) = 0;
    curlshape(0,2) = -z;

    // shape(1,2) = y * l4;   // face1
    curlshape(1,0) = l4-y;
    curlshape(1,1) = y;
    curlshape(1,2) = 0;

    // shape(2,0) = y * z;    // face1
    // shape(2,1) = y * z;
    // shape(2,2) = y * z;
    curlshape(2,0) = z-y;
    curlshape(2,1) = y;
    curlshape(2,2) = -z;

    //shape(3,0) = z * l4;   // face2
    curlshape(3,0) = 0;
    curlshape(3,1) = l4-z;
    curlshape(3,2) = z;

    // shape(4,2) = x * l4;   // face2
    curlshape(4,0) = -x;
    curlshape(4,1) = -l4+x;
    curlshape(4,2) = 0;

    // shape(5,0) = x * z;    // face2
    // shape(5,1) = x * z;
    // shape(5,2) = x * z;
    curlshape(5,0) = -x;
    curlshape(5,1) = x-z;
    curlshape(5,2) = z;
    
    // shape(6,0) = y * l4;   // face3 
    curlshape(6,0) = 0;
    curlshape(6,1) = -y;
    curlshape(6,2) = -l4+y;

    //    shape(7,1) = x * l4;   // face3
    curlshape(7,0) = x;
    curlshape(7,1) = 0;
    curlshape(7,2) = l4-x;

    // shape(8,0) = x * y;    // face3
    // shape(8,1) = x * y;    
    // shape(8,2) = x * y;
    curlshape(8,0) = x;
    curlshape(8,1) = -y;
    curlshape(8,2) = y-x;

    // shape(9 ,0) = y * z;    // face4
    curlshape(9,0) = 0;
    curlshape(9,1) = y;
    curlshape(9,2) = -z;
    
    // shape(10,1) = x * z;    // face4
    curlshape(10,0) = -x;
    curlshape(10,1) = 0;
    curlshape(10,2) = z;

    // shape(11,2) = x * y;    // face4
    curlshape(11,0) = x;
    curlshape(11,1) = -y;
    curlshape(11,2) = 0;
    }







    void FE_NedelecTet3 :: Orthogonalize()
    {
    int i, j;
    
    Mat<NDOF> fiphij;

    Mat<3,NDOF> edgemoments;
    FE_Segm2L2 segm2;
    
    for (i = 0; i < 6; i++)
    {
    ComputeEdgeMoments (i, segm2, edgemoments, 4);
    for (j = 0; j < NDOF; j++)
    {
    fiphij(i, j)    = edgemoments(0, j);
    fiphij(6+i, j)  = edgemoments(1, j);
    fiphij(12+i, j) = edgemoments(2, j);
    }
    }

    Mat<3,NDOF> facemoments;
    FE_RTTrig0 rttrig0;
    
    for (i = 0; i < 4; i++)
    {
    ComputeFaceMoments (i, rttrig0, facemoments, 4);
    for (j = 0; j < NDOF; j++)
    {
    fiphij(18+3*i, j) =  facemoments(1, j);
    fiphij(19+3*i, j) =  facemoments(0, j);
    fiphij(20+3*i, j) =  facemoments(2, j);
    }
    }
    
    CalcInverse (fiphij, trans);


    // curl-free edge shape functions:
    // int nd = NEDGEDOF;
    Mat<NEDGEDOF> fiphij2;
    Matrix<> edgemoments2(3, NEDGEDOF);
    
    for (i = 0; i < 6; i++)
    {
    ComputeEdgeMoments (i, segm2, edgemoments2, 4, 2);

    for (j = 0; j < NEDGEDOF; j++)
    {
    fiphij2(i, j) = edgemoments2(1, j);
    fiphij2(6+i, j) = edgemoments2(2, j);
    }
    }

    CalcInverse (fiphij2, trans2);
    // should be Id


    
    Mat<NFACEDOF> fiphij3;
    Mat<3,NFACEDOF> facemoments3;
    
    for (i = 0; i < 4; i++)
    {
    ComputeFaceMoments (i, rttrig0, facemoments3, 4, 3);
    for (j = 0; j < NFACEDOF; j++)
    {
    fiphij3(3*i  , j) =  facemoments3(1, j);
    fiphij3(3*i+1, j) =  facemoments3(0, j);
    fiphij3(3*i+2, j) =  facemoments3(2, j);
    // for (int k = 0; k < 3; k++)
    // fiphij3(3*i+k, j) = facemoments3(k, j);
    }
    }
    
    CalcInverse (fiphij3, trans3);
    // should be block diagonal
    }
  */

  /* ******************** Hex Elements ************************* */ 

  // Array<HCurlFiniteElement<3>::IPData> FE_NedelecHex1::ipdata;
  
  FE_NedelecHex1 :: FE_NedelecHex1()
    : HCurlFiniteElement<3> (12, 1)
  {
    order++;
  }

  FE_NedelecHex1 :: ~FE_NedelecHex1()
  {
    // if(!ipdatadestructed)
    // ipdata.DeleteAll();
  }

  void FE_NedelecHex1 :: 
  CalcShape (const IntegrationPoint & ip, 
	     SliceMatrix<> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);

    shape = 0;

    shape (0,0) = (1-z) * (1 - y);
    shape (1,0) = (1-z) * (-y);
    shape (2,1) = (1-z) * (-1+x);
    shape (3,1) = (1-z) * x; 

    shape (4,0) = z * (1 - y);
    shape (5,0) = z * (-y);
    shape (6,1) = z * (-1+x);
    shape (7,1) = z * x; 
    
    shape (8,2)  = (1-x)* (1-y); 
    shape (9,2)  = x    * (1-y); 
    shape (10,2) = x    * y; 
    shape (11,2) = (1-x)* y;
  }











  /* **************************** Tet3 without gradients ******************* */


  // Array<HCurlFiniteElement<3>::IPData> FE_NedelecTet3NoGrad::ipdata;
  Mat<FE_NedelecTet3NoGrad::NFACEDOF> FE_NedelecTet3NoGrad::trans3;

  FE_NedelecTet3NoGrad :: FE_NedelecTet3NoGrad()
    : HCurlFiniteElement<3> (NDOF, 2)
  {
    // if (!ipdata.Size())
      Orthogonalize();
      // CalcIPData(ipdata);
      // controldestruction.NewControl(&ipdatadestructed);
  }

  FE_NedelecTet3NoGrad :: ~FE_NedelecTet3NoGrad()
  {
    // if(!ipdatadestructed)
    // ipdata.DeleteAll();
  }


  void FE_NedelecTet3NoGrad :: 
  CalcShape (const IntegrationPoint & ip, 
	     SliceMatrix<> shape) const
  {
    shape = 0.0; //!
    FlatMatrixFixWidth<3> tet1shape(6, &shape(0,0));
    tet1.FE_NedelecTet1::CalcShape (ip, tet1shape);
    
    Mat<NFACEDOF,3> hshape3;
    FE_NedelecTet3NoGrad::CalcShape3 (ip, hshape3);
    for (int i = 0; i < 12; i+=3)
      {
	for (int j = 0; j < 3; j++)
	  for (int k = 0; k < 3; k++)
	    shape(6+i+j,k) = 
	      trans3(i  ,i+j) * hshape3(i  ,k) + 
	      trans3(i+1,i+j) * hshape3(i+1,k) + 
	      trans3(i+2,i+j) * hshape3(i+2,k);
      }
  }


  void FE_NedelecTet3NoGrad :: 
  CalcCurlShape (const IntegrationPoint & ip, 
		 FlatMatrixFixWidth<3> curlshape) const
  {
    FlatMatrixFixWidth<3> tet1curlshape(6, &curlshape(0,0));
    tet1.FE_NedelecTet1::CalcCurlShape (ip, tet1curlshape);

    Mat<NFACEDOF,3> hcurlshape3;
    FE_NedelecTet3NoGrad::CalcCurlShape3 (ip, hcurlshape3);
    for (int i = 0; i < 12; i+=3)
      {
	for (int j = 0; j < 3; j++)
	  for (int k = 0; k < 3; k++)
	    curlshape(6+i+j,k) = 
	      trans3(i  ,i+j) * hcurlshape3(i  ,k) + 
	      trans3(i+1,i+j) * hcurlshape3(i+1,k) + 
	      trans3(i+2,i+j) * hcurlshape3(i+2,k);
      }
  }


  // face shape functions
  void FE_NedelecTet3NoGrad :: 
  CalcShape3 (const IntegrationPoint & ip, 
	      FlatMatrixFixWidth<3> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);
    double l4 = 1 - x - y - z;
    shape = 0.0;

    shape(0,1) = z * l4;   // face1
    shape(1,2) = y * l4;   // face1
    shape(2,0) = y * z;    // face1
    shape(2,1) = y * z;
    shape(2,2) = y * z;

    shape(3,0) = z * l4;   // face2
    shape(4,2) = x * l4;   // face2
    shape(5,0) = x * z;    // face2
    shape(5,1) = x * z;
    shape(5,2) = x * z;

    shape(6,0) = y * l4;   // face3 
    shape(7,1) = x * l4;   // face3
    shape(8,0) = x * y;    // face3
    shape(8,1) = x * y;    
    shape(8,2) = x * y;

    shape(9 ,0) = y * z;    // face4
    shape(10,1) = x * z;    // face4
    shape(11,2) = x * y;    // face4
  }


  // face shape functions
  void FE_NedelecTet3NoGrad :: 
  CalcCurlShape3 (const IntegrationPoint & ip, 
		  FlatMatrixFixWidth<3> curlshape) const
  {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);
    double l4 = 1 - x - y - z;

    // shape(0,1) = z * l4;   // face1
    curlshape(0,0) = -l4 + z;
    curlshape(0,1) = 0;
    curlshape(0,2) = -z;

    // shape(1,2) = y * l4;   // face1
    curlshape(1,0) = l4-y;
    curlshape(1,1) = y;
    curlshape(1,2) = 0;

    // shape(2,0) = y * z;    // face1
    // shape(2,1) = y * z;
    // shape(2,2) = y * z;
    curlshape(2,0) = z-y;
    curlshape(2,1) = y;
    curlshape(2,2) = -z;

    //shape(3,0) = z * l4;   // face2
    curlshape(3,0) = 0;
    curlshape(3,1) = l4-z;
    curlshape(3,2) = z;

    // shape(4,2) = x * l4;   // face2
    curlshape(4,0) = -x;
    curlshape(4,1) = -l4+x;
    curlshape(4,2) = 0;

    // shape(5,0) = x * z;    // face2
    // shape(5,1) = x * z;
    // shape(5,2) = x * z;
    curlshape(5,0) = -x;
    curlshape(5,1) = x-z;
    curlshape(5,2) = z;
    
    // shape(6,0) = y * l4;   // face3 
    curlshape(6,0) = 0;
    curlshape(6,1) = -y;
    curlshape(6,2) = -l4+y;

    //    shape(7,1) = x * l4;   // face3
    curlshape(7,0) = x;
    curlshape(7,1) = 0;
    curlshape(7,2) = l4-x;

    // shape(8,0) = x * y;    // face3
    // shape(8,1) = x * y;    
    // shape(8,2) = x * y;
    curlshape(8,0) = x;
    curlshape(8,1) = -y;
    curlshape(8,2) = y-x;

    // shape(9 ,0) = y * z;    // face4
    curlshape(9,0) = 0;
    curlshape(9,1) = y;
    curlshape(9,2) = -z;
    
    // shape(10,1) = x * z;    // face4
    curlshape(10,0) = -x;
    curlshape(10,1) = 0;
    curlshape(10,2) = z;

    // shape(11,2) = x * y;    // face4
    curlshape(11,0) = x;
    curlshape(11,1) = -y;
    curlshape(11,2) = 0;
  }



  void FE_NedelecTet3NoGrad :: Orthogonalize()
  {
    int i, j;

    FE_RTTrig0 rttrig0;

    /*    
          Mat<NDOF> fiphij;

          Mat<3,NDOF> edgemoments;
          FE_Segm2L2 segm2;
    
          for (i = 0; i < 6; i++)
          {
          ComputeEdgeMoments (i, segm2, edgemoments, 4);
          for (j = 0; j < NDOF; j++)
	  {
          fiphij(i, j)    = edgemoments(0, j);
          fiphij(6+i, j)  = edgemoments(1, j);
          fiphij(12+i, j) = edgemoments(2, j);
	  }
          }

          Mat<3,NDOF> facemoments;
    
          for (i = 0; i < 4; i++)
          {
          ComputeFaceMoments (i, rttrig0, facemoments, 4);
          for (j = 0; j < NDOF; j++)
	  {
          fiphij(18+3*i, j) =  facemoments(1, j);
          fiphij(19+3*i, j) =  facemoments(0, j);
          fiphij(20+3*i, j) =  facemoments(2, j);
	  }
          }
    
          CalcInverse (fiphij, trans);

          // curl-free edge shape functions:
          // int nd = NEDGEDOF;
          Mat<NEDGEDOF> fiphij2;
          Matrix<> edgemoments2(3, NEDGEDOF);
    
          for (i = 0; i < 6; i++)
          {
          ComputeEdgeMoments (i, segm2, edgemoments2, 4, 2);

          for (j = 0; j < NEDGEDOF; j++)
	  {
          fiphij2(i, j) = edgemoments2(1, j);
          fiphij2(6+i, j) = edgemoments2(2, j);
	  }
          }

          CalcInverse (fiphij2, trans2);
          // should be Id
          */


    Mat<NFACEDOF> fiphij3;
    Mat<3,NFACEDOF> facemoments3;
    
    for (i = 0; i < 4; i++)
      {
	ComputeFaceMoments (i, rttrig0, facemoments3, 4, 3);
	for (j = 0; j < NFACEDOF; j++)
	  {
	    fiphij3(3*i  , j) =  facemoments3(1, j);
	    fiphij3(3*i+1, j) =  facemoments3(0, j);
	    fiphij3(3*i+2, j) =  facemoments3(2, j);
	    /*
              for (int k = 0; k < 3; k++)
	      fiphij3(3*i+k, j) = facemoments3(k, j);
	    */
	  }
      }
    
    CalcInverse (fiphij3, trans3);
    // should be block diagonal
  }
  

















  /* ******************** Prism elements *********************** */

  /*
    Array<HCurlFiniteElement<3>::IPData> FE_NedelecPrism1::ipdata;

    FE_NedelecPrism1 :: FE_NedelecPrism1()
    : HCurlFiniteElement<3> (ET_PRISM, 9, 2)
    {
    CalcIPData(ipdata);
    controldestruction.NewControl(&ipdatadestructed);
    }

    FE_NedelecPrism1 :: ~FE_NedelecPrism1()
    {
    if(!ipdatadestructed)
    ipdata.DeleteAll();
    }

    void FE_NedelecPrism1 :: 
    CalcShape (const IntegrationPoint & ip, 
    FlatMatrixFixWidth<3> shape) const
    {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);

    shape = 0;

    shape (0,0) = (1-z) * (1 - y);
    shape (1,0) = (1-z) * (-y);
    shape (2,0) = (1-z) * (y);
    
    shape (3,0) = z * (1 - y);
    shape (4,0) = z * (-y);
    shape (5,0) = z * (y);

    shape (0,1) = (1-z) * x;
    shape (1,1) = (1-z) * x;
    shape (2,1) = (1-z) * (1-x);
    
    shape (3,1) = z * x;
    shape (4,1) = z * x;
    shape (5,1) = z * (1-x);

    shape (6,2) = 1 - x - y;
    shape (7,2) = x;
    shape (8,2) = y;
    }
  */



  // template <int OZ>
  // Array<HCurlFiniteElement<3>::IPData> FE_TNedelecPrism2<OZ>::ipdata;
  template <int OZ>
  Matrix<> FE_TNedelecPrism2<OZ>::trans;
  template <int OZ>
  Matrix<> FE_TNedelecPrism2<OZ>::trans2;
  template <int OZ>
  Matrix<> FE_TNedelecPrism2<OZ>::trans3;

  template <int ZORDER>
  FE_TNedelecPrism2<ZORDER> :: FE_TNedelecPrism2()
    : HCurlFiniteElement<3> (NDOF, MAXORDER)
  {
    // if (!ipdata.Size())
      Orthogonalize();
    
      // CalcIPData(ipdata);
      // controldestruction.NewControl(&ipdatadestructed);
  }

  template <int ZORDER>
  FE_TNedelecPrism2<ZORDER> :: ~FE_TNedelecPrism2()
  {
    // if(!ipdatadestructed)
    // ipdata.DeleteAll();
  }

  template <int ZORDER>
  void FE_TNedelecPrism2<ZORDER> :: 
  CalcShape (const IntegrationPoint & ip, 
	     SliceMatrix<> shape) const
  {
    int i, j;

    shape = 0.0; //!
    /*
      Mat<NDOF,3> hshape;
      CalcShape1 (ip, hshape);
      shape = Trans (trans) * hshape;
    */

    Mat<NEDGEDOF,3> shape2, hshape2;
    CalcShape2 (ip, hshape2);
    shape2 = Trans (trans2) * hshape2;

    for (i = 0; i < NEDGEDOF; i++)
      for (j = 0; j < 3; j++)
	shape(i+9, j) = shape2(i,j);


    Mat<NFACEDOF,3> shape3, hshape3;
    CalcShape3 (ip, hshape3);
    shape3 = Trans (trans3) * hshape3;

    for (i = 0; i < NFACEDOF; i++)
      for (j = 0; j < 3; j++)
	shape(i+9+NEDGEDOF, j) = shape3(i,j);


    //    (*testout) << "shape = " << endl << shape << endl;
    //    (*testout) << "shape3 = " << endl << shape3 << endl;

    Mat<9,3> loshape;
    prism1.CalcShape (ip, loshape);
    for (i = 0; i < 9; i++)
      for (j = 0; j < 3; j++)
	shape(i,j) = loshape(i,j);
  }



  template <int ZORDER>
  void FE_TNedelecPrism2<ZORDER> ::
  CalcShape1 (const IntegrationPoint & ip, 
	      FlatMatrixFixWidth<3> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);

    Vec<ZORDER+1> polz;
    
    int i, j;
    polz(0) = 1;
    for (i = 0; i < ZORDER; i++)
      polz(i+1) = z * polz(i);

    shape = 0;
    int ii = 0;
    for (j = 0; j < ZORDER+1; j++)
      {
	shape(ii++, 0) =     polz(j);
	shape(ii++, 0) = x * polz(j);
	shape(ii++, 0) = y * polz(j);
	shape(ii++, 1) =     polz(j);
	shape(ii++, 1) = x * polz(j);
	shape(ii++, 1) = y * polz(j);
      }

    for (j = 0; j < ZORDER; j++)
      {
	shape(ii++, 2) =     polz(j);
	shape(ii++, 2) = x * polz(j);
	shape(ii++, 2) = y * polz(j);
	shape(ii++, 2) = x*x * polz(j);
	shape(ii++, 2) = x*y * polz(j);
	shape(ii++, 2) = y*y * polz(j);
      }
  }





  template <int ZORDER>
  void FE_TNedelecPrism2<ZORDER> ::
  CalcShape2 (const IntegrationPoint & ip, 
	      FlatMatrixFixWidth<3> shape) const
  {
    // \nabla of edge-bubble

    double x = ip(0);
    double y = ip(1);
    double z = ip(2);

    shape = 0;

    // \nabla of horizontal edges:

    // x y 
    shape(0,0) = y;
    shape(0,1) = x;
    shape(0,2) = 0;

    // x y z
    shape(1,0) = y*z;
    shape(1,1) = x*z;
    shape(1,2) = x*y;

    // x (1-x-y) 
    shape(2,0) = 1-2*x-y;
    shape(2,1) = -x;
    shape(2,2) = 0;

    // x (1-x-y) * z
    shape(3,0) = (1-2*x-y)*z;
    shape(3,1) = -x*z;
    shape(3,2) = x*(1-x-y);

    // y (1-x-y) 
    shape(4,0) = -y;
    shape(4,1) = 1-x-2*y;
    shape(4,2) = 0;

    // y (1-x-y) * z
    shape(5,0) = -y*z;
    shape(5,1) = (1-x-2*y)*z;
    shape(5,2) = y*(1-x-y);

    Vec<ZORDER+1> polz, dpolz;
    
    int i, j;
    polz(0) = 1;
    for (i = 1; i < ZORDER+1; i++)
      polz(i) = z * polz(i-1);
    dpolz(0) = 0;
    for (i = 1; i < ZORDER+1; i++)
      dpolz(i) = i * polz(i-1);
    
    int ii = 6;
    for (j = 0; j < ZORDER-1; j++)
      {
	// \nabla z * (1-z) * polz
	shape(ii, 0) = 0;
	shape(ii, 1) = 0;
	shape(ii, 2) = (dpolz(j) * z*(1-z) + polz(j) * (1-2*z));
	ii++;
	// \nabla x * z * (1-z) * polz
	shape(ii, 0) = z * (1-z) * polz(j);
	shape(ii, 1) = 0;
	shape(ii, 2) = x * (dpolz(j) * z*(1-z) + polz(j) * (1-2*z));
	ii++;
	// \nabla y * (1-z) * polz
	shape(ii, 0) = 0;
	shape(ii, 1) = z * (1-z) * polz(j);
	shape(ii, 2) = y * (dpolz(j) * z*(1-z) + polz(j) * (1-2*z));
	ii++;
      }
  }





  template <int ZORDER>
  void FE_TNedelecPrism2<ZORDER> ::
  CalcShape3 (const IntegrationPoint & ip, 
	      FlatMatrixFixWidth<3> shape) const
  {
    // \nabla of edge-bubble
    int i, j;

    double x = ip(0);
    double y = ip(1);
    double z = ip(2);
    // double l3 = 1-x-y;

    Vec<3> shapev;
    Mat<6,2> shapeh;
    
    shape = 0.0;

    shapev(0) = x*(1-x-y);
    shapev(1) = y*(1-x-y);
    shapev(2) = x*y;

    shapeh = 0.0;
    for (j = 0; j < 2; j++)
      {
	shapeh(3*j  ,j) = 1;
	shapeh(3*j+1,j) = x;
	shapeh(3*j+2,j) = y;
      }
			    
    double pz;
    int ii = 0;

    pz = z*(1-z);
    for (i = 0; i < ZORDER-1; i++)
      {
	for (j = 0; j < 6; j++,ii++)
	  {
	    shape(ii,0) = pz*shapeh(j,0);
	    shape(ii,1) = pz*shapeh(j,1);
	  }
	pz *= (z-0.5);
      }
    pz = 1;
    for (i = 0; i < ZORDER; i++)
      {
	for (j = 0; j < 3; j++,ii++)
	  shape(ii,2) = pz*shapev(j);
	pz *= (z-0.5);
      }
  }










  template <int ZORDER>
  void FE_TNedelecPrism2<ZORDER> :: Orthogonalize()
  {
    int nd = NDOF;
    int i, j, k;

    Matrix<> fiphij(nd);
    FE_TSegmL2<MAXORDER-1> segm;
    Matrix<> edgemoments(MAXORDER+1, nd);
    
    // horizontal+vertical edges
    int base = 9;
    for (i = 0; i < 9; i++)
      {
	int nedge = (i < 6) ? 2 : ZORDER;

	ComputeEdgeMoments (i, segm, edgemoments, 2*MAXORDER);

	nedge--;
	for (j = 0; j < nd; j++)
	  {
	    fiphij(i, j) = edgemoments(0,j);
	    for (k = 0; k < nedge; k++)
	      fiphij(base+k, j) = edgemoments(k+1,j);
	  }
	base += nedge;
      }

    for (i = 2; i < 5; i++)
      {
	FE_TFaceTest<2,ZORDER> facetest;
	int ntest = facetest.GetNDof();
	Matrix<> facemoments(ntest, nd);

	ComputeFaceMoments (i, facetest, facemoments, 2*MAXORDER);
	
	for (j = 0; j < nd; j++)
	  for (k = 0; k < ntest; k++)
	    fiphij(base+k, j) = facemoments(k,j);

	base += ntest;
      }

    trans.SetSize (nd, nd);
    CalcInverse (fiphij, trans);




    nd = NEDGEDOF;
    Matrix<> fiphij2(nd);
    
    // horizontal+vertical edges
    base = 0;
    for (i = 0; i < 9; i++)
      {
	int nedge = (i < 6) ? 2 : ZORDER;

	ComputeEdgeMoments (i, segm, edgemoments, 2*MAXORDER, 2);
	nedge--;
	for (j = 0; j < nd; j++)
	  for (k = 0; k < nedge; k++)
	    fiphij2(base+k, j) = edgemoments(k+1,j);

	base += nedge;
      }

    trans2.SetSize (nd, nd);
    CalcInverse (fiphij2, trans2);



    nd = NFACEDOF;
    Matrix<> fiphij3(nd);
    base = 0;

    for (i = 2; i < 5; i++)
      {
	FE_TFaceTest<2,ZORDER> facetest;
	int ntest = facetest.GetNDof();
	Matrix<> facemoments(ntest, nd);

	ComputeFaceMoments (i, facetest, facemoments, 2*MAXORDER, 3);
	
	for (j = 0; j < nd; j++)
	  for (k = 0; k < ntest; k++)
	    fiphij3(base+k, j) = facemoments(k,j);

	base += ntest;
      }

    (*testout) << "fiphij3 = " << endl << fiphij3 << endl;
    trans3.SetSize (nd, nd);
    CalcInverse (fiphij3, trans3);
  }
  



  template class  FE_TNedelecPrism2<1>;
  template class  FE_TNedelecPrism2<2>;
  template class  FE_TNedelecPrism2<3>;
  template class  FE_TNedelecPrism2<4>;







  
  // RT0 x ZORDER-2  x Q (0, ZORDER-1)
  template <int ZORDER>
  class FE_TVolTest3 : public HDivFiniteElement<3>
  {
  private:
    ///
    // Array<IPData> ipdata;

  public:
    enum { NDOF = 4 * ZORDER - 3 };
    enum { MAXORDER = (1 > ZORDER) ? 1 : ZORDER };

    ///
    FE_TVolTest3()
      : HDivFiniteElement<3> (NDOF, MAXORDER)
    { ; }

    ///
    virtual ~FE_TVolTest3()
    { 
      // ipdata.DeleteAll(); 
    }
    virtual ELEMENT_TYPE ElementType() const { return ET_PRISM; }
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    SliceMatrix<> shape) const
    {
      double x = ip(0);
      double y = ip(1);
      double z = ip(2);

      Mat<3,2> rt0;
      Vec<ZORDER> pz;
      FE_TSegmL2<ZORDER-1> segm;
      
      rt0(0,0) = 1;
      rt0(0,1) = 0;
      rt0(1,0) = 0;
      rt0(1,1) = 1;
      rt0(2,0) = x;
      rt0(2,1) = y;
      
      IntegrationPoint ipz(z, 0, 0, 1);
      segm.CalcShape (ipz, pz);

      shape = 0;
      int i, j, ii = 0;
      for (i = 0; i < 3; i++)
	for (j = 0; j < ZORDER-1; j++)
	  {
	    shape(ii, 0) = rt0(i,0) * pz(j);
	    shape(ii, 1) = rt0(i,1) * pz(j);
	    ii++;
	  }
      for (j = 0; j < ZORDER; j++)
	shape(ii++, 2) = pz(j);
    }
  };
  
  



  /* ************************ Prism3 ******************************* */
  




  FE_Trig3Pot :: FE_Trig3Pot()
    : ScalarFiniteElement<2> (10, 3)
  {
    ;
  }

  FE_Trig3Pot :: ~FE_Trig3Pot()
  {
    ;
  }

  void FE_Trig3Pot :: CalcDShape (const IntegrationPoint & ip, 
				  BareSliceMatrix<> dshape) const
  {
    cerr << "calcdshape not implemnted" << endl;
  }

  void FE_Trig3Pot :: CalcShape (const IntegrationPoint & ip, 
				 BareSliceVector<> shape) const
				 
  {
    double x = ip.Point()[0];
    double y = ip.Point()[1];
    double lam3 = 1-x-y;


    shape(0) = x;
    shape(1) = y;
    shape(2) = lam3;

    /*
      const int edges[3][2] = 
      { { 3, 1 },
      { 3, 2 },
      { 1, 2 } };
    */

    shape(3) = 3 * x * lam3 * (lam3+x);
    shape(4) = 7.5 * x * lam3 * (x-lam3);

    shape(5) = 3 * y * lam3 * (lam3+y);
    shape(6) = 7.5 * y * lam3 * (y-lam3);

    shape(7) = 3 * x * y * (x+y);
    shape(8) = 7.5 * x * y * (y-x);
    shape(9) = 60 * x*y*lam3;  // int_T is 0.5
  }


  // template <int ORDER>
  // Array<HCurlFiniteElement<3>::IPData> FE_TNedelecPrism3<ORDER>::ipdata;
  template <int ORDER>
  Matrix<> FE_TNedelecPrism3<ORDER>::trans;
  template <int ORDER>
  Matrix<> FE_TNedelecPrism3<ORDER>::trans2;
  template <int ORDER>
  Matrix<> FE_TNedelecPrism3<ORDER>::trans_quad;
  template <int ORDER>
  Matrix<> FE_TNedelecPrism3<ORDER>::trans_trig;


  template <int ZORDER>
  FE_TNedelecPrism3<ZORDER> :: FE_TNedelecPrism3()
    : HCurlFiniteElement<3> (NDOF, MAXORDER)
  {
    // if (!ipdata.Size())
      Orthogonalize();

      // CalcIPData(ipdata);
      // controldestruction.NewControl(&ipdatadestructed);
  }

  template <int ZORDER>
  FE_TNedelecPrism3<ZORDER> :: ~FE_TNedelecPrism3()
  {
    // if(!ipdatadestructed)
    // ipdata.DeleteAll();
  }

  template <int ZORDER>
  void FE_TNedelecPrism3<ZORDER> :: 
  CalcShape (const IntegrationPoint & ip, 
	     SliceMatrix<> shape) const
  {

    shape = 0.0; //!
    int i, j;

    /*
      Mat<NDOF,3> hshape;
      CalcShape1 (ip, hshape);
      shape = Trans(trans) * hshape;
    */

    Mat<9,3> loshape;
    prism1.CalcShape (ip, loshape);

    Mat<NEDGEDOF,3> shape2, hshape2;
    CalcShape2 (ip, hshape2);
    shape2 = Trans (trans2) * hshape2;

    Mat<NQUADFACEDOF,3> shape3, hshape3;
    CalcShape3 (ip, hshape3);
    shape3 = Trans (trans_quad) * hshape3;

    Mat<NTRIGFACEDOF+NINNERDOF,3> shape4, hshape4;
    CalcShape4 (ip, hshape4);
    shape4 = Trans (trans_trig) * hshape4;

    for (i = 0; i < 9; i++)
      for (j = 0; j < 3; j++)
	shape(i,j) = loshape(i,j);

    for (i = 0; i < NEDGEDOF; i++)
      for (j = 0; j < 3; j++)
	shape(i+9, j) = shape2(i,j);

    for (i = 0; i < NQUADFACEDOF; i++)
      for (j = 0; j < 3; j++)
	shape(i+9+NEDGEDOF+6, j) = shape3(i,j);

    for (i = 0; i < NTRIGFACEDOF; i++)
      for (j = 0; j < 3; j++)
	shape(i+9+NEDGEDOF, j) = shape4(i,j);

    for (i = 0; i < NINNERDOF; i++)
      for (j = 0; j < 3; j++)
	shape(i+9+NEDGEDOF+6+NQUADFACEDOF, j) = shape4(i+6,j);

    //    (*testout) << "shape = " << endl << shape << endl;
    //    (*testout) << "shape4 = " << endl << shape4 << endl;
  }



  template <int ZORDER>
  void FE_TNedelecPrism3<ZORDER> ::
  CalcShape1 (const IntegrationPoint & ip, 
	      FlatMatrixFixWidth<3> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);

    IntegrationPoint ipxy(x,y,0,1);
    IntegrationPoint ipz(z,0,0,1);

    Vec<6> p2;
    Vec<10> p3;
    Vec<ZORDER+1> pz;

    h1trig2.CalcShape (ipxy, p2);
    h1trig3.CalcShape (ipxy, FlatVector<>(p3));
    segm.CalcShape (ipz, FlatVector<>(pz));
    
    int i, j, ii = 0;
    shape = 0;
    for (i = 0; i < 6; i++)
      for (j = 0; j < ZORDER+1; j++)
	{
	  shape(ii++, 0) = p2(i) * pz(j);
	  shape(ii++, 1) = p2(i) * pz(j);
	}

    for (i = 0; i < 10; i++)
      for (j = 0; j < ZORDER; j++)
	shape(ii++, 2) = p3(i) * pz(j);
  }





  template <int ZORDER>
  void FE_TNedelecPrism3<ZORDER> ::
  CalcShape2 (const IntegrationPoint & ip, 
	      FlatMatrixFixWidth<3> shape) const
  {
    // \nabla of edge-bubble

    double x = ip(0);
    double y = ip(1);
    double z = ip(2);
    shape = 0;

    IntegrationPoint ipxy(x,y,0,1);
    IntegrationPoint ipz(z,0,0,1);
    FE_Trig3EdgeBubble febub;

    // plane edge bubbles:
    Vec<6> edgebub;
    Mat<6,2> gradedgebub;

    febub.CalcShape (ipxy, edgebub);
    febub.CalcDShape (ipxy, gradedgebub);

    int i, j, ii = 0;

    for (i = 0; i < 6; i++)
      {
	shape(ii, 0) = z*gradedgebub(i, 0);
	shape(ii, 1) = z*gradedgebub(i, 1);
	shape(ii, 2) = edgebub(i);  
	ii++;
	shape(ii, 0) = (1-z)*gradedgebub(i, 0);
	shape(ii, 1) = (1-z)*gradedgebub(i, 1);
	shape(ii, 2) = -edgebub(i);  
	ii++;
      }

    
    Vec<ZORDER+1> polz;
    Mat<ZORDER+1,1> dpolz;

    segm.CalcShape (ipz, polz);
    segm.CalcDShape (ipz, dpolz);

    for (j = 0; j < ZORDER-1; j++)
      {
	// \nabla z * (1-z) * polz
	shape(ii, 0) = 0;
	shape(ii, 1) = 0;
	shape(ii, 2) = (dpolz(j,0) * z*(1-z) + polz(j) * (1-2*z));
	ii++;
	// \nabla x * z * (1-z) * polz
	shape(ii, 0) = z * (1-z) * polz(j);
	shape(ii, 1) = 0;
	shape(ii, 2) = x * (dpolz(j,0) * z*(1-z) + polz(j) * (1-2*z));
	ii++;
	// \nabla y * (1-z) * polz
	shape(ii, 0) = 0;
	shape(ii, 1) = z * (1-z) * polz(j);
	shape(ii, 2) = y * (dpolz(j,0) * z*(1-z) + polz(j) * (1-2*z));
	ii++;
      }
    //    (*testout) << "calcshape2, ip = " << ip << endl << "shape: " << endl << shape << endl;
  }




  template <int ZORDER>
  void FE_TNedelecPrism3<ZORDER> ::
  CalcShape3 (const IntegrationPoint & ip, 
	      FlatMatrixFixWidth<3> shape) const
  {
    // \nabla edgebubble * z-bubble

    double x = ip(0);
    double y = ip(1);
    double z = ip(2);

    shape = 0;

    IntegrationPoint ipxy(x,y,0,1);
    IntegrationPoint ipz(z, 0, 0, 1);
    
    FE_Trig3EdgeBubble febub;
    FE_TSegmL2<ZORDER-1> segm;

    Vec<6> edgebub;
    Mat<6,2> gradedgebub;
    Vec<ZORDER> pz;
    Mat<ZORDER,1> dpz;

    febub.CalcShape (ipxy, edgebub);
    febub.CalcDShape (ipxy, gradedgebub);
    segm.CalcShape (ipz, pz);
    segm.CalcDShape (ipz, dpz);


    Mat<3,2> nedelec1;
    nedelec1(0,0) = 1;
    nedelec1(0,1) = 0;
    nedelec1(1,0) = 0;
    nedelec1(1,1) = 1;
    nedelec1(2,0) = y;
    nedelec1(2,1) = -x;


    int i, j, ii = 0;
    for (i = 0; i < 6; i++)
      for (j = 0; j < ZORDER-1; j++)
	{
	  shape(ii, 0) = gradedgebub(i, 0) * pz(j) * z*(z-1);
	  shape(ii, 1) = gradedgebub(i, 1) * pz(j) * z*(z-1);
	  shape(ii, 2) = 0;  
	  ii++;
	}
    for (i = 0; i < 3; i++)
      for (j = 0; j < ZORDER-1; j++)
	{
	  shape(ii, 0) = nedelec1(i, 0) * pz(j) * z*(z-1);
	  shape(ii, 1) = nedelec1(i, 1) * pz(j) * z*(z-1);
	  shape(ii, 2) = 0;  
	  ii++;
	}

    for (i = 0; i < 6; i++)
      for (j = 0; j < ZORDER; j++)
	{
	  shape(ii, 0) = 0;
	  shape(ii, 1) = 0;
	  shape(ii, 2) = edgebub(i) * pz(j); 
	  ii++;
	}
  }




  template <int ZORDER>
  void FE_TNedelecPrism3<ZORDER> ::
  CalcShape4 (const IntegrationPoint & ip, 
	      FlatMatrixFixWidth<3> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);

    shape = 0;
    int i, ii = 0;
    double pz = 1;
    for (i = 0; i < ZORDER+1; i++)
      {
	shape(ii++,0) = pz * y * (1-x-y);
	shape(ii++,1) = pz * x * (1-x-y);
	shape(ii  ,0) = pz * x * y;
	shape(ii++,1) = pz * x * y;
	pz *= z - 0.5;
      }

    pz = 1;
    for (i = 0; i < ZORDER; i++)
      {
	shape(ii++,2) = pz * x * y * (1-x-y);
	pz *= z - 0.5;
      }
  }

  template <int ZORDER>
  void FE_TNedelecPrism3<ZORDER> ::
  CalcInner (const IntegrationPoint & ip, 
	     FlatMatrixFixWidth<3> shape) const
  {
    ;
  }

  template <int ZORDER>
  void FE_TNedelecPrism3<ZORDER> :: GetInternalDofs (Array<int> & idofs) const
  {
    idofs.SetSize (0);
    for (int i = NDOF - NINNERDOF; i < NDOF; i++)
      idofs.Append (i);
  }

  template <int ZORDER>
  void FE_TNedelecPrism3<ZORDER> :: Orthogonalize()
  {
    int nd = NDOF;
    int i, j, k;

    Matrix<> fiphij(nd);
    FE_TSegmL2<MAXORDER-1> segm;
    Matrix<> edgemoments(MAXORDER, nd);
    
    // horizontal+vertical edges
    int base = 9;
    for (i = 0; i < 9; i++)
      {
	int nedge = (i < 6) ? 3 : ZORDER;

	ComputeEdgeMoments (i, segm, edgemoments, 2*MAXORDER);

	nedge--;
	for (j = 0; j < nd; j++)
	  {
	    fiphij(i, j) = edgemoments(0,j);
	    for (k = 0; k < nedge; k++)
	      fiphij(base+k, j) = edgemoments(k+1,j);
	  }
	base += nedge;
      }

    // trig faces:
    Matrix<> trigfacemoments(3,nd);
    FE_RTTrig0 rttrig0;
    for (i = 0; i < 2; i++)
      {
	ComputeFaceMoments (i, rttrig0, trigfacemoments, 6);
	for (j = 0; j < nd; j++)
	  {
	    fiphij(base  , j) =  trigfacemoments(1, j);
	    fiphij(base+1, j) = -trigfacemoments(0, j);
	    fiphij(base+2, j) = -trigfacemoments(2, j);
	  }
	base += 3;
      }

    // quad faces
    for (i = 2; i < 5; i++)
      {
	FE_TFaceTest<3,ZORDER> facetest;
	int ntest = facetest.GetNDof();
	Matrix<> facemoments(ntest, nd);

	ComputeFaceMoments (i, facetest, facemoments, 2*MAXORDER);
	
	for (j = 0; j < nd; j++)
	  for (k = 0; k < ntest; k++)
	    fiphij(base+k, j) = facemoments(k,j);
	base += ntest;
      }

    FE_TVolTest3<ZORDER> voltest;
    int ntest = voltest.GetNDof();
    Matrix<> volmoments(ntest, nd);
    
    ComputeVolMoments (voltest, volmoments, 2*MAXORDER);
    for (j = 0; j < nd; j++)
      for (k = 0; k < ntest; k++)
	fiphij(base+k, j) = volmoments(k,j);
    base += ntest;
    
    trans.SetSize (nd, nd);
    CalcInverse (fiphij, trans);



    nd = NEDGEDOF;
    Matrix<> fiphij2(nd);
    
    // horizontal+vertical edges
    base = 0;
    for (i = 0; i < 9; i++)
      {
	int nedge = (i < 6) ? 3 : ZORDER;

	ComputeEdgeMoments (i, segm, edgemoments, 2*MAXORDER, 2);
	nedge--;
	for (j = 0; j < nd; j++)
	  for (k = 0; k < nedge; k++)
	    fiphij2(base+k, j) = edgemoments(k+1,j);

	base += nedge;
      }

    trans2.SetSize (nd, nd);
    CalcInverse (fiphij2, trans2);



    nd = NQUADFACEDOF;
    Matrix<> fiphij3(nd);
    
    // quad faces
    base = 0;
    for (i = 2; i < 5; i++)
      {
	FE_TFaceTest<3,ZORDER> facetest;
	int ntest = facetest.GetNDof();
	Matrix<> facemoments(ntest, nd);

	ComputeFaceMoments (i, facetest, facemoments, 2*MAXORDER, 3);
	
	for (j = 0; j < nd; j++)
	  for (k = 0; k < ntest; k++)
	    fiphij3(base+k, j) = facemoments(k,j);
	base += ntest;
      }

    trans_quad.SetSize (nd, nd);
    CalcInverse (fiphij3, trans_quad);




    nd = NTRIGFACEDOF+NINNERDOF;
    Matrix<> fiphij4(nd);
    
    // trig faces
    base = 0;
    for (i = 0; i < 2; i++)
      {
	Matrix<> facemoments(3, nd);

	ComputeFaceMoments (i, rttrig0, facemoments, 4, 4);
	for (j = 0; j < nd; j++)
	  {
	    fiphij4(base  , j) = facemoments(1, j);
	    fiphij4(base+1, j) = facemoments(0, j);
	    fiphij4(base+2, j) = facemoments(2, j);
	  }
	base += 3;
      }
    
    {
      FE_TVolTest3<ZORDER> voltest;
      int ntest = voltest.GetNDof();
      Matrix<> volmoments(ntest, nd);
    
      ComputeVolMoments (voltest, volmoments, 2*MAXORDER, 4);
      for (j = 0; j < nd; j++)
	for (k = 0; k < ntest; k++)
	  fiphij4(base+k, j) = volmoments(k,j);
      base += ntest;
    }

    // (*testout) << "fiphij4 = " << endl << fiphij4 << endl;
    trans_trig.SetSize (nd, nd);
    CalcInverse (fiphij4, trans_trig);
  }
  



  template class  FE_TNedelecPrism3<1>;
  template class  FE_TNedelecPrism3<2>;
  template class  FE_TNedelecPrism3<3>;












  
  // RT0 x ZORDER-2  x Q (0, ZORDER-1)
  template <int ZORDER>
  class FE_TVolTest3NoGrad : public HDivFiniteElement<3>
  {
  private:
    ///
    // Array<IPData> ipdata;

  public:
    // enum { NDOF = 4 * ZORDER - 3 };
    enum { NDOF = 3 * ZORDER - 2 };
    enum { MAXORDER = (1 > ZORDER) ? 1 : ZORDER };

    ///
    FE_TVolTest3NoGrad()
      : HDivFiniteElement<3> (NDOF, MAXORDER)
    { ; }

    ///
    virtual ~FE_TVolTest3NoGrad()
    { 
      // ipdata.DeleteAll();
    }

    virtual ELEMENT_TYPE ElementType() const { return ET_PRISM; }
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    SliceMatrix<> shape) const
    {
      double x = ip(0);
      double y = ip(1);
      double z = ip(2);

      Mat<3,2> rt0;
      Vec<ZORDER> pz;
      FE_TSegmL2<ZORDER-1> segm;
      
      rt0(0,0) = 1;
      rt0(0,1) = 0;
      rt0(1,0) = 0;
      rt0(1,1) = 1;
      rt0(2,0) = x;
      rt0(2,1) = y;
      
      IntegrationPoint ipz(z, 0, 0, 1);
      segm.CalcShape (ipz, pz);

      shape = 0;
      int i, j, ii = 0;
      for (i = 0; i < 3; i++)
	for (j = 0; j < ZORDER-1; j++)
	  {
	    shape(ii, 0) = rt0(i,0) * pz(j);
	    shape(ii, 1) = rt0(i,1) * pz(j);
	    ii++;
	  }
      for (j = 0; j < 1; j++)
	shape(ii++, 2) = pz(j);
    }
  };
  
  




  
  /* ************************ Prism3 No Gradients ******************************* */
  
  // template <int ORDER>
  // Array<HCurlFiniteElement<3>::IPData> FE_TNedelecPrism3NoGrad<ORDER>::ipdata;
  template <int ORDER>
  Matrix<> FE_TNedelecPrism3NoGrad<ORDER>::trans_quad;
  template <int ORDER>
  Matrix<> FE_TNedelecPrism3NoGrad<ORDER>::trans_trig;


  template <int ZORDER>
  FE_TNedelecPrism3NoGrad<ZORDER> :: FE_TNedelecPrism3NoGrad()
    : HCurlFiniteElement<3> (NDOF, MAXORDER)
  {
    // if (!ipdata.Size())
      Orthogonalize();

      // CalcIPData(ipdata);
      // controldestruction.NewControl(&ipdatadestructed);
  }

  template <int ZORDER>
  FE_TNedelecPrism3NoGrad<ZORDER> :: ~FE_TNedelecPrism3NoGrad()
  {
    // if(!ipdatadestructed)
    // ipdata.DeleteAll();
  }

  template <int ZORDER>
  void FE_TNedelecPrism3NoGrad<ZORDER> :: 
  CalcShape (const IntegrationPoint & ip, 
	     SliceMatrix<> shape) const
  {
    int i, j;

    shape = 0.0; //!
    Mat<9,3> loshape;
    prism1.CalcShape (ip, loshape);

    /*
      Mat<NEDGEDOF,3> shape2, hshape2;
      CalcShape2 (ip, hshape2);
      shape2 = Trans (trans2) * hshape2;
    */

    Mat<NQUADFACEDOF,3> shape3, hshape3;
    CalcShape3 (ip, hshape3);
    shape3 = Trans (trans_quad) * hshape3;

    Mat<NTRIGFACEDOF+NINNERDOF,3> shape4, hshape4;
    CalcShape4 (ip, hshape4);
    shape4 = Trans (trans_trig) * hshape4;

    for (i = 0; i < 9; i++)
      for (j = 0; j < 3; j++)
	shape(i,j) = loshape(i,j);

    /*
      for (i = 0; i < NEDGEDOF; i++)
      for (j = 0; j < 3; j++)
      shape(i+9, j) = shape2(i,j);
    */

    for (i = 0; i < NQUADFACEDOF; i++)
      for (j = 0; j < 3; j++)
	shape(i+9+6, j) = shape3(i,j);

    for (i = 0; i < NTRIGFACEDOF; i++)
      for (j = 0; j < 3; j++)
	shape(i+9, j) = shape4(i,j);

    for (i = 0; i < NINNERDOF; i++)
      for (j = 0; j < 3; j++)
	shape(i+9+6+NQUADFACEDOF, j) = shape4(i+6,j);

    //    (*testout) << "shape = " << endl << shape << endl;
    //    (*testout) << "shape4 = " << endl << shape4 << endl;
  }



  template <int ZORDER>
  void FE_TNedelecPrism3NoGrad<ZORDER> ::
  CalcShape1 (const IntegrationPoint & ip, 
	      FlatMatrixFixWidth<3> shape) const
  {
    cout << "prism-nograd::calcshape1" << endl;
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);

    IntegrationPoint ipxy(x,y,0,1);
    IntegrationPoint ipz(z,0,0,1);

    Vec<6> p2;
    Vec<10> p3;
    Vec<ZORDER+1> pz;

    h1trig2.CalcShape (ipxy, p2);
    h1trig3.CalcShape (ipxy, FlatVector<>(p3));
    segm.CalcShape (ipz, FlatVector<>(pz));
    
    int i, j, ii = 0;
    shape = 0;
    for (i = 0; i < 6; i++)
      for (j = 0; j < ZORDER+1; j++)
	{
	  shape(ii++, 0) = p2(i) * pz(j);
	  shape(ii++, 1) = p2(i) * pz(j);
	}

    for (i = 0; i < 10; i++)
      for (j = 0; j < ZORDER; j++)
	shape(ii++, 2) = p3(i) * pz(j);
  }





  template <int ZORDER>
  void FE_TNedelecPrism3NoGrad<ZORDER> ::
  CalcShape2 (const IntegrationPoint & ip, 
	      FlatMatrixFixWidth<3> shape) const
  {
    cout << "prism-nograd: calchspae2" << endl;
    // \nabla of edge-bubble

    double x = ip(0);
    double y = ip(1);
    double z = ip(2);
    shape = 0;

    IntegrationPoint ipxy(x,y,0,1);
    IntegrationPoint ipz(z,0,0,1);
    FE_Trig3EdgeBubble febub;

    // plane edge bubbles:
    Vec<6> edgebub;
    Mat<6,2> gradedgebub;

    febub.CalcShape (ipxy, edgebub);
    febub.CalcDShape (ipxy, gradedgebub);

    int i, j, ii = 0;

    for (i = 0; i < 6; i++)
      {
	shape(ii, 0) = z*gradedgebub(i, 0);
	shape(ii, 1) = z*gradedgebub(i, 1);
	shape(ii, 2) = edgebub(i);  
	ii++;
	shape(ii, 0) = (1-z)*gradedgebub(i, 0);
	shape(ii, 1) = (1-z)*gradedgebub(i, 1);
	shape(ii, 2) = -edgebub(i);  
	ii++;
      }
    
    Vec<ZORDER+1> polz;
    Mat<ZORDER+1,1> dpolz;
    
    segm.CalcShape (ipz, polz);
    segm.CalcDShape (ipz, dpolz);
    
    for (j = 0; j < ZORDER-1; j++)
      {
	// \nabla z * (1-z) * polz
	shape(ii, 0) = 0;
	shape(ii, 1) = 0;
	shape(ii, 2) = (dpolz(j,0) * z*(1-z) + polz(j) * (1-2*z));
	ii++;
	// \nabla x * z * (1-z) * polz
	shape(ii, 0) = z * (1-z) * polz(j);
	shape(ii, 1) = 0;
	shape(ii, 2) = x * (dpolz(j,0) * z*(1-z) + polz(j) * (1-2*z));
	ii++;
	// \nabla y * (1-z) * polz
	shape(ii, 0) = 0;
	shape(ii, 1) = z * (1-z) * polz(j);
	shape(ii, 2) = y * (dpolz(j,0) * z*(1-z) + polz(j) * (1-2*z));
	ii++;
      }
    //    (*testout) << "calcshape2, ip = " << ip << endl << "shape: " << endl << shape << endl;
  }




  template <int ZORDER>
  void FE_TNedelecPrism3NoGrad<ZORDER> ::
  CalcShape3 (const IntegrationPoint & ip, 
	      FlatMatrixFixWidth<3> shape) const
  {
    // \nabla edgebubble * z-bubble

    double x = ip(0);
    double y = ip(1);
    double z = ip(2);

    shape = 0;

    IntegrationPoint ipxy(x,y,0,1);
    IntegrationPoint ipz(z, 0, 0, 1);
    
    FE_Trig3EdgeBubble febub;
    FE_TSegmL2<ZORDER-1> segm;

    Vec<6> edgebub;
    Mat<6,2> gradedgebub;
    Vec<ZORDER> pz;
    Mat<ZORDER,1> dpz;

    febub.CalcShape (ipxy, edgebub);
    febub.CalcDShape (ipxy, gradedgebub);
    segm.CalcShape (ipz, pz);
    segm.CalcDShape (ipz, dpz);


    Mat<3,2> nedelec1;
    nedelec1(0,0) = 1;
    nedelec1(0,1) = 0;
    nedelec1(1,0) = 0;
    nedelec1(1,1) = 1;
    nedelec1(2,0) = y;
    nedelec1(2,1) = -x;


    int i, j, ii = 0;
    for (i = 0; i < 6; i++)
      for (j = 0; j < ZORDER-1; j++)
	{
	  shape(ii, 0) = gradedgebub(i, 0) * pz(j) * z*(z-1);
	  shape(ii, 1) = gradedgebub(i, 1) * pz(j) * z*(z-1);
	  shape(ii, 2) = 0;  
	  ii++;
	}
    for (i = 0; i < 3; i++)
      for (j = 0; j < ZORDER-1; j++)
	{
	  shape(ii, 0) = nedelec1(i, 0) * pz(j) * z*(z-1);
	  shape(ii, 1) = nedelec1(i, 1) * pz(j) * z*(z-1);
	  shape(ii, 2) = 0;  
	  ii++;
	}

    for (i = 0; i < 6; i++)
      for (j = 0; j < ZORDER; j++)
	{
	  shape(ii, 0) = 0;
	  shape(ii, 1) = 0;
	  shape(ii, 2) = edgebub(i) * pz(j); 
	  ii++;
	}
  }




  template <int ZORDER>
  void FE_TNedelecPrism3NoGrad<ZORDER> ::
  CalcShape4 (const IntegrationPoint & ip, 
	      FlatMatrixFixWidth<3> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);

    shape = 0;
    int i, ii = 0;
    double pz = 1;
    for (i = 0; i < ZORDER+1; i++)
      {
	shape(ii++,0) = pz * y * (1-x-y);
	shape(ii++,1) = pz * x * (1-x-y);
	shape(ii  ,0) = pz * x * y;
	shape(ii++,1) = pz * x * y;
	pz *= z - 0.5;
      }

    /*
      pz = 1;
      for (i = 0; i < ZORDER; i++)
      {
      shape(ii++,2) = pz * x * y * (1-x-y);
      pz *= z - 0.5;
      }
    */
    shape(ii++,2) = x * y * (1-x-y);
  }

  template <int ZORDER>
  void FE_TNedelecPrism3NoGrad<ZORDER> ::
  CalcInner (const IntegrationPoint & ip, 
	     FlatMatrixFixWidth<3> shape) const
  {
    ;
  }

  template <int ZORDER>
  void FE_TNedelecPrism3NoGrad<ZORDER> :: GetInternalDofs (Array<int> & idofs) const
  {
    idofs.SetSize (0);
    for (int i = NDOF - NINNERDOF; i < NDOF; i++)
      idofs.Append (i);
  }

  template <int ZORDER>
  void FE_TNedelecPrism3NoGrad<ZORDER> :: Orthogonalize()
  {
    int nd = NDOF;
    int i, j, k;

    FE_RTTrig0 rttrig0;
    int base;

    /*
      Matrix<> fiphij(nd);
      FE_TSegmL2<MAXORDER-1> segm;
      Matrix<> edgemoments(MAXORDER, nd);
    
      // horizontal+vertical edges
      int base = 9;
      for (i = 0; i < 9; i++)
      {
      int nedge = (i < 6) ? 3 : ZORDER;

      ComputeEdgeMoments (i, segm, edgemoments, 2*MAXORDER);

      nedge--;
      for (j = 0; j < nd; j++)
      {
      fiphij(i, j) = edgemoments(0,j);
      for (k = 0; k < nedge; k++)
      fiphij(base+k, j) = edgemoments(k+1,j);
      }
      base += nedge;
      }

      // trig faces:
      Matrix<> trigfacemoments(3,nd);
      for (i = 0; i < 2; i++)
      {
      ComputeFaceMoments (i, rttrig0, trigfacemoments, 6);
      for (j = 0; j < nd; j++)
      {
      fiphij(base  , j) =  trigfacemoments(1, j);
      fiphij(base+1, j) = -trigfacemoments(0, j);
      fiphij(base+2, j) = -trigfacemoments(2, j);
      }
      base += 3;
      }

      // quad faces
      for (i = 2; i < 5; i++)
      {
      FE_TFaceTest<3,ZORDER> facetest;
      int ntest = facetest.GetNDof();
      Matrix<> facemoments(ntest, nd);

      ComputeFaceMoments (i, facetest, facemoments, 2*MAXORDER);
	
      for (j = 0; j < nd; j++)
      for (k = 0; k < ntest; k++)
      fiphij(base+k, j) = facemoments(k,j);
      base += ntest;
      }

      FE_TVolTest3<ZORDER> voltest;
      int ntest = voltest.GetNDof();
      Matrix<> volmoments(ntest, nd);
    
      ComputeVolMoments (voltest, volmoments, 2*MAXORDER);
      for (j = 0; j < nd; j++)
      for (k = 0; k < ntest; k++)
      fiphij(base+k, j) = volmoments(k,j);
      base += ntest;
    
      trans.SetSize (nd, nd);
      CalcInverse (fiphij, trans);
    */

    /*
      nd = NEDGEDOF;
      Matrix<> fiphij2(nd);
    
      // horizontal+vertical edges
      base = 0;
      for (i = 0; i < 9; i++)
      {
      int nedge = (i < 6) ? 3 : ZORDER;

      ComputeEdgeMoments (i, segm, edgemoments, 2*MAXORDER, 2);
      nedge--;
      for (j = 0; j < nd; j++)
      for (k = 0; k < nedge; k++)
      fiphij2(base+k, j) = edgemoments(k+1,j);

      base += nedge;
      }

      trans2.SetSize (nd, nd);
      CalcInverse (fiphij2, trans2);
    */


    nd = NQUADFACEDOF;
    Matrix<> fiphij3(nd);
    
    // quad faces
    base = 0;
    for (i = 2; i < 5; i++)
      {
	FE_TFaceTest<3,ZORDER> facetest;
	int ntest = facetest.GetNDof();
	Matrix<> facemoments(ntest, nd);

	ComputeFaceMoments (i, facetest, facemoments, 2*MAXORDER, 3);
	
	for (j = 0; j < nd; j++)
	  for (k = 0; k < ntest; k++)
	    fiphij3(base+k, j) = facemoments(k,j);
	base += ntest;
      }

    trans_quad.SetSize (nd, nd);
    CalcInverse (fiphij3, trans_quad);




    nd = NTRIGFACEDOF+NINNERDOF;
    Matrix<> fiphij4(nd);
    
    // trig faces
    base = 0;
    for (i = 0; i < 2; i++)
      {
	Matrix<> facemoments(3, nd);

	ComputeFaceMoments (i, rttrig0, facemoments, 4, 4);
	for (j = 0; j < nd; j++)
	  {
	    fiphij4(base  , j) = facemoments(1, j);
	    fiphij4(base+1, j) = facemoments(0, j);
	    fiphij4(base+2, j) = facemoments(2, j);
	  }
	base += 3;
      }
    
    {
      FE_TVolTest3NoGrad<ZORDER> voltest;
      int ntest = voltest.GetNDof();
      Matrix<> volmoments(ntest, nd);
    
      ComputeVolMoments (voltest, volmoments, 2*MAXORDER, 4);
      for (j = 0; j < nd; j++)
	for (k = 0; k < ntest; k++)
	  fiphij4(base+k, j) = volmoments(k,j);
      base += ntest;
    }

    // (*testout) << "fiphij4 = " << endl << fiphij4 << endl;
    trans_trig.SetSize (nd, nd);
    CalcInverse (fiphij4, trans_trig);
  }
  



  template class  FE_TNedelecPrism3NoGrad<1>;
  template class  FE_TNedelecPrism3NoGrad<2>;
  template class  FE_TNedelecPrism3NoGrad<3>;









  /* ************************ Pyramid elements ********************** */



  /*
  // Array<HCurlFiniteElement<3>::IPData> FE_NedelecPyramid1::ipdata;
  Matrix<> FE_NedelecPyramid1::trans(8);

  FE_NedelecPyramid1 :: FE_NedelecPyramid1()
    : HCurlFiniteElement<3> (8, 3)
  {
    // if (!ipdata.Size())
    // Orthogonalize();
      // CalcIPData(ipdata);
      // controldestruction.NewControl(&ipdatadestructed);
  }

  FE_NedelecPyramid1 :: ~FE_NedelecPyramid1()
  {
    // if(!ipdatadestructed)
    // ipdata.DeleteAll();
  }


  void FE_NedelecPyramid1 :: CalcShape (const IntegrationPoint & ip, 
					SliceMatrix<> shape) const
  {
    shape = 0.0; //!
    Mat<8,3> hshape;
    CalcShape1 (ip, hshape);
    shape = Trans (trans) * hshape;
  }


  void FE_NedelecPyramid1 :: CalcShape1 (const IntegrationPoint & ip, 
                                         FlatMatrixFixWidth<3> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);
    
    if (z == 1) z = 1-1e-8;

    double xr = x/(1-z);
    double yr = y/(1-z);

    Mat<8,3> hshape; 
    hshape = 0;
    // (1-z)^2 RT0 
    hshape(0,0) =      (1-z)*(1-z);
    hshape(1,0) = yr * (1-z)*(1-z);
    hshape(2,1) =      (1-z)*(1-z);
    hshape(3,1) = xr * (1-z)*(1-z);

    // \nabla (1-z) Q1
    hshape(4,2) = 1;
    hshape(5,0) = -(1-z);
    hshape(5,2) = xr;
    hshape(6,1) = -(1-z);
    hshape(6,2) = yr;
    hshape(7,0) = -yr * (1-z);
    hshape(7,1) = -xr * (1-z);
    hshape(7,2) = xr * yr;

    Mat<3,3> finv;
    finv = 0;
    finv(0,0) = 1/(1-z);
    finv(1,1) = 1/(1-z);
    finv(2,0) = xr/(1-z);
    finv(2,1) = yr/(1-z);
    finv(2,2) = 1;

    shape = hshape * Trans (finv);
  }



  void FE_NedelecPyramid1 :: Orthogonalize()
  {
    int nd = 8;

    Matrix<> fiphij(nd);
    fiphij = 0;

    Matrix<> edgemoments(2, nd);
    FE_Segm1L2 segm1;
    
    for (int i = 0; i < 8; i++)
      {
	ComputeEdgeMoments (i, segm1, edgemoments, 2);

	for (int j = 0; j < nd; j++)
	  fiphij(i, j) = edgemoments(0, j);
      }

    CalcInverse (fiphij, trans);
  }
  */

  










  // Array<HCurlFiniteElement<3>::IPData> FE_NedelecPyramid2::ipdata;
  Matrix<> FE_NedelecPyramid2::trans(NDOF);
  Matrix<> FE_NedelecPyramid2::trans2(NEDGEDOF);

  FE_NedelecPyramid2 :: FE_NedelecPyramid2()
    : HCurlFiniteElement<3> (NDOF, 2)
  {
    // if (!ipdata.Size())
      Orthogonalize();
    
      // CalcIPData(ipdata);
      // controldestruction.NewControl(&ipdatadestructed);
  }

  FE_NedelecPyramid2 :: ~FE_NedelecPyramid2()
  {
    // if(!ipdatadestructed)
    // ipdata.DeleteAll();
  }


  void FE_NedelecPyramid2 :: CalcShape (const IntegrationPoint & ip, 
					SliceMatrix<> shape) const
  {
    Mat<NDOF,3> hshape;
    Mat<8,3> shape1;

    Mat<NEDGEDOF,3> hshape2;
    Mat<NEDGEDOF,3> shape2;

    // Mat<4,3> hshape3;
    // Mat<4,3> shape3;

    shape = 0.0; //!
    CalcShape1 (ip, hshape);
    shape = Trans (trans) * hshape;

    CalcShape2 (ip, hshape2);
    shape2 = Trans (trans2) * hshape2;

    pyramid1.CalcShape (ip, shape1);
    int i, j;
    for (i = 0; i < 8; i++)
      for (j = 0; j < 3; j++)
	shape(i,j) = shape1(i,j);

    for (i = 0; i < 8; i++)
      for (j = 0; j < 3; j++)
	shape(i+8,j) = shape2(i,j);
  }


  void FE_NedelecPyramid2 :: CalcShape1 (const IntegrationPoint & ip, 
					 FlatMatrixFixWidth<3> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);

    if (z == 1) z = 1-1e-10;

    double xr = x/(1-z);
    double yr = y/(1-z);

    Mat<21,3> hshape;
    Vec<4> q1shape;
    Vec<9> q2shape;
    Mat<4,2> q1dshape;
    Mat<9,2> q2dshape;

    IntegrationPoint ipxy(xr, yr, 0, 0);
    quad1().CalcShape (ipxy, q1shape);
    quad2().CalcShape (ipxy, q2shape);

    quad1().CalcDShape (ipxy, q1dshape);
    quad2().CalcDShape (ipxy, q2dshape);

    hshape = 0;
    int ii = 0, i;

    // \nabla (1-z) Q1:  .... 4 dof
    for (i = 0; i < 4; i++)
      {
	hshape(ii, 0) = (1-z) * q1dshape(i,0);
	hshape(ii, 1) = (1-z) * q1dshape(i,1);
	hshape(ii, 2) = -q1shape(i);
	ii++;
      }

    // \nabla (1-z)^2 Q2:  .... 9 dof
    for (i = 0; i < 9; i++)
      {
	hshape(ii, 0) = (1-z) * (1-z) * q2dshape(i,0);
	hshape(ii, 1) = (1-z) * (1-z) * q2dshape(i,1);
	hshape(ii, 2) = 2 * (z-1) * q2shape(i);
	ii++;
      }

    // (1-z)^2 [N1 + N-El2-bubble] .... 8-1 dof
    double z2 = (1-z) * (1-z);
    hshape(ii++,0) = z2;
    hshape(ii++,0) = z2 * yr;
    hshape(ii++,1) = z2;
    hshape(ii++,1) = z2 * xr;
    hshape(ii++,0) = z2 * yr * (1-yr);
    hshape(ii++,1) = z2 * xr * (1-xr);

    hshape(ii  ,0) = z2 * yr * (1-yr) * xr;
    hshape(ii++,1) = -z2 * xr * (1-xr) * yr;

    Mat<3,3> finv;
    finv = 0;
    finv(0,0) = 1/(1-z);
    finv(1,1) = 1/(1-z);
    finv(2,0) = xr/(1-z);
    finv(2,1) = yr/(1-z);
    finv(2,2) = 1;

    // shape = ChangeSize(hshape * Trans (finv),shape.Height(),3);
    shape.Rows(0, hshape.Height()) = hshape * Trans (finv);
  }




  void FE_NedelecPyramid2 :: CalcShape2 (const IntegrationPoint & ip, 
					 FlatMatrixFixWidth<3> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);

    if (z == 1) z = 1-1e-8;

    double xr = x/(1-z);
    double yr = y/(1-z);

    Mat<8,3> hshape;
    Vec<4> q1shape;
    Vec<4> q2bshape;
    Mat<4,2> q1dshape;
    Mat<4,2> q2bdshape;


    q1shape(0) = 1;
    q1shape(1) = xr;
    q1shape(2) = yr;
    q1shape(3) = xr*yr;
    q1dshape = 0;
    q1dshape(1,0) = 1;
    q1dshape(2,1) = 1;
    q1dshape(3,0) = yr;
    q1dshape(3,1) = xr;

    q2bshape(0) = xr*(1-xr) * yr;
    q2bshape(1) = xr*(1-xr) * (1-yr);
    q2bshape(2) = xr        * (1-yr)*yr;
    q2bshape(3) = (1-xr)    * (1-yr)*yr;

    q2bdshape(0,0) = (1-2*xr) * yr;
    q2bdshape(0,1) = xr*(1-xr);
    q2bdshape(1,0) = (1-2*xr) * (1-yr);
    q2bdshape(1,1) = -xr*(1-xr);
    q2bdshape(2,0) = (1-yr)*yr;
    q2bdshape(2,1) = xr * (1-2*yr);
    q2bdshape(3,0) = -(1-yr)*yr;
    q2bdshape(3,1) = (1-xr) * (1-2*yr);


    hshape = 0;
    int ii = 0, i;

    // \nabla z(1-z) Q1:
    for (i = 0; i < 4; i++)
      {
	hshape(ii, 0) = z*(1-z) * q1dshape(i,0);
	hshape(ii, 1) = z*(1-z) * q1dshape(i,1);
	hshape(ii, 2) = (1-2*z) * q1shape(i);
	ii++;
      }

    // \nabla (1-z)^2 Q2b:
    for (i = 0; i < 4; i++)
      {
	hshape(ii, 0) = (1-z) * (1-z) * q2bdshape(i,0);
	hshape(ii, 1) = (1-z) * (1-z) * q2bdshape(i,1);
	hshape(ii, 2) = 2*(z-1) * q2bshape(i);
	ii++;
      }

    Mat<3,3> finv;
    finv = 0;
    finv(0,0) = 1/(1-z);
    finv(1,1) = 1/(1-z);
    finv(2,0) = xr/(1-z);
    finv(2,1) = yr/(1-z);
    finv(2,2) = 1;

    // shape = ChangeSize(hshape * Trans (finv),shape.Height(),3);
    shape.Rows(0, hshape.Height()) = hshape * Trans (finv);
  }






  void FE_NedelecPyramid2 :: Orthogonalize()
  {
    // cout << "compute Nedelec pyramid 2" << endl;

    Mat<NDOF> fiphij;
    fiphij = 0;

    Matrix<> edgemoments(2, NDOF);
    FE_Segm1L2 segm1;
    
    for (int i = 0; i < 8; i++)
      {
	ComputeEdgeMoments (i, segm1, edgemoments, 4);

	for (int j = 0; j < NDOF; j++)
	  {
	    fiphij(i, j) = edgemoments(0, j);
	    fiphij(i+8, j) = edgemoments(1, j);
	  }
      }

    Matrix<> facemoments(4,NDOF);
    FE_RTQuad0 rtquad0;
    ComputeFaceMoments (4, rtquad0, facemoments, 4);

    for (int j = 0; j < NDOF; j++)
      {
	fiphij(16, j) = facemoments(0, j);
	fiphij(17, j) = facemoments(1, j);
	fiphij(18, j) = facemoments(2, j);
	fiphij(19, j) = facemoments(3, j);
      }

    CalcInverse (fiphij, trans);


    // edge dofs:

    Mat<NEDGEDOF> fiphij2;
    fiphij2 = 0;

    for (int i = 0; i < NEDGEDOF; i++)
      {
	ComputeEdgeMoments (i, segm1, edgemoments, 4, 2);

	for (int j = 0; j < NEDGEDOF; j++)
	  fiphij2(i, j) = edgemoments(1, j);
      }

    CalcInverse (fiphij2, trans2);
  }









  class FE_Pyramid3RefEdgeBubble : public ScalarFiniteElement<3>
  {
  public:
    ///
    FE_Pyramid3RefEdgeBubble()
      : ScalarFiniteElement<3> (16, 3) { ; } 

    virtual ELEMENT_TYPE ElementType() const override { return ET_PYRAMID; }
    
    ///
    using  ScalarFiniteElement<3>::CalcShape;
    virtual void CalcShape (const IntegrationPoint & ip, 
			    BareSliceVector<> shape) const override
    {
      double x = ip(0);
      double y = ip(1);
      double z = ip(2);
    
      shape.Range(0,ndof) = 0.0; //!
      int ii = 0;

      double fac = z * (1-z);
      shape(ii++) = fac;
      shape(ii++) = fac * x;
      shape(ii++) = fac * y;
      shape(ii++) = fac * x*y;

      // z(1-z)^2 * computed from 2D edge bubbles
      fac = z * (1-z) * (1-z);
      shape(ii++) = fac * (y-0)*(y+1) * (x-0)*(x+1);
      shape(ii++) = fac * (y-0)*(y+1) * (x-1)*(x-2);
      shape(ii++) = fac * (y-1)*(y-2) * (x-0)*(x+1);
      shape(ii++) = fac * (y-1)*(y-2) * (x-1)*(x-2);

      // (1-z)^2 * Q2-edge bubbles:
      fac = (1-z) * (1-z);
      shape(ii++) = fac * y * (1-y) * x;
      shape(ii++) = fac * y * (1-y) * (1-x);
      shape(ii++) = fac * x * (1-x) * y;
      shape(ii++) = fac * x * (1-x) * (1-y);

      
      fac = (1-z) * (1-z) * (1-z);
      double phi3x = (1-2*x) * x*(1-x);   // = (2 x^3 - 3 x^2 + x ) * yr
      double phi3y = (1-2*y) * y*(1-y);

      shape(ii++) = fac * phi3x * y;  
      shape(ii++) = fac * phi3x * (1-y);
      shape(ii++) = fac * x        * phi3y;
      shape(ii++) = fac * (1-x)    * phi3y;
    }


    virtual void CalcDShape (const IntegrationPoint & ip, 
			     BareSliceMatrix<> dshape) const override
    {
      double x = ip(0);
      double y = ip(1);
      double z = ip(2);
    
      int ii = 0;
      dshape.AddSize(GetNDof(), 3) = 0;
      double fac = z * (1-z);
      double dfac = 1 - 2*z;
      /*
        shape(ii++) = fac;
        shape(ii++) = fac * x;
        shape(ii++) = fac * y;
        shape(ii++) = fac * x*y;
      */
      dshape(ii,2) = dfac;
      ii++;
      dshape(ii,0) = fac;
      dshape(ii,2) = dfac * x;
      ii++;
      dshape(ii,1) = fac;
      dshape(ii,2) = dfac * y;
      ii++;
      dshape(ii,0) = fac*y;
      dshape(ii,1) = fac*x;
      dshape(ii,2) = dfac * x*y;
      ii++;


      // z(1-z)^2 * computed from 2D edge bubbles
      fac = z * (1-z) * (1-z);  // z^3 - 2*z^2 + z
      dfac = 3*z*z - 4 * z+1;
      /*
        shape(ii++) = fac * (y-0)*(y+1) * (x-0)*(x+1);
        shape(ii++) = fac * (y-0)*(y+1) * (x-1)*(x-2);
        shape(ii++) = fac * (y-1)*(y-2) * (x-0)*(x+1);
        shape(ii++) = fac * (y-1)*(y-2) * (x-1)*(x-2);
      */
      dshape(ii,0) = fac * (y-0)*(y+1) * (2*x+1);
      dshape(ii,1) = fac * (2*y+1) * (x-0)*(x+1);
      dshape(ii,2) = dfac * (y-0)*(y+1) * (x-0)*(x+1);
      ii++;
      dshape(ii,0) = fac * (y-0)*(y+1) * (2*x-3);
      dshape(ii,1) = fac * (2*y+1) * (x-1)*(x-2);
      dshape(ii,2) = dfac * (y-0)*(y+1) * (x-1)*(x-2);
      ii++;
      dshape(ii,0) = fac * (y-1)*(y-2) * (2*x+1);
      dshape(ii,1) = fac * (2*y-3) * (x-0)*(x+1);
      dshape(ii,2) = dfac * (y-1)*(y-2) * (x-0)*(x+1);
      ii++;
      dshape(ii,0) = fac * (y-1)*(y-2) * (2*x-3);
      dshape(ii,1) = fac * (2*y-3) * (x-1)*(x-2);
      dshape(ii,2) = dfac * (y-1)*(y-2) * (x-1)*(x-2);
      ii++;



      // (1-z)^2 * Q2-edge bubbles:
      fac = (1-z) * (1-z);
      dfac = 2*z-2;
      /*
        shape(ii++) = fac * y * (1-y) * x;
        shape(ii++) = fac * y * (1-y) * (1-x);
        shape(ii++) = fac * x * (1-x) * y;
        shape(ii++) = fac * x * (1-x) * (1-y);
      */
      dshape(ii,0) = fac * y * (1-y);
      dshape(ii,1) = fac * (1-2*y) * x;
      dshape(ii,2) = dfac * y * (1-y) * x;
      ii++;

      dshape(ii,0) = -fac * y * (1-y);
      dshape(ii,1) = fac * (1-2*y) * (1-x);
      dshape(ii,2) = dfac * y * (1-y) * (1-x);
      ii++;

      dshape(ii,0) = fac * (1-2*x) * y;
      dshape(ii,1) = fac * x * (1-x);
      dshape(ii,2) = dfac * x * (1-x) * y;
      ii++;

      dshape(ii,0) = fac * (1-2*x) * (1-y);
      dshape(ii,1) = -fac * x * (1-x);
      dshape(ii,2) = dfac * x * (1-x) * (1-y);
      ii++;

      fac = (1-z) * (1-z) * (1-z);
      dfac = -3 * (1-z) * (1-z);
      double phi3x = (1-2*x) * x*(1-x);   // = (2 x^3 - 3 x^2 + x ) 
      double phi3y = (1-2*y) * y*(1-y);
      double dphi3x = 6 * x*x - 6*x + 1;
      double dphi3y = 6 * y*y - 6*y + 1;
      /*
        shape(ii++) = fac * phi3x * y;  
        shape(ii++) = fac * phi3x * (1-y);
        shape(ii++) = fac * x        * phi3y;
        shape(ii++) = fac * (1-x)    * phi3y;
      */
      dshape(ii,0) = fac * dphi3x * y;
      dshape(ii,1) = fac * phi3x;
      dshape(ii,2) = dfac * phi3x * y;
      ii++;

      dshape(ii,0) = fac * dphi3x * (1-y);
      dshape(ii,1) = -fac * phi3x;
      dshape(ii,2) = dfac * phi3x * (1-y);
      ii++;

      dshape(ii,0) = fac * phi3y;
      dshape(ii,1) = fac * dphi3y * x;
      dshape(ii,2) = dfac * phi3y * x;
      ii++;

      dshape(ii,0) = -fac * phi3y;
      dshape(ii,1) = fac * dphi3y * (1-x);
      dshape(ii,2) = dfac * phi3y * (1-x);
      ii++;
    }

    ///
    // virtual const Array<IPData> & GetIPData () const { return ipdata; }
  }; 



  class FE_Pyramid3RefFaceBubble : public ScalarFiniteElement<3>
  {
    // Array<IPData> ipdata;
    
  public:
    ///
    FE_Pyramid3RefFaceBubble()
      : ScalarFiniteElement<3> (4, 3) 
    { ; } 
    ///
    virtual ~FE_Pyramid3RefFaceBubble() 
    { 
      // ipdata.DeleteAll(); 
    }
    ///

    virtual ELEMENT_TYPE ElementType() const override { return ET_PYRAMID; }

    using ScalarFiniteElement<3>::CalcShape;
    virtual void CalcShape (const IntegrationPoint & ip, 
			    BareSliceVector<> shape) const override
    {
      double x = ip(0);
      double y = ip(1);
      double z = ip(2);
    
      shape.Range(0,ndof) = 0.0;
      int ii = 0;

      double fac = z * (1-z) * (1-z);
      shape(ii++) = fac * x * (1-x) * y;
      shape(ii++) = fac * x * (1-x) * (1-y);
      shape(ii++) = fac * y * (1-y) * x;
      shape(ii++) = fac * y * (1-y) * (1-x);
    }

    virtual void CalcDShape (const IntegrationPoint & ip, 
			     BareSliceMatrix<> dshape) const override
    {
      cerr << "shape not implemented" << endl;
    }
  };




  FE_Quad3 :: FE_Quad3()
    : ScalarFiniteElement<2> (16, 3)
  {
    ;
  }

  FE_Quad3 :: ~FE_Quad3()
  {
    ;
  }

  void FE_Quad3 :: CalcShape (const IntegrationPoint & ip, 
			      BareSliceVector<> shape) const
    
  {
    double x = ip(0);
    double y = ip(1);

    Vec<4> px, py;
    px(0) = 1-x;
    px(1) = x;
    px(2) = x * (1-x);
    px(3) = px(2) * (1-2*x);

    py(0) = 1-y;
    py(1) = y;
    py(2) = y * (1-y);
    py(3) = py(2) * (1-2*y);

    int ii = 0;
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
	shape(ii++) = px(i) * py(j);
  }



  void FE_Quad3 :: CalcDShape (const IntegrationPoint & ip, 
			       BareSliceMatrix<> dshape) const
    
  {
    double x = ip(0);
    double y = ip(1);

    Vec<4> px, py, dpx, dpy;
    px(0) = 1-x;
    px(1) = x;
    px(2) = x * (1-x);
    px(3) = px(2) * (1-2*x);

    py(0) = 1-y;
    py(1) = y;
    py(2) = y * (1-y);
    py(3) = py(2) * (1-2*y);

    dpx(0) = -1;
    dpx(1) = 1;
    dpx(2) = 1-2*x; 
    dpx(3) = 6*x*x-6*x+1;

    dpy(0) = -1;
    dpy(1) = 1;
    dpy(2) = 1-2*x; 
    dpy(3) = 6*x*x-6*x+1;


    int ii = 0;
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
	{
	  dshape(ii,0) = dpx(i) * py(j);
	  dshape(ii,1) = px(i) * dpy(j);
	  ii++;
	}
  }


#ifdef VERY_OLD_NEDELEC_ELS
  // Array<HCurlFiniteElement<3>::IPData> FE_NedelecPyramid3::ipdata;
  Mat<FE_NedelecPyramid3::NDOF> FE_NedelecPyramid3::trans;
  Mat<FE_NedelecPyramid3::NEDGEDOF> FE_NedelecPyramid3::trans2;
  Mat<FE_NedelecPyramid3::NFACEDOF> FE_NedelecPyramid3::trans3;


  FE_NedelecPyramid3 :: FE_NedelecPyramid3()
    : HCurlFiniteElement<3> (NDOF, 3)
  {
    // if (!ipdata.Size())
      Orthogonalize();

      // CalcIPData(ipdata);
      // controldestruction.NewControl(&ipdatadestructed);
  }

  FE_NedelecPyramid3 :: ~FE_NedelecPyramid3()
  {
    // if(!ipdatadestructed)
    // ipdata.DeleteAll();
  }


  void FE_NedelecPyramid3 :: CalcShape (const IntegrationPoint & ip, 
					SliceMatrix<> shape) const
  {
	  cout << "old nedelec pyramid disabled" << endl;

	  shape = 0.0; //!
	  /*
	  Mat<NDOF,3> hshape;
    Mat<8,3> shape1;
    CalcShape1 (ip, hshape);
    shape = Trans (trans) * hshape;

    Mat<NEDGEDOF,3> hshape2;
    Mat<NEDGEDOF,3> shape2;
    CalcShape2 (ip, hshape2);
    shape2 = Trans (trans2) * hshape2;

    Mat<NFACEDOF,3> hshape3;
    Mat<NFACEDOF,3> shape3;
    CalcShape3 (ip, hshape3);
    shape3 = Trans (trans3) * hshape3;

    pyramid1.CalcShape (ip, shape1);
    int i, j;
    for (i = 0; i < 8; i++)
      for (j = 0; j < 3; j++)
	shape(i,j) = shape1(i,j);

    for (i = 0; i < NEDGEDOF; i++)
      for (j = 0; j < 3; j++)
	shape(i+8,j) = shape2(i,j);

    for (i = 0; i < NFACEDOF; i++)
      for (j = 0; j < 3; j++)
	shape(i+8+NEDGEDOF,j) = shape3(i,j);
	*/
  }


  void FE_NedelecPyramid3 :: CalcShape1 (const IntegrationPoint & ip, 
					 FlatMatrixFixWidth<3> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);

    if (z == 1) z = 1-1e-8;

    double xr = x/(1-z);
    double yr = y/(1-z);

    Mat<NDOF,3> hshape;
    Vec<4> q1shape;
    Vec<9> q2shape;
    Vec<16> q3shape;
    Mat<4,2> q1dshape;
    Mat<9,2> q2dshape;
    Mat<16,2> q3dshape;

    IntegrationPoint ipxy(xr, yr, 0, 0);

    quad1().CalcShape (ipxy, q1shape);
    quad2().CalcShape (ipxy, q2shape);
    quad3.CalcShape (ipxy, q3shape);

    quad1().CalcDShape (ipxy, q1dshape);
    quad2().CalcDShape (ipxy, q2dshape);
    quad3.CalcDShape (ipxy, q3dshape);

    hshape = 0;
    int ii = 0;

    // 29 grad dofs:
    // \nabla (1-z) Q1:  .... 4 dof
    for (int i = 0; i < 4; i++)
      {
	hshape(ii, 0) = (1-z) * q1dshape(i,0);
	hshape(ii, 1) = (1-z) * q1dshape(i,1);
	hshape(ii, 2) = -q1shape(i);
	ii++;
      }

    // \nabla (1-z)^2 Q2:  .... 9 dof
    for (int i = 0; i < 9; i++)
      {
	hshape(ii, 0) = (1-z) * (1-z) * q2dshape(i,0);
	hshape(ii, 1) = (1-z) * (1-z) * q2dshape(i,1);
	hshape(ii, 2) = 2 * (z-1) * q2shape(i);
	ii++;
      }

    // \nabla (1-z)^3 Q3:  .... 16 dof
    for (int i = 0; i < 16; i++)
      {
	hshape(ii, 0) = (1-z) * (1-z) * (1-z) * q3dshape(i,0);
	hshape(ii, 1) = (1-z) * (1-z) * (1-z) * q3dshape(i,1);
	hshape(ii, 2) = -3 * (z-1) * (z-1) * q3shape(i);
	ii++;
      }



    // 8 dof
    double z2 = (1-z) * (1-z);
    hshape(ii++,0) = z2;
    hshape(ii++,0) = z2 * yr;
    hshape(ii++,1) = z2;
    hshape(ii++,1) = z2 * xr;
    hshape(ii++,0) = z2 * yr * (1-yr);
    hshape(ii++,0) = z2 * yr * (1-yr) * xr;
    hshape(ii++,1) = z2 * xr * (1-xr);
    hshape(ii++,1) = z2 * xr * (1-xr) * yr;

    // 20 dof
    double z3 = (1-z) * (1-z) * (1-z);
    hshape(ii++,0) = z3;
    hshape(ii++,0) = z3 * xr;
    hshape(ii++,0) = z3 * yr;
    hshape(ii++,0) = z3 * xr * yr;
    hshape(ii++,1) = z3;
    hshape(ii++,1) = z3 * xr;
    hshape(ii++,1) = z3 * yr;
    hshape(ii++,1) = z3 * xr * yr;

    hshape(ii++,0) = z3 * yr * (1-yr);
    hshape(ii++,0) = z3 * yr * (1-yr) * xr;
    hshape(ii++,0) = z3 * yr * (1-yr) * xr * xr;
    hshape(ii++,0) = z3 * yr * (1-yr) * yr; 
    hshape(ii++,0) = z3 * yr * (1-yr) * yr * xr;
    hshape(ii++,0) = z3 * yr * (1-yr) * yr * xr * xr;

    hshape(ii++,1) = z3 * xr * (1-xr);
    hshape(ii++,1) = z3 * xr * (1-xr) * yr;
    hshape(ii++,1) = z3 * xr * (1-xr) * yr * yr;
    hshape(ii++,1) = z3 * xr * (1-xr) * xr;
    hshape(ii++,1) = z3 * xr * (1-xr) * xr * yr;
    hshape(ii++,1) = z3 * xr * (1-xr) * xr * yr * yr;
  
    Mat<3,3> finv;
    finv = 0;
    finv(0,0) = 1/(1-z);
    finv(1,1) = 1/(1-z);
    finv(2,0) = xr/(1-z);
    finv(2,1) = yr/(1-z);
    finv(2,2) = 1;

    shape = hshape * Trans (finv);
  }




  void FE_NedelecPyramid3 :: CalcShape2 (const IntegrationPoint & ip, 
					 FlatMatrixFixWidth<3> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);

    if (z == 1) z = 1-1e-8;

    double xr = x/(1-z);
    double yr = y/(1-z);

    IntegrationPoint ipr(xr, yr, z, 0);
    Mat<16,3> hshape;

    FE_Pyramid3RefEdgeBubble febub;
    febub.CalcDShape (ipr, FlatMatrix<> (hshape));

    Mat<3,3> finv;
    finv = 0;
    finv(0,0) = 1/(1-z);
    finv(1,1) = 1/(1-z);
    finv(2,0) = xr/(1-z);
    finv(2,1) = yr/(1-z);
    finv(2,2) = 1;

    shape = hshape * Trans (finv);
  }





  void FE_NedelecPyramid3 :: CalcShape3 (const IntegrationPoint & ip, 
					 FlatMatrixFixWidth<3> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);

    if (z == 1) z = 1-1e-8;

    double xr = x/(1-z);
    double yr = y/(1-z);

    IntegrationPoint ipr(xr, yr, z, 0);
    Mat<NFACEDOF,3> hshape;
    hshape = 0;

    int i, j, ii = 0;

    // 12 quad face shape functions:

    // gradient of Q3 bubbles:
    hshape(ii, 0) = (1-z) * (1-z) * (1-z) * (1-2*xr)*yr*(1-yr);
    hshape(ii, 1) = (1-z) * (1-z) * (1-z) * xr*(1-xr)*(1-2*yr);
    hshape(ii, 2) = -3 * (z-1) * (z-1) * xr*(1-xr)*yr*(1-yr);
    ii++;

    hshape(ii, 0) = (1-z) * (1-z) * (1-z) * (2*xr-3*xr*xr)*yr*(1-yr);
    hshape(ii, 1) = (1-z) * (1-z) * (1-z) * xr*xr*(1-xr)*(1-2*yr);
    hshape(ii, 2) = -3 * (z-1) * (z-1) * xr*xr*(1-xr)*yr*(1-yr);
    ii++;

    hshape(ii, 0) = (1-z) * (1-z) * (1-z) * (1-2*xr)*yr*yr*(1-yr);
    hshape(ii, 1) = (1-z) * (1-z) * (1-z) * xr*(1-xr)*(2*yr-3*yr*yr);
    hshape(ii, 2) = -3 * (z-1) * (z-1) * xr*(1-xr)*yr*yr*(1-yr);
    ii++;

    hshape(ii, 0) = (1-z) * (1-z) * (1-z) * (2*xr-3*xr*xr)*yr*yr*(1-yr);
    hshape(ii, 1) = (1-z) * (1-z) * (1-z) * xr*xr*(1-xr)*(2*yr-3*yr*yr);
    hshape(ii, 2) = -3 * (z-1) * (z-1) * xr*xr*(1-xr)*yr*yr*(1-yr);
    ii++;

    
    double z3 = (1-z) * (1-z) * (1-z);

    hshape(ii++,0) = z3 * yr * (1-yr);
    // hshape(ii++,0) = z3 * yr * (1-yr) * xr;
    hshape(ii++,0) = z3 * yr * (1-yr) * xr * xr;
    hshape(ii++,0) = z3 * yr * (1-yr) * (1-2*yr); 
    // hshape(ii++,0) = z3 * yr * (1-yr) * (1-2*yr) * xr;
    hshape(ii++,0) = z3 * yr * (1-yr) * (1-2*yr) * xr * xr;

    hshape(ii++,1) = z3 * xr * (1-xr);
    hshape(ii++,1) = z3 * xr * (1-xr) * yr;
    hshape(ii++,1) = z3 * xr * (1-xr) * yr * yr;
    hshape(ii++,1) = z3 * xr * (1-xr) * (1-2*xr);
    // hshape(ii++,1) = z3 * xr * (1-xr) * (1-2*xr) * yr;
    // hshape(ii++,1) = z3 * xr * (1-xr) * (1-2*xr) * yr * yr;


    
    // triangular face bubbles:

    FE_Pyramid3RefFaceBubble facebub;    
    Mat<4,3> gradients;

    facebub.CalcDShape (ipr, gradients);

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
	hshape(ii+i, j) = gradients(i,j);
    ii+= 4;


    double fac = z * (1-z);
    hshape(ii,0) = fac * yr * xr * (1-z);
    hshape(ii,2) = fac * yr * xr * (1-xr);
    ii++;
    hshape(ii,0) = fac * (1-yr) * xr * (1-z);
    hshape(ii,2) = fac * (1-yr) * xr * (1-xr);
    ii++;
    hshape(ii,0) = -fac * yr * (1-xr) * (1-z);
    hshape(ii,2) = fac * yr * xr * (1-xr);
    ii++;
    hshape(ii,0) = -fac * (1-yr) * (1-xr) * (1-z);
    hshape(ii,2) = fac * (1-yr) * xr * (1-xr);
    ii++;
    
    hshape(ii,1) = fac * xr * yr * (1-z);
    hshape(ii,2) = fac * xr * yr * (1-yr);
    ii++;
    hshape(ii,1) = fac * (1-xr) * yr * (1-z);
    hshape(ii,2) = fac * (1-xr) * yr * (1-yr);
    ii++;
    hshape(ii,1) = -fac * xr * (1-yr) * (1-z);
    hshape(ii,2) = fac * xr * yr * (1-yr);
    ii++;
    hshape(ii,1) = -fac * (1-xr) * (1-yr) * (1-z);
    hshape(ii,2) = fac * (1-xr) * yr * (1-yr);
    ii++;
    

    Mat<3,3> finv;
    finv = 0;
    finv(0,0) = 1/(1-z);
    finv(1,1) = 1/(1-z);
    finv(2,0) = xr/(1-z);
    finv(2,1) = yr/(1-z);
    finv(2,2) = 1;

    shape = hshape * Trans (finv);
  }




  void FE_NedelecPyramid3 :: GetInternalDofs (Array<int> & idofs) const
  {
    idofs.SetSize (0);
    for (int i = NDOF - NINNERDOF; i < NDOF; i++)
      idofs.Append (i);
  }


  void FE_NedelecPyramid3 :: Orthogonalize()
  {
    // cout << "compute Nedelec pyramid 3" << endl;

    Matrix<> fiphij(NDOF);
    fiphij = 0;

    Matrix<> edgemoments(3, NDOF);
    FE_Segm2L2 segm2;
    
    for (int i = 0; i < 8; i++)
      {
	ComputeEdgeMoments (i, segm2, edgemoments, 5);

	for (int j = 0; j < NDOF; j++)
	  {
	    fiphij(i, j) = edgemoments(0, j);
	    fiphij(i+8, j) = edgemoments(1, j);
	    fiphij(i+16, j) = edgemoments(2, j);
	  }
      }

    int ii = 24;

    Matrix<> facemoments3(3,NDOF);
    FE_RTTrig0 rttrig0;

    // 4*3 = 12
    for (int i = 0; i < 4; i++)
      {
	ComputeFaceMoments (i, rttrig0, facemoments3, 4);
	for (int j = 0; j < NDOF; j++)
	  {
	    fiphij(ii, j)   = facemoments3(1, j);
	    fiphij(ii+1, j) = facemoments3(0, j);
	    fiphij(ii+2, j) = facemoments3(2, j);
	  }
	ii+=3;
      }

    FE_TFaceTest<3,3> quadtest;
    int nqtest = quadtest.GetNDof();   // 12
    Matrix<> facemoments4(nqtest, NDOF);

    ComputeFaceMoments (4, quadtest, facemoments4, 6);
	
    for (int j = 0; j < NDOF; j++)
      for (int k = 0; k < nqtest; k++)
	fiphij(ii+k, j) = facemoments4(k,j);
    ii += nqtest;

    Matrix<> f2(ii, NDOF);
    
    for (int i = 0; i < ii; i++)
      for (int j = 0; j < NDOF; j++)
	f2(i,j) = fiphij(i,j);

    // trans.SetSize (NDOF, NDOF);
    Householder (f2, trans);

    /*
      (*testout) << "pyramid3, fiphij = " << endl << fiphij << endl;
      (*testout) << "trans = " << endl << trans << endl;
      (*testout) << "check = " << endl << (fiphij * trans) << endl;
    */

    Mat<NEDGEDOF> fiphij2;
    fiphij2 = 0;

    for (int i = 0; i < 8; i++)
      {
	ComputeEdgeMoments (i, segm2, edgemoments, 4, 2);

	for (int j = 0; j < NEDGEDOF; j++)
	  {
	    fiphij2(i, j) = edgemoments(1, j);
	    fiphij2(i+8, j) = edgemoments(2, j);
	  }
      }
    
    // trans2.SetSize (NEDGEDOF);
    CalcInverse (fiphij2, trans2);

    /*
      (*testout) << "pyramid3, fiphij2 = " << endl << fiphij2 << endl;
      (*testout) << "trans2 = " << endl << trans2 << endl;
      (*testout) << "check = " << endl << (fiphij2 * trans2) << endl;
    */




    Matrix<> fiphij3(NFACEDOF);
    fiphij3 = 0;

    // 4*3 = 12
    ii = 0;
    for (int i = 0; i < 4; i++)
      {
	ComputeFaceMoments (i, rttrig0, facemoments3, 4, 3);
	for (int j = 0; j < NFACEDOF; j++)
	  {
	    fiphij3(ii, j)   = facemoments3(1, j);
	    fiphij3(ii+1, j) = facemoments3(0, j);
	    fiphij3(ii+2, j) = facemoments3(2, j);
	  }
	ii+=3;
      }

    ComputeFaceMoments (4, quadtest, facemoments4, 6, 3);
	
    for (int j = 0; j < NFACEDOF; j++)
      for (int k = 0; k < nqtest; k++)
	fiphij3(ii+k, j) = facemoments4(k,j);

    // trans2.SetSize (NEDGEDOF);
    CalcInverse (fiphij3, FlatMatrix<> (trans3));

    /*
      (*testout) << "pyramid3, fiphij3 = " << endl << fiphij3 << endl;
      (*testout) << "trans3 = " << endl << trans3 << endl;
      (*testout) << "check = " << endl << (fiphij3 * trans3) << endl;
    */
  }
#endif
  

}
