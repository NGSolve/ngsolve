/*********************************************************************/
/* File:   hdivfe.cc                                                 */
/* Author: Joachim Schoeberl                                         */
/* Date:   5. Jul. 2001                                              */
/*********************************************************************/

/* 
   Raviart Thomas finite elements
*/




#include <fem.hpp>
#include "h1lofe.hpp"
#include "hdivlofe.hpp"

namespace ngfem
{

  template <int D>
  string  HDivFiniteElement<D> :: ClassName() const
  {
    return "HDivFiniteElement"; 
  }
  
  template <int D>
  void HDivFiniteElement<D> ::
  CalcDivShape (const IntegrationPoint & ip, 
		SliceVector<> divshape) const
  {
    double eps = 1e-5;
    ArrayMem<double, 200> hm1(DIM*ndof), hm2(DIM*ndof), 
      hm3(DIM*ndof), hm4(DIM*ndof), hmi(DIM*ndof);

    FlatMatrixFixWidth<DIM> shape1(ndof, &hm1[0]);
    FlatMatrixFixWidth<DIM> shape2(ndof, &hm2[0]);
    FlatMatrixFixWidth<DIM> shape3(ndof, &hm3[0]);
    FlatMatrixFixWidth<DIM> shape4(ndof, &hm4[0]);
    
    FlatMatrix<> dshapei(ndof,D, &hmi[0]);
    divshape = 0;

    for (int i = 0; i < D; i++)
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

	 
	for (int j = 0; j < ndof; j++)
	  divshape(j) += dshapei(j,i);
      }
  }


  


  template <int D>
  void HDivFiniteElement<D> ::
  CalcNormalShape (const IntegrationPoint & ip, 
                   SliceVector<> nshape) const
  {
    Array<int> dnums;
    int fnr = ip.FacetNr();
    if (fnr < 0) cerr << "HDivFE::CalcNormalShape: not a facet ip" << endl;
    GetFacetDofs (fnr, dnums);

    NORMAL * normals = ElementTopology::GetNormals(ElementType());
    Vec<D> normal_ref;
    for (int i = 0; i < D; i++)
      normal_ref(i) = normals[fnr][i];

    MatrixFixWidth<D> shape(GetNDof());
    CalcShape (ip, shape);
    for (int i = 0; i < dnums.Size(); i++)
      nshape(i) = InnerProduct (normal_ref, shape.Row(dnums[i]));
  }



  /// compute shape
  template <int D>
  void HDivFiniteElement<D> ::
  CalcMappedShape (const BaseMappedIntegrationPoint & bmip,
                                  SliceMatrix<> shape) const
  {
    auto mip = static_cast<const MappedIntegrationPoint<D,D>&> (bmip);
    CalcShape (mip.IP(), shape);
    Mat<DIM> trans = (1.0/mip.GetJacobiDet()) * mip.GetJacobian();
    for (int i = 0; i < ndof; i++)
      {
        Vec<DIM> hs = shape.Row(i);
        shape.Row(i) = trans * hs;
      }
  }
  

  template <int D>
  void HDivFiniteElement<D> ::
  CalcMappedShape (const BaseMappedIntegrationRule & bmir, SliceMatrix<> shapes) const
  {
    auto mir = static_cast<const MappedIntegrationRule<D,D>&> (bmir);
    for (int i = 0; i < mir.Size(); i++)
      CalcMappedShape (mir[i], shapes.Cols(i*D, (i+1)*D));
  }

  template <int D>
  void HDivFiniteElement<D> ::
  CalcMappedShape (const SIMD<MappedIntegrationPoint<DIM,DIM>> & mip,
                   BareSliceMatrix<SIMD<double>> shape) const
  {
    throw ExceptionNOSIMD(string("SIMD - HDivFE::CalcMappedShape not overloaded, et = ")
                          + typeid(*this).name());
  }
  
  
  template <int D>
  void HDivFiniteElement<D> ::
  CalcMappedShape (const SIMD_BaseMappedIntegrationRule & mir, 
                   BareSliceMatrix<SIMD<double>> shapes) const
  {
    throw ExceptionNOSIMD(string("SIMD - HDivFE::CalcMappedShape not overloaded, et = ")
                          + typeid(*this).name());
  }

  template <int D>
  void HDivFiniteElement<D> ::
  CalcMappedNormalShape (const SIMD_BaseMappedIntegrationRule & mir, 
                         BareSliceMatrix<SIMD<double>> shapes) const
  {
    throw ExceptionNOSIMD(string("SIMD - HDivFE::CalcMappedShape not overloaded, et = ")
                          + typeid(*this).name());
  }
  
  
  /// compute curl of shape
  template <int D>
  void HDivFiniteElement<D> ::
  CalcMappedDivShape (const MappedIntegrationPoint<DIM,DIM> & mip,
                      SliceVector<> divshape) const
  {
    CalcDivShape (mip.IP(), divshape);
    divshape /= mip.GetJacobiDet();
  }

  template <int D>
  void HDivFiniteElement<D> ::
  CalcMappedDivShape (const SIMD_BaseMappedIntegrationRule & mir, 
                      BareSliceMatrix<SIMD<double>> shapes) const
  {
    throw ExceptionNOSIMD("SIMD - HDivFE::CalcMappedDivShape not overloaded");
  }


  template <int D>
  void HDivFiniteElement<D> ::
  Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, 
            FlatMatrixFixWidth<D> vals) const
  {
    MatrixFixWidth<D> shape(ndof);
    for (int i = 0; i < ir.GetNIP(); i++)
      {
        CalcShape (ir[i], shape);
        vals.Row(i) = Trans(shape) * coefs;
      }
  }
  
  template <int D>
  void HDivFiniteElement<D> ::
  EvaluateTrans (const IntegrationRule & ir, 
                 FlatMatrixFixWidth<D> vals,
                 FlatVector<double> coefs) const
  {
    MatrixFixWidth<D> shape(ndof);
    coefs = 0;
    for (int i = 0; i < ir.GetNIP(); i++)
      {
        CalcShape (ir[i], shape);
        coefs += shape * vals.Row(i);
      }
  }
  

  template <int D>
  void HDivFiniteElement<D> ::
  GetFacetDofs(int i, Array<int> & dnums) const
  { 
    cout  << " GetFacetDofs for nothing " << endl; 
    dnums.SetSize(0);
  }; 

  template <int D>
  void HDivFiniteElement<D> ::
  Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const
  {
    throw ExceptionNOSIMD ("HDivFE::Evaluate (simd) not overloaded");
  }
  
  template <int D>
  void HDivFiniteElement<D> ::
  AddTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values,
            BareSliceVector<> coefs) const
  {
    throw ExceptionNOSIMD ("HDivFE::AddTrans (simd) not overloaded");    
  }
  
  template <int D>
  void HDivFiniteElement<D> ::
  EvaluateDiv (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs, BareVector<SIMD<double>> values) const
  {
    throw ExceptionNOSIMD ("HDivFE::EvaluateDiv (simd) not overloaded");    
  }
    
  template <int D>
  void HDivFiniteElement<D> ::
  AddDivTrans (const SIMD_BaseMappedIntegrationRule & ir, BareVector<SIMD<double>> values,
               BareSliceVector<> coefs) const
  {
    throw ExceptionNOSIMD ("HDivFE::EvaluateDivTrans (simd) not overloaded");        
  }
      




  
  template <int D>
  void HDivFiniteElement<D> ::
  ComputeFaceMoments (int fnr, ScalarFiniteElement<DIM-1> & testfe,
		      FlatMatrix<> & moments, int order, int shapenr) const
  {
    int j;
    
    int test_ndof = testfe.GetNDof();
    int dim = D;
    
    MatrixFixWidth<DIM> shape(ndof);
    Matrix<> shapen(ndof, 1);
    Vector<> testshape(test_ndof);

   moments = 0;

   const IntegrationRule & facerule = SelectIntegrationRule (testfe.ElementType(), order);
    
   if (dim == 2)
     {
       const POINT3D * points = ElementTopology::GetVertices (ElementType());
       const EDGE & edge = ElementTopology::GetEdges (ElementType()) [fnr];
       
       Vec<2> p1, p2, tau, nv;

       for (j = 0; j < 2; j++)
	 {
	   p1(j) = points[edge[0]][j];
	   p2(j) = points[edge[1]][j];
	   tau(j) = p2(j) - p1(j);
	 }
       nv(0) = -tau(1); 
       nv(1) = tau(0);

       for (j = 0; j < facerule.GetNIP(); j++)
	 {
	   const IntegrationPoint & ip = facerule[j];
	   Vec<1> p1d;
	   Vec<2> p2d;
	   
	   p1d(0) = ip(0);
	   p2d = p1 + tau * p1d;
	   
	   IntegrationPoint ip2d(p2d, 0);
	   
	   testfe.CalcShape (ip, testshape);
	   
	   if (shapenr == 1)
	     CalcShape1 (ip2d, shape);
	   else
	     CalcShape2 (ip2d, shape);
	   
	   shapen.Col(0) = shape * nv;
	   moments += ip.Weight() * (testshape * Trans (shapen));
	 }
     }
   
   else
     
     {
       const POINT3D * points = ElementTopology::GetVertices (ElementType());
       const FACE & face = ElementTopology::GetFaces (ElementType()) [fnr];

       Vec<3> p1, p2, p3, nv, tau1, tau2;

       if (face[3] == -1) // trig
	 for (j = 0; j < 3; j++)
	   {
	     p1(j) = points[face[0]][j];
	     p2(j) = points[face[1]][j];
	     p3(j) = points[face[2]][j];
	     tau1(j) = p1(j) - p3(j);
	     tau2(j) = p2(j) - p3(j);
	   }
       else  // quad
	 ;
       /*
	for (j = 0; j <= 2; j++)
	  {
	    p1[j] = points[face[1]-1][j];
	    p2[j] = points[face[3]-1][j];
	    p3[j] = points[face[0]-1][j];
	    tau1[j] = p1[j] - p3[j];
	    tau2[j] = p2[j] - p3[j];
	  }
       */	
       
       nv = Cross (tau1, tau2);

       /*
      (*testout) << "face " << fnr << " has points " 
		 << face[0] << " " << face[1] << " " << face[2] << " " 
		 << face[3] << endl;
      (*testout) << "p1 = " << p1 << ", p2 = " << p2 << ", p3 = " << p3 << endl;
      (*testout) << "tau1 = " << tau1 << endl;
      (*testout) << "tau2 = " << tau2 << endl;
       */

      for (j = 0; j < facerule.GetNIP(); j++)
	{
	  const IntegrationPoint & ip = facerule[j];
	  
	  Vec<3> p3d;
	  p3d = p3 + ip(0) * tau1 + ip(1) * tau2;

	  IntegrationPoint ip3d(&p3d(0), 0);
	  
	  testfe.CalcShape (ip, testshape);
	  
	  if (shapenr == 1)
	    CalcShape1 (ip3d, shape);
	  else
	    CalcShape2 (ip3d, shape);
	  
	  shapen.Col(0) = shape * nv;
	  moments += ip.Weight() * (testshape * Trans (shapen));
	}
    }
}



  template<int D>
  list<tuple<string,double>> HDivFiniteElement<D> :: Timing () const
  {
    list<tuple<string,double>>timings;
    IntegrationRule ir(ElementType(), Order());
    SIMD_IntegrationRule simdir(ElementType(), Order());
    Matrix<> shape(GetNDof(), D);
    Vector<> coefs(GetNDof());
    Matrix<> values(ir.Size(), D);
    Vector<> divvalues(ir.Size());
    Vector<SIMD<double>> adivvalues(simdir.Size());
    Matrix<SIMD<double>> avalues(D, simdir.Size());
    FE_ElementTransformation<D,D> trafo(ElementType());
    static LocalHeap lh (100000, "FE - Timing");
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
    timings.push_back(make_tuple("CalcShape", time/steps*1e9/(D*GetNDof())));

    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> Evaluate(ir, coefs, values);
                     }, maxtime);
    timings.push_back(make_tuple("Evaluate",time/steps*1e9/(D*GetNDof()*ir.GetNIP())));

    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> Evaluate(simdmir, coefs, avalues);
                     }, maxtime);
    timings.push_back(make_tuple("Evaluate(SIMD)", time/steps*1e9/(D*GetNDof()*ir.GetNIP())));

    /*
    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> EvaluateDiv(mir, coefs, divvalues);
                     }, maxtime);
    timings.push_back(make_tuple("Evaluate Grad", time/steps*1e9/(D*GetNDof()*ir.GetNIP())));
    */
    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> EvaluateDiv(simdmir, coefs, adivvalues);
                     }, maxtime);
    timings.push_back(make_tuple("Evaluate Grad(SIMD)", time/steps*1e9/(GetNDof()*ir.GetNIP())));

    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> EvaluateTrans(ir, values, coefs);
                     }, maxtime);
    timings.push_back(make_tuple("Evaluate Trans", time/steps*1e9/(D*GetNDof()*ir.GetNIP())));

    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> AddTrans(simdmir, avalues, coefs);
                     }, maxtime);
    timings.push_back(make_tuple("Evaluate Trans (SIMD)", time/steps*1e9/(D*GetNDof()*ir.GetNIP())));

    /*
    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> EvaluateDivTrans(mir, divvalues, coefs);
                     }, maxtime);
    timings.push_back(make_tuple("Evaluate Trans Grad", time/steps*1e9/(D*GetNDof()*ir.GetNIP())));
    */
    
    time = RunTiming([&]() {
                     for (size_t i = 0; i < steps; i++)
                       this -> AddDivTrans(simdmir, adivvalues, coefs);
                     }, maxtime);
    timings.push_back(make_tuple("Evaluate Trans Grad(SIMD)", time/steps*1e9/(GetNDof()*ir.GetNIP())));

    return timings;
  }


  template <int D>
  void HDivFiniteElement<D> ::
  CalcDualShape (const BaseMappedIntegrationPoint & bmip, SliceMatrix<> shape) const
  {
    // throw Exception(string("CalcDualShape not implemented for H(div) element ")+typeid(*this).name());
    static bool first = true;
    if (first)
      cerr << "CalcDualShape not implemented for H(div) element " << typeid(*this).name() << endl;
    first = false;
  }

  template class HDivFiniteElement<0>;
  template class HDivFiniteElement<1>;
  template class HDivFiniteElement<2>;
  template class HDivFiniteElement<3>;

  


  FE_RTTrig0 :: FE_RTTrig0()
    : HDivFiniteElement<2> (3, 1)
  {
    ;
  }

  FE_RTTrig0 :: ~FE_RTTrig0()
  {
    ;
  }

  void FE_RTTrig0 :: 
  CalcShape (const IntegrationPoint & ip, 
	     SliceMatrix<> shape) const
  {
    double x = ip(0);
    double y = ip(1);

    /*
    // old orientation, needed for H(curl)
    shape (0,0) = -x;
    shape (0,1) = 1-y;

    shape (1,0) = x-1;
    shape (1,1) = y;

    shape (2,0) = -x;
    shape (2,1) = -y;
    */
    // new orientation, needed for H(div)
    shape (0,0) = x;
    shape (0,1) = y-1;

    shape (1,0) = x-1;
    shape (1,1) = y;

    shape (2,0) = x;
    shape (2,1) = y;
  }





  // Array<HDivFiniteElement<2>::IPData> FE_RTTrig0plus::ipdata;
  
  FE_RTTrig0plus :: FE_RTTrig0plus()
    : HDivFiniteElement<2> (3, 1)
  {
    // CalcIPData(ipdata);
  }

  FE_RTTrig0plus :: ~FE_RTTrig0plus()
  {
    ;
  }

  void FE_RTTrig0plus :: 
  CalcShape (const IntegrationPoint & ip, 
	     SliceMatrix<> shape) const
  {
    double x = ip.Point()[0];
    double y = ip.Point()[1];

    shape (0,0) = -x;
    shape (0,1) = 1-y;

    shape (1,0) = 1-x;
    shape (1,1) = -y;

    shape (2,0) = -x;
    shape (2,1) = -y;

    shape (3,0) = x-x*x-2*x*y;
    shape (3,1) = -y+y*y+2*x*y;
  }


  // Array<HDivFiniteElement<2>::IPData> FE_BDMTrig1::ipdata;
  Matrix<> FE_BDMTrig1::trans(6);
  
  FE_BDMTrig1 :: FE_BDMTrig1()
    : HDivFiniteElement<2> (6, 1)
  {
    Orthogonalize();
    // CalcIPData(ipdata);
  }
  
  FE_BDMTrig1 :: ~FE_BDMTrig1()
  {
    ;
  }
  
  void FE_BDMTrig1 :: CalcShape (const IntegrationPoint & ip, 
				 SliceMatrix<> shape) const
  {
    Mat<6,2> shape1;
    
    CalcShape1 (ip, shape1);
    shape = Trans (trans) * shape1;
  }
  
  void FE_BDMTrig1 :: CalcShape1 (const IntegrationPoint & ip, 
				  FlatMatrixFixWidth<2> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    shape = 0;
    
    shape(0,0) = 1;
    shape(1,0) = x;
    shape(2,0) = y;
    
    shape(3,1) = 1;
    shape(4,1) = x;
    shape(5,1) = y;
  }
  

  void FE_BDMTrig1 :: Orthogonalize()
  {
    cout << "compute BDM trig 1" << endl;
    
    int nd = 6;
    int i, j;
    Matrix<double> fiphij(nd);
    
    Matrix<> edgemoments(2,nd);
    FE_Segm1 segm1;
    
    for (i = 0; i < 3; i++)
      {
	ComputeFaceMoments (i, segm1, edgemoments, 4);
	for (j = 0; j < nd; j++)
	  {
	    fiphij(2*i, j)   = edgemoments(0, j);
	    fiphij(2*i+1, j) = edgemoments(1, j);
	  }
      }
    CalcInverse (fiphij, trans);
    
    (*testout) << "BDMTrig1" << endl
	       << "fiphij = " << endl << fiphij << endl
	       << "trans = " << endl << trans << endl;
  }
  













  HDivNormalSegm0 :: HDivNormalSegm0()
    : HDivNormalFiniteElement<1> (1, 0)
  {
    ;
  }

  void HDivNormalSegm0 :: 
  CalcShape (const IntegrationPoint & ip, 
	     FlatVector<> shape) const
  {
    shape(0) = 1;
  }














  // Array<HDivFiniteElement<2>::IPData> FE_RTQuad0::ipdata;

FE_RTQuad0 :: FE_RTQuad0()
  : HDivFiniteElement<2> (4, 1) {; }

FE_RTQuad0 :: ~FE_RTQuad0()
{
  ;
}

  void FE_RTQuad0 :: 
  CalcShape (const IntegrationPoint & ip, 
             SliceMatrix<> shape) const
{
  double x = ip(0);
  double y = ip(1);
  shape = 0;
  
  shape (0, 1) = 1-y;
  shape (1, 1) = y;
  shape (2, 0) = 1-x;
  shape (3, 0) = x;
}



  // Array<HDivFiniteElement<3>::IPData> FE_BDMTet1::ipdata;
Matrix<> FE_BDMTet1::trans(12);

FE_BDMTet1 :: FE_BDMTet1()
  : HDivFiniteElement<3> (12, 1)
{
  Orthogonalize();
  // CalcIPData(ipdata);
}

FE_BDMTet1 :: ~FE_BDMTet1()
{
  ;
}

  void FE_BDMTet1 :: CalcShape (const IntegrationPoint & ip, 
                                SliceMatrix<> shape) const
{
  Mat<12,3> shape1;
  CalcShape1 (ip, shape1);
  shape = Trans (trans) * shape1;
}


void FE_BDMTet1 :: CalcShape1 (const IntegrationPoint & ip, 
			       FlatMatrixFixWidth<3> shape) const
{
  double x = ip(0);
  double y = ip(1);
  double z = ip(2);
  shape = 0;

  for (int comp = 0; comp < 3; comp++)
    {
      int base = 4 * comp;
      shape (base  ,comp) = 1;
      shape (base+1,comp) = x;
      shape (base+2,comp) = y;
      shape (base+3,comp) = z;
    }
}




void FE_BDMTet1 :: Orthogonalize()
{
  cout << "compute BDM1 tet" << endl;
  
  int nd = 12;
  
  //  const POINT3D * points = MeshAccess::ME_GetVertices (ET_TET);

  Matrix<> fiphij(nd);
  fiphij = 0;
  
  Matrix<> moments(3,nd);
  ScalarFE<ET_TRIG,1> trig1;
  
  for (int i = 0; i < 4; i++)
    {
      ComputeFaceMoments (i, trig1, moments, 2);
      if (i == 0 || i == 2)
	moments *= -1;  // ???

      (*testout) << "moments = " << moments << endl;

      for (int j = 0; j < nd; j++)
	{
	  fiphij(3*i  , j) = moments(0, j);
	  fiphij(3*i+1, j) = moments(1, j);
	  fiphij(3*i+2, j) = moments(2, j);
	}
    }

  (*testout) << "BDMTet1" << endl
	     << "fiphij = " << endl << fiphij << endl;
  
  CalcInverse (fiphij, trans);

  (*testout) << "BDMTet1" << endl
	     << "fiphij = " << endl << fiphij << endl
	     << "trans = " << endl << trans << endl;
}
  


}


