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
  using namespace ngfem;


#ifdef OLD
  template <int D>
  void HDivFiniteElement<D> :: 
  CalcIPData (Array<IPData> & ipdata)
  {
    if (ipdata.Size() == 0)
      {
	const Array<IntegrationPoint*> & ipts = 
	  GetIntegrationRules().GetIntegrationPoints (eltype);
    
	/*
	(*testout) << "Calc HDiv-IP Data for element type  " << eltype
		   << ", ndof = " << GetNDof() 
		   << ": " << ipts.Size() << endl;
	*/

	ipdata.SetSize (ipts.Size());

	DynamicMem<double> * block = new DynamicMem<double> (ipts.Size() * ndof * (DIM+1));
	block->SetName ("HDiv-FiniteElement IPData");
	double * hp = block->Ptr(); 

	for (int i = 0; i < ipts.Size(); i++)
	  {
	    ipdata[i].shape.AssignMemory (ndof, hp);
	    hp += ndof * DIM;
	    ipdata[i].divshape.AssignMemory (ndof, hp);
	    hp += ndof;
	    
	    CalcShape (*ipts[i], ipdata[i].shape);
	    CalcDivShape (*ipts[i], ipdata[i].divshape);
	  }
      }

    p_ipdata = &ipdata[0];
  }
#endif
  
  template <int D>
  void HDivFiniteElement<D> ::
  CalcDivShape (const IntegrationPoint & ip, 
		FlatVector<> divshape) const
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


  /// compute shape
  template <int D>
  void HDivFiniteElement<D> ::
  CalcMappedShape (const SpecificIntegrationPoint<DIM,DIM> & sip,
                   FlatMatrixFixWidth<DIM> shape) const
  {
    CalcShape (sip.IP(), shape);
    Mat<DIM> trans = (1.0/sip.GetJacobiDet()) * sip.GetJacobian();
    for (int i = 0; i < ndof; i++)
      {
        Vec<DIM> hs = shape.Row(i);
        shape.Row(i) = trans * hs;
      }
  }

  /// compute curl of shape
  template <int D>
  void HDivFiniteElement<D> ::
  CalcMappedDivShape (const SpecificIntegrationPoint<DIM,DIM> & sip,
                      FlatVector<> divshape) const
  {
    CalcDivShape (sip.IP(), divshape);
    divshape /= sip.GetJacobiDet();
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

   const IntegrationRule & facerule = 
     GetIntegrationRules().SelectIntegrationRule (testfe.ElementType(), order);
    
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
	   
	   shapen = shape * nv;
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
	  
	  shapen = shape * nv;
	  moments += ip.Weight() * (testshape * Trans (shapen));
	}
    }








   /*

  int i, j, k, l, m;
  Vector shape, testshape;

  if (SpatialDim() == 2)
    {
      const IntegrationRule & linerule = 
	GetIntegrationRules().SelectIntegrationRule (ET_SEGM, order);
      
      const POINT3D * points = MeshAccess::ME_GetVertices (ElementType());
      const EDGE & edge = MeshAccess::ME_GetEdges (ElementType()) [fnr-1];

      double p1[2], p2[2], nu[2], tau[2];
      for (j = 0; j <= 1; j++)
	{
	  p1[j] = points[edge[0]-1][j];
	  p2[j] = points[edge[1]-1][j];
	  tau[j] = p2[j] - p1[j];
	}
      nu[0] = -tau[1];
      nu[1] = tau[0];

      for (j = 1; j <= linerule.GetNIP(); j++)
	{
	  const IntegrationPoint & ip = linerule.GetIP(j);
	  
	  double p[3];
	  for (k = 0; k <= 1; k++)
	    p[k] = p1[k] + tau[k] * ip.Point()[0];
	  p[2] = 0;

	  IntegrationPoint ip3d(p, 0);
	  
	  testfe.CalcShape (ip, testshape);
	  for (k = 1; k <= 2; k++)
	    {
	      if (shapenr == 1)
		CalcShape1 (ip3d, shape, k);
	      else
		CalcShape2 (ip3d, shape, k);

	      if (j == 1 && k == 1)
		{
		  moments.SetSize(testshape.Size(), shape.Size());
		  moments.SetScalar (0);
		}
	      
	      for (l = 1; l <= shape.Size(); l++)
		for (m = 1; m <= testshape.Size(); m++)
		  moments.Elem(m, l) += 
		    nu[k-1] * ip.Weight() * shape.Get(l) * testshape.Get(m);
	    }
	}
    }
  else
    {
      const POINT3D * points = MeshAccess::ME_GetVertices (ElementType());
      const FACE & face = MeshAccess::ME_GetFaces (ElementType()) [fnr-1];

      const IntegrationRule & rule = 
	GetIntegrationRules().SelectIntegrationRule (face[3] ? ET_QUAD : ET_TRIG, order);
      

      double p1[3], p2[3], p3[3], nu[3], tau1[3], tau2[3];

      if (face[3] == 0)
	for (j = 0; j <= 2; j++)
	  {
	    p1[j] = points[face[0]-1][j];
	    p2[j] = points[face[1]-1][j];
	    p3[j] = points[face[2]-1][j];
	    tau1[j] = p1[j] - p3[j];
	    tau2[j] = p2[j] - p3[j];
	  }
      else
	for (j = 0; j <= 2; j++)
	  {
	    p1[j] = points[face[1]-1][j];
	    p2[j] = points[face[3]-1][j];
	    p3[j] = points[face[0]-1][j];
	    tau1[j] = p1[j] - p3[j];
	    tau2[j] = p2[j] - p3[j];
	  }
	

      nu[0] = tau1[1]*tau2[2]-tau1[2]*tau2[1];
      nu[1] = tau1[2]*tau2[0]-tau1[0]*tau2[2];
      nu[2] = tau1[0]*tau2[1]-tau1[1]*tau2[0];

      (*testout) << "face " << fnr << " has points " 
		 << face[0] << " " << face[1] << " " << face[2] << " " 
		 << face[3] << endl;
      (*testout) << "p3 = " << p3[0] << ", " << p3[1] << ", " << p3[2] << endl;
      (*testout) << "tau1 = " << tau1[0] << ", " << tau1[1] << ", " << tau1[2] << endl;
      (*testout) << "tau2 = " << tau2[0] << ", " << tau2[1] << ", " << tau2[2] << endl;
      for (j = 1; j <= rule.GetNIP(); j++)
	{
	  const IntegrationPoint & ip = rule.GetIP(j);
	  
	  double p[3];
	  for (k = 0; k <= 2; k++)
	    p[k] = p3[k] + tau1[k] * ip.Point()[0]  
	      + tau2[k] * ip.Point()[1];

	  IntegrationPoint ip3d(p, 0);
	  
	  testfe.CalcShape (ip, testshape);
	  (*testout) << "testshape = " << testshape << endl;
	  for (k = 1; k <= 3; k++)
	    {
	      if (shapenr == 1)
		CalcShape1 (ip3d, shape, k);
	      else
		CalcShape2 (ip3d, shape, k);

	      if (j == 1 && k == 1)
		{
		  moments.SetSize(testshape.Size(), shape.Size());
		  moments.SetScalar (0);
		}
	      
	      for (l = 1; l <= shape.Size(); l++)
		for (m = 1; m <= testshape.Size(); m++)
		  moments.Elem(m, l) += 
		    nu[k-1] * ip.Weight() * shape.Get(l) * testshape.Get(m);
	    }
	}
    }
  */
}










  



  // Array<HDivFiniteElement<2>::IPData> FE_RTTrig0::ipdata;

  FE_RTTrig0 :: FE_RTTrig0()
    : HDivFiniteElement<2> (ET_TRIG, 3, 1)
  {
    // CalcIPData(ipdata);
  }

  FE_RTTrig0 :: ~FE_RTTrig0()
  {
    ;
  }

  void FE_RTTrig0 :: 
  CalcShape (const IntegrationPoint & ip, 
	     FlatMatrixFixWidth<2> shape) const
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
    : HDivFiniteElement<2> (ET_TRIG, 3, 1)
  {
    // CalcIPData(ipdata);
  }

  FE_RTTrig0plus :: ~FE_RTTrig0plus()
  {
    ;
  }

  void FE_RTTrig0plus :: 
  CalcShape (const IntegrationPoint & ip, 
	     FlatMatrixFixWidth<2> shape) const
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
    : HDivFiniteElement<2> (ET_TRIG, 6, 1)
  {
    Orthogonalize();
    // CalcIPData(ipdata);
  }
  
  FE_BDMTrig1 :: ~FE_BDMTrig1()
  {
    ;
  }
  
  void FE_BDMTrig1 :: CalcShape (const IntegrationPoint & ip, 
				 FlatMatrixFixWidth<2> shape) const
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
    : HDivNormalFiniteElement<1> (ET_SEGM, 1, 0)
  {
    ;
  }

  void HDivNormalSegm0 :: 
  CalcShape (const IntegrationPoint & ip, 
	     FlatVector<> shape) const
  {
    shape(0) = 1;
  }













#ifdef ABC

#ifdef NONE




  
  // Array<IPDataHDiv*> FE_BDMTrig1plus::ipdata;
  DenseMatrix FE_BDMTrig1plus::trans;

  FE_BDMTrig1plus :: FE_BDMTrig1plus()
  {
    if (!ipdata.Size())
      {
	Orthogonalize();
	CalcIPData(ET_TRIG, ipdata);
      }
  }

  FE_BDMTrig1plus :: ~FE_BDMTrig1plus()
  {
    ;
  }

  void FE_BDMTrig1plus :: CalcShape (const IntegrationPoint & ip, 
				 Vector & shape,
				 int comp) const
  {
    Vector shape1;
    CalcShape1 (ip, shape1, comp);
    trans.MultTrans (shape1, shape);
  }

  void FE_BDMTrig1plus :: CalcShape1 (const IntegrationPoint & ip, 
				  Vector & shape,
				  int comp) const
  {
    shape.SetSize(7);
    double x = ip.Point()[0];
    double y = ip.Point()[1];
    shape.SetScalar (0);
    
    switch (comp)
      {
      case 1:
	{
	  shape (1) = 1;
	  shape (2) = x;
	  shape (3) = y;
	  shape (7) = x-x*x-2*x*y;
	  break;
	}
      case 2:
	{
	  shape (4) = 1;
	  shape (5) = x;
	  shape (6) = y;
	  shape (7) = -y+y*y+2*x*y;
	  break;
	}
      }
  }


  void FE_BDMTrig1plus :: Orthogonalize()
  {
    cout << "compute BDM trig 1" << endl;

    int i, j, k, l;
    int nd = 7;

    const double points[3][3] =
    { { 1, 0, 0 },
      { 0, 1, 0 },
      { 0, 0, 0 } };
    const int edges[3][2] = 
    { { 3, 1 },
      { 3, 2 },
      { 1, 2 } };

    double ipx[3] = { 0, 0.5, 1 };
    double ipw[3] = { 1.0/6.0, 2.0/3.0, 1.0/6.0 };

    DenseMatrix fiphij(nd);
    Vector shape1(nd), sum1(nd), sum2(nd);
    
    fiphij.SetScalar (0);
    for (i = 1; i <= 3; i++)
      {
	double p1[2], p2[2], nu[2], tau[2];
	for (j = 0; j <= 1; j++)
	  {
	    p1[j] = points[edges[i-1][0]-1][j];
	    p2[j] = points[edges[i-1][1]-1][j];
	    tau[j] = p2[j] - p1[j];
	  }
	nu[0] = -tau[1];
	nu[1] = tau[0];
	
	sum1.SetScalar (0);
	sum2.SetScalar (0);

	for (j = 0; j <= 2; j++)
	  {
	    double p[3];
	    for (k = 0; k <= 1; k++)
	      p[k] = p1[k] + tau[k] * ipx[j];
	    p[2] = 0;

	    IntegrationPoint ip(p, ipw[j]);

	    for (k = 1; k <= 2; k++)
	      {
		CalcShape1 (ip, shape1, k);
		sum1.Add (nu[k-1] * ipw[j], shape1);
		sum2.Add ( (2*ipx[j]-1) * nu[k-1] * ipw[j], shape1);
	      }
	  }

	for (j = 1; j <= nd; j++)
	  {
	    fiphij.Elem(i, j) = sum1.Get(j);
	    fiphij.Elem(3+i, j) = sum2.Get(j);
	  }
      }
    fiphij.Elem(7, 7) = 1;

    trans.SetSize (nd);
    CalcInverse (fiphij, trans);
  }
  















  Array<IPDataHDiv*> FE_BDFMTrig2::ipdata;
  DenseMatrix FE_BDFMTrig2::trans;

  FE_BDFMTrig2 :: FE_BDFMTrig2()
  {
    if (!ipdata.Size())
      {
	Orthogonalize();
	CalcIPData(ET_TRIG, ipdata);
      }
  }

  FE_BDFMTrig2 :: ~FE_BDFMTrig2()
  {
    ;
  }

  void FE_BDFMTrig2 :: CalcShape (const IntegrationPoint & ip, 
				 Vector & shape,
				 int comp) const
  {
    Vector shape1;
    CalcShape1 (ip, shape1, comp);
    trans.MultTrans (shape1, shape);
  }

  void FE_BDFMTrig2 :: CalcShape1 (const IntegrationPoint & ip, 
				   Vector & shape,
				   int comp) const
  {
    shape.SetSize(12);
    shape.SetScalar (0);

    double x = ip.Point()[0];
    double y = ip.Point()[1];

    int base = 6*(comp-1);

    shape (base+1) = 1;
    shape (base+2) = x;
    shape (base+3) = y;
    shape (base+4) = x*x;
    shape (base+5) = x*y;
    shape (base+6) = y*y;
  }


  void FE_BDFMTrig2 :: Orthogonalize()
  {
    cout << "compute BDFM2 trig" << endl;

    int i, j, k, l;
    int nd = 12;

    const double points[3][3] =
    { { 1, 0, 0 },
      { 0, 1, 0 },
      { 0, 0, 0 } };

    const IntegrationRule & trigrule = 
      GetIntegrationRules().SelectIntegrationRule (ET_TRIG, 5);

    DenseMatrix fiphij(nd);
    fiphij.SetScalar (0);
    Vector shape1(nd), sum1(nd), sum2(nd), sum3(nd);

    DenseMatrix edgemoments;
    FE_Segm2L2 segm2;

    for (i = 1; i <= 3; i++)
      {
	ComputeFaceMoments (i, segm2, edgemoments, 4);

	for (j = 1; j <= nd; j++)
	  {
	    fiphij.Elem(i, j) = edgemoments.Get(1, j);
	    fiphij.Elem(3+i, j) = edgemoments.Get(2, j);
	    fiphij.Elem(6+i, j) = edgemoments.Get(3, j);
	  }
      }


    // test with nedelec function
    double p1[3], p2[3], p3[3], p0[3], tau1[3], tau2[3];
    for (j = 0; j <= 2; j++)
      {
	p1[j] = points[0][j];
	p2[j] = points[1][j];
	p3[j] = points[2][j];
      }

    for (k = 1; k <= 3; k++) 
      {  

	for (j = 0; j <= 2; j++)
	  {
	    switch (k)
	      {
	      case 1:
		p0[j] = p1[j];
		tau1[j] = p2[j] - p1[j];
		tau2[j] = p3[j] - p1[j];
		break;
	      case 2:
		p0[j] = p2[j];
		tau1[j] = p3[j] - p2[j];
		tau2[j] = p1[j] - p2[j];
		break;
	      case 3:
		p0[j] = p3[j];
		tau1[j] = p1[j] - p3[j];
		tau2[j] = p2[j] - p3[j];
		break;
	      }
	  }
	
	sum1.SetScalar (0);
	

	for (l = 1; l <= trigrule.GetNIP(); l++)
	  { // integration points
		
	    const IntegrationPoint & tip = trigrule.GetIP(l);
	    
	    double p[3];
	    for (j = 0; j <= 2; j++)
	      {
		p[j] = p0[j] + 
		  tip.Point()[0] * tau1[j] + 
		  tip.Point()[1] * tau2[j];
	      }
	    
	    IntegrationPoint ip(p, 1.0/6.0);

	    double rt0[2];
	    rt0[0] = -tip.Point()[1];
	    rt0[1] = tip.Point()[0];
	    
	    for (j = 1; j <= 2; j++)
	      {
		CalcShape1 (ip, shape1, j);
		sum1.Add ((rt0[0] * tau1[j-1] + rt0[1] * tau2[j-1]) * 
			  tip.Weight(), shape1);
	      }
	  }

	for (j = 1; j <= nd; j++)
	  fiphij.Elem(9 + k, j) = sum1.Get(j);	      
      }


    DenseMatrix trans1 (nd);
    CalcInverse (fiphij, trans1);


    trans.SetSize(12, 9);
    for (i = 1; i <= 12; i++)
      {
	for (j = 1; j <= 6; j++)
	  trans.Elem(i, j) = trans1.Elem(i,j);
	for (j = 1; j <= 3; j++)
	  trans.Elem(i, 6+j) = trans1.Elem(i, 9+j);
      }
    cout << "done" << endl;
  }










  
  Array<IPDataHDiv*> FE_BDMTrig2::ipdata;
  DenseMatrix FE_BDMTrig2::trans;

  FE_BDMTrig2 :: FE_BDMTrig2()
  {
    if (!ipdata.Size())
      {
	Orthogonalize();
	CalcIPData(ET_TRIG, ipdata);
      }
  }

  FE_BDMTrig2 :: ~FE_BDMTrig2()
  {
    ;
  }

  void FE_BDMTrig2 :: CalcShape (const IntegrationPoint & ip, 
				 Vector & shape,
				 int comp) const
  {
    Vector shape1;
    CalcShape1 (ip, shape1, comp);
    trans.MultTrans (shape1, shape);
    
    double x = ip.Point()[0];
    double y = ip.Point()[1];

    // lowest order RT:
    switch (comp)
      {
      case 1:
	{
	  shape (1) = -x;
	  shape (2) = x-1;
	  shape (3) = -x;
	  break;
	}
      case 2:
	{
	  shape (1) = 1-y;
	  shape (2) = y;
	  shape (3) = -y;
	  break;
	}
      }
  }

  void FE_BDMTrig2 :: CalcShape1 (const IntegrationPoint & ip, 
				  Vector & shape,
				  int comp) const
  {
    shape.SetSize(12);
    double x = ip.Point()[0];
    double y = ip.Point()[1];
    shape.SetScalar (0);
    
    switch (comp)
      {
      case 1:
	{
	  shape (1) = 1;
	  shape (2) = x;
	  shape (3) = y;
	  shape (4) = x*x;
	  shape (5) = x*y;
	  shape (6) = y*y;
	  break;
	}
      case 2:
	{
	  shape (7) = 1;
	  shape (8) = x;
	  shape (9) = y;
	  shape (10) = x*x;
	  shape (11) = x*y;
	  shape (12) = y*y;
	  break;
	}
      }
  }


  void FE_BDMTrig2 :: Orthogonalize()
  {
    cout << "compute BDM trig" << endl;

    int i, j, k, l;
    int nd = 12;

    const double points[3][3] =
    { { 1, 0, 0 },
      { 0, 1, 0 },
      { 0, 0, 0 } };
    const int edges[3][2] = 
    { { 3, 1 },
      { 3, 2 },
      { 1, 2 } };

    const IntegrationRule & linerule = 
      GetIntegrationRules().SelectIntegrationRule (ET_SEGM, 4);
    const IntegrationRule & trigrule = 
      GetIntegrationRules().SelectIntegrationRule (ET_TRIG, 5);


    DenseMatrix fiphij(nd);
    fiphij.SetScalar (0);
    Vector shape1(nd), sum1(nd), sum2(nd), sum3(nd);

    for (i = 1; i <= 3; i++)
      {
	double p1[2], p2[2], nu[2], tau[2];
	for (j = 0; j <= 1; j++)
	  {
	    p1[j] = points[edges[i-1][0]-1][j];
	    p2[j] = points[edges[i-1][1]-1][j];
	    tau[j] = p2[j] - p1[j];
	  }
	nu[0] = -tau[1];
	nu[1] = tau[0];
	
	sum1.SetScalar (0);
	sum2.SetScalar (0);
	sum3.SetScalar (0);

	for (j = 1; j <= linerule.GetNIP(); j++)
	  {
	    const IntegrationPoint & ipl = linerule.GetIP(j);

	    double p[3];
	    double xi = ipl.Point()[0];

	    for (k = 0; k <= 1; k++)
	      p[k] = p1[k] + tau[k] * xi;
	    p[2] = 0;

	    IntegrationPoint ip(p, 1);

	    for (k = 1; k <= 2; k++)
	      {
		CalcShape1 (ip, shape1, k);
		sum1.Add (nu[k-1] * ipl.Weight(), shape1);
		sum2.Add ( (2*xi-1) * nu[k-1] * ipl.Weight(), shape1);
		sum3.Add ( (sqr(2*xi-1)-1.0/3.0) * nu[k-1] * ipl.Weight(), shape1);
	      }
	  }

	for (j = 1; j <= nd; j++)
	  {
	    fiphij.Elem(i, j) = sum1.Get(j);
	    fiphij.Elem(3+i, j) = sum2.Get(j);
	    fiphij.Elem(6+i, j) = sum3.Get(j);
	  }
      }



    // test with nedelec function
    double p1[3], p2[3], p3[3], p0[3], tau1[3], tau2[3];
    for (j = 0; j <= 2; j++)
      {
	p1[j] = points[0][j];
	p2[j] = points[1][j];
	p3[j] = points[2][j];
      }

    for (k = 1; k <= 3; k++) 
      {  

	for (j = 0; j <= 2; j++)
	  {
	    switch (k)
	      {
	      case 1:
		p0[j] = p1[j];
		tau1[j] = p2[j] - p1[j];
		tau2[j] = p3[j] - p1[j];
		break;
	      case 2:
		p0[j] = p2[j];
		tau1[j] = p3[j] - p2[j];
		tau2[j] = p1[j] - p2[j];
		break;
	      case 3:
		p0[j] = p3[j];
		tau1[j] = p1[j] - p3[j];
		tau2[j] = p2[j] - p3[j];
		break;
	      }
	  }
	
	sum1.SetScalar (0);
	

	for (l = 1; l <= trigrule.GetNIP(); l++)
	  { // integration points
		
	    const IntegrationPoint & tip = trigrule.GetIP(l);
	    
	    double p[3];
	    for (j = 0; j <= 2; j++)
	      {
		p[j] = p0[j] + 
		  tip.Point()[0] * tau1[j] + 
		  tip.Point()[1] * tau2[j];
	      }
	    
	    IntegrationPoint ip(p, 1.0/6.0);

	    double rt0[2];
	    rt0[0] = -tip.Point()[1];
	    rt0[1] = tip.Point()[0];
	    
	    for (j = 1; j <= 2; j++)
	      {
		CalcShape1 (ip, shape1, j);
		sum1.Add ((rt0[0] * tau1[j-1] + rt0[1] * tau2[j-1]) * 
			  tip.Weight(), shape1);
	      }
	  }

	for (j = 1; j <= nd; j++)
	  fiphij.Elem(9 + k, j) = sum1.Get(j);	      
      }


    trans.SetSize (nd);
    CalcInverse (fiphij, trans);
  }
  


  
  Array<IPDataHDiv*> FE_BDMTrig2plus::ipdata;
  DenseMatrix FE_BDMTrig2plus::trans;

  FE_BDMTrig2plus :: FE_BDMTrig2plus()
  {
    if (!ipdata.Size())
      {
	Orthogonalize();
	CalcIPData(ET_TRIG, ipdata);
      }
  }

  FE_BDMTrig2plus :: ~FE_BDMTrig2plus()
  {
    ;
  }

  void FE_BDMTrig2plus :: CalcShape (const IntegrationPoint & ip, 
				 Vector & shape,
				 int comp) const
  {
    Vector shape1;
    CalcShape1 (ip, shape1, comp);
    trans.MultTrans (shape1, shape);
  }

  void FE_BDMTrig2plus :: CalcShape1 (const IntegrationPoint & ip, 
				  Vector & shape,
				  int comp) const
  {
    shape.SetSize(14);
    double x = ip.Point()[0];
    double y = ip.Point()[1];
    shape.SetScalar (0);
    
    switch (comp)
      {
      case 1:
	{
	  shape (1) = 1;
	  shape (2) = x;
	  shape (3) = y;
	  shape (4) = x*x;
	  shape (5) = x*y;
	  shape (6) = y*y;
	  shape (13) = x*x-x*x*x-2*x*x*y;
	  shape (14) = 2*x*y - 2*x*x*y - 3*x*y*y;
	  break;
	}
      case 2:
	{
	  shape (7) = 1;
	  shape (8) = x;
	  shape (9) = y;
	  shape (10) = x*x;
	  shape (11) = x*y;
	  shape (12) = y*y;
	  shape (13) = -2*x*y + 3*x*x*y + 2*x*y*y;
	  shape (14) = -y*y + 2*x*y*y + y*y*y;
	  break;
	}
      }
  }


  void FE_BDMTrig2plus :: Orthogonalize()
  {
    cout << "compute BDM trig 2+" << endl;



    int i, j, k, l;
    int nd = 14;

    const double points[3][3] =
    { { 1, 0, 0 },
      { 0, 1, 0 },
      { 0, 0, 0 } };
    const int edges[3][2] = 
    { { 3, 1 },
      { 3, 2 },
      { 1, 2 } };

    const IntegrationRule & linerule = 
      GetIntegrationRules().SelectIntegrationRule (ET_SEGM, 4);
    const IntegrationRule & trigrule = 
      GetIntegrationRules().SelectIntegrationRule (ET_TRIG, 5);


    DenseMatrix fiphij(nd);
    fiphij.SetScalar (0);
    Vector shape1(nd), sum1(nd), sum2(nd), sum3(nd);

    for (i = 1; i <= 3; i++)
      {
	double p1[2], p2[2], nu[2], tau[2];
	for (j = 0; j <= 1; j++)
	  {
	    p1[j] = points[edges[i-1][0]-1][j];
	    p2[j] = points[edges[i-1][1]-1][j];
	    tau[j] = p2[j] - p1[j];
	  }
	nu[0] = -tau[1];
	nu[1] = tau[0];
	
	sum1.SetScalar (0);
	sum2.SetScalar (0);
	sum3.SetScalar (0);

	for (j = 1; j <= linerule.GetNIP(); j++)
	  {
	    const IntegrationPoint & ipl = linerule.GetIP(j);

	    double p[3];
	    double xi = ipl.Point()[0];

	    for (k = 0; k <= 1; k++)
	      p[k] = p1[k] + tau[k] * xi;
	    p[2] = 0;

	    IntegrationPoint ip(p, 1);

	    for (k = 1; k <= 2; k++)
	      {
		CalcShape1 (ip, shape1, k);
		sum1.Add (nu[k-1] * ipl.Weight(), shape1);
		sum2.Add ( (2*xi-1) * nu[k-1] * ipl.Weight(), shape1);
		sum3.Add ( (sqr(2*xi-1)-1.0/3.0) * nu[k-1] * ipl.Weight(), shape1);
	      }
	  }

	for (j = 1; j <= nd; j++)
	  {
	    fiphij.Elem(i, j) = sum1.Get(j);
	    fiphij.Elem(3+i, j) = sum2.Get(j);
	    fiphij.Elem(6+i, j) = sum3.Get(j);
	  }
      }



    // test with nedelec function
    double p1[3], p2[3], p3[3], p0[3], tau1[3], tau2[3];
    for (j = 0; j <= 2; j++)
      {
	p1[j] = points[0][j];
	p2[j] = points[1][j];
	p3[j] = points[2][j];
      }

    for (k = 1; k <= 3; k++) 
      {  

	for (j = 0; j <= 2; j++)
	  {
	    switch (k)
	      {
	      case 1:
		p0[j] = p1[j];
		tau1[j] = p2[j] - p1[j];
		tau2[j] = p3[j] - p1[j];
		break;
	      case 2:
		p0[j] = p2[j];
		tau1[j] = p3[j] - p2[j];
		tau2[j] = p1[j] - p2[j];
		break;
	      case 3:
		p0[j] = p3[j];
		tau1[j] = p1[j] - p3[j];
		tau2[j] = p2[j] - p3[j];
		break;
	      }
	  }
	
	sum1.SetScalar (0);
	

	for (l = 1; l <= trigrule.GetNIP(); l++)
	  { // integration points
		
	    const IntegrationPoint & tip = trigrule.GetIP(l);
	    
	    double p[3];
	    for (j = 0; j <= 2; j++)
	      {
		p[j] = p0[j] + 
		  tip.Point()[0] * tau1[j] + 
		  tip.Point()[1] * tau2[j];
	      }
	    
	    IntegrationPoint ip(p, 1.0/6.0);

	    double rt0[2];
	    rt0[0] = -tip.Point()[1];
	    rt0[1] = tip.Point()[0];
	    
	    for (j = 1; j <= 2; j++)
	      {
		CalcShape1 (ip, shape1, j);
		sum1.Add ((rt0[0] * tau1[j-1] + rt0[1] * tau2[j-1]) * 
			  tip.Weight(), shape1);
	      }
	  }

	for (j = 1; j <= nd; j++)
	  fiphij.Elem(9 + k, j) = sum1.Get(j);	      
      }

    for (i = 13; i <= 14; i++)
      fiphij.Elem(i,i) = 1;

    trans.SetSize (nd);
    CalcInverse (fiphij, trans);
  }
  




#endif

#endif


  // Array<HDivFiniteElement<2>::IPData> FE_RTQuad0::ipdata;

FE_RTQuad0 :: FE_RTQuad0()
  : HDivFiniteElement<2> (ET_QUAD, 4, 1)
{
  // CalcIPData(ipdata);
}

FE_RTQuad0 :: ~FE_RTQuad0()
{
  ;
}

void FE_RTQuad0 :: 
CalcShape (const IntegrationPoint & ip, 
	   FlatMatrixFixWidth<2> shape) const
{
  double x = ip(0);
  double y = ip(1);
  shape = 0;
  
  shape (0, 1) = 1-y;
  shape (1, 1) = y;
  shape (2, 0) = 1-x;
  shape (3, 0) = x;
}


#ifdef OLD
#ifdef abc


  
  // Array<IPDataHDiv*> FE_BDMQuad1::ipdata;
  DenseMatrix FE_BDMQuad1::trans;

  FE_BDMQuad1 :: FE_BDMQuad1()
  {
    // if (!ipdata.Size())
      {
	Orthogonalize();
	// CalcIPData(ET_QUAD, ipdata);
      }
  }

  FE_BDMQuad1 :: ~FE_BDMQuad1()
  {
    ;
  }

  void FE_BDMQuad1 :: CalcShape (const IntegrationPoint & ip, 
				 Vector & shape,
				 int comp) const
  {
    Vector shape1;
    CalcShape1 (ip, shape1, comp);
    trans.MultTrans (shape1, shape);
  }

  void FE_BDMQuad1 :: CalcShape1 (const IntegrationPoint & ip, 
				  Vector & shape,
				  int comp) const
  {
    shape.SetSize(8);
    double x = ip.Point()[0];
    double y = ip.Point()[1];
    shape.SetScalar (0);
    
    switch (comp)
      {
      case 1:
	{
	  shape (1) = 1;
	  shape (2) = x;
	  shape (3) = y;
	  shape (4) = x*y;
	  break;
	}
      case 2:
	{
	  shape (5) = 1;
	  shape (6) = x;
	  shape (7) = y;
	  shape (8) = x*y;
	  break;
	}
      }
  }


  void FE_BDMQuad1 :: Orthogonalize()
  {
    cout << "compute BDM trig 1" << endl;

    int i, j, k, l;
    int nd = 8;

    const double points[4][3] =
    { { 0, 0, 0 },
      { 1, 0, 0 },
      { 1, 1, 0 },
      { 0, 1, 0 } };
    const int edges[4][2] = 
    { { 1, 2 },
      { 4, 3 },
      { 1, 4 },
      { 2, 3 } };

    double ipx[3] = { 0, 0.5, 1 };
    double ipw[3] = { 1.0/6.0, 2.0/3.0, 1.0/6.0 };

    DenseMatrix fiphij(nd);
    Vector shape1(nd), sum1(nd), sum2(nd);

    for (i = 1; i <= 4; i++)
      {
	double p1[2], p2[2], nu[2], tau[2];
	for (j = 0; j <= 1; j++)
	  {
	    p1[j] = points[edges[i-1][0]-1][j];
	    p2[j] = points[edges[i-1][1]-1][j];
	    tau[j] = p2[j] - p1[j];
	  }
	nu[0] = p1[1]-p2[1];
	nu[1] = p2[0]-p1[0];

	sum1.SetScalar (0);
	sum2.SetScalar (0);

	for (j = 0; j <= 2; j++)
	  {
	    double p[3];
	    for (k = 0; k <= 1; k++)
	      p[k] = p1[k] + tau[k] * ipx[j];
	    p[2] = 0;

	    IntegrationPoint ip(p, ipw[j]);

	    for (k = 1; k <= 2; k++)
	      {
		CalcShape1 (ip, shape1, k);
		sum1.Add (nu[k-1] * ipw[j], shape1);
		sum2.Add ( (2*ipx[j]-1) * nu[k-1] * ipw[j], shape1);
	      }
	  }

	for (j = 1; j <= nd; j++)
	  {
	    fiphij.Elem(i, j) = sum1.Get(j);
	    fiphij.Elem(4+i, j) = sum2.Get(j);
	  }
      }

    trans.SetSize (nd);
    CalcInverse (fiphij, trans);
  }
  















  // Array<IPDataHDiv*> FE_RTSegm0::ipdata;

  FE_RTSegm0 :: FE_RTSegm0()
  {
    // if (!ipdata.Size())
    // CalcIPData(ET_SEGM, ipdata);
  }

  FE_RTSegm0 :: ~FE_RTSegm0()
  {
    ;
  }

  void FE_RTSegm0 :: 
  CalcShape (const IntegrationPoint & ip, 
	     Vector & shape,
	     int comp) const
  {
    shape.SetSize(1);
    shape.SetScalar (0);
    double x = ip.Point()[0];

    shape (1) = 1;
  }


  void FE_RTSegm0 :: CalcDShape (const IntegrationPoint & ip, 
			      DenseMatrix & dshape,
			      int comp) const
  {
    dshape.SetSize (1, 1);
    dshape.SetScalar (0);
  }













  // Array<IPDataHDiv*> FE_RTSegm1::ipdata;

  FE_RTSegm1 :: FE_RTSegm1()
  {
    // if (!ipdata.Size())
    // CalcIPData(ET_SEGM, ipdata);
  }

  FE_RTSegm1 :: ~FE_RTSegm1()
  {
    ;
  }

  void FE_RTSegm1 :: 
  CalcShape (const IntegrationPoint & ip, 
	     Vector & shape,
	     int comp) const
  {
    shape.SetSize(2);
    shape.SetScalar (0);
    double x = ip.Point()[0];

    shape (1) = 1;
    shape (2) = x-0.5;
  }


  void FE_RTSegm1 :: CalcDShape (const IntegrationPoint & ip, 
			      DenseMatrix & dshape,
			      int comp) const
  {
    dshape.SetSize (1, 2);
    dshape.Elem(1,1) = 0;
    dshape.Elem(1,2) = 1;
  }










  // Array<IPDataHDiv*> FE_RTSegm2::ipdata;

  FE_RTSegm2 :: FE_RTSegm2()
  {
    // if (!ipdata.Size())
    // CalcIPData(ET_SEGM, ipdata);
  }

  FE_RTSegm2 :: ~FE_RTSegm2()
  {
    ;
  }

  void FE_RTSegm2 :: 
  CalcShape (const IntegrationPoint & ip, 
	     Vector & shape,
	     int comp) const
  {
    shape.SetSize(3);
    shape.SetScalar (0);
    double x = ip.Point()[0];

    shape (1) = 1;
    shape (2) = x-0.5;
    shape (3) = sqr (x-0.5) - 1.0/12.0;
  }


  void FE_RTSegm2 :: CalcDShape (const IntegrationPoint & ip, 
			      DenseMatrix & dshape,
			      int comp) const
  {
    double x = ip.Point()[0];

    dshape.SetSize (1, 3);
    dshape.Elem(1,1) = 0;
    dshape.Elem(1,2) = 1;
    dshape.Elem(1,3) = 2*x-1;
  }












#endif
#endif


  // Array<HDivFiniteElement<3>::IPData> FE_BDMTet1::ipdata;
Matrix<> FE_BDMTet1::trans(12);

FE_BDMTet1 :: FE_BDMTet1()
  : HDivFiniteElement<3> (ET_TET, 12, 1)
{
  Orthogonalize();
  // CalcIPData(ipdata);
}

FE_BDMTet1 :: ~FE_BDMTet1()
{
  ;
}

void FE_BDMTet1 :: CalcShape (const IntegrationPoint & ip, 
			      FlatMatrixFixWidth<3> shape) const
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
  FE_Trig1 trig1;
  
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
  







#ifdef ABC




#ifdef none

  
  Array<IPDataHDiv*> FE_BDFMTet2::ipdata;
  DenseMatrix FE_BDFMTet2::trans;
  DenseMatrix FE_BDFMTet2::trans2;

  FE_BDFMTet2 :: FE_BDFMTet2()
  {
    if (!ipdata.Size())
      {
	Orthogonalize();
	CalcIPData(ET_TET, ipdata);
      }
  }

  FE_BDFMTet2 :: ~FE_BDFMTet2()
  {
    ;
  }

  void FE_BDFMTet2 :: CalcShape (const IntegrationPoint & ip, 
				 Vector & shape,
				 int comp) const
  {
    Vector shape1, shape2, hshape(12);
    shape.SetSize(18);
    
    CalcShape1 (ip, shape1, comp);
    trans.MultTrans (shape1, hshape);
    shape.SetPart (13, hshape);
    
    CalcShape2 (ip, shape2, comp);
    trans2.MultTrans (shape2, hshape);
    shape.SetPart (1, hshape);
  }


  void FE_BDFMTet2 :: CalcShape1 (const IntegrationPoint & ip, 
				  Vector & shape,
				  int comp) const
  {
    shape.SetSize(30);
    double x = ip.Point()[0];
    double y = ip.Point()[1];
    double z = ip.Point()[2];
    shape.SetScalar (0);

    int base = 10 * (comp-1);

    shape (base+1) = 1;
    shape (base+2) = x;
    shape (base+3) = y;
    shape (base+4) = z;
    shape (base+5) = x*x;
    shape (base+6) = x*y;
    shape (base+7) = x*z;
    shape (base+8) = y*y;
    shape (base+9) = y*z;
    shape (base+10) = z*z;
  }


  void FE_BDFMTet2 :: CalcShape2 (const IntegrationPoint & ip, 
				  Vector & shape,
				  int comp) const
  {
    shape.SetSize(12);
    double x = ip.Point()[0];
    double y = ip.Point()[1];
    double z = ip.Point()[2];
    shape.SetScalar (0);

    int base = 4 * (comp-1);

    shape (base+1) = 1;
    shape (base+2) = x;
    shape (base+3) = y;
    shape (base+4) = z;
  }




  void FE_BDFMTet2 :: Orthogonalize()
  {
    cout << "compute BDFM2 tet" << endl;

    int i, j, k, l;
    int nd = 30;

    const POINT3D * points = MeshAccess::ME_GetVertices (ET_TET);

    const IntegrationRule & tetrule = 
      GetIntegrationRules().SelectIntegrationRule (ET_TET, 5);

    DenseMatrix fiphij(nd);
    fiphij.SetScalar (0);
    Vector shape1(nd), sum1(nd), sum2(nd), sum3(nd);
    Vector sum4(nd), sum5(nd), sum6(nd);

    // all quadratic on faces
    DenseMatrix moments;
    FE_Trig2 trig2;

    for (i = 1; i <= 4; i++)
      {
	ComputeFaceMoments (i, trig2, moments, 4, 1);
	if (i == 1 || i == 3)
	  moments.Scale (-1);

	for (j = 1; j <= nd; j++)
	  {
	    fiphij.Elem(3*(i-1)+1, j) = moments.Get(1, j);
	    fiphij.Elem(3*(i-1)+2, j) = moments.Get(2, j);
	    fiphij.Elem(3*(i-1)+3, j) = moments.Get(3, j);
	    fiphij.Elem(3*(i-1)+13, j) = moments.Get(4, j);
	    fiphij.Elem(3*(i-1)+14, j) = moments.Get(5, j);
	    fiphij.Elem(3*(i-1)+15, j) = moments.Get(6, j);
	  }
      }
    
    // test with volume constants and hdiv-bubbles:
    sum1.SetScalar(0);
    sum2.SetScalar(0);
    sum3.SetScalar(0);
    sum4.SetScalar(0);
    sum5.SetScalar(0);
    sum6.SetScalar(0);

    for (j = 1; j <= tetrule.GetNIP(); j++)
      {
	const IntegrationPoint & ip = tetrule.GetIP(j);
	
	double ipx = ip.Point()[0];
	double ipy = ip.Point()[1];
	double ipz = ip.Point()[1];
	double l4 = 1-ipx-ipy-ipz;

	for (k = 1; k <= 3; k++)
	  {
	    CalcShape1 (ip, shape1, k);
	    // eta1 = (1,0,0)
	    if (k == 1)
	      sum1.Add (ip.Weight(), shape1);

	    if (k == 2)
	      sum2.Add (ip.Weight(), shape1);

	    if (k == 3)
	      sum3.Add (ip.Weight(), shape1);

	    // eta4 = (0, y(l4-z), z(y-l4))
	    if (k == 2)
	      sum4.Add (ipy*(l4-ipz)*ip.Weight(), shape1);
	    if (k == 3)
	      sum4.Add (ipz*(ipy-l4)*ip.Weight(), shape1);

	    if (k == 1)
	      sum5.Add (ipx*(l4-ipz)*ip.Weight(), shape1);
	    if (k == 3)
	      sum5.Add (ipz*(ipx-l4)*ip.Weight(), shape1);

	    if (k == 1)
	      sum6.Add (ipx*(l4-ipy)*ip.Weight(), shape1);
	    if (k == 2)
	      sum6.Add (ipy*(ipx-l4)*ip.Weight(), shape1);
	  }
      }

    for (j = 1; j <= nd; j++)
      {
	fiphij.Elem(25, j) = sum1.Get(j);
	fiphij.Elem(26, j) = sum2.Get(j);
	fiphij.Elem(27, j) = sum3.Get(j);
	fiphij.Elem(28, j) = sum4.Get(j);
	fiphij.Elem(29, j) = sum5.Get(j);
	fiphij.Elem(30, j) = sum6.Get(j);
      }
    


    DenseMatrix trans1 (nd);
    CalcInverse (fiphij, trans1);

    trans.SetSize(30, 6);
    for (i = 1; i <= 30; i++)
      {
	for (j = 1; j <= 6; j++)
	  trans.Elem(i, j) = trans1.Elem(i, 24+j);
      }




    // only linear
    nd = 12;
    DenseMatrix fiphij2(nd);
    fiphij.SetScalar (0);

    FE_Trig1 trig1;

    for (i = 1; i <= 4; i++)
      {
	ComputeFaceMoments (i, trig1, moments, 4, 2);
	if (i == 1 || i == 3)
	  moments.Scale (-1);

	for (j = 1; j <= nd; j++)
	  {
	    fiphij2.Elem(3*(i-1)+1, j) = moments.Get(1, j);
	    fiphij2.Elem(3*(i-1)+2, j) = moments.Get(2, j);
	    fiphij2.Elem(3*(i-1)+3, j) = moments.Get(3, j);
	  }
      }
    
    trans2.SetSize(nd);
    CalcInverse (fiphij2, trans2);
  }
  









  
  Array<IPDataHDiv*> FE_BDMPrism1::ipdata;
  DenseMatrix FE_BDMPrism1::trans;

  FE_BDMPrism1 :: FE_BDMPrism1()
  {
    if (!ipdata.Size())
      {
	Orthogonalize();
	CalcIPData(ET_PRISM, ipdata);
      }
  }

  FE_BDMPrism1 :: ~FE_BDMPrism1()
  {
    ;
  }

  void FE_BDMPrism1 :: CalcShape (const IntegrationPoint & ip, 
				 Vector & shape,
				 int comp) const
  {
    Vector shape1;
    CalcShape1 (ip, shape1, comp);
    trans.MultTrans (shape1, shape);
  }


void FE_BDMPrism1 :: CalcShape1 (const IntegrationPoint & ip, 
				  Vector & shape,
				  int comp) const
{
  int i, j;
  shape.SetSize(18);
  shape.SetScalar (0);

  double x = ip.Point()[0];
  double y = ip.Point()[1];
  double z = ip.Point()[2];
    
  int base = 6*(comp-1);
  
  shape.Elem (base+1) = 1;
  shape.Elem (base+2) = x;
  shape.Elem (base+3) = y;
  shape.Elem (base+4) = z;
  shape.Elem (base+5) = x*z;
  shape.Elem (base+6) = y*z;
}


void FE_BDMPrism1 :: Orthogonalize()
{
  cout << "compute BDM prism1" << endl;
  
  int i, j, k, l;
  int nd = 18;
  
  DenseMatrix fiphij(nd);
  fiphij.SetScalar (0);

  DenseMatrix moments;
  FE_Trig1 trig1;
  FE_Quad1 quad1;

  // 1..6: trig dofs
  for (i = 1; i <= 2; i++)
    {
      ComputeFaceMoments (i, trig1, moments, 5);
      if (i == 1)
	moments.Scale (-1);
      for (j = 1; j <= nd; j++)
	{
	  fiphij.Elem(3*(i-1)+1, j) = moments.Get(1, j);
	  fiphij.Elem(3*(i-1)+2, j) = moments.Get(2, j);
	  fiphij.Elem(3*(i-1)+3, j) = moments.Get(3, j);
	}
    }

  // 7..18 quad dofs
  for (i = 3; i <= 5; i++)
    {
      ComputeFaceMoments (i, quad1, moments, 5);
      for (j = 1; j <= nd; j++)
	{
	  fiphij.Elem(6+4*(i-3)+1, j) = moments.Get(1, j);
	  fiphij.Elem(6+4*(i-3)+2, j) = moments.Get(2, j);
	  fiphij.Elem(6+4*(i-3)+3, j) = moments.Get(3, j);
	  fiphij.Elem(6+4*(i-3)+4, j) = moments.Get(4, j);
	}
    }

  trans.SetSize(nd);
  CalcInverse (fiphij, trans);
}
  







  
  Array<IPDataHDiv*> FE_BDMPrism1p::ipdata;
  DenseMatrix FE_BDMPrism1p::trans;

  FE_BDMPrism1p :: FE_BDMPrism1p()
  {
    if (!ipdata.Size())
      {
	Orthogonalize();
	CalcIPData(ET_PRISM, ipdata);
      }
  }

  FE_BDMPrism1p :: ~FE_BDMPrism1p()
  {
    ;
  }

  void FE_BDMPrism1p :: CalcShape (const IntegrationPoint & ip, 
				 Vector & shape,
				 int comp) const
  {
    Vector shape1;
    CalcShape1 (ip, shape1, comp);
    trans.MultTrans (shape1, shape);
  }


void FE_BDMPrism1p :: CalcShape1 (const IntegrationPoint & ip, 
				  Vector & shape,
				  int comp) const
{
  int i, j;
  shape.SetSize(21);
  shape.SetScalar (0);

  double x = ip.Point()[0];
  double y = ip.Point()[1];
  double z = ip.Point()[2];
    
  int base = 6*(comp-1);
  
  shape.Elem (base+1) = 1;
  shape.Elem (base+2) = x;
  shape.Elem (base+3) = y;
  shape.Elem (base+4) = z;
  shape.Elem (base+5) = x*z;
  shape.Elem (base+6) = y*z;

  if (comp==3)
    {
      shape.Elem(19) = 1 * z * (1-z);
      shape.Elem(20) = x * z * (1-z);
      shape.Elem(21) = y * z * (1-z);
    }
}


void FE_BDMPrism1p :: Orthogonalize()
{
  cout << "compute BDM prism1+" << endl;
  
  int i, j, k, l;
  int nd = 21;
  
  DenseMatrix fiphij(nd);
  fiphij.SetScalar (0);

  DenseMatrix moments;
  FE_Trig1 trig1;
  FE_Quad1 quad1;

  // 1..6: trig dofs
  for (i = 1; i <= 2; i++)
    {
      ComputeFaceMoments (i, trig1, moments, 5);
      if (i == 1)
	moments.Scale (-1);
      for (j = 1; j <= nd; j++)
	{
	  fiphij.Elem(3*(i-1)+1, j) = moments.Get(1, j);
	  fiphij.Elem(3*(i-1)+2, j) = moments.Get(2, j);
	  fiphij.Elem(3*(i-1)+3, j) = moments.Get(3, j);
	}
    }

  // 7..18 quad dofs
  for (i = 3; i <= 5; i++)
    {
      ComputeFaceMoments (i, quad1, moments, 5);
      for (j = 1; j <= nd; j++)
	{
	  fiphij.Elem(6+4*(i-3)+1, j) = moments.Get(1, j);
	  fiphij.Elem(6+4*(i-3)+2, j) = moments.Get(2, j);
	  fiphij.Elem(6+4*(i-3)+3, j) = moments.Get(3, j);
	  fiphij.Elem(6+4*(i-3)+4, j) = moments.Get(4, j);
	}
    }

  fiphij.Elem(19,19) = 1;
  fiphij.Elem(20,20) = 1;
  fiphij.Elem(21,21) = 1;

  trans.SetSize(nd);
  CalcInverse (fiphij, trans);
}
  












  
  Array<IPDataHDiv*> FE_BDFMPrism2::ipdata;
  DenseMatrix FE_BDFMPrism2::trans;

  FE_BDFMPrism2 :: FE_BDFMPrism2()
  {
    if (!ipdata.Size())
      {
	Orthogonalize();
	CalcIPData(ET_PRISM, ipdata);
      }
  }

  FE_BDFMPrism2 :: ~FE_BDFMPrism2()
  {
    ;
  }

  void FE_BDFMPrism2 :: CalcShape (const IntegrationPoint & ip, 
				 Vector & shape,
				 int comp) const
  {
    Vector shape1;
    CalcShape1 (ip, shape1, comp);
    trans.MultTrans (shape1, shape);
  }


void FE_BDFMPrism2 :: CalcShape1 (const IntegrationPoint & ip, 
				  Vector & shape,
				  int comp) const
{
  int i, j;
  shape.SetSize(27);
  shape.SetScalar (0);

  double x = ip.Point()[0];
  double y = ip.Point()[1];
  double z = ip.Point()[2];
    
  FE_BDFMTrig2 bdfmtrig2;
  FE_Trig1 trig1;
  FE_Segm1 segm1;
  FE_Segm2 segm2;

  IntegrationPoint ipxy(x,y,0,0);
  IntegrationPoint ipz(z,0,0,0);

  Vector vxy1(9), vxy2(3), vz1(2), vz2(3);
  
  if (comp <= 2)
    bdfmtrig2.CalcShape (ipxy, vxy1, comp);
  trig1.CalcShape (ipxy, vxy2);
  segm1.CalcShape (ipz, vz1);
  segm2.CalcShape (ipz, vz2);

  if (comp <= 2)
    for (i = 1; i <= vxy1.Size(); i++)
      for (j = 1; j <= vz1.Size(); j++)
	shape.Elem(i+(j-1)*9) = vxy1.Get(i) * vz1.Get(j);
  else
    for (i = 1; i <= vxy2.Size(); i++)
      for (j = 1; j <= vz2.Size(); j++)
	shape.Elem(18+i+(j-1)*3) = vxy2.Get(i) * vz2.Get(j);
}


void FE_BDFMPrism2 :: Orthogonalize()
{
  cout << "compute BDFM2 prism" << endl;
  
  int i, j, k, l;
  int nd = 27;
  
  DenseMatrix fiphij(nd);
  fiphij.SetScalar (0);

  DenseMatrix moments;
  FE_Trig1 trig1;
  FE_Quad1 quad1;

  // 1..6: trig dofs
  for (i = 1; i <= 2; i++)
    {
      ComputeFaceMoments (i, trig1, moments, 5);
      if (i == 1)
	moments.Scale (-1);
      for (j = 1; j <= nd; j++)
	{
	  fiphij.Elem(3*(i-1)+1, j) = moments.Get(1, j);
	  fiphij.Elem(3*(i-1)+2, j) = moments.Get(2, j);
	  fiphij.Elem(3*(i-1)+3, j) = moments.Get(3, j);
	}
    }

  // 7..18 quad dofs
  for (i = 3; i <= 5; i++)
    {
      ComputeFaceMoments (i, quad1, moments, 5);
      for (j = 1; j <= nd; j++)
	{
	  fiphij.Elem(6+4*(i-3)+1, j) = moments.Get(1, j);
	  fiphij.Elem(6+4*(i-3)+2, j) = moments.Get(2, j);
	  fiphij.Elem(6+4*(i-3)+3, j) = moments.Get(3, j);
	  fiphij.Elem(6+4*(i-3)+4, j) = moments.Get(4, j);
	}
    }


  const POINT3D * points = MeshAccess::ME_GetVertices (ET_PRISM);
  const IntegrationRule & prismrule = 
    GetIntegrationRules().SelectIntegrationRule (ET_PRISM, 5);

  Vector shape1(nd), sum1(nd), sum2(nd), sum3(nd);

  // test with nedelec function
  double p1[3], p2[3], p3[3], p0[3], tau1[3], tau2[3];
  for (j = 0; j <= 2; j++)
    {
      p1[j] = points[0][j];
      p2[j] = points[1][j];
      p3[j] = points[2][j];
    }  
  
  for (k = 1; k <= 3; k++) 
    {  
      for (j = 0; j <= 2; j++)
	{
	  switch (k)
	    {
	    case 1:
	      p0[j] = p1[j];
	      tau1[j] = p2[j] - p1[j];
	      tau2[j] = p3[j] - p1[j];
	      break;
	    case 2:
	      p0[j] = p2[j];
	      tau1[j] = p3[j] - p2[j];
	      tau2[j] = p1[j] - p2[j];
	      break;
	    case 3:
	      p0[j] = p3[j];
	      tau1[j] = p1[j] - p3[j];
	      tau2[j] = p2[j] - p3[j];
	      break;
	    }
	}
      
      sum1.SetScalar (0);
      sum2.SetScalar (0);
      sum3.SetScalar (0);
      
      for (l = 1; l <= prismrule.GetNIP(); l++)
	{ // integration points
	  
	  const IntegrationPoint & tip = prismrule.GetIP(l);
	  
	  double p[3];
	  for (j = 0; j <= 1; j++)
	    {
	      p[j] = p0[j] + 
		tip.Point()[0] * tau1[j] + 
		tip.Point()[1] * tau2[j];
	    }
	  p[2] = tip.Point()[2];

	  IntegrationPoint ip(p, 1.0/6.0);
	  
	  double rt0[2];
	  rt0[0] = -tip.Point()[1];
	  rt0[1] = tip.Point()[0];
	  
	  for (j = 1; j <= 2; j++)
	    {
	      CalcShape1 (ip, shape1, j);

	      sum1.Add ((rt0[0] * tau1[j-1] + rt0[1] * tau2[j-1]) * 
			p[2] * tip.Weight(), shape1);
	      sum2.Add ((rt0[0] * tau1[j-1] + rt0[1] * tau2[j-1]) * 
			(1-p[2]) * tip.Weight(), shape1);
	    }
	  CalcShape1 (ip, shape1, 3);
	  sum3.Add (tip.Point()[0] * tip.Weight(), shape1);
	}
      
      (*testout) << "sum1, 2, 3 = " << endl
		 << sum1 << endl
		 << sum2 << endl
		 << sum3 << endl;

      for (j = 1; j <= nd; j++)
	{
	  fiphij.Elem(18+3*(k-1)+1, j) = sum1.Get(j);	      
	  fiphij.Elem(18+3*(k-1)+2, j) = sum2.Get(j);	      
	  fiphij.Elem(18+3*(k-1)+3, j) = sum3.Get(j);	      
	}
    }

  trans.SetSize(nd);
  CalcInverse (fiphij, trans);

  (*testout) << "fiphij =  " << endl << fiphij << endl;
  (*testout) << "trans = " << endl << trans << endl;
}
  

#endif
#endif











}


