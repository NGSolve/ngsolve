/*********************************************************************/
/* File:   scalarfe.cpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


/* 
   Finite Element Definitions
*/




#include <fem.hpp>

namespace ngfem
{
  
  using namespace ngfem;




  template <int D>
  ScalarFiniteElement<D> :: ~ScalarFiniteElement ()
  {
    // delete block;
    ;
  }

  template <int D>
  void ScalarFiniteElement<D> ::
  CalcDShape (const IntegrationPoint & ip, 
	      FlatMatrix<> dshape) const
  {
    int nd = GetNDof();
    int sdim = SpatialDim();

    double eps = 2e-5;
    ArrayMem<double, 100> hm1(nd), hm2(nd), hm3(nd), hm4(nd);
    FlatVector<> 
      shape1(nd, &hm1[0]), 
      shape2(nd, &hm2[0]), 
      shape3(nd, &hm3[0]), 
      shape4(nd, &hm4[0]);

    for (int i = 0; i < sdim; i++)
      {
	IntegrationPoint ip1 = ip;
	IntegrationPoint ip2 = ip;
	((double*) ip1.Point()) [i] -= eps;
	((double*) ip2.Point()) [i] += eps;
	CalcShape (ip1, shape1);
	CalcShape (ip2, shape2);

	((double*) ip1.Point()) [i] -= eps;
	((double*) ip2.Point()) [i] += eps;
	CalcShape (ip1, shape3);
	CalcShape (ip2, shape4);

	for (int j = 0; j < nd; j++)
	  dshape(j, i) = 
	    2/(3*eps) * (shape2(j) - shape1(j)) 
	    -1/(12*eps) * (shape4(j) - shape3(j));
      }
  }

  /*
    a ( eps - (-eps) ) + b ( 2 eps - (-2eps) )  = 1
    a ( eps^3 - (-eps)^3) + b ( 8 eps^3 - -8 eps^3 ) = 0
    
    2 a + 4 b = 1 / eps
    2 a + 16 b = 0  

    b = -1 / 12 eps
    a = 2 / 3 eps
  */


  /// compute dshape, matrix: ndof x spacedim
  template<int D>
  void ScalarFiniteElement<D> :: CalcMappedDShape (const BaseSpecificIntegrationPoint & bsip, 
                                               FlatMatrix<> dshape) const
  {
    const SpecificIntegrationPoint<D,D> & sip = 
      static_cast<const SpecificIntegrationPoint<D,D> &> (bsip);
    Vec<D> hv;
          
    CalcDShape (sip.IP(), dshape);
    for (int i = 0; i < dshape.Height(); i++)
      {
        hv = dshape.Row(i);
        dshape.Row(i) = Trans (sip.GetJacobianInverse ()) * hv;
      }
  }


  template<int D>
  void ScalarFiniteElement<D> :: CalcDDShape (const IntegrationPoint & ip, 
					  FlatMatrix<> ddshape) const
  {
    int i, j, k;
    int nd = GetNDof();
    int sdim = SpatialDim();

    double eps = 1e-7;
    Matrix<> dshape1(nd, sdim), dshape2(nd, sdim);

    for (i = 0; i < sdim; i++)
      {
	IntegrationPoint ip1 = ip;
	IntegrationPoint ip2 = ip;
	((double*) ip1.Point()) [i] -= eps;
	((double*) ip2.Point()) [i] += eps;

	CalcDShape (ip1, dshape1);
	CalcDShape (ip2, dshape2);
	dshape2 -= dshape1;
	dshape2 *= (0.5 / eps);
	for (j = 0; j < nd; j++)
	  for (k = 0; k < sdim; k++)
	    ddshape(j,sdim*i+k) = dshape2(j,k);
      }  
  }




  template<int D>
  void ScalarFiniteElement<D> :: 
  CalcIPData (ELEMENT_TYPE et,
	      IPDataArray & ipdata)
  {
    if (!ipdata.data.Size())
      {
	const Array<IntegrationPoint*> & ipts = 
	  GetIntegrationRules().GetIntegrationPoints (et);

	/*	
          (*testout) << "Calc IP Data for " << typeid(*this).name() 
          << ", dim = " << dimspace
          << ", type = " << ElementTopology::GetElementName (eltype)
          << ", ndof = " << ndof
          << ", ipts = " << ipts.Size() << endl;
	*/

	ipdata.data.SetSize (ipts.Size());
	// block = new DynamicMem<double> (ipts.Size() * ndof * (1+dimspace));
	ipdata.block.Alloc(ipts.Size() * ndof * (1+dimspace));
	ipdata.block.SetName ("FiniteElement Ipdata");
	double * hp = ipdata.block.Ptr();
	for (int i = 0; i < ipts.Size(); i++)
	  {
	    ipdata.data[i].shape.AssignMemory (ndof, hp);  
	    hp += ndof;
	    ipdata.data[i].dshape.AssignMemory (ndof, dimspace, hp);
	    hp += ndof*dimspace;

	    CalcShape (*ipts[i], ipdata.data[i].shape);
	    CalcDShape (*ipts[i], ipdata.data[i].dshape);
	  }
      }
    
    p_ipdata = &ipdata.data[0];
  }
  
  template<int D>
  const IntegrationRule &
  ScalarFiniteElement<D> :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), Order());
  }




  template<int D>
  void ScalarFiniteElement<D> ::
  EvaluateShapeGrid (const IntegrationRuleTP<D> & ir,
		     const FlatVector<double> coefs,
		     FlatVector<double> gridvalues,
		     LocalHeap & lh) const
  {
    gridvalues = 0;
    void * heapp = lh.GetPointer();
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	gridvalues(i) = InnerProduct (coefs, GetShape(ir[i], lh));
	lh.CleanUp (heapp);
      }
    
  }
		
  template<int D>		  
  void ScalarFiniteElement<D> ::
  EvaluateShapeGridTrans (const IntegrationRuleTP<D> & ir,
			  const FlatVector<double> gridvalues,
			  FlatVector<double> coefs,
			  LocalHeap & lh) const
  {
    coefs = 0.0;
    void * heapp = lh.GetPointer();
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	coefs += gridvalues(i) * GetShape(ir[i], lh);
	lh.CleanUp (heapp);
      }
  }
				  


  template<int D>
  void ScalarFiniteElement<D> ::
  EvaluateDShapeGrid (const IntegrationRuleTP<D> & ir,
		      const FlatVector<double> coefs,
		      FlatMatrixFixWidth<D> gridvalues,
		      LocalHeap & lh) const
  {
    void * heapp = lh.GetPointer();
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	Vec<D> v = Trans (GetDShape(ir[i], lh)) * coefs;
        gridvalues.Row(i) = v;

	lh.CleanUp (heapp);
      }
  }
		
  template<int D>		  
  void ScalarFiniteElement<D> ::
  EvaluateDShapeGridTrans (const IntegrationRuleTP<D> & ir,
			   const FlatMatrixFixWidth<D> gridvalues,
			   FlatVector<double> coefs,
			   LocalHeap & lh) const
  {
    coefs = 0.0;
    FlatVector<> v(SpatialDim(), lh);
    void * heapp = lh.GetPointer();
    for (int i = 0; i < ir.GetNIP(); i++)
      {
	for (int j = 0; j < SpatialDim(); j++)
	  v(j) = gridvalues(i,j);
	coefs += GetDShape(ir[i], lh) * v;
	lh.CleanUp (heapp);
      }
  }
				  








  /* Specific finite elements: */
  

  // Array<T_ScalarFiniteElement<FE_Tet0,3,1>::IPData> FE_Tet0::ipdata;
  ScalarFiniteElement<3>::IPDataArray FE_Tet0::ipdata;

  const IntegrationRule &
  FE_Tet0 :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  }



  // Array<T_ScalarFiniteElement<FE_Tet1,3,4>::IPData> FE_Tet1::ipdata;
  ScalarFiniteElement<3>::IPDataArray FE_Tet1::ipdata;

  const IntegrationRule &
  FE_Tet1 :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  }

  void FE_Tet1 :: GetDofs (Array<Dof> & dofs) const
  {
    dofs.SetSize (0);
    for (int i = 0; i < 4; i++)
      dofs.Append (Dof (Node (NT_VERTEX, i), 0));
  }


  // Array<T_ScalarFiniteElement<FE_Tet2,3,10>::IPData> FE_Tet2::ipdata;
  ScalarFiniteElement<3>::IPDataArray FE_Tet2::ipdata;






  ScalarFiniteElement<3>::IPDataArray FE_Tet2HB::ipdata;

  FE_Tet2HB :: FE_Tet2HB()
    : ScalarFiniteElement<3> (ET_TET, 10, 2)

  {
    if (!ipdata.data.Size())
      {
	CalcIPData(ET_TET, ipdata);
      }
  }

  FE_Tet2HB :: ~FE_Tet2HB()
  {
    ;
  }

  void FE_Tet2HB :: CalcShape (const IntegrationPoint & ip, 
			       FlatVector<> shape) const
			     
  {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);
    double lam4 = 1 - x - y - z;

    shape(0) = x;
    shape(1) = y;
    shape(2) = z;
    shape(3) = lam4;

    shape(4) = 4 * x * y;
    shape(5) = 4 * x * z;
    shape(6) = 4 * x * lam4;
    shape(7) = 4 * y * z;
    shape(8) = 4 * y * lam4;
    shape(9) = 4 * z * lam4;
  }




  ScalarFiniteElement<2>::IPDataArray FE_Trig0::ipdata;

  FE_Trig0 :: FE_Trig0()
    : ScalarFiniteElement<2> (ET_TRIG, 1, 0)
  {
    if (!ipdata.Size())
      CalcIPData(ET_TRIG, ipdata);
  }

  FE_Trig0 :: ~FE_Trig0()
  {
    ;
  }

  void FE_Trig0 :: CalcShape (const IntegrationPoint & ip, 
			      FlatVector<> shape) const
			     
  {
    shape(0) = 1;
  }

  void FE_Trig0 :: CalcDShape (const IntegrationPoint & ip, 
                               FlatMatrix<> dshape) const
			      
  {
    dshape = 0;
  }

  const IntegrationRule &
  FE_Trig0 :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  }




  // Array<T_ScalarFiniteElement<FE_Trig1,2,3>::IPData> FE_Trig1::ipdata;
  ScalarFiniteElement<2>::IPDataArray FE_Trig1::ipdata;

  const IntegrationRule &
  FE_Trig1 :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  } 



  // Array<T_ScalarFiniteElement<FE_Trig2,2,6>::IPData> FE_Trig2::ipdata;
  ScalarFiniteElement<2>::IPDataArray FE_Trig2::ipdata;

  const IntegrationRule &
  FE_Trig2 :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 2);
  } 


  // Array<ScalarFiniteElement::IPData> FE_Trig2HB::ipdata;
  ScalarFiniteElement<2>::IPDataArray FE_Trig2HB::ipdata;

  FE_Trig2HB :: FE_Trig2HB()
    : ScalarFiniteElement<2> (ET_TRIG, 6, 2)
  {
    if (!ipdata.Size())
      CalcIPData(ET_TRIG, ipdata);
  }

  FE_Trig2HB :: ~FE_Trig2HB()
  {
    ;
  }

  void FE_Trig2HB :: CalcShape (const IntegrationPoint & ip, 
				FlatVector<> shape) const
			     
  {
    double x = ip(0);
    double y = ip(1);
    double lam3 = 1-x-y;

    shape(0) = x;
    shape(1) = y;
    shape(2) = lam3;
    shape(3) = 4 * y * lam3;
    shape(4) = 4 * x * lam3;
    shape(5) = 4 * x * y;
  }











  // Array<ScalarFiniteElement::IPData> FE_Trig3Pot::ipdata;
  ScalarFiniteElement<2>::IPDataArray FE_Trig3Pot::ipdata;

  FE_Trig3Pot :: FE_Trig3Pot()
    : ScalarFiniteElement<2> (ET_TRIG, 10, 3)
  {
    if (!ipdata.Size())
      CalcIPData(ET_TRIG, ipdata);
  }

  FE_Trig3Pot :: ~FE_Trig3Pot()
  {
    ;
  }

  void FE_Trig3Pot :: CalcShape (const IntegrationPoint & ip, 
				 FlatVector<> shape) const
				 
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


  // Array<T_ScalarFiniteElement<FE_NC_Trig1,2,3>::IPData> FE_NC_Trig1::ipdata;
  ScalarFiniteElement<2>::IPDataArray FE_NC_Trig1::ipdata;

  const IntegrationRule &
  FE_NC_Trig1 :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  } 









  // Array<ScalarFiniteElement::IPData> FE_Tet3Pot::ipdata;
  ScalarFiniteElement<3>::IPDataArray FE_Tet3Pot::ipdata;

  FE_Tet3Pot :: FE_Tet3Pot()
    : ScalarFiniteElement<3> (ET_TET, 20, 3)
  {
    if (!ipdata.Size())
      CalcIPData(ET_TET, ipdata);
  }

  FE_Tet3Pot :: ~FE_Tet3Pot()
  {
    ;
  }

  void FE_Tet3Pot :: CalcShape (const IntegrationPoint & ip, 
				FlatVector<> shape) const
				 
  {
    int i;
    double lami[4];

    lami[0] = ip.Point()[0];
    lami[1] = ip.Point()[1];
    lami[2] = ip.Point()[2];
    lami[3] = 1-lami[0]-lami[1]-lami[2];


    for (i = 0; i < 4; i++)
      shape(i) = lami[i];

    const int edges[6][2] = 
      { { 4, 1 },
        { 4, 2 },
        { 4, 3 }, 
        { 1, 2 },
        { 1, 3 },
        { 2, 3 }};
    const int faces[4][3] =
      { { 4, 2, 3 },
        { 4, 1, 3 },
        { 4, 1, 2 },
        { 1, 2, 3 } };


    for (i = 0; i < 6; i++)
      {
	double l1 = lami[edges[i][0]-1];
	double l2 = lami[edges[i][1]-1];
	shape(2*i+4) = 3 * l1 * l2 * (l1+l2);
	shape(2*i+5) = 7.5 * l1 * l2 * (l2-l1);
      }
    
    for (i = 0; i < 4; i++)
      {
	shape(i+16) = 
	  60 * lami[faces[i][0]-1] * lami[faces[i][1]-1] * lami[faces[i][2]-1];
      }
  }











  // Array<ScalarFiniteElement::IPData> FE_Quad0::ipdata;
  ScalarFiniteElement<2>::IPDataArray FE_Quad0::ipdata;

  FE_Quad0 :: FE_Quad0()
    : ScalarFiniteElement<2> (ET_QUAD, 1, 0)
  {
    if (!ipdata.Size())
      CalcIPData(ET_QUAD, ipdata);
  }

  FE_Quad0 :: ~FE_Quad0()
  {
    ;
  }

  void FE_Quad0 :: CalcShape (const IntegrationPoint & ip, 
			      FlatVector<> shape) const
			     
  {
    shape(0) = 1;
  }

  void FE_Quad0 :: CalcDShape (const IntegrationPoint & ip, 
			       FlatMatrix<> dshape) const
			      
  {
    dshape = 0;
  }

  const IntegrationRule &
  FE_Quad0 :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  }




  // Array<ScalarFiniteElement::IPData> FE_Quad1::ipdata;
  ScalarFiniteElement<2>::IPDataArray FE_Quad1::ipdata;

  FE_Quad1 :: FE_Quad1()
    : ScalarFiniteElement<2> (ET_QUAD, 4, 1)
  {
    if (!ipdata.Size())
      CalcIPData(ET_QUAD, ipdata);
  }

  FE_Quad1 :: ~FE_Quad1()
  {
    ;
  }

  void FE_Quad1 :: CalcShape (const IntegrationPoint & ip, 
			      FlatVector<> shape) const
			      
  {
    double x = ip.Point()[0];
    double y = ip.Point()[1];

    shape(0) = (1-x) * (1-y);
    shape(1) =    x  * (1-y);
    shape(2) =    x  *  y;
    shape(3) = (1-x) *  y;
  }

  void FE_Quad1 :: CalcDShape (const IntegrationPoint & ip, 
			       FlatMatrix<> dshape) const
			      
  {
    double x = ip.Point()[0];
    double y = ip.Point()[1];

    dshape(0,0) = -(1-y);
    dshape(0,1) = -(1-x);

    dshape(1,0) = 1-y;
    dshape(1,1) = -x;

    dshape(2,0) = y;
    dshape(2,1) = x;

    dshape(3,0) = -y;
    dshape(3,1) = (1-x);
  }


  const IntegrationRule &
  FE_Quad1:: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  }






  // Array<ScalarFiniteElement::IPData> FE_Quad2::ipdata;
  ScalarFiniteElement<2>::IPDataArray FE_Quad2::ipdata;

  FE_Quad2 :: FE_Quad2()
    : ScalarFiniteElement<2> (ET_QUAD, 9, 2)
  {
    if (!ipdata.Size())
      CalcIPData(ET_QUAD, ipdata);
  }

  FE_Quad2 :: ~FE_Quad2()
  {
    ;
  }

  void FE_Quad2 :: CalcShape (const IntegrationPoint & ip, 
			      FlatVector<> shape) const
			      
  {
    double x = ip.Point()[0];
    double y = ip.Point()[1];

    Vec<3> px, py;
    px(0) = (1-x) * (1-2*x);
    px(1) = 4 * x * (1-x);
    px(2) = x * (2*x-1);
    py(0) = (1-y) * (1-2*y);
    py(1) = 4 * y * (1-y);
    py(2) = y * (2*y-1);

    int ii = 0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
	shape(ii++) = px(i) * py(j);
  }

  void FE_Quad2 :: CalcDShape (const IntegrationPoint & ip, 
			       FlatMatrix<> dshape) const
    
  {
    double x = ip.Point()[0];
    double y = ip.Point()[1];

    Vec<3> px, py, dpx, dpy;
    px(0) = (1-x) * (1-2*x);
    px(1) = 4 * x * (1-x);
    px(2) = x * (2*x-1);

    py(0) = (1-y) * (1-2*y);
    py(1) = 4 * y * (1-y);
    py(2) = y * (2*y-1);

    dpx(0) = 4*x-3;
    dpx(1) = 8*x-4;
    dpx(2) = 4*x-1;

    dpy(0) = 4*y-3;
    dpy(1) = 8*y-4;
    dpy(2) = 4*y-1;
    
    int ii = 0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
	{
	  dshape(ii,0) = dpx(i) * py(j);
	  dshape(ii,1) = px(i) * dpy(j);
	  ii++;
	}
  }


  const IntegrationRule &
  FE_Quad2:: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  }








  // Array<ScalarFiniteElement::IPData> FE_Quad3::ipdata;
  ScalarFiniteElement<2>::IPDataArray FE_Quad3::ipdata;

  FE_Quad3 :: FE_Quad3()
    : ScalarFiniteElement<2> (ET_QUAD, 16, 3)
  {
    if (!ipdata.Size())
      CalcIPData(ET_QUAD, ipdata);
  }

  FE_Quad3 :: ~FE_Quad3()
  {
    ;
  }

  void FE_Quad3 :: CalcShape (const IntegrationPoint & ip, 
			      FlatVector<> shape) const
    
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
			       FlatMatrix<> dshape) const
    
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












  // Array<ScalarFiniteElement::IPData> FE_Quad2aniso::ipdata;
  ScalarFiniteElement<2>::IPDataArray FE_Quad2aniso::ipdata;

  FE_Quad2aniso :: FE_Quad2aniso()
    : ScalarFiniteElement<2> (ET_QUAD, 6, 2)
  {
    if (!ipdata.Size())
      CalcIPData(ET_QUAD, ipdata);
  }

  FE_Quad2aniso :: ~FE_Quad2aniso()
  {
    ;
  }

  void FE_Quad2aniso :: CalcShape (const IntegrationPoint & ip, 
				   FlatVector<> shape) const
			      
  {
    double x = ip.Point()[0];
    double y = ip.Point()[1];

    shape(0) = (1-x)*(1-2*x) * (1-y);
    shape(1) = x*(2*x-1) * (1-y);
    shape(2) = x*(2*x-1) * y;
    shape(3) = (1-x)*(1-2*x) * y;
    shape(4) = 4*x*(1-x) * (1-y);
    shape(5) = 4*x*(1-x) * y;
  }













  // Array<ScalarFiniteElement::IPData> FE_Pyramid0::ipdata;
  ScalarFiniteElement<3>::IPDataArray FE_Pyramid0::ipdata;

  FE_Pyramid0 :: FE_Pyramid0()
    : ScalarFiniteElement<3> (ET_PYRAMID, 1, 0)
  {
    if (!ipdata.Size())
      CalcIPData(ET_PYRAMID, ipdata);
  }

  FE_Pyramid0 :: ~FE_Pyramid0()
  {
    ;
  }

  void FE_Pyramid0 :: CalcShape (const IntegrationPoint & ip, 
				 FlatVector<> shape) const
			     
  {
    shape(0) = 1;
  }

  void FE_Pyramid0 :: CalcDShape (const IntegrationPoint & ip, 
				  FlatMatrix<> dshape) const
				  
  {
    dshape = 0;
  }

  const IntegrationRule &
  FE_Pyramid0 :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  }









  // Array<ScalarFiniteElement::IPData> FE_Pyramid1::ipdata;
  ScalarFiniteElement<3>::IPDataArray FE_Pyramid1::ipdata;

  FE_Pyramid1 :: FE_Pyramid1()
    : ScalarFiniteElement<3> (ET_PYRAMID, 5, 1)
  {
    if (!ipdata.Size())
      CalcIPData(ET_PYRAMID, ipdata);
  }

  FE_Pyramid1 :: ~FE_Pyramid1()
  {
    ;
  }

  void FE_Pyramid1 :: CalcShape (const IntegrationPoint & ip, 
				 FlatVector<> shape) const
			     
  {
    double x = ip.Point()[0];
    double y = ip.Point()[1];
    double z = ip.Point()[2];
    if (z == 1) z = 1-1e-10;
    shape(0) = (1-z-x)*(1-z-y) / (1-z);
    shape(1) = x*(1-z-y) / (1-z);
    shape(2) = x*y / (1-z);
    shape(3) = (1-z-x)*y / (1-z);
    shape(4) = z;
  }


  void FE_Pyramid1 :: CalcDShape (const IntegrationPoint & ip, 
				  FlatMatrix<> dshape) const
			     
  {
    double x = ip.Point()[0];
    double y = ip.Point()[1];
    double z = ip.Point()[2];
    if (z == 1) z = 1-1e-10;

    // shape(0) = (1-z-x)*(1-z-y) / (1-z);
    dshape(0,0) = -(1-z-y) / (1-z);
    dshape(0,1) = -(1-z-x) / (1-z);
    dshape(0,2) = x*y / (1-z) / (1-z) - 1;

    // shape(1) = x*(1-z-y) / (1-z);
    dshape(1,0) = (1-z-y) / (1-z);
    dshape(1,1) = -x / (1-z);
    dshape(1,2) = -x*y / (1-z) / (1-z);
    
    // shape(2) = x*y / (1-z);
    dshape(2,0) = y / (1-z);
    dshape(2,1) = x / (1-z);
    dshape(2,2) = x * y / (1-z) / (1-z);
    
    // shape(3) = (1-z-x)*y / (1-z);
    dshape(3,0) = -y / (1-z);
    dshape(3,1) = (1-z-x) / (1-z);
    dshape(3,2) = -x*y / (1-z) / (1-z);

    // shape(4) = z;
    dshape(4,0) = 0;
    dshape(4,1) = 0;
    dshape(4,2) = 1;
  }





  const IntegrationRule &
  FE_Pyramid1 :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  }








  // Array<ScalarFiniteElement::IPData> FE_Pyramid2::ipdata;
  ScalarFiniteElement<3>::IPDataArray FE_Pyramid2::ipdata;

  FE_Pyramid2 :: FE_Pyramid2() : ScalarFiniteElement<3> (ET_PYRAMID)
  {
    if (!ipdata.Size())
      CalcIPData(ET_PYRAMID, ipdata);
  }

  FE_Pyramid2 :: ~FE_Pyramid2()
  {
    ;
  }

  void FE_Pyramid2 :: CalcShape (const IntegrationPoint & ip, 
				 FlatVector<> shape) const
			     
  {
    double x = ip.Point()[0];
    double y = ip.Point()[1];
    double z = ip.Point()[2];

    shape = 0;
    shape(0) = (1-z-x)*(1-z-y) / (1-z);
    shape(1) = x*(1-z-y) / (1-z);
    shape(2) = x*y / (1-z);
    shape(3) = (1-z-x)*y / (1-z);
    shape(4) = z;
  }











  // Array<ScalarFiniteElement::IPData> FE_Prism0::ipdata;
  ScalarFiniteElement<3>::IPDataArray FE_Prism0::ipdata;

  FE_Prism0 :: FE_Prism0()
    : ScalarFiniteElement<3> (ET_PRISM, 1, 0)
  {
    if (!ipdata.Size())
      CalcIPData(ET_PRISM, ipdata);
  }

  FE_Prism0 :: ~FE_Prism0()
  {
    ;
  }

  void FE_Prism0 :: CalcShape (const IntegrationPoint & ip, 
			       FlatVector<> shape) const
			     
  {
    shape(0) = 1;
  }

  void FE_Prism0 :: CalcDShape (const IntegrationPoint & ip, 
				FlatMatrix<> dshape) const
			      
  {
    dshape = 0;
  }

  const IntegrationRule &
  FE_Prism0 :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  }




  // Array<ScalarFiniteElement::IPData> FE_Prism1::ipdata;
  ScalarFiniteElement<3>::IPDataArray FE_Prism1::ipdata;

  FE_Prism1 :: FE_Prism1()
    : ScalarFiniteElement<3> (ET_PRISM, 6, 1)
  {
    if (!ipdata.Size())
      CalcIPData(ET_PRISM, ipdata);
  }

  FE_Prism1 :: ~FE_Prism1()
  {
    ;
  }

  void FE_Prism1 :: CalcShape (const IntegrationPoint & ip, 
			       FlatVector<> shape) const
    
  {
    double x = ip.Point()[0];
    double y = ip.Point()[1];
    double z = ip.Point()[2];

    shape(0) = x * (1-z);
    shape(1) = y * (1-z);
    shape(2) = (1-x-y) * (1-z);
    shape(3) = x * z;
    shape(4) = y * z;
    shape(5) = (1-x-y) * z;
    
  }

  void FE_Prism1 :: CalcDShape (const IntegrationPoint & ip, 
				FlatMatrix<> dshape) const
    
  {
    double x = ip.Point()[0];
    double y = ip.Point()[1];
    double z = ip.Point()[2];

    // shape(0) = x * (1-z);
    dshape(0,0) = 1-z;
    dshape(0,1) = 0;
    dshape(0,2) = -x;

    // shape(1) = y * (1-z);
    dshape(1,0) = 0;
    dshape(1,1) = 1-z;
    dshape(1,2) = -y;

    // shape(2) = (1-x-y) * (1-z);
    dshape(2,0) = -(1-z);
    dshape(2,1) = -(1-z);
    dshape(2,2) = -(1-x-y);

    // shape(3) = x * z;
    dshape(3,0) = z;
    dshape(3,1) = 0;
    dshape(3,2) = x;

    // shape(4) = y * z;
    dshape(4,0) = 0;
    dshape(4,1) = z;
    dshape(4,2) = y;

    // shape(5) = (1-x-y) * z;
    dshape(5,0) = -z;
    dshape(5,1) = -z;
    dshape(5,2) = 1-x-y;
  }


  const IntegrationRule &
  FE_Prism1 :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  }






  // Array<ScalarFiniteElement::IPData> FE_Prism2::ipdata;
  ScalarFiniteElement<3>::IPDataArray FE_Prism2::ipdata;

  FE_Prism2 :: FE_Prism2()
    : ScalarFiniteElement<3> (ET_PRISM, 18, 3)
  {
    if (!ipdata.Size())
      CalcIPData(ET_PRISM, ipdata);
  }

  FE_Prism2 :: ~FE_Prism2()
  {
    ;
  }

  void FE_Prism2 :: CalcShape (const IntegrationPoint & ip, 
			       FlatVector<> shape) const
			     
  {
    double x = ip.Point()[0];
    double y = ip.Point()[1];
    double z = ip.Point()[2];

    shape(0) = x * (1-z);
    shape(1) = y * (1-z);
    shape(2) = (1-x-y) * (1-z);
    shape(3) = x * z;
    shape(4) = y * z;
    shape(5) = (1-x-y) * z;

    shape(6) = 4 * x * (1-x-y) * (1-z);
    shape(7) = 4 * x * y       * (1-z);
    shape(8) = 4 * y * (1-x-y) * (1-z);
    shape(9) = 4 * x * (1-x-y) * z;
    shape(10) = 4 * x * y       * z;
    shape(11) = 4 * y * (1-x-y) * z;


    shape(12) = x * (1-z) * z;
    shape(13) = y * (1-z) * z;
    shape(14) = (1-x-y) * (1-z) * z;
    shape(15) = 4 * x * (1-x-y) * (1-z) * z;
    shape(16) = 4 * x * y       * (1-z) * z;
    shape(17) = 4 * y * (1-x-y) * (1-z) * z;
  }



  void FE_Prism2 :: CalcDShape (const IntegrationPoint & ip, 
				FlatMatrix<> dshape) const
  {
    ScalarFiniteElement<3>::CalcDShape (ip, dshape);
  }




  // Array<ScalarFiniteElement::IPData> FE_Prism2aniso::ipdata;
  ScalarFiniteElement<3>::IPDataArray FE_Prism2aniso::ipdata;

  FE_Prism2aniso :: FE_Prism2aniso()
    : ScalarFiniteElement<3> (ET_PRISM, 12, 3)
  {
    if (!ipdata.Size())
      CalcIPData(ET_PRISM, ipdata);
  }

  FE_Prism2aniso :: ~FE_Prism2aniso()
  {
    ;
  }

  void FE_Prism2aniso :: CalcShape (const IntegrationPoint & ip, 
				    FlatVector<> shape) const
			     
  {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);
    double lam3 = 1-x-y;

    shape(0) = x * (2*x-1) * (1-z);
    shape(1) = y * (2*y-1) * (1-z);
    shape(2) = lam3 * (2*lam3-1) * (1-z);
    shape(3) = x * (2*x-1) * z;
    shape(4) = y * (2*y-1) * z;
    shape(5) = lam3 * (2*lam3-1) * z;

    shape(6) = 4 * x * lam3 * (1-z);
    shape(7) = 4 * x * y       * (1-z);
    shape(8) = 4 * y * lam3 * (1-z);
    shape(9) = 4 * x * lam3 * z;
    shape(10) = 4 * x * y       * z;
    shape(11) = 4 * y * lam3 * z;
  }


  void FE_Prism2aniso :: CalcDShape (const IntegrationPoint & ip, 
				     FlatMatrix<> dshape) const
  {
    ScalarFiniteElement<3>::CalcDShape (ip, dshape);
  }





  // Array<ScalarFiniteElement::IPData> FE_Prism2HBaniso::ipdata;
  ScalarFiniteElement<3>::IPDataArray FE_Prism2HBaniso::ipdata;

  FE_Prism2HBaniso :: FE_Prism2HBaniso()
    : ScalarFiniteElement<3> (ET_PRISM, 12, 3)
  {
    if (!ipdata.Size())
      CalcIPData(ET_PRISM, ipdata);
  }

  FE_Prism2HBaniso :: ~FE_Prism2HBaniso()
  {
    ;
  }

  void FE_Prism2HBaniso :: CalcShape (const IntegrationPoint & ip, 
                                      FlatVector<> shape) const
			     
  {
    double x = ip(0);
    double y = ip(1);
    double z = ip(2);

    shape(0) = x * (1-z);
    shape(1) = y * (1-z);
    shape(2) = (1-x-y) * (1-z);
    shape(3) = x * z;
    shape(4) = y * z;
    shape(5) = (1-x-y) * z;

    shape(6) = 4 * x * (1-x-y) * (1-z);
    shape(7) = 4 * x * y       * (1-z);
    shape(8) = 4 * y * (1-x-y) * (1-z);
    shape(9) = 4 * x * (1-x-y) * z;
    shape(10) = 4 * x * y       * z;
    shape(11) = 4 * y * (1-x-y) * z;
  }









  // Array<ScalarFiniteElement::IPData> FE_Prism3aniso::ipdata;
  ScalarFiniteElement<3>::IPDataArray FE_Prism3aniso::ipdata;

  FE_Prism3aniso :: FE_Prism3aniso() 
    : ScalarFiniteElement<3> (ET_PRISM)
  {
    if (!ipdata.Size())
      CalcIPData(ET_PRISM, ipdata);
  }

  FE_Prism3aniso :: ~FE_Prism3aniso()
  {
    ;
  }

  void FE_Prism3aniso :: CalcShape (const IntegrationPoint & ip, 
				    FlatVector<> shape) const
			     
  {
    double x = ip.Point()[0];
    double y = ip.Point()[1];
    double lam3 = 1-x-y;
    double z = ip.Point()[2];
    int base;
    for (base = 0; base <= 10; base += 10)
      {
	double shapez = (base == 0) ? z : 1-z;

	shape(base  ) = shapez * x;
	shape(base+1) = shapez * y;
	shape(base+2) = shapez * lam3;
	shape(base+3) = shapez * x*x*y;
	shape(base+4) = shapez * x*y*y;
	shape(base+5) = shapez * x*x*lam3;
	shape(base+6) = shapez * x*lam3*lam3;
	shape(base+7) = shapez * y*y*lam3;
	shape(base+8) = shapez * y*lam3*lam3;
	shape(base+9) = shapez * x*y*lam3;
      }
  }















  // Array<ScalarFiniteElement::IPData> FE_Hex1::ipdata;
  ScalarFiniteElement<3>::IPDataArray FE_Hex1::ipdata;

  FE_Hex1 :: FE_Hex1()
    : ScalarFiniteElement<3> (ET_HEX, 8, 1)
  {
    if (!ipdata.Size())
      CalcIPData(ET_HEX, ipdata);
  }

  FE_Hex1 :: ~FE_Hex1()
  {
    ;
  }

  void FE_Hex1 :: CalcShape (const IntegrationPoint & ip, 
			     FlatVector<> shape) const
			     
  {
    double x = ip.Point()[0];
    double y = ip.Point()[1];
    double z = ip.Point()[2];

    shape(0) = (1-x) * (1-y) * (1-z);
    shape(1) =    x  * (1-y) * (1-z);
    shape(2) =    x  *    y  * (1-z);
    shape(3) = (1-x) *    y  * (1-z);
    shape(4) = (1-x) * (1-y) *    z ;
    shape(5) =    x  * (1-y) *    z ;
    shape(6) =    x  *    y  *    z ;
    shape(7) = (1-x) *    y  *    z ;
  }






  // Array<ScalarFiniteElement::IPData> FE_Hex0::ipdata;
  ScalarFiniteElement<3>::IPDataArray FE_Hex0::ipdata;

  FE_Hex0 :: FE_Hex0()
    : ScalarFiniteElement<3> (ET_HEX, 1, 0)
  {
    if (!ipdata.Size())
      CalcIPData(ET_HEX, ipdata);
  }

  FE_Hex0 :: ~FE_Hex0()
  {
    ;
  }

  void FE_Hex0 :: CalcShape (const IntegrationPoint & ip, 
			     FlatVector<> shape) const
			     
  {
    double x = ip.Point()[0];
    double y = ip.Point()[1];
    double z = ip.Point()[2];

    shape(0) = 1;
  }


















  // Array<T_ScalarFiniteElement<FE_Segm0,1,1>::IPData> FE_Segm0::ipdata;
  ScalarFiniteElement<1>::IPDataArray FE_Segm0::ipdata;

  ScalarFiniteElement<1>::IPDataArray FE_SegmDummy::ipdata;

  // Array<T_ScalarFiniteElement<FE_Segm1,1,2>::IPData> FE_Segm1::ipdata;
  ScalarFiniteElement<1>::IPDataArray FE_Segm1::ipdata;

  //  Array<T_ScalarFiniteElement<FE_Segm1L2,1,2>::IPData> FE_Segm1L2::ipdata;
  ScalarFiniteElement<1>::IPDataArray FE_Segm1L2::ipdata;
  FE_Segm1L2 ::  FE_Segm1L2() { ; }

  // Array<T_ScalarFiniteElement<FE_Segm2,1,3>::IPData> FE_Segm2::ipdata;
  ScalarFiniteElement<1>::IPDataArray FE_Segm2::ipdata;

  // Array<T_ScalarFiniteElement<FE_Segm2L2,1,3>::IPData> FE_Segm2L2::ipdata;
  ScalarFiniteElement<1>::IPDataArray FE_Segm2L2::ipdata;
  FE_Segm2L2 :: FE_Segm2L2() { ; }

  // Array<T_ScalarFiniteElement<FE_Segm2HB,1,3>::IPData> FE_Segm2HB::ipdata;
  ScalarFiniteElement<1>::IPDataArray FE_Segm2HB::ipdata;

  // Array<T_ScalarFiniteElement<FE_NcSegm1,1,1>::IPData> FE_NcSegm1::ipdata;
  ScalarFiniteElement<1>::IPDataArray FE_NcSegm1::ipdata;

  template <int ORDER>
  ScalarFiniteElement<1>::IPDataArray FE_TSegmL2<ORDER>::ipdata;
  // Array<ScalarFiniteElement::IPData> FE_TSegmL2<ORDER>::ipdata;

  /*
    Array<ScalarFiniteElement::IPData> FE_TSegmL2<0>::ipdata;
    Array<ScalarFiniteElement::IPData> FE_TSegmL2<1>::ipdata;
    Array<ScalarFiniteElement::IPData> FE_TSegmL2<2>::ipdata;
    Array<ScalarFiniteElement::IPData> FE_TSegmL2<3>::ipdata;
  */
  template <int ORDER>
  FE_TSegmL2<ORDER> :: FE_TSegmL2()
    : ScalarFiniteElement<1> (ET_SEGM, ORDER+1, ORDER)
  {
    if (!ipdata.Size())
      CalcIPData(ET_SEGM, ipdata);
  }

  template <int ORDER>
  FE_TSegmL2<ORDER> :: ~FE_TSegmL2()
  {
    ;
  }

  template <int ORDER>
  void FE_TSegmL2<ORDER> :: CalcShape (const IntegrationPoint & ip, 
				       FlatVector<> shape) const
  {
    double x = ip.Point()[0];
    shape = 0;

    if (ORDER >= 0) shape(0) = 1;
    if (ORDER >= 1) shape(1) = 2*x-1;
    if (ORDER >= 2) shape(2) = (2*x-1)*(2*x-1)-1.0/3.0;
    if (ORDER >= 3) shape(3) = (2*x-1)*(2*x-1)*(2*x-1);
    if (ORDER >= 4)
      {
	throw Exception ("TSegmL2: Legendre polynomials not implemented");
      }
  }

  template <int ORDER>
  void FE_TSegmL2<ORDER> :: CalcDShape (const IntegrationPoint & ip, 
					FlatMatrix<> dshape) const
  {
    double x = ip.Point()[0];
    dshape = 0;

    if (ORDER >= 0) dshape(0) = 0;
    if (ORDER >= 1) dshape(1) = 2;
    if (ORDER >= 2) dshape(2) = 8*x-4;
    if (ORDER >= 3) dshape(3) = 6*(2*x-1)*(2*x-1);
    if (ORDER >= 4)
      {
	throw Exception ("TSegmL2: Legendre polynomials not implemented");
      }
  }

  template class  FE_TSegmL2<0>;
  template class  FE_TSegmL2<1>;
  template class  FE_TSegmL2<2>;
  template class  FE_TSegmL2<3>;








  // Array<ScalarFiniteElement::IPData> FE_Segm3Pot::ipdata;
  ScalarFiniteElement<1>::IPDataArray FE_Segm3Pot::ipdata;

  FE_Segm3Pot :: FE_Segm3Pot()
    : ScalarFiniteElement<1> (ET_SEGM, 4, 3)
  {
    if (!ipdata.Size())
      CalcIPData(ET_SEGM, ipdata);
  }

  FE_Segm3Pot :: ~FE_Segm3Pot()
  {
    ;
  }

  void FE_Segm3Pot :: CalcShape (const IntegrationPoint & ip, 
				 FlatVector<> shape) const
				 
  {
    double x = ip.Point()[0];
    double lam2 = 1-x;

    shape(0) = x;
    shape(1) = lam2;

    shape(2) = 3 * x * lam2 * (lam2+x);
    shape(3) = 7.5 * x * lam2 * (x-lam2);
  }





  










  // Array<ScalarFiniteElement::IPData> FE_NcTrig1::ipdata;
  ScalarFiniteElement<2>::IPDataArray FE_NcTrig1::ipdata;


  FE_NcTrig1 :: FE_NcTrig1() : ScalarFiniteElement<2> (ET_TRIG)
  {
    if (!ipdata.Size())
      CalcIPData(ET_TRIG, ipdata);
  }

  FE_NcTrig1 :: ~FE_NcTrig1()
  {
    ;
  }

  void FE_NcTrig1 :: CalcShape (const IntegrationPoint & ip, 
				FlatVector<> shape) const
			     
  {
    double x = ip.Point()[0];
    double y = ip.Point()[1];

    shape(0) = 1-2*y;
    shape(1) = 1-2*x;
    shape(2) = 2*(x+y)-1;
  }






  // Array<ScalarFiniteElement::IPData> FE_NcTet1::ipdata;
  ScalarFiniteElement<3>::IPDataArray FE_NcTet1::ipdata;
  
  FE_NcTet1 :: FE_NcTet1() 
    : ScalarFiniteElement<3> (ET_TET)
  {
    if (!ipdata.Size())
      CalcIPData(ET_TET, ipdata);
  }

  FE_NcTet1 :: ~FE_NcTet1()
  {
    ;
  }

  void FE_NcTet1 :: CalcShape (const IntegrationPoint & ip, 
			       FlatVector<> shape) const
    
  {
    double x = ip.Point()[0];
    double y = ip.Point()[1];
    double z = ip.Point()[2];

    shape(0) = 1-2*x;
    shape(1) = 1-2*y;
    shape(2) = 1-2*z;
    shape(3) = 2*(x+y+z)-1;
  }


  /* Dummy elements */
  ScalarFiniteElement<2>::IPDataArray FE_TrigDummy::ipdata;
  ScalarFiniteElement<2>::IPDataArray FE_QuadDummy::ipdata;

  ScalarFiniteElement<3>::IPDataArray FE_HexDummy::ipdata;
  ScalarFiniteElement<3>::IPDataArray FE_TetDummy::ipdata;
  ScalarFiniteElement<3>::IPDataArray FE_PrismDummy::ipdata;
  ScalarFiniteElement<3>::IPDataArray FE_PyramidDummy::ipdata;








  template class ScalarFiniteElement<1>;
  template class ScalarFiniteElement<2>;
  template class ScalarFiniteElement<3>;

}

