/*********************************************************************/
/* File:   scalarfe.cpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


/* 
   Finite Element Definitions
*/




#include <fem.hpp>
#include <h1lofe.hpp>

namespace ngfem
{
  
  using namespace ngfem;


  template <int D>
  void ScalarFiniteElement<D> ::
  CalcDShape (const IntegrationPoint & ip, 
	      FlatMatrixFixWidth<D> dshape) const
  
  {
    static bool firsttime = true;
    if (firsttime)
      {
	cout << "WARNING: CalcDShape not overloaded for class, using numerical differentiation " << typeid(this).name() << ", ndof = " << ndof << endl;
        firsttime = false;
      }
    
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
  void ScalarFiniteElement<D> :: 
  CalcMappedDShape (const SpecificIntegrationPoint<D,D> & sip, 
                    FlatMatrixFixWidth<D> dshape) const
  {
    CalcDShape (sip.IP(), dshape);
    for (int i = 0; i < dshape.Height(); i++)
      {
        Vec<D> hv = dshape.Row(i);
        dshape.Row(i) = Trans (sip.GetJacobianInverse ()) * hv;
      }
  }


  template<int D>
  void ScalarFiniteElement<D> :: CalcDDShape (const IntegrationPoint & ip, 
					  FlatMatrix<> ddshape) const
  {
    int nd = GetNDof();
    int sdim = SpatialDim();

    double eps = 1e-7;
    Matrix<> dshape1(nd, sdim), dshape2(nd, sdim);

    for (int i = 0; i < sdim; i++)
      {
	IntegrationPoint ip1 = ip;
	IntegrationPoint ip2 = ip;
	((double*) ip1.Point()) [i] -= eps;
	((double*) ip2.Point()) [i] += eps;

	CalcDShape (ip1, dshape1);
	CalcDShape (ip2, dshape2);
	dshape2 -= dshape1;
	dshape2 *= (0.5 / eps);
	for (int j = 0; j < nd; j++)
	  for (int k = 0; k < sdim; k++)
	    ddshape(j,sdim*i+k) = dshape2(j,k);
      }  
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
  
  const IntegrationRule &
  FE_Tet0 :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  }

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



  const IntegrationRule &
  FE_Trig0 :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  }


  const IntegrationRule &
  FE_Trig1 :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  } 


  const IntegrationRule &
  FE_Trig2 :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 2);
  } 



  const IntegrationRule &
  FE_Quad0 :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  }

  const IntegrationRule &
  FE_Quad1:: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  }


  const IntegrationRule &
  FE_Quad2:: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  }


  const IntegrationRule &
  FE_Pyramid0 :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  }



  const IntegrationRule &
  FE_Pyramid1 :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  }



  const IntegrationRule &
  FE_Prism0 :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  }




  const IntegrationRule &
  FE_Prism1 :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), 1);
  }






  // Array<ScalarFiniteElement::IPData> FE_Prism2::ipdata;
  // ScalarFiniteElement<3>::IPDataArray FE_Prism2::ipdata;

  FE_Prism2 :: FE_Prism2()
    : ScalarFiniteElement<3> (ET_PRISM, 18, 3)
  {
    // if (!ipdata.Size()) CalcIPData(ET_PRISM, ipdata);
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
				FlatMatrixFixWidth<3> dshape) const
  {
    ScalarFiniteElement<3>::CalcDShape (ip, dshape);
  }




  // Array<ScalarFiniteElement::IPData> FE_Prism2aniso::ipdata;
  // ScalarFiniteElement<3>::IPDataArray FE_Prism2aniso::ipdata;


  FE_Prism2aniso :: FE_Prism2aniso()
    : ScalarFiniteElement<3> (ET_PRISM, 12, 3)
  {
    // if (!ipdata.Size()) CalcIPData(ET_PRISM, ipdata);
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
				     FlatMatrixFixWidth<3> dshape) const
  {
    ScalarFiniteElement<3>::CalcDShape (ip, dshape);
  }





  // Array<ScalarFiniteElement::IPData> FE_Prism2HBaniso::ipdata;
  // ScalarFiniteElement<3>::IPDataArray FE_Prism2HBaniso::ipdata;

  FE_Prism2HBaniso :: FE_Prism2HBaniso()
    : ScalarFiniteElement<3> (ET_PRISM, 12, 3)
  {
    // if (!ipdata.Size()) CalcIPData(ET_PRISM, ipdata);
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
  // ScalarFiniteElement<3>::IPDataArray FE_Prism3aniso::ipdata;

  FE_Prism3aniso :: FE_Prism3aniso() 
    : ScalarFiniteElement<3> (ET_PRISM)
  {
    // if (!ipdata.Size()) CalcIPData(ET_PRISM, ipdata);
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






  template class  FE_TSegmL2<0>;
  template class  FE_TSegmL2<1>;
  template class  FE_TSegmL2<2>;
  template class  FE_TSegmL2<3>;


  template class ScalarFiniteElement<1>;
  template class ScalarFiniteElement<2>;
  template class ScalarFiniteElement<3>;
}

