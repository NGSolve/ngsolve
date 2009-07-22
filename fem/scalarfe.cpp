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
  double ScalarFiniteElement<D> :: 
  Evaluate (const IntegrationPoint & ip, FlatVector<double> x) const
  {
    VectorMem<20, double> shape(ndof);
    CalcShape (ip, shape);
    return InnerProduct (shape, x);
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



  /*
  template<int D>
  const IntegrationRule &
  ScalarFiniteElement<D> :: NodalIntegrationRule() const
  {
    return GetIntegrationRules().SelectNodalIntegrationRule (ElementType(), Order());
  }
  */

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
				  



  void FE_Tet1 :: GetDofs (Array<Dof> & dofs) const
  {
    /*
    Dof da[] = { Dof (Node (NT_VERTEX, 0), 0),
		 Dof (Node (NT_VERTEX, 1), 0),
		 Dof (Node (NT_VERTEX, 2), 0),
		 Dof (Node (NT_VERTEX, 3), 0) };
    */	 

    dofs.SetSize (0);
    for (int i = 0; i < 4; i++)
      dofs.Append (Dof (Node (NT_VERTEX, i), 0));
  }







  template class  FE_TSegmL2<0>;
  template class  FE_TSegmL2<1>;
  template class  FE_TSegmL2<2>;
  template class  FE_TSegmL2<3>;


  template class ScalarFiniteElement<1>;
  template class ScalarFiniteElement<2>;
  template class ScalarFiniteElement<3>;
}

