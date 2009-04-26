/*********************************************************************/
/* File:   myIntegrator.cpp                                          */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/

/*
  
My own simple integrators for the Poisson Equation

*/



#include <fem.hpp>   
#include "myIntegrator.hpp"

namespace ngfem
{
  
  /*
    Calculates the element matrix.

    Input is:
    the finite element: fel 
    the geometry of the element: eltrans
    
    Output is:
    the element matrix: elmat
    
    Efficient memorymanagement is provided my locheap
  */
  void MyLaplaceIntegrator ::
  AssembleElementMatrix (const FiniteElement & base_fel,
                         const ElementTransformation & eltrans, 
                         FlatMatrix<double> & elmat,
                         LocalHeap & lh) const
  {
    /*
      tell the compiler that we are expecting to get an scalar fe in 2D.
      if not, an exception will be raised
    */
    const ScalarFiniteElement<2> & fel =
      dynamic_cast<const ScalarFiniteElement<2> &> (base_fel);

    // number of element basis functions:
    int ndof = fel.GetNDof();

    elmat.AssignMemory (ndof, ndof, lh);
    elmat = 0;

    Matrix<> dshape_ref(ndof, 2); // gradient on reference element
    Matrix<> dshape(ndof, 2);     // gradient on mapped element

    /*
      get integration rule for element geometry, 
      integration order is 2 times element order
    */
    const IntegrationRule & ir = 
      SelectIntegrationRule (fel.ElementType(), 2*fel.Order());

    // loop over integration points
    for (int i = 0 ; i < ir.GetNIP(); i++)
      {
        
        // calculate Jacobi matrix in the integration point
        SpecificIntegrationPoint<2,2> sip(ir[i], eltrans, lh);

        // lambda(x)
        double lam = coef_lambda -> Evaluate (sip);

        /*
          gradient on reference element
          the i-th row of the matrix is the grad of the i-th shape function
        */
        fel.CalcDShape (ir[i], dshape_ref);
        
        // transform it for the mapped element
        dshape = dshape_ref * sip.GetJacobianInverse();
        
        // integration weight and Jacobi determinant
        double fac = sip.IP().Weight() * fabs (sip.GetJacobiDet());

        // elmat_{i,j} += (fac*lam) * InnerProduct (grad shape_i, grad shape_j)
        elmat += (fac*lam) * dshape * Trans(dshape);
      }     
  }









  /*
    Calculates the element vector.

    Input is:
    the finite element: fel 
    the geometry of the element: eltrans
    
    Output is:
    the element vector: elvec
    
    Efficient memorymanagement is provided my locheap
  */
  void MySourceIntegrator ::
  AssembleElementVector (const FiniteElement & base_fel,
                         const ElementTransformation & eltrans, 
                         FlatVector<double> & elvec,
                         LocalHeap & lh) const
  {
    const ScalarFiniteElement<2> & fel =
      dynamic_cast<const ScalarFiniteElement<2> &> (base_fel);

    int ndof = fel.GetNDof();

    elvec.AssignMemory (ndof, lh);
    elvec = 0;

    Vector<> shape(ndof); 

    const IntegrationRule & ir = 
      SelectIntegrationRule (fel.ElementType(), 2*fel.Order());

    for (int i = 0 ; i < ir.GetNIP(); i++)
      {
        SpecificIntegrationPoint<2,2> sip(ir[i], eltrans, lh);

        double f = coef_f -> Evaluate (sip);

        // calculate shape functions 
        fel.CalcShape (ir[i], shape);
        
        // integration weight and Jacobi determinant
        double fac = sip.IP().Weight() * fabs (sip.GetJacobiDet());

        // elvec_{i} += (fac*f) shape_i
        elvec += (fac*f) * shape;
      }     
  }








  
  namespace init_mylaplace
  {
    class Init
    { 
    public:  
      Init ();
    };        
    
    Init::Init()
    {
      // register the integrator mylaplace for 2D space, requiring one coefficient function
      GetIntegrators().AddBFIntegrator ("mylaplace", 2, 1,
                                        MyLaplaceIntegrator::Create);

      GetIntegrators().AddLFIntegrator ("mysource", 2, 1,
                                        MySourceIntegrator::Create);

    }

    Init init;
  }
}

