#ifndef FILE_SCALARFE
#define FILE_SCALARFE

/*********************************************************************/
/* File:   scalarfe.hpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

namespace ngfem
{


  /**
     Scalar finite element.
     Provides shape functions and derivaties.
  */
  template <int D>
  class NGS_DLL_HEADER ScalarFiniteElement : public FiniteElement
  {
  public:
    virtual string ClassName(void) const {return "ScalarFiniteElement";}
    ///
    ScalarFiniteElement () { ; } // dimspace = D; }
    ///
    ScalarFiniteElement (ELEMENT_TYPE aeltype, 
			 int andof = 0, int aorder = 0)
      : FiniteElement (D, aeltype, andof, aorder) 
    { ; }

    ///
    virtual ~ScalarFiniteElement () { ; }

    /**
       returns shape functions in point ip.
       returns stored values for valid ip.IPNr(), else computes values
    */
    const FlatVector<> GetShape (const IntegrationPoint & ip, 
				 LocalHeap & lh) const
    {
      FlatVector<> shape(ndof, lh);
      CalcShape (ip, shape);
      return shape;
    }

    /**
       returns derivatives in point ip.
       returns stored values for valid ip.IPNr(), else computes values
    */
    const FlatMatrixFixWidth<D> 
    GetDShape (const IntegrationPoint & ip, LocalHeap & lh) const
    {
      FlatMatrixFixWidth<D> dshape(ndof, lh);
      CalcDShape (ip, dshape);
      return dshape;
    }


    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const = 0;
  
    /// compute dshape, matrix: ndof x spacedim
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<D> dshape) const;

    /// compute dshape, matrix: ndof x spacedim
    virtual void CalcMappedDShape (const SpecificIntegrationPoint<D,D> & sip, 
				   FlatMatrixFixWidth<D> dshape) const;




    /**
       returns second derivatives in point ip.
       returns stored values for valid ip.IPNr(), else computes values
    */
    const FlatMatrix<> GetDDShape (const IntegrationPoint & ip, LocalHeap & lh) const
    {
      FlatMatrix<> ddshape(ndof, D*D, lh);
      CalcDDShape (ip, ddshape);
      return ddshape;
    }

    /// compute dshape, matrix: ndof x (spacedim spacedim)
    virtual void CalcDDShape (const IntegrationPoint & ip, 
			      FlatMatrix<> ddshape) const;



    /// evaluate \sum x_i shape_i
    virtual double Evaluate (const IntegrationPoint & ip, FlatVector<> x) const;
    virtual Vec<D> EvaluateGrad (const IntegrationPoint & ip, FlatVector<> x) const;

    
    virtual void Evaluate (const IntegrationRule & ir, FlatVector<> coefs, FlatVector<> values) const;
    virtual void EvaluateGrad (const IntegrationRule & ir, FlatVector<> coefs, FlatMatrixFixWidth<D> values) const;

    virtual void EvaluateTrans (const IntegrationRule & ir, FlatVector<> values, FlatVector<> coefs) const;
    virtual void EvaluateGradTrans (const IntegrationRule & ir, FlatMatrixFixWidth<D> values, FlatVector<> coefs) const;



    virtual void EvaluateShapeGrid (const IntegrationRuleTP<D> & ir,
				    const FlatVector<double> coefs,
				    FlatVector<double> gridvalues,
				    LocalHeap & lh) const;
				  
    virtual void EvaluateShapeGridTrans (const IntegrationRuleTP<D> & ir,
					 const FlatVector<double> gridvalues,
					 FlatVector<double> coefs,
					 LocalHeap & lh) const;
				  
    virtual void EvaluateDShapeGrid (const IntegrationRuleTP<D> & ir,
				     const FlatVector<double> coefs,
				     FlatMatrixFixWidth<D> gridvalues,
				     LocalHeap & lh) const;
				  
    virtual void EvaluateDShapeGridTrans (const IntegrationRuleTP<D> & ir,
					  const FlatMatrixFixWidth<D> gridvalues,
					  FlatVector<double> coefs,
					  LocalHeap & lh) const;
  };







  








}

#endif
