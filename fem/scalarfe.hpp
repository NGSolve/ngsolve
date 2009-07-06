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
  class ScalarFiniteElement : public FiniteElement
  {
  public:
    virtual string ClassName(void) const {return "ScalarFiniteElement";}
    ///
    ScalarFiniteElement () { dimspace = D; }
    ///
    ScalarFiniteElement (ELEMENT_TYPE aeltype, 
			 int andof = 0, int aorder = 0)
      : FiniteElement (D, aeltype, andof, aorder) { ; }
    ///
    virtual ~ScalarFiniteElement () { ; }
    ///
    // virtual const IntegrationRule & NodalIntegrationRule() const;

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

    virtual double
    Evaluate (const IntegrationPoint & ip, 
	      FlatVector<double> x, LocalHeap & lh) const
    {
      return InnerProduct (GetShape(ip, lh), x);
    }  


    /**
       returns second derivatives in point ip.
       returns stored values for valid ip.IPNr(), else computes values
    */
    const FlatMatrix<> GetDDShape (const IntegrationPoint & ip, LocalHeap & lh) const
    {
      FlatMatrix<> ddshape(ndof, dimspace*dimspace, lh);
      CalcDDShape (ip, ddshape);
      return ddshape;
    }

    /// compute dshape, matrix: ndof x (spacedim spacedim)
    virtual void CalcDDShape (const IntegrationPoint & ip, 
			      FlatMatrix<> ddshape) const;



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







  


  
  
  template <int DIM>
  class DShapeElement
  {
    double * data;
  public:
    DShapeElement (double * adata) : data(adata) { ; }
    void operator= (AutoDiff<DIM> ad) 
    { 
      for (int i = 0; i < DIM; i++) 
        data[i] = ad.DValue(i); 
    }
  };



  template <int DIM>
  class DShapeAssign
  {
    double * dshape;
  public:
    DShapeAssign (FlatMatrixFixWidth<DIM> mat)
    { dshape = &mat(0,0); }

    DShapeAssign (double * adshape)
    { dshape = adshape; }

    DShapeElement<DIM> operator[] (int i) const
    { return DShapeElement<DIM> (dshape + i*DIM); }

    const DShapeAssign Addr (int i) const
    { return DShapeAssign (dshape+i*DIM); } 
  };

  








  /**
     Base-element for template polymorphism.
     Barton and Nackman Trick
  */

  template <class FEL, ELEMENT_TYPE ET, int NDOF, int ORDER>
  class T_ScalarFiniteElement : public ScalarFiniteElement<ET_trait<ET>::DIM>
  {

  public:
    
  protected:
    enum { DIM = ET_trait<ET>::DIM };

    T_ScalarFiniteElement ()
      : ScalarFiniteElement<DIM> (ET, NDOF, ORDER) { ; }

    virtual ~T_ScalarFiniteElement() { ; }

  public:

    /*
  const FlatVec<NDOF> & GetShape (const IntegrationPoint & ip,
                                  LocalHeap & lh) const
    {
      ;
    }

    const Mat<NDOF,DIM> & GetDShape (const IntegrationPoint & ip,
				     LocalHeap & lh) const
    {
      ;
    }
    */


    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const
    {
      double pt[DIM];
      for (int i = 0; i < DIM; i++) pt[i] = ip(i);
      FEL::T_CalcShape (pt, shape); 
    }

    static void CalcShapeStat (const IntegrationPoint & ip, 
                               FlatVector<> shape)
    {
      double pt[DIM];
      for (int i = 0; i < DIM; i++) pt[i] = ip(i);
      FEL::T_CalcShape (pt, shape); 
    }
    
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<DIM> dshape) const
    {
      AutoDiff<DIM> adp[DIM];
      for (int i = 0; i < DIM; i++)
        adp[i] = AutoDiff<DIM> (ip(i), i);
      
      DShapeAssign<DIM> ds(dshape); 
      FEL::T_CalcShape (adp, ds);
    }

    static void CalcDShapeStat (const IntegrationPoint & ip, 
				FlatMatrixFixWidth<DIM> dshape)
    {
      AutoDiff<DIM> adp[DIM];
      for (int i = 0; i < DIM; i++)
        adp[i] = AutoDiff<DIM> (ip(i), i);
      
      DShapeAssign<DIM> ds(dshape); 
      FEL::T_CalcShape (adp, ds);
    }

    virtual void 
    CalcMappedDShape (const SpecificIntegrationPoint<DIM,DIM> & sip, 
                      FlatMatrixFixWidth<DIM> dshape) const
    {
      AutoDiff<DIM> adp[DIM];
      
      for (int i = 0; i < DIM; i++)
        adp[i].Value() = sip.IP()(i);
      
      for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
          adp[i].DValue(j) = sip.GetJacobianInverse()(i,j);
      
      DShapeAssign<DIM> ds(dshape); 
      FEL::T_CalcShape (adp, ds);
    }
  };











}

#endif
