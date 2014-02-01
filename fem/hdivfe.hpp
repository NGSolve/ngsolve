#ifndef FILE_HDIVFE
#define FILE_HDIVFE

/*********************************************************************/
/* File:   hdivfe.hpp                                                */
/* Author: Joachim Schoeberl                                         */
/* Date:   5. Jul. 2001                                              */
/*********************************************************************/

namespace ngfem
{

  /**
     Finite Elements for H(div)
     Raviart-Thomas, BDM, BDFM
  */
  template <int D>
  class NGS_DLL_HEADER HDivFiniteElement : public FiniteElement
  {
  public:
    enum { DIM = D };

  public:
    ///
    INLINE HDivFiniteElement (int andof, int aorder)
      : FiniteElement (andof, aorder) { ; } 

    ///
    INLINE HDivFiniteElement () { ; }

    ///
    virtual ~HDivFiniteElement () { ; }

    /// 
    virtual string ClassName() const;

    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip,
			    SliceMatrix<> shape) const = 0;

    /// compute div of shape
    virtual void CalcDivShape (const IntegrationPoint & ip,
			       SliceVector<> divshape) const;

    /// compute shape
    virtual void CalcMappedShape (const MappedIntegrationPoint<DIM,DIM> & sip,
				  SliceMatrix<> shape) const;

    /// compute div of shape
    virtual void CalcMappedDivShape (const MappedIntegrationPoint<DIM,DIM> & sip,
				     SliceVector<> divshape) const;



    INLINE const FlatMatrixFixWidth<DIM> GetShape (const IntegrationPoint & ip,
                                                   LocalHeap & lh) const
    {
      FlatMatrixFixWidth<DIM> shape(ndof, lh);
      CalcShape (ip, shape);
      return shape;
    }

    INLINE const FlatVector<> GetDivShape (const IntegrationPoint & ip,
                                           LocalHeap & lh) const
    {
      FlatVector<> divshape(ndof, lh);
      CalcDivShape (ip, divshape);
      return divshape;
    }

    
    virtual void Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, 
			   FlatMatrixFixWidth<D> vals) const
    {
      MatrixFixWidth<D> shape(ndof);
      for (int i = 0; i < ir.GetNIP(); i++)
	{
	  CalcShape (ir[i], shape);
	  vals.Row(i) = Trans(shape) * coefs;
	}
    }

    virtual void GetFacetDofs(int i, Array<int> & dnums) const
    { cout  << " GetFacetDofs for nothing " << endl; dnums.SetSize(0);}; 


  protected:

    /// compute basis, will be orthogonalized
    virtual void CalcShape1 (const IntegrationPoint & ip,
			     FlatMatrixFixWidth<DIM> shape) const { ; }

    ///
    virtual void CalcShape2 (const IntegrationPoint & ip,
			     FlatMatrixFixWidth<DIM> shape) const { ; }

    ///
    void ComputeFaceMoments (int fnr, ScalarFiniteElement<DIM-1> & testfe,
			     FlatMatrix<> & moments,
			     int order, int shape = 1) const;
  };




  /**
    HDivNormalFiniteElement
  */

  template <int D>
  class NGS_DLL_HEADER HDivNormalFiniteElement : public FiniteElement
  {
  public:
    enum { DIM = D };

  public:
    ///
    HDivNormalFiniteElement (int andof, int aorder)
      : FiniteElement (andof, aorder){;}

    ///
    virtual ~HDivNormalFiniteElement () { ; }

    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip,
			    FlatVector<> shape) const = 0;

    ///
    const FlatVector<> GetShape (const IntegrationPoint & ip,
				 LocalHeap & lh) const
    {
      FlatVector<> shape(ndof, lh);
      CalcShape (ip, shape);
      return shape;
    }

  };



}

#endif



