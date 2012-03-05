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
    HDivFiniteElement (ELEMENT_TYPE aeltype, int andof, int aorder)
      : FiniteElement (DIM, aeltype, andof, aorder) { ; } 

    HDivFiniteElement () 
    {
      ; // dimspace = D; 
    }

    ///
    virtual ~HDivFiniteElement () { ; }


    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip,
			    FlatMatrixFixWidth<DIM> shape) const = 0;

    /// compute div of shape
    virtual void CalcDivShape (const IntegrationPoint & ip,
			       FlatVector<> divshape) const;

    /// compute shape
    virtual void CalcMappedShape (const MappedIntegrationPoint<DIM,DIM> & sip,
				  FlatMatrixFixWidth<DIM> shape) const;

    /// compute div of shape
    virtual void CalcMappedDivShape (const MappedIntegrationPoint<DIM,DIM> & sip,
				     FlatVector<> divshape) const;



    const FlatMatrixFixWidth<DIM> GetShape (const IntegrationPoint & ip,
					    LocalHeap & lh) const
    {
      /*
      if (ip.IPNr() >= 0 && p_ipdata)
	{
	  return p_ipdata[ip.IPNr()].shape;
	}
      else
      */
	{
	  FlatMatrixFixWidth<DIM> shape(ndof, lh);
	  CalcShape (ip, shape);
	  return shape;
	}
    }

    const FlatVector<> GetDivShape (const IntegrationPoint & ip,
				    LocalHeap & lh) const
    {
      /*
      if (ip.IPNr() >= 0 && p_ipdata)
	{
	  return p_ipdata[ip.IPNr()].divshape;
	}
      else
      */
	{
	  FlatVector<> divshape(ndof, lh);
	  CalcDivShape (ip, divshape);
	  return divshape;
	}
    }

    
    virtual void Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, FlatMatrixFixWidth<D> vals) const
    {
      MatrixFixWidth<D> shape(ndof);
      for (int i = 0; i < ir.GetNIP(); i++)
	{
	  CalcShape (ir[i], shape);
	  vals.Row(i) = Trans(shape) * coefs;
	}
    }


  protected:
    ///
    // void CalcIPData (Array<IPData> & ipdata);

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




  /*
    HDivNormalFiniteElement
  */


  ///
  template <int D>
  class NGS_DLL_HEADER HDivNormalFiniteElement : public FiniteElement
  {
  public:
    enum { DIM = D };

  public:
    ///
    HDivNormalFiniteElement (ELEMENT_TYPE aeltype, int andof, int aorder)
      : FiniteElement (DIM, aeltype, andof, aorder){;}

    ///
    virtual ~HDivNormalFiniteElement () { ; }


    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip,
			    FlatVector<> shape) const = 0;

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



