#ifndef FILE_ELEMENTTRANSFORMATION
#define FILE_ELEMENTTRANSFORMATION

/*********************************************************************/
/* File:   elementtransformation.hpp                                 */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/



/*
  Transformation from reference element to actual element
*/



/*
  extern "C"
  {
  DLL_HEADER void Ng_GetElementTransformation (int ei, const double * xi, 
  double * x, double * dxdxi);

  DLL_HEADER void Ng_GetSurfaceElementTransformation (int sei, const double * xi, 
  double * x, double * dxdxi);
  }

  namespace netgen
  {
  template <int DIM_EL, int DIM_SPACE> 
  DLL_HEADER void Ng_MultiElementTransformation (int elnr, int npts,
  const double * xi, size_t sxi,
  double * x, size_t sx,
  double * dxdxi, size_t sdxdxi);
  }
*/


namespace ngfem
{

  /**
     Transformation from reference element to physical element.
     Uses Netgen-meshaccess to perform transformation
  */
  class NGS_DLL_HEADER ElementTransformation
  {
  protected:
    /// element number
    int elnr;
    /// material property
    int elindex;
    /// element in R^dim
    // int dim;
    /// is it a boundary element ?
    // bool boundary;
    /// geometry of element
    ELEMENT_TYPE eltype;

    bool higher_integration_order;
    bool iscurved;

  public:
    ElementTransformation * specific;
    ///
    ElementTransformation () { higher_integration_order = false; specific = NULL; }
    ///
    virtual ~ElementTransformation() { delete specific; }
    ///
    virtual void SetElement (bool aboundary, int aelnr, int aelindex)
    {
      // boundary = Boundary();  // aboundary;
      elnr = aelnr;
      elindex = aelindex;
      if (specific) specific -> SetElement (aboundary, aelnr, aelindex);
      // iscurved = false;
    }
    ///
    void SetElementType (ELEMENT_TYPE aet) { eltype = aet; }
    ///
    int GetElementNr () const { return elnr; }
    ///
    int GetElementIndex () const { return elindex; }
    ///
    ELEMENT_TYPE GetElementType () const { return eltype; }
    ///
    virtual void CalcJacobian (const IntegrationPoint & ip,
			       FlatMatrix<> dxdxi) const
    {
      specific -> CalcJacobian (ip, dxdxi);
    }

    ///
    void CalcJacobian (const IntegrationPoint & ip,
		       FlatMatrix<> dxdxi,
		       LocalHeap & lh) const
    {
      CalcJacobian (ip, dxdxi);
    }

    ///
    virtual void CalcPoint (const IntegrationPoint & ip,
			    FlatVector<> point) const
    {
      specific -> CalcPoint (ip, point);
    }

    void CalcPoint (const IntegrationPoint & ip,
		    FlatVector<> point,
		    LocalHeap & lh) const
    {
      CalcPoint (ip, point);
    }

    ///
    virtual void CalcPointJacobian (const IntegrationPoint & ip,
				    FlatVector<> point, FlatMatrix<> dxdxi) const
    {
      specific -> CalcPointJacobian (ip, point, dxdxi);
    }

    void CalcPointJacobian (const IntegrationPoint & ip,
			    FlatVector<> point, FlatMatrix<> dxdxi,
			    LocalHeap & lh) const
    {
      CalcPointJacobian (ip, point, dxdxi);
    }


    /*
      template <int DIMS, int DIMR>
      void CalcMultiPointJacobian (const IntegrationRule & ipts,
      FlatArray<Vec<DIMR> > point, 
      FlatArray<Mat<DIMR,DIMS> > dxdxi,
      LocalHeap & lh) const
      {
      netgen::Ng_MultiElementTransformation <DIMS,DIMR> (elnr, ipts.Size(),
      &ipts[0](0), &ipts[1](0)-&ipts[0](0),
      &point[0](0), &point[1](0)-&point[0](0),
      &dxdxi[0](0,0), &dxdxi[1](0,0)-&dxdxi[0](0,0));
      }
    */



    virtual void CalcMultiPointJacobian (const IntegrationRule & ir,
					 BaseMappedIntegrationRule & mir) const
    {
      specific -> CalcMultiPointJacobian (ir, mir);
    }
    

    ///
    void CalcNormalVector (const IntegrationPoint & ip,
			   FlatVector<> nv,
			   LocalHeap & lh) const
    {
      if (Boundary())
	{
	  if (SpaceDim() == 2)
	    {
	      Mat<2,1> dxdxi;
	      CalcJacobian (ip, dxdxi);
	      // Ng_GetSurfaceElementTransformation (elnr+1, &ip(0), 0, &dxdxi(0));
	      double len = sqrt (sqr (dxdxi(0,0)) + sqr(dxdxi(1,0)));
	      nv(0) = -dxdxi(1,0) / len; //SZ 
	      nv(1) = dxdxi(0,0) / len;
	    }
	  else
	    {
	      Mat<3,2> dxdxi;
	      CalcJacobian (ip, dxdxi);
	      // Ng_GetSurfaceElementTransformation (elnr+1, &ip(0), 0, &dxdxi(0));
	      nv(0) = dxdxi(1,0) * dxdxi(2,1) - dxdxi(2,0) * dxdxi(1,1);
	      nv(1) = dxdxi(2,0) * dxdxi(0,1) - dxdxi(0,0) * dxdxi(2,1);
	      nv(2) = dxdxi(0,0) * dxdxi(1,1) - dxdxi(1,0) * dxdxi(0,1);
	      nv /= L2Norm (nv);
	    }
	}
    }
  
  
    ///
    virtual int SpaceDim () const
    {
      return specific->SpaceDim();
    }
    
    ///
    virtual bool Boundary() const
    {
      return specific->Boundary();
    }


    void SetHigherIntegrationOrder(void) {higher_integration_order = true;}
    void UnSetHigherIntegrationOrder(void) {higher_integration_order = false;}
    bool HigherIntegrationOrderSet(void) const 
    {
      return higher_integration_order;
    }
  
    virtual bool IsCurvedElement() const
    {
      return specific -> IsCurvedElement();
      // return iscurved;
    }

    virtual void GetSort (FlatArray<int> sort) const
    {
      if (specific) specific -> GetSort (sort);
    }

    virtual BaseSpecificIntegrationPoint & operator() (const IntegrationPoint & ip, LocalHeap & lh) const
    {
      return (*specific) (ip, lh);
    }


    virtual BaseMappedIntegrationRule & operator() (const IntegrationRule & ir, LocalHeap & lh) const
    {
      return (*specific) (ir, lh);
    }

  };











  /*
    Transformation from reference element to physical element.
    Uses finite element fel to describe mapping
  */
  template <int DIMS, int DIMR>
  class NGS_DLL_HEADER FE_ElementTransformation : public ElementTransformation
  {
    /// finite element defining transformation.
    const ScalarFiniteElement<DIMS> * fel;

    /// matrix with points, dim * np
    FlatMatrix<> pointmat;
    ///
    bool pointmat_ownmem;

    /// normal vectors (only surfelements)
    FlatMatrix<> nvmat;
  public:
    ///
    FE_ElementTransformation ()
      : pointmat(0,0,0), pointmat_ownmem(false), nvmat(0,0,0)
    { ; }

    ~FE_ElementTransformation ()
    {
      if (pointmat_ownmem) delete [] &pointmat(0,0); 
    }

    ///
    virtual void SetElement (const FiniteElement * afel, int aelnr, int aelindex)
    {
      fel = static_cast<const ScalarFiniteElement<DIMS>*> (afel); 
      elnr = aelnr; 
      elindex = aelindex;
    }

    ///
    const FiniteElement & GetElement () const { return *fel; }

  
    ELEMENT_TYPE GetElementType () const { return fel->ElementType(); }
  
    ///
    virtual void CalcJacobian (const IntegrationPoint & ip,
			       FlatMatrix<> dxdxi) const
    {
      for (int i = 0; i < DIMR; i++)
	dxdxi.Row(i) = fel->EvaluateGrad (ip, pointmat.Row(i));
    }


    ///
    virtual void CalcPoint (const IntegrationPoint & ip, 
			    FlatVector<> point) const
    {
      for (int i = 0; i < DIMR; i++)
	point(i) = fel->Evaluate (ip, pointmat.Row(i));
    }

    ///
    virtual void CalcPointJacobian (const IntegrationPoint & ip,
				    FlatVector<> point, 
				    FlatMatrix<> dxdxi) const
    {
      CalcPoint (ip, point);
      CalcJacobian (ip, dxdxi);
    }


    template <typename T0, typename T1, typename T2>
    void CalcMultiPointJacobian (FlatArray<T0> ipts,
				 FlatArray<T1> point, 
				 FlatArray<T2> dxdxi,
				 LocalHeap & lh) const
    {
      for (int i = 0; i < ipts.Size(); i++)
	CalcPointJacobian (IntegrationPoint (ipts[i]), point[i], dxdxi[i], lh);
    }


    virtual void CalcMultiPointJacobian (const IntegrationRule & ir,
					 BaseMappedIntegrationRule & bmir) const
    {
      MappedIntegrationRule<DIMS,DIMR> & mir = static_cast<MappedIntegrationRule<DIMS,DIMR> &> (bmir);
      for (int i = 0; i < ir.Size(); i++)
	{
	  CalcPointJacobian (ir[i], mir[i].Point(), mir[i].Jacobian());
	  mir[i].Compute();
	}
    }

    ///
    const FlatMatrix<> & PointMatrix () const { return pointmat; }
    ///
    FlatMatrix<> & PointMatrix () { return pointmat; }
    ///
    void AllocPointMatrix (int spacedim, int vertices)
    {
      if (pointmat_ownmem) delete [] &pointmat(0,0);
      pointmat.AssignMemory (spacedim, vertices, new double[spacedim*vertices]);
      pointmat_ownmem = 1;
    }
    ///
    const FlatMatrix<> & NVMatrix () const { return nvmat; }

    template <typename T>
    void CalcNormalVector (const IntegrationPoint & ip,
			   T & nv,
			   LocalHeap & lh) const
    {
      for (int i = 0; i < nvmat.Height(); i++)
	nv(i) = nvmat(i,0) ;
    }

    ///
    int SpaceDim () const
    {
      return pointmat.Height(); 
    }

    bool Boundary(void) const
    {
      return pointmat.Height() != ElementTopology::GetSpaceDim (fel->ElementType());
    }

    void GetSort (FlatArray<int> sort) const
    { ; }


    virtual BaseSpecificIntegrationPoint & operator() (const IntegrationPoint & ip, LocalHeap & lh) const
    {
      return *new (lh) SpecificIntegrationPoint<DIMS,DIMR> (ip, *this, lh);
    }

    virtual BaseMappedIntegrationRule & operator() (const IntegrationRule & ir, LocalHeap & lh) const
    {
      return *new (lh) MappedIntegrationRule<DIMS,DIMR> (ir, *this, lh);
    }

  };






}




#endif






