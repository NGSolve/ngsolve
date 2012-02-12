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



namespace ngfem
{

  /**
     Transformation from reference element to physical element.
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

    /// is the element curved ?
    bool iscurved;

  public:
    /// polymorphism: if specific is set, use it.
    ElementTransformation * specific;
    ///
    ElementTransformation () { higher_integration_order = false; specific = NULL; }
    ///
    virtual ~ElementTransformation() { delete specific; }
    /// set data: is it a boundary, element number, and element index
    virtual void SetElement (bool aboundary, int aelnr, int aelindex)
    {
      // boundary = Boundary();  // aboundary;
      elnr = aelnr;
      elindex = aelindex;
      if (specific) specific -> SetElement (aboundary, aelnr, aelindex);
      // iscurved = false;
    }
    /// set geometric type of element
    void SetElementType (ELEMENT_TYPE aet) { eltype = aet; }
    /// return element number
    int GetElementNr () const { return elnr; }
    /// return element index
    int GetElementIndex () const { return elindex; }
    /// return element geometry type 
    ELEMENT_TYPE GetElementType () const { return eltype; }
    /// calculates the Jacobi matrix
    virtual void CalcJacobian (const IntegrationPoint & ip,
			       FlatMatrix<> dxdxi) const
    {
      specific -> CalcJacobian (ip, dxdxi);
    }

    /*
    /// Calculates the Jacobi matrix
    void CalcJacobian (const IntegrationPoint & ip,
		       FlatMatrix<> dxdxi,
		       LocalHeap & lh) const
    {
      CalcJacobian (ip, dxdxi);
    }
    */
    /// calculate the mapped point
    virtual void CalcPoint (const IntegrationPoint & ip,
			    FlatVector<> point) const
    {
      specific -> CalcPoint (ip, point);
    }

    /*
    /// calculate the mapped point
    void CalcPoint (const IntegrationPoint & ip,
		    FlatVector<> point,
		    LocalHeap & lh) const
    {
      CalcPoint (ip, point);
    }
    */

    /// calculate point and Jacobi matrix
    virtual void CalcPointJacobian (const IntegrationPoint & ip,
				    FlatVector<> point, FlatMatrix<> dxdxi) const
    {
      specific -> CalcPointJacobian (ip, point, dxdxi);
    }

    /*
    /// calculate point and Jacobi matrix
    void CalcPointJacobian (const IntegrationPoint & ip,
			    FlatVector<> point, FlatMatrix<> dxdxi,
			    LocalHeap & lh) const
    {
      CalcPointJacobian (ip, point, dxdxi);
    }
    */

    /// Calculate points and Jacobimatrices in all points of integrationrule
    virtual void CalcMultiPointJacobian (const IntegrationRule & ir,
					 BaseMappedIntegrationRule & mir) const
    {
      specific -> CalcMultiPointJacobian (ir, mir);
    }
    

    /// Calcs the normal vector in ip
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
  
  
    /// returns space dimension of physical elements
    virtual int SpaceDim () const
    {
      return specific->SpaceDim();
    }
    
    /// is it a mapping for boundary elements ?
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

    /// has the element non-constant Jacobian ?
    virtual bool IsCurvedElement() const
    {
      return specific -> IsCurvedElement();
    }

    virtual void GetSort (FlatArray<int> sort) const
    {
      if (specific) specific -> GetSort (sort);
    }

    /// return a mapped integration point on localheap
    virtual BaseMappedIntegrationPoint & operator() (const IntegrationPoint & ip, LocalHeap & lh) const
    {
      return (*specific) (ip, lh);
    }


    /// return a mapped integration rule on localheap
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
    FE_ElementTransformation ();

    ///
    ~FE_ElementTransformation ();

    ///
    virtual void SetElement (const FiniteElement * afel, int aelnr, int aelindex)
    {
      fel = static_cast<const ScalarFiniteElement<DIMS>*> (afel); 
      elnr = aelnr; 
      elindex = aelindex;
      SetElementType (fel->ElementType());
    }

    ///
    const FiniteElement & GetElement () const { return *fel; }

  
    ELEMENT_TYPE GetElementType () const { return fel->ElementType(); }
  
    ///
    virtual void CalcJacobian (const IntegrationPoint & ip,
			       FlatMatrix<> dxdxi) const;

    ///
    virtual void CalcPoint (const IntegrationPoint & ip, 
			    FlatVector<> point) const;

    ///
    virtual void CalcPointJacobian (const IntegrationPoint & ip,
				    FlatVector<> point, 
				    FlatMatrix<> dxdxi) const;




    virtual void CalcMultiPointJacobian (const IntegrationRule & ir,
					 BaseMappedIntegrationRule & bmir) const;

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
      return *new (lh) SpecificIntegrationPoint<DIMS,DIMR> (ip, *this);
    }

    virtual BaseMappedIntegrationRule & operator() (const IntegrationRule & ir, LocalHeap & lh) const
    {
      return *new (lh) MappedIntegrationRule<DIMS,DIMR> (ir, *this, lh);
    }

  };






}




#endif






