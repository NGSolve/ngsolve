#ifndef FILE_ELEMENTTRANSFORMATION
#define FILE_ELEMENTTRANSFORMATION

/*********************************************************************/
/* File:   elementtransformation.hpp                                 */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

namespace ngfem
{

  /*
    Transformation from reference element to actual element
  */


#ifndef NETGEN_ELTRANS

  /*
    Transformation from reference element to physical element.
    Uses finite element fel to describe mapping
  */
  class NGS_DLL_HEADER ElementTransformation
  {
    /// finite element defining transformation.
    const FiniteElement * fel;

    /// matrix with points, dim * np
    FlatMatrix<> pointmat;
    ///
    bool pointmat_ownmem;

    /// normal vectors (only surfelements)
    FlatMatrix<> nvmat;

    /// element number
    int elnr;

    /// material property
    int elindex;

    bool higher_integration_order;
  public:
    ///
    ElementTransformation ()
      : pointmat(0,0,0), pointmat_ownmem(false), nvmat(0,0,0), higher_integration_order(false)
    { ; }

    ~ElementTransformation ()
    {
      if (pointmat_ownmem) delete [] &pointmat(0,0); 
    }

    ///
    void SetElement (const FiniteElement * afel, int aelnr, int aelindex)
    {
      fel = afel; 
      elnr = aelnr; 
      elindex = aelindex;
    }


    ///
    const FiniteElement & GetElement () const { return *fel; }
    ///
    int GetElementNr () const { return elnr; }
    ///
    int GetElementIndex () const { return elindex; }

    void SetHigherIntegrationOrder(void) {higher_integration_order = true;}
    void UnSetHigherIntegrationOrder(void) {higher_integration_order = false;}

    bool HigherIntegrationOrderSet(void) const
    {
      // (*testout) << "eltrans hios " << higher_integration_order << endl;
      return higher_integration_order;
    }
  
  


    ELEMENT_TYPE GetElementType () const
    { return fel->ElementType(); }
  
    ///
    template <typename T>
    void CalcJacobian (const IntegrationPoint & ip,
		       MatExpr<T> & dxdxi,
		       LocalHeap & lh) const
    {
      switch (ElementTopology::GetSpaceDim (fel->ElementType()))
	{
	case 1:
	  dxdxi = pointmat * static_cast<const ScalarFiniteElement<1>*> (fel)->GetDShape(ip, lh);
	  break;
	case 2:
	  dxdxi = pointmat * static_cast<const ScalarFiniteElement<2>*> (fel)->GetDShape(ip, lh);
	  break;
	case 3:
	  dxdxi = pointmat * static_cast<const ScalarFiniteElement<3>*> (fel)->GetDShape(ip, lh);
	  break;
	}

    }

    ///
    template <typename T>
    void CalcJacobian (const IntegrationPoint & ip,
		       MatExpr<T> & dxdxi) const
    {
      LocalHeap lh(1000, "calc jac");
      switch (ElementTopology::GetSpaceDim (fel->ElementType())) 
	{
	case 1:
	  dxdxi = pointmat * static_cast<const ScalarFiniteElement<1>*> (fel)->GetDShape(ip, lh);
	  break;
	case 2:
	  dxdxi = pointmat * static_cast<const ScalarFiniteElement<2>*> (fel)->GetDShape(ip, lh);
	  break;
	case 3:
	  dxdxi = pointmat * static_cast<const ScalarFiniteElement<3>*> (fel)->GetDShape(ip, lh);
	  break;
	}
    }


    ///
    template <typename T>
    void CalcPoint (const IntegrationPoint & ip, 
		    // MatExpr<T> & point,
		    T & point,
		    LocalHeap & lh) const
    {
      switch (ElementTopology::GetSpaceDim (fel->ElementType())) 
	{
	case 1:
	  point = pointmat * static_cast<const ScalarFiniteElement<1>*> (fel)->GetShape(ip, lh);
	  break;
	case 2:
	  point = pointmat * static_cast<const ScalarFiniteElement<2>*> (fel)->GetShape(ip, lh);
	  break;
	case 3:
	  point = pointmat * static_cast<const ScalarFiniteElement<3>*> (fel)->GetShape(ip, lh);
	  break;
	}

      // point = pointmat * fel->GetShape(ip, lh);
    }

    /*
      template <int D, typename T1, typename T2>
      void CalcPointJacobian (const IntegrationPoint & ip,
      T1 & point, T2 & dxdxi,
      LocalHeap & lh) const
    */
    template <int S, int R>
    void CalcPointJacobian (const IntegrationPoint & ip,
			    // T1 & point, T2 & dxdxi,
			    Vec<R> & point, Mat<R,S> & dxdxi,
			    LocalHeap & lh) const
    {
      dxdxi = pointmat * static_cast<const ScalarFiniteElement<S>*> (fel)->GetDShape(ip, lh);
      point = pointmat * static_cast<const ScalarFiniteElement<S>*> (fel)->GetShape(ip, lh);
    
      /*
      // switch (fel->SpatialDim())
      {
      case 1:
      dxdxi = pointmat * static_cast<const ScalarFiniteElement<1>*> (fel)->GetDShape(ip, lh);
      point = pointmat * static_cast<const ScalarFiniteElement<1>*> (fel)->GetShape(ip, lh);
      break;
      case 2:
      dxdxi = pointmat * static_cast<const ScalarFiniteElement<2>*> (fel)->GetDShape(ip, lh);
      point = pointmat * static_cast<const ScalarFiniteElement<2>*> (fel)->GetShape(ip, lh);
      break;
      case 3:
      dxdxi = pointmat * static_cast<const ScalarFiniteElement<3>*> (fel)->GetDShape(ip, lh);
      point = pointmat * static_cast<const ScalarFiniteElement<3>*> (fel)->GetShape(ip, lh);
      break;
      }
      */
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



    BaseSpecificIntegrationPoint & operator() (const IntegrationPoint & ip, LocalHeap & lh) const
    {
      if (Boundary())
	{
	  if (SpaceDim() == 2)
	    return *new SpecificIntegrationPoint<1,2> (ip, *this, lh);
	  else
	    return *new SpecificIntegrationPoint<2,3> (ip, *this, lh);
	}
      else
	{
	  if (SpaceDim() == 2)
	    return *new SpecificIntegrationPoint<2,2> (ip, *this, lh);
	  else
	    return *new SpecificIntegrationPoint<3,3> (ip, *this, lh);
	}
    }

    BaseMappedIntegrationRule & operator() (const IntegrationRule & ir, LocalHeap & lh) const
    {
      if (Boundary())
	{
	  if (SpaceDim() == 2)
	    return *new MappedIntegrationRule<1,2> (ir, *this, lh);
	  else
	    return *new MappedIntegrationRule<2,3> (ir, *this, lh);
	}
      else
	{
	  if (SpaceDim() == 2)
	    return *new MappedIntegrationRule<2,2> (ir, *this, lh);
	  else
	    return *new MappedIntegrationRule<3,3> (ir, *this, lh);
	}
    }


  };




#else




  /**
     Transformation from reference element to physical element.
     Uses Netgen-meshaccess to perform transformation
  */
  class NGS_DLL_HEADER ElementTransformation
  {
    /// element number
    int elnr;

    /// material property
    int elindex;
  
    /// element in R^dim
    int dim;
    /// is it a boundary element ?
    bool boundary;
    /// geometry of element
    ELEMENT_TYPE eltype;

    /*
    const Array<Vec<3> > * pts;
    const Array<Mat<3,3> > * dxdxis;
    const Array<int> * first_of_element;
    */

    bool higher_integration_order;

    bool iscurved;

    /*
    double buffer[100];
    bool buffervalid;
    */

  public:
    ///
    ElementTransformation () { /* pts = 0; */ higher_integration_order = false; /* buffervalid = false; */ }
    /*
      ElementTransformation (Array<Vec<3> > * apts,
      Array<Mat<3,3> > * adxdxis,
      Array<int> * afirst_of_element)
      : pts(apts), dxdxis(adxdxis), first_of_element(afirst_of_element) { ; }
    */

    /*
    void SetGeometryData (const Array<Vec<3> > * apts,
			  const Array<Mat<3,3> > * adxdxis,
			  const Array<int> * afirst_of_element)
    {
      pts = apts;
      dxdxis = adxdxis;
      first_of_element = afirst_of_element;
    }
    */

    ///
    void SetElement (bool aboundary, int aelnr, int aelindex)
    {
      boundary = aboundary;
      elnr = aelnr;
      elindex = aelindex;
      dim = Ng_GetDimension();
      if (boundary)
	{
	  iscurved = Ng_IsSurfaceElementCurved (elnr+1);
	}
      else
	{
	  iscurved = Ng_IsElementCurved (elnr+1);
	}
      // buffervalid = 0;
    }

    void SetElementType (ELEMENT_TYPE aet) { eltype = aet; }
    ///
    int GetElementNr () const { return elnr; }
    ///
    int GetElementIndex () const { return elindex; }
    ///
    ELEMENT_TYPE GetElementType () const { return eltype; }
  
    ///
    template <typename T>
    void CalcJacobian (const IntegrationPoint & ip,
		       T & dxdxi) const
    {
      if (boundary)
	Ng_GetSurfaceElementTransformation (elnr+1, &ip(0), 0, &dxdxi(0,0));
      else
	{
	  Ng_GetElementTransformation (elnr+1, &ip(0), 0, &dxdxi(0,0));
	  /*
	  Ng_GetBufferedElementTransformation (elnr+1, &ip(0),
					       0, &dxdxi(0,0), 
					       const_cast<double*> (&buffer[0]), buffervalid);
	  const_cast<bool&> (buffervalid) = true;
	  */
	}
    }

    template <typename T>
    void CalcJacobian (const IntegrationPoint & ip,
		       T & dxdxi,
		       LocalHeap & lh) const
    {
      if (boundary)
	Ng_GetSurfaceElementTransformation (elnr+1, &ip(0),
					    0, &dxdxi(0,0));
					  
      else
	{
	  Ng_GetElementTransformation (elnr+1, &ip(0), 0, &dxdxi(0,0));
	  /*
	  Ng_GetBufferedElementTransformation (elnr+1, &ip(0),
					       0, &dxdxi(0,0), 
					       const_cast<double*> (&buffer[0]), buffervalid);
	  const_cast<bool&> (buffervalid) = true;
	  */
	}
    }




    ///
    template <typename T>
    void CalcPoint (const IntegrationPoint & ip,
		    T & point,
		    LocalHeap & lh) const
    {
      if (boundary)
	Ng_GetSurfaceElementTransformation (elnr+1, &ip(0), &point(0), 0);
      else
	{
	  Ng_GetElementTransformation (elnr+1, &ip(0), &point(0), 0);

	  /*
	  Ng_GetBufferedElementTransformation (elnr+1, &ip(0), &point(0), 0, 
					       const_cast<double*> (&buffer[0]), buffervalid);
	  const_cast<bool&> (buffervalid) = true;
	  */
	}
    }


    // template <int D>
    // template <typename T1, typename T2>
    template <int S, int R>
    void CalcPointJacobian (const IntegrationPoint & ip,
			    // T1 & point, T2 & dxdxi,
			    Vec<R> & point, Mat<R,S> & dxdxi,
			    LocalHeap & lh) const
    {
      if (boundary)
	Ng_GetSurfaceElementTransformation (elnr+1, &ip(0), &point(0), &dxdxi(0,0));
      else
	{
	  /*
	  if (ip.precomputed_geometry && pts)
	    {
	      point = (*pts)[ (*first_of_element)[elnr] + ip.Nr()];
	      dxdxi = (*dxdxis)[ (*first_of_element)[elnr] + ip.Nr()];
	    }
	  else
	  */
	    {
	      Ng_GetElementTransformation (elnr+1, &ip(0), &point(0), &dxdxi(0,0));
	      /*
	      Ng_GetBufferedElementTransformation (elnr+1, &ip(0), &point(0), &dxdxi(0,0),
						   const_cast<double*> (&buffer[0]), buffervalid);
	      const_cast<bool&> (buffervalid) = true;
	      */
	    }
	}
    }


    template <typename T0, int DIMS, int DIMR>
    void CalcMultiPointJacobian (FlatArray<T0> ipts,
				 FlatArray<Vec<DIMR> > point, 
				 FlatArray<Mat<DIMR,DIMS> > dxdxi,
				 LocalHeap & lh) const
    {
      netgen::Ng_MultiElementTransformation <DIMS,DIMR> (elnr, ipts.Size(),
							 &ipts[0](0), &ipts[1](0)-&ipts[0](0),
							 &point[0](0), &point[1](0)-&point[0](0),
							 &dxdxi[0](0,0), &dxdxi[1](0,0)-&dxdxi[0](0,0));
    }










    ///
    template <typename T>
    void CalcNormalVector (const IntegrationPoint & ip,
			   T & nv,
			   LocalHeap & lh) const
    {
      //cout << "calc normal vec" << endl;
      //cout << "bound = " << boundary << ", dim = " << dim << endl;
      if (boundary)
	{
	  if (dim == 2)
	    {
	      Mat<2,1> dxdxi;
	      Ng_GetSurfaceElementTransformation (elnr+1, &ip(0), 0, &dxdxi(0));
	      //  cout << "dxdxi = " << dxdxi << endl;
	      double len = sqrt (sqr (dxdxi(0,0)) + sqr(dxdxi(1,0)));
	      //nv(0) = dxdxi(1,0) / len;  // original 
	      //nv(1) = -dxdxi(0,0) / len;
	      nv(0) = -dxdxi(1,0) / len; //SZ 
	      nv(1) = dxdxi(0,0) / len;
	      //cout << "nv = " << nv << endl;
	      // (*testout)<<"NormalVector="<<nv<<endl;
	    }
	  else
	    {
	      // berechne aus dxdxi

	      Mat<3,2> dxdxi;
	      Ng_GetSurfaceElementTransformation (elnr+1, &ip(0), 0, &dxdxi(0));

	      nv(0) = dxdxi(1,0) * dxdxi(2,1) - dxdxi(2,0) * dxdxi(1,1);
	      nv(1) = dxdxi(2,0) * dxdxi(0,1) - dxdxi(0,0) * dxdxi(2,1);
	      nv(2) = dxdxi(0,0) * dxdxi(1,1) - dxdxi(1,0) * dxdxi(0,1);
	      nv /= L2Norm (nv);
	    }

	}
    }
  
  
    ///
    int SpaceDim () const
    {
      return dim;
    }

    bool Boundary(void) const
    {
      return boundary;
    }


    void SetHigherIntegrationOrder(void) {higher_integration_order = true;}
    void UnSetHigherIntegrationOrder(void) {higher_integration_order = false;}

    bool HigherIntegrationOrderSet(void) const 
    {
      return higher_integration_order;
    }
  
    bool IsCurvedElement() const
    {
      return iscurved;
    }


    void GetSort (FlatArray<int> sort) const
    {
      int vnums[12];
      if (boundary)
	Ng_GetSurfaceElement (elnr+1, vnums);
      else
	Ng_GetElement (elnr+1, vnums);

      switch (eltype)
	{
	case ET_TRIG:
	  for (int i = 0; i < 3; i++) sort[i] = i;
	  if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
	  if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);
	  if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
	  // vnums[sort[0]] < vnums[sort[1]] < vnums[sort[2]]
	  break; 

	case ET_TET:
	  for (int i = 0; i < 4; i++) sort[i] = i;
	  if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
	  if (vnums[sort[2]] > vnums[sort[3]]) Swap (sort[2], sort[3]);
	  if (vnums[sort[0]] > vnums[sort[2]]) Swap (sort[0], sort[2]);
	  if (vnums[sort[1]] > vnums[sort[3]]) Swap (sort[1], sort[3]);
	  if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);

	  // vnums[sort[0]] < vnums[sort[1]] < vnums[sort[2]] < vnums[sort[3]]
	  break; 

	case ET_PRISM:
	  for (int i = 0; i < 6; i++) sort[i] = i;

	  if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
	  if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);
	  if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
        
	  if (vnums[sort[3]] > vnums[sort[4]]) Swap (sort[3], sort[4]);
	  if (vnums[sort[4]] > vnums[sort[5]]) Swap (sort[4], sort[5]);
	  if (vnums[sort[3]] > vnums[sort[4]]) Swap (sort[3], sort[4]);
	  break;

	default:
	  throw Exception ("undefined eltype in ElementTransformation::GetSort()\n");
	}
    }



    BaseSpecificIntegrationPoint & operator() (const IntegrationPoint & ip, LocalHeap & lh) const
    {
      if (boundary)
	{
	  if (dim == 2)
	    return *new SpecificIntegrationPoint<1,2> (ip, *this, lh);
	  else
	    return *new SpecificIntegrationPoint<2,3> (ip, *this, lh);
	}
      else
	{
	  if (dim == 2)
	    return *new SpecificIntegrationPoint<2,2> (ip, *this, lh);
	  else
	    return *new SpecificIntegrationPoint<3,3> (ip, *this, lh);
	}
    }


    BaseMappedIntegrationRule & operator() (const IntegrationRule & ir, LocalHeap & lh) const
    {
      if (Boundary())
	{
	  if (SpaceDim() == 2)
	    return *new MappedIntegrationRule<1,2> (ir, *this, lh);
	  else
	    return *new MappedIntegrationRule<2,3> (ir, *this, lh);
	}
      else
	{
	  if (SpaceDim() == 2)
	    return *new MappedIntegrationRule<2,2> (ir, *this, lh);
	  else
	    return *new MappedIntegrationRule<3,3> (ir, *this, lh);
	}
    }

  };



#endif



}


#endif






