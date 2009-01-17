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


#ifndef NETGEN_ELTRANS

/*
    Transformation from reference element to physical element.
    Uses finite element fel to describe mapping
*/
class ElementTransformation
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
    : pointmat(0,0,0), nvmat(0,0,0), pointmat_ownmem(false), higher_integration_order(false)
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
    switch (fel->SpatialDim())
      {
      case 1:
        dxdxi = pointmat * static_cast<const NodalFiniteElement<1>*> (fel)->GetDShape(ip, lh);
        break;
      case 2:
        dxdxi = pointmat * static_cast<const NodalFiniteElement<2>*> (fel)->GetDShape(ip, lh);
        break;
      case 3:
        dxdxi = pointmat * static_cast<const NodalFiniteElement<3>*> (fel)->GetDShape(ip, lh);
        break;
      }

  }

  ///
  template <typename T>
  void CalcJacobian (const IntegrationPoint & ip,
		     MatExpr<T> & dxdxi) const
  {
    LocalHeap lh(1000);
    switch (fel->SpatialDim())
      {
      case 1:
        dxdxi = pointmat * static_cast<const NodalFiniteElement<1>*> (fel)->GetDShape(ip, lh);
        break;
      case 2:
        dxdxi = pointmat * static_cast<const NodalFiniteElement<2>*> (fel)->GetDShape(ip, lh);
        break;
      case 3:
        dxdxi = pointmat * static_cast<const NodalFiniteElement<3>*> (fel)->GetDShape(ip, lh);
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
    switch (fel->SpatialDim())
      {
      case 1:
        point = pointmat * static_cast<const NodalFiniteElement<1>*> (fel)->GetShape(ip, lh);
        break;
      case 2:
        point = pointmat * static_cast<const NodalFiniteElement<2>*> (fel)->GetShape(ip, lh);
        break;
      case 3:
        point = pointmat * static_cast<const NodalFiniteElement<3>*> (fel)->GetShape(ip, lh);
        break;
      }

        // point = pointmat * fel->GetShape(ip, lh);
  }

  template <typename T1, typename T2>
  void CalcPointJacobian (const IntegrationPoint & ip,
			  T1 & point, T2 & dxdxi,
			  LocalHeap & lh) const
  {
    switch (fel->SpatialDim())
      {
      case 1:
        dxdxi = pointmat * static_cast<const NodalFiniteElement<1>*> (fel)->GetDShape(ip, lh);
        point = pointmat * static_cast<const NodalFiniteElement<1>*> (fel)->GetShape(ip, lh);
        break;
      case 2:
        dxdxi = pointmat * static_cast<const NodalFiniteElement<2>*> (fel)->GetDShape(ip, lh);
        point = pointmat * static_cast<const NodalFiniteElement<2>*> (fel)->GetShape(ip, lh);
        break;
      case 3:
        dxdxi = pointmat * static_cast<const NodalFiniteElement<3>*> (fel)->GetDShape(ip, lh);
        point = pointmat * static_cast<const NodalFiniteElement<3>*> (fel)->GetShape(ip, lh);
        break;
      }
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
    return pointmat.Height() != fel->SpatialDim();
  }

  void GetSort (FlatArray<int> sort) const
  { ; }
};




#else




/**
    Transformation from reference element to physical element.
    Uses Netgen-meshaccess to perform transformation
*/
class ElementTransformation
{
  /// element number
  int elnr;

  /// material property
  int elindex;
  
  /// elemente in R^dim
  int dim;
  /// is it a boundary element ?
  bool boundary;
  /// geometry of element
  ELEMENT_TYPE eltype;

  const ARRAY<Vec<3> > * pts;
  const ARRAY<Mat<3,3> > * dxdxis;
  const ARRAY<int> * first_of_element;

  bool higher_integration_order;

  bool iscurved;

  double buffer[100];
  bool buffervalid;

public:
  ///
  ElementTransformation () { pts = 0; higher_integration_order = false; buffervalid = false; }
  /*
  ElementTransformation (ARRAY<Vec<3> > * apts,
			 ARRAY<Mat<3,3> > * adxdxis,
			 ARRAY<int> * afirst_of_element)
    : pts(apts), dxdxis(adxdxis), first_of_element(afirst_of_element) { ; }
  */

  void SetGeometryData (const ARRAY<Vec<3> > * apts,
			const ARRAY<Mat<3,3> > * adxdxis,
			const ARRAY<int> * afirst_of_element)
  {
    pts = apts;
    dxdxis = adxdxis;
    first_of_element = afirst_of_element;
  }

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
    buffervalid = 0;
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
      Ng_GetSurfaceElementTransformation (elnr+1, &ip(0),
					  0, &dxdxi(0,0));
    else
      {
        Ng_GetBufferedElementTransformation (elnr+1, &ip(0),
                                             0, &dxdxi(0,0), 
                                             const_cast<double*> (&buffer[0]), buffervalid);
        const_cast<bool&> (buffervalid) = true;
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
        Ng_GetBufferedElementTransformation (elnr+1, &ip(0),
                                             0, &dxdxi(0,0), 
                                             const_cast<double*> (&buffer[0]), buffervalid);
        const_cast<bool&> (buffervalid) = true;
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
        Ng_GetBufferedElementTransformation (elnr+1, &ip(0), &point(0), 0, 
                                             const_cast<double*> (&buffer[0]), buffervalid);
        const_cast<bool&> (buffervalid) = true;
      }
  }



  template <typename T1, typename T2>
  void CalcPointJacobian (const IntegrationPoint & ip,
			  T1 & point, T2 & dxdxi,
			  LocalHeap & lh) const
  {
    if (boundary)
      Ng_GetSurfaceElementTransformation (elnr+1, &ip(0), &point(0), &dxdxi(0,0));
    else
      {
	if (ip.precomputed_geometry && pts)
	  {
	    point = (*pts)[ (*first_of_element)[elnr] + ip.Nr()];
	    dxdxi = (*dxdxis)[ (*first_of_element)[elnr] + ip.Nr()];
	  }
	else
          {
            Ng_GetBufferedElementTransformation (elnr+1, &ip(0), &point(0), &dxdxi(0,0),
                                                 const_cast<double*> (&buffer[0]), buffervalid);
            const_cast<bool&> (buffervalid) = true;
          }
      }
  }


  template <typename T0, typename T1, typename T2>
  void CalcMultiPointJacobian (FlatArray<T0> ipts,
                               FlatArray<T1> point, 
                               FlatArray<T2> dxdxi,
                               LocalHeap & lh) const
  {
    if (boundary)
      ;
      // Ng_GetSurfaceElementTransformation (elnr+1, &ip(0), &point(0), &dxdxi(0,0));
    else
      {
        /*
        for (int i = 0; i < ipts.Size(); i++)
          Ng_GetElementTransformation (elnr+1, &ipts[i](0), &point[i](0), &dxdxi[i](0,0));
        */

        Ng_GetMultiElementTransformation (elnr+1, ipts.Size(),
                                          &ipts[0](0), &ipts[1](0)-&ipts[0](0),
                                          &point[0](0), &point[1](0)-&point[0](0),
                                          &dxdxi[0](0,0), &dxdxi[1](0,0)-&dxdxi[0](0,0));
      }
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
      }
  }
};



#endif
#endif






