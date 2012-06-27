#ifndef FILE_INTRULE
#define FILE_INTRULE

/*********************************************************************/
/* File:   intrule.hpp                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

namespace ngfem
{

  /// An integration point 
  class IntegrationPoint
  {
  private:
    /// number within intergration Rule
    int nr; 
    /// coordinates (empty values for 1D and 2D)
    double pi[3];
    /// weight of integration point
    double weight;
    /// point is on facetnr, -1 for volume
    int facetnr;

  public:
    ///
    bool precomputed_geometry;
  public:
    ///
    IntegrationPoint & operator=(const IntegrationPoint & aip)
    {
      nr = aip.Nr();
      pi[0] = aip(0);
      pi[1] = aip(1);
      pi[2] = aip(2);
      weight = aip.Weight();
      precomputed_geometry = aip.precomputed_geometry;
      facetnr = -1;
      return *this;
    }

    ///
    IntegrationPoint (const double api[3], double aw)
    {
      nr = -1;
      pi[0] = api[0];
      pi[1] = api[1];
      pi[2] = api[2];
      weight = aw;
      facetnr = -1;
      precomputed_geometry = 0;
    }

    ///
    IntegrationPoint (double p1 = 0, double p2 = 0, double p3 = 0, double aw = 0)
    {
      nr = -1;
      pi[0] = p1;
      pi[1] = p2;
      pi[2] = p3;
      weight = aw;
      facetnr = -1;
      precomputed_geometry = 0;
    }

    ///
    IntegrationPoint (const FlatVector<double> & ap, double aw)
    {
      nr = -1;
      pi[0] = (ap.Size() >= 1) ? ap(0) : 0;
      pi[1] = (ap.Size() >= 2) ? ap(1) : 0;
      pi[2] = (ap.Size() >= 3) ? ap(2) : 0;
      weight = aw;
      facetnr = -1;
      precomputed_geometry = 0;
    }


    template <int D>
    IntegrationPoint (const Vec<D> & ap, double aw = -1)
    {
      nr = -1;
      for (int j = 0; j < D; j++)
	pi[j] = ap(j);
      weight = aw;
      facetnr = -1;
      precomputed_geometry = 0;
    }

    ///
    IntegrationPoint (const IntegrationPoint & aip)
    { *this = aip; }

    ///
    void SetNr (int anr) { nr = anr; }
    ///
    const double * Point () const { return pi; }
    ///
    int Size() const { return 3; }
    ///
    const double & operator() (int i) const { return pi[i]; }
    ///
    double & operator() (int i) { return pi[i]; }
    ///
    double Weight () const { return weight; }
    ///
    void SetWeight (double w) { weight = w; }

    /// number within integration rule
    int Nr () const { return nr; }

    /// global number of ip
    int IPNr () const { return -1; }
  public:
    int & FacetNr() { return facetnr; }
    int FacetNr() const { return facetnr; }

    ///
    friend NGS_DLL_HEADER ostream & operator<< (ostream & ost, const IntegrationPoint & ip);
  };



  class NGS_DLL_HEADER ElementTransformation;

  /**
     Base class for MappedIntegrationPoint.
     A specific integration point is the mapped point, and stores
     point coordinates, the Jacobimatrix, and the determinant of the Jacobimatrix.
  */
  class BaseMappedIntegrationPoint
  {
  protected:
    /// IP on the reference element
    const IntegrationPoint * ip;
    /// computed by the transformation
    const ElementTransformation * eltrans;
    /// fabs(det)
    double measure; 
  public:
    ///
    BaseMappedIntegrationPoint (const IntegrationPoint & aip,
				const ElementTransformation & aeltrans)
      : ip(&aip), eltrans(&aeltrans)  { ; }
    /// 
    const IntegrationPoint & IP () const { return *ip; }
    ///
    const ElementTransformation & GetTransformation () const { return *eltrans; }
    ///
    int GetIPNr() const { return ip->Nr(); }
    ///
    double GetMeasure() const { return measure; }
    ///
    double GetWeight() const { return measure * ip->Weight(); }
  };


  template <int R, typename SCAL = double>
  class DimMappedIntegrationPoint : public BaseMappedIntegrationPoint
  {
  protected:
    ///
    Vec<R,SCAL> point;
  public:
    ///
    DimMappedIntegrationPoint (const IntegrationPoint & aip,
				 const ElementTransformation & aeltrans)
      : BaseMappedIntegrationPoint (aip, aeltrans)
    { ; }
    ///
    const Vec<R,SCAL> & GetPoint () const { return point; }
    Vec<R,SCAL> & Point () { return point; }
    ///
    SCAL operator() (int i) const { return point(i); }
  };


  /// ip, dimension source, dimension range
  template <int DIMS, int DIMR, typename SCAL = double> 
  class MappedIntegrationPoint : public DimMappedIntegrationPoint<DIMR,SCAL>
  {
  private:
    /// Jacobi matrix
    Mat<DIMR,DIMS,SCAL> dxdxi;
    /// (pseudo)inverse of Jacobi matrix
    Mat<DIMS,DIMR,SCAL> dxidx;
    /// Jacobian 
    SCAL det;
    /// for boundary points
    Vec<DIMR,SCAL> normalvec;
    Vec<DIMR,SCAL> tangentialvec;  // for independent integrator
  
 
  public:
    typedef SCAL TSCAL;
    ///
    NGS_DLL_HEADER MappedIntegrationPoint (const IntegrationPoint & aip,
					   const ElementTransformation & aeltrans);

    MappedIntegrationPoint (const IntegrationPoint & aip,
			    const ElementTransformation & aeltrans,
			    int /* dummy */)
      : DimMappedIntegrationPoint<DIMR,SCAL> (aip, aeltrans)
    { ; }

    ///
    MappedIntegrationPoint (const IntegrationPoint & aip,
			    const ElementTransformation & aeltrans,
			    const FlatVec<DIMR, SCAL> ax,
			    const Mat<DIMR, DIMS, SCAL> & adxdxi)
      : DimMappedIntegrationPoint<DIMR,SCAL> (aip, aeltrans)
    {
      this->point = ax;
      dxdxi = adxdxi;
      Compute();
    }

    void Compute ()
    {
      if (DIMS == DIMR)
	{
	  det = Det (dxdxi);
	  dxidx = Inv (dxdxi);
	  normalvec = TSCAL(0.0);
	  tangentialvec = TSCAL(0.0);
	}
      else
	{
	  if (DIMR == 3)
	    {
	      normalvec = Cross (Vec<3,SCAL> (dxdxi.Col(0)),
				 Vec<3,SCAL> (dxdxi.Col(1)));
	      det = L2Norm (normalvec);
	      normalvec /= det;
	    }
	  else if (DIMR == 2)
	    {
	      det = sqrt ( sqr (dxdxi(0,0)) + sqr (dxdxi(1,0)));

	      normalvec(0) = -dxdxi(1,0) / det;
	      normalvec(1) = dxdxi(0,0) / det;
	    }
	  else
	    {
	      det = 1.0;
	      normalvec = 1.0;
	    }
	
	  Mat<DIMS,DIMS,SCAL> ata, iata;
	
	  ata = Trans (dxdxi) * dxdxi;
	  iata = Inv (ata);
	  dxidx = iata * Trans (dxdxi);
	  tangentialvec = TSCAL(0.0);
	}
      this->measure = fabs (det);
    }
  
    ///
    const Mat<DIMR,DIMS,SCAL> & GetJacobian() const { return dxdxi; }
    Mat<DIMR,DIMS,SCAL> & Jacobian() { return dxdxi; }
    ///
    SCAL GetJacobiDet() const { return det; }
    ///
    const Mat<DIMS,DIMR,SCAL> & GetJacobianInverse() const { return dxidx; }
    ///
    const Vec<DIMR,SCAL> GetNV () const { return normalvec; }
    /// 
    void SetNV ( const Vec<DIMR,SCAL> & vec) { normalvec = vec; }
    ///
    void SetTV ( const Vec<DIMR,SCAL> & vec) { tangentialvec = vec; }
    ///
    const Vec<DIMR,SCAL> GetTV () const { return tangentialvec; }
    ///
    int IsBoundary () const { return DIMS != DIMR; }

    ///
    void CalcHesse (Mat<2> & ddx1, Mat<2> & ddx2) const;
    void CalcHesse (Mat<2> & ddx1, Mat<2> & ddx2, Mat<2> & ddx3) const; 
    void CalcHesse (Mat<3> & ddx1, Mat<3> & ddx2, Mat<3> & ddx3) const;
  };


  template <int DIMS, int DIMR, typename SCAL> 
  inline ostream & operator<< (ostream & ost, const MappedIntegrationPoint<DIMS,DIMR,SCAL> & sip)
  {
    ost << sip.GetPoint() << ", dxdxi = " << sip.GetJacobian();
    return ost;
  }





  /**
     An integration rule.
     Contains array of integration points
  */

  class IntegrationRule : public Array<IntegrationPoint> 
  {
  public:
    ///
    IntegrationRule () { ; }

    /**
       An integration rule for element type of certain order.
       Obtains a reference to the precomputed integration rule
    */
    NGS_DLL_HEADER IntegrationRule (ELEMENT_TYPE eltype, int order);


    NGS_DLL_HEADER IntegrationRule (const IntegrationRule & ir2)
    {
      throw Exception ("should not call copy-constructor of ir\n");
    }

    NGS_DLL_HEADER IntegrationRule (int asize, LocalHeap & lh)
      : Array<IntegrationPoint> (asize, lh)
    {
      ;
    }

    NGS_DLL_HEADER IntegrationRule (int asize, double (*pts)[3], double * weights);

    // make it polymorphic
    virtual ~IntegrationRule() { ; }

    ///
    void AddIntegrationPoint (const IntegrationPoint & ip)
    { 
      Append (ip);
    }

    /// number of integration points
    int GetNIP() const { return Size(); }
  };

  inline ostream & operator<< (ostream & ost, const IntegrationRule & ir)
  {
    for (int i = 0; i < ir.GetNIP(); i++)
      ost << ir[i] << endl;
    return ost;
  }




  template <int D>
  class IntegrationRuleTP : public IntegrationRule
  {
    const IntegrationRule *irx, *iry, *irz;
    // int sort[8];

    // ArrayMem<Vec<D>, 100> xi;
    // ArrayMem<double, 100> weight;
    ArrayMem<Vec<D>, 100> x;
    ArrayMem<Mat<D,D>, 100> dxdxi;
    ArrayMem<Mat<D,D>, 100> dxdxi_duffy;

  public:
    NGS_DLL_HEADER IntegrationRuleTP (const ElementTransformation & eltrans,
		       int order, bool compute_mapping, LocalHeap & lh);

    // tensor product rule for a facet
    NGS_DLL_HEADER IntegrationRuleTP (ELEMENT_TYPE eltype, FlatArray<int> sort, 
		       NODE_TYPE nt, int nodenr, int order, LocalHeap & lh);

    const IntegrationRule & GetIRX() const { return *irx; }
    const IntegrationRule & GetIRY() const { return *iry; }
    const IntegrationRule & GetIRZ() const { return *irz; }

    int GetNIP() const { return Size(); }
    double GetWeight (int i) const { return (*this)[i].Weight(); }  // weight[i]; }
    const Vec<D> GetXi(int i) const 
    { 
      Vec<D> ret;
      for (int j = 0; j < D; j++)
	ret(j) = (*this)[i](j);
      return ret;  // return xi[i]; 
    }
    Vec<D> & GetPoint(int i) const { return x[i]; }
    Mat<D,D> & GetJacobian(int i) const { return dxdxi[i]; }
    Mat<D,D> & GetDuffyJacobian(int i) const { return dxdxi_duffy[i]; }
  };




  /**
     Computes Gaussean integration.
     Integration rule on (0,1), contains n points,
     Exact for polynomials up to order 2n-1
  */ 
  NGS_DLL_HEADER extern void ComputeGaussRule (int n, 
				Array<double> & xi, 
				Array<double> & wi);




  /**
     Computes Gauss-Jacobi integration rule.
     Integration rule on (0,1), contains n points.
     Assumes that the polynomial has alpha roots in 1, and beta roots in 0.
     Exact for polynomials up to order 2n-1 + alpha + beta
  */ 
  NGS_DLL_HEADER extern void ComputeGaussJacobiRule (int n, 
				      Array<double> & xi, 
				      Array<double> & wi,
				      double alf,
				      double bet);


  /**
     Computes Gauss-Hermite integration rule.
     Integration rule on R, contains n points.
     Exact for exp{-x*x} * p(x),  with polynomials p(x) up to order 2n-1
  */ 
  NGS_DLL_HEADER extern void ComputeHermiteRule (int n, 
				  Array<double> & x,
				  Array<double> & w);
  



  /// Get a reference to the integration-rules container
  // extern NGS_DLL_HEADER const IntegrationRules & GetIntegrationRules ();

  extern NGS_DLL_HEADER const IntegrationRule & SelectIntegrationRule (ELEMENT_TYPE eltype, int order);
  extern NGS_DLL_HEADER const IntegrationRule & SelectIntegrationRuleJacobi10 (int order);
  extern NGS_DLL_HEADER const IntegrationRule & SelectIntegrationRuleJacobi20 (int order);


  // transformation of (d-1) dimensional integration points on facets to 
  // d-dimensional point in volumes
  class Facet2ElementTrafo
  {
  protected:
    // mutable IntegrationPoint elpoint;  
    ELEMENT_TYPE eltype;
    // const POINT3D * points;
    FlatVector<Vec<3> > points;
    const EDGE * edges;
    const FACE * faces;
    EDGE hedges[4];
    FACE hfaces[6];
  public:
    Facet2ElementTrafo(ELEMENT_TYPE aeltype) 
      : eltype(aeltype),
	points(99,(double*)ElementTopology::GetVertices (aeltype))
    {
      // points = ElementTopology::GetVertices (eltype);
      edges = ElementTopology::GetEdges (eltype);
      faces = ElementTopology::GetFaces (eltype);
    }
  
    Facet2ElementTrafo(ELEMENT_TYPE aeltype, const FlatArray<int> & vnums) 
      : eltype(aeltype),
	points(99,(double*)ElementTopology::GetVertices (aeltype))
    {
      // points = ElementTopology::GetVertices (eltype);
      edges = ElementTopology::GetEdges (eltype);
      faces = ElementTopology::GetFaces (eltype);

      if (eltype == ET_TRIG)
	{
	  for (int i = 0; i < 3; i++)
	    {
	      hedges[i][0] = edges[i][0];
	      hedges[i][1] = edges[i][1];
	      if (vnums[hedges[i][0]] > vnums[hedges[i][1]])
		swap (hedges[i][0], hedges[i][1]);
	    }
	  edges = &hedges[0];
	}

      if (eltype == ET_QUAD)
	{
	  for (int i = 0; i < 4; i++)
	    {
	      hedges[i][0] = edges[i][0];
	      hedges[i][1] = edges[i][1];
	      if (vnums[hedges[i][0]] > vnums[hedges[i][1]])
		swap (hedges[i][0], hedges[i][1]);
	    }
	  edges = &hedges[0];
	}

      if (eltype == ET_TET)
	{
	  for (int i = 0; i < 4; i++)
	    {
	      hfaces[i][0] = faces[i][0];
	      hfaces[i][1] = faces[i][1];
	      hfaces[i][2] = faces[i][2];
	      if (vnums[hfaces[i][0]] > vnums[hfaces[i][1]]) swap (hfaces[i][0], hfaces[i][1]);
	      if (vnums[hfaces[i][1]] > vnums[hfaces[i][2]]) swap (hfaces[i][1], hfaces[i][2]);
	      if (vnums[hfaces[i][0]] > vnums[hfaces[i][1]]) swap (hfaces[i][0], hfaces[i][1]);
	    }
	  faces = &hfaces[0];
	}


      if (eltype == ET_PRISM)
	{
	  for (int i = 0; i < 2; i++)
	    {
	      hfaces[i][0] = faces[i][0];
	      hfaces[i][1] = faces[i][1];
	      hfaces[i][2] = faces[i][2];
	      if (vnums[hfaces[i][0]] > vnums[hfaces[i][1]]) swap (hfaces[i][0], hfaces[i][1]);
	      if (vnums[hfaces[i][1]] > vnums[hfaces[i][2]]) swap (hfaces[i][1], hfaces[i][2]);
	      if (vnums[hfaces[i][0]] > vnums[hfaces[i][1]]) swap (hfaces[i][0], hfaces[i][1]);
	    }
	  for (int i = 2; i < 5; i++)
	    {
	      int jmin = 0;
	      for (int j = 1; j < 4; j++)
		if (vnums[faces[i][j]] < vnums[faces[i][jmin]]) jmin = j;
	      int j1 = (jmin+1)%4;
	      int j2 = (jmin+2)%4;
	      int j3 = (jmin+3)%4;
	      if (vnums[faces[i][j3]] < vnums[faces[i][j1]]) swap (j1, j3);

	      hfaces[i][0] = faces[i][jmin];
	      hfaces[i][1] = faces[i][j1];
	      hfaces[i][2] = faces[i][j2];
	      hfaces[i][3] = faces[i][j3];
	    }
	  faces = &hfaces[0];
	}

    }


    void operator()(int fnr, const IntegrationPoint &ipfac, IntegrationPoint & ipvol) const 
    {
      ELEMENT_TYPE facettype = ElementTopology::GetFacetType(eltype, fnr);

      switch (facettype)
	{
	case ET_SEGM:
	  {
	    FlatVec<3> p1 = points (edges[fnr][0]);
	    FlatVec<3> p2 = points (edges[fnr][1]);

	    // Vec<3> volp = p2 + ipfac(0) * (p1-p2);
	    // ipvol = volp;
	    ipvol = Vec<3> (p2 + ipfac(0) * (p1-p2));
	    break;
	  }
	case ET_TRIG:
	  {
	    FlatVec<3> p0 = points(faces[fnr][0]);
	    FlatVec<3> p1 = points(faces[fnr][1]);
	    FlatVec<3> p2 = points(faces[fnr][2]);

            // for (int j = 0; j < 3; j++)
	    //   ipvol(j) = p2[j] + ipfac(0)*(p0[j]-p2[j]) + ipfac(1)*(p1[j]-p2[j]);
	    ipvol = Vec<3> (p2 + ipfac(0) * (p0-p2) + ipfac(1)*(p1-p2));
	    break;
	  }
	case ET_QUAD:
	  {
	    FlatVec<3> p0 = points(faces[fnr][0]);
	    FlatVec<3> p1 = points(faces[fnr][1]);
	    FlatVec<3> p2 = points(faces[fnr][3]);

            // for (int j = 0; j < 3; j++)
	    //  ipvol(j) = p0[j] + ipfac(0)*(p1[j]-p0[j]) + ipfac(1)*(p2[j]-p0[j]);
	    ipvol = Vec<3> (p0 + ipfac(0) * (p1-p0) + ipfac(1)*(p2-p0));
	    break;
	  }
	default:
	  throw Exception ("undefined facet type in Facet2ElementTrafo()\n");
	} 
      ipvol.FacetNr() = fnr;
    }

    const IntegrationPoint operator()(int fnr, const IntegrationPoint &ip1d) const 
    {
      IntegrationPoint elpoint;
      operator()(fnr, ip1d, elpoint);
      return elpoint;
    }

    IntegrationRule & operator() (int fnr, const IntegrationRule & irfacet, LocalHeap & lh)
    {
      IntegrationRule & irvol = *new (lh) IntegrationRule (irfacet.GetNIP(), lh);

      switch (ElementTopology::GetFacetType(eltype, fnr))
	{
	case ET_SEGM:
	  {
	    FlatVec<3> p1 = points (edges[fnr][0]);
	    FlatVec<3> p2 = points (edges[fnr][1]);

	    for (int i = 0; i < irfacet.GetNIP(); i++)
	      irvol[i] = Vec<3> (p2 + irfacet[i](0) * (p1-p2));
	    break;
	  }
	case ET_TRIG:
	  {
	    FlatVec<3> p0 = points(faces[fnr][0]);
	    FlatVec<3> p1 = points(faces[fnr][1]);
	    FlatVec<3> p2 = points(faces[fnr][2]);

	    for (int i = 0; i < irfacet.GetNIP(); i++)
	      irvol[i] = Vec<3> (p2 + irfacet[i](0) * (p0-p2) + irfacet[i](1)*(p1-p2));
	    break;
	  }
	case ET_QUAD:
	  {
	    FlatVec<3> p0 = points(faces[fnr][0]);
	    FlatVec<3> p1 = points(faces[fnr][1]);
	    FlatVec<3> p2 = points(faces[fnr][3]);

	    for (int i = 0; i < irfacet.GetNIP(); i++)
	      irvol[i] = Vec<3> (p0 + irfacet[i](0) * (p1-p0) + irfacet[i](1)*(p2-p0));
	    break;
	  }
	default:
	  throw Exception ("undefined facet type in Facet2ElementTrafo()\n");
	} 

      /*
      for (int i = 0; i < irfacet.GetNIP(); i++)
	(*this) (fnr, irfacet[i], irvol[i]);
      */
      for (int i = 0; i < irfacet.GetNIP(); i++)
	irvol[i].FacetNr() = fnr;

      return irvol;
    }
  };







  // transformation of (d-1) dimensional integration points on facets to 
  // d-dimensional point in volumes
  class Facet2SurfaceElementTrafo
  {
  protected:
    // mutable IntegrationPoint elpoint;  
    ELEMENT_TYPE eltype;
    const POINT3D * points;
    const EDGE * edges;
    const FACE * faces;
    EDGE hedges[4];
    FACE hfaces[6];
  public:
    Facet2SurfaceElementTrafo(ELEMENT_TYPE aeltype) : eltype(aeltype) 
    {
      points = ElementTopology::GetVertices (eltype);
      edges = ElementTopology::GetEdges (eltype);
      faces = ElementTopology::GetFaces (eltype);
    }
    
    Facet2SurfaceElementTrafo(ELEMENT_TYPE aeltype, FlatArray<int> & vnums) : eltype(aeltype) 
    {
      points = ElementTopology::GetVertices (eltype);
      edges = ElementTopology::GetEdges (eltype);
      faces = ElementTopology::GetFaces (eltype);

      if (eltype == ET_SEGM)
	{
	  hedges[0][0] = edges[0][0];
	  hedges[0][1] = edges[0][1];
	  if (vnums[hedges[0][0]] > vnums[hedges[0][1]])
	    swap (hedges[0][0], hedges[0][1]);
	  edges = &hedges[0];
	}

      if (eltype == ET_TRIG)
	{
	  hfaces[0][0] = faces[0][0];
	  hfaces[0][1] = faces[0][1];
	  hfaces[0][2] = faces[0][2];
	  if (vnums[hfaces[0][0]] > vnums[hfaces[0][1]]) swap (hfaces[0][0], hfaces[0][1]);
	  if (vnums[hfaces[0][1]] > vnums[hfaces[0][2]]) swap (hfaces[0][1], hfaces[0][2]);
	  if (vnums[hfaces[0][0]] > vnums[hfaces[0][1]]) swap (hfaces[0][0], hfaces[0][1]);
	  
	  faces = &hfaces[0];
	}

      if (eltype == ET_QUAD)
	{
	  int jmin = 0;
	  for (int j = 1; j < 4; j++)
	    if (vnums[faces[0][j]] < vnums[faces[0][jmin]]) jmin = j;
	  int j1 = (jmin+1)%4;
	  int j2 = (jmin+2)%4;
	  int j3 = (jmin+3)%4;
	  if (vnums[faces[0][j3]] < vnums[faces[0][j1]]) swap (j1, j3);
	  
	  hfaces[0][0] = faces[0][jmin];
	  hfaces[0][1] = faces[0][j1];
	  hfaces[0][2] = faces[0][j2];
	  hfaces[0][3] = faces[0][j3];
	  faces = &hfaces[0];
	}
    }


    void operator()(const IntegrationPoint &ipfac, IntegrationPoint & ipvol) const 
    {
      int fnr = 0;
      switch (eltype)
	{
	case ET_SEGM:
	  {
	    const POINT3D & p1 = points[edges[fnr][0]];
	    const POINT3D & p2 = points[edges[fnr][1]];
            for (int j = 0; j < 3; j++)
              ipvol(j) = p2[j] + (ipfac(0))*(p1[j]-p2[j]);
	    break;
	  }
	case ET_TRIG:
	  {
	    const POINT3D & p0 = points[faces[fnr][0]];
	    const POINT3D & p1 = points[faces[fnr][1]];
	    const POINT3D & p2 = points[faces[fnr][2]];
	    ipvol(0) = p2[0] + ipfac(0)*(p0[0]-p2[0]) + ipfac(1)*(p1[0]-p2[0]);
	    ipvol(1) = p2[1] + ipfac(0)*(p0[1]-p2[1]) + ipfac(1)*(p1[1]-p2[1]);
	    ipvol(2) = p2[2] + ipfac(0)*(p0[2]-p2[2]) + ipfac(1)*(p1[2]-p2[2]);
	    break;
	  }
	case ET_QUAD:
	  {
	    const POINT3D & p0 = points[faces[fnr][0]];
	    const POINT3D & p1 = points[faces[fnr][1]];
	    const POINT3D & p2 = points[faces[fnr][3]];
	    ipvol(0) = p0[0] + ipfac(0)*(p1[0]-p0[0]) + ipfac(1)*(p2[0]-p0[0]);
	    ipvol(1) = p0[1] + ipfac(0)*(p1[1]-p0[1]) + ipfac(1)*(p2[1]-p0[1]);
	    ipvol(2) = p0[2] + ipfac(0)*(p1[2]-p0[2]) + ipfac(1)*(p2[2]-p0[2]);
	    break;
	  }
	default:
	  throw Exception ("undefined facet type in Facet2ElementTrafo()\n");

	} 
      /*      cerr << "*** mapping integrationpoint for element " << eltype << " and facel " << fnr << " of type " << facettype << endl;
	      cerr << "  * ipfac = " << ipfac;
	      cerr << "  * ipvol = " << ipvol;*/
    }
    const IntegrationPoint operator()(const IntegrationPoint &ip1d) const 
    {
      IntegrationPoint elpoint;
      operator()(ip1d, elpoint);
      return elpoint;
    }
  };













  class NGS_DLL_HEADER BaseMappedIntegrationRule 
  { 
  protected:
    const IntegrationRule & ir;
    const ElementTransformation & eltrans;
    char * baseip;
    int incr;

  public:
    BaseMappedIntegrationRule (const IntegrationRule & air,
			       const ElementTransformation & aeltrans)
      : ir(air), eltrans(aeltrans) { ; }

    int Size() const { return ir.Size(); }
    const IntegrationRule & IR() const { return ir; }
    const ElementTransformation & GetTransformation () const { return eltrans; }

    const BaseMappedIntegrationPoint & operator[] (int i) const
    { return *static_cast<const BaseMappedIntegrationPoint*> ((void*)(baseip+i*incr)); }
  };

  template <int DIM_ELEMENT, int DIM_SPACE>
  class NGS_DLL_HEADER MappedIntegrationRule : public BaseMappedIntegrationRule
  {
    FlatArray< MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE> > sips;
  public:
    MappedIntegrationRule (const IntegrationRule & ir, 
			   const ElementTransformation & aeltrans, 
			   LocalHeap & lh);

    MappedIntegrationRule (const IntegrationRule & ir, 
			   const ElementTransformation & eltrans, 
			   int dummy,
			   LocalHeap & lh)
      : BaseMappedIntegrationRule (ir, eltrans), sips(ir.GetNIP(), lh)
    {
      baseip = (char*)(void*)(BaseMappedIntegrationPoint*)(&sips[0]);
      incr = (char*)(void*)(&sips[1]) - (char*)(void*)(&sips[0]);
    }
    
    MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE> & operator[] (int i) const
    { 
      return sips[i]; 
    }
  };

#define SpecificIntegrationPoint MappedIntegrationPoint 
  /*
#define BaseSpecificIntegrationPoint BaseMappedIntegrationPoint 
#define DimSpecificIntegrationPoint DimMappedIntegrationPoint 
  */



}


#endif
