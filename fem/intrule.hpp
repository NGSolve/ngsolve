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
  class /* alignas(8) */ IntegrationPoint
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
    INLINE IntegrationPoint & operator=(const IntegrationPoint & aip)
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
    INLINE IntegrationPoint (const double api[3], double aw)
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
    INLINE IntegrationPoint (double p1 = 0, double p2 = 0, double p3 = 0, double aw = 0)
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
    INLINE IntegrationPoint (const FlatVector<double> & ap, double aw)
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
    INLINE IntegrationPoint (const Vec<D> & ap, double aw = -1)
    {
      nr = -1;
      for (int j = 0; j < D; j++)
	pi[j] = ap(j);
      weight = aw;
      facetnr = -1;
      precomputed_geometry = 0;
    }

    ///
    INLINE IntegrationPoint (const IntegrationPoint & aip)
    { *this = aip; }

    ///
    INLINE void SetNr (int anr) { nr = anr; }
    ///
    // const double * Point () const { return pi; }
    ///
    INLINE int Size() const { return 3; }
    ///
    INLINE const double & operator() (int i) const { return pi[i]; }
    ///
    INLINE double & operator() (int i) { return pi[i]; }
    ///
    INLINE double Weight () const { return weight; }
    ///
    INLINE void SetWeight (double w) { weight = w; }

    /// number within integration rule
    INLINE int Nr () const { return nr; }

    /// global number of ip
    INLINE int IPNr () const { return -1; }

    INLINE FlatVec<3, double> Point() { return &pi[0]; }
    INLINE FlatVec<3, const double> Point() const { return &pi[0]; }

    INLINE int & FacetNr() { return facetnr; }
    INLINE int FacetNr() const { return facetnr; }


    template <int DIM> 
    INLINE operator Vec<DIM, AutoDiff<DIM>> () const
    {
      Vec<DIM, AutoDiff<DIM> > adp;
      for (int i = 0; i < DIM; i++)
        adp[i] = AutoDiff<DIM> (pi[i], i);
      return adp;
    }
    


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
    INLINE BaseMappedIntegrationPoint (const IntegrationPoint & aip,
				const ElementTransformation & aeltrans)
      : ip(&aip), eltrans(&aeltrans)  { ; }
    /// 
    INLINE const IntegrationPoint & IP () const { return *ip; }
    ///
    INLINE const ElementTransformation & GetTransformation () const { return *eltrans; }
    ///
    INLINE int GetIPNr() const { return ip->Nr(); }
    ///
    INLINE double GetMeasure() const { return measure; }
    ///
    INLINE double GetWeight() const { return measure * ip->Weight(); }
  };


  template <int R, typename SCAL = double>
  class DimMappedIntegrationPoint : public BaseMappedIntegrationPoint
  {
  protected:
    ///
    Vec<R,SCAL> point;
  public:
    ///
    INLINE DimMappedIntegrationPoint (const IntegrationPoint & aip,
				 const ElementTransformation & aeltrans)
      : BaseMappedIntegrationPoint (aip, aeltrans)
    { ; }
    ///
    INLINE const Vec<R,SCAL> & GetPoint () const { return point; }
    INLINE Vec<R,SCAL> & Point () { return point; }
    ///
    INLINE SCAL operator() (int i) const { return point(i); }
  };


  /// ip, dimension source, dimension range
  template <int DIMS, int DIMR, typename SCAL = double> 
  class MappedIntegrationPoint : public DimMappedIntegrationPoint<DIMR,SCAL>
  {
  private:
    /// Jacobi matrix
    Mat<DIMR,DIMS,SCAL> dxdxi;
    /// (pseudo)inverse of Jacobi matrix
    SCAL det;
    /// for boundary points
    Vec<DIMR,SCAL> normalvec;
    Vec<DIMR,SCAL> tangentialvec;  // for independent integrator
  
 
  public:
    typedef SCAL TSCAL;
    ///
    NGS_DLL_HEADER MappedIntegrationPoint (const IntegrationPoint & aip,
					   const ElementTransformation & aeltrans);

    INLINE MappedIntegrationPoint (const IntegrationPoint & aip,
			    const ElementTransformation & aeltrans,
			    int /* dummy */)
      : DimMappedIntegrationPoint<DIMR,SCAL> (aip, aeltrans)
    { ; }

    ///
    INLINE MappedIntegrationPoint (const IntegrationPoint & aip,
			    const ElementTransformation & aeltrans,
			    const FlatVec<DIMR, SCAL> ax,
			    const Mat<DIMR, DIMS, SCAL> & adxdxi)
      : DimMappedIntegrationPoint<DIMR,SCAL> (aip, aeltrans)
    {
      this->point = ax;
      dxdxi = adxdxi;
      Compute();
    }

    INLINE void Compute ()
    {
      if (DIMS == DIMR)
	{
	  det = Det (dxdxi);
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
	  tangentialvec = TSCAL(0.0);
	}
      this->measure = fabs (det);
    }
  
    ///
    INLINE const Mat<DIMR,DIMS,SCAL> & GetJacobian() const { return dxdxi; }
    INLINE Mat<DIMR,DIMS,SCAL> & Jacobian() { return dxdxi; }
    ///
    INLINE SCAL GetJacobiDet() const { return det; }
    ///
    // const Mat<DIMS,DIMR,SCAL> & GetJacobianInverse() const { return dxidx; }
    INLINE const Mat<DIMS,DIMR,SCAL> GetJacobianInverse() const 
    { 
      if (DIMS == DIMR)
        // return Inv (dxdxi);
        return 1.0/det * Trans (Cof (dxdxi));
      else
        {
	  Mat<DIMS,DIMS,SCAL> ata, iata;
	  ata = Trans (dxdxi) * dxdxi;
	  iata = Inv (ata);
	  return (iata * Trans (dxdxi));
        }
    }


    INLINE operator Vec<DIMS, AutoDiff<DIMR,TSCAL>> () const
    {
      Vec<DIMS, AutoDiff<DIMR, TSCAL> > adp;
      
      /*
      for (int i = 0; i < DIMS; i++)
        adp[i].Value() = this->IP()(i);
      for (int i = 0; i < DIMS; i++)
        for (int j = 0; j < DIMR; j++)
          adp[i].DValue(j) = GetJacobianInverse()(i,j);
      */
      Mat<DIMS,DIMR,TSCAL> ijac = GetJacobianInverse();
      for (int i = 0; i < DIMS; i++)
        adp[i] = AutoDiff<DIMR,TSCAL> (this->IP()(i), &ijac(i,0));
      return adp;
    }



    ///
    INLINE const Vec<DIMR,SCAL> GetNV () const { return normalvec; }
    /// 
    INLINE void SetNV ( const Vec<DIMR,SCAL> & vec) { normalvec = vec; }
    ///
    INLINE void SetTV ( const Vec<DIMR,SCAL> & vec) { tangentialvec = vec; }
    ///
    INLINE const Vec<DIMR,SCAL> GetTV () const { return tangentialvec; }
    ///
    INLINE int IsBoundary () const { return DIMS != DIMR; }

    ///
    void CalcHesse (Mat<2> & ddx1, Mat<2> & ddx2) const;
    void CalcHesse (Mat<2> & ddx1, Mat<2> & ddx2, Mat<2> & ddx3) const; 
    void CalcHesse (Mat<3> & ddx1, Mat<3> & ddx2, Mat<3> & ddx3) const;
  };


  /*
  template <int DIM> 
  INLINE Vec<DIM, AutoDiff<DIM>> Mip2Ad (const MappedIntegrationPoint<DIM,DIM> & mip)
  {
    Vec<DIM, AutoDiff<DIM> > adp;

    for (int i = 0; i < DIM; i++)
      adp[i].Value() = mip.IP()(i);
    for (int i = 0; i < DIM; i++)
      for (int j = 0; j < DIM; j++)
        adp[i].DValue(j) = mip.GetJacobianInverse()(i,j);
    
    return adp;
  }
  */


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
    INLINE IntegrationRule () { ; }

    /**
       An integration rule for element type of certain order.
       Obtains a reference to the precomputed integration rule
    */
    NGS_DLL_HEADER IntegrationRule (ELEMENT_TYPE eltype, int order);
    
    INLINE IntegrationRule (int asize, IntegrationPoint * pip)
      : Array<IntegrationPoint> (asize, pip) { ; }


    INLINE NGS_DLL_HEADER IntegrationRule (const IntegrationRule & ir2)
      : Array<IntegrationPoint> (ir2.Size(), &ir2[0])
    { ; }

    INLINE NGS_DLL_HEADER IntegrationRule (int asize, LocalHeap & lh)
      : Array<IntegrationPoint> (asize, lh)
    { ; }

    INLINE NGS_DLL_HEADER IntegrationRule (int asize, double (*pts)[3], double * weights);

    // make it polymorphic
    HD virtual ~IntegrationRule() { ; }

    ///
    INLINE void AddIntegrationPoint (const IntegrationPoint & ip)
    { 
      Append (ip);
    }

    /// number of integration points
    INLINE int GetNIP() const { return Size(); }

    INLINE IntegrationRule Range (int first, int next) const
    {
      return IntegrationRule (next-first, &(*this)[first]);
    }
  };

  DLL_HEADER ostream & operator<< (ostream & ost, const IntegrationRule & ir);

  /*
    DG Integration rule:
    volume IR contains points of facet IR.
    facet IR are Gauss-rules
    boundary_volume_factor is maximal ration of boundary-weight to volume-weight
   */
  class DGIntegrationRule : public IntegrationRule
  {
    Array<IntegrationRule*> facetrules;
    double boundary_volume_factor;
  public:
    NGS_DLL_HEADER DGIntegrationRule (ELEMENT_TYPE eltype, int order);
    int GetNFacets () const { return facetrules.Size(); }
    const IntegrationRule & GetFacetIntegrationRule (int fnr) const 
    { return *facetrules[fnr]; }
    double BoundaryVolumeFactor () const { return boundary_volume_factor; }
  };
  
  DLL_HEADER ostream & operator<< (ostream & ost, const DGIntegrationRule & ir);

  template <int D>
  class IntegrationRuleTP : public IntegrationRule
  {
    const IntegrationRule *irx, *iry, *irz;

    // ArrayMem<Vec<D>, 100> x;
    // ArrayMem<Mat<D,D>, 100> dxdxi;
    ArrayMem<Mat<D,D>, 100> dxdxi_duffy;
    Mat<D,D> dxdxi_permute;

  public:
    NGS_DLL_HEADER IntegrationRuleTP (const ElementTransformation & eltrans,
                                      int order, LocalHeap * lh = NULL);

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
    // Vec<D> & GetPoint(int i) const { return x[i]; }
    // Mat<D,D> & GetJacobian(int i) const { return dxdxi[i]; }
    Mat<D,D> & GetDuffyJacobian(int i) const { return dxdxi_duffy[i]; }
    Mat<D,D> GetPermutationJacobian() const { return dxdxi_permute; }
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
     Gauss-Lobatto Rule with n points on interval [0,1].
     Contains 0 and 1 as points
     exact for polynomials up to order 2n-3
   */
  NGS_DLL_HEADER extern 
  void ComputeGaussLobattoRule (int n, Array<double>& xi, Array<double>& wi);


  /**
     Computes Gauss-Hermite integration rule.
     Integration rule on R, contains n points.
     Exact for exp{-x*x} * p(x),  with polynomials p(x) up to order 2n-1
  */ 
  NGS_DLL_HEADER extern void ComputeHermiteRule (int n, 
				  Array<double> & x,
				  Array<double> & w);
  

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

      if (eltype == ET_HEX)
        {
          for (int i = 0; i < 6; i++)
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
	case ET_POINT:
	  {
	    ipvol = Vec<3> (points (fnr) );
	    break;
	  }

	case ET_SEGM:
	  {
	    FlatVec<3> p1 = points (edges[fnr][0]);
	    FlatVec<3> p2 = points (edges[fnr][1]);

	    ipvol = Vec<3> (p2 + ipfac(0) * (p1-p2));
	    break;
	  }
	case ET_TRIG:
	  {
	    FlatVec<3> p0 = points(faces[fnr][0]);
	    FlatVec<3> p1 = points(faces[fnr][1]);
	    FlatVec<3> p2 = points(faces[fnr][2]);

	    ipvol = Vec<3> (p2 + ipfac(0) * (p0-p2) + ipfac(1)*(p1-p2));
	    break;
	  }
	case ET_QUAD:
	  {
	    FlatVec<3> p0 = points(faces[fnr][0]);
	    FlatVec<3> p1 = points(faces[fnr][1]);
	    FlatVec<3> p2 = points(faces[fnr][3]);

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
    IntegrationRule ir;
    const ElementTransformation & eltrans;
    char * baseip;
    int incr;

  public:
    INLINE BaseMappedIntegrationRule (const IntegrationRule & air,
			       const ElementTransformation & aeltrans)
      : ir(air.Size(),&air[0]), eltrans(aeltrans) { ; }

    INLINE int Size() const { return ir.Size(); }
    INLINE const IntegrationRule & IR() const { return ir; }
    INLINE const ElementTransformation & GetTransformation () const { return eltrans; }

    INLINE const BaseMappedIntegrationPoint & operator[] (int i) const
    { return *static_cast<const BaseMappedIntegrationPoint*> ((void*)(baseip+i*incr)); }
  };

  template <int DIM_ELEMENT, int DIM_SPACE>
  class NGS_DLL_HEADER MappedIntegrationRule : public BaseMappedIntegrationRule
  {
    FlatArray< MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE> > mips;
  public:
    MappedIntegrationRule (const IntegrationRule & ir, 
			   const ElementTransformation & aeltrans, 
			   LocalHeap & lh);

    INLINE MappedIntegrationRule (const IntegrationRule & ir, 
			   const ElementTransformation & eltrans, 
			   int dummy,
			   LocalHeap & lh)
      : BaseMappedIntegrationRule (ir, eltrans), mips(ir.GetNIP(), lh)
    {
      baseip = (char*)(void*)(BaseMappedIntegrationPoint*)(&mips[0]);
      incr = (char*)(void*)(&mips[1]) - (char*)(void*)(&mips[0]);
    }

    INLINE MappedIntegrationRule (const IntegrationRule & air, 
			   const ElementTransformation & aeltrans, 
                           FlatArray< MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE> > amips)
      : BaseMappedIntegrationRule (air, aeltrans), mips(amips)
    {
      baseip = (char*)(void*)(BaseMappedIntegrationPoint*)(&mips[0]);
      incr = (char*)(void*)(&mips[1]) - (char*)(void*)(&mips[0]);
    }
    
    INLINE MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE> & operator[] (int i) const
    { 
      return mips[i]; 
    }

    INLINE MappedIntegrationRule Range(int first, int next)
    {
      return MappedIntegrationRule (ir.Range(first,next), eltrans, mips.Range(first,next));
    }
  };



  template <int DIM_ELEMENT, int DIM_SPACE>
  inline ostream & operator<< (ostream & ost, const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> & ir)
  {
    for (int i = 0; i < ir.Size(); i++)
      ost << ir[i] << endl;
    return ost;
  }




#define SpecificIntegrationPoint MappedIntegrationPoint 

}


#endif
