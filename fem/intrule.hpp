#ifndef FILE_INTRULE
#define FILE_INTRULE

/*********************************************************************/
/* File:   intrule.hpp                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


#include "elementtopology.hpp"  // for VorB

namespace ngcomp { class MeshAccess; }

namespace ngfem
{
  template <int DIM, typename T>
  class TIP;
  
  template <typename T>
  class TIP<0,T>
  {
  public:
    int8_t facetnr = -1;
    VorB vb = VOL;

    [[deprecated("Use TIP(facetnr, vb) instead")]]    
    TIP () = default;
    TIP (const TIP &) = default;
    TIP (TIP &&) = default;
    template <typename T2>
    TIP (const TIP<0,T2> & tip)
      : facetnr(tip.facetnr), vb(tip.vb) { }

    TIP (int8_t afacetnr, VorB avb)
      : facetnr(afacetnr), vb(avb) { } 
    
    TIP & operator= (const TIP &) = default;
    TIP & operator= (TIP &&) = default;
    
    explicit TIP (Vec<0,T> v, int8_t afacetnr=-1, VorB avb=VOL)
      : facetnr(afacetnr), vb(avb) { ; }    
    template <typename T1, typename T2>
    TIP (TIP<0,T1> ip1, TIP<0,T2> ip2) { ; } 
  };

  
  template <typename T>
  class TIP<1,T>
  {
  public:
    T x;
    int8_t facetnr = -1;
    VorB vb = VOL;

    [[deprecated("Use TIP(facetnr, vb) instead")]]    
    TIP () = default;
    TIP (int8_t afacetnr, VorB avb)
      : facetnr(afacetnr), vb(avb) { } 
    
    TIP (const TIP &) = default;
    TIP (TIP &&) = default;
    template <typename T2>
    TIP (const TIP<1,T2> & tip)
      : x(tip.x), facetnr(tip.facetnr), vb(tip.vb) { } 

    TIP & operator= (const TIP &) = default;
    TIP & operator= (TIP &&) = default;
    
    TIP (T _x, int8_t afacetnr, VorB avb)
      : x(_x), facetnr(afacetnr), vb(avb) { ; }
    
    explicit TIP (Vec<1,T> v, int8_t afacetnr = -1, VorB avb = VOL)
      : x(v(0)), facetnr(afacetnr), vb(avb) { ; }
    
    template <typename T1, typename T2>
    TIP (TIP<1,T1> ip1, TIP<1,T2> ip2)
      : x(ip1.x, ip2.x), facetnr(ip1.facetnr), vb(ip1.vb) { ; } 
  };

  template <typename T>
  class TIP<2,T>
  {
  public:
    T x, y;
    int8_t facetnr = -1;
    VorB vb = VOL;

    [[deprecated("Use TIP(facetnr, vb) instead")]]        
    TIP () = default;
    TIP (int8_t afacetnr, VorB avb)
      : facetnr(afacetnr), vb(avb) { } 
    
    TIP (const TIP &) = default;
    TIP (TIP &&) = default;
    template <typename T2>
    TIP (const TIP<2,T2> & tip)
      : x(tip.x), y(tip.y), facetnr(tip.facetnr), vb(tip.vb) { }

    TIP & operator= (const TIP &) = default;
    TIP & operator= (TIP &&) = default;
    
    TIP (T _x, T _y, int8_t afacetnr, VorB avb)
      : x(_x), y(_y), facetnr(afacetnr), vb(avb) { ; }
    explicit TIP (Vec<2,T> v, int8_t afacetnr = -1, VorB avb = VOL)
      : x(v(0)), y(v(1)), facetnr(afacetnr), vb(avb) { ; }
    
    template <typename T1, typename T2>
    TIP (TIP<2,T1> ip1, TIP<2,T2> ip2)
      : x(ip1.x, ip2.x), y(ip1.y, ip2.y), facetnr(ip1.facetnr), vb(ip1.vb) { ; } 
  };
  template <typename T>
  class TIP<3,T>
  {
  public:
    T x, y, z;
    int8_t facetnr = -1;
    VorB vb = VOL;

    [[deprecated("Use TIP(facetnr, vb) instead")]]            
    TIP () = default;
    TIP (int8_t afacetnr, VorB avb)
      : facetnr(afacetnr), vb(avb) { } 
    
    TIP (const TIP &) = default;
    TIP (TIP &&) = default;
    template <typename T2>
    TIP (const TIP<3,T2> & tip)
      : x(tip.x), y(tip.y), z(tip.z), facetnr(tip.facetnr), vb(tip.vb) { }
    
    TIP & operator= (const TIP &) = default;
    TIP & operator= (TIP &&) = default;
    
    TIP (T _x, T _y, T _z, int8_t afacetnr, VorB avb)
      : x(_x), y(_y), z(_z), facetnr(afacetnr), vb(avb) { ; }
    
    explicit TIP (Vec<3,T> v, int8_t afacetnr = -1, VorB avb = VOL)
      : x(v(0)), y(v(1)), z(v(2)), facetnr(afacetnr), vb(avb) { ; }
    
    template <typename T1, typename T2>    
    TIP (TIP<3,T1> ip1, TIP<3,T2> ip2)
      : x(ip1.x, ip2.x), y(ip1.y, ip2.y), z(ip1.z, ip2.z),
        facetnr(ip1.facetnr), vb(ip1.vb) { ; } 
  };

  template <typename T>
  inline ostream & operator<< (ostream & ost, TIP<0,T> tip)
  {
    return ost;
  }
  template <typename T>
  inline ostream & operator<< (ostream & ost, TIP<1,T> tip)
  {
    return ost << "x = " << tip.x;
  }
  template <typename T>
  inline ostream & operator<< (ostream & ost, TIP<2,T> tip)
  {
    return ost << "x = " << tip.x << ", y = " << tip.y;
  }
  template <typename T>
  inline ostream & operator<< (ostream & ost, TIP<3,T> tip)
  {
    return ost << "x = " << tip.x << ", y = " << tip.y << ", z = " << tip.z;
  }


  
  /// An integration point 
  class IntegrationPoint
  {
  private:
    /// number within integration Rule
    int nr; 
    /// coordinates (empty values for 1D and 2D)
    double pi[3];
    /// weight of integration point
    double weight;
    /// point is on facetnr, -1 for volume
    int8_t facetnr = -1;
    /// co-dimension of point (0..vol, 1..bnd, 2..bbnd, 3..bbbnd=vertex)
    VorB vb = VOL;
    ///
    // bool precomputed_geometry;
  public:
    ///
    INLINE IntegrationPoint & operator=(const IntegrationPoint & aip)
    {
      nr = aip.Nr();
      pi[0] = aip(0);
      pi[1] = aip(1);
      pi[2] = aip(2);
      weight = aip.Weight();
      // precomputed_geometry = aip.precomputed_geometry;
      facetnr = aip.facetnr;
      vb = aip.vb;
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
      // precomputed_geometry = 0;
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
      // precomputed_geometry = 0;
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
      // precomputed_geometry = 0;
    }


    template <int D>
    INLINE IntegrationPoint (const Vec<D> & ap, double aw = -1)
    {
      nr = -1;
      for (int j = 0; j < D; j++)
	pi[j] = ap(j);
      weight = aw;
      facetnr = -1;
      // precomputed_geometry = 0;
    }

    ///
    INLINE IntegrationPoint (const IntegrationPoint & aip)
    { *this = aip; }

    // void SetPrecomputedGeometry(bool value) { precomputed_geometry = value; }
    // bool GetPrecomputedGeometry() const { return precomputed_geometry; }
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
    // INLINE int IPNr () const { return -1; }

    INLINE auto Point() { return FlatVec<3> (&pi[0]); }
    INLINE auto Point() const { return FlatVec<3, const double> (&pi[0]); }

    void SetFacetNr (int afacetnr, VorB avb = BND)
    { facetnr = afacetnr; vb = avb; }
    INLINE auto FacetNr() const { return facetnr; }
    INLINE VorB VB() const { return vb; } 
    
    template <int DIM> 
    INLINE operator Vec<DIM, AutoDiff<DIM>> () const
    {
      Vec<DIM, AutoDiff<DIM> > adp;
      for (int i = 0; i < DIM; i++)
        adp[i] = AutoDiff<DIM> (pi[i], i);
      return adp;
    }

    INLINE operator TIP<0,double> () const { return TIP<0,double>(facetnr, vb); }
    INLINE operator TIP<1,double> () const { return TIP<1,double>(pi[0], facetnr, vb); }
    INLINE operator TIP<2,double> () const { return TIP<2,double>(pi[0], pi[1], facetnr, vb); }
    INLINE operator TIP<3,double> () const { return TIP<3,double>(pi[0], pi[1], pi[2], facetnr, vb); } 

    INLINE operator TIP<0,AutoDiff<0>> () const
    { return TIP<0,AutoDiff<0>>(facetnr, vb); } 
    INLINE operator TIP<1,AutoDiff<1>> () const
    { return TIP<1,AutoDiff<1>>(AutoDiff<1> (pi[0],0), facetnr, vb); }
    INLINE operator TIP<2,AutoDiff<2>> () const
    { return TIP<2,AutoDiff<2>>(AutoDiff<2> (pi[0],0), AutoDiff<2> (pi[1],1), facetnr, vb); }
    INLINE operator TIP<3,AutoDiff<3>> () const
    { return TIP<3,AutoDiff<3>>(AutoDiff<3> (pi[0],0), AutoDiff<3> (pi[1],1), AutoDiff<3> (pi[2],2), facetnr, vb); } 
    
    template <int D>
    INLINE ngfem::TIP<D,double> TIp() const;
    
    friend NGS_DLL_HEADER ostream & operator<< (ostream & ost, const IntegrationPoint & ip);
  };

  template <int D> INLINE ngfem::TIP<D,double> IntegrationPoint :: TIp() const
  { return ngfem::TIP<D,double> (*this); }
  

  struct MeshPoint
  {
    double x,y,z;
    ngcomp::MeshAccess* mesh;
    VorB vb;
    int nr;
  };


  class NGS_DLL_HEADER ElementTransformation;
  class NGS_DLL_HEADER BaseMappedIntegrationRule;
  /**
     Base class for MappedIntegrationPoint.
     A specific integration point is the mapped point, and stores
     point coordinates, the Jacobimatrix, and the determinant of the Jacobimatrix.
  */
  class BaseMappedIntegrationPoint
  {
  protected:
    /// IP on the reference element
    IntegrationPoint ip;
    /// computed by the transformation
    const ElementTransformation * eltrans;
    ///
    bool owns_trafo = false;
    ///
    bool is_complex;
    /// fabs(det)
    double measure;
    
    BaseMappedIntegrationPoint (const BaseMappedIntegrationPoint&) = default;
  public:
    ///
    INLINE BaseMappedIntegrationPoint () = default;
    ///
    INLINE BaseMappedIntegrationPoint (const IntegrationPoint & aip,
                                       const ElementTransformation & aeltrans)
      : ip(aip), eltrans(&aeltrans)  { ; }
    ///
    NGS_DLL_HEADER virtual ~BaseMappedIntegrationPoint ();
    /// 
    INLINE const IntegrationPoint & IP () const { return ip; }
    ///
    INLINE const ElementTransformation & GetTransformation () const { return *eltrans; }
    ///
    INLINE int GetIPNr() const { return ip.Nr(); }
    ///
    INLINE double GetMeasure() const { return measure; }
    void SetMeasure (double _measure) { measure = _measure; }    
    ///
    INLINE double GetWeight() const { return measure * ip.Weight(); }

    NGS_DLL_HEADER FlatVector<> GetPoint() const;
    FlatMatrix<> GetJacobian() const;

    // implemented in elementtransformation.hpp
    INLINE int DimElement() const; // { return eltrans->ElementDim(); }
    INLINE int DimSpace() const; // { return eltrans->SpaceDim(); } 
    
    FlatVector<Complex> GetPointComplex() const;
    FlatMatrix<Complex> GetJacobianComplex() const;
    // dimension of range
    // [[deprecated("Use DimSpace instead")]]
    // NGS_DLL_HEADER int Dim() const;  
    // VorB ElementVB() const; 
    bool IsComplex() const { return is_complex; }
    void SetOwnsTrafo (bool aowns_trafo = true) { owns_trafo = aowns_trafo; }
    virtual void IntegrationRuleFromPoint(std::function<void(const BaseMappedIntegrationRule&)> func) const { ; } 
  };

  template <typename SCAL = double>
  class alignas (sizeof(SCAL)) ScalMappedIntegrationPoint : public BaseMappedIntegrationPoint
  {
  protected:
    SCAL det;
  public:
    using BaseMappedIntegrationPoint :: BaseMappedIntegrationPoint;
    ScalMappedIntegrationPoint() { is_complex = false; }

    ScalMappedIntegrationPoint(const IntegrationPoint & aip,
                                       const ElementTransformation & aeltrans) 
      : BaseMappedIntegrationPoint(aip,aeltrans){ is_complex = false; }
    ///
    INLINE SCAL GetJacobiDet() const { return det; }
  };
  
  template<> INLINE ScalMappedIntegrationPoint<Complex> :: ScalMappedIntegrationPoint()
  { is_complex = true; }
  
  template<> INLINE ScalMappedIntegrationPoint<Complex> 
    :: ScalMappedIntegrationPoint(const IntegrationPoint & aip,
                                       const ElementTransformation & aeltrans)
      : BaseMappedIntegrationPoint(aip,aeltrans) { 
    is_complex = true; }

  
  template <int R, typename SCAL = double>
  class DimMappedIntegrationPoint : public ScalMappedIntegrationPoint<SCAL>
  {
  protected:
    ///
    Vec<R,SCAL> point;
    Vec<R,SCAL> normalvec;
    Vec<R,SCAL> tangentialvec;  // for independent integrator
    using ScalMappedIntegrationPoint<SCAL>::det;
  public:
    ///
    using ScalMappedIntegrationPoint<SCAL>::ScalMappedIntegrationPoint;
    /*
    INLINE DimMappedIntegrationPoint () = default;
    ///
    INLINE DimMappedIntegrationPoint (const IntegrationPoint & aip,
				 const ElementTransformation & aeltrans)
      : BaseMappedIntegrationPoint (aip, aeltrans)
    { ; }
    */
    ///
    INLINE const Vec<R,SCAL> & GetPoint () const { return point; }
    INLINE Vec<R,SCAL> & Point () { return point; }
    ///
    INLINE SCAL operator() (int i) const { return point(i); }

    ///
    INLINE const Vec<R,SCAL> & GetNV () const { return normalvec; }
    /// 
    INLINE void SetNV ( Vec<R,SCAL> vec) { normalvec = vec; }
    ///
    INLINE void SetTV ( Vec<R,SCAL> vec) { tangentialvec = vec; }
    ///
    INLINE const Vec<R,SCAL> GetTV () const { return tangentialvec; }
  };


  /// ip, dimension source, dimension range
  template <int DIMS, int DIMR, typename SCAL = double> 
  class MappedIntegrationPoint : public DimMappedIntegrationPoint<DIMR,SCAL>
  {
    static_assert(DIMS <= DIMR, "DIM-source > DIM-range !!");
  private:
    /// Jacobi matrix
    Mat<DIMR,DIMS,SCAL> dxdxi;
    /// for boundary points
    using DimMappedIntegrationPoint<DIMR,SCAL>::normalvec;
    using DimMappedIntegrationPoint<DIMR,SCAL>::tangentialvec;
    using DimMappedIntegrationPoint<DIMR,SCAL>::det;

  public:
    typedef SCAL TSCAL;
    ///
    NGS_DLL_HEADER MappedIntegrationPoint () = default;
    ///
    NGS_DLL_HEADER MappedIntegrationPoint (const IntegrationPoint & aip,
					   const ElementTransformation & aeltrans)
      : DimMappedIntegrationPoint<DIMR,SCAL> (aip, aeltrans)
    {
      this->eltrans->CalcPointJacobian(this->IP(), this->point, dxdxi);
      this->Compute();
    }

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

    INLINE void CheckDims() const
    {
      static_assert(DIMR >= DIMS, "DIMR >= DIMS not satisfied");
    }
    
    INLINE void Compute ()
    {
      if constexpr (DIMS == DIMR)
	{
	  det = Det (dxdxi);
	  normalvec = TSCAL(0.0);
	  tangentialvec = TSCAL(0.0);
	}
      else
	{
	  if (DIMR == 3)
	    {
	      if(DIMS == 2)
		{
		  normalvec = Cross (Vec<3,SCAL> (dxdxi.Col(0)),
				     Vec<3,SCAL> (dxdxi.Col(1)));
		  det = L2Norm (normalvec);
		  normalvec /= det;
                  tangentialvec = TSCAL(0.0);
		}
	      else if (DIMS == 1)
		{
		  // CHECK!
		  normalvec = TSCAL(0.0);
		  tangentialvec = Vec<3,SCAL>(dxdxi.Col(0));
		  det = L2Norm(tangentialvec);
		  tangentialvec /= det;
		}
              else
                {
                  det = 1;
                }
	    }
	  else if (DIMR == 2)
	    {
              if (DIMS == 1)
                {
                  det = sqrt ( ngstd::sqr (dxdxi(0,0)) + ngstd::sqr (dxdxi(1,0)));
                  
                  normalvec(0) = -dxdxi(1,0) / det;
                  normalvec(1) = dxdxi(0,0) / det;
                  //tangentialvec = TSCAL(0.0);
                  tangentialvec(0) = -normalvec(1);
                  tangentialvec(1) = normalvec(0);
                }
              else
                {
                  det = 1;
                }
            }
	  else
	    {
	      det = 1.0;
	      normalvec = 1.0;
              tangentialvec = TSCAL(0.0);
	    }
	}
      this->measure = fabs (det);
    }
  
    ///
    INLINE const Mat<DIMR,DIMS,SCAL> & GetJacobian() const { return dxdxi; }
    INLINE Mat<DIMR,DIMS,SCAL> & Jacobian() { return dxdxi; }
    ///
    // const Mat<DIMS,DIMR,SCAL> & GetJacobianInverse() const { return dxidx; }
    INLINE const Mat<DIMS,DIMR,SCAL> GetJacobianInverse() const 
    { 
      if (DIMS == DIMR)
        // return Inv (dxdxi);
        return 1.0/det * Trans (Cof (dxdxi));
      else
        {
	  // Mat<DIMS,DIMS,SCAL> ata, iata;
	  // ata = Trans (dxdxi) * dxdxi;
	  // iata = Inv (ata);
          // auto iata = Inv (Trans (dxdxi) * dxdxi);
	  // return iata * Trans (dxdxi);
          return Inv (Trans (dxdxi) * dxdxi) * Trans (dxdxi);
        }
    }

    INLINE const Mat<DIMR,DIMS,SCAL> GetJacobianCofactor() const 
    { 
      if (DIMS == DIMR)
        return Cof (dxdxi);
      else
        return det * Trans(GetJacobianInverse());
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
    // INLINE VorB VB() const { return VorB(DIMR-DIMS); }
    //INLINE int IsBoundary () const { return DIMS != DIMR; }

    ///
    void CalcHesse (Mat<1> & ddx1, Mat<1> & ddx2) const;
    void CalcHesse (Mat<2> & ddx1, Mat<2> & ddx2) const;
    void CalcHesse (Mat<1> & ddx1, Mat<1> & ddx2, Mat<1> & ddx3) const;
    void CalcHesse (Mat<2> & ddx1, Mat<2> & ddx2, Mat<2> & ddx3) const; 
    void CalcHesse (Mat<3> & ddx1, Mat<3> & ddx2, Mat<3> & ddx3) const;

    void CalcHesse (Vec<DIMR,Mat<DIMS,DIMS>> & ddx1) const;
    Vec<DIMR,Mat<DIMS,DIMS>> CalcHesse() const
    {
      if constexpr(std::is_same<SCAL,double>::value)
        {
          Vec<DIMR,Mat<DIMS,DIMS>> hesse;
          CalcHesse (hesse);
          return hesse;
        }
      throw Exception ("CalcHesse not available for complex mips");
    }



    void IntegrationRuleFromPoint(std::function<void(const BaseMappedIntegrationRule&)> func) const override;
  };



  inline ostream & operator<< (ostream & ost, const BaseMappedIntegrationPoint & mip)
  {
    ost << mip.GetPoint(); //  << ", dxdxi = " << mip.GetJacobian();
    return ost;
  }

  template <int DIMS, int DIMR, typename SCAL> 
  inline ostream & operator<< (ostream & ost, const MappedIntegrationPoint<DIMS,DIMR,SCAL> & mip)
  {
    ost << mip.GetPoint() << ", dxdxi = " << mip.GetJacobian();
    return ost;
  }





  /**
     An integration rule.
     Contains array of integration points
  */

  class IntegrationRule : public Array<IntegrationPoint> 
  {
    int dimension = -1;
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

    // copies pointers
    INLINE NGS_DLL_HEADER IntegrationRule (const IntegrationRule & ir2)
      : Array<IntegrationPoint> (ir2.Size(), ir2.data), dimension(ir2.dimension)
    { ; }

    INLINE NGS_DLL_HEADER IntegrationRule (IntegrationRule && ir2) = default;

    INLINE NGS_DLL_HEADER IntegrationRule (size_t asize, LocalHeap & lh)
      : Array<IntegrationPoint> (asize, lh)
    { ; }

    INLINE NGS_DLL_HEADER IntegrationRule (size_t asize, double (*pts)[3], double * weights);
    
    // make it polymorphic
    HD virtual ~IntegrationRule() { ; }

    IntegrationRule & operator= (IntegrationRule && ir2) = default;
    IntegrationRule & operator= (const IntegrationRule & ir2) = delete;
    
    ///
    INLINE void AddIntegrationPoint (const IntegrationPoint & ip)
    { 
      Append (ip);
    }

    IntegrationRule Copy() const;
    
    /// number of integration points
    INLINE size_t GetNIP() const { return Size(); }

    INLINE IntegrationRule Range (size_t first, size_t next) const
    {
      return IntegrationRule (next-first, &(*this)[first]);
    }
    int Dim() const { return dimension; }
    void SetDim (int dim) { dimension = dim; }

    struct IntegrationRuleSplitArray
    {
        const IntegrationRule & ir;
        static constexpr size_t BS=128;
        size_t size;

        IntegrationRuleSplitArray(const IntegrationRule & ir_)
        : ir(ir_)
        {
            size = (ir.Size()+BS-1)/BS;
        }

        IntegrationRule operator[](size_t i) const
        {
            size_t first = i*BS;
            size_t next = min((i+1)*BS, ir.Size());
            return ir.Range(first, next);
        }

        auto begin() const { return AOWrapperIterator(*this, 0); }
        auto end() const { return AOWrapperIterator(*this, size); }
    };

    IntegrationRuleSplitArray Split() { return *this; }
  };

  NGS_DLL_HEADER ostream & operator<< (ostream & ost, const IntegrationRule & ir);

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
  
  NGS_DLL_HEADER ostream & operator<< (ostream & ost, const DGIntegrationRule & ir);

  template <int D>
  class IntegrationRuleTP : public IntegrationRule
  {
    const IntegrationRule *irx, *iry, *irz;

    ArrayMem<Mat<D,D>, 100> dxdxi_duffy;
    Mat<D,D> dxdxi_permute;

  public:
    NGS_DLL_HEADER IntegrationRuleTP (const ElementTransformation & eltrans,
                                      IVec<D> order, bool compute_duffy = true, bool compute_points = true); 

    // tensor product rule for a facet
    NGS_DLL_HEADER IntegrationRuleTP (ELEMENT_TYPE eltype, FlatArray<int> sort, 
		       NODE_TYPE nt, int nodenr, int order, LocalHeap & lh);

    const IntegrationRule & GetIRX() const { return *irx; }
    const IntegrationRule & GetIRY() const { return *iry; }
    const IntegrationRule & GetIRZ() const { return *irz; }

    int GetNIP() const
    {
      switch (D)
        {
        case 1: return irx->GetNIP();
        case 2: return irx->GetNIP()*iry->GetNIP();
        case 3: return irx->GetNIP()*iry->GetNIP()*irz->GetNIP();
        }
      return 0;
    }
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

    template <ELEMENT_TYPE ET>
    Mat<D,D> GetDuffyJacobian(int ix, int iy, int iz) const 
    { 
      if (dxdxi_duffy.Size())
        return dxdxi_duffy[(ix*iry->GetNIP()+iy)*irz->GetNIP()+iz]; 
      else
        {
          if (ET == ET_TET)
            {
              double x = (*irx)[ix](0);
              double invx = 1.0 / (1-x);
                  
              double y = (*iry)[iy](0);
              double invxy = 1.0 / ( (1-x) * (1-y) );
                      
              double z = (*irz)[iz](0);

              Mat<3> trans3; 
              
              // hex -> tet transform
              trans3(0,0) = 1;
              trans3(0,1) = 0;
              trans3(0,2) = 0;
              
              trans3(1,0) = y*invx;
              trans3(1,1) = 1*invx;
              trans3(1,2) = 0;
              
              trans3(2,0) = z*invxy;
              trans3(2,1) = z*invxy;
              trans3(2,2) = 1*invxy;
              
              return trans3 * dxdxi_permute;
            }
          else
            {
              cout << "need to compute duffy for ET " << ET << endl;
              return Mat<D,D> (0.0);
            }
        }
    }

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
     Gauss-Lobatto Rule with n points on interval [0,1].
     Contains 0 as points
     exact for polynomials up to order 2n-2
   */
  NGS_DLL_HEADER extern 
  void ComputeGaussRadauRule (int n, Array<double>& xi, Array<double>& wi);


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

  INLINE IntegrationRule :: IntegrationRule (ELEMENT_TYPE eltype, int order)
  { 
    const IntegrationRule & ir = SelectIntegrationRule (eltype, order);
    size = ir.Size();
    data = &ir[0];
    mem_to_delete = NULL;
    dimension = ElementTopology::GetSpaceDim(eltype);
  }



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
    bool swapped = false; // new orientation with 2 tet-classes
    VorB vb = BND;        // facet codimension
  public:
    Facet2ElementTrafo(ELEMENT_TYPE aeltype, VorB _vb = BND) 
      : eltype(aeltype),
	points(99,(Vec<3>*)ElementTopology::GetVertices (aeltype)),
        vb(_vb)
    {
      // points = ElementTopology::GetVertices (eltype);
      edges = ElementTopology::GetEdges (eltype);
      faces = ElementTopology::GetFaces (eltype);
    }
  
    // Facet2ElementTrafo(ELEMENT_TYPE aeltype, const FlatArray<int> & vnums)

    template <typename T>
    Facet2ElementTrafo(ELEMENT_TYPE aeltype, const BaseArrayObject<T> & vnums) 
      : eltype(aeltype),
	points(99,(Vec<3>*)ElementTopology::GetVertices (aeltype))
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
          swapped = vnums[2] > vnums[3];
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

    ELEMENT_TYPE FacetType (int fnr) const
    {
      if (vb == VOL)
        return eltype;
      if (vb == BND)
        return ElementTopology::GetFacetType(eltype, fnr);
      else
        {
          if (Dim(eltype)-int(vb) == 1)
            return ET_SEGM;
          else
            return ET_POINT;
        }
    }

    size_t GetNFacets () const
    {
      if (vb == VOL)
        return 1;
      else if (vb == BND)
        return ElementTopology::GetNFacets(eltype);
      else
        {
          if (Dim(eltype)-int(vb) == 1)
            return ElementTopology::GetNEdges(eltype);
          else
            return ElementTopology::GetNVertices(eltype); // points
        }
    }
    
    void operator()(int fnr, const IntegrationPoint &ipfac, IntegrationPoint & ipvol) const 
    {
      if (vb == VOL)
        {
          ipvol = ipfac;
          return;
        }
      
      switch (FacetType(fnr))
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
      ipvol.SetWeight(ipfac.Weight());
      // ipvol.FacetNr() = fnr;
      ipvol.SetFacetNr(fnr, vb);
    }

    FlatMatrix<> GetJacobian(int fnr, LocalHeap & lh) const
    {
      // ELEMENT_TYPE facettype = ElementTopology::GetFacetType(eltype, fnr);
      // switch (facettype)
      switch (FacetType(fnr))
	{
	case ET_SEGM:
	  {
            FlatMatrix<> mat(2,1,lh);
	    FlatVec<3> p1 = points (edges[fnr][0]);
	    FlatVec<3> p2 = points (edges[fnr][1]);
            mat.Col(0) = (p1 - p2).Range(0,2);
            return mat;
	    break;
	  }
	case ET_TRIG:
	  {
            FlatMatrix<> mat(3,2,lh);
	    FlatVec<3> p0 = points(faces[fnr][0]);
	    FlatVec<3> p1 = points(faces[fnr][1]);
	    FlatVec<3> p2 = points(faces[fnr][2]);
            mat.Col(0) = p0 - p2;
            mat.Col(1) = p1 - p2;
            return mat;
	    break;
	  }
	case ET_QUAD:
	  {
            FlatMatrix<> mat(3,2,lh);
	    FlatVec<3> p0 = points(faces[fnr][0]);
	    FlatVec<3> p1 = points(faces[fnr][1]);
	    FlatVec<3> p2 = points(faces[fnr][3]);
            mat.Col(0) = p1 - p0;
            mat.Col(1) = p2 - p0;
            return mat;
	    break;
	  }
	default:
	  throw Exception ("undefined facet type in Facet2ElementTrafo::GetJacobian(..)\n");
	}
    }

    
    const IntegrationPoint operator()(int fnr, const IntegrationPoint &ip1d) const 
    {
      IntegrationPoint elpoint;
      operator()(fnr, ip1d, elpoint);
      return elpoint;
    }

    IntegrationRule & operator() (int fnr, const IntegrationRule & irfacet, LocalHeap & lh)
    {
      if (vb == VOL) return const_cast<IntegrationRule&> (irfacet);
      
      IntegrationRule & irvol = *new (lh) IntegrationRule (irfacet.GetNIP(), lh);

      // switch (ElementTopology::GetFacetType(eltype, fnr))
      switch (FacetType(fnr))        
	{
	case ET_POINT:
	  {
	    irvol[0] = Vec<3> (points (fnr) );
	    break;
	  }
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

      for (int i = 0; i < irfacet.GetNIP(); i++)
        {
          irvol[i].SetFacetNr(fnr, vb);
          irvol[i].SetWeight(irfacet[i].Weight());
        }

      return irvol;
    }


    NGS_DLL_HEADER const class SIMD_IntegrationRule & operator() (int fnr, const class SIMD_IntegrationRule & irfacet, LocalHeap & lh);

  };







  // transformation of (d-1) dimensional integration points on facets to 
  // d-dimensional point in volumes
  class Facet2SurfaceElementTrafo
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
    Facet2SurfaceElementTrafo(ELEMENT_TYPE aeltype) :
      eltype(aeltype),
      points(99,(Vec<3>*)ElementTopology::GetVertices (aeltype))
    {
      // points = ElementTopology::GetVertices (eltype);
      edges = ElementTopology::GetEdges (eltype);
      faces = ElementTopology::GetFaces (eltype);
    }
    
    Facet2SurfaceElementTrafo(ELEMENT_TYPE aeltype, FlatArray<int> & vnums)
      : eltype(aeltype),
        points(99,(Vec<3>*)ElementTopology::GetVertices (aeltype))
    {
      // points = ElementTopology::GetVertices (eltype);
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
	case ET_POINT:
	  {
	    ipvol = Vec<3> (points (fnr) );
	    break;
	  }
          
	case ET_SEGM:
	  {
            /*
	    const POINT3D & p1 = points[edges[fnr][0]];
	    const POINT3D & p2 = points[edges[fnr][1]];
            for (int j = 0; j < 3; j++)
              ipvol(j) = p2[j] + (ipfac(0))*(p1[j]-p2[j]);
            */
	    FlatVec<3> p1 = points (edges[fnr][0]);
	    FlatVec<3> p2 = points (edges[fnr][1]);
	    ipvol = Vec<3> (p2 + ipfac(0) * (p1-p2));
	    break;
	  }
	case ET_TRIG:
	  {
            /*
	    const POINT3D & p0 = points[faces[fnr][0]];
	    const POINT3D & p1 = points[faces[fnr][1]];
	    const POINT3D & p2 = points[faces[fnr][2]];
	    ipvol(0) = p2[0] + ipfac(0)*(p0[0]-p2[0]) + ipfac(1)*(p1[0]-p2[0]);
	    ipvol(1) = p2[1] + ipfac(0)*(p0[1]-p2[1]) + ipfac(1)*(p1[1]-p2[1]);
	    ipvol(2) = p2[2] + ipfac(0)*(p0[2]-p2[2]) + ipfac(1)*(p1[2]-p2[2]);
            */
	    FlatVec<3> p0 = points(faces[fnr][0]);
	    FlatVec<3> p1 = points(faces[fnr][1]);
	    FlatVec<3> p2 = points(faces[fnr][2]);
	    ipvol = Vec<3> (p2 + ipfac(0) * (p0-p2) + ipfac(1)*(p1-p2));
	    break;
	  }
	case ET_QUAD:
	  {
            /*
	    const POINT3D & p0 = points[faces[fnr][0]];
	    const POINT3D & p1 = points[faces[fnr][1]];
	    const POINT3D & p2 = points[faces[fnr][3]];
	    ipvol(0) = p0[0] + ipfac(0)*(p1[0]-p0[0]) + ipfac(1)*(p2[0]-p0[0]);
	    ipvol(1) = p0[1] + ipfac(0)*(p1[1]-p0[1]) + ipfac(1)*(p2[1]-p0[1]);
	    ipvol(2) = p0[2] + ipfac(0)*(p1[2]-p0[2]) + ipfac(1)*(p2[2]-p0[2]);
            */
	    FlatVec<3> p0 = points(faces[fnr][0]);
	    FlatVec<3> p1 = points(faces[fnr][1]);
	    FlatVec<3> p2 = points(faces[fnr][3]);
	    ipvol = Vec<3> (p0 + ipfac(0) * (p1-p0) + ipfac(1)*(p2-p0));
	    break;
	  }
	default:
	  throw Exception ("undefined facet type in Facet2ElementTrafo()\n");

	} 
      /*      cerr << "*** mapping integrationpoint for element " << eltype << " and facel " << fnr << " of type " << facettype << endl;
	      cerr << "  * ipfac = " << ipfac;
	      cerr << "  * ipvol = " << ipvol;*/      
    }
      
    IntegrationRule & operator() (const IntegrationRule & irfacet, LocalHeap & lh)
    {
      IntegrationRule & irvol = *new (lh) IntegrationRule (irfacet.GetNIP(), lh);
      int fnr = 0;
      switch (eltype)
        {
        case ET_POINT:
          {
            irvol[0] = Vec<3> (points (fnr) );
            break;
          }
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
	  throw Exception ("undefined facet type in Facet2SurfaceElementTrafo()\n");
        } 
      
      for (int i = 0; i < irfacet.GetNIP(); i++)
        {
          // irvol[i].SetFacetNr(fnr);
          irvol[i].SetWeight(irfacet[i].Weight());
        }
      
      return irvol;
    }

    const IntegrationPoint operator()(const IntegrationPoint &ip1d) const 
    {
      IntegrationPoint elpoint;
      operator()(ip1d, elpoint);
      return elpoint;
    }

    NGS_DLL_HEADER class SIMD_IntegrationRule & operator() (const class SIMD_IntegrationRule & irfacet, LocalHeap & lh);

    IntegrationRule & Inverse (const IntegrationRule & ir, LocalHeap & lh);    
    class SIMD_IntegrationRule & Inverse (const class SIMD_IntegrationRule & ir, LocalHeap & lh);
  };













  class NGS_DLL_HEADER BaseMappedIntegrationRule 
  { 
  protected:
    IntegrationRule ir;
    const ElementTransformation & eltrans;
    char * baseip;
    size_t incr;
    // mir on other element as needed for evaluating DG jump terms
    const BaseMappedIntegrationRule * other_mir = nullptr;
    
  public:    
    INLINE BaseMappedIntegrationRule (const IntegrationRule & air,
                                      const ElementTransformation & aeltrans)
      : ir(air.Size(),&air[0]), eltrans(aeltrans) { ; }
    INLINE ~BaseMappedIntegrationRule ()
    {
      ir.NothingToDelete();
    }
    INLINE size_t Size() const { return ir.Size(); }
    INLINE const IntegrationRule & IR() const { return ir; }
    INLINE const ElementTransformation & GetTransformation () const { return eltrans; }

    INLINE BaseMappedIntegrationPoint & operator[] (size_t i) const
    { return *static_cast<BaseMappedIntegrationPoint*> ((void*)(baseip+i*incr)); }

    virtual BaseMappedIntegrationRule & Range(size_t first, size_t next, LocalHeap & lh) const = 0;

    auto begin () const { return AOWrapperIterator<BaseMappedIntegrationRule> (*this, 0); }
    auto end () const { return AOWrapperIterator<BaseMappedIntegrationRule> (*this, Size()); }

    int DimElement() const;
    int DimSpace() const;
    
    virtual SliceMatrix<> GetPoints() const = 0;
    virtual SliceMatrix<> GetNormals() const = 0;    
    virtual SliceMatrix<Complex> GetPointsComplex() const
    { throw Exception("don't have complex ir"); }
    virtual void ComputeNormalsAndMeasure (ELEMENT_TYPE et) { throw Exception ("ComputeNormalsAndMeasure(ET) not overloaded"); }
    virtual void ComputeNormalsAndMeasure (ELEMENT_TYPE et, int facetnr) = 0;
    virtual bool IsComplex() const = 0;

    // for DG jump terms
    void SetOtherMIR (const BaseMappedIntegrationRule * other) { other_mir = other; }
    auto GetOtherMIR () const { return other_mir; }
  };

  template <int DIM_ELEMENT, int DIM_SPACE, typename SCAL = double>
  class NGS_DLL_HEADER MappedIntegrationRule : public BaseMappedIntegrationRule
  {
    static_assert(DIM_ELEMENT <= DIM_SPACE, "DIM-source > DIM-range !!");    
    FlatArray< MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE, SCAL> > mips;
  public:
    typedef SCAL TSCAL;    
    MappedIntegrationRule (const IntegrationRule & ir, 
			   const ElementTransformation & aeltrans, 
			   Allocator & lh);

    INLINE MappedIntegrationRule (const IntegrationRule & ir, 
                                  const ElementTransformation & eltrans, 
                                  int dummy,
                                  Allocator & lh)
      : BaseMappedIntegrationRule (ir, eltrans), mips(ir.Size(), lh)
    {
      baseip = (char*)(void*)(BaseMappedIntegrationPoint*)(&mips[0]);
      incr = sizeof (MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE,SCAL>);
    }

    INLINE MappedIntegrationRule (const IntegrationRule & air, 
                                  const ElementTransformation & aeltrans, 
                                  FlatArray< MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE,SCAL> > amips)
      : BaseMappedIntegrationRule (air, aeltrans), mips(amips)
    {
      baseip = (char*)(void*)(BaseMappedIntegrationPoint*)(&mips[0]);
      if (mips.Size() > 1)
        incr = (char*)(void*)(&mips[1]) - (char*)(void*)(&mips[0]);
      else
        incr = 0;
    }
    
    INLINE MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE> & operator[] (size_t i) const
    { 
      // return mips[i];
      return static_cast<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE> &> (BaseMappedIntegrationRule::operator[] (i));
    }

    INLINE MappedIntegrationRule Range(size_t first, size_t next) const
    {
      return MappedIntegrationRule (ir.Range(first,next), eltrans, mips.Range(first,next));
    }

    virtual BaseMappedIntegrationRule & Range(size_t first, size_t next, LocalHeap & lh) const
    {
      return *new (lh) MappedIntegrationRule (ir.Range(first,next), eltrans, mips.Range(first,next));
    }

    auto begin() { return mips.begin(); }
    auto end() { return mips.end(); }
    
    virtual SliceMatrix<> GetPoints() const
    {
      return SliceMatrix<> (mips.Size(), DIM_SPACE*sizeof(SCAL)/sizeof(double),
                            sizeof(MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE, SCAL>) / sizeof(double),
                            const_cast<double*> (&mips[0].GetPoint()(0)));
    }

    virtual SliceMatrix<> GetNormals() const
    {
      return SliceMatrix<> (mips.Size(), DIM_SPACE*sizeof(SCAL)/sizeof(double),
                            sizeof(MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE, SCAL>) / sizeof(double),
                            const_cast<double*> (&mips[0].GetNV()(0)));
    }

    virtual void ComputeNormalsAndMeasure (ELEMENT_TYPE et);
    virtual void ComputeNormalsAndMeasure (ELEMENT_TYPE et, int facetnr);
    virtual bool IsComplex() const { return false; } 
  };

  template <int DIM_ELEMENT, int DIM_SPACE>
  class NGS_DLL_HEADER MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE,Complex> : public BaseMappedIntegrationRule
  {
    using SCAL = Complex;
    FlatArray< MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE, Complex> > mips;
  public:
    MappedIntegrationRule (const IntegrationRule & ir, 
			   const ElementTransformation & aeltrans, 
			   Allocator & lh);

    INLINE MappedIntegrationRule (const IntegrationRule & ir, 
                                  const ElementTransformation & eltrans, 
                                  int dummy,
                                  Allocator & lh)
      : BaseMappedIntegrationRule (ir, eltrans), mips(ir.Size(), lh)
    {
      baseip = (char*)(void*)(BaseMappedIntegrationPoint*)(&mips[0]);
      incr = sizeof (MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE,SCAL>);
    }

    INLINE MappedIntegrationRule (const IntegrationRule & air, 
                                  const ElementTransformation & aeltrans, 
                                  FlatArray< MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE,SCAL> > amips)
      : BaseMappedIntegrationRule (air, aeltrans), mips(amips)
    {
      baseip = (char*)(void*)(BaseMappedIntegrationPoint*)(&mips[0]);
      if (mips.Size() > 1)
        incr = (char*)(void*)(&mips[1]) - (char*)(void*)(&mips[0]);
      else
        incr = 0;
    }
    
    INLINE MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE,SCAL> & operator[] (int i) const
    { 
      return mips[i]; 
    }

    INLINE MappedIntegrationRule Range(size_t first, size_t next) const
    {
      return MappedIntegrationRule (ir.Range(first,next), eltrans, mips.Range(first,next));
    }

    virtual BaseMappedIntegrationRule & Range(size_t first, size_t next, LocalHeap & lh) const
    {
      return *new (lh) MappedIntegrationRule (ir.Range(first,next), eltrans, mips.Range(first,next));
    }
    
    virtual SliceMatrix<> GetPoints() const
    {
      throw Exception("don't have real points for complex ir");
    }

    virtual SliceMatrix<Complex> GetPointsComplex() const
    {
      return SliceMatrix<Complex> (mips.Size(), DIM_SPACE,
                                   //&mips[1].GetPoint()(0) - &mips[0].GetPoint()(0),
                                   sizeof(MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE, SCAL>) / sizeof(Complex),
                                   const_cast<Complex*> (&mips[0].GetPointComplex()(0)));
    }

    virtual SliceMatrix<> GetNormals() const
    {
      throw Exception("never tested");
    }
    

    virtual void ComputeNormalsAndMeasure (ELEMENT_TYPE et, int facetnr);
    virtual bool IsComplex() const { return true; }     
  };

  
  ostream & operator<< (ostream & ost, const BaseMappedIntegrationRule & mir);

  
  template <int DIM_ELEMENT, int DIM_SPACE, typename SCAL>
  inline ostream & operator<< (ostream & ost, const MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE,SCAL> & mir)
  {
    for (auto & mip : mir)
      ost << mip << endl;
    return ost;
  }

  template <int DIMS, int DIMR, typename SCAL>
  void MappedIntegrationPoint<DIMS,DIMR,SCAL>::
  IntegrationRuleFromPoint(std::function<void(const BaseMappedIntegrationRule&)> func) const
  {
    if constexpr (std::is_same_v<SCAL,double> || std::is_same_v<SCAL,Complex>)
      {
        FlatArray<MappedIntegrationPoint<DIMS,DIMR,SCAL>> ia(1, const_cast<MappedIntegrationPoint*>(this));
        IntegrationRule ir(1, const_cast<IntegrationPoint*>(&this->IP()));
        MappedIntegrationRule<DIMS,DIMR,SCAL> mir(ir, this->GetTransformation(), ia);
        func (mir);
      }
  }
  



  
  // deprecated, don't use SpecificIntegrationPoint anymore
  template <int DIMS, int DIMR, typename SCAL>   
  using SpecificIntegrationPoint = MappedIntegrationPoint<DIMS,DIMR,SCAL>;
}

namespace ngcore
{
  using ngbla::Vec;
  using ngbla::Mat;
  using ngbla::FlatVector;
  
  template<>
  class alignas(sizeof(SIMD<double>)) SIMD<ngfem::IntegrationPoint>
  {
    SIMD<double> x[3], weight;
    int facetnr = -1;
    ngfem::VorB vb = ngfem::VOL;
  public:
    static constexpr size_t Size() { return SIMD<double>::Size(); }

    SIMD() = default;
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;
    
    template <typename Function>
    SIMD (const Function & func)
    {
      ngfem::IntegrationPoint ip[Size()];
      for (int i = 0; i < Size(); i++) ip[i] = func(i);
      for (int j = 0; j < 3; j++)
        x[j] = [&ip,j] (int i) { return ip[i](j); };
      weight = [&ip] (int i) { return ip[i].Weight(); };
    }

    template <int DIM>
    operator Vec<DIM,SIMD<double>>() const
    {
      Vec<DIM,SIMD<double>> hp;
      for (int i = 0; i < DIM; i++)
        hp(i) = x[i];
      return hp;
    }
    
    const SIMD<double> & operator() (size_t i) const { return x[i]; }
    SIMD<double> & operator() (size_t i) { return x[i]; }
    SIMD<double> Weight() const { return weight; }
    SIMD<double> & Weight() { return weight; }
    ngfem::IntegrationPoint operator[] (size_t i) const
    { return ngfem::IntegrationPoint(x[0][i], x[1][i], x[2][i], weight[i]); }

    template<int I>
    auto Get()
    {
      static_assert(I>=0 && I<SIMD<double>::Size(), "Index out of range");
      return (*this)[I];
    }

    int FacetNr() const { return facetnr; }
    void SetFacetNr (int afacetnr, ngfem::VorB avb = ngfem::BND)
    { facetnr = afacetnr; vb = avb; }      
    INLINE ngfem::VorB VB() const { return vb; } 

    template <int DIM> 
    INLINE operator Vec<DIM, ngstd::AutoDiff<DIM,SIMD<double>>> () const
    {
      Vec<DIM, ngstd::AutoDiff<DIM,SIMD<double>> > adp;
      for (int i = 0; i < DIM; i++)
        adp[i] = ngstd::AutoDiff<DIM,SIMD<double>> (x[i], i);
      return adp;
    }

    template <int D>
    INLINE ngfem::TIP<D,SIMD<double>> TIp() const;
    INLINE operator ngfem::TIP<0,ngcore::SIMD<double>> () const { return ngfem::TIP<0,ngcore::SIMD<double>>(facetnr, vb); }
    INLINE operator ngfem::TIP<1,ngcore::SIMD<double>> () const { return ngfem::TIP<1,ngcore::SIMD<double>>(x[0], facetnr, vb); }
    INLINE operator ngfem::TIP<2,ngcore::SIMD<double>> () const { return ngfem::TIP<2,ngcore::SIMD<double>>(x[0], x[1], facetnr, vb); }
    INLINE operator ngfem::TIP<3,ngcore::SIMD<double>> () const { return ngfem::TIP<3,ngcore::SIMD<double>>(x[0], x[1], x[2], facetnr, vb); }

    /*
    template <int DIM> 
    INLINE operator Vec<DIM, AutoDiff<DIM,SIMD<double>>> () const
    {
      Vec<DIM, AutoDiff<DIM,SIMD<double>> > adp;
      for (int i = 0; i < DIM; i++)
        adp[i] = AutoDiff<DIM,SIMD<double>> (x[i], i);
      return adp;
    }
    */
  };

  /*
  template <> INLINE ngfem::TIP<0,SIMD<double>> SIMD<ngfem::IntegrationPoint> :: TIp<0>() const
  { return ngfem::TIP<0,ngstd::SIMD<double>>(); }
  template <> INLINE ngfem::TIP<1,SIMD<double>> SIMD<ngfem::IntegrationPoint> :: TIp<1>() const
  { return ngfem::TIP<1,ngstd::SIMD<double>>(x[0]); }
  template <> INLINE ngfem::TIP<2,SIMD<double>> SIMD<ngfem::IntegrationPoint> :: TIp<2>() const
  { return ngfem::TIP<2,ngstd::SIMD<double>>(x[0], x[1]); }
  template <> INLINE ngfem::TIP<3,SIMD<double>> SIMD<ngfem::IntegrationPoint> :: TIp<3>() const
  { return ngfem::TIP<3,ngstd::SIMD<double>>(x[0], x[1], x[2]); }
  */
  template <int D> INLINE ngfem::TIP<D,SIMD<double>> SIMD<ngfem::IntegrationPoint> :: TIp() const
  {
    // return ngfem::TIP<D,ngstd::SIMD<double>> (*this);
    ngfem::TIP<D,ngcore::SIMD<double>> tip = *this;
    return tip;
  }
  
  
  template <>
  class SIMD<ngfem::BaseMappedIntegrationPoint>
  {
  protected:
    SIMD<ngfem::IntegrationPoint> ip;
    const ngfem::ElementTransformation * eltrans;
    SIMD<double> measure;
    SIMD<double> det;
    // bool owns_trafo = false;
  public:
    SIMD() = default;
    SIMD (const SIMD<ngfem::IntegrationPoint> & aip,
          const ngfem::ElementTransformation & aeltrans)
      : ip(aip), eltrans(&aeltrans) { ; }
    // ~SIMD();
    const SIMD<ngfem::IntegrationPoint> & IP () const { return ip; }
    const ngfem::ElementTransformation & GetTransformation () const { return *eltrans; }
    // int GetIPNr() const { return ip.Nr(); }

    int DimElement() const;
    int DimSpace() const;

    void SetMeasure (SIMD<double> _measure) { measure = _measure; }
    SIMD<double> GetMeasure() const { return measure; }
    SIMD<double> GetWeight() const { return measure * ip.Weight(); }
    SIMD<double> GetJacobiDet() const { return det; }
    NGS_DLL_HEADER FlatVector<SIMD<double>> GetPoint() const;
    // virtual void Print (ostream & ost) const = 0;
  };

  template <int R>
  class SIMD<ngfem::DimMappedIntegrationPoint<R>> : public SIMD<ngfem::BaseMappedIntegrationPoint>
  {
  protected:
    ngbla::Vec<R,SIMD<double>> point, normalvec, tangentialvec;
  public:
    SIMD() = default;
    SIMD (const SIMD<ngfem::IntegrationPoint> & aip,
          const ngfem::ElementTransformation & eltrans)
      : SIMD<ngfem::BaseMappedIntegrationPoint> (aip, eltrans) { ; }
    const ngbla::Vec<R,SIMD<double>> & GetPoint() const { return point; }
    ngbla::Vec<R,SIMD<double>> & Point() { return point; }

    INLINE const Vec<R,SIMD<double>> GetNV () const { return normalvec; }
    INLINE Vec<R,SIMD<double>> & NV () { return normalvec; }    
    INLINE void SetNV (Vec<R,SIMD<double>> vec) { normalvec = vec; }

    INLINE const Vec<R,SIMD<double>> & GetTV () const { return tangentialvec; }
    INLINE void SetTV (Vec<R,SIMD<double>> vec) { tangentialvec = vec; }
  };

  template <int DIMS, int DIMR>
  class SIMD<ngfem::MappedIntegrationPoint<DIMS,DIMR>> : public SIMD<ngfem::DimMappedIntegrationPoint<DIMR>>
  {
    static_assert(DIMS <= DIMR, "DIM-source > DIM-range !!");    
  protected:
    using SIMD<ngfem::DimMappedIntegrationPoint<DIMR>>::measure;
    using SIMD<ngfem::DimMappedIntegrationPoint<DIMR>>::det;
    using SIMD<ngfem::DimMappedIntegrationPoint<DIMR>>::normalvec;
    using SIMD<ngfem::DimMappedIntegrationPoint<DIMR>>::tangentialvec;
    ngbla::Mat<DIMR,DIMS,SIMD<double>> dxdxi;
  public:
    SIMD () = default;
    SIMD (const SIMD<ngfem::IntegrationPoint> & aip,
          const ngfem::ElementTransformation & aeltrans,
          int /* dummy */)
      : SIMD<ngfem::DimMappedIntegrationPoint<DIMR>> (aip, aeltrans)
      { ; }

    SIMD (const SIMD<ngfem::IntegrationPoint> & aip,
          const ngfem::ElementTransformation & aeltrans,
          const Vec<DIMR, SIMD<double>> ax,
          const Mat<DIMR, DIMS, SIMD<double> > & adxdxi)
      : SIMD<ngfem::DimMappedIntegrationPoint<DIMR>> (aip, aeltrans)
    {
      this->point = ax;
      dxdxi = adxdxi;
      Compute();
    }
    
    const Mat<DIMR,DIMS,SIMD<double>> & GetJacobian() const { return dxdxi; }
    Mat<DIMR,DIMS,SIMD<double>> & Jacobian() { return dxdxi; }
    
    int Dim() const { return DIMR; }

    INLINE void Compute ()
    {
      if constexpr (DIMS == DIMR)
	{
	  det = Det (dxdxi);
	  normalvec = SIMD<double>(0.0);
	  tangentialvec = SIMD<double>(0.0);
	}
      else
	{
	  if (DIMR == 3)
	    {
              if (DIMS == 2)
                {
                  normalvec = Cross (Vec<3,SIMD<double>> (dxdxi.Col(0)),
                                     Vec<3,SIMD<double>> (dxdxi.Col(1)));
                  det = L2Norm (normalvec);
                  normalvec /= det;
                  tangentialvec = SIMD<double>(0.0);
                }
              else if (DIMS == 1)
                {
                  normalvec = SIMD<double>(0.0);
		  tangentialvec = Vec<3,SIMD<double>>(dxdxi.Col(0));
		  det = L2Norm(tangentialvec);
		  tangentialvec /= det;
                }
              else
                {
                  det = 1;
                }
            }
	  else if (DIMR == 2)
	    {
              if (DIMS == 1)
                {
                  det = sqrt ( ngstd::sqr (dxdxi(0,0)) + ngstd::sqr (dxdxi(1,0)));
                  
                  normalvec(0) = -dxdxi(1,0) / det;
                  normalvec(1) = dxdxi(0,0) / det;
                  //tangentialvec = SIMD<double>(0.0);
                  tangentialvec(0) = -normalvec(1);
                  tangentialvec(1) = normalvec(0);
                }
              else
                {
                  det = 1;
                }
	    }
	  else
	    {
	      det = 1.0;
	      normalvec = 1.0;
              tangentialvec = SIMD<double>(0.0);
	    }
	}
      measure = fabs (det);
    }
  
    // SIMD<double> GetJacobiDet() const { return det; }
    ///
    // const Mat<DIMS,DIMR,SCAL> & GetJacobianInverse() const { return dxidx; }
    INLINE const Mat<DIMS,DIMR,SIMD<double>> GetJacobianInverse() const 
    { 
      if (DIMS == DIMR)
        return 1.0/det * Trans (Cof (dxdxi));
      else
        {
	  Mat<DIMS,DIMS,SIMD<double>> ata, iata;
	  ata = Trans (dxdxi) * dxdxi;
	  iata = Inv (ata);
	  return (iata * Trans (dxdxi));
        }
    }

    INLINE const Mat<DIMR,DIMS,SIMD<double>> GetJacobianCofactor() const 
    { 
      if (DIMS == DIMR)
        return Cof (dxdxi);
      else
        return det * Trans(GetJacobianInverse());
    }

    
    INLINE operator Vec<DIMS, ngstd::AutoDiff<DIMR,SIMD<double>>> () const
    {
      Vec<DIMS, ngstd::AutoDiff<DIMR, SIMD<double>> > adp;

      Mat<DIMS,DIMR,SIMD<double>> ijac = GetJacobianInverse();
      for (int i = 0; i < DIMS; i++)
        adp[i].Value() = this->IP()(i);
      for (int i = 0; i < DIMS; i++)
        for (int j = 0; j < DIMR; j++)
          adp[i].DValue(j) = ijac(i,j);
      /*
      Mat<DIMS,DIMR,SIMD<double>> ijac = GetJacobianInverse();
      for (int i = 0; i < DIMS; i++)
        adp[i] = AutoDiff<DIMR,SIMD<double>> (this->IP()(i), &ijac(i,0));
      */
      return adp;
    }

    void CalcHesse (Vec<DIMR,Mat<DIMS,DIMS,SIMD<double>>> & ddx1) const;
   
    void Print (ostream & ost) const
    {
      ost << "ip = " << this->ip << std::endl;
      ost << "Point = " << this->point << std::endl;
      ost << "Jacobian = " << dxdxi << std::endl;
      ost << "normal = " << this->GetNV() << std::endl;
    }
  };
}

namespace ngfem
{
  template <int D>
  INLINE auto GetTIP (const IntegrationPoint & ip)
  {
    return ngfem::TIP<D,double> (ip);
  }

  template <int D>
  INLINE auto GetTIPGrad (const IntegrationPoint & ip)
  {
    TIP<D,AutoDiff<D>> tip = ip;
    return tip;
    // return TIP<D,AutoDiff<D>>(ip);
  }

  template <int D>
  INLINE auto GetTIP (const SIMD<IntegrationPoint> & ip)
  {
    // return ngfem::TIP<D,SIMD<double> (ip);
    return ip.TIp<D>();
  }
  
  template <int D>
  INLINE auto GetTIPGrad (const SIMD<IntegrationPoint> & ip)
  {
    Vec<D, AutoDiff<D,SIMD<double>>> adp = ip;
    return TIP<D,AutoDiff<D,SIMD<double>>> (adp, ip.FacetNr(), ip.VB());
    // TIP<D,AutoDiff<D,SIMD<double>>> tip = ip;
    // return tip;
    // return TIP<D,AutoDiff<D>>(ip);
  }

  
  template<int DIMS, int DIMR>
  INLINE void GetTIP1( const SIMD<MappedIntegrationPoint<DIMS,DIMR>> & mip, TIP<DIMS,AutoDiff<DIMR,SIMD<double>>> & adp);

  template<int DIMR>
  INLINE void GetTIP1( const SIMD<MappedIntegrationPoint<0,DIMR>> & mip, TIP<0,AutoDiff<DIMR,SIMD<double>>> & adp) 
  {
    adp.facetnr = mip.IP().FacetNr();
    adp.vb = mip.IP().VB();
  }
  
  template<int DIMR>
  INLINE void GetTIP1 (const SIMD<MappedIntegrationPoint<1,DIMR>> & mip, TIP<1,AutoDiff<DIMR,SIMD<double>>> & adp) 
  {
    Mat<1,DIMR,SIMD<double>> ijac = mip.GetJacobianInverse();
    const auto &ip = mip.IP();
    adp.x.Value() = ip(0);
    for (int j = 0; j < DIMR; j++)
      adp.x.DValue(j) = ijac(0,j);
    adp.facetnr = mip.IP().FacetNr();
    adp.vb = mip.IP().VB();
  }
  
  
  template<int DIMR>
  INLINE void GetTIP1 (const SIMD<MappedIntegrationPoint<2,DIMR>> & mip, TIP<2,AutoDiff<DIMR,SIMD<double>>> & adp) 
  {
    Mat<2,DIMR,SIMD<double>> ijac = mip.GetJacobianInverse();
    const auto &ip = mip.IP();
    adp.x.Value() = ip(0);
    adp.y.Value() = ip(1);
    for (int j = 0; j < DIMR; j++)
      adp.x.DValue(j) = ijac(0,j);
    for (int j = 0; j < DIMR; j++)
      adp.y.DValue(j) = ijac(1,j);
    adp.facetnr = mip.IP().FacetNr();
    adp.vb = mip.IP().VB();
  }
  
  
    template<int DIMR>
    INLINE void GetTIP1 (const SIMD<MappedIntegrationPoint<3,DIMR>> & mip, TIP<3, AutoDiff<DIMR,SIMD<double>>> & adp) 
    {
      Mat<3,DIMR,SIMD<double>> ijac = mip.GetJacobianInverse();
      const auto &ip = mip.IP();
      adp.x.Value() = ip(0);
      adp.y.Value() = ip(1);
      adp.z.Value() = ip(2);
      for (int j = 0; j < DIMR; j++)
        adp.x.DValue(j) = ijac(0,j);
      for (int j = 0; j < DIMR; j++)
        adp.y.DValue(j) = ijac(1,j);
      for (int j = 0; j < DIMR; j++)
        adp.z.DValue(j) = ijac(2,j);
      adp.facetnr = mip.IP().FacetNr();
      adp.vb = mip.IP().VB();
    }






  
  template<int DIMS, int DIMR>
  INLINE auto GetTIP (const SIMD<MappedIntegrationPoint<DIMS,DIMR>> & mip)
  {
    TIP<DIMS,AutoDiff<DIMR,SIMD<double>>> tip(mip.IP().FacetNr(), mip.IP().VB());
    GetTIP1 (mip, tip);
    return tip;
  }

  /*
    template<int DIMR>
    INLINE auto GetTIP( const SIMD<MappedIntegrationPoint<0,DIMR>> & mip) -> TIP<0,AutoDiff<DIMR,SIMD<double>>>
    {
      TIP<0,AutoDiff<DIMR,SIMD<double>>> adp;
      adp.facetnr = mip.IP().FacetNr();
      adp.vb = mip.IP().VB();
      return adp;
    }

    template<int DIMR>
    INLINE auto GetTIP( const SIMD<MappedIntegrationPoint<1,DIMR>> & mip) -> TIP<1,AutoDiff<DIMR,SIMD<double>>>
    {
      TIP<1,AutoDiff<DIMR,SIMD<double>>> adp;
      Mat<1,DIMR,SIMD<double>> ijac = mip.GetJacobianInverse();
      const auto &ip = mip.IP();
      adp.x.Value() = ip(0);
      for (int j = 0; j < DIMR; j++)
        adp.x.DValue(j) = ijac(0,j);
      adp.facetnr = mip.IP().FacetNr();
      adp.vb = mip.IP().VB();
      return adp;
    }

  
    template<int DIMR>
    INLINE auto GetTIP( const SIMD<MappedIntegrationPoint<2,DIMR>> & mip) -> TIP<2,AutoDiff<DIMR,SIMD<double>>>
    {
      TIP<2,AutoDiff<DIMR,SIMD<double>>> adp;      
      Mat<2,DIMR,SIMD<double>> ijac = mip.GetJacobianInverse();
      const auto &ip = mip.IP();
      adp.x.Value() = ip(0);
      adp.y.Value() = ip(1);
      for (int j = 0; j < DIMR; j++)
        adp.x.DValue(j) = ijac(0,j);
      for (int j = 0; j < DIMR; j++)
        adp.y.DValue(j) = ijac(1,j);
      adp.facetnr = mip.IP().FacetNr();
      adp.vb = mip.IP().VB();
      return adp;
    }

  
    template<int DIMR>
    INLINE auto GetTIP(const SIMD<MappedIntegrationPoint<3,DIMR>> & mip) -> TIP<3, AutoDiff<DIMR,SIMD<double>>>
    {
      Mat<3,DIMR,SIMD<double>> ijac = mip.GetJacobianInverse();
      TIP<3, AutoDiff<DIMR,SIMD<double>>> adp;
      const auto &ip = mip.IP();
      adp.x.Value() = ip(0);
      adp.y.Value() = ip(1);
      adp.z.Value() = ip(2);
      for (int j = 0; j < DIMR; j++)
        adp.x.DValue(j) = ijac(0,j);
      for (int j = 0; j < DIMR; j++)
        adp.y.DValue(j) = ijac(1,j);
      for (int j = 0; j < DIMR; j++)
        adp.z.DValue(j) = ijac(2,j);
      adp.facetnr = mip.IP().FacetNr();
      adp.vb = mip.IP().VB();
      return adp;
    }
  */


  template <int DIMS, int DIMR>
  INLINE void GetTIP1( const MappedIntegrationPoint<DIMS,DIMR> & mip);
  
  template<int DIMR>
  INLINE void GetTIP1( const MappedIntegrationPoint<0,DIMR> & mip, TIP<0,AutoDiff<DIMR>> & adp)
  {
    adp.facetnr = mip.IP().FacetNr();
    adp.vb = mip.IP().VB();
  }

  template<int DIMR>
  INLINE void GetTIP1( const MappedIntegrationPoint<1,DIMR> & mip, TIP<1,AutoDiff<DIMR>> & adp)
  {
    Mat<1,DIMR> ijac = mip.GetJacobianInverse();
    const auto &ip = mip.IP();
    adp.x.Value() = ip(0);
    for (int j = 0; j < DIMR; j++)
      adp.x.DValue(j) = ijac(0,j);
    adp.facetnr = mip.IP().FacetNr();
    adp.vb = mip.IP().VB();
  }
  
  template<int DIMR>
  INLINE void GetTIP1( const MappedIntegrationPoint<2,DIMR> & mip, TIP<2,AutoDiff<DIMR>> &adp)
  {
    Mat<2,DIMR> ijac = mip.GetJacobianInverse();
    const auto &ip = mip.IP();
    adp.x.Value() = ip(0);
    adp.y.Value() = ip(1);
    for (int j = 0; j < DIMR; j++)
      adp.x.DValue(j) = ijac(0,j);
    for (int j = 0; j < DIMR; j++)
      adp.y.DValue(j) = ijac(1,j);
    adp.facetnr = mip.IP().FacetNr();
    adp.vb = mip.IP().VB();
  }

  
  template<int DIMR>
  INLINE void GetTIP1(const MappedIntegrationPoint<3,DIMR> & mip, TIP<3, AutoDiff<DIMR>> & adp)
  {
    Mat<3,DIMR> ijac = mip.GetJacobianInverse();
    const auto &ip = mip.IP();
    adp.x.Value() = ip(0);
    adp.y.Value() = ip(1);
    adp.z.Value() = ip(2);
    for (int j = 0; j < DIMR; j++)
      adp.x.DValue(j) = ijac(0,j);
    for (int j = 0; j < DIMR; j++)
      adp.y.DValue(j) = ijac(1,j);
    for (int j = 0; j < DIMR; j++)
      adp.z.DValue(j) = ijac(2,j);
    adp.facetnr = mip.IP().FacetNr();
    adp.vb = mip.IP().VB();
  }

  

  template<int DIMS, int DIMR>
  INLINE auto GetTIP (const MappedIntegrationPoint<DIMS,DIMR> & mip) //  -> TIP<DIMS,AutoDiff<DIMR>>;
  {
    TIP<DIMS,AutoDiff<DIMR>> tip(mip.IP().FacetNr(), mip.IP().VB());
    GetTIP1(mip, tip);
    return tip;
  }


  template <int DIMS, int DIMR>
  auto GetTIPHesse (const MappedIntegrationPoint<DIMS, DIMR> & mip)
  {
    Vec<DIMR, Mat<DIMS,DIMS>> hesse;
    mip.CalcHesse(hesse);
    Mat<DIMS,DIMR> ijac = mip.GetJacobianInverse();
    
    Vec<DIMR, Mat<DIMR,DIMR>> hessemapped;
    for (int i = 0; i < DIMR; i++)
      hessemapped(i) =  Trans(ijac) * hesse(i) * ijac; 

    Vec<DIMS, Mat<DIMR,DIMR>> hessemapped2 = Mat<DIMR,DIMR>(0.0);
    for (int i =  0; i < DIMR; i++)
      for (int j = 0; j < DIMS; j++)
        hessemapped2(j) += ijac(j,i) * hessemapped(i);

    TIP<DIMS, AutoDiffDiff<DIMR>> tip = GetTIP(mip);
    
    if constexpr(DIMS >= 1)
      for (int j = 0; j < DIMR; j++)
        for (int k = 0; k < DIMR; k++)
          tip.x.DDValue(j,k) = -hessemapped2(0)(j,k);
    if constexpr(DIMS >= 2)
      for (int j = 0; j < DIMR; j++)
        for (int k = 0; k < DIMR; k++)
          tip.y.DDValue(j,k) = -hessemapped2(1)(j,k);
    if constexpr(DIMS >= 3)
      for (int j = 0; j < DIMR; j++)
        for (int k = 0; k < DIMR; k++)
          tip.z.DDValue(j,k) = -hessemapped2(2)(j,k);
    
    return tip;
  }


  template <int DIMS, int DIMR>
  auto GetTIPHesse (const SIMD<MappedIntegrationPoint<DIMS, DIMR>> & mip)
  {
    Vec<DIMR, Mat<DIMS,DIMS, SIMD<double>>> hesse;
    mip.CalcHesse(hesse);
    Mat<DIMS,DIMR, SIMD<double>> ijac = mip.GetJacobianInverse();
    
    Vec<DIMR, Mat<DIMR,DIMR, SIMD<double>>> hessemapped;
    for (int i = 0; i < DIMR; i++)
      hessemapped(i) =  Trans(ijac) * hesse(i) * ijac; 

    Vec<DIMS, Mat<DIMR,DIMR, SIMD<double>>> hessemapped2 = Mat<DIMR,DIMR,SIMD<double>>(0.0);
    for (int i =  0; i < DIMR; i++)
      for (int j = 0; j < DIMS; j++)
        hessemapped2(j) += ijac(j,i) * hessemapped(i);

    TIP<DIMS, AutoDiffDiff<DIMR, SIMD<double>>> tip = GetTIP(mip);
    
    if constexpr(DIMS >= 1)
      for (int j = 0; j < DIMR; j++)
        for (int k = 0; k < DIMR; k++)
          tip.x.DDValue(j,k) = -hessemapped2(0)(j,k);
    if constexpr(DIMS >= 2)
      for (int j = 0; j < DIMR; j++)
        for (int k = 0; k < DIMR; k++)
          tip.y.DDValue(j,k) = -hessemapped2(1)(j,k);
    if constexpr(DIMS >= 3)
      for (int j = 0; j < DIMR; j++)
        for (int k = 0; k < DIMR; k++)
          tip.z.DDValue(j,k) = -hessemapped2(2)(j,k);
    
    return tip;
  }


  
  /*
  template<int DIMR>
  INLINE auto GetTIP( const MappedIntegrationPoint<0,DIMR> & mip) -> TIP<0,AutoDiff<DIMR>>
  {
    TIP<0,AutoDiff<DIMR>> adp;
    adp.facetnr = mip.IP().FacetNr();
    adp.vb = mip.IP().VB();
    return adp;
  }

    template<int DIMR>
    INLINE auto GetTIP( const MappedIntegrationPoint<1,DIMR> & mip) -> TIP<1,AutoDiff<DIMR>>
    {
      TIP<1,AutoDiff<DIMR>> adp;
      Mat<1,DIMR> ijac = mip.GetJacobianInverse();
      const auto &ip = mip.IP();
      adp.x.Value() = ip(0);
      for (int j = 0; j < DIMR; j++)
        adp.x.DValue(j) = ijac(0,j);
      adp.facetnr = mip.IP().FacetNr();
      adp.vb = mip.IP().VB();
      return adp;
    }

  
    template<int DIMR>
    INLINE auto GetTIP( const MappedIntegrationPoint<2,DIMR> & mip) -> TIP<2,AutoDiff<DIMR>>
    {
      TIP<2,AutoDiff<DIMR>> adp;      
      Mat<2,DIMR> ijac = mip.GetJacobianInverse();
      const auto &ip = mip.IP();
      adp.x.Value() = ip(0);
      adp.y.Value() = ip(1);
      for (int j = 0; j < DIMR; j++)
        adp.x.DValue(j) = ijac(0,j);
      for (int j = 0; j < DIMR; j++)
        adp.y.DValue(j) = ijac(1,j);
      adp.facetnr = mip.IP().FacetNr();
      adp.vb = mip.IP().VB();
      return adp;
    }

  
    template<int DIMR>
    INLINE auto GetTIP(const MappedIntegrationPoint<3,DIMR> & mip) -> TIP<3, AutoDiff<DIMR>>
    {
      Mat<3,DIMR> ijac = mip.GetJacobianInverse();
      TIP<3, AutoDiff<DIMR>> adp;
      const auto &ip = mip.IP();
      adp.x.Value() = ip(0);
      adp.y.Value() = ip(1);
      adp.z.Value() = ip(2);
      for (int j = 0; j < DIMR; j++)
        adp.x.DValue(j) = ijac(0,j);
      for (int j = 0; j < DIMR; j++)
        adp.y.DValue(j) = ijac(1,j);
      for (int j = 0; j < DIMR; j++)
        adp.z.DValue(j) = ijac(2,j);
      adp.facetnr = mip.IP().FacetNr();
      adp.vb = mip.IP().VB();
      return adp;
    }
  */

  



  
  class SIMD_IntegrationRule : public Array<SIMD<IntegrationPoint>>
  {
    int dimension = -1;
    size_t nip = -47;
    const SIMD_IntegrationRule *irx = nullptr, *iry = nullptr, *irz = nullptr; // for tensor product IR
  public:
    SIMD_IntegrationRule () = default;
    inline SIMD_IntegrationRule (ELEMENT_TYPE eltype, int order);
    SIMD_IntegrationRule (const SIMD_IntegrationRule & ir) = delete;
    SIMD_IntegrationRule (SIMD_IntegrationRule && ir) = default;
    NGS_DLL_HEADER SIMD_IntegrationRule (const IntegrationRule & ir);
    NGS_DLL_HEADER SIMD_IntegrationRule (const IntegrationRule & ir, LocalHeap & lh);
    NGS_DLL_HEADER SIMD_IntegrationRule (int nip, LocalHeap & lh);
    NGS_DLL_HEADER ~SIMD_IntegrationRule ()
    {
      delete [] mem_to_delete;
      mem_to_delete = nullptr;
    }

    SIMD_IntegrationRule & operator= (const SIMD_IntegrationRule &) = delete;
    SIMD_IntegrationRule & operator= (SIMD_IntegrationRule &&) = default;
    
    SIMD_IntegrationRule (size_t asize, SIMD<IntegrationPoint> * pip)
      : Array<SIMD<IntegrationPoint>> (asize, pip), nip(asize*SIMD<IntegrationPoint>::Size()) { }

    INLINE SIMD_IntegrationRule Range (size_t first, size_t next) const
    {
      return SIMD_IntegrationRule (next-first, &(*this)[first]);
    }

    size_t GetNIP() const { return nip; } // Size()*SIMD<double>::Size(); }
    void SetNIP(size_t _nip) { nip = _nip; }

    SIMD_IntegrationRule Clone() const
    {
      SIMD_IntegrationRule ir2(Size(), &(*this)[0]);
      ir2.dimension = dimension;
      ir2.nip = nip;
      ir2.irx = irx;
      ir2.iry = iry;
      ir2.irz = irz;
      return ir2;
    }


    
    bool IsTP() const { return irx != nullptr; } 
    const SIMD_IntegrationRule & GetIRX() const { return *irx; }
    const SIMD_IntegrationRule & GetIRY() const { return *iry; }
    const SIMD_IntegrationRule & GetIRZ() const { return *irz; }
    void SetIRX(const SIMD_IntegrationRule * ir) { irx = ir; }
    void SetIRY(const SIMD_IntegrationRule * ir) { iry = ir; }
    void SetIRZ(const SIMD_IntegrationRule * ir) { irz = ir; }
  };

  extern NGS_DLL_HEADER const SIMD_IntegrationRule & SIMD_SelectIntegrationRule (ELEMENT_TYPE eltype, int order);

  inline SIMD_IntegrationRule :: SIMD_IntegrationRule (ELEMENT_TYPE eltype, int order)
  { 
    const SIMD_IntegrationRule & ir = SIMD_SelectIntegrationRule (eltype, order);
    size = ir.Size();
    data = &ir[0];
    mem_to_delete = NULL;
    dimension = ElementTopology::GetSpaceDim(eltype);
    nip = ir.nip;
    irx = ir.irx;
    iry = ir.iry;
    irz = ir.irz;
  }



  
  class NGS_DLL_HEADER SIMD_BaseMappedIntegrationRule 
  { 
  protected:
    SIMD_IntegrationRule ir;
    const ElementTransformation & eltrans;
    char * baseip;
    size_t incr;
    int dim_element, dim_space;

    // mir on other element as needed for evaluating DG jump terms
    const SIMD_BaseMappedIntegrationRule * other_mir = nullptr;

    BareSliceMatrix<SIMD<double>> points{0,0,0, nullptr};
    BareSliceMatrix<SIMD<double>> normals{0,0,0, nullptr};
  public:
    SIMD_BaseMappedIntegrationRule (const SIMD_IntegrationRule & air,
                                    const ElementTransformation & aeltrans)
      : ir(air.Size(),&air[0]), eltrans(aeltrans)
    {
      ir.SetIRX(&air.GetIRX());
      ir.SetIRY(&air.GetIRY());
      ir.SetIRZ(&air.GetIRZ());
      ir.SetNIP(air.GetNIP());
    }
    ~SIMD_BaseMappedIntegrationRule ()
      { ir.NothingToDelete(); }

    INLINE size_t Size() const { return ir.Size(); }
    INLINE const SIMD_IntegrationRule & IR() const { return ir; }
    INLINE const ElementTransformation & GetTransformation () const { return eltrans; }
    virtual void ComputeNormalsAndMeasure (ELEMENT_TYPE et, int facetnr) = 0;    
    INLINE const SIMD<BaseMappedIntegrationPoint> & operator[] (size_t i) const
    { return *static_cast<const SIMD<BaseMappedIntegrationPoint>*> ((void*)(baseip+i*incr)); }
    INLINE int DimElement() const { return dim_element; }
    INLINE int DimSpace() const { return dim_space; }
    // virtual ABareMatrix<double> GetPoints() const = 0;
    // virtual BareSliceMatrix<SIMD<double>> GetPoints() const = 0;
    SliceMatrix<SIMD<double>> GetPoints() const { return points.AddSize(ir.Size(), dim_space) ; }
    SliceMatrix<SIMD<double>> GetNormals() const { return normals.AddSize(ir.Size(), dim_space) ; }
    virtual void Print (ostream & ost) const = 0;
    bool IsComplex() const { return false; }
    BareSliceMatrix<SIMD<Complex>> GetPointsComplex() const { throw ExceptionNOSIMD("Not implemented"); }

    // for DG jump terms
    void SetOtherMIR (const SIMD_BaseMappedIntegrationRule * other) { other_mir = other; }
    auto GetOtherMIR () const { return other_mir; }

    virtual void TransformGradient (BareSliceMatrix<SIMD<double>> grad) const = 0; // covariant transformation
    virtual void TransformGradientTrans (BareSliceMatrix<SIMD<double>> grad) const = 0; // covariant transformation transpose
  };

  inline ostream & operator<< (ostream & ost, const SIMD_BaseMappedIntegrationRule & mir)
  {
    mir.Print (ost);
    return ost;
  }

  template <int DIM_ELEMENT, int DIM_SPACE>
  class NGS_DLL_HEADER SIMD_MappedIntegrationRule : public SIMD_BaseMappedIntegrationRule
  {
    FlatArray<SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>> mips;
  public:
    SIMD_MappedIntegrationRule (const SIMD_IntegrationRule & ir, 
                                const ElementTransformation & aeltrans, 
                                Allocator & lh);
    SIMD_MappedIntegrationRule (const SIMD_IntegrationRule & ir, 
                                const ElementTransformation & aeltrans,
                                int dummy, 
                                Allocator & lh)
      : SIMD_BaseMappedIntegrationRule (ir, aeltrans), mips(ir.Size(), lh)
      {
        dim_element = DIM_ELEMENT;
        dim_space = DIM_SPACE;
        baseip = (char*)(void*)(SIMD<BaseMappedIntegrationPoint>*)(&mips[0]);
        incr = sizeof (SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>);

        for (size_t i = 0; i < ir.Size(); i++)
          new (&mips[i]) SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>> (ir[i], eltrans, -1);

        /*
        new (&points) BareSliceMatrix<SIMD<double>> (sizeof(SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>)/sizeof(SIMD<double>),
                                                     &mips[0].Point()(0),
                                                     DummySize(mips.Size(), DIM_SPACE));
        
        new (&normals) BareSliceMatrix<SIMD<double>> (sizeof(SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>)/sizeof(SIMD<double>),
                                                      &mips[0].NV()(0),
                                                      DummySize(mips.Size(), DIM_SPACE));
        */
        new (&points) BareSliceMatrix<SIMD<double>> (mips.Size(), DIM_SPACE,
                                                     sizeof(SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>)/sizeof(SIMD<double>),
                                                     &mips[0].Point()(0));
        
        new (&normals) BareSliceMatrix<SIMD<double>> (mips.Size(), DIM_SPACE,
                                                      sizeof(SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>)/sizeof(SIMD<double>),
                                                      &mips[0].NV()(0));
      }

    virtual void ComputeNormalsAndMeasure (ELEMENT_TYPE et, int facetnr) override;
    SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>> & operator[] (size_t i) const 
    { 
      return mips[i]; 
    }
    /*
    virtual BareSliceMatrix<SIMD<double>> GetPoints() const
    {
      return BareSliceMatrix<SIMD<double>> (sizeof(SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>)/sizeof(SIMD<double>),
                                            &mips[0].Point()(0),
                                            DummySize(mips.Size(), DIM_SPACE));
    }
    */
    virtual void Print (ostream & ost) const override;

    virtual void TransformGradient (BareSliceMatrix<SIMD<double>> grad) const override;
    virtual void TransformGradientTrans (BareSliceMatrix<SIMD<double>> grad) const override;
  };
}


#endif
