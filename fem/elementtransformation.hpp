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

#include "intrule.hpp"

namespace ngcomp { class GridFunction; }

namespace ngfem
{
  class FiniteElement;
  template <int D> class ScalarFiniteElement;
  
  /**
     Transformation from reference element to physical element.
  */
  class NGS_DLL_HEADER ElementTransformation
  {
  protected:
    /// geometry of element
    ELEMENT_TYPE eltype;
    /// element number
    int elnr;
    /// material property
    int elindex;
    ///
    bool higher_integration_order;
    /// is the element curved ?
    bool iscurved = false;
    bool is_complex = false;
  public:
    ///
    // ElementTransformation () { higher_integration_order = false; } 
    ElementTransformation (ELEMENT_TYPE et, VorB /* aboundary */, int aelnr, int aelindex)
      : eltype(et), elnr(aelnr), elindex(aelindex)
    {
      higher_integration_order = false; 
    }

    ElementTransformation(ELEMENT_TYPE et, ElementId ei, int aelindex)
      : eltype(et), elnr(ei.Nr()), elindex(aelindex)
    {
      higher_integration_order = false;
    }
    
    ///
    virtual ~ElementTransformation() { ; } 
    /*
    /// set data: is it a boundary, element number, and element index
    virtual void SetElement (bool aboundary, int aelnr, int aelindex)
    {
      elnr = aelnr;
      elindex = aelindex;
    }
    */
    /// set geometric type of element
    // void SetElementType (ELEMENT_TYPE aet) { eltype = aet; }
    /// return element number
    int GetElementNr () const { return elnr; }
    void SetElementNr (int _elnr) { elnr = _elnr; }
    ///
    ElementId GetElementId () const { return ElementId(VB(), elnr); }
    /// return element index
    int GetElementIndex () const { return elindex; }
    /// return element geometry type 
    ELEMENT_TYPE GetElementType () const { return eltype; }
    /// calculates the Jacobi matrix
    virtual void CalcJacobian (const IntegrationPoint & ip,
			       FlatMatrix<> dxdxi) const = 0;

    /// calculate the mapped point
    virtual void CalcPoint (const IntegrationPoint & ip,
			    FlatVector<> point) const = 0;

    /// calculate point and Jacobi matrix
    virtual void CalcPointJacobian (const IntegrationPoint & ip,
				    FlatVector<> point, FlatMatrix<> dxdxi) const = 0;

    /// Calculate points and Jacobimatrices in all points of integrationrule
    virtual void CalcMultiPointJacobian (const IntegrationRule & ir,
					 BaseMappedIntegrationRule & mir) const = 0;

    /// Calculate points and Jacobimatrices in all points of integrationrule
    virtual void CalcMultiPointJacobian (const SIMD_IntegrationRule & ir,
					 SIMD_BaseMappedIntegrationRule & mir) const;
    
    /// Calcs the normal vector in ip
    void CalcNormalVector (const IntegrationPoint & ip,
			   FlatVector<> nv,
			   LocalHeap & lh) const
    {
      if(VB()==BND)
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
    virtual int SpaceDim () const = 0;

    /// is it a mapping for boundary or codim 2 elements ?
    virtual VorB VB() const = 0;

    virtual int ElementDim () const
    { return SpaceDim() - int(VB()); }
    
    void SetHigherIntegrationOrder(void) {higher_integration_order = true;}
    void UnSetHigherIntegrationOrder(void) {higher_integration_order = false;}
    bool HigherIntegrationOrderSet(void) const 
    {
      return higher_integration_order;
    }

    /// has the element non-constant Jacobian ?
    virtual bool IsCurvedElement() const 
    {
      return iscurved;
    }

    bool IsComplex() const { return is_complex; }
    
    virtual void GetSort (FlatArray<int> sort) const
    { ; }

    /// return a mapped integration point on localheap
    virtual BaseMappedIntegrationPoint & operator() (const IntegrationPoint & ip, Allocator & lh) const = 0;

    /// return a mapped integration rule on localheap
    virtual BaseMappedIntegrationRule & operator() (const IntegrationRule & ir, Allocator & lh) const = 0;

    /// return a mapped integration rule on localheap
    virtual SIMD_BaseMappedIntegrationRule & operator() (const SIMD_IntegrationRule & ir, Allocator & lh) const;

    template <int DIMS, int DIMR> 
      void CalcHesse (const SIMD<ngfem::IntegrationPoint> & ip, Vec<DIMR, Mat<DIMS,DIMS,SIMD<double>>> & hesse) const
    {
      VCalcHesse (ip, &hesse(0)(0,0));
    }
    virtual void VCalcHesse (const SIMD<ngfem::IntegrationPoint> & ip, SIMD<double> * hesse) const;
        
    virtual bool BelongsToMesh (const void * mesh) const { return true; }
    virtual const void * GetMesh () const { return NULL; }

    const ElementTransformation & AddDeformation (const ngcomp::GridFunction * gf, LocalHeap & lh) const
    {
      if (!gf) return *this;
      return VAddDeformation(gf,lh);
    }
    virtual const ElementTransformation & VAddDeformation (const ngcomp::GridFunction * gf, LocalHeap & lh) const
    {
      throw Exception ("don't know how to deform");
    }
    
    void * userdata = nullptr;
    class UserDataStack {
      ElementTransformation & trafo;
      void * save;
    public:
      UserDataStack (ElementTransformation & atrafo)
        : trafo(atrafo)
      {
        save = trafo.userdata;
        trafo.userdata = nullptr;
      }
      ~UserDataStack ()
      {
        trafo.userdata = save;
      }
    };
    [[nodiscard]]
      auto PushUserData() const
      {
        return UserDataStack(const_cast<ElementTransformation&>(*this));
      }
    
  private:
    ElementTransformation (const ElementTransformation & eltrans2) { ; }
    ElementTransformation & operator= (const ElementTransformation & eltrans2) 
    { return *this; }
  };



  
  int BaseMappedIntegrationPoint :: DimElement() const { return eltrans->ElementDim(); }
  int BaseMappedIntegrationPoint :: DimSpace() const { return eltrans->SpaceDim(); } 



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
    Matrix<> pointmat;
    /// normal vectors (only surfelements)
    FlatMatrix<> nvmat;
  public:
    /// type of element, np x dim point-matrix
    FE_ElementTransformation (ELEMENT_TYPE type, SliceMatrix<> pmat);
    /// trafo of reference-element
    FE_ElementTransformation (ELEMENT_TYPE type);
    ///
    FE_ElementTransformation ();

    ///
    ~FE_ElementTransformation ();

    ///
    virtual void SetElement (const ScalarFiniteElement<DIMS> * afel, int aelnr, int aelindex)
    {
      // fel = static_cast<const ScalarFiniteElement<DIMS>*> (afel);
      fel = afel; 
      elnr = aelnr; 
      elindex = aelindex;
      eltype = fel->ElementType();
      pointmat.SetSize (DIMR, fel->GetNDof());
    }

    ///
    const FiniteElement & GetElement () const { return *fel; }

  
    ELEMENT_TYPE GetElementType () const { return fel->ElementType(); }
  
    ///
    virtual void CalcJacobian (const IntegrationPoint & ip,
			       FlatMatrix<> dxdxi) const override;

    ///
    virtual void CalcPoint (const IntegrationPoint & ip, 
			    FlatVector<> point) const override;

    ///
    virtual void CalcPointJacobian (const IntegrationPoint & ip,
				    FlatVector<> point, 
				    FlatMatrix<> dxdxi) const override;

    virtual void CalcMultiPointJacobian (const IntegrationRule & ir,
					 BaseMappedIntegrationRule & bmir) const override;

    virtual void CalcMultiPointJacobian (const SIMD_IntegrationRule & ir,
					 SIMD_BaseMappedIntegrationRule & mir) const override;

    ///
    FlatMatrix<> PointMatrix () const { return pointmat; }
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
    int SpaceDim () const override
    {
      return pointmat.Height(); 
    }

    VorB VB(void) const override
    {
      return pointmat.Height()==ElementTopology::GetSpaceDim(fel->ElementType()) ? VOL :
	(pointmat.Height()==ElementTopology::GetSpaceDim(fel->ElementType())-1 ? BND : BBND);
    }

    void GetSort (FlatArray<int> sort) const override
    { ; }


    virtual BaseMappedIntegrationPoint & operator() (const IntegrationPoint & ip, Allocator & lh) const override
    {
      return *new (lh) MappedIntegrationPoint<DIMS,DIMR> (ip, *this);
    }

    virtual BaseMappedIntegrationRule & operator() (const IntegrationRule & ir, Allocator & lh) const override
    {
      return *new (lh) MappedIntegrationRule<DIMS,DIMR> (ir, *this, lh);
    }

    virtual SIMD_BaseMappedIntegrationRule & operator() (const SIMD_IntegrationRule & ir, Allocator & lh) const override
    {
      return *new (lh) SIMD_MappedIntegrationRule<DIMS,DIMR> (ir, *this, lh);
    }

  };

  extern NGS_DLL_HEADER const ElementTransformation & GetFEElementTransformation (ELEMENT_TYPE et);

}


#endif






