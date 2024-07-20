#ifndef FILE_INTEGRATOR
#define FILE_INTEGRATOR

/*********************************************************************/
/* File:   integrator.hpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


// #include "elementtransformation.hpp"
#include "finiteelement.hpp"
#include "differentialoperator.hpp"

namespace ngfem
{

  /*
    bilinear-form and linear-form integrators
  */


  // enum FORM_TYPE { VOLUME_FORM, BOUNDARY_FORM, SKELETON_FORM, CURVE_FORM };


  /**
     Base class for linear-form and bilinear-form integrators.
     Provides integration order, restriction to subdomains
  */
  class NGS_DLL_HEADER Integrator
  {
  protected:
    /// define only on some sub-domains
    BitArray definedon;

    /// if >= 0, use this order of integration
    int integration_order;

    // lower bound for the integration order
    int higher_integration_order;

    /// if >= 0, use this order of integration for all terms
    static int common_integration_order;

    /// plane element and constant coefficients 
    bool const_coef;

    ///
    string name;

    /// integration only along curve
    bool is_curve_integrator = false;
    Array < FlatVector < double > * > curve_ips;
    Array < FlatVector < double > * > curve_ip_tangents;
    Array <int> continuous_curveparts;
  
    int cachecomp;
    int bonus_intorder = 0;

    /// define only on some sub-domains
    shared_ptr<BitArray> definedon_element = nullptr;
    std::array<unique_ptr<IntegrationRule>,25> userdefined_intrules;
    std::array<unique_ptr<SIMD_IntegrationRule>,25> userdefined_simd_intrules;

    mutable bool simd_evaluate = true;

    shared_ptr<ngcomp::GridFunction> deformation; // ALE for this integrator
    
  protected:
    void DeleteCurveIPs ( void );

  public:
    /// constructor
    Integrator() throw ();

    /// destructor
    virtual ~Integrator();

    /// integrates on the boundary, or on the domain ?
    //[[deprecated("Use VB() instead")]]
    virtual bool BoundaryForm () const
      { return VB() == BND; }

    virtual VorB VB() const = 0;

    class DGFormulation
    {
    public:
      bool neighbor_testfunction = true; // trivially parallel if no neighbor testfunction
      bool element_boundary = false;   // loop over edges, or loop over element-boundaries
      DGFormulation(bool nbtest=true, bool eb=false)
        : neighbor_testfunction(nbtest), element_boundary(eb) { ; }
    };
    virtual DGFormulation GetDGFormulation() const { return DGFormulation(); }
      
    /// integrates just on the skeleton, standard is NO
    virtual bool SkeletonForm () const { return false; }  

    virtual bool VolumeForm () const
    {
      if ( VB() != VOL || SkeletonForm() || IntegrationAlongCurve()) return false;
      return true;
    }

    /// Is Integrator defined on this sub-domain ?
    bool DefinedOn (int mat) const;

    /// defined only on some subdomains
    void SetDefinedOn (const BitArray & adefinedon);

    /// defined only on elements (elements/boundary elements/facets/..)
    void SetDefinedOnElements (shared_ptr<BitArray> adefinedonelem)
    {
      definedon_element = adefinedonelem;
    }

    void SetIntegrationRule(ELEMENT_TYPE et, const IntegrationRule& ir)
    {
      userdefined_intrules[int(et)] = make_unique<IntegrationRule>(ir.Copy());
      userdefined_simd_intrules[int(et)] = make_unique<SIMD_IntegrationRule>(*userdefined_intrules[int(et)]);
    }

    void SetIntegrationRule(const IntegrationRule& ir)
    {
      for(auto i : Range(25))
        {
          userdefined_intrules[i] = make_unique<IntegrationRule>(ir.Copy());
          userdefined_simd_intrules[i] = make_unique<SIMD_IntegrationRule>(*userdefined_intrules[i]);
        }
    }

    inline const IntegrationRule& GetIntegrationRule(ELEMENT_TYPE et, int order) const
    {
      return userdefined_intrules[et] ? *userdefined_intrules[et] : SelectIntegrationRule(et,order);
    }
    inline const SIMD_IntegrationRule& GetSIMDIntegrationRule(ELEMENT_TYPE et, int order) const
    {
      return userdefined_simd_intrules[et] ? *userdefined_simd_intrules[et] : SIMD_SelectIntegrationRule(et,order);
    }

    /// defined only on some elements/facets/boundary elements
    shared_ptr<BitArray> GetDefinedOnElements () const { return definedon_element; } 
    
    /// Is Integrator defined on this element ?
    bool DefinedOnElement (int elem) const
    {
      if (definedon_element != nullptr && !definedon_element->Test(elem))
        return false;
      else
        return true;
    }
    
    /// defined only on some subdomains
    const BitArray & GetDefinedOn () const { return definedon; } 

    /// defined only on some subdomains (0-based)
    void SetDefinedOn (const Array<int> & regions);

    bool DefinedOnSubdomainsOnly() const
    { return definedon.Size() != 0; }



    /// use exactly this integration order for all integrals
    static void SetCommonIntegrationOrder (int cio)
    { 
      common_integration_order = cio; 
    }

    static int GetCommonIntegrationOrder ()
    { 
      return common_integration_order;
    }

    /// set minimal integration order
    void SetHigherIntegrationOrder (int io)
    {
      higher_integration_order = io; 
    }

    /// set integration order
    void SetIntegrationOrder (int io)
    {
      integration_order = io; 
    }

    /// returns integration order
    int GetIntegrationOrder (void) const
    {
      return integration_order;
    }

    void SetBonusIntegrationOrder (int bo) { bonus_intorder = bo; }
    int GetBonusIntegrationOrder() const   { return bonus_intorder; }
    
    /// benefit from constant coefficient
    void SetConstantCoefficient (bool acc = 1)
    { const_coef = acc; }

    /// dimension of element
    virtual int DimElement () const { return -1; }

    /// dimension of space
    virtual int DimSpace () const { return -1; }
    
    /// how many copies of the same element (e.g. elasticity)
    virtual int GetDimension () const { return 1; }
    
    /// 
    void SetName (const string & aname);
    ///
    virtual string Name () const;

    /// does element match integrator ?
    virtual void CheckElement (const FiniteElement & el) const { ; }



    // special hacks by Markus
    bool IntegrationAlongCurve (void) const
    { return is_curve_integrator; }

    void SetIntegrationAlongCurve ( const int npoints );

    void UnSetIntegrationAlongCurve ( void );

    virtual int NumCurvePoints() const
    { return curve_ips.Size(); }

    virtual FlatVector<double> CurvePoint(int i)
    { return *(curve_ips[i]); }

    virtual FlatVector<double> CurvePointTangent(int i)
    { return *(curve_ip_tangents[i]); }

    virtual int GetNumCurveParts() const;
    virtual int GetStartOfCurve(int i) const;
    virtual int GetEndOfCurve(int i) const;

    virtual void AppendCurvePoint(const FlatVector<double> & point);
    virtual void AppendCurvePoint(const FlatVector<double> & point, const FlatVector<double> & tangent);
    virtual void SetCurveClearance();


  
    virtual void SetCacheComp(const int comp)
    { cachecomp = comp; }

    virtual int CacheComp(void) const
    { return cachecomp; }

    virtual void SetFileName(const string & filename);

    virtual void SetFlags (const Flags & flags) { ; }

    bool SimdEvaluate () const { return simd_evaluate; }
    void SetSimdEvaluate (bool b = true) { simd_evaluate = b; }

    void SetDeformation (shared_ptr<ngcomp::GridFunction> adeform) { deformation = adeform; } 
    const shared_ptr<ngcomp::GridFunction> & GetDeformation() const { return deformation; }
  };


  ostream & operator << (ostream & ost, const Integrator & igt);

  /**
     A BilinearFormIntegrator computes the element matrices. Different
     equations are provided by derived classes. An Integrator can be defined
     in the domain or at the boundary.
  */
  class NGS_DLL_HEADER BilinearFormIntegrator : public Integrator
  {
  protected:
    // evaluate something, e.g. energy, ...
    SymbolTable<shared_ptr<DifferentialOperator>> evaluators;

  public:
    // typedef double TSCAL;
    ///
    BilinearFormIntegrator () throw () { ; }
    ///
    virtual ~BilinearFormIntegrator ();

    /// generates symmetric matrix ? 
    virtual xbool IsSymmetric () const = 0;

    /// components of flux
    virtual int DimFlux () const { return -1; }

    /**
       Computes the element matrix.
    */
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
		       const ElementTransformation & eltrans, 
		       FlatMatrix<double> elmat,
		       LocalHeap & lh) const = 0;

    /**
       Computes the element matrix.
       Complex version
    */
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
		       const ElementTransformation & eltrans, 
		       FlatMatrix<Complex> elmat,
		       LocalHeap & lh) const;

    /**
       Computes the element matrix.
       Add the element to elmat
    */
    virtual void
      CalcElementMatrixAdd (const FiniteElement & fel,
                            const ElementTransformation & eltrans, 
                            FlatMatrix<double> elmat,
                            bool & symmetric_so_far,
                            LocalHeap & lh) const;

    /**
       Computes the element matrix.
       Complex version
       Add the element to elmat
    */
    virtual void
      CalcElementMatrixAdd (const FiniteElement & fel,
                            const ElementTransformation & eltrans, 
                            FlatMatrix<Complex> elmat,
                            bool & symmetric_so_far,                            
                            LocalHeap & lh) const;
    

    
    virtual void
    CalcElementMatrixIndependent (const FiniteElement & bfel_master,
                                  const FiniteElement & bfel_master_element,				    
                                  const FiniteElement & bfel_other,
                                  const ElementTransformation & eltrans_master, 
                                  const ElementTransformation & eltrans_master_element, 
                                  const ElementTransformation & eltrans_other,
                                  const IntegrationPoint & ip_master,
                                  const IntegrationPoint & ip_master_element,
                                  const IntegrationPoint & ip_other,
                                  FlatMatrix<double> & elmat,
                                  LocalHeap & lh) const
    {;}
    virtual void
    ApplyElementMatrixIndependent (const FiniteElement & bfel_master,
				   const FiniteElement & bfel_master_element,				    
				   const FiniteElement & bfel_other,
				   const ElementTransformation & eltrans_master, 
				   const ElementTransformation & eltrans_master_element, 
				   const ElementTransformation & eltrans_other,
				   const IntegrationPoint & ip_master,
				   const IntegrationPoint & ip_master_element,
				   const IntegrationPoint & ip_other,
				   const FlatVector<double> & elx,
				   Vector<double> & result,
				   LocalHeap & lh) const
    {;}
    virtual void
    CalcElementMatrixIndependent (const FiniteElement & bfel_master,
                                  const FiniteElement & bfel_master_element,				    
                                  const FiniteElement & bfel_other,
                                  const ElementTransformation & eltrans_master, 
                                  const ElementTransformation & eltrans_master_element, 
                                  const ElementTransformation & eltrans_other,
                                  const IntegrationPoint & ip_master,
                                  const IntegrationPoint & ip_master_element,
                                  const IntegrationPoint & ip_other,
                                  FlatMatrix<Complex> & elmat,
                                  LocalHeap & lh) const
    {
      FlatMatrix<double> rmat;
      CalcElementMatrixIndependent(bfel_master,bfel_master_element,bfel_other,
                                   eltrans_master, eltrans_master_element, eltrans_other,
                                   ip_master, ip_master_element, ip_other,
                                   rmat, lh);
      elmat.AssignMemory(rmat.Height(), rmat.Width(), lh);
      elmat = rmat;
    }

    virtual void
    CalcElementMatrixIndependent (const FiniteElement & bfel_master,
				      const FiniteElement & bfel_other,
				      const ElementTransformation & eltrans_master, 
				      const ElementTransformation & eltrans_other,
				      const IntegrationPoint & ip_master,
				      const IntegrationPoint & ip_other,
				      FlatMatrix<double> & elmat,
				      LocalHeap & lh) const
    {;}
    virtual void
    CalcElementMatrixIndependent (const FiniteElement & bfel_master,
				      const FiniteElement & bfel_other,
				      const ElementTransformation & eltrans_master, 
				      const ElementTransformation & eltrans_other,
				      const IntegrationPoint & ip_master,
				      const IntegrationPoint & ip_other,
				      FlatMatrix<Complex> & elmat,
				      LocalHeap & lh) const
    {
      FlatMatrix<double> rmat;
      CalcElementMatrixIndependent(bfel_master,bfel_other,
				       eltrans_master, eltrans_other,
				       ip_master, ip_other,
				       rmat, lh);
      elmat.AssignMemory(rmat.Height(), rmat.Width(), lh);
      elmat = rmat;
    }


    virtual void
    CalcElementMatrixDiag (const FiniteElement & fel,
			   const ElementTransformation & eltrans, 
			   FlatVector<double> diag,
			   LocalHeap & lh) const;




    virtual void
    CalcLinearizedElementMatrix (const FiniteElement & fel, 
				 const ElementTransformation & eltrans, 
				 FlatVector<double> elveclin,
				 FlatMatrix<double> elmat,
				 LocalHeap & lh) const;

    virtual void
    CalcLinearizedElementMatrix (const FiniteElement & fel, 
				 const ElementTransformation & eltrans, 
				 FlatVector<Complex> elveclin,
				 FlatMatrix<Complex> elmat,
				 LocalHeap & lh) const;


    virtual void *  
    PrecomputeData (const FiniteElement & fel, 
		    const ElementTransformation & eltrans, 
		    LocalHeap & lh) const { return 0; }
  

    virtual void 
    ApplyElementMatrix (const FiniteElement & fel, 
			const ElementTransformation & eltrans, 
			const FlatVector<double> elx, 
			FlatVector<double> ely,
			void * precomputed,
			LocalHeap & lh) const;

    virtual void 
    ApplyElementMatrix (const FiniteElement & fel, 
			const ElementTransformation & eltrans, 
			const FlatVector<Complex> elx, 
			FlatVector<Complex> ely,
			void * precomputed,
			LocalHeap & lh) const;

    virtual void 
    ApplyElementMatrixTrans (const FiniteElement & fel, 
                             const ElementTransformation & eltrans, 
                             const FlatVector<double> elx, 
                             FlatVector<double> ely,
                             void * precomputed,
                             LocalHeap & lh) const;
    
    virtual void 
    ApplyElementMatrixTrans (const FiniteElement & fel, 
                             const ElementTransformation & eltrans, 
                             const FlatVector<Complex> elx, 
                             FlatVector<Complex> ely,
                             void * precomputed,
                             LocalHeap & lh) const;
    


    
    /*
    template < int S, class T>
    void ApplyElementMatrix (const FiniteElement & fel, 
			     const ElementTransformation & eltrans, 
			     const FlatVector< Vec<S,T> > & elx, 
			     FlatVector< Vec<S,T> > & ely,
			     void * precomputed,
			     LocalHeap & lh) const
    {
      //cout << "call baseclass ApplyElementMatrix, type = " << typeid(*this).name() << endl;
      FlatMatrix<T> mat;
      CalcElementMatrix (fel, eltrans, mat, lh);
      ely = mat * elx;
    }
    */


    virtual void 
    ApplyLinearizedElementMatrix (const FiniteElement & fel, 
				  const ElementTransformation & eltrans, 
                                  FlatVector<double> ellin,
                                  FlatVector<double> elx, 
				  FlatVector<double> ely,
				  LocalHeap & lh) const;

    virtual void 
    ApplyLinearizedElementMatrix (const FiniteElement & fel, 
				  const ElementTransformation & eltrans, 
                                  FlatVector<Complex> ellin, 
                                  FlatVector<Complex> elx, 
				  FlatVector<Complex> ely,
				  LocalHeap & lh) const;



    virtual double Energy (const FiniteElement & fel, 
			   const ElementTransformation & eltrans, 
                           FlatVector<double> elx, 
			   LocalHeap & lh) const;

    virtual double Energy (const FiniteElement & fel, 
			   const ElementTransformation & eltrans, 
                           FlatVector<Complex> elx, 
			   LocalHeap & lh) const;




    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseMappedIntegrationPoint & bmip,
	      BareSliceVector<double> elx, 
	      FlatVector<double> flux,
	      bool applyd,
	      LocalHeap & lh) const;


    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseMappedIntegrationPoint & bmip,
              BareSliceVector<Complex> elx, 
	      FlatVector<Complex> flux,
	      bool applyd,
	      LocalHeap & lh) const;


    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseMappedIntegrationRule & mir,
              BareSliceVector<double> elx, 
	      BareSliceMatrix<double> flux,
	      bool applyd,
	      LocalHeap & lh) const;


    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseMappedIntegrationRule & mir,
              BareSliceVector<Complex> elx, 
	      BareSliceMatrix<Complex> flux,
	      bool applyd,
	      LocalHeap & lh) const;


    virtual void
    CalcFlux (const FiniteElement & fel,
	      const FiniteElement & felflux,
	      const ElementTransformation & eltrans,
              BareSliceVector<> elx, 
	      FlatVector<> flux,
	      bool applyd,
	      LocalHeap & lh) const;


    virtual void
    CalcFluxMulti (const FiniteElement & fel,
		   const BaseMappedIntegrationPoint & bmip,
		   int m,
		   FlatVector<double> elx, 
		   FlatVector<double> flux,
		   bool applyd,
		   LocalHeap & lh) const;


    virtual void
    CalcFluxMulti (const FiniteElement & fel,
		   const BaseMappedIntegrationPoint & bmip,
		   int m,
                   FlatVector<Complex> elx, 
		   FlatVector<Complex> flux,
		   bool applyd,
		   LocalHeap & lh) const;


    virtual void
    ApplyBTrans (const FiniteElement & fel,
		 const BaseMappedIntegrationPoint & bmip,
                 FlatVector<double> elx, 
		 FlatVector<double> ely,
		 LocalHeap & lh) const;

    virtual void
    ApplyBTrans (const FiniteElement & fel,
		 const BaseMappedIntegrationPoint & bmip,
                 FlatVector<Complex> elx, 
		 FlatVector<Complex> ely,
		 LocalHeap & lh) const;

    
    virtual void
    ApplyBTrans (const FiniteElement & fel,
		 const BaseMappedIntegrationRule & mir,
                 FlatMatrix<double> elx, 
		 FlatVector<double> ely,
		 LocalHeap & lh) const;

    virtual void
    ApplyBTrans (const FiniteElement & fel,
		 const BaseMappedIntegrationRule & mir,
                 FlatMatrix<Complex> elx, 
		 FlatVector<Complex> ely,
		 LocalHeap & lh) const;
    

    virtual void ApplyDMat (const FiniteElement & bfel,
			    const BaseMappedIntegrationPoint & bmip,
                            FlatVector<double> elx, 
			    FlatVector<double> eldx,
			    LocalHeap & lh) const;

    virtual void ApplyDMat (const FiniteElement & bfel,
			    const BaseMappedIntegrationPoint & bmip,
                            FlatVector<Complex> elx, 
			    FlatVector<Complex> eldx,
			    LocalHeap & lh) const;
  
    virtual void ApplyDMat (const FiniteElement & bfel,
			    const BaseMappedIntegrationRule & mir,
                            FlatMatrix<double> elx, 
			    FlatMatrix<double> eldx,
			    LocalHeap & lh) const;

    virtual void ApplyDMat (const FiniteElement & bfel,
			    const BaseMappedIntegrationRule & mir,
                            FlatMatrix<Complex> elx, 
			    FlatMatrix<Complex> eldx,
			    LocalHeap & lh) const;
  
    virtual void ApplyDMatInv (const FiniteElement & bfel,
			       const BaseMappedIntegrationPoint & bmip,
                               FlatVector<double> elx, 
			       FlatVector<double> eldx,
			       LocalHeap & lh) const;

    virtual void ApplyDMatInv (const FiniteElement & bfel,
			       const BaseMappedIntegrationPoint & bmip,
                               FlatVector<Complex> elx, 
			       FlatVector<Complex> eldx,
			       LocalHeap & lh) const;
  
    virtual void ApplyDMatInv (const FiniteElement & bfel,
			       const BaseMappedIntegrationRule & mir,
                               FlatMatrix<double> elx, 
			       FlatMatrix<double> eldx,
			       LocalHeap & lh) const;

    virtual void ApplyDMatInv (const FiniteElement & bfel,
			       const BaseMappedIntegrationRule & mir,
                               FlatMatrix<Complex> elx, 
			       FlatMatrix<Complex> eldx,
			       LocalHeap & lh) const;


    shared_ptr<DifferentialOperator> GetEvaluator(string name) const
    {
      return evaluators[name];
    }

    const auto & GetEvaluators() const { return evaluators; }
    /*
    virtual const IntegrationRule & GetIntegrationRule (const FiniteElement & fel,
							const bool use_higher_integration_order = false) const;
    */
    bool geom_free = false;
  };

  
/*
  class FacetNeighbourElInfo{
    public:
      //finite Element of the neighbour element      
      const FiniteElement & volumefel; 
      //local Facet Number from volumeElements view
      int LocalFacetNr;
      //Transformation of the neighbouring element
      const ElementTransformation & eltrans; 
      //Vertices of the Element
      FlatArray<int> & ElVertices;
      bool nonempty;
      FacetNeighbourElInfo(const FiniteElement & vvolumefel, int lLocalFacetNr, 
			   const ElementTransformation & eeltrans,
			   FlatArray<int> & eElVertices):volumefel(vvolumefel),
			   LocalFacetNr(lLocalFacetNr), eltrans(eeltrans), 
			   ElVertices(eElVertices), nonempty(true)
			   {;}
      FacetNeighbourElInfo():nonempty(false){;};
      bool IsEmpty(){return !nonempty;}
  };
*/

  class NGS_DLL_HEADER FacetBilinearFormIntegrator : public BilinearFormIntegrator
  {
    public:
      
    FacetBilinearFormIntegrator() // const Array<shared_ptr<CoefficientFunction>> & coeffs) 
      : BilinearFormIntegrator() { ; }

    ~FacetBilinearFormIntegrator() { ; }

      
      
      
    virtual VorB VB () const 
    { return VOL; }

    virtual bool SkeletonForm () const 
    { return 1; }
    
    virtual void CalcElementMatrix (const FiniteElement & fel,
				    const ElementTransformation & eltrans, 
				    FlatMatrix<double> elmat,
				    LocalHeap & lh) const {
      throw Exception ("FacetBilinearFormIntegrator can not assemble volumetric element matrices!");
    }
    

    virtual void
    CalcFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
			 const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
			 const FiniteElement & volumefel2, int LocalFacetNr2,
			 const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
			 FlatMatrix<double> elmat,
			 LocalHeap & lh) const{ 
      throw Exception ("FacetBilinearFormIntegrator::CalcFacetMatrix for inner facets not implemented!");
    }
    
    virtual void 
    CalcFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
			 const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
			 const FiniteElement & volumefel2, int LocalFacetNr2,
			 const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,	 
			 FlatMatrix<Complex> elmat,
			 LocalHeap & lh) const{ 
      throw Exception ("FacetBilinearFormIntegrator::CalcFacetMatrix<Complex> for inner facets not implemented!");
    }

    virtual void
    CalcLinearizedFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
			 const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
			 const FiniteElement & volumefel2, int LocalFacetNr2,
			 const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
			 FlatVector<double> elvec,
			 FlatMatrix<double> elmat,
			 LocalHeap & lh) const{
      CalcFacetMatrix (volumefel1, LocalFacetNr1, eltrans1, ElVertices1, volumefel2, LocalFacetNr2, eltrans2, ElVertices2, elmat, lh);
      // throw Exception ("FacetBilinearFormIntegrator::CalcLinearizedFacetMatrix for inner facets not implemented!");
    }
    
    virtual void 
    CalcLinearizedFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
			 const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
			 const FiniteElement & volumefel2, int LocalFacetNr2,
			 const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,	 
			 FlatVector<Complex> elvec,
			 FlatMatrix<Complex> elmat,
			 LocalHeap & lh) const{ 
      throw Exception ("FacetBilinearFormIntegrator::CalcLinearizedFacetMatrix<Complex> for inner facets not implemented!");
    }

    
    virtual void
      ApplyFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                        const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                        const FiniteElement & volumefel2, int LocalFacetNr2,
                        const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                        FlatVector<double> elx, FlatVector<double> ely,
                        LocalHeap & lh) const
    { 
      throw Exception ("FacetBilinearFormIntegrator::ApplyFacetMatrix for inner facets not implemented!");
    }
    virtual void
      ApplyFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                        const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                        const FiniteElement & volumefel2, int LocalFacetNr2,
                        const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                        FlatVector<Complex> elx, FlatVector<Complex> ely,
                        LocalHeap & lh) const
    { 
      throw Exception ("FacetBilinearFormIntegrator::ApplyFacetMatrix for inner facets not implemented!");
    }


    virtual void
    CalcFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
			 const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
			 const ElementTransformation & seltrans, FlatArray<int> & SElVertices,  
			 FlatMatrix<double> elmat,
			 LocalHeap & lh) const{ 
      throw Exception ("FacetBilinearFormIntegrator::CalcFacetMatrix for boundary facets not implemented!");
    }
    virtual void 
    CalcFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
			 const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
			 const ElementTransformation & seltrans, FlatArray<int> & SElVertices,  
			 FlatMatrix<Complex> elmat,
			 LocalHeap & lh) const{ 
      throw Exception ("FacetBilinearFormIntegrator::CalcFacetMatrix<Complex> for boundary facets not implemented!");
    }

    virtual void
      CalcLinearizedFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                                 const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                                 const ElementTransformation & seltrans, FlatArray<int> & SElVertices,  
                                 FlatVector<double> vec, FlatMatrix<double> elmat,
                                 LocalHeap & lh) const
    {
      CalcFacetMatrix (volumefel, LocalFacetNr,
                       eltrans, ElVertices, seltrans, SElVertices, elmat, lh);
    }
      
    virtual void
      CalcLinearizedFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                                 const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                                 const ElementTransformation & seltrans, FlatArray<int> & SElVertices,  
                                 FlatVector<Complex> vec, FlatMatrix<Complex> elmat,
                                 LocalHeap & lh) const
    {
      throw Exception ("CalcLinearizedFacetMatrix<Complex> not available");
    }
      
    
    virtual void
      ApplyFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                        const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                        const ElementTransformation & seltrans, FlatArray<int> & SElVertices,
                        FlatVector<double> elx, FlatVector<double> ely,
                        LocalHeap & lh) const
    { 
      throw Exception ("FacetBilinearFormIntegrator::ApplyFacetMatrix for boundary facets not implemented!");
    }
    virtual void
      ApplyFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                        const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                        const ElementTransformation & seltrans, FlatArray<int> & SElVertices,
                        FlatVector<Complex> elx, FlatVector<Complex> ely,
                        LocalHeap & lh) const
    { 
      throw Exception ("FacetBilinearFormIntegrator::ApplyFacetMatrix for boundary facets not implemented!");
    }


    // calculate traces in integration points
    virtual void
      CalcTraceValues (const FiniteElement & volumefel, int LocalFacetNr,
		       const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
		       FlatVector<double> & trace, FlatVector<double> elx, LocalHeap & lh) const
    { 
      throw Exception ("FacetBilinearFormIntegrator::ApplyFacetMatrix for boundary facets not implemented!");
    }

    virtual void
      CalcTraceValues (const FiniteElement & volumefel, int LocalFacetNr,
                       const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
		       FlatVector<Complex> & trace, FlatVector<Complex> elx, LocalHeap & lh) const
    { 
      throw Exception ("FacetBilinearFormIntegrator::ApplyFacetMatrix for boundary facets not implemented!");
    }
    
    virtual void
      ApplyFromTraceValues (const FiniteElement & volumefel, int LocalFacetNr,
			    const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
			    FlatVector<double> trace,
			    FlatVector<double> elx, FlatVector<double> ely, 
			    LocalHeap & lh) const
    { 
      throw Exception ("FacetBilinearFormIntegrator::ApplyFromTraceValues not implemented!");
    }
    
    virtual void
      ApplyFromTraceValues (const FiniteElement & volumefel, int LocalFacetNr,
			    const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
			    FlatVector<Complex> trace,
			    FlatVector<Complex> elx, FlatVector<Complex> ely, 
			    LocalHeap & lh) const
    { 
      throw Exception ("FacetBilinearFormIntegrator::ApplyFromTraceValues not implemented!");
    }


    
  };





  class NGS_DLL_HEADER BlockBilinearFormIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<BilinearFormIntegrator> bfi;
    int dim;
    int comp;
  public:
    BlockBilinearFormIntegrator (shared_ptr<BilinearFormIntegrator> abfi, int adim, int acomp);
    BlockBilinearFormIntegrator (shared_ptr<BilinearFormIntegrator> abfi, int adim);
    virtual ~BlockBilinearFormIntegrator ();

    virtual VorB VB () const override
    { return bfi->VB(); }
    virtual xbool IsSymmetric () const override { return bfi->IsSymmetric(); }
    virtual int DimFlux () const override 
    { return (comp == -1) ? dim * bfi->DimFlux() : bfi->DimFlux(); }
    int GetDim() const { return dim; }
    int GetComp() const { return comp; } 

    const BilinearFormIntegrator & Block () const { return *bfi; }
    shared_ptr<BilinearFormIntegrator> BlockPtr () const { return bfi; }

    virtual void
    CalcElementMatrix (const FiniteElement & bfel, 
		       const ElementTransformation & eltrans, 
		       FlatMatrix<double> elmat,
		       LocalHeap & lh) const override;

    virtual void
    CalcElementMatrix (const FiniteElement & bfel, 
		       const ElementTransformation & eltrans, 
		       FlatMatrix<Complex> elmat,
		       LocalHeap & lh) const override;
    
    virtual void 
    ApplyElementMatrix (const FiniteElement & bfel, 
			const ElementTransformation & eltrans, 
			const FlatVector<double> elx, 
			FlatVector<double> ely,
			void * precomputed,
			LocalHeap & lh) const override;

    virtual void 
    ApplyElementMatrix (const FiniteElement & bfel, 
			const ElementTransformation & eltrans, 
			const FlatVector<Complex> elx, 
			FlatVector<Complex> ely,
			void * precomputed,
			LocalHeap & lh) const override;

    virtual void
    CalcLinearizedElementMatrix (const FiniteElement & bfel,
				 const ElementTransformation & eltrans,
				 FlatVector<double> elveclin,
				 FlatMatrix<double> elmat,
				 LocalHeap & lh) const override;
    virtual void
    CalcLinearizedElementMatrix (const FiniteElement & bfel, 
				 const ElementTransformation & eltrans, 
				 FlatVector<Complex> elveclin,
				 FlatMatrix<Complex> elmat,
				 LocalHeap & lh) const override;

    /*
    virtual void
    CalcFlux (const FiniteElement & fel,
	      const ElementTransformation & eltrans,
	      const IntegrationPoint & ip,
              FlatVector<double> elx, 
	      FlatVector<double> flux,
	      bool applyd,
	      LocalHeap & lh) const;

    virtual void
    CalcFlux (const FiniteElement & fel,
	      const ElementTransformation & eltrans,
	      const IntegrationPoint & ip,
              FlatVector<Complex> elx, 
	      FlatVector<Complex> flux,
	      bool applyd,
	      LocalHeap & lh) const;
    */


    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseMappedIntegrationPoint & bmip,
              BareSliceVector<double> elx, 
	      FlatVector<double> flux,
	      bool applyd,
	      LocalHeap & lh) const override;

    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseMappedIntegrationPoint & bmip,
              BareSliceVector<Complex> elx, 
	      FlatVector<Complex> flux,
	      bool applyd,
	      LocalHeap & lh) const override;

    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseMappedIntegrationRule & mir,
              BareSliceVector<double> elx, 
	      BareSliceMatrix<double> flux,
	      bool applyd,
	      LocalHeap & lh) const override;



    virtual void
    ApplyBTrans (const FiniteElement & bfel,
		 const BaseMappedIntegrationPoint & bmip,
                 FlatVector<double> elx, 
		 FlatVector<double> ely,
		 LocalHeap & lh) const override;

    virtual void
    ApplyBTrans (const FiniteElement & bfel,
		 const BaseMappedIntegrationPoint & bmip,
                 FlatVector<Complex> elx, 
		 FlatVector<Complex> ely,
		 LocalHeap & lh) const override;

    virtual void
    ApplyBTrans (const FiniteElement & fel,
		 const BaseMappedIntegrationRule & mir,
                 FlatMatrix<double> elx, 
		 FlatVector<double> ely,
		 LocalHeap & lh) const override;

    virtual double Energy (const FiniteElement & fel, 
			   const ElementTransformation & eltrans, 
                           FlatVector<double> elx, 
			   LocalHeap & lh) const override;

    virtual string Name () const override;
  };







  class NGS_DLL_HEADER ComplexBilinearFormIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<BilinearFormIntegrator> bfi;
    Complex factor;
  public:
    ComplexBilinearFormIntegrator (shared_ptr<BilinearFormIntegrator> abfi,
				   Complex afactor);

    virtual VorB VB () const override
    { return bfi->VB(); }

    virtual int DimFlux () const override
    { return bfi->DimFlux(); }
    virtual int DimElement () const override
    { return bfi->DimElement(); }
    virtual int DimSpace () const override
    { return bfi->DimSpace(); }
    virtual xbool IsSymmetric () const override
    { return bfi->IsSymmetric(); }


    virtual void GetFactor(Complex & fac) const {fac = factor;}
    virtual void GetFactor(double & fac) const {fac = factor.real();}
  
    virtual shared_ptr<BilinearFormIntegrator> GetBFI(void) const {return bfi;}

    virtual void CheckElement (const FiniteElement & el) const override { bfi->CheckElement(el); }


    virtual void
    CalcElementMatrix (const FiniteElement & fel, 
		       const ElementTransformation & eltrans, 
		       FlatMatrix<double> elmat,
		       LocalHeap & lh) const override;

    virtual void
    CalcElementMatrix (const FiniteElement & fel, 
		       const ElementTransformation & eltrans, 
		       FlatMatrix<Complex> elmat,
		       LocalHeap & lh) const override;
    
    virtual void
    CalcElementMatrixIndependent (const FiniteElement & bfel_master,
				      const FiniteElement & bfel_master_element,				    
				      const FiniteElement & bfel_other,
				      const ElementTransformation & eltrans_master, 
				      const ElementTransformation & eltrans_master_element, 
				      const ElementTransformation & eltrans_other,
				      const IntegrationPoint & ip_master,
				      const IntegrationPoint & ip_master_element,
				      const IntegrationPoint & ip_other,
				      FlatMatrix<double> & elmat,
				      LocalHeap & lh) const override;

    virtual void
    CalcElementMatrixIndependent (const FiniteElement & bfel_master,
				      const FiniteElement & bfel_master_element,				    
				      const FiniteElement & bfel_other,
				      const ElementTransformation & eltrans_master, 
				      const ElementTransformation & eltrans_master_element, 
				      const ElementTransformation & eltrans_other,
				      const IntegrationPoint & ip_master,
				      const IntegrationPoint & ip_master_element,
				      const IntegrationPoint & ip_other,
				      FlatMatrix<Complex> & elmat,
				      LocalHeap & lh) const override;

    virtual void
    CalcElementMatrixIndependent (const FiniteElement & bfel_master,
				      const FiniteElement & bfel_other,
				      const ElementTransformation & eltrans_master, 
				      const ElementTransformation & eltrans_other,
				      const IntegrationPoint & ip_master,
				      const IntegrationPoint & ip_other,
				      FlatMatrix<double> & elmat,
				      LocalHeap & lh) const override;

    virtual void
    CalcElementMatrixIndependent (const FiniteElement & bfel_master,
				      const FiniteElement & bfel_other,
				      const ElementTransformation & eltrans_master, 
				      const ElementTransformation & eltrans_other,
				      const IntegrationPoint & ip_master,
				      const IntegrationPoint & ip_other,
				      FlatMatrix<Complex> & elmat,
				      LocalHeap & lh) const override;

  

    virtual void 
    ApplyElementMatrix (const FiniteElement & fel, 
			const ElementTransformation & eltrans, 
			const FlatVector<Complex> elx, 
			FlatVector<Complex> ely,
			void * precomputed,
			LocalHeap & lh) const override;
    /*
    virtual void
    CalcFlux (const FiniteElement & fel,
	      const ElementTransformation & eltrans,
	      const IntegrationPoint & ip,
              FlatVector<Complex> elx, 
	      FlatVector<Complex> flux,
	      bool applyd,
	      LocalHeap & lh) const;
    */

    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseMappedIntegrationPoint & bmip,
              BareSliceVector<Complex> elx, 
	      FlatVector<Complex> flux,
	      bool applyd,
	      LocalHeap & lh) const override;

    virtual string Name () const override;

    /*
    virtual const IntegrationRule & GetIntegrationRule (const FiniteElement & fel,
							const bool use_higher_integration_order = false) const;
    */
  };





  class NGS_DLL_HEADER TransposeBilinearFormIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<BilinearFormIntegrator> bfi;
  public:
    TransposeBilinearFormIntegrator (shared_ptr<BilinearFormIntegrator> abfi)
      : bfi(abfi) { ; }
    virtual ~TransposeBilinearFormIntegrator () { ; }
  
    shared_ptr<BilinearFormIntegrator> GetBFI(void) const {return bfi;}

    virtual VorB VB () const
    { return bfi->VB(); }

    virtual int DimFlux () const 
    { return bfi->DimFlux(); }
    virtual int DimElement () const
    { return bfi->DimElement(); }
    virtual int DimSpace () const
    { return bfi->DimSpace(); }
    virtual xbool IsSymmetric () const
    { return bfi->IsSymmetric(); }

    virtual void CheckElement (const FiniteElement & el) const
    {
      return bfi->CheckElement (el);
    }

    virtual void
    CalcElementMatrix (const FiniteElement & bfel, 
		       const ElementTransformation & eltrans, 
		       FlatMatrix<double> elmat,
		       LocalHeap & lh) const;
  };



  class NGS_DLL_HEADER BilinearFormIntegratorAnyDim : public BilinearFormIntegrator
  {
    shared_ptr<BilinearFormIntegrator> bfi[4]; // dim 0 ... dim 3
    shared_ptr<BilinearFormIntegrator> any_dim;
  public:
    BilinearFormIntegratorAnyDim (shared_ptr<BilinearFormIntegrator> abfi[4])
    { 
      for (int i = 0; i < 4; i++)
        {
          bfi[i] = abfi[i];
          if (bfi[i]) any_dim = bfi[i];
        }
    }
  
    shared_ptr<BilinearFormIntegrator> GetBFI(int dim) const 
    { 
      if (!bfi[dim]) 
        throw Exception (string("BFI for dimension") + ToString(dim)+"not available");
      bfi[dim]->SetDefinedOn(definedon);
      return bfi[dim];
    }

    virtual VorB VB () const
    { return any_dim->VB(); }

    virtual int DimFlux () const 
    { throw Exception("BFI AnyDim - DimFlux not available"); }
    virtual int DimElement () const
    { throw Exception("BFI AnyDim - DimElement not available"); }
    virtual int DimSpace () const
    { throw Exception("BFI AnyDim - DimSpace not available"); }
    virtual xbool IsSymmetric () const
    { return any_dim->IsSymmetric(); }

    virtual void CheckElement (const FiniteElement & el) const;

    virtual void
    CalcElementMatrix (const FiniteElement & bfel, 
		       const ElementTransformation & eltrans, 
		       FlatMatrix<double> elmat,
		       LocalHeap & lh) const;
    virtual void
    CalcElementMatrix (const FiniteElement & bfel, 
		       const ElementTransformation & eltrans, 
		       FlatMatrix<Complex> elmat,
		       LocalHeap & lh) const;
  };





  class NGS_DLL_HEADER CompoundBilinearFormIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<BilinearFormIntegrator> bfi;
    int comp;
  public:
    CompoundBilinearFormIntegrator (shared_ptr<BilinearFormIntegrator> abfi, int acomp);
  
    shared_ptr<BilinearFormIntegrator> GetBFI(void) const {return bfi;}
    int GetComponent() const {return comp;}
    virtual VorB VB () const override
    { return bfi->VB(); }

    virtual int DimFlux () const override
    { return bfi->DimFlux(); }
    virtual int DimElement () const override
    { return bfi->DimElement(); }
    virtual int DimSpace () const override
    { return bfi->DimSpace(); }
    virtual xbool IsSymmetric () const override
    { return bfi->IsSymmetric(); }
    virtual bool SkeletonForm () const override
    { return bfi->SkeletonForm(); }
    virtual void CheckElement (const FiniteElement & el) const override
    {
      return bfi->CheckElement (dynamic_cast<const CompoundFiniteElement&>(el)[comp]);
    }

    virtual void
    CalcElementMatrix (const FiniteElement & bfel, 
		       const ElementTransformation & eltrans, 
		       FlatMatrix<double> elmat,
		       LocalHeap & lh) const override;

    virtual void
    CalcElementMatrix (const FiniteElement & bfel, 
		       const ElementTransformation & eltrans, 
		       FlatMatrix<Complex> elmat,
		       LocalHeap & lh) const override;


    virtual void
    CalcLinearizedElementMatrix (const FiniteElement & fel, 
				 const ElementTransformation & eltrans, 
				 FlatVector<double> elveclin,
				 FlatMatrix<double> elmat,
				 LocalHeap & lh) const override;
    
    virtual void
    CalcLinearizedElementMatrix (const FiniteElement & fel, 
                                 const ElementTransformation & eltrans, 
                                 FlatVector<Complex> elveclin,
                                 FlatMatrix<Complex> elmat,
                                 LocalHeap & lh) const override;

    virtual void
    ApplyElementMatrix (const FiniteElement & bfel, 
			const ElementTransformation & eltrans, 
                        FlatVector<double> elx,
			FlatVector<double> ely,
			void * precomputed,
			LocalHeap & lh) const override;

    virtual void
    ApplyElementMatrix (const FiniteElement & bfel, 
			const ElementTransformation & eltrans, 
                        FlatVector<Complex> elx,
			FlatVector<Complex> ely,
			void * precomputed,
			LocalHeap & lh) const override;

    virtual void
    ApplyLinearizedElementMatrix (const FiniteElement & bfel, 
				  const ElementTransformation & eltrans, 
                                  FlatVector<double> ellin,
                                  FlatVector<double> elx,
				  FlatVector<double> ely,
				  LocalHeap & lh) const override;

    virtual void
    ApplyLinearizedElementMatrix (const FiniteElement & bfel, 
				  const ElementTransformation & eltrans, 
                                  FlatVector<Complex> ellin,
                                  FlatVector<Complex> elx,
				  FlatVector<Complex> ely,
				  LocalHeap & lh) const override;


    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseMappedIntegrationPoint & bmip,
              BareSliceVector<double> elx, 
	      FlatVector<double> flux,
	      bool applyd,
	      LocalHeap & lh) const override;


    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseMappedIntegrationPoint & bmip,
              BareSliceVector<Complex> elx, 
	      FlatVector<Complex> flux,
	      bool applyd,
	      LocalHeap & lh) const override;


    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseMappedIntegrationRule & mir,
              BareSliceVector<double> elx, 
	      BareSliceMatrix<double> flux,
	      bool applyd,
	      LocalHeap & lh) const override;


    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseMappedIntegrationRule & mir,
              BareSliceVector<Complex> elx, 
	      BareSliceMatrix<Complex> flux,
	      bool applyd,
	      LocalHeap & lh) const override;


    virtual void
    ApplyBTrans (const FiniteElement & fel,
		 const BaseMappedIntegrationPoint & bmip,
                 FlatVector<double> elx, 
		 FlatVector<double> ely,
		 LocalHeap & lh) const override;

    virtual void
    ApplyBTrans (const FiniteElement & fel,
		 const BaseMappedIntegrationPoint & bmip,
                 FlatVector<Complex> elx, 
		 FlatVector<Complex> ely,
		 LocalHeap & lh) const override;

    virtual string Name () const override;
  };









  /**
     Integrator for element vector.
  */
  class NGS_DLL_HEADER LinearFormIntegrator : public Integrator

  {
  public:
    ///
    LinearFormIntegrator () { ; }
    ///
    virtual ~LinearFormIntegrator () { ; }


    /**
       Computes the element vector.
    */
    virtual void 
    CalcElementVector (const FiniteElement & fel,
		       const ElementTransformation & eltrans, 
		       FlatVector<double> elvec,
		       LocalHeap & lh) const;

    virtual void 
    CalcElementVector (const FiniteElement & fel,
		       const ElementTransformation & eltrans, 
		       FlatVector<Complex> elvec,
		       LocalHeap & lh) const;

    
    virtual void
    CalcElementVectorIndependent (const FiniteElement & gfel,
				      const BaseMappedIntegrationPoint & s_mip,
				      const BaseMappedIntegrationPoint & g_mip,
				      FlatVector<double> & elvec,
				      LocalHeap & lh,
				      const bool curveint = false) const
    {
      cerr << "CalcElementVectorIndependent called for base-class!" << endl;
      exit(10);
    }
  
    virtual void
    CalcElementVectorIndependent (const FiniteElement & gfel,
				      const BaseMappedIntegrationPoint & s_mip,
				      const BaseMappedIntegrationPoint & g_mip,
				      FlatVector<Complex> & elvec,
				      LocalHeap & lh,
				      const bool curveint = false) const
    {
      FlatVector<double> rvec(elvec.Size(), lh);
      CalcElementVectorIndependent (gfel, s_mip, g_mip, rvec, lh,curveint);
      elvec = rvec;
    }

  
  };















 class NGS_DLL_HEADER FacetLinearFormIntegrator : public LinearFormIntegrator
  {
    public:
      
    FacetLinearFormIntegrator( /* const Array<shared_ptr<CoefficientFunction>> & coeffs */) 
      : LinearFormIntegrator() { ; }

    ~FacetLinearFormIntegrator() { ; }

    virtual VorB VB () const 
    { return BND; }

    virtual bool SkeletonForm () const 
    { return 1; }
    
    virtual void 
    CalcElementVector (const FiniteElement & bfel, 
			   const ElementTransformation & eltrans, 
			   FlatVector<double> elvec,
			   LocalHeap & lh) const{
      throw Exception ("FacetLinearFormIntegrator can not assemble volumetric element matrices!");
    }
    
    virtual void
    CalcFacetVector (const FiniteElement & volumefel1, int LocalFacetNr1,
			 const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
			 const FiniteElement & volumefel2, int LocalFacetNr2,
			 const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
			 FlatVector<double> elvec,
			 LocalHeap & lh) const{ 
      throw Exception ("FacetLinearFormIntegrator::CalcFacetVector for inner facets not implemented!");
    }
    
    virtual void 
    CalcFacetVector (const FiniteElement & volumefel1, int LocalFacetNr1,
			 const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
			 const FiniteElement & volumefel2, int LocalFacetNr2,
			 const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,	 
			 FlatVector<Complex> elvec,
			 LocalHeap & lh) const{ 
      throw Exception ("FacetLinearFormIntegrator::CalcFacetVector<Complex> for inner facets not implemented!");
    }

    virtual void
    CalcFacetVector (const FiniteElement & volumefel, int LocalFacetNr,
			 const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
			 const ElementTransformation & seltrans,
			 FlatVector<double> elvec,
			 LocalHeap & lh) const{ 
      throw Exception ("FacetLinearFormIntegrator::CalcFacetVector not implemented!");
    }

    virtual void 
    CalcFacetVector (const FiniteElement & volumefel, int LocalFacetNr,
			 const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
			 const ElementTransformation & seltrans,
			 FlatVector<Complex> elvec,
			 LocalHeap & lh) const { 
      throw Exception ("FacetLinearFormIntegrator::CalcFacetVector<complex> not implemented!");
    }
			 
  };  





  class NGS_DLL_HEADER BlockLinearFormIntegrator : public LinearFormIntegrator
  {
    shared_ptr<LinearFormIntegrator> lfi;
    int dim;
    int comp;
  public:
    BlockLinearFormIntegrator (shared_ptr<LinearFormIntegrator> alfi, int adim, int acomp);

    virtual VorB VB () const
    { return lfi->VB(); }


    virtual void 
    CalcElementVector (const FiniteElement & bfel, 
		       const ElementTransformation & eltrans, 
		       FlatVector<double> elvec,
		       LocalHeap & lh) const;
  };





  class NGS_DLL_HEADER ComplexLinearFormIntegrator : public LinearFormIntegrator
  {
    shared_ptr<LinearFormIntegrator> lfi;
    Complex factor;
  public:
    ComplexLinearFormIntegrator (shared_ptr<LinearFormIntegrator> alfi, 
				 Complex afactor);
    virtual ~ComplexLinearFormIntegrator();
    
    virtual VorB VB () const override;
    virtual void CheckElement (const FiniteElement & el) const override;


    virtual void
    CalcElementVector (const FiniteElement & fel, 
		       const ElementTransformation & eltrans, 
		       FlatVector<double> elvec,
		       LocalHeap & lh) const override;

    virtual void
    CalcElementVector (const FiniteElement & fel, 
		       const ElementTransformation & eltrans, 
		       FlatVector<Complex> elvec,
		       LocalHeap & lh) const override;

    virtual void
    CalcElementVectorIndependent (const FiniteElement & gfel, 
                                  const BaseMappedIntegrationPoint & s_mip,
                                  const BaseMappedIntegrationPoint & g_mip,
                                  FlatVector<double> & elvec,
                                  LocalHeap & lh,
                                  const bool curveint = false) const override;
  

    virtual void
    CalcElementVectorIndependent (const FiniteElement & gfel, 
				      const BaseMappedIntegrationPoint & s_mip,
				      const BaseMappedIntegrationPoint & g_mip,
				      FlatVector<Complex> & elvec,
				      LocalHeap & lh,
                                  const bool curveint = false) const override;

    virtual string Name () const override;
  };



  
  class NGS_DLL_HEADER CompoundLinearFormIntegrator : public LinearFormIntegrator
  {
    shared_ptr<LinearFormIntegrator> lfi;
    int comp;
  public:
    CompoundLinearFormIntegrator (shared_ptr<LinearFormIntegrator> alfi, int acomp)
      : lfi(alfi), comp(acomp)
    {
        is_curve_integrator = lfi->IntegrationAlongCurve();
    }

    VorB VB () const override
    { return lfi->VB(); }

    void
    CalcElementVector (const FiniteElement & bfel, 
		       const ElementTransformation & eltrans, 
		       FlatVector<double> elvec,
		       LocalHeap & lh) const override;

    void
    CalcElementVector (const FiniteElement & bfel, 
		       const ElementTransformation & eltrans, 
		       FlatVector<Complex> elvec,
		       LocalHeap & lh) const override;

    void
    CalcElementVectorIndependent (const FiniteElement & gfel,
				      const BaseMappedIntegrationPoint & s_mip,
				      const BaseMappedIntegrationPoint & g_mip,
				      FlatVector<double> & elvec,
				      LocalHeap & lh,
				      const bool curveint = false) const override;

    void
    CalcElementVectorIndependent (const FiniteElement & gfel,
				      const BaseMappedIntegrationPoint & s_mip,
				      const BaseMappedIntegrationPoint & g_mip,
				      FlatVector<Complex> & elvec,
				      LocalHeap & lh,
				      const bool curveint = false) const override;

    string Name () const override
    {
      return string ("CompoundIntegrator (") + lfi->Name() + ")";
    }

    int NumCurvePoints() const override { return lfi->NumCurvePoints(); }
    FlatVector<double> CurvePoint(int i) override
    { return lfi->CurvePoint(i); }
    FlatVector<double> CurvePointTangent(int i) override
    { return lfi->CurvePointTangent(i); }
    int GetNumCurveParts() const override
    { return lfi->GetNumCurveParts(); }
    int GetStartOfCurve(int i) const override
    { return lfi->GetStartOfCurve(i); }
    int GetEndOfCurve(int i) const override
    { return lfi->GetEndOfCurve(i); }

    void AppendCurvePoint(const FlatVector<double> & point) override
    { lfi->AppendCurvePoint(point); }
    void AppendCurvePoint(const FlatVector<double> & point,
                          const FlatVector<double> & tangent) override
    { lfi->AppendCurvePoint(point, tangent); }
    void SetCurveClearance() override
    { lfi->SetCurveClearance(); }
    void SetCacheComp(const int comp) override
    { cachecomp = comp; }
    int CacheComp() const override
    { return lfi->CacheComp(); }

  };





  class NGS_DLL_HEADER LinearFormIntegratorAnyDim : public LinearFormIntegrator
  {
    shared_ptr<LinearFormIntegrator> lfi[4]; // dim 0 ... dim 3
    shared_ptr<LinearFormIntegrator> any_dim;
  public:
    LinearFormIntegratorAnyDim (shared_ptr<LinearFormIntegrator> alfi[4])
    { 
      for (int i = 0; i < 4; i++)
        {
          lfi[i] = alfi[i];
          if (lfi[i]) any_dim = lfi[i];
        }
    }
  
    shared_ptr<LinearFormIntegrator> GetLFI(int dim) const 
    { 
      if (!lfi[dim]) 
        throw Exception (string("LFI for dimension") + ToString(dim)+"not available");
      lfi[dim]->SetDefinedOn(definedon);
      return lfi[dim];
    }

    virtual VorB VB () const
    { return any_dim->VB(); }
    virtual int DimElement () const
    { throw Exception("BFI AnyDim - DimElement not available"); }
    virtual int DimSpace () const
    { throw Exception("BFI AnyDim - DimSpace not available"); }

    virtual void CheckElement (const FiniteElement & el) const;

    virtual void
    CalcElementVector (const FiniteElement & bfel, 
		       const ElementTransformation & eltrans, 
		       FlatVector<double> elvec,
		       LocalHeap & lh) const;

    virtual void
    CalcElementVector (const FiniteElement & bfel, 
		       const ElementTransformation & eltrans, 
		       FlatVector<Complex> elvec,
		       LocalHeap & lh) const;
  };


  class CalcFluxDifferentialOperator : public DifferentialOperator
  {
    shared_ptr<BilinearFormIntegrator> bfi;
    bool applyd;
  public:
    CalcFluxDifferentialOperator (shared_ptr<BilinearFormIntegrator> _bfi, bool _applyd)
      : DifferentialOperator(_bfi->DimFlux(), 1, _bfi->VB(), 0), bfi(_bfi)
    { ; }

    virtual shared_ptr<DifferentialOperator> GetTrace() const override { return nullptr; }

    
    virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationPoint & mip,
		BareSliceMatrix<double,ColMajor> mat,   
		LocalHeap & lh) const override
    {
      throw Exception ("CalcFluxDifferentialOperator::CalcMatrix not available");
    }

    virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationPoint & bmip,
		BareSliceMatrix<Complex,ColMajor> mat, 
		LocalHeap & lh) const override
    {
      throw Exception ("CalcFluxDifferentialOperator::CalcMatrix not available");
    }
      
    virtual void
    CalcMatrix (const FiniteElement & fel,
                const BaseMappedIntegrationRule & mir,
		BareSliceMatrix<double,ColMajor> mat,   
		LocalHeap & lh) const override
    {
      throw Exception ("CalcFluxDifferentialOperator::CalcMatrix not available");
    }

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationRule & mir,
		BareSliceMatrix<Complex,ColMajor> mat,   
		LocalHeap & lh) const override
    {
      throw Exception ("CalcFluxDifferentialOperator::CalcMatrix not available");
    }
      
    
    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const SIMD_BaseMappedIntegrationRule & mir,
		BareSliceMatrix<SIMD<double>> mat) const override
    {
      throw Exception ("CalcFluxDifferentialOperator::CalcMatrix not available");
    }

    using DifferentialOperator::Apply;
    virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationPoint & mip,
	   BareSliceVector<double> x, 
	   FlatVector<double> flux,
	   LocalHeap & lh) const override
    {
      bfi->CalcFlux(fel, mip, x, flux, applyd, lh);
    }

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationPoint & mip,
	   BareSliceVector<Complex> x, 
	   FlatVector<Complex> flux,
	   LocalHeap & lh) const override
    {
      bfi->CalcFlux(fel, mip, x, flux, applyd, lh);
    }

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationRule & mir,
	   BareSliceVector<double> x, 
	   BareSliceMatrix<double> flux,
	   LocalHeap & lh) const override
    {
      bfi->CalcFlux(fel, mir, x, flux, applyd, lh);
    }

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationRule & mir,
	   BareSliceVector<Complex> x, 
	   BareSliceMatrix<Complex> flux,
	   LocalHeap & lh) const override
    {
      bfi->CalcFlux(fel, mir, x, flux, applyd, lh);
    }


    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
	   const SIMD_BaseMappedIntegrationRule & mir,
	   BareSliceVector<double> x, 
	   BareSliceMatrix<SIMD<double>> flux) const override
    {
      // bfi->CalcFlux(fel, mir, x, flux, applyd);
      throw ExceptionNOSIMD (string("CalcFluxDiffop: simd is not supported"));
    }

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
	   const SIMD_BaseMappedIntegrationRule & mir,
	   BareSliceVector<Complex> x, 
	   BareSliceMatrix<SIMD<Complex>> flux) const override
    {
      // bfi->CalcFlux(fel, mir, x, flux, applyd);
      throw ExceptionNOSIMD (string("CalcFluxDiffop: simd is not supported"));
    }
    
  };


  /// container for all integrators
  class NGS_DLL_HEADER Integrators
  {
  public:

    /// description of integrator
    template<typename T>
    class IntegratorInfo
    {
    public:
      string name;
      int spacedim;
      int numcoeffs;
      shared_ptr<T> (*creator)(const Array<shared_ptr<CoefficientFunction>> &);
    
      IntegratorInfo (const string & aname,
		      int aspacedim,
		      int anumcoffs,		      
		      shared_ptr<T> (*acreator)(const Array<shared_ptr<CoefficientFunction>> &));
    };
  

    Array<IntegratorInfo<BilinearFormIntegrator>*> bfis;
    Array<IntegratorInfo<LinearFormIntegrator>*> lfis;
  
  public:
    ///
    Integrators();
    ///
    ~Integrators();  
    ///
    void AddBFIntegrator (const string & aname, int aspacedim, int anumcoeffs,
			  shared_ptr<BilinearFormIntegrator> (*acreator)(const Array<shared_ptr<CoefficientFunction>> &));
    ///
    void AddLFIntegrator (const string & aname, int aspacedim, int anumcoeffs,
			  shared_ptr<LinearFormIntegrator> (*acreator)(const Array<shared_ptr<CoefficientFunction>> &));
    
    ///
    const Array<IntegratorInfo<BilinearFormIntegrator>*> & GetBFIs() const { return bfis; }
    ///
    const IntegratorInfo<BilinearFormIntegrator> * GetBFI(const string & name, int dim) const;
    ///
    shared_ptr<BilinearFormIntegrator> CreateBFI(const string & name, int dim,
                                                 const Array<shared_ptr<CoefficientFunction>> & coeffs) const;
    ///
    shared_ptr<BilinearFormIntegrator> CreateBFI(const string & name, int dim,
                                                 shared_ptr<CoefficientFunction> coef) const;
    shared_ptr<BilinearFormIntegrator> CreateBFI(const string & name, int dim,
                                                 const CoefficientFunction* coef) const;
    
    ///
    const Array<IntegratorInfo<LinearFormIntegrator>*> & GetLFIs() const { return lfis; }
    ///
    const IntegratorInfo<LinearFormIntegrator> * GetLFI(const string & name, int dim) const;
    ///
    shared_ptr<LinearFormIntegrator> CreateLFI(const string & name, int dim,
                                               const Array<shared_ptr<CoefficientFunction>> & coeffs) const;

    shared_ptr<LinearFormIntegrator> CreateLFI(const string & name, int dim,
                                               shared_ptr<CoefficientFunction> coef) const;

    ///
    void Print (ostream & ost) const;
  };

  /// 
  extern NGS_DLL_HEADER Integrators & GetIntegrators ();

  template <typename ... ARGS>
  inline shared_ptr<BilinearFormIntegrator> CreateBFI (ARGS ... args)
  {
    return GetIntegrators().CreateBFI (args...);
  }

  template <typename ... ARGS>
  inline shared_ptr<LinearFormIntegrator> CreateLFI (ARGS ... args)
  {
    return GetIntegrators().CreateLFI (args...);
  }

  /*
  class ConvertCoefs
  {
    Array<shared_ptr<CoefficientFunction> > coefs;
    Array<CoefficientFunction*> pcoefs;
  public:

    ConvertCoefs (const Array<shared_ptr<CoefficientFunction>> & acoefs)
      : coefs (acoefs)  // explicit copy !
    { 
      for (int i = 0; i < acoefs.Size(); i++)
        pcoefs.Append (acoefs[i].get());
    }

    ConvertCoefs (const Array<CoefficientFunction*> & acoefs)
      : pcoefs(acoefs)
    {
      for (int i = 0; i < acoefs.Size(); i++)
        coefs.Append (shared_ptr<CoefficientFunction> (acoefs[i], NOOP_Deleter));
    }

    operator Array<shared_ptr<CoefficientFunction>> () const 
    {
      return Array<shared_ptr<CoefficientFunction>> (coefs); 
    }

    operator Array<CoefficientFunction*> () const 
    {
      return Array<CoefficientFunction*> (pcoefs);
    }
  };
  */

  template <typename BFI>
  class RegisterBilinearFormIntegrator
  {
  public:
    RegisterBilinearFormIntegrator (string label, int dim, int numcoeffs)
    {
      GetIntegrators().AddBFIntegrator (label, dim, numcoeffs, Create);
      // cout << "register bf-integrator '" << label << "'" << endl;
    }
    
    static shared_ptr<BilinearFormIntegrator> Create (const Array<shared_ptr<CoefficientFunction>> & coefs)
    {
      // return shared_ptr<BilinearFormIntegrator> (new BFI (ConvertCoefs (coefs)));
      // return make_shared<BFI>(ConvertCoefs (coefs));
      return make_shared<BFI>(coefs);
    }

    // static shared_ptr<Integrator> Create (Array<CoefficientFunction*> & coefs)
    // {
    // return new BFI (ConvertCoefs (coefs));
    // }
  };



  template <typename LFI>
  class RegisterLinearFormIntegrator
  {
  public:
    RegisterLinearFormIntegrator (string label, int dim, int numcoeffs)
    {
      GetIntegrators().AddLFIntegrator (label, dim, numcoeffs, Create);
      // cout << "register lf-integrator '" << label << "'" << endl;
    }
    
    static shared_ptr<LinearFormIntegrator> Create (const Array<shared_ptr<CoefficientFunction>> & coefs)
    {
      // return shared_ptr<LinearFormIntegrator> (new LFI (ConvertCoefs(coefs)));
      return make_shared<LFI> (coefs);
      // return new LFI (coefs);
    }
  };




}
#endif
