#ifndef FILE_INTEGRATOR
#define FILE_INTEGRATOR

/*********************************************************************/
/* File:   integrator.hpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

namespace ngfem
{

  /*
    bilinear-form and linear-form integrators
  */


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
    Array < FlatVector < double > * > curve_ips;
    Array < FlatVector < double > * > curve_ip_tangents;
    Array <int> continuous_curveparts;
  
    int cachecomp;

  
  protected:
    void DeleteCurveIPs ( void );

  public:
    /// constructor
    Integrator() throw ();

    /// destructor
    virtual ~Integrator();

    /// integrates on the boundary, or on the domain ?
    virtual bool BoundaryForm () const = 0;

    /// integrates just on the skeleton, standard is NO
    virtual bool SkeletonForm () const {return 0;} 

    /// Is Integrator defined on this sub-domain ?
    bool DefinedOn (int mat) const;

    /// defined only on some subdomains
    void SetDefinedOn (const BitArray & adefinedon);

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

    /// benefit from constant coefficient
    void SetConstantCoefficient (bool acc = 1)
    { const_coef = acc; }

    /// dimension of element
    virtual int DimElement () const { return -1; }

    /// dimension of space
    virtual int DimSpace () const { return -1; }

    /// 
    void SetName (const string & aname);
    ///
    virtual string Name () const;

    /// does element match integrator ?
    virtual void CheckElement (const FiniteElement & el) const { ; }



    // special hacks by Markus
    bool IntegrationAlongCurve (void) const
    { return curve_ips.Size() > 0; }

    void SetIntegrationAlongCurve ( const int npoints );

    void UnSetIntegrationAlongCurve ( void );

    int NumCurvePoints(void) const
    { return curve_ips.Size(); }

    FlatVector<double> & CurvePoint(const int i)
    { return *(curve_ips[i]); }

    const FlatVector<double> & CurvePoint(const int i) const
    { return *(curve_ips[i]); }

    FlatVector<double> & CurvePointTangent(const int i)
    { return *(curve_ip_tangents[i]); }

    const FlatVector<double> & CurvePointTangent(const int i) const
    { return *(curve_ip_tangents[i]); }

    int GetNumCurveParts(void) const;
    int GetStartOfCurve(const int i) const;
    int GetEndOfCurve(const int i) const;

    void AppendCurvePoint(const FlatVector<double> & point);
    void AppendCurvePoint(const FlatVector<double> & point, const FlatVector<double> & tangent);
    void SetCurveClearance(void);


  
    virtual void SetCacheComp(const int comp)
    { cachecomp = comp; }

    virtual int CacheComp(void) const
    { return cachecomp; }

    virtual void SetFileName(const string & filename)
    {
      cerr << "SetFileName not defined for Integrator base class" << endl;
    }
  };




  /**
     A BilinearFormIntegrator computes the element matrices. Different
     equations are provided by derived classes. An Integrator can be defined
     in the domain or at the boundary.
  */
  class NGS_DLL_HEADER BilinearFormIntegrator : public Integrator
  {
  public:
    // typedef double TSCAL;
    ///
    BilinearFormIntegrator () throw ();
    ///
    virtual ~BilinearFormIntegrator ();

    /// components of flux
    virtual int DimFlux () const { return -1; }

    /**
       Computes the element matrix.
    */
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
		       const ElementTransformation & eltrans, 
		       FlatMatrix<double> & elmat,
		       LocalHeap & locheap) const
    {
      AssembleElementMatrix (fel, eltrans, elmat, locheap);
    }

    /**
       Computes the element matrix.
       Complex version
    */
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
		       const ElementTransformation & eltrans, 
		       FlatMatrix<Complex> & elmat,
		       LocalHeap & locheap) const;


    /*
       Computes the element matrix. 
       should be renamed to  CalcElementMatrix
    */
    virtual void
    AssembleElementMatrix (const FiniteElement & fel,
			   const ElementTransformation & eltrans, 
			   FlatMatrix<double> & elmat,
			   LocalHeap & locheap) const;

    /*
    virtual void
    AssembleElementMatrix (const FiniteElement & fel,
			   const ElementTransformation & eltrans, 
			   FlatMatrix<Complex> & elmat,
			   LocalHeap & locheap) const;
    */



    virtual void
    AssembleElementMatrixIndependent (const FiniteElement & bfel_master,
				      const FiniteElement & bfel_master_element,				    
				      const FiniteElement & bfel_slave,
				      const ElementTransformation & eltrans_master, 
				      const ElementTransformation & eltrans_master_element, 
				      const ElementTransformation & eltrans_slave,
				      const IntegrationPoint & ip_master,
				      const IntegrationPoint & ip_master_element,
				      const IntegrationPoint & ip_slave,
				      FlatMatrix<double> & elmat,
				      LocalHeap & locheap) const
    {;}
    virtual void
    ApplyElementMatrixIndependent (const FiniteElement & bfel_master,
				   const FiniteElement & bfel_master_element,				    
				   const FiniteElement & bfel_slave,
				   const ElementTransformation & eltrans_master, 
				   const ElementTransformation & eltrans_master_element, 
				   const ElementTransformation & eltrans_slave,
				   const IntegrationPoint & ip_master,
				   const IntegrationPoint & ip_master_element,
				   const IntegrationPoint & ip_slave,
				   const FlatVector<double> & elx,
				   Vector<double> & result,
				   LocalHeap & locheap) const
    {;}
    virtual void
    AssembleElementMatrixIndependent (const FiniteElement & bfel_master,
				      const FiniteElement & bfel_master_element,				    
				      const FiniteElement & bfel_slave,
				      const ElementTransformation & eltrans_master, 
				      const ElementTransformation & eltrans_master_element, 
				      const ElementTransformation & eltrans_slave,
				      const IntegrationPoint & ip_master,
				      const IntegrationPoint & ip_master_element,
				      const IntegrationPoint & ip_slave,
				      FlatMatrix<Complex> & elmat,
				      LocalHeap & locheap) const
    {
      FlatMatrix<double> rmat;
      AssembleElementMatrixIndependent(bfel_master,bfel_master_element,bfel_slave,
				       eltrans_master, eltrans_master_element, eltrans_slave,
				       ip_master, ip_master_element, ip_slave,
				       rmat, locheap);
      elmat.AssignMemory(rmat.Height(), rmat.Width(), locheap);
      elmat = rmat;
    }

    virtual void
    AssembleElementMatrixIndependent (const FiniteElement & bfel_master,
				      const FiniteElement & bfel_slave,
				      const ElementTransformation & eltrans_master, 
				      const ElementTransformation & eltrans_slave,
				      const IntegrationPoint & ip_master,
				      const IntegrationPoint & ip_slave,
				      FlatMatrix<double> & elmat,
				      LocalHeap & locheap) const
    {;}
    virtual void
    AssembleElementMatrixIndependent (const FiniteElement & bfel_master,
				      const FiniteElement & bfel_slave,
				      const ElementTransformation & eltrans_master, 
				      const ElementTransformation & eltrans_slave,
				      const IntegrationPoint & ip_master,
				      const IntegrationPoint & ip_slave,
				      FlatMatrix<Complex> & elmat,
				      LocalHeap & locheap) const
    {
      FlatMatrix<double> rmat;
      AssembleElementMatrixIndependent(bfel_master,bfel_slave,
				       eltrans_master, eltrans_slave,
				       ip_master, ip_slave,
				       rmat, locheap);
      elmat.AssignMemory(rmat.Height(), rmat.Width(), locheap);
      elmat = rmat;
    }


    virtual void
    CalcElementMatrixDiag (const FiniteElement & fel,
			   const ElementTransformation & eltrans, 
			   FlatVector<double> & diag,
			   LocalHeap & lh) const;




    virtual void
    CalcLinearizedElementMatrix (const FiniteElement & fel, 
				 const ElementTransformation & eltrans, 
				 FlatVector<double> & elveclin,
				 FlatMatrix<double> & elmat,
				 LocalHeap & locheap) const;

    virtual void
    CalcLinearizedElementMatrix (const FiniteElement & fel, 
				 const ElementTransformation & eltrans, 
				 FlatVector<Complex> & elveclin,
				 FlatMatrix<Complex> & elmat,
				 LocalHeap & locheap) const;


    virtual void *  
    PrecomputeData (const FiniteElement & fel, 
		    const ElementTransformation & eltrans, 
		    LocalHeap & locheap) const { return 0; }
  

    virtual void 
    ApplyElementMatrix (const FiniteElement & fel, 
			const ElementTransformation & eltrans, 
			const FlatVector<double> & elx, 
			FlatVector<double> & ely,
			void * precomputed,
			LocalHeap & locheap) const;

    virtual void 
    ApplyElementMatrix (const FiniteElement & fel, 
			const ElementTransformation & eltrans, 
			const FlatVector<Complex> & elx, 
			FlatVector<Complex> & ely,
			void * precomputed,
			LocalHeap & locheap) const;


    template < int S, class T>
    void ApplyElementMatrix (const FiniteElement & fel, 
			     const ElementTransformation & eltrans, 
			     const FlatVector< Vec<S,T> > & elx, 
			     FlatVector< Vec<S,T> > & ely,
			     void * precomputed,
			     LocalHeap & locheap) const
    {
      //cout << "call baseclass ApplyElementMatrix, type = " << typeid(*this).name() << endl;
      FlatMatrix<T> mat;
      CalcElementMatrix (fel, eltrans, mat, locheap);
      ely = mat * elx;
    }


    virtual void 
    ApplyLinearizedElementMatrix (const FiniteElement & fel, 
				  const ElementTransformation & eltrans, 
				  const FlatVector<double> & ellin,
				  const FlatVector<double> & elx, 
				  FlatVector<double> & ely,
				  LocalHeap & locheap) const;

    virtual void 
    ApplyLinearizedElementMatrix (const FiniteElement & fel, 
				  const ElementTransformation & eltrans, 
				  const FlatVector<Complex> & ellin, 
				  const FlatVector<Complex> & elx, 
				  FlatVector<Complex> & ely,
				  LocalHeap & locheap) const;



    virtual double Energy (const FiniteElement & fel, 
			   const ElementTransformation & eltrans, 
			   const FlatVector<double> & elx, 
			   LocalHeap & locheap) const;

    virtual double Energy (const FiniteElement & fel, 
			   const ElementTransformation & eltrans, 
			   const FlatVector<Complex> & elx, 
			   LocalHeap & locheap) const;




    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseSpecificIntegrationPoint & bsip,
	      const FlatVector<double> & elx, 
	      FlatVector<double> & flux,
	      bool applyd,
	      LocalHeap & lh) const;


    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseSpecificIntegrationPoint & bsip,
	      const FlatVector<Complex> & elx, 
	      FlatVector<Complex> & flux,
	      bool applyd,
	      LocalHeap & lh) const;


    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseMappedIntegrationRule & mir,
	      const FlatVector<double> & elx, 
	      FlatMatrix<double> & flux,
	      bool applyd,
	      LocalHeap & lh) const;


    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseMappedIntegrationRule & mir,
	      const FlatVector<Complex> & elx, 
	      FlatMatrix<Complex> & flux,
	      bool applyd,
	      LocalHeap & lh) const;



    virtual void
    CalcFluxMulti (const FiniteElement & fel,
		   const BaseSpecificIntegrationPoint & bsip,
		   int m,
		   const FlatVector<double> & elx, 
		   FlatVector<double> & flux,
		   bool applyd,
		   LocalHeap & lh) const;


    virtual void
    CalcFluxMulti (const FiniteElement & fel,
		   const BaseSpecificIntegrationPoint & bsip,
		   int m,
		   const FlatVector<Complex> & elx, 
		   FlatVector<Complex> & flux,
		   bool applyd,
		   LocalHeap & lh) const;


    virtual void
    ApplyBTrans (const FiniteElement & fel,
		 const BaseSpecificIntegrationPoint & bsip,
		 const FlatVector<double> & elx, 
		 FlatVector<double> & ely,
		 LocalHeap & lh) const;

    virtual void
    ApplyBTrans (const FiniteElement & fel,
		 const BaseSpecificIntegrationPoint & bsip,
		 const FlatVector<Complex> & elx, 
		 FlatVector<Complex> & ely,
		 LocalHeap & lh) const;

    
    virtual void
    ApplyBTrans (const FiniteElement & fel,
		 const BaseMappedIntegrationRule & mir,
		 const FlatMatrix<double> & elx, 
		 FlatVector<double> & ely,
		 LocalHeap & lh) const;

    virtual void
    ApplyBTrans (const FiniteElement & fel,
		 const BaseMappedIntegrationRule & mir,
		 const FlatMatrix<Complex> & elx, 
		 FlatVector<Complex> & ely,
		 LocalHeap & lh) const;
    


    virtual void ApplyDMat (const FiniteElement & bfel,
			    const BaseSpecificIntegrationPoint & bsip,
			    const FlatVector<double> & elx, 
			    FlatVector<double> & eldx,
			    LocalHeap & lh) const;

    virtual void ApplyDMat (const FiniteElement & bfel,
			    const BaseSpecificIntegrationPoint & bsip,
			    const FlatVector<Complex> & elx, 
			    FlatVector<Complex> & eldx,
			    LocalHeap & lh) const;
  
    virtual void ApplyDMat (const FiniteElement & bfel,
			    const BaseMappedIntegrationRule & mir,
			    const FlatMatrix<double> & elx, 
			    FlatMatrix<double> & eldx,
			    LocalHeap & lh) const;

    virtual void ApplyDMat (const FiniteElement & bfel,
			    const BaseMappedIntegrationRule & mir,
			    const FlatMatrix<Complex> & elx, 
			    FlatMatrix<Complex> & eldx,
			    LocalHeap & lh) const;
  
    virtual void ApplyDMatInv (const FiniteElement & bfel,
			       const BaseSpecificIntegrationPoint & bsip,
			       const FlatVector<double> & elx, 
			       FlatVector<double> & eldx,
			       LocalHeap & lh) const;

    virtual void ApplyDMatInv (const FiniteElement & bfel,
			       const BaseSpecificIntegrationPoint & bsip,
			       const FlatVector<Complex> & elx, 
			       FlatVector<Complex> & eldx,
			       LocalHeap & lh) const;
  
    virtual void ApplyDMatInv (const FiniteElement & bfel,
			       const BaseMappedIntegrationRule & mir,
			       const FlatMatrix<double> & elx, 
			       FlatMatrix<double> & eldx,
			       LocalHeap & lh) const;

    virtual void ApplyDMatInv (const FiniteElement & bfel,
			       const BaseMappedIntegrationRule & mir,
			       const FlatMatrix<Complex> & elx, 
			       FlatMatrix<Complex> & eldx,
			       LocalHeap & lh) const;




    virtual const IntegrationRule & GetIntegrationRule (const FiniteElement & fel,
							const bool use_higher_integration_order = false) const;
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
      
    FacetBilinearFormIntegrator(Array<CoefficientFunction*> & coeffs) 
      : BilinearFormIntegrator() { ; }

    ~FacetBilinearFormIntegrator() { ; }

      
      
      
    virtual bool BoundaryForm () const 
    { return 0; }

    virtual bool SkeletonForm () const 
    { return 1; }
    
    virtual void CalcElementMatrix (const FiniteElement & fel,
				    const ElementTransformation & eltrans, 
				    FlatMatrix<double> & elmat,
				    LocalHeap & locheap) {
      throw Exception ("FacetBilinearFormIntegrator can not assemble volumetric element matrices!");
    }
    

    virtual void
    AssembleFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
			 const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
			 const FiniteElement & volumefel2, int LocalFacetNr2,
			 const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
			 FlatMatrix<double> & elmat,
			 LocalHeap & locheap) const{ 
      throw Exception ("FacetBilinearFormIntegrator::AssembleFacetMatrix for inner facets not implemented!");
    }
    virtual void 
    AssembleFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
			 const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
			 const FiniteElement & volumefel2, int LocalFacetNr2,
			 const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,	 
			 FlatMatrix<Complex> & elmat,
			 LocalHeap & locheap) const{ 
      throw Exception ("FacetBilinearFormIntegrator::AssembleFacetMatrix<Complex> for inner facets not implemented!");
    }


    virtual void
    AssembleFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
			 const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
			 const ElementTransformation & seltrans,  
			 FlatMatrix<double> & elmat,
			 LocalHeap & locheap) const{ 
      throw Exception ("FacetBilinearFormIntegrator::AssembleFacetMatrix for boundary facets not implemented!");
    }
    virtual void 
    AssembleFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
			 const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
			 const ElementTransformation & seltrans,  
			 FlatMatrix<Complex> & elmat,
			 LocalHeap & locheap) const{ 
      throw Exception ("FacetBilinearFormIntegrator::AssembleFacetMatrix<Complex> for boundary facets not implemented!");
    }


  };





  // returns coef * Identity-matrix
  class RegularizationBilinearFormIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * coef;
  public:
    RegularizationBilinearFormIntegrator (CoefficientFunction * acoef)
      : coef(acoef)
    { ; }

    virtual bool BoundaryForm () const { return 0; }
    virtual int DimFlux () const { return 1; }
    virtual int DimElement () const { return 3; }
    virtual int DimSpace () const { return 3; }

    virtual void
    CalcElementMatrix (const FiniteElement & fel,
		       const ElementTransformation & eltrans, 
		       FlatMatrix<double> & elmat,
		       LocalHeap & locheap) const
    {
      int ndof = fel.GetNDof();

      // elmat.AssignMemory (ndof, ndof, locheap);
      elmat = 0;
      SpecificIntegrationPoint<3,3> sip (GetIntegrationRules().SelectIntegrationRule (fel.ElementType(), 0)[0], 
					 eltrans, locheap);
      double val = Evaluate (*coef, sip);
      elmat = 0;
      for (int i = 0; i < ndof; i++)
	elmat(i,i) = val;
    }

    virtual string Name () const { return "Regularization"; }
  };
  




  class NGS_DLL_HEADER BlockBilinearFormIntegrator : public BilinearFormIntegrator
  {
    BilinearFormIntegrator & bfi;
    int dim;
    int comp;
  public:
    BlockBilinearFormIntegrator (BilinearFormIntegrator & abfi, int adim, int acomp);
    BlockBilinearFormIntegrator (BilinearFormIntegrator & abfi, int adim);
    virtual ~BlockBilinearFormIntegrator ();

    virtual bool BoundaryForm () const
    { return bfi.BoundaryForm(); }

    virtual int DimFlux () const 
    { return (comp == -1) ? dim * bfi.DimFlux() : bfi.DimFlux(); }

    const BilinearFormIntegrator & Block () const { return bfi; }

    virtual void
    CalcElementMatrix (const FiniteElement & bfel, 
		       const ElementTransformation & eltrans, 
		       FlatMatrix<double> & elmat,
		       LocalHeap & locheap) const;

    virtual void
    CalcElementMatrix (const FiniteElement & bfel, 
		       const ElementTransformation & eltrans, 
		       FlatMatrix<Complex> & elmat,
		       LocalHeap & locheap) const;
    
    virtual void 
    ApplyElementMatrix (const FiniteElement & bfel, 
			const ElementTransformation & eltrans, 
			const FlatVector<double> & elx, 
			FlatVector<double> & ely,
			void * precomputed,
			LocalHeap & locheap) const;

    virtual void 
    ApplyElementMatrix (const FiniteElement & bfel, 
			const ElementTransformation & eltrans, 
			const FlatVector<Complex> & elx, 
			FlatVector<Complex> & ely,
			void * precomputed,
			LocalHeap & locheap) const;

    virtual void
    CalcLinearizedElementMatrix (const FiniteElement & bfel,
				 const ElementTransformation & eltrans,
				 FlatVector<double> & elveclin,
				 FlatMatrix<double> & elmat,
				 LocalHeap & locheap) const;
    virtual void
    CalcLinearizedElementMatrix (const FiniteElement & bfel, 
				 const ElementTransformation & eltrans, 
				 FlatVector<Complex> & elveclin,
				 FlatMatrix<Complex> & elmat,
				 LocalHeap & locheap) const;


    virtual void
    CalcFlux (const FiniteElement & fel,
	      const ElementTransformation & eltrans,
	      const IntegrationPoint & ip,
	      const FlatVector<double> & elx, 
	      FlatVector<double> & flux,
	      bool applyd,
	      LocalHeap & lh) const;

    virtual void
    CalcFlux (const FiniteElement & fel,
	      const ElementTransformation & eltrans,
	      const IntegrationPoint & ip,
	      const FlatVector<Complex> & elx, 
	      FlatVector<Complex> & flux,
	      bool applyd,
	      LocalHeap & lh) const;



    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseSpecificIntegrationPoint & bsip,
	      const FlatVector<double> & elx, 
	      FlatVector<double> & flux,
	      bool applyd,
	      LocalHeap & lh) const;

    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseSpecificIntegrationPoint & bsip,
	      const FlatVector<Complex> & elx, 
	      FlatVector<Complex> & flux,
	      bool applyd,
	      LocalHeap & lh) const;

    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseMappedIntegrationRule & mir,
	      const FlatVector<double> & elx, 
	      FlatMatrix<double> & flux,
	      bool applyd,
	      LocalHeap & lh) const;



    virtual void
    ApplyBTrans (const FiniteElement & bfel,
		 const BaseSpecificIntegrationPoint & bsip,
		 const FlatVector<double> & elx, 
		 FlatVector<double> & ely,
		 LocalHeap & lh) const;

    virtual void
    ApplyBTrans (const FiniteElement & bfel,
		 const BaseSpecificIntegrationPoint & bsip,
		 const FlatVector<Complex> & elx, 
		 FlatVector<Complex> & ely,
		 LocalHeap & lh) const;

    virtual void
    ApplyBTrans (const FiniteElement & fel,
		 const BaseMappedIntegrationRule & mir,
		 const FlatMatrix<double> & elx, 
		 FlatVector<double> & ely,
		 LocalHeap & lh) const;

    virtual double Energy (const FiniteElement & fel, 
			   const ElementTransformation & eltrans, 
			   const FlatVector<double> & elx, 
			   LocalHeap & locheap) const;

    virtual string Name () const;
  };







  class NGS_DLL_HEADER ComplexBilinearFormIntegrator : public BilinearFormIntegrator
  {
    const BilinearFormIntegrator & bfi;
    Complex factor;
  public:
    ComplexBilinearFormIntegrator (const BilinearFormIntegrator & abfi,
				   Complex afactor);

    virtual bool BoundaryForm () const
    { return bfi.BoundaryForm(); }

    virtual int DimFlux () const
    { return bfi.DimFlux(); }
    virtual int DimElement () const
    { return bfi.DimElement(); }
    virtual int DimSpace () const
    { return bfi.DimSpace(); }

    virtual void GetFactor(Complex & fac) const {fac = factor;}
    virtual void GetFactor(double & fac) const {fac = factor.real();}
  
    virtual const BilinearFormIntegrator & GetBFI(void) const {return bfi;}

    virtual void CheckElement (const FiniteElement & el) const { bfi.CheckElement(el); }


    virtual void
    CalcElementMatrix (const FiniteElement & fel, 
		       const ElementTransformation & eltrans, 
		       FlatMatrix<double> & elmat,
		       LocalHeap & locheap) const;

    virtual void
    CalcElementMatrix (const FiniteElement & fel, 
		       const ElementTransformation & eltrans, 
		       FlatMatrix<Complex> & elmat,
		       LocalHeap & locheap) const;
    
    virtual void
    AssembleElementMatrixIndependent (const FiniteElement & bfel_master,
				      const FiniteElement & bfel_master_element,				    
				      const FiniteElement & bfel_slave,
				      const ElementTransformation & eltrans_master, 
				      const ElementTransformation & eltrans_master_element, 
				      const ElementTransformation & eltrans_slave,
				      const IntegrationPoint & ip_master,
				      const IntegrationPoint & ip_master_element,
				      const IntegrationPoint & ip_slave,
				      FlatMatrix<double> & elmat,
				      LocalHeap & locheap) const;

    virtual void
    AssembleElementMatrixIndependent (const FiniteElement & bfel_master,
				      const FiniteElement & bfel_master_element,				    
				      const FiniteElement & bfel_slave,
				      const ElementTransformation & eltrans_master, 
				      const ElementTransformation & eltrans_master_element, 
				      const ElementTransformation & eltrans_slave,
				      const IntegrationPoint & ip_master,
				      const IntegrationPoint & ip_master_element,
				      const IntegrationPoint & ip_slave,
				      FlatMatrix<Complex> & elmat,
				      LocalHeap & locheap) const;

    virtual void
    AssembleElementMatrixIndependent (const FiniteElement & bfel_master,
				      const FiniteElement & bfel_slave,
				      const ElementTransformation & eltrans_master, 
				      const ElementTransformation & eltrans_slave,
				      const IntegrationPoint & ip_master,
				      const IntegrationPoint & ip_slave,
				      FlatMatrix<double> & elmat,
				      LocalHeap & locheap) const;

    virtual void
    AssembleElementMatrixIndependent (const FiniteElement & bfel_master,
				      const FiniteElement & bfel_slave,
				      const ElementTransformation & eltrans_master, 
				      const ElementTransformation & eltrans_slave,
				      const IntegrationPoint & ip_master,
				      const IntegrationPoint & ip_slave,
				      FlatMatrix<Complex> & elmat,
				      LocalHeap & locheap) const;

  

    virtual void 
    ApplyElementMatrix (const FiniteElement & fel, 
			const ElementTransformation & eltrans, 
			const FlatVector<Complex> & elx, 
			FlatVector<Complex> & ely,
			void * precomputed,
			LocalHeap & locheap) const;

    virtual void
    CalcFlux (const FiniteElement & fel,
	      const ElementTransformation & eltrans,
	      const IntegrationPoint & ip,
	      const FlatVector<Complex> & elx, 
	      FlatVector<Complex> & flux,
	      bool applyd,
	      LocalHeap & lh) const;


    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseSpecificIntegrationPoint & bsip,
	      const FlatVector<Complex> & elx, 
	      FlatVector<Complex> & flux,
	      bool applyd,
	      LocalHeap & lh) const;

    virtual string Name () const;

  
    virtual const IntegrationRule & GetIntegrationRule (const FiniteElement & fel,
							const bool use_higher_integration_order = false) const;
  };






  class NGS_DLL_HEADER CompoundBilinearFormIntegrator : public BilinearFormIntegrator
  {
    const BilinearFormIntegrator & bfi;
    int comp;
  public:
    CompoundBilinearFormIntegrator (const BilinearFormIntegrator & abfi, int acomp);
  
    const BilinearFormIntegrator * GetBFI(void) const {return &bfi;}

    virtual bool BoundaryForm () const
    { return bfi.BoundaryForm(); }

    virtual int DimFlux () const 
    { return bfi.DimFlux(); }
    virtual int DimElement () const
    { return bfi.DimElement(); }
    virtual int DimSpace () const
    { return bfi.DimSpace(); }


    virtual void
    CalcElementMatrix (const FiniteElement & bfel, 
		       const ElementTransformation & eltrans, 
		       FlatMatrix<double> & elmat,
		       LocalHeap & locheap) const;

    virtual void
    CalcElementMatrix (const FiniteElement & bfel, 
		       const ElementTransformation & eltrans, 
		       FlatMatrix<Complex> & elmat,
		       LocalHeap & locheap) const;


    virtual void
    CalcLinearizedElementMatrix (const FiniteElement & fel, 
				 const ElementTransformation & eltrans, 
				 FlatVector<double> & elveclin,
				 FlatMatrix<double> & elmat,
				 LocalHeap & locheap) const;
    
    virtual void
    CalcLinearizedElementMatrix (const FiniteElement & fel, 
				     const ElementTransformation & eltrans, 
				     FlatVector<Complex> & elveclin,
				     FlatMatrix<Complex> & elmat,
				     LocalHeap & locheap) const;

    virtual void
    ApplyElementMatrix (const FiniteElement & bfel, 
			const ElementTransformation & eltrans, 
			const FlatVector<double> & elx,
			FlatVector<double> & ely,
			void * precomputed,
			LocalHeap & locheap) const;

    virtual void
    ApplyElementMatrix (const FiniteElement & bfel, 
			const ElementTransformation & eltrans, 
			const FlatVector<Complex> & elx,
			FlatVector<Complex> & ely,
			void * precomputed,
			LocalHeap & locheap) const;

    virtual void
    ApplyLinearizedElementMatrix (const FiniteElement & bfel, 
				  const ElementTransformation & eltrans, 
				  const FlatVector<double> & ellin,
				  const FlatVector<double> & elx,
				  FlatVector<double> & ely,
				  LocalHeap & locheap) const;

    virtual void
    ApplyLinearizedElementMatrix (const FiniteElement & bfel, 
				  const ElementTransformation & eltrans, 
				  const FlatVector<Complex> & ellin,
				  const FlatVector<Complex> & elx,
				  FlatVector<Complex> & ely,
				  LocalHeap & locheap) const;

    /*
    virtual void
    CalcFlux (const FiniteElement & bfel,
	      const ElementTransformation & eltrans,
	      const IntegrationPoint & ip,
	      const FlatVector<double> & elx, 
	      FlatVector<double> & flux,
	      bool applyd,
	      LocalHeap & lh) const;

    virtual void
    CalcFlux (const FiniteElement & bfel,
	      const ElementTransformation & eltrans,
	      const IntegrationPoint & ip,
	      const FlatVector<Complex> & elx, 
	      FlatVector<Complex> & flux,
	      bool applyd,
	      LocalHeap & lh) const;
    */

    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseSpecificIntegrationPoint & bsip,
	      const FlatVector<double> & elx, 
	      FlatVector<double> & flux,
	      bool applyd,
	      LocalHeap & lh) const;


    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseSpecificIntegrationPoint & bsip,
	      const FlatVector<Complex> & elx, 
	      FlatVector<Complex> & flux,
	      bool applyd,
	      LocalHeap & lh) const;


    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseMappedIntegrationRule & mir,
	      const FlatVector<double> & elx, 
	      FlatMatrix<double> & flux,
	      bool applyd,
	      LocalHeap & lh) const;


    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseMappedIntegrationRule & mir,
	      const FlatVector<Complex> & elx, 
	      FlatMatrix<Complex> & flux,
	      bool applyd,
	      LocalHeap & lh) const;


    virtual string Name () const;
  };






  /*

  template <int DIM_SPACE>
  class NGS_DLL_HEADER DirichletPenaltyIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * penalty;
  public:
    DirichletPenaltyIntegrator (CoefficientFunction * apenalty)
      : penalty(apenalty) { ; }

    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new DirichletPenaltyIntegrator (coeffs[0]);
    }

    virtual bool BoundaryForm () const
    { return 1; }

    virtual string Name () const
    { return "DirichletPenalty"; }

    virtual void
    CalcElementMatrix (const FiniteElement & fel, 
			   const ElementTransformation & eltrans, 
			   FlatMatrix<double> & elmat,
			   LocalHeap & locheap) const
    {
      int ndof = fel.GetNDof();

      // elmat.AssignMemory (ndof, ndof, locheap);
    
      const IntegrationRule & ir = GetIntegrationRules().SelectIntegrationRule (fel.ElementType(), 0);
      SpecificIntegrationPoint<DIM_SPACE-1,DIM_SPACE> sip (ir[0], eltrans, locheap);

      double val = Evaluate (*penalty, sip);

      elmat = 0;
      for (int i = 0; i < ndof; i++)
	elmat(i,i) = val;
    }


    virtual void 
    ApplyElementMatrix (const FiniteElement & fel, 
			const ElementTransformation & eltrans, 
			const FlatVector<Complex> & elx, 
			FlatVector<Complex> & ely,
			void * precomputed,
			LocalHeap & locheap) const
    {
      const IntegrationRule & ir = GetIntegrationRules().SelectIntegrationRule (fel.ElementType(), 0);
      SpecificIntegrationPoint<DIM_SPACE-1,DIM_SPACE> sip (ir[0], eltrans, locheap);
      double val = Evaluate (*penalty, sip);

      ely = val * elx;
    }
  };

  */



















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
		       FlatVector<double> & elvec,
		       LocalHeap & locheap) const
    {
      AssembleElementVector (fel, eltrans, elvec, locheap);
    }


    virtual void 
    CalcElementVector (const FiniteElement & fel,
		       const ElementTransformation & eltrans, 
		       FlatVector<Complex> & elvec,
		       LocalHeap & locheap) const;


    // old version:
    virtual void 
    AssembleElementVector (const FiniteElement & fel,
			   const ElementTransformation & eltrans, 
			   FlatVector<double> & elvec,
			   LocalHeap & locheap) const;

    virtual void 
    AssembleElementVector (const FiniteElement & fel,
			   const ElementTransformation & eltrans, 
			   FlatVector<Complex> & elvec,
			   LocalHeap & locheap) const;

    virtual void
    AssembleElementVectorIndependent (const FiniteElement & gfel,
				      const BaseSpecificIntegrationPoint & s_sip,
				      const BaseSpecificIntegrationPoint & g_sip,
				      FlatVector<double> & elvec,
				      LocalHeap & locheap,
				      const bool curveint = false) const
    {
      cerr << "AssembleElementVectorIndependent called for base-class!" << endl;
      exit(10);
    }
  
    virtual void
    AssembleElementVectorIndependent (const FiniteElement & gfel,
				      const BaseSpecificIntegrationPoint & s_sip,
				      const BaseSpecificIntegrationPoint & g_sip,
				      FlatVector<Complex> & elvec,
				      LocalHeap & locheap,
				      const bool curveint = false) const
    {
      FlatVector<double> rvec(elvec.Size(), locheap);
      AssembleElementVectorIndependent (gfel, s_sip, g_sip, rvec, locheap,curveint);
      elvec = rvec;
    }

  
  };















 class NGS_DLL_HEADER FacetLinearFormIntegrator : public LinearFormIntegrator
  {
    public:
      
    FacetLinearFormIntegrator(Array<CoefficientFunction*> & coeffs) 
      : LinearFormIntegrator() { ; }

    ~FacetLinearFormIntegrator() { ; }

    virtual bool BoundaryForm () const 
    { return 1; }

    virtual bool SkeletonForm () const 
    { return 1; }
    
    virtual void 
    CalcElementVector (const FiniteElement & bfel, 
			   const ElementTransformation & eltrans, 
			   FlatVector<double> & elvec,
			   LocalHeap & locheap) const{
      throw Exception ("FacetLinearFormIntegrator can not assemble volumetric element matrices!");
    }
    

    virtual void
    AssembleFacetVector (const FiniteElement & volumefel, int LocalFacetNr,
			 const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
			 const ElementTransformation & seltrans,
			 FlatVector<double> & elvec,
			 LocalHeap & locheap) const{ 
      throw Exception ("FacetLinearFormIntegrator::AssembleFacetVector not implemented!");
    }

    virtual void 
    AssembleFacetVector (const FiniteElement & volumefel, int LocalFacetNr,
			 const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
			 const ElementTransformation & seltrans,
			 FlatVector<Complex> & elvec,
			 LocalHeap & locheap) const { 
      throw Exception ("FacetLinearFormIntegrator::AssembleFacetVector<complex> not implemented!");
    }
			 
  };  





  class NGS_DLL_HEADER BlockLinearFormIntegrator : public LinearFormIntegrator
  {
    const LinearFormIntegrator & lfi;
    int dim;
    int comp;
  public:
    BlockLinearFormIntegrator (const LinearFormIntegrator & alfi, int adim, int acomp);

    virtual bool BoundaryForm () const
    { return lfi.BoundaryForm(); }


    virtual void 
    CalcElementVector (const FiniteElement & bfel, 
		       const ElementTransformation & eltrans, 
		       FlatVector<double> & elvec,
		       LocalHeap & locheap) const;
  };





  class NGS_DLL_HEADER ComplexLinearFormIntegrator : public LinearFormIntegrator
  {
    const LinearFormIntegrator & lfi;
    Complex factor;
  public:
    ComplexLinearFormIntegrator (const LinearFormIntegrator & alfi, 
				 Complex afactor)
      : lfi(alfi), factor(afactor)
    { ; }

    virtual bool BoundaryForm () const { return lfi.BoundaryForm(); } 

    virtual void CheckElement (const FiniteElement & el) const { lfi.CheckElement(el); }


    virtual void
    CalcElementVector (const FiniteElement & fel, 
		       const ElementTransformation & eltrans, 
		       FlatVector<double> & elvec,
		       LocalHeap & locheap) const
    {
      throw Exception ("ComplexLinearFormIntegrator: cannot assemble double vector");
    }

    virtual void
    CalcElementVector (const FiniteElement & fel, 
		       const ElementTransformation & eltrans, 
		       FlatVector<Complex> & elvec,
		       LocalHeap & locheap) const
    {
      FlatVector<Complex> rvec(elvec.Size(), locheap);
      lfi.CalcElementVector (fel, eltrans, rvec, locheap);
      elvec = factor * rvec;
    }  


    virtual void
    AssembleElementVectorIndependent (const FiniteElement & gfel, 
				      const BaseSpecificIntegrationPoint & s_sip,
				      const BaseSpecificIntegrationPoint & g_sip,
				      FlatVector<double> & elvec,
				      LocalHeap & locheap,
				      const bool curveint = false) const
    {
      throw Exception ("ComplexLinearFormIntegrator: cannot assemble double vector");
    }
  

    virtual void
    AssembleElementVectorIndependent (const FiniteElement & gfel, 
				      const BaseSpecificIntegrationPoint & s_sip,
				      const BaseSpecificIntegrationPoint & g_sip,
				      FlatVector<Complex> & elvec,
				      LocalHeap & locheap,
				      const bool curveint = false) const
    { 
      FlatVector<double> rvec;

      lfi.AssembleElementVectorIndependent (gfel, s_sip, g_sip,
					    rvec, locheap, curveint);
      elvec.AssignMemory (rvec.Size(), locheap);
      elvec = factor * rvec;
    }




    virtual string Name () const
    {
      return
	string ("ComplexIntegrator (") + lfi.Name() +
	string (")");
    }

  };



  
  class NGS_DLL_HEADER CompoundLinearFormIntegrator : public LinearFormIntegrator
  {
    const LinearFormIntegrator & lfi;
    int comp;
  public:
    CompoundLinearFormIntegrator (const LinearFormIntegrator & alfi, int acomp)
      : lfi(alfi), comp(acomp) { ; }

    virtual bool BoundaryForm () const
    { return lfi.BoundaryForm(); }


    virtual void 
    CalcElementVector (const FiniteElement & bfel, 
		       const ElementTransformation & eltrans, 
		       FlatVector<double> & elvec,
		       LocalHeap & locheap) const
    {
      const CompoundFiniteElement & fel =
	dynamic_cast<const CompoundFiniteElement&> (bfel);

      FlatVector<double> vec1(fel[comp].GetNDof(), locheap);
      lfi.CalcElementVector (fel[comp], eltrans, vec1, locheap);
    
      elvec.AssignMemory (fel.GetNDof(), locheap);
      elvec = 0;

      int base = 0;
      for (int i = 0; i < comp; i++)
	base += fel[i].GetNDof();

      for (int i = 0; i < vec1.Size(); i++)
	elvec(base+i) = vec1(i);
    }  


    virtual void 
    CalcElementVector (const FiniteElement & bfel, 
		       const ElementTransformation & eltrans, 
		       FlatVector<Complex> & elvec,
		       LocalHeap & locheap) const
    {
      const CompoundFiniteElement & fel =
	dynamic_cast<const CompoundFiniteElement&> (bfel);

      FlatVector<Complex> vec1(fel[comp].GetNDof(), locheap);
      lfi.CalcElementVector (fel[comp], eltrans, vec1, locheap);
    
      elvec.AssignMemory (fel.GetNDof(), locheap);
      elvec = 0;

      int base = 0;
      for (int i = 0; i < comp; i++)
	base += fel[i].GetNDof();

      for (int i = 0; i < vec1.Size(); i++)
	elvec(base+i) = vec1(i);
    }  

    virtual void
    AssembleElementVectorIndependent (const FiniteElement & gfel,
				      const BaseSpecificIntegrationPoint & s_sip,
				      const BaseSpecificIntegrationPoint & g_sip,
				      FlatVector<double> & elvec,
				      LocalHeap & locheap,
				      const bool curveint = false) const
    {
      const CompoundFiniteElement & fel =
	dynamic_cast<const CompoundFiniteElement&> (gfel);

      int i;
      FlatVector<double> vec1;
      lfi.AssembleElementVectorIndependent (fel[comp], s_sip, g_sip, vec1, locheap, curveint);
    
      elvec.AssignMemory (fel.GetNDof(), locheap);
      elvec = 0;

      int base = 0;
      for (i = 0; i < comp; i++)
	base += fel[i].GetNDof();

      for (i = 0; i < vec1.Size(); i++)
	elvec(base+i) = vec1(i);
    }
    virtual void
    AssembleElementVectorIndependent (const FiniteElement & gfel,
				      const BaseSpecificIntegrationPoint & s_sip,
				      const BaseSpecificIntegrationPoint & g_sip,
				      FlatVector<Complex> & elvec,
				      LocalHeap & locheap,
				      const bool curveint = false) const
    {
      const CompoundFiniteElement & fel =
	dynamic_cast<const CompoundFiniteElement&> (gfel);

      int i;
      FlatVector<Complex> vec1;
      lfi.AssembleElementVectorIndependent (fel[comp], s_sip, g_sip, vec1, locheap, curveint);
    
      elvec.AssignMemory (fel.GetNDof(), locheap);
      elvec = 0;

      int base = 0;
      for (i = 0; i < comp; i++)
	base += fel[i].GetNDof();

      for (i = 0; i < vec1.Size(); i++)
	elvec(base+i) = vec1(i);
    }

  


    virtual string Name () const
    {
      return
	string ("CompoundIntegrator (") + lfi.Name() +
	string (")");
    }
  };









  /// container for all integrators
  class NGS_DLL_HEADER Integrators
  {
  public:

    /// description of integrator
    class IntegratorInfo
    {
    public:
      string name;
      int spacedim;
      int numcoeffs;
      Integrator* (*creator)(Array<CoefficientFunction*> &);
    
      IntegratorInfo (const string & aname,
		      int aspacedim,
		      int anumcoffs,		      
		      Integrator* (*acreator)(Array<CoefficientFunction*>&));
    };
  

    Array<IntegratorInfo*> bfis;
    Array<IntegratorInfo*> lfis;
  
  public:
    ///
    Integrators();
    ///
    ~Integrators();  
    ///
    void AddBFIntegrator (const string & aname, int aspacedim, int anumcoeffs,
			  Integrator* (*acreator)(Array<CoefficientFunction*>&));
    ///
    void AddLFIntegrator (const string & aname, int aspacedim, int anumcoeffs,
			  Integrator* (*acreator)(Array<CoefficientFunction*>&));
  
    ///
    const Array<IntegratorInfo*> & GetBFIs() const { return bfis; }
    ///
    const IntegratorInfo * GetBFI(const string & name, int dim) const;
    ///
    BilinearFormIntegrator * CreateBFI(const string & name, int dim,
				       Array<CoefficientFunction*> & coeffs) const;
    ///
    BilinearFormIntegrator * CreateBFI(const string & name, int dim,
				       CoefficientFunction * coef) const;

    ///
    const Array<IntegratorInfo*> & GetLFIs() const { return lfis; }
    ///
    const IntegratorInfo * GetLFI(const string & name, int dim) const;
    ///
    LinearFormIntegrator * CreateLFI(const string & name, int dim,
				     Array<CoefficientFunction*> & coeffs) const;

    LinearFormIntegrator * CreateLFI(const string & name, int dim,
				     CoefficientFunction * coef) const;

    ///
    void Print (ostream & ost) const;
  };

  /// 
  extern NGS_DLL_HEADER Integrators & GetIntegrators ();






  template <typename BFI>
  class RegisterBilinearFormIntegrator
  {
  public:
    RegisterBilinearFormIntegrator (string label, int dim, int numcoeffs)
    {
      GetIntegrators().AddBFIntegrator (label, dim, numcoeffs, Create);
      cout << "register bf-integrator '" << label << "'" << endl;
    }
    
    static Integrator * Create (Array<CoefficientFunction*> & coefs)
    {
      return new BFI (coefs);
    }
  };



  template <typename LFI>
  class RegisterLinearFormIntegrator
  {
  public:
    RegisterLinearFormIntegrator (string label, int dim, int numcoeffs)
    {
      GetIntegrators().AddLFIntegrator (label, dim, numcoeffs, Create);
      cout << "register lf-integrator '" << label << "'" << endl;
    }
    
    static Integrator * Create (Array<CoefficientFunction*> & coefs)
    {
      return new LFI (coefs);
    }
  };




}

#endif
