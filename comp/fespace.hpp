#ifndef FILE_FESPACE
#define FILE_FESPACE

/*********************************************************************/
/* File:   fespace.hpp                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


namespace ngcomp
{

  /*
    Finite Element Space
  */


  /*
    transformation from local to global orientation
    used for low order Nedelec elements
  */
  enum TRANSFORM_TYPE { TRANSFORM_MAT_LEFT = 1,
			TRANSFORM_MAT_RIGHT = 2,
			TRANSFORM_MAT_LEFT_RIGHT = 3,
			TRANSFORM_RHS = 4,
			TRANSFORM_SOL = 8 };
		     

  /**
     Base class for finite element space.
     Provides finite elements, global degrees of freedom, 
     and transformations of element-matrices and element-vectors
  */
  class NGS_DLL_HEADER FESpace : public NGS_Object
  {
  protected:
    /// order of finite elements
    int order;
    /// how many components
    int dimension;
    /// complex space
    bool iscomplex;

    /// eliminate element-internal dofs ?
    bool eliminate_internal;

    /// couple (all) neighbouring degrees of freedom (like for jump terms of dg-methods)?
    bool dgjumps;

    /// prolongation operators between multigrid levels
    ngmg::Prolongation *prol;

    /// on which subdomains is the space defined ?
    Array<int> definedon;
    /// on which boundaries is the space defined ?
    Array<int> definedonbound;

    /// prototype: what are the Dirichlet boundaries ?
    BitArray dirichlet_boundaries;

    /// dofs on Dirichlet boundary
    BitArray dirichlet_dofs;
    BitArray free_dofs;


    Array<bool> dirichlet_vertex;
    Array<bool> dirichlet_edge;
    Array<bool> dirichlet_face;

  
    /// Reference - element (low order only)
    FiniteElement * tet;
    /// Reference - element (low order only)
    FiniteElement * prism;
    /// Reference - element (low order only) 
    FiniteElement * pyramid;
    /// Reference - element (low order only)
    FiniteElement * hex;
    /// Reference - element (low order only)
    FiniteElement * trig;
    /// Reference - element (low order only)
    FiniteElement * quad;
    /// Reference - element (low order only)
    FiniteElement * segm;

    /// Evaluator for visualization
    BilinearFormIntegrator * evaluator;
    /// Evaluator for visualization of boundary data
    BilinearFormIntegrator * boundary_evaluator;


    /// if non-zero, pointer to low order space
    FESpace * low_order_space;


    /// new node concept. primarily started for parallelization 
    /// how many low-order dofs per vertex, edge, face, cell ...
    int lodofs_per_node[4];
    /// vertex, edge, face, cell dofs start (and end) here ...
    int first_lodof[5];
    /// the high order vertex, edge, face, cell dofs on each node start (and end) here ...
    Array<int> first_hodofs[4];



    /// if directsolverclustered[i] is true, then the unknowns of domain i are clustered
    Array<bool> directsolverclustered;

    Array<string> directsolvermaterials;

    mutable Array<int> adddirectsolverdofs;

    Array<int> directvertexclusters;
    Array<int> directedgeclusters;
    Array<int> directfaceclusters;
    Array<int> directelementclusters;

    Table<int> * element_coloring;

#ifdef PARALLEL
    class ngparallel::ParallelDofs * paralleldofs;
#endif
  
  public:
    /**
       Constructor.
       Used flags are: \\
       -order=<int>:  finite element order \\
       -dim=<int>:    number of components \\
       -complex:      complex space \\
       -eliminate_internal:  eliminate internal dofs \\
       -dirichlet=<int-list>: dirichlet boundaries, 1-based \\
    */
    FESpace (const MeshAccess & ama, const Flags & flags, 
             bool checkflags = false);
    ///
    virtual ~FESpace ();
  
    /// update dof-tables, old style
    // virtual void Update();

    /// update dof-table
    virtual void Update(LocalHeap & lh);

    /// update element coloring
    virtual void FinalizeUpdate(LocalHeap & lh);
    const Table<int> & ElementColoring() const { return *element_coloring; }

    /// print report to stream
    virtual void PrintReport (ostream & ost);

    /// order of finite elements
    int GetOrder () const { return order; }

    /// how many components
    int GetDimension () const { return dimension; }

    /// complex space ?
    bool IsComplex () const { return iscomplex; }


    /// number of global dofs
    virtual int GetNDof () const = 0;
    /// number of global dofs on the level
    virtual int GetNDofLevel (int level) const;
  
    /// returns finite element. attention: should be thread-safe, but is not always
    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
    /// get dof-nrs of the element
    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;
    /// get remaining dofs after static condensation
    virtual void GetExternalDofNrs (int elnr, Array<int> & dnums) const;

    /// experiments with new preconditioners
    virtual void GetWireBasketDofNrs (int vnr, Array<int> & dnums) const;

    /// get dofs on nr'th node of type nt.
    virtual void GetNodeDofNrs (NODE_TYPE nt, int nr, Array<int> & dnums) const;
    /// get number of low-order dofs for node of type nt
    virtual int GetNLowOrderNodeDofs ( NODE_TYPE nt ) const
    { return lodofs_per_node[nt]; }

    /// get dofs on vertex vnr
    virtual void GetVertexDofNrs (int vnr, Array<int> & dnums) const;
    /// get dofs on edge enr
    virtual void GetEdgeDofNrs (int ednr, Array<int> & dnums) const;
    /// get dofs on face fnr
    virtual void GetFaceDofNrs (int fanr, Array<int> & dnums) const;
    /// get dofs on element (=cell) elnr
    virtual void GetInnerDofNrs (int elnr, Array<int> & dnums) const;

    virtual bool UsesDGCoupling () const {return dgjumps;};

    /// returns surface element for boundary interals
    virtual const FiniteElement & GetSFE (int selnr, LocalHeap & lh) const;
    /// returns dofs of sourface element
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const = 0;

    /// is the FESpace defined for this sub-domain nr ?
    bool DefinedOn (int domnr) const
    { return !definedon.Size() || definedon[domnr]; }
    /// is the FESpace defined for this boundary nr ?
    bool DefinedOnBoundary (int bnr) const
    {return !definedonbound.Size() || definedonbound[bnr]; }

    ///
    void SetDefinedOn (const BitArray & defon);
    ///
    void SetDefinedOnBoundary (const BitArray & defon);

    ///
    void SetDirichletBoundaries (const BitArray & dirbnds);
    /// Get reference element for tet, prism, trig, etc ..
    const FiniteElement & GetFE (ELEMENT_TYPE type) const;

    /// according low-order FESpace (if available)
    FESpace & LowOrderFESpace () { return *low_order_space; }
    /// according low-order FESpace (if available)
    const FESpace & LowOrderFESpace () const { return *low_order_space; }


    /// non Dirichlet dofs
    virtual const BitArray * GetFreeDofs () const;
    ///
    bool IsDirichletDof (int i) const
    { return dirichlet_dofs.Size() && dirichlet_dofs[i]; }

    bool IsDirichletBoundary (int i) const
    { return dirichlet_boundaries.Size() && dirichlet_boundaries[i]; }

    /// is vertex on Dirichlet boundary ?
    bool IsDirichletVertex (int i) const { return dirichlet_vertex.Size() && dirichlet_vertex[i]; }
    /// is edge on Dirichlet boundary ?
    bool IsDirichletEdge (int i) const { return dirichlet_edge.Size() && dirichlet_edge[i]; }
    /// is face on Dirichlet boundary ?
    bool IsDirichletFace (int i) const { return dirichlet_face.Size() && dirichlet_face[i]; }

    /// 
    virtual Table<int> * CreateSmoothingBlocks (const Flags & flags) const;
    // { return CreateSmoothingBlocks(0); }

    /// for anisotropic plane smoothing:
    virtual Array<int> * CreateDirectSolverClusters (const Flags & flags) const
    { return 0; }

    virtual void AddDirectSolverClusterDof(int dn) const
    { adddirectsolverdofs.Append(dn); }

    virtual Array<int> & DirectVertexClusters(void)
    { return directvertexclusters; }
    virtual Array<int> & DirectEdgeClusters(void)
    { return directedgeclusters; }
    virtual Array<int> & DirectFaceClusters(void)
    { return directfaceclusters; }
    virtual Array<int> & DirectElementClusters(void)
    { return directelementclusters; }


    void TransformMat (int elnr, bool boundary,
		       FlatMatrix<double> & mat, TRANSFORM_TYPE type) const
    {
      VTransformMR (elnr, boundary, mat, type);
    }
  
    void TransformMat (int elnr, bool boundary,
		       FlatMatrix<Complex> & mat, TRANSFORM_TYPE type) const
    {
      VTransformMC (elnr, boundary, mat, type);
    }
  


    void TransformVec (int elnr, bool boundary,
		       FlatVector<double> & vec, TRANSFORM_TYPE type) const
    {
      VTransformVR (elnr, boundary, vec, type);
    }
  
    void TransformVec (int elnr, bool boundary,
		       FlatVector<Complex> & vec, TRANSFORM_TYPE type) const
    {
      VTransformVC (elnr, boundary, vec, type);
    }


    template < int S, class T >
    void TransformVec (int elnr, bool boundary,
		       FlatVector< Vec<S,T> >& vec, TRANSFORM_TYPE type) const;

  

    virtual void VTransformMR (int elnr, bool boundary,
			       FlatMatrix<double> & mat, TRANSFORM_TYPE type) const
    { ; }
    virtual void VTransformMC (int elnr, bool boundary,
			       FlatMatrix<Complex> & mat, TRANSFORM_TYPE type) const
    { ; }
    virtual void VTransformVR (int elnr, bool boundary,
			       FlatVector<double> & vec, TRANSFORM_TYPE type) const
    { ; }
    virtual void VTransformVC (int elnr, bool boundary,
			       FlatVector<Complex> & vec, TRANSFORM_TYPE type) const
    { ; }
  
  
    /// Returns multigrid-prolongation
    virtual const ngmg::Prolongation * GetProlongation () const
    { return prol; }
    /// Set multigrid prolongation
    // void SetProlongation (ngmg::Prolongation * aprol)
    // { prol = aprol; }


    /// returns function-evaluator
    const BilinearFormIntegrator * GetEvaluator () const
    { return evaluator; }
    /// returns function-evaluator for boundary values
    const BilinearFormIntegrator * GetBoundaryEvaluator () const
    { return boundary_evaluator; }


    /// generates matrix graph
    virtual MatrixGraph * GetGraph (int level, bool symmetric);

    /// special elements for hacks (used for contact, periodic-boundary-penalty-constraints, ...
    Array<SpecialElement*> specialelements;

    void AppendSpecialElement (SpecialElement * spel)
    { specialelements.Append (spel); }



#ifdef PARALLEL
    virtual void UpdateParallelDofs ();
    virtual void UpdateParallelDofs ( LocalHeap & lh );
    virtual ngparallel::ParallelDofs & GetParallelDofs () const
    { return *paralleldofs; }

    virtual void ResetParallelDofs ();
  
    virtual void UpdateParallelDofs_hoproc();
    virtual void UpdateParallelDofs_loproc();
  
    MatrixGraph * GetConsistentGraph (int level, bool symmetric);
#endif


    virtual int GetRelOrder() const
    { 
      cout << "virtual GetRelOrder called for FiniteElementSpace, not available ! " << endl; 
      return 0; 
    } 

    virtual bool VarOrder() const { return 0; }

    bool timing;
    void Timing () const;
  };






  /**
     A space of continuous finite elements.
     Supports first and second order finite elements.
  */
  class NGS_DLL_HEADER NodalFESpace : public FESpace
  {
    ///
    Array<int> ndlevel;

  public:

    ///
    NodalFESpace (const MeshAccess & ama, const Flags & flags, bool parseflags=false);
    ///
    virtual ~NodalFESpace ();

    static FESpace * Create (const MeshAccess & ma, const Flags & flags);

    ///
    virtual string GetClassName () const
    {
      return "NodalFESpace";
    }

    ///
    virtual void Update (LocalHeap & lh);
    ///
    virtual int GetNDof () const;
    ///
    virtual int GetNDofLevel (int level) const;
    ///
    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;
    ///
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;
  
  
    virtual void GetVertexDofNrs (int vnr, Array<int> & dnums) const;
    virtual void GetEdgeDofNrs (int ednr, Array<int> & dnums) const;
    virtual void GetFaceDofNrs (int fanr, Array<int> & dnums) const;
    virtual void GetInnerDofNrs (int elnr, Array<int> & dnums) const;

    //virtual Array<int> * CreateDirectSolverClusters (int type = 0) const;
    virtual Array<int> * CreateDirectSolverClusters (const Flags & flags) const;

#ifdef PARALLEL
    virtual void UpdateParallelDofs_hoproc();
    virtual void UpdateParallelDofs_loproc();
#endif
  };






  ///
  class NGS_DLL_HEADER NonconformingFESpace : public FESpace
  {
    ///
    Array<int> ndlevel;

  public:
    NonconformingFESpace (const MeshAccess & ama, const Flags & flags, bool parseflags=false);
    virtual ~NonconformingFESpace ();

    virtual string GetClassName () const
    { return "Nonconforming FESpace"; }

    static FESpace * Create (const MeshAccess & ma, const Flags & flags)
    { return new NonconformingFESpace (ma, flags, true); }

    ///
    virtual void Update(LocalHeap & lh);
    ///
    virtual int GetNDof () const;
    ///
    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;
    ///
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;
  };







  ///
  class NGS_DLL_HEADER ElementFESpace : public FESpace
  {
    ///  Array<int> startelement;
    Array<int> ndlevel;
    int n_el_dofs;
  public:
    ///
    ElementFESpace (const MeshAccess & ama, const Flags& flags, bool parseflags=false);

    ///
    ~ElementFESpace ();

    virtual string GetClassName () const
    {
      return "ElementFESpace";
    }

    ///
    virtual void Update(LocalHeap & lh);

    ///
    virtual int GetNDof () const { return ndlevel.Last(); }
  
    ///
    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;

    ///
    virtual int GetNDofLevel (int level) const;

    ///
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;

#ifdef PARALLEL
    virtual void UpdateParallelDofs_hoproc();
    virtual void UpdateParallelDofs_loproc();
#endif
  };





  /// Non-continous fe space on boundary
  class NGS_DLL_HEADER SurfaceElementFESpace : public FESpace
  {
    ///
    Array<int> ndlevel;
    int n_el_dofs;
  public:
    ///
    SurfaceElementFESpace (const MeshAccess & ama, const Flags& flags, 
                           bool checkflags = false);

    ///
    ~SurfaceElementFESpace ();

    ///
    virtual string GetClassName() const
    { return "SurfaceElement"; }

    ///
    virtual void Update(LocalHeap & lh);

    ///
    virtual int GetNDof () const { return ndlevel.Last(); } 

    ///
    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;

    ///
    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;

    ///
    virtual int GetNDofLevel (int level) const;

    ///
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;
  };





  /*


  ///
  class NonConformingFESpace : public FESpace
  {
    ///
    FE_NcSegm1 segm1;
    ///
    FE_NcTrig1 trig1;
    ///
    FE_NcTet1 tet1;

    ///
    HashTable<ngstd::INT<2>,int> *node2face2d;
    ///
    HashTable<ngstd::INT<3>,int> *node2face3d;
    ///
    Array<ngstd::INT<2> > faces;
    ///
    Array<int[4]> elementfaces;
    ///
    Array<int> surfelementfaces;

    ///
    Array<int[5]> parentfaces;

    ///
    Array<short int> finelevelofedge;
    ///
    Array<int> nflevel;
  
  public:
    ///
    NonConformingFESpace (const MeshAccess & ama, const Flags& flags, bool parseflags=false);

    ///
    ~NonConformingFESpace ();

    ///
    virtual string GetClassName() const
    { return "Non-conforming"; }

    ///
    virtual void Update(LocalHeap & lh);

    ///
    virtual int GetNDof () const;
    ///
    virtual int GetNDofLevel (int level) const;

    ///
    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
    ///
    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;

    ///
    virtual const FiniteElement & GetSFE (int selnr, LocalHeap & lh) const;
    ///
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;

    ///
    int GetFacePoint1 (int fnr) const { return faces[fnr][0]; }
    ///
    int GetFacePoint2 (int fnr) const { return faces[fnr][1]; }

    ///
    int GetParentFace1 (int fnr) const { return parentfaces[fnr][0]; }
    ///
    int GetParentFace2 (int fnr) const { return parentfaces[fnr][1]; }
    ///
    int GetParentFace3 (int fnr) const { return parentfaces[fnr][2]; }
    ///
    int GetParentFace4 (int fnr) const { return parentfaces[fnr][3]; }
    ///
    int GetParentFace5 (int fnr) const { return parentfaces[fnr][4]; }
    ///
    int GetFineLevelOfFace (int ednr) const { return finelevelofedge[ednr]; }
  };

  */










  /// A combination of fe-spaces
  class NGS_DLL_HEADER CompoundFESpace : public FESpace
  {
  protected:
    /// pointer to components
    Array<const FESpace*> spaces;
    /// cummlated #dofs of components
    Array<int> cummulative_nd;
    /// 
    Array<int> ndlevel;
  public:
    ///
    CompoundFESpace (const MeshAccess & ama,
		     const Array<const FESpace*> & aspaces,
		     const Flags & flags, bool parseflags=false);
    ///
    virtual ~CompoundFESpace ();

    ///
    virtual string GetClassName () const
    {
      return "CompoundFESpace";
    }

    ///
    virtual void Update(LocalHeap & lh);
    ///
    virtual int GetNDof () const
    { return ndlevel.Last(); }
    ///
    virtual int GetNDofLevel (int level) const
    { return ndlevel[level]; }

    // returns start and end points of dofs corresponding to space "spacenr"
    // first space: spacenr = 0
    int GetStorageStart(int spacenr) const
    { return cummulative_nd[spacenr]; }
    
    ///
    int GetStorageEnd(int spacenr) const
    { return cummulative_nd[spacenr+1]; }


    ///
    const FESpace * operator[] (int i) const { return spaces[i]; }
    FESpace * operator[] (int i) { return const_cast<FESpace*> (spaces[i]); }

    ///
    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
    ///
    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;
    virtual void GetExternalDofNrs (int elnr, Array<int> & dnums) const;


    virtual void GetWireBasketDofNrs (int vnr, Array<int> & dnums) const;
    virtual void GetVertexDofNrs (int vnr, Array<int> & dnums) const;
    virtual void GetEdgeDofNrs (int ednr, Array<int> & dnums) const;
    virtual void GetFaceDofNrs (int fanr, Array<int> & dnums) const;
    virtual void GetInnerDofNrs (int elnr, Array<int> & dnums) const;

    ///
    virtual const FiniteElement & GetSFE (int selnr, LocalHeap & lh) const;
    ///
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;


    template <class MAT>
    void TransformMat (int elnr, bool boundary,
		       MAT & mat, TRANSFORM_TYPE tt) const;

    template <class VEC>
    void TransformVec (int elnr, bool boundary,
		       VEC & vec, TRANSFORM_TYPE tt) const;


    virtual void VTransformMR (int elnr, bool boundary,
			       FlatMatrix<double> & mat, TRANSFORM_TYPE tt) const 
    {
      TransformMat (elnr, boundary, mat, tt);
    }

    virtual void VTransformMC (int elnr, bool boundary,
			       FlatMatrix<Complex> & mat, TRANSFORM_TYPE tt) const
    {
      TransformMat (elnr, boundary, mat, tt);
    }

    virtual void VTransformVR (int elnr, bool boundary,
			       FlatVector<double> & vec, TRANSFORM_TYPE tt) const 
    {
      TransformVec (elnr, boundary, vec, tt);
    }

    virtual void VTransformVC (int elnr, bool boundary,
			       FlatVector<Complex> & vec, TRANSFORM_TYPE tt) const 
    {
      TransformVec (elnr, boundary, vec, tt);
    }

    inline int GetNSpaces () const { return spaces.Size(); } 
#ifdef PARALLEL

    virtual void UpdateParallelDofs_hoproc();
    virtual void UpdateParallelDofs_loproc();

#endif
  };





#ifdef PARALLEL
  class ParallelElementFESpace : public ElementFESpace
  {
  public:

    ///
    ParallelElementFESpace (const MeshAccess & ama,
			    int aorder, int adim, bool acomplex);
    ///
    ParallelElementFESpace (const MeshAccess & ama, const Flags& flags, bool parseflags=false);

    ///
    ~ParallelElementFESpace ()
    {
      ;
    }

    virtual string GetClassName () const
    {
      return "ParallelElementFESpace";
    }


    virtual void UpdateParallelDofs_hoproc();
    virtual void UpdateParallelDofs_loproc();
  };


  class ParallelNodalFESpace : public NodalFESpace
  {

  public:

    ///
    ParallelNodalFESpace (const MeshAccess & ama, const Flags & flags, bool parseflags=false);
    ///
    virtual ~ParallelNodalFESpace ()
    {
      ;
    }

    ///
    virtual string GetClassName () const
    {
      return "ParallelNodalFESpace";
    }


    virtual void UpdateParallelDofs_hoproc();
    virtual void UpdateParallelDofs_loproc();

  };

#endif


  /// Registered FESpace classes
  class NGS_DLL_HEADER FESpaceClasses
  {
  public:
    struct FESpaceInfo
    {
      string name;
      FESpace* (*creator)(const MeshAccess & ma, const Flags & flags);
      FESpaceInfo (const string & aname,
		   FESpace* (*acreator)(const MeshAccess & ma, const Flags & flags));
    };

    Array<FESpaceInfo*> fesa;
  public:
    FESpaceClasses();
    ~FESpaceClasses();  
    void AddFESpace (const string & aname, 
		     FESpace* (*acreator)(const MeshAccess & ma, const Flags & flags));
  
    const Array<FESpaceInfo*> & GetFESpaces() { return fesa; }
    const FESpaceInfo * GetFESpace(const string & name);

    void Print (ostream & ost) const;
  };
 
  extern NGS_DLL_HEADER FESpaceClasses & GetFESpaceClasses ();


}



#endif
