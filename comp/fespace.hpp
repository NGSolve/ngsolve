#ifndef FILE_FESPACE
#define FILE_FESPACE

/*********************************************************************/
/* File:   fespace.hpp                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

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
class FESpace : public NGS_Object
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

  /// prolongation operators between multigrid levels
  ngmg::Prolongation *prol;

  /// on which subdomains is the space defined ?
  ARRAY<int> definedon;
  /// on which boundaries is the space defined ?
  ARRAY<int> definedonbound;

  // BEM is not supported
  // ARRAY<int> BEMboundary;

  /// prototype: what are the (homogeneous) Dirichlet boundaries ?
  BitArray dirichlet_boundaries;

  /// dofs on Dirichlet boundary
  BitArray dirichlet_dofs;

  
  
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
  ARRAY<int> first_hodofs[4];




  /// if directsolverclustered[i] is true, then the unknowns of domain i are clustered
  ARRAY<bool> directsolverclustered;

  ARRAY<string> directsolvermaterials;

  mutable ARRAY<int> adddirectsolverdofs;

  ARRAY<int> directvertexclusters;
  ARRAY<int> directedgeclusters;
  ARRAY<int> directfaceclusters;
  ARRAY<int> directelementclusters;

#ifdef PARALLEL
  class ngparallel::ParallelDofs * paralleldofs;
#endif
  
public:
  /*
  /// old constructor type, eliminated by JS Jan 2008
  FESpace (const MeshAccess & ama, int aorder,
	   int adim, bool acomplex, bool parseflags=false);
  */

  /**
     Constructor.
     Used flags are: \\
     -order=<int>:  finite element order \\
     -dim=<int>:    number of components \\
     -complex:      complex space \\
     -eliminate_internal:  eliminate internal dofs \\
     -dirichlet=<int-list>: dirichlet boundaries, 1-based \\
     
  */
  FESpace (const MeshAccess & ama, const Flags & flags, bool parseflags=false);
  ///
  virtual ~FESpace ();
  
  /// update dof-tables
  virtual void Update();

  /// update dof-table, preferred form
  virtual void Update(LocalHeap & lh);

  /// print report to stream
  virtual void PrintReport (ostream & ost);

  /// order of finite elements
  int GetOrder () const { return order; }

  /// how many components
  int GetDimension () const { return dimension; }

  /// complex space ?
  bool IsComplex () const { return iscomplex; }

  //
  // void SetBEM (bool abem);

  /// will be replaced by getclassname !
  virtual const char * GetType() 
  { return GetClassName().c_str(); }


  /// number of global dofs
  virtual int GetNDof () const = 0;
  /// number of global dofs on the level
  virtual int GetNDofLevel (int level) const;
  
  /// returns finite element. attention: should be thread-safe, but is not always
  virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
  /// get dof-nrs of the element
  virtual void GetDofNrs (int elnr, ARRAY<int> & dnums) const;
  /// get remaining dofs after static condensation
  virtual void GetExternalDofNrs (int elnr, ARRAY<int> & dnums) const;

  /// experiments with new preconditioners
  virtual void GetWireBasketDofNrs (int vnr, ARRAY<int> & dnums) const;

  /// get dofs on nr'th node of type nt.
  virtual void GetNodeDofNrs (NODE_TYPE nt, int nr, ARRAY<int> & dnums) const;
  /// get number of low-order dofs for node of type nt
  virtual int GetNLowOrderNodeDofs ( NODE_TYPE nt ) const
  { return lodofs_per_node[nt]; }

  /// get dofs on vertex vnr
  virtual void GetVertexDofNrs (int vnr, ARRAY<int> & dnums) const;
  /// get dofs on edge enr
  virtual void GetEdgeDofNrs (int ednr, ARRAY<int> & dnums) const;
  /// get dofs on face fnr
  virtual void GetFaceDofNrs (int fanr, ARRAY<int> & dnums) const;
  /// get dofs on element (=cell) elnr
  virtual void GetInnerDofNrs (int elnr, ARRAY<int> & dnums) const;


  /// returns surface element for boundary interals
  virtual const FiniteElement & GetSFE (int selnr, LocalHeap & lh) const;
  /// returns dofs of sourface element
  virtual void GetSDofNrs (int selnr, ARRAY<int> & dnums) const = 0;
  //
  // virtual void GetBEMDofNrs (ARRAY<int> & dnums) const;

  /// is the FESpace defined for this sub-domain nr ?
  bool DefinedOn (int domnr) const
  { return !definedon.Size() || definedon[domnr]; }
  /// 
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



  ///
  virtual void LockSomeDofs (BaseMatrix & mat) const { };

  /// old style
  virtual Table<int> * CreateSmoothingBlocks (int type = 0) const;

  /// 
  virtual Table<int> * CreateSmoothingBlocks (const Flags & flags) const
  { return CreateSmoothingBlocks(0); }

  /// for anisotropic plane smoothing, old style
  virtual BitArray * CreateIntermediatePlanes (int type = 0) const
  { return 0; }
  //virtual ARRAY<int> * CreateDirectSolverClusters (int type = 0) const
  //{ return 0; }
  virtual ARRAY<int> * CreateDirectSolverClusters (const Flags & flags) const
  { return 0; }
  //{ return CreateDirectSolverClusters(0); }

  virtual void AddDirectSolverClusterDof(int dn) const
  { adddirectsolverdofs.Append(dn); }

  virtual ARRAY<int> & DirectVertexClusters(void)
  { return directvertexclusters; }
  virtual ARRAY<int> & DirectEdgeClusters(void)
  { return directedgeclusters; }
  virtual ARRAY<int> & DirectFaceClusters(void)
  { return directfaceclusters; }
  virtual ARRAY<int> & DirectElementClusters(void)
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
  void SetProlongation (ngmg::Prolongation * aprol)
  { prol = aprol; }


  /// returns function-evaluator
  const BilinearFormIntegrator * GetEvaluator () const
  { return evaluator; }
  /// returns function-evaluator for boundary values
  const BilinearFormIntegrator * GetBoundaryEvaluator () const
  { return boundary_evaluator; }


  /// generates matrix graph
  virtual MatrixGraph * GetGraph (int level, bool symmetric);

  /// special elements for hacks (used for contact, periodic-boundary-penalty-constraints, ...
  ARRAY<SpecialElement*> specialelements;

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
};






/**
   A space of continuous finite elements.
   Supports first and second order finite elements.
 */
class NodalFESpace : public FESpace
{
  ///
  ARRAY<int> ndlevel;

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
  virtual void Update(LocalHeap & lh);
  ///
  virtual int GetNDof () const;
  ///
  virtual int GetNDofLevel (int level) const;
  ///
  virtual void GetDofNrs (int elnr, ARRAY<int> & dnums) const;
  ///
  virtual void GetSDofNrs (int selnr, ARRAY<int> & dnums) const;

  
  
  virtual void GetVertexDofNrs (int vnr, ARRAY<int> & dnums) const;
  virtual void GetEdgeDofNrs (int ednr, ARRAY<int> & dnums) const;
  virtual void GetFaceDofNrs (int fanr, ARRAY<int> & dnums) const;
  virtual void GetInnerDofNrs (int elnr, ARRAY<int> & dnums) const;

  //virtual ARRAY<int> * CreateDirectSolverClusters (int type = 0) const;
  virtual ARRAY<int> * CreateDirectSolverClusters (const Flags & flags) const;
#ifdef PARALLEL

   virtual void UpdateParallelDofs_hoproc();
   virtual void UpdateParallelDofs_loproc();

#endif

};






///
class NonconformingFESpace : public FESpace
{
  ///
  ARRAY<int> ndlevel;

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
  virtual void GetDofNrs (int elnr, ARRAY<int> & dnums) const;
  ///
  virtual void GetSDofNrs (int selnr, ARRAY<int> & dnums) const;
};




#ifdef OLD
/// use edge and face tables
class NodalFESpaceAlt : public FESpace
{
  ///
  ARRAY<int> ndlevel;
  ///
  int nv, ned, nfa;
public:

  ///
  NodalFESpaceAlt (const MeshAccess & ama,
		   int aorder, int adim, bool acomplex);

  ///
  ~NodalFESpaceAlt ();


  virtual string GetClassName () const
  {
    return "NodalFESpaceAlt";
  }

  ///
  virtual void Update();
  ///
  virtual int GetNDof () const;
  ///
  virtual int GetNDofLevel (int level) const;

  /*
  ///
  virtual const NodalFiniteElement & GetFE (int elnr, LocalHeap & lh) const
  {
    return static_cast<const NodalFiniteElement&> (FESpace::GetFE(elnr));
  }

  ///
  virtual const FiniteElement & GetSFE (int selnr) const
  {
    return static_cast<const NodalFiniteElement&> (FESpace::GetSFE(selnr));
  }
  */
  ///
  virtual void GetDofNrs (int elnr, ARRAY<int> & dnums) const;
  ///
  virtual void GetSDofNrs (int selnr, ARRAY<int> & dnums) const;


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

  /*
  ///
  virtual void TransformMatrix (int elnr, DenseMatrix & mat) const;
  ///
  virtual void TransformSurfMatrix (int elnr, DenseMatrix & mat) const;
  */
  ///
  virtual void LockSomeDofs (BaseMatrix & mat) const;
  ///
  virtual Table<int> * CreateSmoothingBlocks (int type = 0) const;
};
#endif



///
class ElementFESpace : public FESpace
{
  ///  ARRAY<int> startelement;
  ARRAY<int> ndlevel;
  int n_el_dofs;
public:

  ///
  /*
  ElementFESpace (const MeshAccess & ama,
		  int aorder, int adim, bool acomplex);
  */
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
  virtual int GetNDof () const
    { return ndlevel.Last(); }
  
  ///
  //  virtual const FiniteElement & GetFE (int elnr) const;
  ///
  virtual void GetDofNrs (int elnr, ARRAY<int> & dnums) const;

  ///
  virtual int GetNDofLevel (int level) const;

  ///
  // virtual const FiniteElement & GetSFE (int selnr) const;
  ///
  virtual void GetSDofNrs (int selnr, ARRAY<int> & dnums) const;


#ifdef PARALLEL
   virtual void UpdateParallelDofs_hoproc();
   virtual void UpdateParallelDofs_loproc();
#endif
};





/// Non-continous fe space on boundary
class SurfaceElementFESpace : public FESpace
{
  ///
  ARRAY<int> ndlevel;
  int n_el_dofs;
public:

  ///
  /*
  SurfaceElementFESpace (const MeshAccess & ama,
			 int aorder, int adim, bool acomplex);
  */
  ///
  SurfaceElementFESpace (const MeshAccess & ama, const Flags& flags, bool parseflags=false);

  ///
  ~SurfaceElementFESpace ();

  ///
  virtual string GetClassName() const
    { return "SurfaceElement"; }

  ///
  virtual void Update(LocalHeap & lh);

  ///
  virtual int GetNDof () const
    { return ndlevel.Last(); }

  ///
  virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
  ///
  virtual void GetDofNrs (int elnr, ARRAY<int> & dnums) const;

  ///
  virtual int GetNDofLevel (int level) const;

  ///
  virtual void GetSDofNrs (int selnr, ARRAY<int> & dnums) const;
};








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
  ARRAY<ngstd::INT<2> > faces;
  ///
  ARRAY<int[4]> elementfaces;
  ///
  ARRAY<int> surfelementfaces;

  ///
  ARRAY<int[5]> parentfaces;

  ///
  ARRAY<short int> finelevelofedge;
  ///
  ARRAY<int> nflevel;
  
public:
  ///
  //   NonConformingFESpace (const MeshAccess & ama,
  //			int aorder, int adim, bool acomplex);
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
  virtual void GetDofNrs (int elnr, ARRAY<int> & dnums) const;

  ///
  virtual const FiniteElement & GetSFE (int selnr, LocalHeap & lh) const;
  ///
  virtual void GetSDofNrs (int selnr, ARRAY<int> & dnums) const;

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












/// A combination of fe-spaces
class CompoundFESpace : public FESpace
{
protected:
  /// pointer to components
  ARRAY<const FESpace*> spaces;
  /// cummlated #dofs of components
  ARRAY<int> cummulative_nd;
  /// 
  ARRAY<int> ndlevel;
public:
  /*
  CompoundFESpace (const MeshAccess & ama,
		   const ARRAY<const FESpace*> & aspaces);
  */
  CompoundFESpace (const MeshAccess & ama,
		   const ARRAY<const FESpace*> & aspaces,
                   const Flags & flags, bool parseflags=false);
  ///
  virtual ~CompoundFESpace ();

  virtual string GetClassName () const
  {
    return "CompoundFESpace";
  }

  ///
  // virtual void Update();
  virtual void Update(LocalHeap & lh);
  ///
  virtual int GetNDof () const
  {return ndlevel.Last(); }
  ///
  virtual int GetNDofLevel (int level) const
  { return ndlevel[level]; }

  // returns start and end points of dofs corresponding to space "spacenr"
  // first space: spacenr = 0
  int GetStorageStart(int spacenr) const
  { return cummulative_nd[spacenr]; }

  int GetStorageEnd(int spacenr) const
  { return cummulative_nd[spacenr+1]; }



  const FESpace * operator[] (int i) const { return spaces[i]; }
  FESpace * operator[] (int i) { return const_cast<FESpace*> (spaces[i]); }

  virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
  ///
  virtual void GetDofNrs (int elnr, ARRAY<int> & dnums) const;
  virtual void GetExternalDofNrs (int elnr, ARRAY<int> & dnums) const;


  virtual void GetWireBasketDofNrs (int vnr, ARRAY<int> & dnums) const;
  virtual void GetVertexDofNrs (int vnr, ARRAY<int> & dnums) const;
  virtual void GetEdgeDofNrs (int ednr, ARRAY<int> & dnums) const;
  virtual void GetFaceDofNrs (int fanr, ARRAY<int> & dnums) const;
  virtual void GetInnerDofNrs (int elnr, ARRAY<int> & dnums) const;

  ///
  virtual const FiniteElement & GetSFE (int selnr, LocalHeap & lh) const;
  ///
  virtual void GetSDofNrs (int selnr, ARRAY<int> & dnums) const;


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
class FESpaceClasses
{
public:
  struct FESpaceInfo
  {
    string name;
    FESpace* (*creator)(const MeshAccess & ma, const Flags & flags);
    FESpaceInfo (const string & aname,
		 FESpace* (*acreator)(const MeshAccess & ma, const Flags & flags));
  };

  ARRAY<FESpaceInfo*> fesa;
public:
  FESpaceClasses();
  ~FESpaceClasses();  
  void AddFESpace (const string & aname, 
		   FESpace* (*acreator)(const MeshAccess & ma, const Flags & flags));
  
  const ARRAY<FESpaceInfo*> & GetFESpaces() { return fesa; }
  const FESpaceInfo * GetFESpace(const string & name);

  void Print (ostream & ost) const;
};
 
extern FESpaceClasses & GetFESpaceClasses ();






#endif
