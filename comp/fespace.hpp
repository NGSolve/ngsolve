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


  /**
    transformation from local to global orientation
    used for low order Nedelec elements
  */
  enum TRANSFORM_TYPE { TRANSFORM_MAT_LEFT = 1,
			TRANSFORM_MAT_RIGHT = 2,
			TRANSFORM_MAT_LEFT_RIGHT = 3,
			TRANSFORM_RHS = 4,
			TRANSFORM_SOL = 8 };
  /**
    coupling types: Each degree of freedom is either
     - a local degree of freedom 
     - an interface degree of freedom 
     or
     - a wirebasket degree of freedom
  */
  enum COUPLING_TYPE {  UNUSED_DOF = 0,
			LOCAL_DOF = 1,
			INTERFACE_DOF = 2,
			NONWIREBASKET_DOF = 3,
			WIREBASKET_DOF = 4,
			EXTERNAL_DOF = 6,
			ANY_DOF = 7
		      };
		     

  /**
     constant_order .... one order for everything
     node_type_order ... order for edges, or faces, gradients or curls, but all over the mesh
     variable_order .... a different, anisotropic order for every mesh node
   */
  enum ORDER_POLICY { CONSTANT_ORDER = 0, NODE_TYPE_ORDER = 1, VARIABLE_ORDER = 2, OLDSTYLE_ORDER = 3 };

  
  
  NGS_DLL_HEADER ostream & operator<< (ostream & ost, COUPLING_TYPE ct);


  class FESpace;

  // will be size_t some day 
  typedef int DofId;

  using ngmg::Prolongation;

  /**
     Base class for finite element space.
     Provides finite elements, global degrees of freedom, 
     and transformations of element-matrices and element-vectors
  */
  class NGS_DLL_HEADER FESpace : public NGS_Object
  {
  protected:
    /// global order of finite elements
    int order;
    /// how many components
    int dimension;
    /// complex space
    bool iscomplex;

    /// couple (all) neighbouring degrees of freedom (like for jump terms of dg-methods)?
    bool dgjumps;

    /// debug output to testout
    bool print; 

    /// prolongation operators between multigrid levels
    shared_ptr<Prolongation> prol;// = NULL;
    /// highest multigrid-level for which Update was called (memory allocation)
    int level_updated;

    /// on which subdomains is the space defined ?
    Array<bool> definedon[3];

    /// prototype: what are the Dirichlet boundaries ?
    BitArray dirichlet_boundaries;

    /// dofs on Dirichlet boundary
    BitArray dirichlet_dofs;
    shared_ptr<BitArray> free_dofs;
    shared_ptr<BitArray> external_free_dofs;


    Array<bool> dirichlet_vertex;
    Array<bool> dirichlet_edge;
    Array<bool> dirichlet_face;

  
    /// Reference - element (low order only)
    FiniteElement * tet;  // = NULL;
    /// Reference - element (low order only)
    FiniteElement * prism; // = NULL;
    /// Reference - element (low order only) 
    FiniteElement * pyramid;  // = NULL;
    /// Reference - element (low order only)
    FiniteElement * hex; //  = NULL;
    /// Reference - element (low order only)
    FiniteElement * trig; // = NULL;
    /// Reference - element (low order only)
    FiniteElement * quad;// = NULL;
    /// Reference - element (low order only)
    FiniteElement * segm;// = NULL;
    /// Reference - element (low order only)
    FiniteElement * point;// = NULL;


    FiniteElement * dummy_tet; // = new <DummyFE<ET_TET>();
    FiniteElement * dummy_pyramid; // = new DummyFE<ET_PYRAMID>();
    FiniteElement * dummy_prism; // = new DummyFE<ET_PRISM>();
    FiniteElement * dummy_hex; //  = new DummyFE<ET_HEX>();
    FiniteElement * dummy_trig; // = new DummyFE<ET_TRIG>();
    FiniteElement * dummy_quad; // = new DummyFE<ET_QUAD>();
    FiniteElement * dummy_segm; // = new DummyFE<ET_SEGM>();
    FiniteElement * dummy_point; // = new DummyFE<ET_POINT>();

    /// Evaluator for visualization (new style)
    shared_ptr<DifferentialOperator> evaluator[3];
    /// Evaluator for flux
    shared_ptr<DifferentialOperator> flux_evaluator[3];
    /// Evaluator for visualization (old style)
    shared_ptr<BilinearFormIntegrator> integrator[3];

    /// if non-zero, pointer to low order space
    shared_ptr<FESpace> low_order_space; 

    /// if directsolverclustered[i] is true, then the unknowns of domain i are clustered
    Array<bool> directsolverclustered;

    Array<string> directsolvermaterials;

    mutable Array<int> adddirectsolverdofs;

    Array<int> directvertexclusters;
    Array<int> directedgeclusters;
    Array<int> directfaceclusters;
    Array<int> directelementclusters;

    
    Table<int> element_coloring[3]; 
    Table<int> facet_coloring;  // elements on facet in own colors (DG)
    Array<COUPLING_TYPE> ctofdof;

    ParallelDofs * paralleldofs; // = NULL;

    bool no_low_order_space;

    int et_bonus_order[30]; // order increase for element-type

    typedef int8_t TORDER;

    ORDER_POLICY order_policy = OLDSTYLE_ORDER;
    
    /*
      the function space H(curl) has high order basis funcitons which 
      are gradients, and additional ones which span the domain of the curl.
      In general, for function spaces of the de Rham sequence we refer to 
      functions in the range of the differential operator from the left to 
      the left sub-space, and functions spanning the domain of the right 
      differential operator as the right sub-space. 
      
      We can give different polynomial orders to the left and the
      right sub-space.  This allows to define Raviart-Thomas vs BDM
      elements, or Nedelec-type 1 vs Nedelec-type 2 elements, and
      more: We can skip all high order gradients in H(curl), or also
      define a (high-order) div-free H(div) space.

      We can give different orders for the left and right space for different
      element-types (like trig or quad)
     */
    int et_order_left[30];  // order for range of diff-op from the left 
    int et_order_right[30]; // order for domain of diff-op to the right
    
    Array<TORDER> order_edge; 
    Array<INT<2,TORDER>> order_face_left;
    Array<INT<2,TORDER>> order_face_right; 
    Array<INT<3,TORDER>> order_cell_left;
    Array<INT<3,TORDER>> order_cell_right;
    size_t order_timestamp = 0;
    BitArray is_atomic_dof;

    
    // move ndof and ndof_level to FESpace base class
  private:
    size_t ndof;
    Array<size_t> ndof_level;
  protected:
    void SetNDof (size_t _ndof);
    
  public:
    virtual int GetSpacialDimension() const { return ma->GetDimension();}
    string type;

    /**
       Constructor.
       Used flags are: \\
       -order=<int>:  finite element order \\
       -dim=<int>:    number of components \\
       -complex:      complex space \\
       -dirichlet=<int-list>: dirichlet boundaries, 1-based \\
    */
    FESpace (shared_ptr<MeshAccess> ama, const Flags & flags, 
             bool checkflags = false);
    /// cleanup
    virtual ~FESpace ();

    /// update dof-table
    virtual void Update(LocalHeap & lh);

    /// update element coloring
    virtual void FinalizeUpdate(LocalHeap & lh);

    /// highest level where update/finalize was called
    int GetLevelUpdated() const { return level_updated; }

    const Table<int> & ElementColoring(VorB vb = VOL) const 
    { return element_coloring[vb]; }

    const Table<int> & FacetColoring() const;
    
    /// print report to stream
    virtual void PrintReport (ostream & ost) const;

    /// Dump/restore fespace
    virtual void DoArchive (Archive & archive);

    /// order of finite elements
    int GetOrder () const { return order; }

    void SetBonusOrder (ELEMENT_TYPE et, int bonus) 
    { et_bonus_order[et] = bonus; }

    void SetOrderPolicy (ORDER_POLICY op)
    {
      order_policy = op;
    }
    
    void SetOrder (ELEMENT_TYPE et, TORDER order)
    {
      if (order_policy == CONSTANT_ORDER || order_policy == OLDSTYLE_ORDER)
        order_policy = NODE_TYPE_ORDER;
      et_order_left[et] = et_order_right[et] = order;
    }
    void SetOrderLeft (ELEMENT_TYPE et, TORDER order)
    {
      if (order_policy == CONSTANT_ORDER || order_policy == OLDSTYLE_ORDER)
        order_policy = NODE_TYPE_ORDER;      
      et_order_left[et] = order;
    }
    void SetOrderRight (ELEMENT_TYPE et, TORDER order)
    {
      if (order_policy == CONSTANT_ORDER || order_policy == OLDSTYLE_ORDER)
        order_policy = NODE_TYPE_ORDER;
      et_order_right[et] = order;
    }

    void SetOrder (NodeId ni, TORDER order)
    {
      switch (ni.GetType())
        {
        case NT_VERTEX:
          break;
        case NT_EDGE:
          if (ni.GetNr() < order_edge.Size())
            order_edge[ni.GetNr()] = order;
          break;
        case NT_FACE:
          if (ni.GetNr() < order_face_left.Size())
            order_face_left[ni.GetNr()] = order;
          if (ni.GetNr() < order_face_right.Size())
            order_face_right[ni.GetNr()] = order;
          break;
        case NT_CELL:
          // not yet 
          break;
        case NT_ELEMENT: case NT_FACET:
          break;
        }
    }
    /// how many components
    int GetDimension () const { return dimension; }

    /// complex space ?
    bool IsComplex () const { return iscomplex; }


    /// number of (process-local) dofs
    virtual size_t GetNDof () const { return ndof; } 
    /// number of dofs on the level
    virtual size_t GetNDofLevel (int level) const { return ndof_level[level]; } 

    SymbolTable<shared_ptr<DifferentialOperator>> additional_evaluators;
       
    class Element : public Ngs_Element
    {
      const FESpace & fes;
      Array<DofId> & temp_dnums;
      LocalHeap & lh;
      mutable bool dofs_set = false;
    public:     
      INLINE Element (const FESpace & afes, ElementId id, Array<DofId> & atemp_dnums,
                      LocalHeap & alh)
        : Ngs_Element ((*afes.GetMeshAccess())[id] ), fes(afes), 
          temp_dnums(atemp_dnums), lh(alh) 
      { ; }

      INLINE Element (const Element & el) = default;
      INLINE Element (Element && el) = default;
      // ElementId operator() const { return ei; }

      INLINE FlatArray<DofId> GetDofs() const
      {
        if (!dofs_set)
          fes.GetDofNrs (*this, temp_dnums);
        dofs_set = true;
        return temp_dnums;
      }

      INLINE const ElementTransformation & GetTrafo() const
      {
        return fes.GetMeshAccess()->GetTrafo (ElementId(*this), lh);
      }

      INLINE const FiniteElement & GetFE() const
      {
        return fes.GetFE (ElementId(*this), lh);
      }

      INLINE LocalHeap & GetLH() const
      {
        return lh;
      }
    };

    class ElementIterator
    {
      const FESpace & fes;
      ElementId ei;
      const FlatArray<bool> defined_on;
      Array<DofId> & temp_dnums;      
      LocalHeap & lh;
      void * heappointer;
    public:
      INLINE ElementIterator (const FESpace & afes, ElementId aei, 
                              const FlatArray<bool> adefined_on,
                              Array<DofId> & atemp_dnums, LocalHeap & alh)
        : fes(afes), ei(aei), defined_on(adefined_on), 
          temp_dnums(atemp_dnums), lh(alh), heappointer(lh.GetPointer()) { ; }
      INLINE ElementIterator & operator++ ()
      {
        lh.CleanUp(heappointer);
        ++ei;
        while (ei.Nr() < fes.GetMeshAccess()->GetNE(VorB(ei)) && 
               (defined_on.Size() && 
                !defined_on[fes.GetMeshAccess()->GetElIndex(ei)])
               ) ++ei;
        return *this;
      }
      INLINE Element operator*() const { return Element (fes, ei, temp_dnums, lh); }          
      INLINE bool operator!=(const ElementIterator & id2) const { return ei != id2.ei; }
      INLINE bool operator==(const ElementIterator & id2) const { return ei == id2.ei; }
    };
    
    class ElementRange : public IntRange
    {
      const FESpace & fes;
      Array<bool> definedon;
      const VorB vb;
      mutable Array<DofId> temp_dnums;
      mutable LocalHeap mylh;
      LocalHeap & lh;
    public:
      INLINE ElementRange (const FESpace & afes, VorB avb, IntRange ar, LocalHeap && lh2) 
        : IntRange(ar), fes(afes),
          definedon(fes.definedon[avb].Size(),fes.definedon[avb].Addr(0)),
          // FlatArray<bool>(fes.definedon) : FlatArray<bool>(fes.definedonbound)), 
          vb(avb), mylh(move(lh2)), lh(mylh)
      { ; }

      INLINE ElementRange (const FESpace & afes, VorB avb, IntRange ar, LocalHeap & lh2) 
        : IntRange(ar), fes(afes), 
          definedon(fes.definedon[avb].Size(),fes.definedon[avb].Addr(0)),
          // definedon( (avb==VOL) ? FlatArray<bool> (fes.definedon) : FlatArray<bool> (fes.definedonbound)), 
          vb(avb), mylh(), lh(lh2)
      { ; }

      ElementRange (const ElementRange & r2) = delete;

      INLINE ElementRange (ElementRange && r2) 
        : IntRange(r2), fes(r2.fes), definedon(move(r2.definedon)), vb(r2.vb), 
          temp_dnums(move(r2.temp_dnums)), mylh(move(r2.mylh)), 
          lh( (&r2.mylh == &r2.lh) ? mylh : r2.lh)
      { ; }

      INLINE ~ElementRange () { ; }

      ElementRange & operator= (const ElementRange & r2) = delete;
      
      INLINE ElementIterator begin () const 
      {
        ElementId ei = ElementId(vb,First());
        while ((ei.Nr() < IntRange::end()) && 
               (definedon.Size() && !definedon[fes.GetMeshAccess()->GetElIndex(ei)]))
          ++ei;
        return ElementIterator(fes, ei, definedon, temp_dnums, lh); 
      }

      INLINE ElementIterator end () const 
      {
        return ElementIterator(fes, ElementId(vb,Next()), definedon, temp_dnums, lh); 
      }
    };

    /*
    ElementRange Elements (VorB vb = VOL, int heapsize = 10000) const
    {
      return ElementRange (*this, vb, IntRange (0, ma->GetNE(vb)), LocalHeap(heapsize));
    }
    */

    ElementRange Elements (VorB vb = VOL, LocalHeap && lh = 10000) const
    {
      // cout << "C++ FESpace::Elements with lh rvalue, name = " << lh.name << endl;
      return ElementRange (*this, vb, IntRange (0, ma->GetNE(vb)), move(lh));
    }

    ElementRange Elements (VorB vb, LocalHeap & lh) const
    {
      return ElementRange (*this, vb, IntRange (0, ma->GetNE(vb)), lh);
    }
        

    /// returns finite element. 
    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const = 0;

    [[deprecated("Use GetFE with element-id instead of elnr!")]]    
    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const final;
    [[deprecated("Use GetFE(ElementId(BND,elnr)) instead!")]]    
    virtual const FiniteElement & GetSFE (int elnr, LocalHeap & lh) const final;
    [[deprecated("Use GetFE(ElementId(BBND,elnr)) instead!")]]        
    virtual const FiniteElement & GetCD2FE (int cd2elnr, LocalHeap & lh) const final;

    /// get dof-nrs of the element
    [[deprecated("Use GetDofNrs with element-id instead of elnr!")]]
    void GetDofNrs (int elnr, Array<DofId> & dnums) const
      { GetDofNrs(ElementId(VOL,elnr),dnums); }

    /// get dof-nrs of domain or boundary element elnr
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const = 0;
    virtual void GetDofNrs (NodeId ni, Array<DofId> & dnums) const;
    
    Table<int> CreateDofTable (VorB vorb) const;

    // FlatArray<int> GetDofNrs (ElementId ei, LocalHeap & lh) const;
    
    /// get coupling types of dofs
    virtual void GetDofCouplingTypes (int elnr, Array<COUPLING_TYPE> & dnums) const;
    
    /// get coupling types of dof
    virtual COUPLING_TYPE GetDofCouplingType (DofId dof) const;
    virtual void SetDofCouplingType (DofId dof, COUPLING_TYPE ct) const;
    
    void CheckCouplingTypes() const;

    /// get dof-nrs of the element of certain coupling type
    void GetDofNrs (int elnr, Array<DofId> & dnums, COUPLING_TYPE ctype) const;    

    /// get dofs on nr'th node of type nt.
    [[deprecated("Use GetDofNrs with NodeId instead of nt/nr")]]    
    virtual void GetNodeDofNrs (NODE_TYPE nt, int nr, Array<int> & dnums) const final;
    /// get number of low-order dofs for node of type nt
    // virtual int GetNLowOrderNodeDofs ( NODE_TYPE nt ) const;
    // { return lodofs_per_node[nt]; }

    /// get dofs on vertex vnr
    [[deprecated("Use GetDofNrs(NODE_TYPE(NT_VERTEX,nr) instead")]]
    virtual void GetVertexDofNrs (int vnr, Array<DofId> & dnums) const;
    /// get dofs on edge enr
    [[deprecated("Use GetDofNrs(NODE_TYPE(NT_EDGE,nr) instead")]]    
    virtual void GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const;
    /// get dofs on face fnr
    virtual void GetFaceDofNrs (int fanr, Array<DofId> & dnums) const;
    /// get dofs on element (=cell) elnr
    virtual void GetInnerDofNrs (int elnr, Array<DofId> & dnums) const;

    virtual bool UsesDGCoupling () const throw() { return dgjumps; };

    /// returns dofs of sourface element
    [[deprecated("Use GetDofNrs(ElementId(BND,elnr)) instead!")]]
    void GetSDofNrs (int selnr, Array<DofId> & dnums) const
      { GetDofNrs(ElementId(BND,selnr),dnums); }

    bool DefinedOn(VorB vb, int domnr) const
    { return !definedon[vb].Size() || definedon[vb][domnr]; }

    /// is the FESpace defined for this sub-domain nr ?
    [[deprecated("Use Definedon(VorB,int) instead")]]
    bool DefinedOn (int domnr) const
    { return !definedon[VOL].Size() || definedon[VOL][domnr]; }
    /// is the FESpace defined for this boundary nr ?
    [[deprecated("Use Definedon(VorB,int) instead")]]
    bool DefinedOnBoundary (int bnr) const
    {return !definedon[BND].Size() || definedon[BND][bnr]; }

    /// is the FESpace defined for this sub-domain / boundary nr ?
    [[deprecated("Use DefinedOn(VorB, int) instead")]]
    bool DefinedOn (int index, bool bound) const
    {
      if (bound)
        return !definedon[BND].Size() || definedon[BND][index];
      else
        return !definedon[VOL].Size() || definedon[VOL][index];
    }

    bool DefinedOn (ElementId id) const
    {
      if(!definedon[id.VB()].Size()) return true;
      return definedon[id.VB()][ma->GetElement(id).GetIndex()];
    }

    bool DefinedOn (Ngs_Element el) const
    {
      if(!definedon[el.VB()].Size()) return true;
      return definedon[el.VB()][el.GetIndex()];
    }

    void SetDefinedOn (VorB vb, const BitArray& defon);
    ///
    //[[deprecated("Use SetDefinedOn(VorB, const Bitarray&)")]]
     void SetDefinedOn (const BitArray & defon)
     { SetDefinedOn(VOL,defon); }
    ///
    //[[deprecated("Use SetDefinedOn(VorB, const Bitarray&)")]]
    void SetDefinedOnBoundary (const BitArray & defon)
     { SetDefinedOn(BND,defon); }

    ///
    void SetDirichletBoundaries (const BitArray & dirbnds);
    /// Get reference element for tet, prism, trig, etc ..
    const FiniteElement & GetFE (ELEMENT_TYPE type) const;

    /// according low-order FESpace (if available)
    FESpace & LowOrderFESpace () { return *low_order_space; }
    /// according low-order FESpace (if available)
    const FESpace & LowOrderFESpace () const { return *low_order_space; }
    shared_ptr<FESpace> LowOrderFESpacePtr () const { return low_order_space; }
    ///
    // void SetLowOrderSpace (bool los) { is_low_order_space = los; }
    ///
    // bool IsLowOrderSpace () const { return is_low_order_space; }

    /// non Dirichlet dofs
    virtual shared_ptr<BitArray> GetFreeDofs (bool external = false) const;
    bool IsFreeDof (DofId dof, bool external = false) const
    {
      if (external)
        return external_free_dofs->Test(dof);
      else
        return free_dofs->Test(dof);
    }
    ///
    bool IsDirichletDof (int i) const
    { return dirichlet_dofs.Size() && dirichlet_dofs[i]; }

    bool IsDirichletBoundary (int i) const
    { return dirichlet_boundaries.Size() && dirichlet_boundaries[i]; }

    /// is vertex on Dirichlet boundary ?
    bool IsDirichletVertex (size_t i) const { return dirichlet_vertex.Size() && dirichlet_vertex[i]; }
    /// is edge on Dirichlet boundary ?
    bool IsDirichletEdge (size_t i) const { return dirichlet_edge.Size() && dirichlet_edge[i]; }
    /// is face on Dirichlet boundary ?
    bool IsDirichletFace (size_t i) const { return dirichlet_face.Size() && dirichlet_face[i]; }

    void GetFilteredDofs(COUPLING_TYPE doffilter, BitArray & output, bool freedofsonly=true) const;
    /// 
    virtual shared_ptr<Table<int>> CreateSmoothingBlocks (const Flags & flags) const;
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

    bool IsAtomicDof (size_t nr) const { return (is_atomic_dof.Size() != 0) && is_atomic_dof[nr]; }
    bool HasAtomicDofs () const { return is_atomic_dof.Size() != 0; }

    [[deprecated("Use TransformMat with VorB  instead of bool")]]
    void TransformMat (int elnr, bool boundary,
		       const SliceMatrix<double> & mat, TRANSFORM_TYPE type) const
    {
      TransformMat(ElementId(boundary ? BND : VOL, elnr), mat, type);
    }
  
    [[deprecated("Use TransformMat with VorB  instead of bool")]]
    void TransformMat (int elnr, bool boundary,
		       const SliceMatrix<Complex> & mat, TRANSFORM_TYPE type) const
    {
      TransformMat(ElementId(boundary ? BND : VOL, elnr), mat, type);
    }
  
    [[deprecated("Use TransformVec with VorB  instead of bool")]]
    void TransformVec (int elnr, bool boundary,
		       const FlatVector<double> & vec, TRANSFORM_TYPE type) const
    {
      // VTransformVR (elnr, boundary ? BND : VOL, vec, type);
      VTransformVR (ElementId(boundary ? BND : VOL, elnr), vec, type);
    }
  
    [[deprecated("Use TransformVec with VorB  instead of bool")]]
    void TransformVec (int elnr, bool boundary,
		       const FlatVector<Complex> & vec, TRANSFORM_TYPE type) const
    {
      // VTransformVC (elnr, boundary ? BND : VOL, vec, type);
      VTransformVC (ElementId(boundary ? BND : VOL, elnr), vec, type);      
    }

    [[deprecated("Use TransformMat with VorB  instead of bool")]]    
    void TransformMat (int elnr, VorB vb,
                       const SliceMatrix<double> & mat, TRANSFORM_TYPE type) const
    {
      // VTransformMR (elnr, vb, mat, type);
      VTransformMR (ElementId(vb, elnr), mat, type);            
    }

    [[deprecated("Use TransformMat with VorB  instead of bool")]]    
    void TransformMat (int elnr, VorB vb,
		       const SliceMatrix<Complex> & mat, TRANSFORM_TYPE type) const
    {
      // VTransformMC (elnr, vb, mat, type);
      VTransformMC (ElementId(vb, elnr), mat, type);                  
    }

    [[deprecated("Use TransformVec with VorB  instead of bool")]]        
    void TransformVec (int elnr, VorB vb,
		       const FlatVector<double> & vec, TRANSFORM_TYPE type) const
    {
      // VTransformVR (elnr, vb, vec, type);
      VTransformVR (ElementId(vb, elnr), vec, type);            
    }

    [[deprecated("Use TransformVec with VorB  instead of bool")]]            
    void TransformVec (int elnr, VorB vb,
		       const FlatVector<Complex> & vec, TRANSFORM_TYPE type) const
    {
      // VTransformVC (elnr, vb, vec, type);
      VTransformVC (ElementId(vb, elnr), vec, type);                  
    }

    
    void TransformMat (ElementId ei, 
                       SliceMatrix<double> mat, TRANSFORM_TYPE type) const
    {
      VTransformMR (ei, mat, type);
    }
    void TransformMat (ElementId ei, 
		       SliceMatrix<Complex> mat, TRANSFORM_TYPE type) const
    {
      VTransformMC (ei, mat, type);
    }		
    void TransformVec (ElementId ei, 
		       SliceVector<double> vec, TRANSFORM_TYPE type) const
    {
      VTransformVR (ei, vec, type);
    }
    void TransformVec (ElementId ei, 
		       SliceVector<Complex> vec, TRANSFORM_TYPE type) const
    {
      VTransformVC (ei, vec, type);
    }

    
    template < int S, class T >
    void TransformVec (int elnr, VorB vb,
		       const FlatVector< Vec<S,T> >& vec, TRANSFORM_TYPE type) const;

    /*
    template < class T >
    void TransformVec (ElementId ei,
		       const T & vec, TRANSFORM_TYPE type) const
    {
      TransformVec (ei, vec, type);
    }
    */

    virtual void VTransformMR (ElementId ei,
			       const SliceMatrix<double> mat, TRANSFORM_TYPE type) const
    { ; }
    virtual void VTransformMC (ElementId ei, 
			       const SliceMatrix<Complex> mat, TRANSFORM_TYPE type) const
    { ; }


    virtual void VTransformVR (ElementId ei,
			       const SliceVector<double> vec, TRANSFORM_TYPE type) const
    { ; }
    virtual void VTransformVC (ElementId ei, 
			       const SliceVector<Complex> vec, TRANSFORM_TYPE type) const
    { ; }
  
  
    /// Returns multigrid-prolongation
    virtual shared_ptr<Prolongation> GetProlongation () const { return prol; }
    /// Set multigrid prolongation
    // void SetProlongation (ngmg::Prolongation * aprol)
    // { prol = aprol; }


    /// returns function-evaluator
    shared_ptr<DifferentialOperator> GetEvaluator (VorB vb = VOL) const
    {
      return evaluator[vb];
    }

    [[deprecated("Use GetEvaluator(VorB) instead of GetEvaluator(bool)!")]]
    shared_ptr<DifferentialOperator> GetEvaluator (bool boundary) const
    {
      if(boundary)
	return evaluator[BND];
      else
	return evaluator[VOL];
    }

    shared_ptr<DifferentialOperator> GetFluxEvaluator (VorB vb=VOL) const
    {
      return flux_evaluator[vb];
    }

    [[deprecated("Use GetFluxEvaluator(VorB) instead of GetFluxEvaluator(bool)!")]]
    shared_ptr<DifferentialOperator> GetFluxEvaluator (bool boundary) const
    {
      if(boundary)
	return flux_evaluator[BND];
      else
	return flux_evaluator[VOL];
    }

    virtual SymbolTable<shared_ptr<DifferentialOperator>> GetAdditionalEvaluators () const
    { return additional_evaluators; } 

    /// returns function-evaluator
    [[deprecated("Use GetIntegrator(VorB) instead of GetIntegrator(bool)!")]]    
    shared_ptr<BilinearFormIntegrator> GetIntegrator (bool vb = VOL) const
    {
      return integrator[vb];
    }

    shared_ptr<BilinearFormIntegrator> GetIntegrator (VorB vb = VOL) const
    {
      return integrator[vb];
    }

    /// special elements for hacks (used for contact, periodic-boundary-penalty-constraints, ...
    Array<SpecialElement*> specialelements;

    void AppendSpecialElement (SpecialElement * spel)
    { specialelements.Append (spel); }

    const Array<SpecialElement*> & GetSpecialElements() const {return specialelements;}

    virtual void SolveM(CoefficientFunction & rho, BaseVector & vec,
                        LocalHeap & lh) const;
      
    ParallelDofs & GetParallelDofs () const { return *paralleldofs; }
    virtual void UpdateParallelDofs ();

    //// is FESpace mpi-distributed ?
    bool IsParallel() const;

    /// ndof over all mpi-partitions
    size_t GetNDofGlobal() const;

    virtual int GetRelOrder() const
    { 
      cout << "virtual GetRelOrder called for FiniteElementSpace, not available ! " << endl; 
      return 0; 
    } 

    virtual bool VarOrder() const { return 0; }

    bool timing;
    std::list<std::tuple<std::string,double>> Timing () const;

  protected:
    template <template <ELEMENT_TYPE ET> class FE>
    void SetDummyFE ()
    {
      delete dummy_tet;
      delete dummy_pyramid;
      delete dummy_prism;
      delete dummy_hex;
      delete dummy_trig;
      delete dummy_quad;
      delete dummy_segm;
      delete dummy_point;
      dummy_tet = new FE<ET_TET>();
      dummy_pyramid = new FE<ET_PYRAMID>();
      dummy_prism = new FE<ET_PRISM>();
      dummy_hex = new FE<ET_HEX>();
      dummy_trig = new FE<ET_TRIG>();
      dummy_quad = new FE<ET_QUAD>();
      dummy_segm = new FE<ET_SEGM>();
      dummy_point = new FE<ET_POINT>();
    }
  };



  extern NGS_DLL_HEADER void IterateElements (const FESpace & fes,
			       VorB vb, 
			       LocalHeap & clh, 
			       const function<void(FESpace::Element,LocalHeap&)> & func);
  /*
  template <typename TFUNC>
  inline void IterateElements (const FESpace & fes, 
                               VorB vb, 
                               LocalHeap & clh, 
                               const TFUNC & func)
  {
    IterateElements1 (fes, vb, clh, func);
  }
  */
  
  /*
  template <typename TFUNC>
  inline void IterateElements (const FESpace & fes, 
                               VorB vb, 
                               LocalHeap & clh, 
                               const TFUNC & func)
  {
    
#pragma omp parallel 
    {

#pragma omp single
      {
        const Table<int> & element_coloring = fes.ElementColoring(vb);

        for (FlatArray<int> els_of_col : element_coloring)
          {

            for (int i = 0; i < els_of_col.Size(); i++)
              {
#pragma omp task
                {
                  LocalHeap lh = clh.Split();
                  Array<int> temp_dnums;
                  FESpace::Element el(fes, ElementId (vb, els_of_col[i]), temp_dnums);
                  func (el, lh);
                }
              }

#pragma omp taskwait
          }
      }

    }

  }
  */









#ifdef OLD_REMOVED_FOR_CLANG
  template <typename TFUNC>
  inline void IterateElementsInsideParallel (const FESpace & fes, 
                                             VorB vb, 
                                             LocalHeap & lh, 
                                             const TFUNC & func)
  {
    const Table<int> & element_coloring = fes.ElementColoring(vb);
    
    Array<int> temp_dnums;

    // lh.ClearValues();
    
    for (FlatArray<int> els_of_col : element_coloring)
      
#pragma omp for schedule(dynamic)
      for (int i = 0; i < els_of_col.Size(); i++)
        {
          HeapReset hr(lh);
          FESpace::Element el(fes, ElementId (vb, els_of_col[i]), temp_dnums);
          func (el, lh);
        }
    // cout << "lh, used size = " << lh.UsedSize() << endl;
  }
#endif





  /**
     A space of continuous finite elements.
     Supports first and second order finite elements.
  */
  class NGS_DLL_HEADER NodalFESpace : public FESpace
  {
    ///
    // Array<int> ndlevel;
    bool hb_defined;

  public:

    ///
    NodalFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags=false);
    ///
    virtual ~NodalFESpace ();

    ///
    virtual string GetClassName () const override
    {
      return "NodalFESpace";
    }

    ///
    virtual void Update (LocalHeap & lh) override;
    
    virtual void DoArchive (Archive & archive) override;

    virtual FiniteElement & GetFE(ElementId ei, Allocator & lh) const override;
    ///
    // virtual size_t GetNDof () const throw() override;
    ///
    // virtual size_t GetNDofLevel (int level) const override;
    ///
    // using FESpace::GetDofNrs;
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    ///

    virtual void GetVertexDofNrs (int vnr, Array<DofId> & dnums) const override;
    virtual void GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const override;
    virtual void GetFaceDofNrs (int fanr, Array<DofId> & dnums) const override;
    virtual void GetInnerDofNrs (int elnr, Array<DofId> & dnums) const override;

    virtual Array<int> * CreateDirectSolverClusters (const Flags & flags) const override;
  };






  ///
  class NGS_DLL_HEADER NonconformingFESpace : public FESpace
  {
    ///
    Array<int> ndlevel;

  public:
    NonconformingFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags=false);
    virtual ~NonconformingFESpace ();

    virtual string GetClassName () const override
    { return "Nonconforming FESpace"; }

    ///
    virtual void Update(LocalHeap & lh) override;

    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;
    ///
    virtual size_t GetNDof () const throw() override;
    ///
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
  };







  ///
  class NGS_DLL_HEADER ElementFESpace : public FESpace
  {
    ///  Array<int> startelement;
    // Array<int> ndlevel;
    int n_el_dofs;
  public:
    ///
    ElementFESpace (shared_ptr<MeshAccess> ama, const Flags& flags, bool parseflags=false);

    ///
    ~ElementFESpace ();

    virtual string GetClassName () const override
    {
      return "ElementFESpace";
    }

    ///
    virtual void Update(LocalHeap & lh) override;
    /// 
    virtual void DoArchive (Archive & archive) override;

    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;
    ///
    // virtual size_t GetNDof () const throw() override { return ndlevel.Last(); }
  
    ///
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;

    ///
    // virtual size_t GetNDofLevel (int level) const override;


    virtual void GetVertexDofNrs (int vnr, Array<DofId> & dnums) const override
    { dnums.SetSize (0); }
    virtual void GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const override
    { dnums.SetSize (0); }
    virtual void GetFaceDofNrs (int fanr, Array<DofId> & dnums) const override
    { dnums.SetSize (0); }
    virtual void GetInnerDofNrs (int elnr, Array<DofId> & dnums) const override
    { GetDofNrs (elnr, dnums); }
  };





  /// Non-continous fe space on boundary
  class NGS_DLL_HEADER SurfaceElementFESpace : public FESpace
  {
    ///
    Array<int> ndlevel;
    int n_el_dofs;
  public:
    ///
    SurfaceElementFESpace (shared_ptr<MeshAccess> ama, const Flags& flags, 
                           bool checkflags = false);

    ///
    ~SurfaceElementFESpace ();

    ///
    virtual string GetClassName() const
    { return "SurfaceElement"; }

    ///
    virtual void Update(LocalHeap & lh);

    ///
    virtual size_t GetNDof () const throw() { return ndlevel.Last(); } 

    ///
    virtual const FiniteElement & GetFE (ElementId ei, LocalHeap & lh) const;

    ///
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const;

    ///
    virtual size_t GetNDofLevel (int level) const;

  };




  


  /// A combination of fe-spaces
  class NGS_DLL_HEADER CompoundFESpace : public FESpace
  {
  protected:
    /// pointers to components
    Array<shared_ptr<FESpace>> spaces;
    /// cummlated number of dofs of components
    Array<int> cummulative_nd;
    /// dofs on each multigrid level
    /// Array<int> ndlevel;
  public:
    /// generates a compound space.
    /// components will be added later
    CompoundFESpace (shared_ptr<MeshAccess> ama,
		     const Flags & flags, bool parseflags = false);
    /// generates a compound space 
    /// components are provided in aspaces
    CompoundFESpace (shared_ptr<MeshAccess> ama,
		     const Array<shared_ptr<FESpace>> & aspaces,
		     const Flags & flags, bool parseflags = false);

    /// not much to do.
    /// components will not be deleted
    virtual ~CompoundFESpace ();

    /// add an additional component space
    void AddSpace (shared_ptr<FESpace> fes);

    ///
    virtual string GetClassName () const
    {
      return "CompoundFESpace";
    }

    /// updates also components
    virtual void Update(LocalHeap & lh);
    /// updates also components
    virtual void FinalizeUpdate(LocalHeap & lh);

    /// copies dofcoupling from components
    virtual void UpdateCouplingDofArray();

    /// 
    // virtual size_t GetNDof () const throw() { return cummulative_nd.Last(); } 
    ///
    // virtual size_t GetNDofLevel (int level) const { return ndlevel[level]; }

    IntRange GetRange (int spacenr) const
    { 
      return IntRange(cummulative_nd[spacenr], cummulative_nd[spacenr+1]);
    }

    /// get component space
    shared_ptr<FESpace> operator[] (int i) const { return spaces[i]; }

    /// returns a compound finite element
    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const;
    ///
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const;
    virtual void GetDofNrs (NodeId ni, Array<DofId> & dnums) const;
    ///
    [[deprecated("Use GetDofNrs(NODE_TYPE(NT_VERTEX,nr) instead")]]    
    virtual void GetVertexDofNrs (int vnr, Array<DofId> & dnums) const;
    [[deprecated("Use GetDofNrs(NODE_TYPE(NT_EDGE,nr) instead")]]    
    virtual void GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const;
    virtual void GetFaceDofNrs (int fanr, Array<DofId> & dnums) const;
    virtual void GetInnerDofNrs (int elnr, Array<DofId> & dnums) const;

    
    template <class T> NGS_DLL_HEADER
      void T_TransformMat (ElementId ei, 
                           SliceMatrix<T> mat, TRANSFORM_TYPE tt) const;
    
    template <class T> NGS_DLL_HEADER
      void T_TransformVec (ElementId ei, 
                         SliceVector<T> vec, TRANSFORM_TYPE tt) const;

    virtual void VTransformMR (ElementId ei,
			       SliceMatrix<double> mat, TRANSFORM_TYPE tt) const;
    virtual void VTransformMC (ElementId ei,
                               SliceMatrix<Complex> mat, TRANSFORM_TYPE tt) const;
    virtual void VTransformVR (ElementId ei,
                               SliceVector<double> vec, TRANSFORM_TYPE tt) const;
    virtual void VTransformVC (ElementId ei, 
                               SliceVector<Complex> vec, TRANSFORM_TYPE tt) const;

    /// number of component spaces
    inline int GetNSpaces () const { return spaces.Size(); } 
  };





  /// Registered FESpace classes
  class NGS_DLL_HEADER FESpaceClasses
  {
  public:
    /// descriptor for register fespaces. 
    /// function pointer to create function.
    struct FESpaceInfo
    {
      /// the name
      string name;
      /// function pointer to creator function
      shared_ptr<FESpace> (*creator)(shared_ptr<MeshAccess> ma, const Flags & flags);
      /// creates a descriptor
      FESpaceInfo (const string & aname,
		   shared_ptr<FESpace> (*acreator)(shared_ptr<MeshAccess> ma, const Flags & flags))
	: name(aname), creator(acreator) {;}
    };
  private:
    Array<shared_ptr<FESpaceInfo>> fesa;

  public:
    /// initialize 
    FESpaceClasses() { ; }
    /// cleans up
    ~FESpaceClasses();  

    /// add a descriptor
    void AddFESpace (const string & aname, 
		     shared_ptr<FESpace> (*acreator)(shared_ptr<MeshAccess> ma, const Flags & flags));
  
    /// returns all creators
    const Array<shared_ptr<FESpaceInfo>> & GetFESpaces() { return fesa; }

    /// returns a creator structure
    const shared_ptr<FESpaceInfo> GetFESpace(const string & name);

    /// print available fespaces to stream
    void Print (ostream & ost) const;
  };
 
  /// returns createion object
  extern NGS_DLL_HEADER FESpaceClasses & GetFESpaceClasses ();

  /// creates a fespace of that type
  extern NGS_DLL_HEADER shared_ptr<FESpace> CreateFESpace (const string & type,
                                                           shared_ptr<MeshAccess> ma,
                                                           const Flags & flags);


  /**
     template for registration of finite element spaces.
     provides static Create - function
   */
  template <typename FES>
  class RegisterFESpace
  {
  public:
    /// constructor registers fespace
    RegisterFESpace (string label)
    {
      GetFESpaceClasses().AddFESpace (label, Create);
      // cout << "register fespace '" << label << "'" << endl;
    }
    
    /// creates an fespace of type FES
    static shared_ptr<FESpace> Create (shared_ptr<MeshAccess> ma, const Flags & flags)
    {
      return make_shared<FES> (ma, flags);
    }
  };















#ifdef PARALLEL

  class ParallelMeshDofs : public ParallelDofs
  {
    shared_ptr<MeshAccess> ma;
    Array<Node> dofnodes;
  public:
    ParallelMeshDofs (shared_ptr<MeshAccess> ama, const Array<Node> & adofnodes, 
		      int dim = 1, bool iscomplex = false);

    shared_ptr<MeshAccess> GetMeshAccess() const { return ma; }
    const Array<Node> & GetDofNodes() const { return dofnodes; }
  };
  
#else


  class ParallelMeshDofs : public ParallelDofs 
  {
  public:
    ParallelMeshDofs (shared_ptr<MeshAccess> ama, const Array<NodeId> & adofnodes, 
		      int dim = 1, bool iscomplex = false)
    { ndof = adofnodes.Size(); }
  };

#endif

}





#ifdef PARALLEL
namespace ngstd
{
  template<>
  class MPI_Traits<ngcomp::COUPLING_TYPE>
  {
  public:
    /// returns MPI-type 
    static MPI_Datatype MPIType () 
    { 
      if (sizeof(ngcomp::COUPLING_TYPE) == sizeof(int)) return MPI_INT;
      cout << "please provide MPI_Datatype for COUPLING_TYPE" << endl;
      exit(1);
    }
  };
}
#endif

#endif
