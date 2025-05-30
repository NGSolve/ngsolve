#ifndef FILE_FESPACE
#define FILE_FESPACE

/*********************************************************************/
/* File:   fespace.hpp                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

#include <core/register_archive.hpp>

#include <finiteelement.hpp>
#include <diffop.hpp>
#include <symbolicintegrator.hpp>   // for ProxyFunction

#include <basevector.hpp>
#include <basematrix.hpp>

// #include <paralleldofs.hpp>


#include "ngsobject.hpp"

namespace ngmg
{
  class Prolongation;
}


namespace ngcomp
{
  using namespace ngla;
  
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
			TRANSFORM_SOL = 8,
                        TRANSFORM_SOL_INVERSE = 16};
  /**
    coupling types: Each degree of freedom is either
     - an unused or hidden dof (invisible)
     - a local degree of freedom 
     - an interface degree of freedom 
     or
     - a wirebasket degree of freedom
  */
  enum COUPLING_TYPE : uint8_t {  UNUSED_DOF = 0,
			HIDDEN_DOF = 1,
			LOCAL_DOF = 2,
			CONDENSABLE_DOF = 3,
			INTERFACE_DOF = 4,
			NONWIREBASKET_DOF = 6,
			WIREBASKET_DOF = 8,
			EXTERNAL_DOF = 12,
			VISIBLE_DOF = 14,
			ANY_DOF = 15
		      };
/*
Bit encoding:          
UNUSED                     0 |  0
HIDDEN                     1 |  1
LOCAL                    1 0 |  2
CONDENSABLE              1 1 |  3
INTERFACE              1 0 0 |  4
NONWIREBASKET          1 1 0 |  6
WIREBASKET           1 0 0 0 |  8
EXTERNAL             1 1 0 0 | 12
VISIBLE              1 1 1 0 | 14
ANY                  1 1 1 1 | 15
*/

  /**
     constant_order .... one order for everything
     node_type_order ... order for edges, or faces, gradients or curls, but all over the mesh
     variable_order .... a different, anisotropic order for every mesh node
   */
  enum ORDER_POLICY { CONSTANT_ORDER = 0, NODE_TYPE_ORDER = 1, VARIABLE_ORDER = 2, OLDSTYLE_ORDER = 3 };

  
  
  NGS_DLL_HEADER ostream & operator<< (ostream & ost, COUPLING_TYPE ct);


  class FESpace;
  class BilinearForm;

  // will be size_t some day 
  typedef int DofId;
  enum IRREGULAR_DOF_NR
    {
      NO_DOF_NR = -1,            // don't assemble this dof (it has no regular number)
      NO_DOF_NR_CONDENSE = -2    // condense out this dof, don't assemble to global system
    };
  INLINE bool IsRegularDof (DofId dof) { return dof >= 0; } // ATTENTION for size_t 

  
  using ngmg::Prolongation;


  struct ProxyNode : public shared_ptr<ProxyFunction>
  {
    // shared_ptr<ProxyFunction> proxy;
    std::vector<ProxyNode> list;
    
    ProxyNode (shared_ptr<ProxyFunction> _proxy) : shared_ptr<ProxyFunction>(_proxy) { }
    ProxyNode (std::vector<ProxyNode> _list) : list(_list) { }
    auto operator* () const { return shared_ptr<ProxyFunction> (*this); }
    auto operator[] (int i) const { return list[i]; }

    void SetFESpace (shared_ptr<FESpace> fespace)
    {
      if (*this)
        shared_ptr<ProxyFunction> (*this) -> SetFESpace(fespace);
      else
        for (auto node : list)
          node.SetFESpace(fespace);
    }
  };



  
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

    /// whether the space should update itself on changes to the mesh
    bool autoupdate;

    /// debug output to testout
    bool print; 

    /// prolongation operators between multigrid levels
    shared_ptr<Prolongation> prol;// = NULL;
    /// highest multigrid-level for which Update was called (memory allocation)
    int level_updated;

    /// on which subdomains is the space defined ?
    Array<bool> definedon[4];

    /// prototype: what are the Dirichlet boundaries ?
    BitArray dirichlet_constraints[4];
    BitArray & dirichlet_boundaries = dirichlet_constraints[1];

    /// dofs on Dirichlet boundary
    BitArray dirichlet_dofs;
    shared_ptr<BitArray> free_dofs;
    shared_ptr<BitArray> external_free_dofs;


    Array<bool> dirichlet_vertex;
    Array<bool> dirichlet_edge;
    Array<bool> dirichlet_face;
    
    
    /// Evaluator for visualization (new style)
    shared_ptr<DifferentialOperator> evaluator[4];
    /// Evaluator for flux
    shared_ptr<DifferentialOperator> flux_evaluator[4];

    SymbolTable<shared_ptr<DifferentialOperator>> additional_evaluators;

    /// Evaluator for visualization (old style)
    shared_ptr<BilinearFormIntegrator> integrator[4];

    /// if non-zero, pointer to low order space
    shared_ptr<FESpace> low_order_space; 
    shared_ptr<BaseMatrix> low_order_embedding;
      
    /// if directsolverclustered[i] is true, then the unknowns of domain i are clustered
    Array<bool> directsolverclustered;

    Array<string> directsolvermaterials;

    mutable Array<int> adddirectsolverdofs;

    Array<int> directvertexclusters;
    Array<int> directedgeclusters;
    Array<int> directfaceclusters;
    Array<int> directelementclusters;

    
    Table<int> element_coloring[4]; 
    Table<int> facet_coloring;  // elements on facet in own colors (DG)
    Array<COUPLING_TYPE> ctofdof;

    shared_ptr<ParallelDofs> paralleldofs;

    bool no_low_order_space;

    int et_bonus_order[30]; // order increase for element-type

    typedef int8_t TORDER;

    ORDER_POLICY order_policy = OLDSTYLE_ORDER;
  
    // size_t order_timestamp = 0;
    BitArray is_atomic_dof;

    // only a few spaces (lowest order Nedelec) need the transformation
    // of element vectors
    bool needs_transform_vec = true;

    
    // move ndof and ndof_level to FESpace base class
    size_t ndof;
    Array<size_t> ndof_level;
  protected:
    void SetNDof (size_t _ndof);
    
  public:
    string type;
    SimpleSignal updateSignal;
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
    FESpace (const FESpace &) = delete;
    /// cleanup
    virtual ~FESpace ();

    static DocInfo GetDocu ();
    
    /// update dof-table
    virtual void Update();

    virtual void UpdateDofTables() { ; } 
    virtual void UpdateCouplingDofArray() { ; } 
    virtual void UpdateFreeDofs();
    /// update element coloring
    virtual void FinalizeUpdate();

    /// highest level where update/finalize was called
    int GetLevelUpdated() const { return level_updated; }

    const Table<int> & ElementColoring(VorB vb = VOL) const 
    { return element_coloring[vb]; }

    const Table<int> & FacetColoring() const;
    
    /// print report to stream
    virtual void PrintReport (ostream & ost) const override;

    /// Dump/restore fespace
    virtual void DoArchive (Archive & archive);
    std::tuple<Shallow<shared_ptr<MeshAccess>>, Flags> GetCArgs();

    Array<MemoryUsage> GetMemoryUsage () const override;
    
    /// order of finite elements
    int GetOrder () const { return order; }

    /*
    void SetBonusOrder (ELEMENT_TYPE et, int bonus) 
    { et_bonus_order[et] = bonus; }
    */
    void SetOrderPolicy (ORDER_POLICY op)
    {
      order_policy = op;
    }
    
    virtual void SetOrder (ELEMENT_TYPE et, TORDER order);

    virtual void SetOrder (NodeId ni, int order); 
    virtual int GetOrder (NodeId ni) const; 

    /// how many components
    int GetDimension () const { return dimension; }

    /// complex space ?
    bool IsComplex () const { return iscomplex; }

    virtual int GetSpatialDimension() const { return ma->GetDimension();}

    /// number of (process-local) dofs
    virtual size_t GetNDof () const { return ndof; } 
    /// number of dofs on the level
    virtual size_t GetNDofLevel (int level) const
    { return (level < ndof_level.Size()) ? ndof_level[level] : GetNDof(); } 

    virtual FlatArray<VorB> GetDualShapeNodes (VorB vb) const;
    
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
      auto & GetFESpace() const { return fes; } 
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
          vb(avb), mylh(std::move(lh2)), lh(mylh)
      { ; }

      INLINE ElementRange (const FESpace & afes, VorB avb, IntRange ar, LocalHeap & lh2) 
        : IntRange(ar), fes(afes), 
          definedon(fes.definedon[avb].Size(),fes.definedon[avb].Addr(0)),
          vb(avb), mylh(), lh(lh2)
      { ; }

      ElementRange (const ElementRange & r2) = delete;

      INLINE ElementRange (ElementRange && r2) 
        : IntRange(r2), fes(r2.fes), definedon(std::move(r2.definedon)), vb(r2.vb), 
          temp_dnums(std::move(r2.temp_dnums)), mylh(std::move(r2.mylh)), 
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

    ElementRange Elements (VorB vb = VOL, LocalHeap && lh = 10000) const
    {
      // cout << "C++ FESpace::Elements with lh rvalue, name = " << lh.name << endl;
      return ElementRange (*this, vb, IntRange (0, ma->GetNE(vb)), std::move(lh));
    }

    ElementRange Elements (VorB vb, LocalHeap & lh) const
    {
      return ElementRange (*this, vb, IntRange (0, ma->GetNE(vb)), lh);
    }

    /// returns finite element. 
    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const = 0;


    /// get dof-nrs of domain or boundary element elnr
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const = 0;
    
    virtual void GetDofNrs (NodeId ni, Array<DofId> & dnums) const;
    BitArray GetDofs (const Region & reg) const;
    Table<int> CreateDofTable (VorB vorb) const;

    /// get coupling types of dofs
    virtual void GetDofCouplingTypes (int elnr, Array<COUPLING_TYPE> & dnums) const;
    
    /// get coupling types of dof
    // virtual COUPLING_TYPE GetDofCouplingType (DofId dof) const;
    // make sure we have it, otherwise throw exception
    bool CouplingTypeArrayAvailable() const { return ctofdof.Size() == GetNDof(); }
    COUPLING_TYPE GetDofCouplingType (DofId dof) const
    { return IsRegularDof(dof)
        ? ( (ctofdof.Size()==0) ? WIREBASKET_DOF : ctofdof[dof])  // would like to rely on the ctarray
        : ( (dof == NO_DOF_NR) ? UNUSED_DOF : HIDDEN_DOF ); }
    
    virtual void SetDofCouplingType (DofId dof, COUPLING_TYPE ct) const;
    auto & CouplingTypes() { return ctofdof; }
    void CheckCouplingTypes() const;
      
    /// get dof-nrs of the element of certain coupling type
    void GetDofNrs (ElementId ei, Array<DofId> & dnums, COUPLING_TYPE ctype) const;

    /// get dofs (local numbering) of a certain type
    virtual void GetElementDofsOfType (ElementId ei, Array<DofId> & dnums, COUPLING_TYPE ctype) const;



    /// get dofs on vertex vnr
    // [[deprecated("Use GetDofNrs(NODE_TYPE(NT_VERTEX,nr) instead")]]
    virtual void GetVertexDofNrs (int vnr, Array<DofId> & dnums) const;
    /// get dofs on edge enr
    // [[deprecated("Use GetDofNrs(NODE_TYPE(NT_EDGE,nr) instead")]]    
    virtual void GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const;
    /// get dofs on face fnr
    virtual void GetFaceDofNrs (int fanr, Array<DofId> & dnums) const;
    /// get dofs on element (=cell) elnr
    virtual void GetInnerDofNrs (int elnr, Array<DofId> & dnums) const;
    /// get dofs that are globally defined
    virtual void GetGlobalDofNrs (int gnr, Array<DofId> & dnums) const;

    virtual bool UsesDGCoupling () const throw() { return dgjumps; };

    bool DoesAutoUpdate () const { return autoupdate; };
    void ConnectAutoUpdate();
    

    auto & DefinedOn(VorB vb) const { return definedon[vb]; }
    
    bool DefinedOn(VorB vb, int domnr) const
    { return !definedon[vb].Size() || definedon[vb][domnr]; }


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

    xbool DefinedOnX (Ngs_Element el) const
    {
      // a temporary workaround,
      // clean solution will be to set definedon[BND] correctly
      if (el.VB() <= BND) return DefinedOn(el);
      
      if (!definedon[el.VB()].Size()) return maybe;
      return definedon[el.VB()][el.GetIndex()];
    }

    
    virtual void SetDefinedOn (VorB vb, const BitArray& defon);
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
    // const FiniteElement & GetFE (ELEMENT_TYPE type) const;

    /// according low-order FESpace (if available)
    [[deprecated("Use LowOrderFESpacePtr instead!")]]            
    FESpace & LowOrderFESpace () { return *low_order_space; }
    /// according low-order FESpace (if available)
    [[deprecated("Use LowOrderFESpacePtr instead!")]]                
    const FESpace & LowOrderFESpace () const { return *low_order_space; }
    shared_ptr<FESpace> LowOrderFESpacePtr () const { return low_order_space; }
    shared_ptr<BaseMatrix> LowOrderEmbedding () const { return low_order_embedding; }
    
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
    virtual shared_ptr<Array<int>> CreateDirectSolverClusters (const Flags & flags) const
    { return nullptr; }

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


    bool NeedsTransformVec() const { return needs_transform_vec; }

    void TransformMat (ElementId ei, 
                       SliceMatrix<double> mat, TRANSFORM_TYPE type) const
    {
      if (needs_transform_vec)      
        VTransformMR (ei, mat, type);
    }
    void TransformMat (ElementId ei, 
		       SliceMatrix<Complex> mat, TRANSFORM_TYPE type) const
    {
      if (needs_transform_vec)      
        VTransformMC (ei, mat, type);
    }		
    void TransformVec (ElementId ei, 
		       SliceVector<double> vec, TRANSFORM_TYPE type) const
    {
      if (needs_transform_vec)
        VTransformVR (ei, vec, type);
    }
    void TransformVec (ElementId ei, 
		       SliceVector<Complex> vec, TRANSFORM_TYPE type) const
    {
      if (needs_transform_vec)
        VTransformVC (ei, vec, type);
    }

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
    virtual shared_ptr<Prolongation> GetProlongation () const;
    // { return prol; }
    /// Set multigrid prolongation
    // void SetProlongation (ngmg::Prolongation * aprol)
    // { prol = aprol; }
    virtual void SetHarmonicProlongation (shared_ptr<BilinearForm> bfa, string inverse);

    /// returns function-evaluator
    shared_ptr<DifferentialOperator> GetEvaluator (VorB vb = VOL) const
    {
      return evaluator[vb];
    }


    shared_ptr<DifferentialOperator> GetFluxEvaluator (VorB vb = VOL) const
    {
      return flux_evaluator[vb];
    }


    virtual SymbolTable<shared_ptr<DifferentialOperator>> GetAdditionalEvaluators () const;
      // { return additional_evaluators; } 


    shared_ptr<BilinearFormIntegrator> GetIntegrator (VorB vb = VOL) const;
    /*
    {
      return integrator[vb];
    }
    */

    ProxyNode GetProxyFunction(bool testfunction) const;
    virtual ProxyNode MakeProxyFunction (bool testfunction,
                                         const function<shared_ptr<ProxyFunction>(shared_ptr<ProxyFunction>)> & addblock) const;
    
    auto GetTrialFunction() const { return GetProxyFunction(false); }
    auto GetTestFunction() const { return GetProxyFunction(true); }
    

    virtual shared_ptr<BaseMatrix> GetMassOperator (shared_ptr<CoefficientFunction> rho,
                                                    shared_ptr<Region> defon,
                                                    LocalHeap & lh) const;

    
    virtual shared_ptr<BaseMatrix> CreateMassOperator (shared_ptr<CoefficientFunction> rho,
                                                       shared_ptr<Region> defon,
                                                       bool inverse,
                                                       LocalHeap & lh) const;
    
    virtual void SolveM(CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                        LocalHeap & lh) const;
    virtual void ApplyM(CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                        LocalHeap & lh) const;

    virtual shared_ptr<BaseMatrix> GetTraceOperator (shared_ptr<FESpace> tracespace, bool avg) const;

    virtual shared_ptr<BaseMatrix> ConvertL2Operator (shared_ptr<FESpace> l2space) const;
    
    virtual void GetTrace (const FESpace & tracespace, const BaseVector & in, BaseVector & out, bool avg,
                           LocalHeap & lh) const;
    
    virtual void GetTraceTrans (const FESpace & tracespace, const BaseVector & in, BaseVector & out, bool avg,
                                LocalHeap & lh) const;
    
    shared_ptr<ParallelDofs> GetParallelDofs () const { return paralleldofs; }
    virtual void UpdateParallelDofs ();

    //// is FESpace mpi-distributed ?
    bool IsParallel() const;

    /// ndof over all mpi-partitions
    size_t GetNDofGlobal() const;

    virtual int GetRelOrder() const
    { 
      return 0; 
    } 

    virtual bool VarOrder() const { return 0; }

    bool timing;
    std::list<std::tuple<std::string,double>> Timing () const;




      /*
    [[deprecated("Use GetFE with element-id instead of elnr!")]]    
    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const final;
    [[deprecated("Use GetFE(ElementId(BND,elnr)) instead!")]]    
    virtual const FiniteElement & GetSFE (int elnr, LocalHeap & lh) const final;
    [[deprecated("Use GetFE(ElementId(BBND,elnr)) instead!")]]        
    virtual const FiniteElement & GetCD2FE (int cd2elnr, LocalHeap & lh) const final;
*/
    /// get dof-nrs of the element
    [[deprecated("Use GetDofNrs with element-id instead of elnr!")]]
    void GetDofNrs (int elnr, Array<DofId> & dnums) const
      { GetDofNrs(ElementId(VOL,elnr),dnums); }

    [[deprecated("Use GetDofNrs with element-id instead of elnr!")]]
    void GetDofNrs (int elnr, Array<DofId> & dnums, COUPLING_TYPE ctype) const;

    /// get dofs on nr'th node of type nt.
    [[deprecated("Use GetDofNrs with NodeId instead of nt/nr")]]    
    virtual void GetNodeDofNrs (NODE_TYPE nt, int nr, Array<int> & dnums) const final;
    /// get number of low-order dofs for node of type nt
    // virtual int GetNLowOrderNodeDofs ( NODE_TYPE nt ) const;
    // { return lodofs_per_node[nt]; }

    /// returns dofs of sourface element
    [[deprecated("Use GetDofNrs(ElementId(BND,elnr)) instead!")]]
    void GetSDofNrs (int selnr, Array<DofId> & dnums) const
      { GetDofNrs(ElementId(BND,selnr),dnums); }

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

    [[deprecated("Use GetEvaluator(VorB) instead of GetEvaluator(bool)!")]]
    shared_ptr<DifferentialOperator> GetEvaluator (bool boundary) const
    {
      if(boundary)
	return evaluator[BND];
      else
	return evaluator[VOL];
    }

    [[deprecated("Use GetFluxEvaluator(VorB) instead of GetFluxEvaluator(bool)!")]]
    shared_ptr<DifferentialOperator> GetFluxEvaluator (bool boundary) const
    {
      if(boundary)
	return flux_evaluator[BND];
      else
	return flux_evaluator[VOL];
    }

    /// returns function-evaluator
    [[deprecated("Use GetIntegrator(VorB) instead of GetIntegrator(bool)!")]]    
    shared_ptr<BilinearFormIntegrator> GetIntegrator (bool vb = VOL) const
    {
      return integrator[vb];
    }
  };



  extern NGS_DLL_HEADER void IterateElements (const FESpace & fes,
			       VorB vb, 
			       LocalHeap & clh, 
			       const function<void(FESpace::Element,LocalHeap&)> & func);
 


  /**
     A space of continuous finite elements.
     Supports first and second order finite elements.
  */
  class NGS_DLL_HEADER NodalFESpace : public FESpace
  {
    ///
    // Array<int> ndlevel;
    bool hb_defined;
    Array<bool> used_vertex;
    Array<bool> used_edge;

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
    void Update () override;
    void UpdateCouplingDofArray() override;
    
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

    virtual shared_ptr<Array<int>> CreateDirectSolverClusters (const Flags & flags) const override;
  };






  ///
  class NGS_DLL_HEADER NonconformingFESpace : public FESpace
  {
    ///
    // Array<int> ndlevel;

  public:
    NonconformingFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags=false);
    virtual ~NonconformingFESpace ();

    virtual string GetClassName () const override
    { return "Nonconforming FESpace"; }

    ///
    void Update() override;

    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;
    ///
    // virtual size_t GetNDof () const throw() override;
    ///
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
  };



  ///
  class NGS_DLL_HEADER NonconformingSurfaceFESpace : public FESpace
  {
    ///
    // Array<int> ndlevel;

  public:
    NonconformingSurfaceFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags=false);
    virtual ~NonconformingSurfaceFESpace ();

    virtual string GetClassName () const override
    { return "Nonconforming surface FESpace"; }

    ///
    void Update() override;

    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;
    ///
    // virtual size_t GetNDof () const throw() override;
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
    void Update() override;
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
    // Array<int> ndlevel;
    int n_el_dofs;
  public:
    ///
    SurfaceElementFESpace (shared_ptr<MeshAccess> ama, const Flags& flags, 
                           bool checkflags = false);

    ///
    ~SurfaceElementFESpace ();

    ///
    virtual string GetClassName() const override
    { return "SurfaceElement"; }

    ///
    void Update() override;

    ///
    // virtual size_t GetNDof () const throw() { return ndlevel.Last(); }

    ///
    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;

    ///
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;

    ///
    // virtual size_t GetNDofLevel (int level) const;

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
    bool all_the_same;
    bool do_subspace_update = true;
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
    auto & Spaces() const { return spaces; }

    ///
    string GetClassName () const override
    {
      return "CompoundFESpace";
    }

    /// updates also components
    void Update() override;
    /// updates also components
    void FinalizeUpdate() override;

    ProxyNode MakeProxyFunction (bool testfunction,
                                 const function<shared_ptr<ProxyFunction>(shared_ptr<ProxyFunction>)> & addblock) const override;

    /// copies dofcoupling from components
    void UpdateCouplingDofArray() override;
    virtual void UpdateFreeDofs() override;
    
    void SetDefinedOn (VorB vb, const BitArray& defon) override;

    DofRange GetRange (int spacenr) const
    {
      if (spacenr+1 >= cummulative_nd.Size())
        throw Exception("spacenr >= cummulative_nd.Size() in CompoundFESpace!");
      
      return DofRange(IntRange(cummulative_nd[spacenr], cummulative_nd[spacenr+1]),
                      spaces[spacenr]->GetParallelDofs());
    }

    shared_ptr<BaseMatrix> EmbeddingOperator (int spacenr) const;
    shared_ptr<BaseMatrix> RestrictionOperator (int spacenr) const;
    
    /// get component space
    shared_ptr<FESpace> operator[] (int i) const { return spaces[i]; }

    /// returns a compound finite element
    FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;
    ///
    void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    void GetDofNrs (NodeId ni, Array<DofId> & dnums) const override;
    void GetElementDofsOfType (ElementId ei, Array<DofId> & dnums, COUPLING_TYPE ctype) const override;
    ///
    [[deprecated("Use GetDofNrs(NODE_TYPE(NT_VERTEX,nr) instead")]]    
    void GetVertexDofNrs (int vnr, Array<DofId> & dnums) const override;
    [[deprecated("Use GetDofNrs(NODE_TYPE(NT_EDGE,nr) instead")]]    
    void GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const override;
    void GetFaceDofNrs (int fanr, Array<DofId> & dnums) const override;
    void GetInnerDofNrs (int elnr, Array<DofId> & dnums) const override;

    void SolveM(CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                        LocalHeap & lh) const override;
    void ApplyM(CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                        LocalHeap & lh) const override;
    
    template <class T> NGS_DLL_HEADER
      void T_TransformMat (ElementId ei, 
                           SliceMatrix<T> mat, TRANSFORM_TYPE tt) const;
    
    template <class T> NGS_DLL_HEADER
      void T_TransformVec (ElementId ei, 
                         SliceVector<T> vec, TRANSFORM_TYPE tt) const;

    void VTransformMR (ElementId ei,
                       SliceMatrix<double> mat, TRANSFORM_TYPE tt) const override;
    void VTransformMC (ElementId ei,
                       SliceMatrix<Complex> mat, TRANSFORM_TYPE tt) const override;
    void VTransformVR (ElementId ei,
                       SliceVector<double> vec, TRANSFORM_TYPE tt) const override;
    void VTransformVC (ElementId ei, 
                       SliceVector<Complex> vec, TRANSFORM_TYPE tt) const override;

    /// number of component spaces
    inline int GetNSpaces () const { return spaces.Size(); }

    void SetDoSubspaceUpdate(bool _do_subspace_update)
    { do_subspace_update = _do_subspace_update; }
  };



  class NGS_DLL_HEADER CompoundFESpaceAllSame : public CompoundFESpace
  {
    bool interleaved;
  public:
    CompoundFESpaceAllSame (shared_ptr<FESpace> space, int dim, const Flags & flags,
                            bool checkflags = false);
    virtual string GetClassName () const override;

    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override {
      return *new (alloc) VectorFiniteElement (spaces[0]->GetFE(ei, alloc), spaces.Size());
    }

    virtual FlatArray<VorB> GetDualShapeNodes (VorB vb) const override
    {
      return spaces[0]->GetDualShapeNodes(vb);
    }

    virtual shared_ptr<BaseMatrix> GetTraceOperator (shared_ptr<FESpace> tracespace, bool avg) const override;
  };

  class NGS_DLL_HEADER MatrixFESpace : public CompoundFESpace
  {
    bool symmetric;
    bool deviatoric;
    bool skewsymmetric;
    int vdim;
  public:
    MatrixFESpace (shared_ptr<FESpace> space, int avdim, const Flags & flags,
                   bool checkflags = false);
    virtual string GetClassName () const override;
    FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;    
  };


  
  template <typename BASESPACE>
  class VectorFESpace : public CompoundFESpace
  {
  public:
    VectorFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, 
                   bool checkflags = false)
      : CompoundFESpace (ama, flags)
    {
      Array<string> dirichlet_comp;
      string dirnames[] = { "dirichletx", "dirichlety", "dirichletz" };
      for (int i = 0; i <  ma->GetDimension(); i++)
        {
          Flags tmpflags = flags;
          if (flags.StringFlagDefined(dirnames[i]))
            tmpflags.SetFlag ("dirichlet", flags.GetStringFlag(dirnames[i]));
          if (flags.StringFlagDefined(dirnames[i]+"_bbnd"))
            tmpflags.SetFlag ("dirichlet_bbnd", flags.GetStringFlag(dirnames[i]+"_bbnd"));
          AddSpace (make_shared<BASESPACE> (ama, tmpflags));
        }

      for (auto vb : { VOL, BND, BBND, BBBND })
        {
          if (auto eval = spaces[0] -> GetEvaluator(vb))
            evaluator[vb] = make_shared<VectorDifferentialOperator> (eval, ma->GetDimension());
          if (auto fluxeval = spaces[0] -> GetFluxEvaluator(vb))
            flux_evaluator[vb] = make_shared<VectorDifferentialOperator> (fluxeval, ma->GetDimension());
        }

      auto additional = spaces[0]->GetAdditionalEvaluators();
      for (int i = 0; i < additional.Size(); i++)
        additional_evaluators.Set (additional.GetName(i),
                                   make_shared<VectorDifferentialOperator>(additional[i], ma->GetDimension()));

      type = "Vector"+(*this)[0]->type;
    }

    virtual string GetClassName () const override
    {
      return "Vector"+ (*this)[0]->GetClassName();
    }

    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override {
        return *new (alloc) VectorFiniteElement (spaces[0]->GetFE(ei, alloc), spaces.Size());
    }

    virtual FlatArray<VorB> GetDualShapeNodes (VorB vb) const override
    {
      return spaces[0]->GetDualShapeNodes(vb);
    }

  };


  
  class NGS_DLL_HEADER ApplyMass : public BaseMatrix
  {
  protected:
    shared_ptr<FESpace> fes;
    shared_ptr<CoefficientFunction> rho;
    bool inverse;
    shared_ptr<Region> definedon;
    LocalHeap & lh;
  public:
    ///
    ApplyMass (shared_ptr<FESpace> afes,
               shared_ptr<CoefficientFunction> arho,
               bool ainverse,
               shared_ptr<Region> adefinedon,
               LocalHeap & alh);
    virtual ~ApplyMass();
    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = nullptr) const override;

    virtual bool IsComplex() const override
    {
      return fes->IsComplex();
    }

    shared_ptr<BaseMatrix> CreateDeviceMatrix() const override
    {
      return fes->CreateMassOperator(rho, definedon, inverse, lh)->CreateDeviceMatrix();
    }
    
    
    virtual void Mult (const BaseVector & v, BaseVector & prod) const override;
    virtual void MultAdd (double val, const BaseVector & v, BaseVector & prod) const override;
    virtual void MultAdd (Complex val, const BaseVector & v, BaseVector & prod) const override;
    virtual void MultTransAdd (double val, const BaseVector & v, BaseVector & prod) const override;
    
    virtual AutoVector CreateVector () const override;
    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;
    
    virtual int VHeight() const override
    {
      return fes->GetNDof(); 
    }
    ///
    virtual int VWidth() const override
    {
      return fes->GetNDof();       
    }
  };



  class NGS_DLL_HEADER ApplyTrace : public BaseMatrix
  {
  protected:
    shared_ptr<FESpace> fes;
    shared_ptr<FESpace> festrace;
    bool average;
    LocalHeap & lh;
  public:
    ApplyTrace (shared_ptr<FESpace> afes,
                shared_ptr<FESpace> afestrace,
                bool aaverage,
                LocalHeap & alh);
    virtual ~ApplyTrace();

    virtual int VHeight() const override { return festrace->GetNDof(); }
    virtual int VWidth() const override { return fes->GetNDof(); }
    virtual bool IsComplex() const override { return fes->IsComplex(); }
    
    virtual void Mult (const BaseVector & v, BaseVector & prod) const override;
    virtual void MultAdd (double val, const BaseVector & v, BaseVector & prod) const override;
    virtual void MultAdd (Complex val, const BaseVector & v, BaseVector & prod) const override;
    virtual void MultTransAdd (double val, const BaseVector & v, BaseVector & prod) const override;
    
    virtual AutoVector CreateVector () const override;
    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;
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
      /// function pointer to docu function
      DocInfo (*getdocu)();
      /// creates a descriptor
      FESpaceInfo (const string & aname,
		   shared_ptr<FESpace> (*acreator)(shared_ptr<MeshAccess> ma, const Flags & flags),
                   DocInfo (*agetdocu)())
	: name(aname), creator(acreator), getdocu(agetdocu) {;}
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
		     shared_ptr<FESpace> (*acreator)(shared_ptr<MeshAccess> ma, const Flags & flags),
                     DocInfo (*getdocu)() = FESpace::GetDocu);
  
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
    RegisterClassForArchive<FES, FESpace> regclass;
    /// constructor registers fespace
    RegisterFESpace (string label)
    {
      GetFESpaceClasses().AddFESpace (label, Create, FES::GetDocu);
      // cout << "register fespace '" << label << "'" << endl;
    }
    
    /// creates an fespace of type FES
    static shared_ptr<FESpace> Create (shared_ptr<MeshAccess> ma, const Flags & flags)
    {
      return make_shared<FES> (ma, flags);
    }
  };

  
}





namespace ngcore
{
  template<typename T> struct MPI_typetrait;
  
  template<>
  struct MPI_typetrait<ngcomp::COUPLING_TYPE>
  {
    static auto MPIType () 
    { 
      static_assert ( (sizeof(ngcomp::COUPLING_TYPE) == sizeof(char)) ||
                      (sizeof(ngcomp::COUPLING_TYPE) == sizeof(int)) );
      if constexpr (sizeof(ngcomp::COUPLING_TYPE) == sizeof(char)) return MPI_typetrait<char>::MPIType();
      if constexpr (sizeof(ngcomp::COUPLING_TYPE) == sizeof(int)) return MPI_typetrait<int>::MPIType();
    }
  };
}


#endif
