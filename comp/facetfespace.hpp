/**
   High Order Finite Element Space for Facets
*/
/* facetfespace.*pp
 * An fespace for functions living on the facets (=edges in 2D, faces in 3D)
 *   - the functions on different facets are independent (no continuity accross
 *     vertices ( or edges in 3D))
 *   - the functions on the facets are accessed via the elements
 *     (=FacetVolumeElements), which then can access their facets;
 *     see facetfe.hpp for details
 *   - the ordering of dofs is facet'wise
 *   - no low-order space
 *   - no prolongation so far (will probably never be needed)
 *   - some additional utility functions 
 */
#ifndef FACET_FESPACE_HPP
#define FACET_FESPACE_HPP

namespace ngcomp
{

  class NGS_DLL_HEADER FacetFESpace : public FESpace 
  {
  protected:  
    // Level
    int level;
    // Number of Facets
    int nfa;
    // Number of coarse facets, number of fine facets;
    int ncfa;
    // Number of Elements
    int nel;
  
    Array<int> first_facet_dof;
    Array<int> first_inner_dof;  // for highest_order_dc
  
    // relative order to mesh-order
    int rel_order; 
  
    Array<INT<2> > order_facet;
    Array<bool> fine_facet;
  
    int ndof;
    Array<int> ndlevel;
    bool var_order; 
    bool print; 
    bool highest_order_dc;
    
  public:
    ///
    FacetFESpace (const MeshAccess & ama, const Flags & flags, bool parseflags=false);
    ///
    virtual ~FacetFESpace ();
    ///
    virtual string GetClassName () const
    {
      return "FacetFESpace";
    }
  
    ///
    virtual void Update(LocalHeap & lh);
  
    //  virtual void UpdateDofTables();
    virtual void UpdateCouplingDofArray();    
    ///
    virtual int GetNDof () const;
    ///
    virtual int GetNDofLevel (int level) const;
    ///
    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
    ///
    virtual const FiniteElement & GetSFE (int selnr, LocalHeap & lh) const; 
    ///
    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;
    ///

    IntRange GetFacetDofs (int nr) const
    { 
      return IntRange (first_facet_dof[nr], first_facet_dof[nr+1]); 
    }
  
    virtual void GetFacetDofNrs (int nr, Array<int> & dnums) const
    {
      dnums.SetSize(0);
      dnums += nr;
      dnums += GetFacetDofs(nr);
    }

    ///
    virtual int GetNFacetDofs (int felnr) const 
    { return (first_facet_dof[felnr+1]-first_facet_dof[felnr] + 1); }
    ///
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;
    ///
    virtual Table<int> * CreateSmoothingBlocks (const Flags & precflags) const;
    ///
    virtual Array<int> * CreateDirectSolverClusters (const Flags & precflags) const;
     
    // some utility functions for convenience
    ///
    virtual void GetVertexNumbers(int elnr, Array<int>& vnums) 
    { ma.GetElVertices(elnr, vnums); };
    ///
    virtual INT<2> GetFacetOrder(int fnr) 
    { return order_facet[fnr]; };
  

    virtual int GetFirstFacetDof(int fanr) const 
    {
      return first_facet_dof[fanr];
    }

    virtual void GetVertexDofNrs ( int nr, Array<int> & dnums ) const
    {
      dnums.SetSize(0);
    }

    virtual void GetEdgeDofNrs ( int nr, Array<int> & dnums ) const
    {
      dnums.SetSize(0);
      if ( ma.GetDimension() == 3 ) return;

      dnums += nr;
      dnums += GetFacetDofs(nr);
    }

    virtual void GetFaceDofNrs (int nr, Array<int> & dnums) const
    {
      dnums.SetSize(0);
      if (ma.GetDimension() == 2) return;

      dnums += nr;
      dnums += GetFacetDofs(nr);
    }
  
    virtual void GetInnerDofNrs (int elnr, Array<int> & dnums) const
    {
      dnums.SetSize(0);
    }


#ifdef PARALLEL_NOT_JS_20110614
    // virtual void UpdateParallelDofs_hoproc();
    virtual void UpdateParallelDofs_loproc();
#endif
  
  };









  class NGS_DLL_HEADER EdgeFESpace : public FESpace 
  {
  protected:  
    int ned;
  
    Array<int> first_edge_dof;
  
    Array<int> order_edge;
    Array<bool> fine_edge;
  
    int ndof;
    bool print; 
  
  public:
    ///
    EdgeFESpace (const MeshAccess & ama, const Flags & flags, bool parseflags=false);
    ///
    virtual ~EdgeFESpace ();
    ///
    // static FESpace * Create (const MeshAccess & ma, const Flags & flags);
    ///
    virtual string GetClassName () const
    {
      return "EdgeFESpace";
    }
  
    ///
    virtual void Update(LocalHeap & lh);
    ///
    virtual void UpdateCouplingDofArray();    
    ///
    virtual int GetNDof () const { return first_edge_dof[ned]; }
    ///
    virtual int GetNDofLevel (int level) const { return GetNDof(); }
    ///
    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
    ///
    virtual const FiniteElement & GetSFE (int selnr, LocalHeap & lh) const; 
    ///
    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;

    virtual void GetVertexDofNrs ( int elnum, Array<int> & dnums ) const  { dnums.SetSize(0); }

    virtual void GetEdgeDofNrs ( int elnum, Array<int> & dnums ) const;

    virtual void GetFaceDofNrs ( int elnum, Array<int> & dnums ) const  { dnums.SetSize(0); }
    virtual void GetInnerDofNrs (int elnr, Array<int> & dnums) const { dnums.SetSize(0); }

    ///
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;
    ///
    virtual Table<int> * CreateSmoothingBlocks (const Flags & precflags) const;
    ///
    virtual Array<int> * CreateDirectSolverClusters (const Flags & precflags) const;
  };     







}



#endif
