#ifndef VECTOR_FACET_FESPACE_HPP
#define VECTOR_FACET_FESPACE_HPP

#include <fem.hpp>
#include <comp.hpp>

class VectorFacetFESpace : public FESpace
{
protected:
  /// Level
  int level;
  /// Number of Facets
  int nfacets;
  /// 
  int ncfacets;
  ///
  int nel;

  ARRAY<int> first_facet_dof;
  int ndof_lo;

  int rel_order;

  ARRAY<INT<2> > order_facet;
  ARRAY<bool> fine_facet;

  int ndof;
  ARRAY<int> ndlevel;
  bool var_order;
  bool print;



public:
  ///
  /*
  VectorFacetFESpace () :
    FESpace (),
    level(0),
    nfacets(0),
    ncfacets(0),
    nel(0),
    ndof(0),
    ndof_lo(0),
    var_order(0),
    print(0)
  {
    first_facet_dof.SetSize(0);
    order_facet.SetSize(0);
    fine_facet.SetSize(0);
    ndlevel.SetSize(0);
  }
  */

  VectorFacetFESpace ( const MeshAccess & ama, const Flags & flags, 
		       bool parseflags = false );

  virtual ~VectorFacetFESpace ()
  {
    ;
  }

  static FESpace * Create ( const MeshAccess & ma, const Flags & flags );

  virtual string GetClassName () const 
  {
    return "VectorFacetFESpace";
  }

  virtual void Update(LocalHeap& lh);

  virtual int GetNDof() const { return ndof; }

  virtual int GetNDofLevel ( int i ) const { return ndlevel[i]; }

  virtual int GetNDofLowOrder () const
  { return ndof_lo; }

  virtual const FiniteElement & GetFE ( int elnr, LocalHeap & lh ) const;
  virtual const FiniteElement & GetSFE ( int selnr, LocalHeap & lh ) const;

  virtual void GetFacetDofNrs ( int felnr, ARRAY<int> & dnums ) const;

  virtual int GetNFacetDofs ( int felnr ) const;

  virtual void GetDofNrs ( int elnr, ARRAY<int> & dnums ) const;

  virtual void GetWireBasketDofNrs(int elnr, ARRAY<int> & dnums) const;
  ///
  virtual void GetExternalDofNrs (int elnr, ARRAY<int> & dnums) const;
  ///
  virtual void GetSDofNrs (int selnr, ARRAY<int> & dnums) const;
  ///
  virtual Table<int> * CreateSmoothingBlocks (const Flags & precflags) const;
  ///
  virtual ARRAY<int> * CreateDirectSolverClusters (const Flags & precflags) const;
  
  // some utility functions for convenience
  ///
  virtual void GetVertexNumbers(int elnr, ARRAY<int>& vnums) ;
  ///
  virtual INT<2> GetFacetOrder(int fnr) ;

    virtual int GetFirstFacetDof(int fanr) const;


  virtual void GetVertexDofNrs ( int elnum, ARRAY<int> & dnums ) const;

  virtual void GetEdgeDofNrs ( int elnum, ARRAY<int> & dnums ) const;

  virtual void GetFaceDofNrs (int felnr, ARRAY<int> & dnums) const;
  
#ifdef PARALLEL
   virtual void UpdateParallelDofs_hoproc();
   virtual void UpdateParallelDofs_loproc();
#endif



};

#endif

