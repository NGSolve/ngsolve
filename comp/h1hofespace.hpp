#ifndef FILE_H1HOFESPACE
#define FILE_H1HOFESPACE

/*********************************************************************/
/* File:   h1hofespace.hpp                                           */
/* Author: Start                                                     */
/* Date:   10. Feb. 2003                                             */
/*********************************************************************/

/**
   High Order Finite Element Space
*/


class H1HighOrderFESpace : public FESpace
{
protected:

  // Level
  int level;

  // Number of Edges
  int ned;
  // Number of Faces
  int nfa;
  // Number of Elements
  int nel;
  // Number of Vertex
  int nv;

  ARRAY<int> first_edge_dof;
  ARRAY<int> first_face_dof;
  ARRAY<int> first_element_dof;

  ARRAY<bool> dirichlet_vertex;
  ARRAY<bool> dirichlet_edge;
  ARRAY<bool> dirichlet_face;

  /// relative order to mesh-order
  int rel_order; 
  bool var_order; 
  

  ARRAY<int> order_edge;
  ARRAY<INT<2> > order_face;
  ARRAY<INT<3> > order_inner;
  ARRAY<int> order_avertex; 
  ARRAY<bool> fine_edge; 
  ARRAY<bool> fine_face; 

  int ndof;
  int uniform_order_inner;
  int uniform_order_face;
  int uniform_order_edge;
  int uniform_order_quad;
  int uniform_order_trig;
  ARRAY<INT<3> > dom_order_min; 
  ARRAY<INT<3> > dom_order_max;
  int smoother; 
  
  ARRAY<int> ndlevel;

  int augmented;

  bool level_adapted_order; 


  int loworderindex; 

  bool minext, optext;
  bool fast_pfem;
  bool fast_pfem_sz;
  bool plate;  // thin plate: prism trig dofs are internal
  int testf1,testf2,testi1,testi2,testi3;

  ARRAY<INT<2> > defined_on_one_side_of_bounding_curve;

  bool print; 

public:

  H1HighOrderFESpace (const MeshAccess & ama, const Flags & flags, bool parseflags=false);
  ///
  virtual ~H1HighOrderFESpace ();

  static FESpace * Create (const MeshAccess & ma, const Flags & flags);

  virtual string GetClassName () const
  {
    return "H1HighOrderFESpace";
  }

  ///
  virtual void Update(LocalHeap & lh);
  ///
  virtual void PrintReport (ostream & ost);

  ///
  virtual int GetNDof () const;
  ///
  virtual int GetNDofLevel (int alevel) const;
  ///
  virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
  ///
  virtual const FiniteElement & GetSFE (int elnr, LocalHeap & lh) const;
  ///
  virtual void GetDofNrs (int elnr, ARRAY<int> & dnums) const;

  virtual void GetWireBasketDofNrs (int vnr, ARRAY<int> & dnums) const;
  virtual void GetVertexDofNrs (int vnr, ARRAY<int> & dnums) const;
  virtual void GetEdgeDofNrs (int ednr, ARRAY<int> & dnums) const;
  virtual void GetFaceDofNrs (int fanr, ARRAY<int> & dnums) const;
  virtual void GetInnerDofNrs (int elnr, ARRAY<int> & dnums) const;

  ///
  virtual void GetExternalDofNrs (int elnr, ARRAY<int> & dnums) const;
  ///
  virtual void GetSDofNrs (int selnr, ARRAY<int> & dnums) const;
  
  virtual Table<int> * CreateSmoothingBlocks (const Flags & precflags) const; 
  ///
  virtual ARRAY<int> * CreateDirectSolverClusters (const Flags & flags) const;

  ///
  int GetFirstFaceDof(int i) const {return(first_face_dof[i]);} ;  
  ///
  int GetFirstEdgeDof(int i) const {return(first_edge_dof[i]);} ; 
  ///
  int GetFirstElementDof(int i) const {return(first_element_dof[i]);} ; 

  bool IsDirichletVertex (int i) const { return dirichlet_vertex.Size() && dirichlet_vertex[i]; }
  bool IsDirichletEdge (int i) const { return dirichlet_edge.Size() && dirichlet_edge[i]; }
  bool IsDirichletFace (int i) const { return dirichlet_face.Size() && dirichlet_face[i]; }

  void UpdateDofTables ();
  void SetEdgeOrder (int enr, int eo) { order_edge[enr] = eo; }
  void SetFaceOrder (int fnr, int fo) { order_face[fnr] = INT<2> (fo, fo); }
  void SetFaceOrder (int fnr, int ox, int oy) { order_face[fnr] = INT<2> (ox, oy); }
  void SetElementOrder (int elnr, int elo) 
  { order_inner[elnr] = INT<3> (elo, elo, elo); }
  void SetElementOrder (int elnr, int ox, int oy, int oz) 
  { order_inner[elnr] = INT<3> (ox, oy, oz); }

  int GetAugmented() const { return augmented; }

  /// get relative (to mesh) order of finite elements
  virtual int GetRelOrder() const { return rel_order; }
  virtual bool VarOrder() const { return var_order; }

  void RestrictToOneSideOfBoundingCurve(int index1, int index2);
  void DeleteOneSideRestrictions(void);

#ifdef PARALLEL
#ifdef OLD_PARALLEL_UPDATE
  virtual void UpdateParallelDofs_hoproc();
#endif
  virtual void UpdateParallelDofs_loproc();
#endif
};

#endif

