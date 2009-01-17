#ifndef FILE_HDIVHYBIRDHOFESPACE

#define FILE_HDIVHYBRIDHOFESPACE

/*********************************************************************/
/* File:   hdivhybridhofespace.hpp                                   */
/* Author: Sabine Zaglmayr, JS                                       */
/* Date:   15. Feb. 2003, May 2006                                   */
/*********************************************************************/

/**
   HDiv High Order Finite Element Space
*/


class HDivHybridHighOrderFESpace : public FESpace
{
private:
  
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
  ARRAY<int> first_inner_dof;

  /// relative order to mesh-order
  int rel_order; 

  ARRAY<int> order_edge;
  ARRAY<int> order_inner;
  ARRAY<int> order_face;

  int ndof;
  ARRAY<int> ndlevel;

public:

  HDivHybridHighOrderFESpace (const MeshAccess & ama, const Flags & flags, bool parseflags=false);
  ///
  virtual ~HDivHybridHighOrderFESpace ();

  static FESpace * Create (const MeshAccess & ma, const Flags & flags);

  virtual string GetClassName () const
  {
    return "HDivHybridHighOrderFESpace";
  }

  ///
  virtual void Update();
  ///
  virtual int GetNDof () const;
  ///
  virtual int GetNDofLevel (int level) const;
  ///
  virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
  ///
  virtual const FiniteElement & GetSFE (int selnr, LocalHeap & lh) const; // 2D: array =0.;
  ///
  virtual void GetDofNrs (int elnr, ARRAY<int> & dnums) const;
  ///
  virtual void GetExternalDofNrs (int elnr, ARRAY<int> & dnums) const;
  ///
  virtual void GetSDofNrs (int selnr, ARRAY<int> & dnums) const;
  ///
  virtual Table<int> * CreateSmoothingBlocks (const Flags & precflags) const;


  /// 
  virtual void GetVertexDofNrs (int vnr, ARRAY<int> & dnums) const;
  /// 
  virtual void GetEdgeDofNrs (int ednr, ARRAY<int> & dnums) const;
  /// 
  virtual void GetFaceDofNrs (int fanr, ARRAY<int> & dnums) const;
  /// 
  virtual void GetInnerDofNrs (int elnr, ARRAY<int> & dnums) const;

};

#endif





