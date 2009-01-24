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

  Array<int> first_edge_dof;
  Array<int> first_face_dof;
  Array<int> first_inner_dof;

  /// relative order to mesh-order
  int rel_order; 

  Array<int> order_edge;
  Array<int> order_inner;
  Array<int> order_face;

  int ndof;
  Array<int> ndlevel;

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
  virtual void GetDofNrs (int elnr, Array<int> & dnums) const;
  ///
  virtual void GetExternalDofNrs (int elnr, Array<int> & dnums) const;
  ///
  virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;
  ///
  virtual Table<int> * CreateSmoothingBlocks (const Flags & precflags) const;


  /// 
  virtual void GetVertexDofNrs (int vnr, Array<int> & dnums) const;
  /// 
  virtual void GetEdgeDofNrs (int ednr, Array<int> & dnums) const;
  /// 
  virtual void GetFaceDofNrs (int fanr, Array<int> & dnums) const;
  /// 
  virtual void GetInnerDofNrs (int elnr, Array<int> & dnums) const;

};

#endif





