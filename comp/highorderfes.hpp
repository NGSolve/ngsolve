#ifdef OLD_HIGHORDERFES

#ifndef FILE_HIGHORDERFES
#define FILE_HIGHORDERFES

/*********************************************************************/
/* File:   highorderfes.hh                                           */
/* Author: Joachim Schoeberl                                         */
/* Date:   28. Oct. 2000                                             */
/*********************************************************************/

/** 
   High Order Finite Element Space
*/
class NodalFESpaceP : public FESpace
{
  ///
  int p;
  ///
  ARRAY<int> eldofs;
  ///
  int ndof;

  int nv, nedge, nface, ne;

  int n_edge_dofs;
  int n_face_dofs;
  int n_el_dofs;
  int augmented;
public:

  ///
  NodalFESpaceP (const MeshAccess & ama, 	       
		 int ap, int adim, bool acomplex);

  ///
  ~NodalFESpaceP ();

  static FESpace * Create (const MeshAccess & ma, const Flags & flags);

  virtual string GetClassName () const
  {
    return "NodalFESpaceP";
  }

  ///
  virtual void Update();

  ///
  virtual int GetNDof () const
  { return ndof; }

  ///
  virtual void GetDofNrs (int elnr, ARRAY<int> & dnums) const;

  ///
  virtual void GetSDofNrs (int selnr, ARRAY<int> & dnums) const;

  ///
  int GetEdgeDof (int enr, int i1, int i2, int lam1, int lam2) const;
  ///
  int GetFaceDof (int fnr, int i1, int i2, int i3, int lam1, int lam2, int lam3) const;
  ///
  int GetQuadFaceDof (int fnr, int i1, int i2, int i3, int i4, 
		      int lam1, int lam2) const;
  ///
  int GetElementDof (int elnr, int lam1, int lam2, int lam3, int lam4) const;




  template <class MAT>
  void TransformMat (int elnr, bool boundary,
		     MAT & mat, TRANSFORM_TYPE tt) const;

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


  ///
  virtual Table<int> * CreateSmoothingBlocks (int type = 0) const;
};






#endif


#endif
