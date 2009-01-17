#ifndef FILE_SPARSECHOLESKY
#define FILE_SPARSECHOLESKY

/* *************************************************************************/
/* File:   sparsecholesky.hpp                                              */
/* Author: Joachim Schoeberl                                               */
/* Date:   18. Jun. 97                                                     */
/* *************************************************************************/

/*
  sparse cholesky factorization
*/


class BaseSparseCholesky : public BaseMatrix
{ 
public:
  virtual void Smooth (BaseVector & u, const BaseVector & f, BaseVector & y) const = 0;
};


/**
   A sparse cholesky factorization.
   The unknowns are reordered by the minimum degree
   ordering algorithm
*/
template<class TM, class TV_ROW, class TV_COL>
class SparseCholesky : public BaseSparseCholesky
{
  int height, nze;

  MoveableMem<int> order, firstinrow, firstinrow_ri, rowindex2, blocknrs;
  ///
  MoveableMem<TM> lfact;
  MoveableMem<TM> diag;
  ///
  BitArray * inner;
  ARRAY<int> * cluster;
  ///
  MinimumDegreeOrdering * mdo;
  int maxrow;
  const SparseMatrix<TM,TV_ROW,TV_COL> & mat;

public:
  typedef TV_COL TV;
  typedef TV_ROW TVX;
  typedef typename mat_traits<TV_ROW>::TSCAL TSCAL_VEC;
  typedef typename mat_traits<TM>::TSCAL TSCAL_MAT;

  ///
  SparseCholesky (const SparseMatrix<TM,TV_ROW,TV_COL> & a, 
		  BitArray * ainner = NULL,
		  ARRAY<int> * acluster = NULL,
		  bool allow_refactor = 0);
  ///
  ~SparseCholesky ();
  ///
  int VHeight() const { return height; }
  ///
  int VWidth() const { return height; }
  ///
  void Allocate (const ARRAY<int> & aorder, 
		 const ARRAY<MDOVertex> & vertices,
		 const int * blocknr);
  ///
  void Factor (); // const int * blocknr);
  ///
  void FactorNew (const SparseMatrix<TM,TV_ROW,TV_COL> & a);
  ///
  virtual void Mult (const BaseVector & x, BaseVector & y) const;
  ///
  virtual void MultAdd (TSCAL_VEC s, const BaseVector & x, BaseVector & y) const;
  /**
     A = L+D+L^T
     y = f - (L+D)^T u
     w = C^{-1} (y - L u)
     u += w
     y -= (L+D)^T w
   **/
  virtual void Smooth (BaseVector & u, const BaseVector & f, BaseVector & y) const;
  ///
  virtual ostream & Print (ostream & ost) const;

  virtual void MemoryUsage (ARRAY<MemoryUsageStruct*> & mu) const
  {
    mu.Append (new MemoryUsageStruct ("SparseChol", nze*sizeof(TM), 1));
  }


  ///
  void Set (int i, int j, const TM & val);
  ///
  const TM & Get (int i, int j) const;
  ///
  void SetOrig (int i, int j, const TM & val)
  { Set (order[i], order[j], val); }

  virtual BaseVector * CreateVector () const
  {
    return new VVector<TV> (height);
  }
};



#endif
