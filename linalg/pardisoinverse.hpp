#ifndef FILE_PARDISOINVERSE
#define FILE_PARDISOINVERSE

/* *************************************************************************/
/* File:   pardisoinverse.hpp                                              */
/* Author: Florian Bachinger                                               */
/* Date:   Feb. 04                                                         */
/* *************************************************************************/

/*
  interface to the PARDISO - package
*/





////////////////////////////////////////////////////////////////////////////////


template<class TM, class TV_ROW, class TV_COL>
class PardisoInverse : public BaseMatrix
{
  int height, nze, entrysize;

  int pt[64];
  int hparams[64];
  int * rowstart, * indices;
  typename mat_traits<TM>::TSCAL * matrix;

  int symmetric, matrixtype, spd;

  const BitArray * inner;
  const ARRAY<int> * cluster;

  //
  void SetMatrixType(TM entry);
  
public:
  typedef typename mat_traits<TM>::TV_COL TV;
  typedef typename mat_traits<TM>::TV_ROW TVX;
  typedef typename mat_traits<TM>::TSCAL TSCAL;

  ///
  PardisoInverse (const SparseMatrix<TM,TV_ROW,TV_COL> & a, 
		  const BitArray * ainner = NULL,
		  const ARRAY<int> * acluster = NULL,
		  int symmetric = 0);
  
  ///
  PardisoInverse (const ARRAY<int> & aorder, 
		  const ARRAY<CliqueEl*> & cliques,
		  const ARRAY<MDOVertex> & vertices,
		  int symmetric = 0);		  
  ///
  virtual ~PardisoInverse ();
  ///
  int VHeight() const { return height; }
  ///
  int VWidth() const { return height; }
  ///
  void Allocate (const ARRAY<int> & aorder, 
		 const ARRAY<CliqueEl*> & cliques,
		 const ARRAY<MDOVertex> & vertices);
  ///
  void Factor (const int * blocknr);
  ///
  void FactorNew (const SparseMatrix<TM,TV_ROW,TV_COL> & a);
  ///
  virtual void Mult (const BaseVector & x, BaseVector & y) const;

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
  //  void SetOrig (int i, int j, const TM & val)
  //{ ; }

  virtual BaseVector * CreateVector () const
  {
    return new VVector<TV> (height);
  }
};



#endif
