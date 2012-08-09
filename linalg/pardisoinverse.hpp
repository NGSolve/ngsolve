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


namespace ngla
{
  using ngbla::integer;

  template<class TM, 
	   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class PardisoInverse : public SparseFactorization
  {
    integer height;             // matrix size in scalars
    integer compressed_height;  // matrix size after compression in scalars
    integer nze, entrysize;
    bool print;

    //CL: is this also working on 32-bit architectures?
    integer pt[128];
    integer hparams[64];

    Array<integer> rowstart, indices; 
    typename mat_traits<TM>::TSCAL * matrix;

    integer matrixtype;
    int symmetric, spd;

    void SetMatrixType(); // TM entry

    bool compressed;
    Array<int> compress;
    BitArray used;
  
  public:
    typedef TV_COL TV;
    typedef TV_ROW TVX;
    typedef typename mat_traits<TM>::TSCAL TSCAL;

    ///
    PardisoInverse (const SparseMatrix<TM,TV_ROW,TV_COL> & a, 
		    const BitArray * ainner = NULL,
		    const Array<int> * acluster = NULL,
		    int symmetric = 0);
    ///
    virtual ~PardisoInverse ();
    ///
    int VHeight() const { return height/entrysize; }
    ///
    int VWidth() const { return height/entrysize; }
    ///
    void Factor (const int * blocknr);
    ///
    void FactorNew (const SparseMatrix<TM,TV_ROW,TV_COL> & a);
    ///
    virtual void Mult (const BaseVector & x, BaseVector & y) const;


    ///
    virtual ostream & Print (ostream & ost) const;

    virtual void MemoryUsage (Array<MemoryUsageStruct*> & mu) const
    {
      mu.Append (new MemoryUsageStruct ("SparseChol", nze*sizeof(TM), 1));
    }
    virtual BaseVector * CreateVector () const
    {
      return new VVector<TV> (height/entrysize);
    }
  };

}

#endif
