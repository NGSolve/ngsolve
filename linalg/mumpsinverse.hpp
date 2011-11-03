#ifndef FILE_MUMPSINVERSE
#define FILE_MUMPSINVERSE

/* *************************************************************************/
/* File:   mumpsinverse.hpp                                               */
/* Author: Joachin Schoeberl                                               */
/* Date:   Apr. 09                                                         */
/* *************************************************************************/

namespace ngparallel
{
  class ParallelDofs;
}

namespace ngla
{

  /*
    interface to the sparse direct solver Mumps 
  */

#ifdef USE_MUMPS
#include "dmumps_c.h"
#include "zmumps_c.h"


  template <class TSCAL>
  class mumps_trait
  {
  public:
    typedef DMUMPS_STRUC_C MUMPS_STRUC_C;
    typedef double MUMPS_TSCAL;
    static void MumpsFunction (DMUMPS_STRUC_C * id)
    { dmumps_c (id); }
  };

  template <>
  class mumps_trait<Complex>
  {
  public:
    typedef ZMUMPS_STRUC_C MUMPS_STRUC_C;
    typedef mumps_double_complex MUMPS_TSCAL;
    static void MumpsFunction (ZMUMPS_STRUC_C * id)
    { zmumps_c (id); }
  };
#endif



  template<class TM, 
	   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class MumpsInverse : public BaseMatrix
  {
    typedef typename mat_traits<TM>::TV_COL TV;
    typedef typename mat_traits<TM>::TV_ROW TVX;
    typedef typename mat_traits<TM>::TSCAL TSCAL;

#ifdef USE_MUMPS
    typedef typename mumps_trait<TSCAL>::MUMPS_STRUC_C MUMPS_STRUC_C;
    MUMPS_STRUC_C mumps_id;
#endif
    int height, nze, entrysize;

    bool symmetric, iscomplex;

    const BitArray * inner;
    const Array<int> * cluster;

  public:
    ///
    MumpsInverse (const SparseMatrix<TM,TV_ROW,TV_COL> & a, 
                  const BitArray * ainner = NULL,
                  const Array<int> * acluster = NULL,
                  bool symmetric = 0);

    ///
    ~MumpsInverse ();

    ///
    int VHeight() const { return height; }
    
    ///
    int VWidth() const { return height; }

    ///
    virtual void Mult (const BaseVector & x, BaseVector & y) const;

    ///
    virtual BaseVector * CreateVector () const
    {
      return new VVector<TV> (height);
    }
  };







  template<class TM, 
	   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class ParallelMumpsInverse : public BaseMatrix
  {
    typedef typename mat_traits<TM>::TV_COL TV;
    typedef typename mat_traits<TM>::TV_ROW TVX;
    typedef typename mat_traits<TM>::TSCAL TSCAL;

#ifdef USE_MUMPS
    typedef typename mumps_trait<TSCAL>::MUMPS_STRUC_C MUMPS_STRUC_C;
    MUMPS_STRUC_C mumps_id;
#endif
    int height, nze, entrysize;

    bool symmetric, iscomplex;

    const BitArray * inner;
    const Array<int> * cluster;
    Array<int> select;
    Array<int> loc2glob;
    int num_globdofs;
  public:
    ///
    ParallelMumpsInverse (const SparseMatrixTM<TM> & a, 
			  const BitArray * ainner,
			  const Array<int> * acluster,
			  const ngparallel::ParallelDofs * pardofs,
			  bool symmetric = 0);

    ///
    ~ParallelMumpsInverse ();

    ///
    int VHeight() const { return height; }
    
    ///
    int VWidth() const { return height; }

    ///
    virtual void Mult (const BaseVector & x, BaseVector & y) const;

    ///
    virtual BaseVector * CreateVector () const
    {
      return new VVector<TV> (height);
    }
  };










}

#endif
