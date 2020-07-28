#ifndef FILE_MUMPSINVERSE
#define FILE_MUMPSINVERSE

/* *************************************************************************/
/* File:   mumpsinverse.hpp                                               */
/* Author: Joachin Schoeberl                                               */
/* Date:   Apr. 09                                                         */
/* *************************************************************************/

/*
namespace ngparallel
{
  class ParallelDofs;
}
*/

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



  template<class TM, 
	   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class MumpsInverse : public BaseMatrix
  {
    typedef typename mat_traits<TM>::TV_COL TV;
    typedef typename mat_traits<TM>::TV_ROW TVX;
    typedef typename mat_traits<TM>::TSCAL TSCAL;


    typedef typename mumps_trait<TSCAL>::MUMPS_STRUC_C MUMPS_STRUC_C;
    MUMPS_STRUC_C mumps_id;

    int height, nze, entrysize;

    bool symmetric, iscomplex;

    shared_ptr<BitArray> inner;
    shared_ptr<const Array<int>> cluster;

    NgMPI_Comm comm;

  public:
    ///
    MumpsInverse (const SparseMatrix<TM,TV_ROW,TV_COL> & a, 
                  shared_ptr<BitArray> ainner = nullptr,
                  shared_ptr<const Array<int>> acluster = nullptr,
                  bool symmetric = false);

    ///
    ~MumpsInverse ();

    ///
    int VHeight() const { return height/entrysize; }
    
    ///
    int VWidth() const { return height/entrysize; }
    ///
    virtual bool IsComplex() const { return iscomplex; }
    ///
    virtual void Mult (const BaseVector & x, BaseVector & y) const;

    ///
    virtual AutoVector CreateVector () const
    {
      return make_unique<VVector<TV>> (height);
    }

    virtual AutoVector CreateRowVector () const
    {
      return make_unique<VVector<TV>> (height);
    }

    virtual AutoVector CreateColVector () const
    {
      return make_unique<VVector<TV>> (height);
    }

  };





  template<class TM, 
	   class TVEC = typename mat_traits<TM>::TV_ROW> 
  class ParallelMumpsInverse : public BaseMatrix
  {
    typedef typename mat_traits<TM>::TV_COL TV;
    typedef typename mat_traits<TM>::TSCAL TSCAL;


    typedef typename mumps_trait<TSCAL>::MUMPS_STRUC_C MUMPS_STRUC_C;
    MUMPS_STRUC_C mumps_id;
    int height, nze, entrysize;

    bool symmetric, iscomplex;

    // const BitArray * inner;
    // const Array<int> * cluster;
    Array<int> select;
    Array<int> loc2glob;
    int num_globdofs;
  public:
    ///
    ParallelMumpsInverse (const BaseSparseMatrix & a, 
			  shared_ptr<BitArray> ainner,
			  shared_ptr<const Array<int>> acluster,
			  shared_ptr<ParallelDofs> pardofs,
			  bool symmetric = 0);

    ///
    ~ParallelMumpsInverse ();

    ///
    int VHeight() const { return height; }
    
    ///
    int VWidth() const { return height; }
    ///
    virtual bool IsComplex() const { return iscomplex; }
    ///
    virtual void Mult (const BaseVector & x, BaseVector & y) const;

    ///
    virtual AutoVector CreateVector () const;
    virtual AutoVector CreateRowVector () const;
    virtual AutoVector CreateColVector () const;
    /*
    {
      return new ParallelVVector<TV> (height, paralleldofs);
      // return new VVector<TV> (height);
    }
    */

  private:
    void MumpsFunction (MUMPS_STRUC_C & mumps_id)
    {
      mumps_trait<TSCAL>::MumpsFunction (&mumps_id);
    }
  };





#endif




}

#endif
