#ifndef FILE_MUMPSINVERSE
#define FILE_MUMPSINVERSE

/* *************************************************************************/
/* File:   mumpsinverse.hpp                                               */
/* Author: Joachin Schoeberl                                               */
/* Date:   Apr. 09                                                         */
/* *************************************************************************/

namespace ngla
{

  /*
    interface to the sparse direct solver Mumps 
  */

#ifdef USE_MUMPS
#include "dmumps_c.h"
#include "cmumps_c.h"


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
    typedef CMUMPS_STRUC_C MUMPS_STRUC_C;
    typedef mumps_complex MUMPS_TSCAL;
    static void MumpsFunction (CMUMPS_STRUC_C * id)
    { cmumps_c (id); }
  };
#endif



  template<class TM, class TV_ROW, class TV_COL>
  class MumpsInverse : public BaseMatrix
  {
    typedef typename mat_traits<TM>::TV_COL TV;
    typedef typename mat_traits<TM>::TV_ROW TVX;
    typedef typename mat_traits<TM>::TSCAL TSCAL;


#ifdef USE_MUMPS
    typedef typename mumps_trait<TSCAL>::MUMPS_STRUC_C MUMPS_STRUC_C;
    MUMPS_STRUC_C id;
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
    // virtual ostream & Print (ostream & ost) const;
    ///
    virtual BaseVector * CreateVector () const
    {
      return new VVector<TV> (height);
    }
  };

}

#endif
