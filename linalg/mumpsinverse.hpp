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


not available yet



  namespace double_mumps
  {
    //#include "mumps/dsp_defs.h"
    //#include "mumps/util.h"
    extern "C"
    {
    
#include <slu_ddefs.h>
#include <slu_util.h>
#include <supermatrix.h>
    }
  
  }
  namespace complex_mumps
  {
    extern "C"
    {
      using namespace double_mumps;
#include <slu_zdefs.h>
#include <slu_util.h>
#include <supermatrix.h>
#include <slu_scomplex.h>
    }

  }

  using namespace double_mumps;
  using namespace complex_mumps;

#endif

  ////////////////////////////////////////////////////////////////////////////////


  template<class TM, class TV_ROW, class TV_COL>
  class MumpsInverse : public BaseMatrix
  {
#ifdef USE_MUMPS
    SuperMatrix A, B, L, U;

    MumpsStat_t     stat;
#endif
    int height, nze, entrysize;

    int * perm_c, * perm_r;
    int * colstart, * indices;
    typename mat_traits<TM>::TSCAL * matrix;

    int symmetric, iscomplex;

    const BitArray * inner;
    const Array<int> * cluster;

  public:
    typedef typename mat_traits<TM>::TV_COL TV;
    typedef typename mat_traits<TM>::TV_ROW TVX;
    typedef typename mat_traits<TM>::TSCAL TSCAL;

    ///
    MumpsInverse (const SparseMatrix<TM,TV_ROW,TV_COL> & a, 
		    const BitArray * ainner = NULL,
		    const Array<int> * acluster = NULL,
		    int symmetric = 0);
    //   ///
    //   MumpsInverse (const SparseMatrix<TM> & a, 
    // 		  const BitArray * ainner = NULL,
    // 		  const Array<int> * acluster = NULL,
    // 		  int symmetric = 0);
  
    ///
    MumpsInverse (const Array<int> & aorder, 
		    const Array<CliqueEl*> & cliques,
		    const Array<MDOVertex> & vertices,
		    int symmetric = 0);		  
    ///
    ~MumpsInverse ();
    ///
    int VHeight() const { return height; }
    ///
    int VWidth() const { return height; }
    ///
    void Allocate (const Array<int> & aorder, 
		   const Array<CliqueEl*> & cliques,
		   const Array<MDOVertex> & vertices);
    ///
    void Factor (const int * blocknr);
    ///
    void FactorNew (const SparseMatrix<TM> & a);
    ///
    virtual void Mult (const BaseVector & x, BaseVector & y) const;

    ///
    virtual ostream & Print (ostream & ost) const;

    virtual void MemoryUsage (Array<MemoryUsageStruct*> & mu) const
    {
      mu.Append (new MemoryUsageStruct ("MumpsInverse", nze*sizeof(TM), 1));
    }

    ///
    void Set (int i, int j, const TM & val);
    ///
    const TM & Get (int i, int j) const;

    virtual BaseVector * CreateVector () const
    {
      return new VVector<TV> (height);
    }
  };

}

#endif
