#ifndef FILE_SUPERLUINVERSE
#define FILE_SUPERLUINVERSE

/* *************************************************************************/
/* File:   superluinverse.hpp                                              */
/* Author: Florian Bachinger                                               */
/* Date:   Feb. 04                                                         */
/* *************************************************************************/

namespace ngla
{

  /*
    interface to the sparse direct solver SuperLU 
  */


#ifdef USE_SUPERLU

  namespace double_superlu
  {
    //#include "superlu/dsp_defs.h"
    //#include "superlu/util.h"
    extern "C"
    {
    
#include <slu_ddefs.h>
#include <slu_util.h>
#include <supermatrix.h>
    }
  
  }
  namespace complex_superlu
  {
    extern "C"
    {
      using namespace double_superlu;
#include <slu_zdefs.h>
#include <slu_util.h>
#include <supermatrix.h>
#include <slu_scomplex.h>
    }

  }

  using namespace double_superlu;
  using namespace complex_superlu;

#endif

  ////////////////////////////////////////////////////////////////////////////////


template<class TM, 
   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
   class TV_COL = typename mat_traits<TM>::TV_COL>
  class SuperLUInverse : public BaseMatrix
  {
#ifdef USE_SUPERLU
    SuperMatrix A, B, L, U;

    SuperLUStat_t     stat;
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
    SuperLUInverse (const SparseMatrix<TM,TV_ROW,TV_COL> & a, 
		    const BitArray * ainner = NULL,
		    const Array<int> * acluster = NULL,
		    int symmetric = 0);
    //   ///
    //   SuperLUInverse (const SparseMatrix<TM> & a, 
    // 		  const BitArray * ainner = NULL,
    // 		  const Array<int> * acluster = NULL,
    // 		  int symmetric = 0);
  
    ///
    SuperLUInverse (const Array<int> & aorder, 
		    const Array<CliqueEl*> & cliques,
		    const Array<MDOVertex> & vertices,
		    int symmetric = 0);		  
    ///
    ~SuperLUInverse ();
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

    virtual Array<MemoryUsage> GetMemoryUsage () const
    {
      return { MemoryUsage ("SuperLUInverse", nze*sizeof(TM), 1) };
    }

    

/*    ///
    void Set (int i, int j, const TM & val);
    ///
    const TM & Get (int i, int j) const;
*/
    virtual AutoVector CreateVector () const
    {
      return new VVector<TV> (height);
    }
  };

}

#endif
