#ifndef FILE_SUPERLU_DIST
#define FILE_SUPERLU_DIST

/* *************************************************************************/
/* File:   superlu_dist.hpp                                                */
/* Author: Astrid                                                          */
/* Date:   2007                                                            */
/* *************************************************************************/

/*
  interface to the sparse direct solver SuperLU_DIST
*/


#ifdef PARALLEL

#include <parallelngs.hpp>
#include <la.hpp>

#ifdef USE_SUPERLU_DIST
namespace double_superlu_dist
{
  extern "C" {
#include <superlu_ddefs.h>
#include <superlu_defs.h>
#include <supermatrix.h>
#include <util_dist.h>
  }

}
namespace complex_superlu_dist
{
  extern "C" {
    using namespace double_superlu_dist;
    #include <superlu_zdefs.h>
    //#include <superlu_defs.h>
    //#include <supermatrix.h>
    //#include <util_dist.h>
  }
}



#endif
//using namespace double_superlu_dist;
//using namespace complex_superlu_dist;
using namespace ngparallel;


////////////////////////////////////////////////////////////////////////////////


template<class TM, class TV_ROW, class TV_COL>
class SuperLU_DIST_Inverse : public ParallelBaseMatrix
{
#ifdef USE_SUPERLU_DIST
  double_superlu_dist::SuperMatrix A;  // matrix which contains complete rows belonging to masterdofs
  // whatever 
  double_superlu_dist::ScalePermstruct_t scalepermstruct;
  // stores lu factorization
  double_superlu_dist::LUstruct_t double_lustruct;
  complex_superlu_dist::LUstruct_t complex_lustruct;
  // 
  double_superlu_dist::SOLVEstruct_t double_solvestruct;
  complex_superlu_dist::SOLVEstruct_t complex_solvestruct;
  double_superlu_dist::gridinfo_t gridinfo;

  double_superlu_dist::SuperLUStat_t     stat;
  double_superlu_dist::superlu_options_t options;

#endif

  int nrows_gridp, ncols_gridp;

  int height_local_tm, height_local_scal, nnz_local_tm, nnz_local_scal, entrysize;
  int height_global_tm, height_global_scal;

  // matrix graph, consists of
  int * rowptr, *colind;
  typename mat_traits<TM>::TSCAL * matrix;

  //  int * colstart, * indices;
  //  typename mat_traits<TM>::TSCAL * matrix;

  int symmetric, iscomplex;
  int first_row_tm, first_row_scal;
  const BitArray * inner;
  const ARRAY<int> * cluster;

  ARRAY<int> * index_tm;

public:
  typedef typename mat_traits<TM>::TV_COL TV;
  typedef typename mat_traits<TM>::TV_ROW TVX;
  typedef typename mat_traits<TM>::TSCAL TSCAL;

  ///
  SuperLU_DIST_Inverse (const ParallelSparseMatrix<TM,TV_ROW,TV_COL> & a, 
		  const BitArray * ainner = NULL,
		  const ARRAY<int> * acluster = NULL,
		  int symmetric = 0);
//   ///
//   SuperLUInverse (const SparseMatrix<TM> & a, 
// 		  const BitArray * ainner = NULL,
// 		  const ARRAY<int> * acluster = NULL,
// 		  int symmetric = 0);
  
  ///
//   SuperLU_DIST_Inverse (const ARRAY<int> & aorder, 
// 		  const ARRAY<CliqueEl*> & cliques,
// 		  const ARRAY<MDOVertex> & vertices,
// 		  int symmetric = 0);		  
  ///
  ~SuperLU_DIST_Inverse ();
  ///
  int VHeight_Local() const { return height_local_tm; }
  int VHeight_Global() const { return height_global_tm; }
  ///
  int VWidth_Local() const { return height_global_tm; }
  int VWidth_Global() const { return height_global_tm; }
  ///
  void Allocate (const ARRAY<int> & aorder, 
		 const ARRAY<CliqueEl*> & cliques,
		 const ARRAY<MDOVertex> & vertices);
  ///
  void Factor (const int * blocknr);
  ///
  void FactorNew (const SparseMatrix<TM> & a);
  ///
  virtual void Mult (const BaseVector & x, BaseVector & y) const;

  ///
  virtual ostream & Print (ostream & ost) const;

  virtual void MemoryUsage (ARRAY<MemoryUsageStruct*> & mu) const
  {
    mu.Append (new MemoryUsageStruct ("SuperLU_DIST_Inverse", nnz_local_tm*sizeof(TM), 1));
  }

  ///
  void Set (int i, int j, const TM & val);
  ///
  const TM & Get (int i, int j) const;

  virtual BaseVector * CreateVector () const
  {
    ParallelDofs * aparalleldofs = const_cast<ParallelDofs * > (this->paralleldofs);
    return new ParallelVVector<TV> (aparalleldofs->GetNDof() , aparalleldofs );
  }

  bool IsSuperLUDof ( int dof ) const;
  bool IsSuperLUMasterDof ( int dof ) const;
  int SuperLUIndex ( int dof ) const;
  int SuperLUMasterIndex ( int dof ) const;



};



#endif
#endif
