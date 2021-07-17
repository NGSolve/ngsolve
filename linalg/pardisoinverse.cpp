/* *************************************************************************/
/* File:   pardisoinverse.cpp                                              */
/* Author: Florian Bachinger                                               */
/* Date:   Feb. 2004                                                       */
/* *************************************************************************/

#include <la.hpp>
#include "pardisoinverse.hpp"

#ifdef USE_MKL
#include <mkl_service.h>
#endif // USE_MKL

// #include "/opt/intel/mkl/include/mkl_service.h"

#define F77_FUNC(func)  func ## _

using ngbla::integer;

extern "C" 
{
  /* PARDISO prototype. */

#ifdef USE_PARDISO

#ifdef USE_MKL
#define MKL_PARDISO

  void mkl_free_buffers (void);
  void F77_FUNC(pardiso)
    (void * pt, integer * maxfct, integer * mnum, integer * mtype, integer * phase, integer * n, 
     double * a, integer * ia, integer * ja, integer * perm, integer * nrhs, integer * iparam, 
     integer * msglvl, double * b, double * x, integer * error);

#else // USE_MKL

#ifdef USE_PARDISO400
extern  integer F77_FUNC(pardisoinit)
    (void *, integer *, integer *, integer *, double *, integer *);
#else // USE_PARDISO400
extern  integer F77_FUNC(pardisoinit)
    (void *, integer *, integer *);
#endif // USE_PARDISO400
  integer F77_FUNC(pardiso)
    (void * pt, integer * maxfct, integer * mnum, integer * mtype, integer * phase, integer * n, 
     double * a, integer * ia, integer * ja, integer * perm, integer * nrhs, integer * iparam, 
     integer * msglvl, double * b, double * x, integer * error);

#endif

#else // USE_PARDISO
  // Neither MKL nor PARDISO linked at compile-time
  // check for MKL at run-time and set function pointers if available
#define MKL_PARDISO

  void (*mkl_free_buffers) (void) = nullptr;
  void (*F77_FUNC(pardiso))
    (void * pt, integer * maxfct, integer * mnum, integer * mtype, integer * phase, integer * n,
     double * a, integer * ia, integer * ja, integer * perm, integer * nrhs, integer * iparam,
     integer * msglvl, double * b, double * x, integer * error) = nullptr;

#endif // USE_PARDISO
}


#ifdef USE_MKL
namespace ngstd
{
    extern int mkl_max_threads;
}
#endif // USE_MKL



namespace ngla
{
#ifdef USE_PARDISO
  bool is_pardiso_available = true;
#else
  static SharedLibrary libmkl;
  static bool LoadMKL()
  {
      try
      {
#ifdef WIN32
          libmkl.Load("mkl_rt.dll");
#else
          libmkl.Load("libmkl_rt.so");
#endif
          mkl_free_buffers = libmkl.GetFunction<decltype(mkl_free_buffers)>("mkl_free_buffers");
          F77_FUNC(pardiso) = libmkl.GetFunction<decltype(pardiso_)>("pardiso_");
          return mkl_free_buffers && F77_FUNC(pardiso);
      }
      catch(const std::runtime_error &)
      {
          return false;
      }
  };
  bool is_pardiso_available = LoadMKL();
#endif

  int pardiso_msg = 0;


  class SubsetAll 
  {
  public:
    bool Used(int i) { return true; }
    bool Used(int i, int j) { return true; }
  };

  class SubsetFree
  {
    const BitArray & free;
  public:
    SubsetFree (const BitArray & afree) : free(afree) { ; }
    bool Used(int i) { return free.Test(i); }
    bool Used(int i, int j) { return free.Test(i) && free.Test(j); }
  };

  class SubsetCluster
  {
    const Array<int> & cluster;
  public:
    SubsetCluster (const Array<int> & acluster) : cluster(acluster) { ; }
    bool Used(int i) { return cluster[i]; }
    bool Used(int i, int j) { return cluster[i] && (cluster[i] == cluster[j]); }
  };


  template<class TM>
  PardisoInverseTM<TM> :: 
  PardisoInverseTM (shared_ptr<const SparseMatrixTM<TM>> a,
		    shared_ptr<BitArray> ainner,
		    shared_ptr<const Array<int>> acluster,
		    int asymmetric)
    : SparseFactorization (a, ainner, acluster)
  { 
    static Timer timer("Pardiso Inverse");
    RegionTimer reg (timer);
    GetMemoryTracer().SetName("PardisoInverseTM<" + Demangle(typeid(TM).name()) + ">");
    GetMemoryTracer().Track(rowstart, "rowstart",
                            indices, "indices",
                            matrix, "matrix");


    if (getenv ("PARDISOMSG"))
      pardiso_msg = 1;

    /*
    cout << "set mkl_rt interface" << endl;
#ifdef MKL_ILP64
    MKL_Set_Interface_Layer(MKL_INTERFACE_ILP64);
#else
    MKL_Set_Interface_Layer(MKL_INTERFACE_LP64);
#endif
    MKL_Set_Threading_Layer (MKL_THREADING_GNU);
    */

    print = bool (pardiso_msg); // false;

    symmetric = asymmetric;
    compressed = false;

    (*testout) << "Pardiso, symmetric = " << symmetric << endl;

    if (inner && cluster)
      throw Exception("PardisoInverse: Cannot use inner and cluster");

    if ( (inner && inner->Size() < a->Height()) ||
	 (cluster && cluster->Size() < a->Height() ) )
      {
	cout << "PardisoInverse: Size of inner/cluster does not match matrix size!" << endl;
	throw Exception("Invalid parameters inner/cluster. Thrown by PardisoInverse.");
      }

    if ( int( mat_traits<TM>::WIDTH) != int(mat_traits<TM>::HEIGHT) )
      {
	cout << "PardisoInverse: Each entry in the square matrix has to be a square matrix!" << endl;
	throw Exception("No Square Matrix. Thrown by PardisoInverse.");
      }


    entrysize = mat_traits<TM>::HEIGHT; 
    height = a->Height() * entrysize;

    *testout << "matrix.InverseTpye = " <<  a->GetInverseType() << endl;
    spd = ( a->GetInverseType() == PARDISOSPD ) ? 1 : 0;

    integer maxfct = 1, mnum = 1, phase = 12, nrhs = 1, msglevel = print, error = 0;
    integer * params = const_cast <integer*> (&hparams[0]);

    for (int i = 0; i < 64; i++)
      params[i] = 0;

    params[0] = 1; // no pardiso defaults
    // params[2] = 8; // 1 processor (?)
    params[2] = TaskManager::GetMaxThreads(); 

    params[1] = 0; // fill in 0..MDO, 2..metis
    params[3] = params[4] = params[5] = params[7] = params[8] = 
      params[11] = params[12] = params[18] = 0;
    params[6] = 16;
    params[9] = 13;  // perturbation 1E-10
    params[10] = 1;
    params[12] = (symmetric) ? 0 : 1;  // scling + matching
    // JS
    params[6] = 0;
    params[17] = -1;
    params[20] = 1;  // 1x1 and 2x2 bunc and Kaufman pivoting
    
    params[26] = 1; // check input matrix
    params[59] = 0; // 0..incore, 1..depending on MKL_PARDISO_MAX_CORE_SIZE, 2..ooc

    for (int i = 0; i < 128; i++) pt[i] = 0;

#ifdef MKL_PARDISO
    //    no init in MKL PARDISO
#else

    integer retvalue;
#ifdef USE_PARDISO400
    double dparm[64]; 
    integer solver = 0;
    F77_FUNC(pardisoinit) (pt,  &matrixtype, &solver, params, dparm, &retvalue); 
    cout << "pardisoinit, retval = " << retvalue << endl;
#else
    retvalue = F77_FUNC(pardisoinit) (pt,  &matrixtype, params); 
#endif
    // cout << "init success" << endl;
    // cout << "retvalue = " << retvalue << endl;
#endif
    
    SetMatrixType();

    if (inner)
      GetPardisoMatrix (*a, SubsetFree (*inner));
    else if (cluster)
      GetPardisoMatrix (*a, SubsetCluster (*cluster));
    else
      GetPardisoMatrix (*a, SubsetAll());

    nze = rowstart[compressed_height];


    // call pardiso for factorization:
    // time1 = clock();
    cout << IM(3) << "call pardiso ..." << flush;

    if (task_manager) task_manager -> StopWorkers();

#ifdef USE_MKL
    mkl_set_num_threads(mkl_max_threads);
#endif // USE_MKL

    // retvalue =
    if (matrix.Size() > 0)
      F77_FUNC(pardiso) ( pt, &maxfct, &mnum, &matrixtype, &phase, &compressed_height, 
			reinterpret_cast<double *>(matrix.Data()),
			rowstart.Data(), indices.Data(), NULL, &nrhs, params, &msglevel,
			NULL, NULL, &error );
#ifdef USE_MKL
    mkl_set_num_threads(1);
#endif // USE_MKL
    
    if (task_manager) task_manager -> StartWorkers();

    cout << IM(3) << " done" << endl;

    if ( error != 0 )
      {
	cout << IM(1) << "Setup and Factorization: PARDISO returned error " << error << "!" << endl;
	
	string errmsg;
	switch (error)
	  {
	  case -1: errmsg = "input inconsistent"; break;
	  case -2: errmsg = "not enough memory"; break;
	  case -3: errmsg = "reordering problem"; break;
	  case -4: errmsg = "zero pivot, numerical factorization or iterative refinement problem"; break;
	  case -5: errmsg = "unclassified (internal) error"; break;
	  case -6: errmsg = "preordering failed"; break;
	  default: ;
	  }
	
	cout << "err = " << errmsg << endl;

	switch (error)
	  {
	  case -4: 
	    {
	      cout << "iparam(20) = " << params[19] << endl; break;
	    }
	  default: ;
	  }
	
	cout << "symmetric = " << symmetric << endl;
	cout << "spd = " << spd << endl;
	cout << "compressed = " << compressed << endl;
	cout << "inner = " << inner << endl;
	cout << "cluster = " << cluster << endl;
	
	if (compressed_height < 1000)
	  {
	    ofstream err("pardiso.err");
	    err << "ngsolve-matrix = " << endl << a << endl;
	    err << "pardiso matrix = " << endl;
	    for (int i = 0; i < compressed_height; i++)
	      {
		err << "Row " << i << " start " << rowstart[i] << ": ";
		if ( inner ) err << " free=" << inner->Test(i) << " ";
		if ( cluster ) err << " cluster=" << (*cluster)[i] << " ";
		for (int j = rowstart[i]; j < rowstart[i+1]; j++)
		  err << "c=" << indices[j-1]-1 << ", v=" << matrix[j-1] << "   ";
		err << "\n";
	      }
	    cout << "wrote matrix to file 'pardiso.err', please check" << endl;
	  }
	throw Exception("PardisoInverse: Setup and Factorization failed.");
      }

    memory_allocated_in_pardiso_lib = 1024*params[15];
    GetMemoryTracer().Alloc(memory_allocated_in_pardiso_lib);

    /*
    (*testout) << endl << "Direct Solver: PARDISO by Schenk/Gaertner." << endl;
    (*testout) << "Matrix prepared for PARDISO in " <<
      double(time1 - starttime)/CLOCKS_PER_SEC << " sec." << endl;
    (*testout) << "Factorization by PARDISO done in " << 
      double(time2 - time1)/CLOCKS_PER_SEC << " sec." << endl << endl;
    */
  }
  
  
  template<class TM> template <typename TSUBSET>
  void PardisoInverseTM<TM> :: 
  GetPardisoMatrix (const SparseMatrixTM<TM> & a, TSUBSET subset)
  {
    Array<int> icompress (a.Height());
    icompress = -1;
    
    compress.SetSize(0);

    int cnt = 0;
    for (int i = 0; i < a.Height(); i++)
      if (subset.Used(i))
	{
	  icompress[i] = cnt++;
	  compress.Append(i);
	}
    compressed = true;
    compressed_height = cnt * entrysize;

    rowstart.SetSize (compressed_height+1);
    rowstart = 0;
    if (symmetric)
      {
	// --- transform lower left to upper right triangular matrix ---
	// 1.) build array 'rowstart':
	// (a) get nr. of entries for each row

	for (int i = 0; i < a.Height(); i++)
	  {
	    if (!subset.Used(i)) continue;

	    FlatArray<int> ind = a.GetRowIndices(i);
	    for (int j = 0; j < ind.Size(); j++ )
	      {
		if (!subset.Used(i,ind[j])) continue;
		int ccol = entrysize*icompress[ind[j]];

		if (i != ind[j])
		  for (int k = 0; k < entrysize; k++ )
		    rowstart[ccol+k+1] += entrysize;
		else
		  for (int k = 0; k < entrysize; k++ )
		    rowstart[ccol+k+1] += entrysize-k;
	      }
	  }
	
	// (b) accumulate
	rowstart[0] = 0;
	for (int i = 1; i <= compressed_height; i++) 
	  rowstart[i] += rowstart[i-1];
	

	indices.SetSize(rowstart[compressed_height]);
	matrix.SetSize (rowstart[compressed_height]);

	// 2.) build whole matrix:
	Array<int> counter(compressed_height);
	counter = 0;
	for (int i = 0; i < a.Height(); i++ )
	  {
	    if (!subset.Used(i)) continue;
	    int ci = entrysize*icompress[i];

	    FlatArray<int> ind = a.GetRowIndices(i);
	    FlatVector<TM> values = a.GetRowValues(i);

	    for (int j = 0; j < ind.Size(); j++)
	      {
		if (!subset.Used (i,ind[j])) continue;

		int ccol = entrysize*icompress[ind[j]];
		for (int k = 0; k < entrysize; k++)
		  for (int l = 0; l < entrysize; l++ )
		    if ( (i != ind[j]) || (k <= l) )
		      {
			indices[rowstart[ccol+k]+counter[ccol+k]] = ci+l+1;
			matrix[rowstart[ccol+k]+counter[ccol+k]] = Access(values[j],l,k);
			counter[ccol+k]++;
		      }
	      }	
	  }

	for (int i = 0; i <= compressed_height; i++)
	  rowstart[i]++;
      }
    else
      {
	int counter = 0;
	for (int i = 0; i < a.Height(); i++ )
	  {
	    if (!subset.Used(i)) continue;
	    FlatArray<int> ind = a.GetRowIndices(i);

	    int rowelems = 0;
	    for (int j = 0; j < ind.Size(); j++)
	      if (subset.Used(i,ind[j]))
		rowelems += entrysize;
	    
	    int ci = entrysize * icompress[i];
	    for (int k = 0; k < entrysize; k++, counter += rowelems )
	      rowstart[ci+k] = counter+1;
	  }
	rowstart[compressed_height] = counter+1;
	
	indices.SetSize(counter);
	matrix.SetSize (counter);

	for (int i = 0; i < a.Height(); i++ )
	  {
	    if (!subset.Used(i)) continue;
	    FlatArray<int> ind = a.GetRowIndices(i);
	    FlatVector<TM> values = a.GetRowValues(i);

	    int ci = entrysize*icompress[i];
	    
	    int counter = 0;
	    for (int j = 0; j < ind.Size(); j++ )
	      if (subset.Used(i,ind[j]))
		{
		  int ccol = entrysize*icompress[ind[j]];
		  for (int k=0; k<entrysize; k++ )
		    for (int l=0; l<entrysize; l++ )
		      {
			indices[rowstart[ci+k]+counter+l-1] = ccol+l+1;
			matrix[rowstart[ci+k]+counter+l-1] = Access(values[j],k,l);
		      }
		  counter+=entrysize;
		}
	  }
      }
  }




  template<class TM, class TV_ROW, class TV_COL>
  void PardisoInverse<TM,TV_ROW,TV_COL> :: 
  Mult (const BaseVector & x, BaseVector & y) const
  {
    static Timer timer(string("Pardiso Solve, mat = ") + typeid(TM).name() + ", vec = " + typeid(TV_ROW).name());
    RegionTimer reg (timer);

    FlatVector<TVX> fx = x.FV<TVX> ();
    FlatVector<TVX> fy = y.FV<TVX> ();

    integer maxfct = 1, mnum = 1, phase = 33, msglevel = 0, error = 0;
    integer nrhs = fx.Size() / (height/entrysize);
    
    if (fx.Size() != fy.Size()) 
      {
	cout << "PardisoInverse::Mult .. sizes don't match" << endl;
	cout << "type<TVX> = " << typeid(TVX).name() << endl;
	cout << "type<TM> = " << typeid(TM).name() << endl;
	cout << "fx.size = " << fx.Size() << endl;
	cout << "fy.size = " << fy.Size() << endl;
	cout << "size(x) = " << x.Size() << endl;
	cout << "size(y) = " << y.Size() << endl;
	cout << "height = " << height/entrysize << endl;
      }

    FlatMatrix<TVX> mx(nrhs, height/entrysize, (TVX*)fx.Data());
    FlatMatrix<TVX> my(nrhs, height/entrysize, (TVX*)fy.Data());

    integer * params = const_cast <integer*> (&hparams[0]);

    if(task_manager)
      task_manager->SuspendWorkers(1000);

#ifdef USE_MKL
    mkl_set_num_threads(mkl_max_threads);
#endif // USE_MKL

    if (matrix.Size() > 0)
      {
        if (compressed)
          {
            Matrix<TVX> hx(nrhs, compress.Size());
            Matrix<TVX> hy(nrhs, compress.Size());
            hx = mx.Cols(compress);
            F77_FUNC(pardiso) ( const_cast<integer*>(pt), &maxfct, &mnum, 
                                const_cast<integer*>(&matrixtype),
                                &phase, const_cast<integer*>(&compressed_height), 
                                reinterpret_cast<double *>(matrix.Data()),
                                rowstart.Data(), indices.Data(),
                                NULL, &nrhs, params, &msglevel,
                                reinterpret_cast<double *>(hx.Data()), 
                                reinterpret_cast<double *>(hy.Data()), &error );
            my = 0; 
            my.Cols(compress) = hy;
          }
        else
          {
            F77_FUNC(pardiso) ( const_cast<integer *>(pt), &maxfct, &mnum, 
                                const_cast<integer *>(&matrixtype),
                                &phase, const_cast<integer *>(&compressed_height), 
                                reinterpret_cast<double *>(matrix.Data()),
                                rowstart.Data(), indices.Data(),
                                NULL, &nrhs, params, &msglevel,
                                static_cast<double *>((void*)fx.Data()), 
                                static_cast<double *>((void*)fy.Data()), &error );
          }
      }

#ifdef USE_MKL
    mkl_set_num_threads(1);
#endif // USE_MKL

    if(task_manager)
      task_manager->ResumeWorkers();

    if ( error != 0 )
      cout << "Apply Inverse: PARDISO returned error " << error << "!" << endl;
  }
  
  
  template<class TM, class TV_ROW, class TV_COL>
  void PardisoInverse<TM,TV_ROW,TV_COL> ::
  MultTrans (const BaseVector & x, BaseVector & y) const
  {
    const_cast<integer&>(hparams[11]) = 2; // Solve transposed matrix

    // See https://software.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/sparse-solver-routines/onemkl-pardiso-parallel-direct-sparse-solver-interface/pardiso-iparm-parameter.html
    
    Mult(x,y);
    const_cast<integer&>(hparams[11]) = 0;
  }




  template<>
  void PardisoInverse<double,Complex,Complex> :: 
  Mult (const BaseVector & x, BaseVector & y) const
  {
    static Timer timer ("Pardiso Solve, mat = d, vec = c");
    RegionTimer reg (timer);


    FlatMatrixFixWidth<2> hx(x.Size(), (double*)x.Memory());
    Matrix<> tx(2, compressed_height);
    if (compressed)
      tx = Trans (hx.Rows(compress));
      // tx = Trans (hx).Cols(compress));
    // for (int i = 0; i < compress.Size(); i++)
    // tx.Col(i) = hx.Row(compress[i]);
    else
      tx = Trans (hx);

    FlatMatrixFixWidth<2> hy(y.Size(), (double*)y.Memory());
    Matrix<> ty(2, compressed_height);

      
    integer maxfct = 1, mnum = 1, phase = 33, msglevel = 0, error = 0;
    integer nrhs = 2;

    integer * params = const_cast <integer*> (&hparams[0]);

    if(task_manager)
      task_manager->SuspendWorkers(1000);

#ifdef USE_MKL
    mkl_set_num_threads(mkl_max_threads);
#endif // USE_MKL

    if (matrix.Size() > 0)
      F77_FUNC(pardiso) ( const_cast<integer *>(pt), 
                          &maxfct, &mnum, const_cast<integer *>(&matrixtype),
                          &phase, const_cast<integer *>(&compressed_height), 
                          reinterpret_cast<double *>(matrix.Data()),
                          rowstart.Data(), indices.Data(),
                          NULL, &nrhs, params, &msglevel, tx.Data(), ty.Data(),
                          &error );
#ifdef USE_MKL
      mkl_set_num_threads(1);
#endif // USE_MKL

    if(task_manager)
      task_manager->ResumeWorkers();

    if ( error != 0 )
      cout << "Apply Inverse: PARDISO returned error " << error << "!" << endl;

    if (compressed)
      {
	hy = 0;
	hy.Rows(compress) = Trans(ty);
      }
    else
      hy = Trans (ty);
  }
  
  






  template<class TM>
  ostream & PardisoInverseTM<TM> :: Print (ostream & ost) const
  {
    cout << "PardisoInverse::Print not implemented!" << endl;
    return ost; 
  }


  template<class TM>
  PardisoInverseTM<TM> :: ~PardisoInverseTM()
  {
    integer maxfct = 1, mnum = 1, phase = -1, nrhs = 1, msglevel = 1, error;
    integer * params = const_cast <integer*> (&hparams[0]);

    //    cout << "call pardiso (clean up) ..." << endl;
    if (task_manager) task_manager -> StopWorkers();
    F77_FUNC(pardiso) ( pt, &maxfct, &mnum, &matrixtype, &phase, &compressed_height, NULL,
			rowstart.Data(), indices.Data(), NULL, &nrhs, params, &msglevel,
			NULL, NULL, &error );
#ifdef MKL_PARDISO
    mkl_free_buffers();
#endif // MKL_PARDISO
    GetMemoryTracer().Free(memory_allocated_in_pardiso_lib);
    memory_allocated_in_pardiso_lib = 0;
    if (task_manager) task_manager -> StartWorkers();
    if (error != 0)
      cout << "Clean Up: PARDISO returned error " << error << "!" << endl;
  }






  template<class TM>
  void PardisoInverseTM<TM> :: SetMatrixType() // TM entry)
  {
    if (mat_traits<TM>::IS_COMPLEX)
      {
	if ( symmetric ) 
	  matrixtype = 6;   // complex symmetric
	else
	  matrixtype = 13;  // complex genral
      }
    else
      {
	if ( symmetric )
	  {
	    if ( spd ) 
	      matrixtype = 2;    // pos def
	    else
	      matrixtype = -2;   // symmetric indef
	  }
	else 
	  matrixtype = 11;       // general non-sym
      }

    if (print)
      cout << "spd = " << int(spd) << ", sym = " << int(symmetric) 
	   << ", complex = " << int(mat_traits<TM>::IS_COMPLEX)
	   << ", matrixtype = " << matrixtype << endl;
    *testout << "pardiso matrixtype = " << matrixtype << endl;
  }







  template class PardisoInverseTM<double>;
  template class PardisoInverseTM<Complex>;
  template class PardisoInverse<double>;
  template class PardisoInverse<Complex>;
  template class PardisoInverse<double,Complex,Complex>;
#if MAX_SYS_DIM >= 1
  template class PardisoInverseTM<Mat<1,1,double> >;
  template class PardisoInverseTM<Mat<1,1,Complex> >;
  template class PardisoInverse<Mat<1,1,double> >;
  template class PardisoInverse<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class PardisoInverseTM<Mat<2,2,double> >;
  template class PardisoInverseTM<Mat<2,2,Complex> >;
  template class PardisoInverse<Mat<2,2,double> >;
  template class PardisoInverse<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class PardisoInverseTM<Mat<3,3,double> >;
  template class PardisoInverseTM<Mat<3,3,Complex> >;
  template class PardisoInverse<Mat<3,3,double> >;
  template class PardisoInverse<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class PardisoInverseTM<Mat<4,4,double> >;
  template class PardisoInverseTM<Mat<4,4,Complex> >;
  template class PardisoInverse<Mat<4,4,double> >;
  template class PardisoInverse<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class PardisoInverse<Mat<5,5,double> >;
  template class PardisoInverse<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class PardisoInverse<Mat<6,6,double> >;
  template class PardisoInverse<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class PardisoInverse<Mat<7,7,double> >;
  template class PardisoInverse<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class PardisoInverse<Mat<8,8,double> >;
  template class PardisoInverse<Mat<8,8,Complex> >;
#endif





#if MAX_CACHEBLOCKS >= 2
  template class PardisoInverse<double, Vec<2,double>, Vec<2,double> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class PardisoInverse<double, Vec<3,double>, Vec<3,double> >;
  template class PardisoInverse<double, Vec<4,double>, Vec<4,double> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class PardisoInverse<double, Vec<5,double>, Vec<5,double> >;
  template class PardisoInverse<double, Vec<6,double>, Vec<6,double> >;
  template class PardisoInverse<double, Vec<7,double>, Vec<7,double> >;
  template class PardisoInverse<double, Vec<8,double>, Vec<8,double> >;
  template class PardisoInverse<double, Vec<9,double>, Vec<9,double> >;
  template class PardisoInverse<double, Vec<10,double>, Vec<10,double> >;
  template class PardisoInverse<double, Vec<11,double>, Vec<11,double> >;
  template class PardisoInverse<double, Vec<12,double>, Vec<12,double> >;
  template class PardisoInverse<double, Vec<13,double>, Vec<13,double> >;
  template class PardisoInverse<double, Vec<14,double>, Vec<14,double> >;
  template class PardisoInverse<double, Vec<15,double>, Vec<15,double> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class PardisoInverse<double, Vec<2,Complex>, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class PardisoInverse<double, Vec<3,Complex>, Vec<3,Complex> >;
  template class PardisoInverse<double, Vec<4,Complex>, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class PardisoInverse<double, Vec<5,Complex>, Vec<5,Complex> >;
  template class PardisoInverse<double, Vec<6,Complex>, Vec<6,Complex> >;
  template class PardisoInverse<double, Vec<7,Complex>, Vec<7,Complex> >;
  template class PardisoInverse<double, Vec<8,Complex>, Vec<8,Complex> >;
  template class PardisoInverse<double, Vec<9,Complex>, Vec<9,Complex> >;
  template class PardisoInverse<double, Vec<10,Complex>, Vec<10,Complex> >;
  template class PardisoInverse<double, Vec<11,Complex>, Vec<11,Complex> >;
  template class PardisoInverse<double, Vec<12,Complex>, Vec<12,Complex> >;
  template class PardisoInverse<double, Vec<13,Complex>, Vec<13,Complex> >;
  template class PardisoInverse<double, Vec<14,Complex>, Vec<14,Complex> >;
  template class PardisoInverse<double, Vec<15,Complex>, Vec<15,Complex> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class PardisoInverse<Complex, Vec<2,Complex>, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class PardisoInverse<Complex, Vec<3,Complex>, Vec<3,Complex> >;
  template class PardisoInverse<Complex, Vec<4,Complex>, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class PardisoInverse<Complex, Vec<5,Complex>, Vec<5,Complex> >;
  template class PardisoInverse<Complex, Vec<6,Complex>, Vec<6,Complex> >;
  template class PardisoInverse<Complex, Vec<7,Complex>, Vec<7,Complex> >;
  template class PardisoInverse<Complex, Vec<8,Complex>, Vec<8,Complex> >;
  template class PardisoInverse<Complex, Vec<9,Complex>, Vec<9,Complex> >;
  template class PardisoInverse<Complex, Vec<10,Complex>, Vec<10,Complex> >;
  template class PardisoInverse<Complex, Vec<11,Complex>, Vec<11,Complex> >;
  template class PardisoInverse<Complex, Vec<12,Complex>, Vec<12,Complex> >;
  template class PardisoInverse<Complex, Vec<13,Complex>, Vec<13,Complex> >;
  template class PardisoInverse<Complex, Vec<14,Complex>, Vec<14,Complex> >;
  template class PardisoInverse<Complex, Vec<15,Complex>, Vec<15,Complex> >;
#endif



}

