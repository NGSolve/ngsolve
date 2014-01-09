/* *************************************************************************/
/* File:   superluinverse.cpp                                              */
/* Author: Florian Bachinger                                               */
/* Date:   Feb. 2004                                                       */
/* *************************************************************************/

#include <la.hpp>


#ifdef USE_SUPERLU




namespace ngla
{
  using namespace ngla;
  using namespace ngstd;
  
  using namespace double_superlu;
  using namespace complex_superlu;

namespace superlufunc

{
  using namespace superlufunc;
  template<class TM>
  int NumRows(TM entry)
  {
    return entry.Height();
  }
  int NumRows(double entry) { return 1; }
  int NumRows(Complex entry) { return 1; }
  
  
  template<class TM>
  int NumCols(TM entry)
  {
    return entry.Width();
  }
  int NumCols(double entry) { return 1; }
  int NumCols(Complex entry) { return 1; }
  
  
  template<class TM>
  typename TM::TSCAL Elem(TM entry, int i, int j)
  {
    return entry(i,j);
  }
  double Elem (double entry, int i, int j) { return entry; }
  Complex Elem (Complex entry, int i, int j) { return entry; }
  
  
  
  template<class TM>
  int IsComplex(TM entry)
  {
    return ( sizeof( Elem(entry,0,0) ) == sizeof( Complex ) );
  }
  
}

  using namespace superlufunc;


  template <class TM, class TV_ROW, class TV_COL>
  SuperLUInverse<TM,TV_ROW,TV_COL> :: 
  SuperLUInverse (const SparseMatrix<TM,TV_ROW,TV_COL> & a, 
		  const BitArray * ainner,
		  const Array<int> * acluster,
		  int asymmetric)
  { 
    //        (*testout) << "matrix = " << a << endl;

    symmetric = asymmetric;
    inner = ainner;
    cluster = acluster;
    cout << "SuperLU called ...";

    if ( (inner && (inner->Size() < a.Height())) ||
         (cluster && (cluster->Size() < a.Height())) )
      {
	cout << "SuperLU: Size of inner/cluster does not match matrix size!" << endl;
	throw Exception("Invalid parameters inner/cluster. Thrown by PardisoInverse.");
      }

    clock_t starttime, time1, time2;
    starttime = clock();


    // prepare matrix and parameters for SuperLU
    TM entry = a(0,0);
    if ( NumRows(entry) != NumCols(entry) )
      {
	cout << "SuperLU: Each entry in the square matrix has to be a square matrix!" << endl;
	throw Exception("No Square Matrix. Thrown by SuperLUInverse.");
      }
    entrysize = NumRows(entry);
    iscomplex = IsComplex<TM>(entry);

    height = a.Height() * entrysize;
    colstart = new int[height+1];
    int * counter = new int[height];

    int i, j, k, l;
    int col;
 
    for ( i=0; i<height; i++ ) counter[i] = colstart[i+1] = 0;

    
    if ( symmetric )
      {
	// --- transform lower left to full matrix (in compressed column storage format) ---
	// 1.) build array 'colstart':
	// (a) get nr. of entries for each col
	for ( i=1; i<=a.Height(); i++ )
	  {
	    for ( j=0; j<a.GetRowIndices(i-1).Size(); j++ )
	      {
		col = a.GetRowIndices(i-1)[j];

		if (  (!inner && !cluster) ||
		      (inner && (inner->Test(i-1) && inner->Test(col) ) ) ||
		      (!inner && cluster && 
		           ((*cluster)[i-1] == (*cluster)[col] 
			    && (*cluster)[i-1] ))  )
		  {
		    for ( k=0; k<entrysize; k++ )
		      {
			colstart[col*entrysize+k+1] += entrysize;
			if (i-1 != col) colstart[(i-1)*entrysize+k+1] += entrysize;
		      }
		  }
		else if ( i-1 == col )
		  {
		    for ( k=0; k<entrysize; k++ )
		      colstart[col*entrysize+k+1] ++;
		  }
	      }
	  }
	// (b) accumulate
	colstart[0] = 0;
	for ( i=1; i<=height; i++ ) colstart[i] += colstart[i-1];
	nze = colstart[height];
	


	// 2.) build whole matrix:
 	indices = new int[nze];
 	matrix = new TSCAL[nze];      

	for ( i=0; i<a.Height(); i++ )
	  {
	    for ( j=0; j<a.GetRowIndices(i).Size(); j++ )
	      {
		col = a.GetRowIndices(i)[j];

		if (  (!inner && !cluster) ||
		      (inner && (inner->Test(i) && inner->Test(col) ) ) ||
		      (!inner && cluster && 
		           ((*cluster)[i] == (*cluster)[col] 
			    && (*cluster)[i] ))  )
		  {
		    entry = a(i,col);
		    for ( k=0; k<entrysize; k++)
		      for ( l=0; l<entrysize; l++ )
			{
			  indices[ colstart[col*entrysize+k]+
				   counter[col*entrysize+k] ] = i*entrysize+l;
			  matrix[ colstart[col*entrysize+k]+
				  counter[col*entrysize+k] ] = Elem(entry,l,k);
			  counter[col*entrysize+k]++;

			  if ( i != col )
			    {
			      indices[ colstart[i*entrysize+l]+
				       counter[i*entrysize+l] ] = col*entrysize+k;
			      matrix[ colstart[i*entrysize+l]+
				      counter[i*entrysize+l] ] = Elem(entry,k,l);
			      counter[i*entrysize+l]++;
			    }
			}
		  }
		else if (i == col)
		  {
		    // in the case of 'inner' or 'cluster': 1 on the diagonal for
		    // unused dofs.
		    for ( l=0; l<entrysize; l++ )
		      {
			indices[ colstart[col*entrysize+l]+
				 counter[col*entrysize+l] ] = col*entrysize+l;
			matrix[ colstart[col*entrysize+l]+
				counter[col*entrysize+l] ] = 1;
			counter[col*entrysize+l]++;
		      }
		  }
	      }
	  }
	

	
// 	(*testout) << endl << "col, colstart / indices, matrix-entries" << endl;
// 	for ( i=0; i<height; i++ )
// 	  {
// 	    (*testout) << endl << i+1 << ", " << colstart[i] << ":   ";
// 	    for ( j=colstart[i]; j<colstart[i+1]; j++ )
// 	      (*testout) << indices[j] << ", " << matrix[j] << "      ";
// 	  }
	
      }
    else
      {
	// --- transform matrix to compressed column storage format ---
	// 1.) build array 'colstart':
	// (a) get nr. of entries for each col
	for ( i=1; i<=a.Height(); i++ )
	  {
	    for ( j=0; j<a.GetRowIndices(i-1).Size(); j++ )
	      {
		col = a.GetRowIndices(i-1)[j];

		if (  (!inner && !cluster) ||
		      (inner && (inner->Test(i-1) && inner->Test(col) ) ) ||
		      (!inner && cluster && 
		           ((*cluster)[i-1] == (*cluster)[col] 
			    && (*cluster)[i-1] ))  )
		  {
		    for ( k=0; k<entrysize; k++ )
			colstart[col*entrysize+k+1] += entrysize;
		  }
		else if ( i-1 == col )
		  {
		    for ( k=0; k<entrysize; k++ )
		      colstart[col*entrysize+k+1] ++;
		  }
	      }
	  }
	// (b) accumulate
	colstart[0] = 0;
	for ( i=1; i<=height; i++ ) colstart[i] += colstart[i-1];
	nze = colstart[height];


	// 2.) build whole matrix:
	indices = new int[a.NZE() * entrysize * entrysize ];
	matrix = new TSCAL[a.NZE() * entrysize * entrysize ];      


	for ( i=0; i<a.Height(); i++ )
	  {
	    for ( j=0; j<a.GetRowIndices(i).Size(); j++ )
	      {
		col = a.GetRowIndices(i)[j];

		if (  (!inner && !cluster) ||
		      (inner && (inner->Test(i) && inner->Test(col) ) ) ||
		      (!inner && cluster && 
		           ((*cluster)[i] == (*cluster)[col] 
			    && (*cluster)[i] ))  )
		  {
		    entry = a(i,col);
		    for ( k=0; k<entrysize; k++)
		      for ( l=0; l<entrysize; l++ )
			{
			  indices[ colstart[col*entrysize+k]+
				   counter[col*entrysize+k] ] = i*entrysize+l;
			  matrix[ colstart[col*entrysize+k]+
				  counter[col*entrysize+k] ] = Elem(entry,l,k);
			  counter[col*entrysize+k]++;
			}
		  }
		else if (i == col)
		  {
		    // in the case of 'inner' or 'cluster': 1 on the diagonal for
		    // unused dofs.
		    for ( l=0; l<entrysize; l++ )
		      {
			indices[ colstart[col*entrysize+l]+
				 counter[col*entrysize+l] ] = col*entrysize+l;
			matrix[ colstart[col*entrysize+l]+
				counter[col*entrysize+l] ] = 1;
			counter[col*entrysize+l]++;
		      }
		  }
	      }
	  }

	
// 	   (*testout) << endl << "col, colstart / indices, matrix-entries" << endl;
// 	   for ( i=0; i<height; i++ )
// 	   {
// 	   (*testout) << endl << i << ", " << colstart[i] << ":   ";
// 	   for ( j=colstart[i]; j<colstart[i+1]; j++ )
// 	   (*testout) << indices[j] << ", " << matrix[j] << "      ";
// 	   }
// 	   (*testout) << endl;
	
      }

    


    // todo: rather only factor, not factor and solve
    TSCAL *rhs = new TSCAL[height];
    for ( i=0; i<height; i++ ) rhs[i] = 0;

    superlu_options_t options;
    int               error;

    perm_c = new int[height];
    perm_r = new int[height];

    set_default_options(&options);
    StatInit(&stat);
    

    // call superlu for factorization:
    time1 = clock();

    
    if ( iscomplex )
      {
	zCreate_CompCol_Matrix(&A, height, height, nze, 
			       reinterpret_cast<doublecomplex *>(matrix), indices, colstart, 
			       SLU_NC, SLU_Z, SLU_GE);
	zCreate_Dense_Matrix(&B, height, 1, 
			     reinterpret_cast<doublecomplex *>(rhs), height, SLU_DN, SLU_Z, SLU_GE);

	zgssv( &options, &A, perm_c, perm_r, &L, &U, &B, &stat, &error );
      }
    else
      {
	dCreate_CompCol_Matrix(&A, height, height, nze, 
			       reinterpret_cast<double *>(matrix), indices, colstart, 
			       SLU_NC, SLU_D, SLU_GE);

	dCreate_Dense_Matrix(&B, height, 1, 
			     reinterpret_cast<double *>(rhs), height, SLU_DN, SLU_D, SLU_GE);
	
	dgssv( &options, &A, perm_c, perm_r, &L, &U, &B, &stat, &error );
      }

    time2 = clock();


    if ( error != 0 )
      {
	cout << "Setup and Factorization: SuperLU returned error " << error << "!" << endl;
	throw Exception("SuperLUInverse: Setup and Factorization failed.");
      }

    (*testout) << endl << "Direct Solver: SuperLU by Lawrence Berkeley National Laboratory." << endl;
    (*testout) << "Matrix prepared for SuperLU in " <<
      double(time1 - starttime)/CLOCKS_PER_SEC << " sec." << endl;
    (*testout) << "Factorization by SuperLU done in " << 
      double(time2 - time1)/CLOCKS_PER_SEC << " sec." << endl << endl;
    
    cout << " done " << endl;
    delete [] counter;
    delete [] rhs;    
 }
  
  

  template <class TM, class TV_ROW, class TV_COL>
  SuperLUInverse<TM,TV_ROW,TV_COL> :: 

  SuperLUInverse (const Array<int> & aorder, 
		  const Array<CliqueEl*> & cliques,
		  const Array<MDOVertex> & vertices,
		  int symmetric)
  {
    Allocate (aorder, cliques, vertices);
  }
  

  

  template <class TM, class TV_ROW, class TV_COL>
  void SuperLUInverse<TM, TV_ROW,TV_COL> :: 
  Allocate (const Array<int> & aorder, 
	    const Array<CliqueEl*> & cliques,
	    const Array<MDOVertex> & vertices)
  {
    cout << "SuperLUInverse::Allocate not implemented!" << endl;
  }
  


  template <class TM, class TV_ROW, class TV_COL>
  void SuperLUInverse<TM,TV_ROW,TV_COL> :: 
  FactorNew (const SparseMatrix<TM> & a)
  {
    throw Exception ("SuperLUInverse::FactorNew not implemented");
  }






  template <class TM, class TV_ROW, class TV_COL>
  void SuperLUInverse<TM,TV_ROW, TV_COL> :: Factor (const int * blocknr)
  {
    cout << "SuperLUInverse::Factor not implemented!" << endl;
  }
  
  



  template <class TM, class TV_ROW, class TV_COL>
  void SuperLUInverse<TM,TV_ROW,TV_COL> :: 
  Mult (const BaseVector & x, BaseVector & y) const
  {

    FlatVector<TVX> fx = x.FV<TVX> ();
    // dynamic_cast<T_BaseVector<TVX> &> (const_cast<BaseVector &> (x)).FV();
    FlatVector<TVX> fy = y.FV<TVX> ();
    // dynamic_cast<T_BaseVector<TVX> &> (y).FV();
    
    int error;
    fy = fx;
    
    if ( iscomplex )
      {
	zCreate_Dense_Matrix(const_cast<SuperMatrix *>(&B), height, 1, 
			     static_cast<doublecomplex *>(fy.Data()), 
			     height, SLU_DN, SLU_Z, SLU_GE);

	zgstrs( NOTRANS, const_cast<SuperMatrix *>(&L), const_cast<SuperMatrix *>(&U), 
		perm_c, perm_r, const_cast<SuperMatrix *>(&B), 
		const_cast<SuperLUStat_t *>(&stat), &error );
      }
    else
      {
	dCreate_Dense_Matrix(const_cast<SuperMatrix *>(&B), height, 1, 
			     static_cast<double *>(fy.Data()), 
			     height, SLU_DN, SLU_D, SLU_GE);

	dgstrs( NOTRANS, const_cast<SuperMatrix *>(&L), const_cast<SuperMatrix *>(&U), 
		perm_c, perm_r, const_cast<SuperMatrix *>(&B), 
		const_cast<SuperLUStat_t *>(&stat), &error );
      }

    if ( error != 0 )
      cout << "Apply Inverse: SuperLU returned error " << error << "!" << endl;


    if (inner)
      {
	for (int i=0; i<height/entrysize; i++)
	  if (!inner->Test(i)) 
	    for (int j=0; j<entrysize; j++ ) fy(i*entrysize+j) = 0.0;
      }
    else if (cluster)
      {
	for (int i=0; i<height/entrysize; i++)
	  if (!(*cluster)[i]) 
	    for (int j=0; j<entrysize; j++ ) fy(i*entrysize+j) = 0.0;
      }
  }
  
  




/*
  template <class TM, class TV_ROW, class TV_COL>
  void SuperLUInverse<TM,TV_ROW,TV_COL> :: Set (int i, int j, const TM & val)
  {
    cout << "SuperLUInverse::Set not implemented!" << endl;
  }



  template <class TM, class TV_ROW, class TV_COL>
  const TM & SuperLUInverse<TM,TV_ROW,TV_COL> :: Get (int i, int j) const
  {
    cout << "SuperLUInverse::Get not implemented!" << endl;
  }
*/

  template <class TM, class TV_ROW, class TV_COL>
  ostream & SuperLUInverse<TM,TV_ROW,TV_COL> :: Print (ostream & ost) const
  {
    cout << "SuperLUInverse::Print not implemented!" << endl;
    return ost; 
  }



  template <class TM, class TV_ROW, class TV_COL>
  SuperLUInverse<TM,TV_ROW,TV_COL> :: ~SuperLUInverse()
  {
//     Destroy_CompCol_Matrix(&A);
//     Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree(&stat);

    delete [] perm_c;
    delete [] perm_r;
    delete [] colstart;
    delete [] indices;
    delete [] matrix;
  }







  template class SuperLUInverse<double>;
  template class SuperLUInverse<Complex>;
  template class SuperLUInverse<double,Complex,Complex>;
#if MAX_SYS_DIM >= 1
  template class SuperLUInverse<Mat<1,1,double> >;
  template class SuperLUInverse<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class SuperLUInverse<Mat<2,2,double> >;
  template class SuperLUInverse<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class SuperLUInverse<Mat<3,3,double> >;
  template class SuperLUInverse<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class SuperLUInverse<Mat<4,4,double> >;
  template class SuperLUInverse<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class SuperLUInverse<Mat<5,5,double> >;
  template class SuperLUInverse<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class SuperLUInverse<Mat<6,6,double> >;
  template class SuperLUInverse<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class SuperLUInverse<Mat<7,7,double> >;
  template class SuperLUInverse<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class SuperLUInverse<Mat<8,8,double> >;
  template class SuperLUInverse<Mat<8,8,Complex> >;
#endif






  template class SuperLUInverse<double, Vec<2,double>, Vec<2,double> >;
  template class SuperLUInverse<double, Vec<3,double>, Vec<3,double> >;
  template class SuperLUInverse<double, Vec<4,double>, Vec<4,double> >;
  template class SuperLUInverse<double, Vec<5,double>, Vec<5,double> >;
  template class SuperLUInverse<double, Vec<6,double>, Vec<6,double> >;
  template class SuperLUInverse<double, Vec<7,double>, Vec<7,double> >;
  template class SuperLUInverse<double, Vec<8,double>, Vec<8,double> >;
  template class SuperLUInverse<double, Vec<9,double>, Vec<9,double> >;
  template class SuperLUInverse<double, Vec<10,double>, Vec<10,double> >;
  template class SuperLUInverse<double, Vec<11,double>, Vec<11,double> >;
  template class SuperLUInverse<double, Vec<12,double>, Vec<12,double> >;
  template class SuperLUInverse<double, Vec<13,double>, Vec<13,double> >;
  template class SuperLUInverse<double, Vec<14,double>, Vec<14,double> >;
  template class SuperLUInverse<double, Vec<15,double>, Vec<15,double> >;
  
  template class SuperLUInverse<double, Vec<2,Complex>, Vec<2,Complex> >;
  template class SuperLUInverse<double, Vec<3,Complex>, Vec<3,Complex> >;
  template class SuperLUInverse<double, Vec<4,Complex>, Vec<4,Complex> >;
  template class SuperLUInverse<double, Vec<5,Complex>, Vec<5,Complex> >;
  template class SuperLUInverse<double, Vec<6,Complex>, Vec<6,Complex> >;
  template class SuperLUInverse<double, Vec<7,Complex>, Vec<7,Complex> >;
  template class SuperLUInverse<double, Vec<8,Complex>, Vec<8,Complex> >;
  template class SuperLUInverse<double, Vec<9,Complex>, Vec<9,Complex> >;
  template class SuperLUInverse<double, Vec<10,Complex>, Vec<10,Complex> >;
  template class SuperLUInverse<double, Vec<11,Complex>, Vec<11,Complex> >;
  template class SuperLUInverse<double, Vec<12,Complex>, Vec<12,Complex> >;
  template class SuperLUInverse<double, Vec<13,Complex>, Vec<13,Complex> >;
  template class SuperLUInverse<double, Vec<14,Complex>, Vec<14,Complex> >;
  template class SuperLUInverse<double, Vec<15,Complex>, Vec<15,Complex> >;

  template class SuperLUInverse<Complex, Vec<2,Complex>, Vec<2,Complex> >;
  template class SuperLUInverse<Complex, Vec<3,Complex>, Vec<3,Complex> >;
  template class SuperLUInverse<Complex, Vec<4,Complex>, Vec<4,Complex> >;
  template class SuperLUInverse<Complex, Vec<5,Complex>, Vec<5,Complex> >;
  template class SuperLUInverse<Complex, Vec<6,Complex>, Vec<6,Complex> >;
  template class SuperLUInverse<Complex, Vec<7,Complex>, Vec<7,Complex> >;
  template class SuperLUInverse<Complex, Vec<8,Complex>, Vec<8,Complex> >;
  template class SuperLUInverse<Complex, Vec<9,Complex>, Vec<9,Complex> >;
  template class SuperLUInverse<Complex, Vec<10,Complex>, Vec<10,Complex> >;
  template class SuperLUInverse<Complex, Vec<11,Complex>, Vec<11,Complex> >;
  template class SuperLUInverse<Complex, Vec<12,Complex>, Vec<12,Complex> >;
  template class SuperLUInverse<Complex, Vec<13,Complex>, Vec<13,Complex> >;
  template class SuperLUInverse<Complex, Vec<14,Complex>, Vec<14,Complex> >;
  template class SuperLUInverse<Complex, Vec<15,Complex>, Vec<15,Complex> >;
  
  /*



  template class SuperLUInverse<double, Vec<2,double>, Vec<2,double> >;
  template class SuperLUInverse<double, Vec<3,double>, Vec<3,double> >;
  template class SuperLUInverse<double, Vec<4,double>, Vec<4,double> >;
  template class SuperLUInverse<double, Vec<5,double>, Vec<5,double> >;
  template class SuperLUInverse<double, Vec<6,double>, Vec<6,double> >;
  template class SuperLUInverse<double, Vec<7,double>, Vec<7,double> >;
  template class SuperLUInverse<double, Vec<8,double>, Vec<8,double> >;
  template class SuperLUInverse<double, Vec<9,double>, Vec<9,double> >;
  template class SuperLUInverse<double, Vec<10,double>, Vec<10,double> >;
  template class SuperLUInverse<double, Vec<11,double>, Vec<11,double> >;
  template class SuperLUInverse<double, Vec<12,double>, Vec<12,double> >;
  template class SuperLUInverse<double, Vec<13,double>, Vec<13,double> >;
  template class SuperLUInverse<double, Vec<14,double>, Vec<14,double> >;
  template class SuperLUInverse<double, Vec<15,double>, Vec<15,double> >;

  template class SuperLUInverse<double, Vec<2,Complex>, Vec<2,Complex> >;
  template class SuperLUInverse<double, Vec<3,Complex>, Vec<3,Complex> >;
  template class SuperLUInverse<double, Vec<4,Complex>, Vec<4,Complex> >;
  template class SuperLUInverse<double, Vec<5,Complex>, Vec<5,Complex> >;
  template class SuperLUInverse<double, Vec<6,Complex>, Vec<6,Complex> >;
  template class SuperLUInverse<double, Vec<7,Complex>, Vec<7,Complex> >;
  template class SuperLUInverse<double, Vec<8,Complex>, Vec<8,Complex> >;
  template class SuperLUInverse<double, Vec<9,Complex>, Vec<9,Complex> >;
  template class SuperLUInverse<double, Vec<10,Complex>, Vec<10,Complex> >;
  template class SuperLUInverse<double, Vec<11,Complex>, Vec<11,Complex> >;
  template class SuperLUInverse<double, Vec<12,Complex>, Vec<12,Complex> >;
  template class SuperLUInverse<double, Vec<13,Complex>, Vec<13,Complex> >;
  template class SuperLUInverse<double, Vec<14,Complex>, Vec<14,Complex> >;
  template class SuperLUInverse<double, Vec<15,Complex>, Vec<15,Complex> >;

  template class SuperLUInverse<Complex, Vec<2,Complex>, Vec<2,Complex> >;
  template class SuperLUInverse<Complex, Vec<3,Complex>, Vec<3,Complex> >;
  template class SuperLUInverse<Complex, Vec<4,Complex>, Vec<4,Complex> >;
  template class SuperLUInverse<Complex, Vec<5,Complex>, Vec<5,Complex> >;
  template class SuperLUInverse<Complex, Vec<6,Complex>, Vec<6,Complex> >;
  template class SuperLUInverse<Complex, Vec<7,Complex>, Vec<7,Complex> >;
  template class SuperLUInverse<Complex, Vec<8,Complex>, Vec<8,Complex> >;
  template class SuperLUInverse<Complex, Vec<9,Complex>, Vec<9,Complex> >;
  template class SuperLUInverse<Complex, Vec<10,Complex>, Vec<10,Complex> >;
  template class SuperLUInverse<Complex, Vec<11,Complex>, Vec<11,Complex> >;
  template class SuperLUInverse<Complex, Vec<12,Complex>, Vec<12,Complex> >;
  template class SuperLUInverse<Complex, Vec<13,Complex>, Vec<13,Complex> >;
  template class SuperLUInverse<Complex, Vec<14,Complex>, Vec<14,Complex> >;
  template class SuperLUInverse<Complex, Vec<15,Complex>, Vec<15,Complex> >;
  */


  /*

  template class SuperLUInverse<double>;
  template class SuperLUInverse<Complex>;
#if MAX_SYS_DIM >= 1
  template class SuperLUInverse<Mat<1,1,double> >;
  template class SuperLUInverse<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class SuperLUInverse<Mat<2,2,double> >;
  template class SuperLUInverse<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class SuperLUInverse<Mat<3,3,double> >;
  template class SuperLUInverse<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class SuperLUInverse<Mat<4,4,double> >;
  template class SuperLUInverse<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class SuperLUInverse<Mat<5,5,double> >;
  template class SuperLUInverse<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class SuperLUInverse<Mat<6,6,double> >;
  template class SuperLUInverse<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class SuperLUInverse<Mat<7,7,double> >;
  template class SuperLUInverse<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class SuperLUInverse<Mat<8,8,double> >;
  template class SuperLUInverse<Mat<8,8,Complex> >;
#endif
  */


}









#endif



