/* *************************************************************************/
/* File:   mumpsinverse.cpp                                                */
/* Author: Joachim Schoeberl                                               */
/* Date:   May. 2009                                                       */
/* *************************************************************************/

#ifdef USE_MUMPS

#include <la.hpp>



namespace ngla
{
  using namespace ngla;
  using namespace ngstd;


  
#define JOB_INIT -1
#define JOB_END -2

#define JOB_ANALYSIS 1
#define JOB_FACTOR 2
#define JOB_SOLVE 3

#define USE_COMM_WORLD -987654

  

  template <class TM, class TV_ROW, class TV_COL>
  MumpsInverse<TM,TV_ROW,TV_COL> :: 
  MumpsInverse (const SparseMatrix<TM,TV_ROW,TV_COL> & a, 
                const BitArray * ainner,
                const Array<int> * acluster,
                bool asymmetric)
  { 
    symmetric = asymmetric;
    inner = ainner;
    cluster = acluster;


    cout << "Mumps called ..." << endl;

    if ( ( inner && inner->Size() < a.Height() ) ||
	 ( cluster && cluster->Size() < a.Height() ) )
      {
	cout << "Mumps: Size of inner/cluster does not match matrix size!" << endl;
	throw Exception("Invalid parameters inner/cluster. Thrown by MumpsInverse.");
      }

    clock_t starttime, time1 = 0, time2 = 0; 
    starttime = clock();

    // prepare matrix and parameters for Mumps
    if ( int( mat_traits<TM>::WIDTH) != int(mat_traits<TM>::HEIGHT) )
      {
	cout << "Mumps: Each entry in the square matrix has to be a square matrix!" << endl;
	throw Exception("No Square Matrix. Thrown by MumpsInverse.");
      }

    entrysize = mat_traits<TM>::HEIGHT;   // NumRows(entry);
    iscomplex = mat_traits<TM>::IS_COMPLEX;

    height = a.Height() * entrysize;
    int * colstart = new int[height+1];
    int * counter = new int[height];

    int * col_indices = 0, * row_indices = 0;
    TSCAL * matrix = 0;

 
    for ( int i = 0; i < height; i++ ) counter[i] = colstart[i+1] = 0;
    
    if ( symmetric )
      {
        /* later 
	// --- transform lower left to full matrix (in compressed column storage format) ---
	// 1.) build array 'colstart':
	// (a) get nr. of entries for each col
	for ( i=1; i<=a.Height(); i++ )
        {
        for ( j=0; j<a.GetRowIndices(i-1).Size(); j++ )
        {
        int col = a.GetRowIndices(i-1)[j];

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
 	col_indices = new int[nze];
 	row_indices = new int[nze];
 	matrix = new TSCAL[nze];      

	for ( i=0; i<a.Height(); i++ )
        {
        for ( j=0; j<a.GetRowIndices(i).Size(); j++ )
        {
        int col = a.GetRowIndices(i)[j];

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
        col_indices[ colstart[col*entrysize+k]+
        counter[col*entrysize+k] ] = i*entrysize+l;
        matrix[ colstart[col*entrysize+k]+
        counter[col*entrysize+k] ] = Elem(entry,l,k);
        counter[col*entrysize+k]++;

        if ( i != col )
        {
        col_indices[ colstart[i*entrysize+l]+
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
        col_indices[ colstart[col*entrysize+l]+
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
        */
	
      }
    else
      {
	// --- transform matrix to compressed column storage format ---

	// 1.) build array 'colstart':
	// (a) get nr. of entries for each col
	for (int i = 0; i < a.Height(); i++ )
	  {
	    for (int j = 0; j < a.GetRowIndices(i).Size(); j++ )
	      {
		int col = a.GetRowIndices(i)[j];
                
		if (  (!inner && !cluster) ||
		      (inner && (inner->Test(i) && inner->Test(col) ) ) ||
		      (!inner && cluster && 
                       ((*cluster)[i] == (*cluster)[col] 
                        && (*cluster)[i] ))  )
		  {
		    for (int k=0; k<entrysize; k++ )
                      colstart[col*entrysize+k+1] += entrysize;
		  }
		else if ( i == col )
		  {
		    for (int k=0; k<entrysize; k++ )
		      colstart[col*entrysize+k+1] ++;
		  }
	      }
	  }

	// (b) accumulate
	colstart[0] = 0;
	for (int i = 1; i <= height; i++ ) colstart[i] += colstart[i-1];
	nze = colstart[height];


	// 2.) build whole matrix:
	col_indices = new int[a.NZE() * entrysize * entrysize ];
	row_indices = new int[a.NZE() * entrysize * entrysize ];
	matrix = new TSCAL[a.NZE() * entrysize * entrysize ];      

	for (int i = 0; i < a.Height(); i++ )
	  {
	    for (int j = 0; j<a.GetRowIndices(i).Size(); j++ )
	      {
		int col = a.GetRowIndices(i)[j];

		if (  (!inner && !cluster) ||
		      (inner && (inner->Test(i) && inner->Test(col) ) ) ||
		      (!inner && cluster &&
                       ((*cluster)[i] == (*cluster)[col] 
                        && (*cluster)[i] ))  )
		  {
		    TM entry = a(i,col);
		    for (int k = 0; k < entrysize; k++)
		      for (int l = 0; l < entrysize; l++ )
			{
			  col_indices[ colstart[col*entrysize+k]+
                                       counter[col*entrysize+k] ] = i*entrysize+l + 1;
			  row_indices[ colstart[col*entrysize+k]+
                                       counter[col*entrysize+k] ] = col*entrysize+k + 1;
			  matrix[ colstart[col*entrysize+k]+
				  counter[col*entrysize+k] ] = Access(entry,l,k);
			  counter[col*entrysize+k]++;
			}
		  }
		else if (i == col)
		  {
		    // in the case of 'inner' or 'cluster': 1 on the diagonal for
		    // unused dofs.
		    for (int l=0; l<entrysize; l++ )
		      {
			col_indices[ colstart[col*entrysize+l]+
                                     counter[col*entrysize+l] ] = col*entrysize+l + 1;
			row_indices[ colstart[col*entrysize+l]+
                                     counter[col*entrysize+l] ] = col*entrysize+l + 1;
			matrix[ colstart[col*entrysize+l]+
				counter[col*entrysize+l] ] = 1;
			counter[col*entrysize+l]++;
		      }
		  }
	      }
	  }
      }


    id.job=JOB_INIT; 
    id.par=1; 
    id.sym=0;
    id.comm_fortran=USE_COMM_WORLD;
    // dmumps_c(&id);
    mumps_trait<TSCAL>::MumpsFunction (&id);

    cout << "MUMPS version number is " << id.version_number << endl;


    /* Define the problem on the host */
    id.n   = height; 
    id.nz  = nze;
    id.irn = col_indices;
    id.jcn = row_indices;
    // id.rhs = rhs;

    if (symmetric)
      id.sym = 1;   // spd
    else
      id.sym = 0;   // non-symmetric
    

    id.icntl[0]=-1; 
    id.icntl[1]=-1; 
    id.icntl[2]=-1; 
    id.icntl[3]=0;

    id.job = JOB_ANALYSIS;
    // dmumps_c(&id);
    mumps_trait<TSCAL>::MumpsFunction (&id);

    cout << "analysis done, now factor" << endl;

    // id.a   = (double*)(void*)matrix; 
    id.a   = (typename mumps_trait<TSCAL>::MUMPS_TSCAL*)matrix; 
    // id.a   = (typename mumps_trait<TSCAL>::MUMPS_TSCAL*)(void*)matrix; 


    id.job = JOB_FACTOR;
    // dmumps_c(&id);
    mumps_trait<TSCAL>::MumpsFunction (&id);

    cout << "factorization done" << endl;

    time2 = clock();




    /*
    if ( error != 0 )
      {
        cout << "Setup and Factorization: Mumps returned error " << error << "!" << endl;
        throw Exception("MumpsInverse: Setup and Factorization failed.");
      }
    */


    (*testout) << endl << "Direct Solver: Mumps by Lawrence Berkeley National Laboratory." << endl;
    (*testout) << "Matrix prepared for Mumps in " <<
      double(time1 - starttime)/CLOCKS_PER_SEC << " sec." << endl;
    (*testout) << "Factorization by Mumps done in " << 
      double(time2 - time1)/CLOCKS_PER_SEC << " sec." << endl << endl;

    
    cout << " done " << endl;
    delete [] colstart;
    delete [] counter;
    // delete [] rhs;    
  }
  
  
  /*
    template <class TM, class TV_ROW, class TV_COL>
    MumpsInverse<TM,TV_ROW,TV_COL> :: 
    MumpsInverse (const Array<int> & aorder, 
    const Array<CliqueEl*> & cliques,
    const Array<MDOVertex> & vertices,
    int symmetric)
    {
    Allocate (aorder, cliques, vertices);
    }
  */

  
  /*
    template <class TM, class TV_ROW, class TV_COL>
    void MumpsInverse<TM, TV_ROW,TV_COL> :: 
    Allocate (const Array<int> & aorder, 
    const Array<CliqueEl*> & cliques,
    const Array<MDOVertex> & vertices)
    {
    cout << "MumpsInverse::Allocate not implemented!" << endl;
    }
  


    template <class TM, class TV_ROW, class TV_COL>
    void MumpsInverse<TM,TV_ROW,TV_COL> :: 
    FactorNew (const SparseMatrix<TM> & a)
    {
    throw Exception ("MumpsInverse::FactorNew not implemented");
    }



    template <class TM, class TV_ROW, class TV_COL>
    void MumpsInverse<TM,TV_ROW, TV_COL> :: Factor (const int * blocknr)
    {
    cout << "MumpsInverse::Factor not implemented!" << endl;
    }
  */  




  template <class TM, class TV_ROW, class TV_COL>
  void MumpsInverse<TM,TV_ROW,TV_COL> :: 
  Mult (const BaseVector & x, BaseVector & y) const
  {

    FlatVector<TVX> fx = 
      dynamic_cast<T_BaseVector<TVX> &> (const_cast<BaseVector &> (x)).FV();
    FlatVector<TVX> fy = 
      dynamic_cast<T_BaseVector<TVX> &> (y).FV();
    
    fy = fx;

    MUMPS_STRUC_C & ncid = const_cast<MUMPS_STRUC_C&> (id);

    // ncid.rhs = &y.FVDouble()(0);
    ncid.rhs = (typename mumps_trait<TSCAL>::MUMPS_TSCAL*)&fy(0);

    ncid.job = JOB_SOLVE;
    mumps_trait<TSCAL>::MumpsFunction (&ncid);


    if (inner)
      {
	for (int i = 0; i < height/entrysize; i++)
	  if (!inner->Test(i)) 
	    for (int j = 0; j < entrysize; j++ ) fy(i*entrysize+j) = 0;
      }
    else if (cluster)
      {
	for (int i = 0; i < height/entrysize; i++)
	  if (!(*cluster)[i]) 
	    for (int j = 0; j < entrysize; j++ ) fy(i*entrysize+j) = 0;
      }
  }
  
  




  /*
    template <class TM, class TV_ROW, class TV_COL>
    void MumpsInverse<TM,TV_ROW,TV_COL> :: Set (int i, int j, const TM & val)
    {
    cout << "MumpsInverse::Set not implemented!" << endl;
    }



    template <class TM, class TV_ROW, class TV_COL>
    const TM & MumpsInverse<TM,TV_ROW,TV_COL> :: Get (int i, int j) const
    {
    cout << "MumpsInverse::Get not implemented!" << endl;
    }


    template <class TM, class TV_ROW, class TV_COL>
    ostream & MumpsInverse<TM,TV_ROW,TV_COL> :: Print (ostream & ost) const
    {
    cout << "MumpsInverse::Print not implemented!" << endl;
    return ost; 
    }
  */


  template <class TM, class TV_ROW, class TV_COL>
  MumpsInverse<TM,TV_ROW,TV_COL> :: ~MumpsInverse()
  {
    /*
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
    */
  }







  template class MumpsInverse<double>;
  template class MumpsInverse<Complex>;
  template class MumpsInverse<double,Complex,Complex>;
#if MAX_SYS_DIM >= 1
  template class MumpsInverse<Mat<1,1,double> >;
  template class MumpsInverse<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class MumpsInverse<Mat<2,2,double> >;
  template class MumpsInverse<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class MumpsInverse<Mat<3,3,double> >;
  template class MumpsInverse<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class MumpsInverse<Mat<4,4,double> >;
  template class MumpsInverse<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class MumpsInverse<Mat<5,5,double> >;
  template class MumpsInverse<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class MumpsInverse<Mat<6,6,double> >;
  template class MumpsInverse<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class MumpsInverse<Mat<7,7,double> >;
  template class MumpsInverse<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class MumpsInverse<Mat<8,8,double> >;
  template class MumpsInverse<Mat<8,8,Complex> >;
#endif



  template class MumpsInverse<double, Vec<2,double>, Vec<2,double> >;
  template class MumpsInverse<double, Vec<3,double>, Vec<3,double> >;
  template class MumpsInverse<double, Vec<4,double>, Vec<4,double> >;
  template class MumpsInverse<double, Vec<5,double>, Vec<5,double> >;
  template class MumpsInverse<double, Vec<6,double>, Vec<6,double> >;
  template class MumpsInverse<double, Vec<7,double>, Vec<7,double> >;
  template class MumpsInverse<double, Vec<8,double>, Vec<8,double> >;
  template class MumpsInverse<double, Vec<9,double>, Vec<9,double> >;
  template class MumpsInverse<double, Vec<10,double>, Vec<10,double> >;
  template class MumpsInverse<double, Vec<11,double>, Vec<11,double> >;
  template class MumpsInverse<double, Vec<12,double>, Vec<12,double> >;
  template class MumpsInverse<double, Vec<13,double>, Vec<13,double> >;
  template class MumpsInverse<double, Vec<14,double>, Vec<14,double> >;
  template class MumpsInverse<double, Vec<15,double>, Vec<15,double> >;
  
  template class MumpsInverse<double, Vec<2,Complex>, Vec<2,Complex> >;
  template class MumpsInverse<double, Vec<3,Complex>, Vec<3,Complex> >;
  template class MumpsInverse<double, Vec<4,Complex>, Vec<4,Complex> >;
  template class MumpsInverse<double, Vec<5,Complex>, Vec<5,Complex> >;
  template class MumpsInverse<double, Vec<6,Complex>, Vec<6,Complex> >;
  template class MumpsInverse<double, Vec<7,Complex>, Vec<7,Complex> >;
  template class MumpsInverse<double, Vec<8,Complex>, Vec<8,Complex> >;
  template class MumpsInverse<double, Vec<9,Complex>, Vec<9,Complex> >;
  template class MumpsInverse<double, Vec<10,Complex>, Vec<10,Complex> >;
  template class MumpsInverse<double, Vec<11,Complex>, Vec<11,Complex> >;
  template class MumpsInverse<double, Vec<12,Complex>, Vec<12,Complex> >;
  template class MumpsInverse<double, Vec<13,Complex>, Vec<13,Complex> >;
  template class MumpsInverse<double, Vec<14,Complex>, Vec<14,Complex> >;
  template class MumpsInverse<double, Vec<15,Complex>, Vec<15,Complex> >;

  template class MumpsInverse<Complex, Vec<2,Complex>, Vec<2,Complex> >;
  template class MumpsInverse<Complex, Vec<3,Complex>, Vec<3,Complex> >;
  template class MumpsInverse<Complex, Vec<4,Complex>, Vec<4,Complex> >;
  template class MumpsInverse<Complex, Vec<5,Complex>, Vec<5,Complex> >;
  template class MumpsInverse<Complex, Vec<6,Complex>, Vec<6,Complex> >;
  template class MumpsInverse<Complex, Vec<7,Complex>, Vec<7,Complex> >;
  template class MumpsInverse<Complex, Vec<8,Complex>, Vec<8,Complex> >;
  template class MumpsInverse<Complex, Vec<9,Complex>, Vec<9,Complex> >;
  template class MumpsInverse<Complex, Vec<10,Complex>, Vec<10,Complex> >;
  template class MumpsInverse<Complex, Vec<11,Complex>, Vec<11,Complex> >;
  template class MumpsInverse<Complex, Vec<12,Complex>, Vec<12,Complex> >;
  template class MumpsInverse<Complex, Vec<13,Complex>, Vec<13,Complex> >;
  template class MumpsInverse<Complex, Vec<14,Complex>, Vec<14,Complex> >;
  template class MumpsInverse<Complex, Vec<15,Complex>, Vec<15,Complex> >;
}









#endif



