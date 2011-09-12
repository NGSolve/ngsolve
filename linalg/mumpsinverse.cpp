/* *************************************************************************/
/* File:   mumpsinverse.cpp                                                */
/* Author: Joachim Schoeberl                                               */
/* Date:   May. 2009                                                       */
/* *************************************************************************/

#ifdef USE_MUMPS

#include <la.hpp>
#include <parallelngs.hpp>


namespace netgen {
  extern int id, ntasks;
}

namespace ngla
{
  using namespace ngstd;
  using namespace netgen;

  
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
    static int timer = NgProfiler::CreateTimer ("Mumps Inverse");
    NgProfiler::RegionTimer reg (timer);

    VT_OFF();

    symmetric = asymmetric;
    inner = ainner;
    cluster = acluster;



    if (id == 0)
      {
	if ( ( inner && inner->Size() < a.Height() ) ||
	     ( cluster && cluster->Size() < a.Height() ) )
	  {
	    cout << "Mumps: Size of inner/cluster does not match matrix size!" << endl;
	    throw Exception("Invalid parameters inner/cluster. Thrown by MumpsInverse.");
	  }
	
	if ( int( mat_traits<TM>::WIDTH) != int(mat_traits<TM>::HEIGHT) )
	  {
	    cout << "Mumps: Each entry in the square matrix has to be a square matrix!" << endl;
	    throw Exception("No Square Matrix. Thrown by MumpsInverse.");
	  }
      }



    clock_t starttime, time1 = 0, time2 = 0; 
    starttime = clock();

    entrysize = mat_traits<TM>::HEIGHT; 
    iscomplex = mat_traits<TM>::IS_COMPLEX;


    int * colstart = 0;
    int * counter = 0;
    int * col_indices = 0, * row_indices = 0;
    TSCAL * matrix = 0;

    if (id == 0)
      {
	height = a.Height() * entrysize;
	
	int * colstart = new int[height+1];
	int * counter = new int[height];
	
	for ( int i = 0; i < height; i++ ) 
	  counter[i] = 0;

	for ( int i = 0; i < height; i++ ) 
	  colstart[i+1] = 0;
	
	if ( symmetric )
	  {
	    cout << "copy matrix symmetric" << endl;
	    
	    col_indices = new int[a.NZE() * entrysize * entrysize ];
	    row_indices = new int[a.NZE() * entrysize * entrysize ];
	    matrix = new TSCAL[a.NZE() * entrysize * entrysize ];      
	    
	    int ii = 0;
	    for (int i = 0; i < a.Height(); i++ )
	      {
		FlatArray<const int> rowind = a.GetRowIndices(i);

		for (int j = 0; j < rowind.Size(); j++ )
		  {
		    int col = rowind[j];
		    
		    if (  (!inner && !cluster) ||
			  (inner && (inner->Test(i) && inner->Test(col) ) ) ||
			  (!inner && cluster &&
			   ((*cluster)[i] == (*cluster)[col] 
			    && (*cluster)[i] ))  )
		      {
			TM entry = a(i,col);
			for (int l = 0; l < entrysize; l++ )
			  for (int k = 0; k < entrysize; k++)
			    {
			      int rowi = i*entrysize+l+1;
			      int coli = col*entrysize+k+1;
			      TSCAL val = Access(entry,l,k);

			      if (rowi >= coli)
				{
				  col_indices[ii] = coli;
				  row_indices[ii] = rowi;
				  matrix[ii] = val;
				  ii++;
				}
			    }
		      }
		    else if (i == col)
		      {
			// in the case of 'inner' or 'cluster': 1 on the diagonal for
			// unused dofs.
			for (int l=0; l<entrysize; l++ )
			  {
			    col_indices[ii] = col*entrysize+l+1;
			    row_indices[ii] = col*entrysize+l+1;
			    matrix[ii] = 1;
			    ii++;
			  }
		      }
		  }
	      }
	    nze = ii;
	  }
	else
	  {
	    cout << "copy matrix non-symmetric" << endl;
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
			      row_indices[ colstart[col*entrysize+k]+
					   counter[col*entrysize+k] ] = i*entrysize+l + 1;
			      col_indices[ colstart[col*entrysize+k]+
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
      }



    mumps_id.job=JOB_INIT; 
    mumps_id.par= (ntasks == 1) ? 1 : 0;
    mumps_id.sym= symmetric ? 1 : 0;
    // mumps_id.comm_fortran=USE_COMM_WORLD;
    mumps_id.comm_fortran = MPI_Comm_c2f (ngparallel::ngs_comm);
    mumps_trait<TSCAL>::MumpsFunction (&mumps_id);

    if (id == 0)
      cout << "MUMPS version number is " << mumps_id.version_number << endl;


    /* Define the problem on the host */
    mumps_id.n   = height; 
    mumps_id.nz  = nze;
    mumps_id.irn = row_indices;
    mumps_id.jcn = col_indices;


    mumps_id.icntl[0]=-1; 
    mumps_id.icntl[1]=-1; 
    mumps_id.icntl[2]=-1; 
    mumps_id.icntl[3]=0;
    // mumps_id.icntl[12]=1;  // root schur complement sequential
    mumps_id.icntl[13]=30; // increased due to error -9
    // cout << "icntl(7) = " << mumps_id.icntl[6] << endl;
    // cout << "icntl(22) = " << mumps_id.icntl[21] << endl;

    // mumps_id.comm_fortran=USE_COMM_WORLD;
    mumps_id.comm_fortran = MPI_Comm_c2f (ngparallel::ngs_comm);
    mumps_id.job = JOB_ANALYSIS;

    if (id == 0)
      cout << "analysis ... " << flush;

    mumps_trait<TSCAL>::MumpsFunction (&mumps_id);

    // cout << "num floating-point ops = " << mumps_id.rinfog[0] << endl;
    if (mumps_id.infog[0])
      {
	cout << "analysis done" << endl;
	cout << "error-code = " << mumps_id.infog[0] << flush;
      }



    mumps_id.a   = (typename mumps_trait<TSCAL>::MUMPS_TSCAL*)matrix; 

    mumps_id.job = JOB_FACTOR;
    
    if (id == 0)
      cout << "factor ... " << flush;

    MPI_Barrier (ngparallel::ngs_comm);


    mumps_trait<TSCAL>::MumpsFunction (&mumps_id);

    if (mumps_id.infog[0] != 0)
      {
	cout << " factorization done" << endl;
	cout << "error-code = " << mumps_id.infog[0] << endl;
	cout << "info(1) = " << mumps_id.info[0] << endl;
	cout << "info(2) = " << mumps_id.info[1] << endl;
      }
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

    
    if (id == 0)
      cout << " done " << endl;
    delete [] colstart;
    delete [] counter;
    // delete [] rhs;    


    delete [] col_indices;
    delete [] row_indices;
    delete [] matrix;
    VT_ON();
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
    static int timer = NgProfiler::CreateTimer ("Mumps mult inverse");
    NgProfiler::RegionTimer reg (timer);

    VT_OFF();

    if (id == 0)
      {
	FlatVector<TVX> fx = x.FV<TVX>();
	FlatVector<TVX> fy = y.FV<TVX>();
	
	fy = fx;
	
	MUMPS_STRUC_C & ncid = const_cast<MUMPS_STRUC_C&> (mumps_id);
	
	ncid.rhs = (typename mumps_trait<TSCAL>::MUMPS_TSCAL*)& (y.FV<TSCAL>()(0));
	
	ncid.job = JOB_SOLVE;
	mumps_trait<TSCAL>::MumpsFunction (&ncid);
	
	if (inner)
	  {
	    for (int i = 0; i < height/entrysize; i++)
	      if (!inner->Test(i)) 
		for (int j = 0; j < entrysize; j++ ) fy(i*entrysize+j) = 0.0;
	  }
	else if (cluster)
	  {
	    for (int i = 0; i < height/entrysize; i++)
	      if (!(*cluster)[i]) 
		for (int j = 0; j < entrysize; j++ ) fy(i*entrysize+j) = 0.0;
	  }
      }
    else
      {
	MUMPS_STRUC_C & ncid = const_cast<MUMPS_STRUC_C&> (mumps_id);
	ncid.rhs = (typename mumps_trait<TSCAL>::MUMPS_TSCAL*)& (y.FV<TSCAL>()(0));
	
	ncid.job = JOB_SOLVE;
	mumps_trait<TSCAL>::MumpsFunction (&ncid);
      }

    VT_ON();
  }
  
  



  template <class TM, class TV_ROW, class TV_COL>
  MumpsInverse<TM,TV_ROW,TV_COL> :: ~MumpsInverse()
  {
    mumps_id.job=JOB_END; 
    /*
    mumps_id.par= (ntasks == 1) ? 1 : 0;
    mumps_id.sym = symmetric ? 1 : 0;
    // mumps_id.comm_fortran=USE_COMM_WORLD;
    mumps_id.comm_fortran = MPI_Comm_c2f (ngparallel::ngs_comm);
    */
    mumps_trait<TSCAL>::MumpsFunction (&mumps_id);
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
  /*
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
  */

  /*
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
  */
}









#endif



