/* *************************************************************************/
/* File:   pardisoinverse.cpp                                              */
/* Author: Florian Bachinger                                               */
/* Date:   Feb. 2004                                                       */
/* *************************************************************************/

#include <la.hpp>

#ifdef USE_PARDISO

#define F77_FUNC(func)  func ## _

extern "C"
{
  /* PARDISO prototype. */

#ifdef USE_MKL

  void F77_FUNC(pardiso)
    (void * pt, int * maxfct, int * mnum, int * mtype, int * phase, int * n, 
     double * a, int * ia, int * ja, int * perm, int * nrhs, int * iparam, 
     int * msglvl, double * b, double * x, int * error);

#else

#ifdef USE_PARDISO400
extern  int F77_FUNC(pardisoinit)
    (void *, int *, int *, int *, double *, int *);
#else
extern  int F77_FUNC(pardisoinit)
    (void *, int *, int *);
#endif
  int F77_FUNC(pardiso)
    (void * pt, int * maxfct, int * mnum, int * mtype, int * phase, int * n, 
     double * a, int * ia, int * ja, int * perm, int * nrhs, int * iparam, 
     int * msglvl, double * b, double * x, int * error);

#endif
}






namespace ngla
{
  using namespace ngla;
  using namespace ngstd;


  int pardiso_msg = 0;


  template<class TM, class TV_ROW, class TV_COL>
  PardisoInverse<TM,TV_ROW,TV_COL> :: 
  PardisoInverse (const SparseMatrix<TM,TV_ROW,TV_COL> & a, 
		  const BitArray * ainner,
		  const Array<int> * acluster,
		  int asymmetric)
    : SparseFactorization (a)
  { 
    static int timer = NgProfiler::CreateTimer ("Pardiso Inverse");
    NgProfiler::RegionTimer reg (timer);


    if (getenv ("PARDISOMSG"))
      pardiso_msg = 1;



    print = bool (pardiso_msg); // false;

    symmetric = asymmetric;
    inner = ainner;
    cluster = acluster;

    (*testout) << "Pardiso, symmetric = " << symmetric << endl;

    if (inner && cluster)
      throw Exception("PardisoInverse: Cannot use inner and cluster");

    if ( (inner && inner->Size() < a.Height()) ||
	 (cluster && cluster->Size() < a.Height() ) )
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
    height = a.Height() * entrysize;

    rowstart = new int[height+1];
    indices = new int[a.NZE() * entrysize * entrysize ];
    matrix = new TSCAL[a.NZE() * entrysize * entrysize ];     

    *testout << "matrix.InverseTpye = " <<  a.GetInverseType() << endl;
    spd = ( a.GetInverseType() == PARDISOSPD ) ? 1 : 0;

    int maxfct = 1, mnum = 1, phase = 12, nrhs = 1, msglevel = print, error;
    int * params = const_cast <int*> (&hparams[0]);

    for (int i = 0; i < 64; i++)
      params[i] = 0;

    params[0] = 1; // no pardiso defaults
    params[2] = 8; // 1 processor (?)

    params[1] = 0; // fill in 0..MDO, 2..metis
    params[3] = params[4] = params[5] = params[7] = params[8] = 
      params[11] = params[12] = params[18] = 0;
    params[6] = 16;
    params[9] = 10;  // perturbation 1E-10
    params[10] = 1;

    // JS
    params[6] = 0;
    params[17] = -1;
    params[20] = 1;  // 1x1 and 2x2 bunc and Kaufman pivoting
    
    params[26] = 1; // check input matrix
    params[59] = 0; // 0..incore, 1..depending on MKL_PARDISO_MAX_CORE_SIZE, 2..ooc

    for (int i = 0; i < 128; i++) pt[i] = 0;

#ifdef USE_MKL
    //    no init in MKL PARDISO
    //     retvalue = F77_FUNC(pardiso) ( pt, &maxfct, &mnum, &matrixtype, &phase, &height, 
    // 				   reinterpret_cast<double *>(matrix),
    // 				   rowstart, indices, NULL, &nrhs, params, &msglevel,
    // 				   NULL, NULL, &error );
#else
    int retvalue;
#ifdef USE_PARDISO400
    double dparm[64];
    int solver = 0;
    F77_FUNC(pardisoinit) (pt,  &matrixtype, &solver, params, dparm, &retvalue); 
#else
    retvalue = 
      F77_FUNC(pardisoinit) (pt,  &matrixtype, params); 
#endif
    // cout << "init success" << endl;
    // cout << "retvalue = " << retvalue << endl;
#endif
    
    SetMatrixType();

    if ( symmetric )
      {
	// --- transform lower left to upper right triangular matrix ---
	// 1.) build array 'rowstart':
	// (a) get nr. of entries for each row

	for (int i=0; i<=height; i++ ) rowstart[i] = 0;

	for (int i=0; i < a.Height(); i++ )
	  {
	    for (int j = 0; j < a.GetRowIndices(i).Size(); j++ )
	      {
		int col = a.GetRowIndices(i)[j];
		if ( i != col )
		  {
		    if (  (!inner && !cluster) ||
			  (inner && (inner->Test(i) && inner->Test(col) ) ) ||
			  (!inner && cluster && 
		           ((*cluster)[i] == (*cluster)[col] 
			    && (*cluster)[i] ))  )
		      {
			for (int k=0; k<entrysize; k++ )
			  rowstart[col*entrysize+k+1] += entrysize;
		      }
		  }
		else if ( (!inner && !cluster) || 
			  (inner && inner->Test(i)) ||
			  (!inner && cluster && (*cluster)[i]) )
		  {
		    for (int k=0; k<entrysize; k++ )
		      rowstart[col*entrysize+k+1] += entrysize-k;
		  }
		else
		  {
		    for (int k=0; k<entrysize; k++ )
		      rowstart[col*entrysize+k+1] ++;
		  }
	      }
	  }
	// (b) accumulate
	rowstart[0] = 1;
	for (int i=1; i<=height; i++ ) rowstart[i] += rowstart[i-1];
	
	
	// 2.) build whole matrix:
	int * counter = new int[height];
	for (int i=0; i<a.Height(); i++ )
	  {
	    // int rowsize = a.GetRowIndices(i).Size();
	    for (int k=0; k<entrysize; k++ ) counter[i*entrysize+k]=0;

	    for (int j=0; j<a.GetRowIndices(i).Size(); j++ )
	      {
		int col = a.GetRowIndices(i)[j];

		if ( i != col )
		  {
		    if ( (!inner && !cluster) ||
			 (inner && (inner->Test(i) && inner->Test(col)) ) ||
			 (!inner && cluster && 
			  ((*cluster)[i] == (*cluster)[col] 
			   && (*cluster)[i] ))  )
		      {
			TM entry = a(i,col);
			for (int k=0; k<entrysize; k++)
			  for (int l=0; l<entrysize; l++ )
			    {
			      indices[ rowstart[col*entrysize+k]+
				       counter[col*entrysize+k]-1 ] = i*entrysize+l+1;
			      matrix[ rowstart[col*entrysize+k]+
				      counter[col*entrysize+k]-1 ] = Access(entry,l,k);
			      counter[col*entrysize+k]++;
			    }
		      }
		  }
		else if ( (!inner && !cluster) || 
			  (inner && inner->Test(i)) ||
			  (!inner && cluster && (*cluster)[i]) )
		  {
		    TM entry = a(i,col);
		    for (int l=0; l<entrysize; l++ )
		      for (int k=0; k<=l; k++)
			{
			  indices[ rowstart[col*entrysize+k]+
				   counter[col*entrysize+k]-1 ] = i*entrysize+l+1;
			  matrix[ rowstart[col*entrysize+k]+
				  counter[col*entrysize+k]-1 ] = Access(entry,l,k);
			  counter[col*entrysize+k]++;
			}
		  }
		else
		  {
		    // in the case of 'inner' or 'cluster': 1 on the diagonal for
		    // unused dofs.
		    for (int l=0; l<entrysize; l++ )
		      {
			indices[ rowstart[col*entrysize+l]+
				 counter[col*entrysize+l]-1 ] = i*entrysize+l+1;
			matrix[ rowstart[col*entrysize+l]+
				counter[col*entrysize+l]-1 ] = 1;
			counter[col*entrysize+l]++;
		      }
		  }
	      }
	  }
	

	/*
	  (*testout) << endl << "row, rowstart / indices, matrix-entries" << endl;
	  for ( i=0; i<height; i++ )
	  {
	  (*testout) << endl << i+1 << ", " << rowstart[i] << ":   ";
	  for ( j=rowstart[i]-1; j<rowstart[i+1]-1; j++ )
	  (*testout) << indices[j] << ", " << matrix[j] << "      ";
	  }
	  (*testout) << endl;
	*/

	delete [] counter;
      }
    else
      {
	if (print)
	  cout << "non-symmetric pardiso, cluster = " << cluster << ", entrysize = " << entrysize << endl;

	// for non-symmetric matrices:
	int counter = 0, rowelems;

	for (int i = 0; i < a.Height(); i++ )
	  {
	    rowelems = 0;
	    int rowsize = a.GetRowIndices(i).Size();

	    // get effective number of elements in the row

	    if ( inner )
	      {
		for (int j = 0; j < rowsize; j++)
		  {
		    int col = a.GetRowIndices(i)[j];
		    if (inner->Test(i) && inner->Test(col))
		      rowelems += entrysize;
		    else if ( i == col )
		      rowelems++;
		  }
	      }
            else if (cluster)
	      {
		for (int j = 0; j < rowsize; j++ )
		  {
		    int col = a.GetRowIndices(i)[j];
		    if ( (*cluster)[i] == (*cluster)[col] && (*cluster)[i])
		      rowelems+=entrysize;
		    else if ( i == col )
		      rowelems++;
		  }
	      }
	    else 
	      rowelems = rowsize * entrysize;


	    for (int k = 0; k < entrysize; k++ )
	      rowstart[i*entrysize+k] = counter+rowelems*k+1;
	    
	    //	    (*testout) << "rowelems: " << rowelems << endl;


            /*
	    for ( j=0; j<rowsize; j++ )
	      {
		col = a.GetRowIndices(i)[j];

		if ( (!inner && !cluster) ||
		     (inner && ( inner->Test(i) && inner->Test(col) )) ||
		     (!inner && cluster &&
		      ((*cluster)[i] == (*cluster)[col] 
		       && (*cluster)[i] )) )
		  {
		    entry = a(i,col);
		    for ( k=0; k<entrysize; k++ )
		      for ( l=0; l<entrysize; l++ )
			{
			  indices[counter+rowelems*k+l] = col*entrysize+l+1;
			  matrix[counter+rowelems*k+l] = Access(entry,k,l);
			}

		    counter+=entrysize;
		  }
		else if ( i == col )
		  {
		    // in the case of 'inner' or 'cluster': 1 on the diagonal for
		    // unused dofs.
		    for ( l=0; l<entrysize; l++ )
		      {
			indices[counter+rowelems*l+l] = col*entrysize+l+1;
			matrix[counter+rowelems*l+l] = 1;
		      }
		    counter++;
		  }
	      }
            */

            if (inner)
              {
                for (int j=0; j<rowsize; j++ )
                  {
                    int col = a.GetRowIndices(i)[j];
                    
                    if ( (!inner && !cluster) ||
                         (inner && ( inner->Test(i) && inner->Test(col) )) ||
                         (!inner && cluster &&
		      ((*cluster)[i] == (*cluster)[col] 
		       && (*cluster)[i] )) )
                      {
                        TM entry = a(i,col);
                        for (int k=0; k<entrysize; k++ )
                          for (int l=0; l<entrysize; l++ )
                            {
                              indices[counter+rowelems*k+l] = col*entrysize+l+1;
                              matrix[counter+rowelems*k+l] = Access(entry,k,l);
                            }
                        
                        counter+=entrysize;
                      }
                    else if ( i == col )
                      {
		    // in the case of 'inner' or 'cluster': 1 on the diagonal for
		    // unused dofs.
                        for (int l=0; l<entrysize; l++ )
                          {
                            indices[counter+rowelems*l+l] = col*entrysize+l+1;
                            matrix[counter+rowelems*l+l] = 1;
                          }
                        counter++;
                      }
                  }
              }
            else if (cluster)
              {
                for (int j=0; j<rowsize; j++ )
                  {
                    int col = a.GetRowIndices(i)[j];
                    
                    if ( (*cluster)[i] == (*cluster)[col] && (*cluster)[i] )
                      {
                        TM entry = a(i,col);
                        for (int k=0; k<entrysize; k++ )
                          for (int l=0; l<entrysize; l++ )
                            {
                              indices[counter+rowelems*k+l] = col*entrysize+l+1;
                              matrix[counter+rowelems*k+l] = Access(entry,k,l);
                            }
                        counter+=entrysize;
                      }
                    else if ( i == col )
                      {
                        for (int l=0; l<entrysize; l++ )
                          {
                            indices[counter+rowelems*l+l] = col*entrysize+l+1;
                            matrix[counter+rowelems*l+l] = 1;
                          }
                        counter+=entrysize;
                      }
                  }
              }
            else
              for (int j = 0; j < rowsize; j++)
                {
                  int col = a.GetRowIndices(i)[j];
                  // TM entry = a(i,col);
		  const TM & entry = a.GetRowValue(i,j);
                  for (int k = 0; k < entrysize; k++)
                    for (int l = 0; l < entrysize; l++)
                      {
                        indices[counter+rowelems*k+l] = col*entrysize+l+1;
                        matrix[counter+rowelems*k+l] = Access(entry,k,l);
                      }
                  
                  counter+=entrysize;
                }


	    // counter += rowelems * ( entrysize-1 );    should not be here ???? (JS, Oct 2009)
	  }
	rowstart[height] = counter+1;
	
	/*
	  (*testout) << endl << "row, rowstart / indices, matrix-entries" << endl;
	  for ( i=0; i<height; i++ )
	  {
	  (*testout) << endl << i+1 << ", " << rowstart[i] << ":   ";
	  for ( j=rowstart[i]-1; j<rowstart[i+1]-1; j++ )
	  (*testout) << indices[j] << ", " << matrix[j] << "      ";
	  }
	  (*testout) << endl;
	*/

      }

    nze = rowstart[height];

    
    {
      int n = a.Height();
      used.SetSize (n);
      used.Set();
      
      for (int i = 0; i < n; i++)
	if (a.GetPositionTest (i,i) == -1)
	  used.Clear(i);
      
      if (inner)
	for (int i = 0; i < n; i++)
	  if (!inner->Test(i))
	    used.Clear(i);
      
      if (cluster)
	for (int i = 0; i < n; i++)
	  if (!(*cluster)[i])
	    used.Clear(i);
    }

    // call pardiso for factorization:
    // time1 = clock();
    cout << "call pardiso ..." << flush;


    


    // retvalue = 
    F77_FUNC(pardiso) ( pt, &maxfct, &mnum, &matrixtype, &phase, &height, 
			reinterpret_cast<double *>(matrix),
			rowstart, indices, NULL, &nrhs, params, &msglevel,
			NULL, NULL, &error );
    
    cout << " done" << endl;

    if ( error != 0 )
      {
	cout << "Setup and Factorization: PARDISO returned error " << error << "!" << endl;
	
	string errmsg;
	switch (error)
	  {
	  case -1: errmsg = "input inconsistent"; break;
	  case -2: errmsg = "not enough memory"; break;
	  case -3: errmsg = "reordering problem"; break;
	  case -4: errmsg = "zero pivot, numerical factorization or iterative refinement problem"; break;
	  case -5: errmsg = "unclassified (internal) error"; break;
	  case -6: errmsg = "preordering failed"; break;
	  default: 
	    ;
	  }
	
	cout << "err = " << errmsg << endl;

	switch (error)
	  {
	  case -4: 
	    {
	      cout << "iparam(20) = " << params[19] << endl; break;
	    }
	  default:
	    ;
	  }

	ofstream err("pardiso.err");
	err << "ngsolve-matrix = " << endl << a << endl;
	err << "pardiso matrix = " << endl;
	for (int i = 0; i < height; i++)
	  {
	    err << "Row " << i << " start " << rowstart[i] << ": ";
	    if ( inner ) err << " free=" << inner->Test(i) << " ";
	    if ( cluster ) err << " cluster=" << (*cluster)[i] << " ";
	    for (int j = rowstart[i]; j < rowstart[i+1]; j++)
	      err << "c=" << indices[j-1]-1 << ", v=" << matrix[j-1] << "   ";
	    err << "\n";
	  }
	
	cout << "wrote matrix to file 'pardiso.err', please check" << endl;
	throw Exception("PardisoInverse: Setup and Factorization failed.");
      }

    /*
    (*testout) << endl << "Direct Solver: PARDISO by Schenk/Gaertner." << endl;
    (*testout) << "Matrix prepared for PARDISO in " <<
      double(time1 - starttime)/CLOCKS_PER_SEC << " sec." << endl;
    (*testout) << "Factorization by PARDISO done in " << 
      double(time2 - time1)/CLOCKS_PER_SEC << " sec." << endl << endl;
    */
  }
  
  /*
  template<class TM, class TV_ROW, class TV_COL>
  PardisoInverse<TM,TV_ROW,TV_COL> :: 
  PardisoInverse (const Array<int> & aorder, 
		  const Array<CliqueEl*> & cliques,
		  const Array<MDOVertex> & vertices,
		  int symmetric)
  {
    Allocate (aorder, cliques, vertices);
  }

  
  template<class TM, class TV_ROW, class TV_COL>
  void PardisoInverse<TM,TV_ROW,TV_COL> :: 
  Allocate (const Array<int> & aorder, 
	    const Array<CliqueEl*> & cliques,
	    const Array<MDOVertex> & vertices)
  {
    cout << "PardisoInverse::Allocate not implemented!" << endl;
  }
  */
  

  template<class TM, class TV_ROW, class TV_COL>
  void PardisoInverse<TM,TV_ROW,TV_COL> :: 
  FactorNew (const SparseMatrix<TM,TV_ROW,TV_COL> & a)
  {
    throw Exception ("PardisoInverse::FactorNew not implemented");
  }





  template<class TM, class TV_ROW, class TV_COL>
  void PardisoInverse<TM,TV_ROW,TV_COL> :: Factor (const int * blocknr)
  {
    cout << "PardisoInverse::Factor not implemented!" << endl;
  }
  
  


  template<class TM, class TV_ROW, class TV_COL>
  void PardisoInverse<TM,TV_ROW,TV_COL> :: 
  Mult (const BaseVector & x, BaseVector & y) const
  {
    static int timer = NgProfiler::CreateTimer ("Pardiso Solve");
    NgProfiler::RegionTimer reg (timer);

    FlatVector<TVX> fx = x.FV<TVX> ();
    // dynamic_cast<T_BaseVector<TVX> &> (const_cast<BaseVector &> (x)).FV();
    FlatVector<TVX> fy = y.FV<TVX> ();
    // dynamic_cast<T_BaseVector<TVX> &> (y).FV();
    

    int maxfct = 1, mnum = 1, phase = 33, nrhs = 1, msglevel = 0, error;
    // int params[64];
    int * params = const_cast <int*> (&hparams[0]);

    /*
      params[0] = 1; // no pardiso defaults
      params[2] = 1; // 1 processor

      params[1] = 2;
      params[3] = params[4] = params[5] = params[7] = params[8] = 
      params[11] = params[12] = params[18] = 0;
      params[6] = 16;
      params[9] = 20;
      params[10] = 1;
    */

    /*
    params[0] = 1; // no pardiso defaults
    params[1] = 2;
    params[2] = 1; // 1 processor

    params[3] = params[4] = params[5] = params[7] = params[8] = 
      params[11] = params[12] = params[18] = 0;
    params[6] = 16;
    params[9] = 20;
    params[10] = 1;
    */

    // JS
    // params[6] = 0;
    // params[17] = -1;
    // params[20] = 0;


    /*
    params[0] = 0; // 
    params[2] = 8; 
    */

    


#ifdef USE_PARDISO400
    F77_FUNC(pardiso) ( const_cast<long int *>(pt), &maxfct, &mnum, const_cast<int *>(&matrixtype),
			&phase, const_cast<int *>(&height), 
			reinterpret_cast<double *>(matrix),
			rowstart, indices,
			NULL, &nrhs, params, &msglevel,
			static_cast<double *>(fx.Data()), 
			static_cast<double *>(fy.Data()), &error );
#else
    F77_FUNC(pardiso) ( const_cast<int *>(pt), &maxfct, &mnum, const_cast<int *>(&matrixtype),
			&phase, const_cast<int *>(&height), 
			reinterpret_cast<double *>(matrix),
			rowstart, indices,
			NULL, &nrhs, params, &msglevel,
			static_cast<double *>(fx.Data()), 
			static_cast<double *>(fy.Data()), &error );
#endif
    if ( error != 0 )
      cout << "Apply Inverse: PARDISO returned error " << error << "!" << endl;

    for (int i=0; i<height/entrysize; i++)
      if (!used.Test(i))
	for (int j=0; j<entrysize; j++ ) fy(i*entrysize+j) = 0.0;
    

    /*
    if (inner)
      {
	for (int i=0; i<height/entrysize; i++)
	  if (!inner->Test(i)) 
	    for (int j=0; j<entrysize; j++ ) fy(i*entrysize+j) = 0;
      }
    else if (cluster)
      {
	for (int i=0; i<height/entrysize; i++)
	  if (!(*cluster)[i]) 
	    for (int j=0; j<entrysize; j++ ) fy(i*entrysize+j) = 0;
      }
    */

  }
  
  






  /*
  template<class TM, class TV_ROW, class TV_COL>
  void PardisoInverse<TM,TV_ROW,TV_COL> :: Set (int i, int j, const TM & val)
  {
    cout << "PardisoInverse::Set not implemented!" << endl;
  }


  template<class TM, class TV_ROW, class TV_COL>
  const TM & PardisoInverse<TM,TV_ROW,TV_COL> :: Get (int i, int j) const
  {
    cout << "PardisoInverse::Get not implemented!" << endl;
    static TM dummy(0.0);
    return dummy;
  }
  */


  template<class TM, class TV_ROW, class TV_COL>
  ostream & PardisoInverse<TM,TV_ROW,TV_COL> :: Print (ostream & ost) const
  {
    cout << "PardisoInverse::Print not implemented!" << endl;
    return ost; 
  }


  template<class TM, class TV_ROW, class TV_COL>
  PardisoInverse<TM,TV_ROW,TV_COL> :: ~PardisoInverse()
  {
    int maxfct = 1, mnum = 1, phase = -1, nrhs = 1, msglevel = 1, error;
    int * params = const_cast <int*> (&hparams[0]);

    //    cout << "call pardiso (clean up) ..." << endl;
    F77_FUNC(pardiso) ( pt, &maxfct, &mnum, &matrixtype, &phase, &height, NULL,
			rowstart, indices, NULL, &nrhs, params, &msglevel,
			NULL, NULL, &error );
    if ( error != 0 )
      cout << "Clean Up: PARDISO returned error " << error << "!" << endl;

    delete [] rowstart;
    delete [] indices;
    delete [] matrix;
  }






  template<class TM, class TV_ROW, class TV_COL>
  void PardisoInverse<TM,TV_ROW,TV_COL> :: SetMatrixType() // TM entry)
  {
    // if ( IsComplex(entry) )
    if (mat_traits<TM>::IS_COMPLEX)
      {
	if ( symmetric ) matrixtype = 6;
	else
	  {
	    matrixtype = 13;
	    (*testout) << "PARDISO: Assume matrix type to be complex non-symmetric." << endl;
	    (*testout) << "Warning: Works, but is not optimal for Hermitian matrices." << endl;
	  }
      }
    else
      {
	if ( symmetric )
	  {
	    /*
	    matrixtype = 2;    // +2 .. pos def., -2  .. symmetric indef
	    (*testout) << "PARDISO: Assume matrix type to be symmetric positive definite." << endl;
	    (*testout) << "Warning: This might cause errors for symmetric indefinite matrices!!" << endl;
	    */

	    if ( spd ) 
	      matrixtype = 2;
	    else
	      matrixtype = -2;
	  }
	else matrixtype = 11;
      }

    if (print)
      cout << "spd = " << int(spd) << ", sym = " << int(symmetric) 
	   << ", complex = " << int(mat_traits<TM>::IS_COMPLEX)
	   << ", matrixtype = " << matrixtype << endl;
    *testout << "pardiso matrixtype = " << matrixtype << endl;
  }







  template class PardisoInverse<double>;
  template class PardisoInverse<Complex>;
  template class PardisoInverse<double,Complex,Complex>;
#if MAX_SYS_DIM >= 1
  template class PardisoInverse<Mat<1,1,double> >;
  template class PardisoInverse<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class PardisoInverse<Mat<2,2,double> >;
  template class PardisoInverse<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class PardisoInverse<Mat<3,3,double> >;
  template class PardisoInverse<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
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









#endif

