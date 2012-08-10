/* *************************************************************************/
/* File:   pardisoinverse.cpp                                              */
/* Author: Florian Bachinger                                               */
/* Date:   Feb. 2004                                                       */
/* *************************************************************************/

#include <la.hpp>

#ifdef USE_PARDISO

#define F77_FUNC(func)  func ## _

using ngbla::integer;

extern "C"
{
  /* PARDISO prototype. */

#ifdef USE_MKL

  void F77_FUNC(pardiso)
    (void * pt, integer * maxfct, integer * mnum, integer * mtype, integer * phase, integer * n, 
     double * a, integer * ia, integer * ja, integer * perm, integer * nrhs, integer * iparam, 
     integer * msglvl, double * b, double * x, integer * error);

#else

#ifdef USE_PARDISO400
extern  integer F77_FUNC(pardisoinit)
    (void *, integer *, integer *, integer *, double *, integer *);
#else
extern  integer F77_FUNC(pardisoinit)
    (void *, integer *, integer *);
#endif
  integer F77_FUNC(pardiso)
    (void * pt, integer * maxfct, integer * mnum, integer * mtype, integer * phase, integer * n, 
     double * a, integer * ia, integer * ja, integer * perm, integer * nrhs, integer * iparam, 
     integer * msglvl, double * b, double * x, integer * error);

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
    : SparseFactorization (a, ainner, acluster)
  { 
    static Timer timer("Pardiso Inverse");
    RegionTimer reg (timer);


    if (getenv ("PARDISOMSG"))
      pardiso_msg = 1;


    print = bool (pardiso_msg); // false;

    symmetric = asymmetric;
    compressed = false;

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
    compressed_height = height;
    cout << "NZE = " << a.NZE() << endl;
    rowstart.SetSize(compressed_height+1);
    indices.SetSize(a.NZE() * sqr(entrysize));
    matrix.SetSize (a.NZE() * sqr(entrysize));
    // matrix = new TSCAL[a.NZE() * entrysize * entrysize ];     

    *testout << "matrix.InverseTpye = " <<  a.GetInverseType() << endl;
    spd = ( a.GetInverseType() == PARDISOSPD ) ? 1 : 0;

    integer maxfct = 1, mnum = 1, phase = 12, nrhs = 1, msglevel = print, error;
    integer * params = const_cast <integer*> (&hparams[0]);

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
    //     retvalue = F77_FUNC(pardiso) ( pt, &maxfct, &mnum, &matrixtype, &phase, &compressed_height, 
    // 				   reinterpret_cast<double *>(matrix),
    // 				   rowstart, indices, NULL, &nrhs, params, &msglevel,
    // 				   NULL, NULL, &error );
#else
    int retvalue;
#ifdef USE_PARDISO400
    double dparm[64]; 
    integer solver = 0;
    F77_FUNC(pardisoinit) (pt,  &matrixtype, &solver, params, dparm, &retvalue); 
#else
    retvalue = 
      F77_FUNC(pardisoinit) (pt,  &matrixtype, params); 
#endif
    // cout << "init success" << endl;
    // cout << "retvalue = " << retvalue << endl;
#endif
    
    SetMatrixType();


    Array<int> icompress(a.Height());
    compress.SetSize(0);
    if (inner)
      {
	compressed = true;
	icompress = -1;

	for (int i = 0; i < a.Height(); i++)
	  if (inner->Test(i))
	    {
	      icompress[i] = compress.Size();
	      compress.Append (i);
	    }
	
	compressed_height = compress.Size() * entrysize;
	rowstart.SetSize(compressed_height+1);
      }


    if ( symmetric )
      {
	// --- transform lower left to upper right triangular matrix ---
	// 1.) build array 'rowstart':
	// (a) get nr. of entries for each row


	
	rowstart = 0;
	
	if (inner)
	  for (int i = 0; i < a.Height(); i++ )
	    {
	      if (inner->Test(i))
		for (int j = 0; j < a.GetRowIndices(i).Size(); j++ )
		  {
		    int col = a.GetRowIndices(i)[j];
		    int ccol = icompress[col];
		    
		    if ( i != col )
		      {
			if (inner->Test(col))
			  for (int k = 0; k < entrysize; k++ )
			    rowstart[ccol*entrysize+k+1] += entrysize;
		      }
		    else
		      {
			for (int k=0; k<entrysize; k++ )
			  rowstart[ccol*entrysize+k+1] += entrysize-k;
		      }
		  }
	    }

	else
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
	rowstart[0] = 0;
	for (int i = 1; i <= compressed_height; i++) 
	  rowstart[i] += rowstart[i-1];
	
	
	// 2.) build whole matrix:
	Array<int> counter(compressed_height);

	if (inner)
	  {
	    for (int i = 0; i < a.Height(); i++ )
	      {
		int ci = icompress[i];
		if (ci < 0) continue;
		for (int k = 0; k < entrysize; k++) counter[ci*entrysize+k]=0;

		for (int j = 0; j < a.GetRowIndices(i).Size(); j++)
		  {
		    int col = a.GetRowIndices(i)[j];
		    int ccol = icompress[col];

		    if (i != col)
		      {
			if (inner->Test(i) && inner->Test(col))
			  {
			    TM entry = a(i,col);
			    for (int k = 0; k < entrysize; k++)
			      for (int l = 0; l < entrysize; l++ )
				{
				  indices[rowstart[ccol*entrysize+k]+
					  counter[ccol*entrysize+k]] = ci*entrysize+l+1;
				  matrix[rowstart[ccol*entrysize+k]+
					 counter[ccol*entrysize+k]] = Access(entry,l,k);
				  counter[ccol*entrysize+k]++;
				}
			  }
		      }
		    else 
		      {
			if (inner->Test(i))
			  {
			    *testout << "i = " << i << ", ci = " << ci << " ccol = " << ccol << endl;
			    *testout << "rowstart = " << rowstart[ccol*entrysize] << "; counter = " << counter[ccol*entrysize] << endl;
			    TM entry = a(i,col);
			    for (int l = 0; l < entrysize; l++ )
			      for (int k = 0; k <= l; k++)
				{
				  indices[rowstart[ccol*entrysize+k]+
					  counter[ccol*entrysize+k]] = ci*entrysize+l+1;
				  matrix[rowstart[ccol*entrysize+k]+
					 counter[ccol*entrysize+k]] = Access(entry,l,k);
				  counter[ccol*entrysize+k]++;
				}
			  }
			/*
			else
			  {
			    // in the case of 'inner' or 'cluster': 1 on the diagonal for
			    // unused dofs.
			    for (int l=0; l<entrysize; l++ )
			      {
				indices[rowstart[col*entrysize+l]+
					counter[col*entrysize+l]] = i*entrysize+l+1;
				matrix[rowstart[col*entrysize+l]+
				       counter[col*entrysize+l]] = 1;
				counter[col*entrysize+l]++;
			      }
			  }
			*/
		      }
		  }
	      }
	  }	
	else
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
				       counter[col*entrysize+k]] = i*entrysize+l+1;
			      matrix[ rowstart[col*entrysize+k]+
				      counter[col*entrysize+k]] = Access(entry,l,k);
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
				   counter[col*entrysize+k]] = i*entrysize+l+1;
			  matrix[ rowstart[col*entrysize+k]+
				  counter[col*entrysize+k]] = Access(entry,l,k);
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
				 counter[col*entrysize+l]] = i*entrysize+l+1;
			matrix[ rowstart[col*entrysize+l]+
				counter[col*entrysize+l]] = 1;
			counter[col*entrysize+l]++;
		      }
		  }
	      }
	  }
	    
	// pardiso is 1-based
	for (int i = 0; i <= compressed_height; i++)
	  rowstart[i]++;
	    
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

	// delete [] counter;
      }
    else // non-symmetric
      {
	if (print)
	  cout << "non-symmetric pardiso, cluster = " << cluster << ", entrysize = " << entrysize << endl;

	// for non-symmetric matrices:
	int counter = 0;

	for (int i = 0; i < a.Height(); i++ )
	  {
	    int rowelems = 0;
	    int rowsize = a.GetRowIndices(i).Size();
	    int ci = inner ? icompress[i] : i;

	    // get effective number of elements in the row

	    if ( inner )
	      {
		if (inner -> Test(i))
		  for (int j = 0; j < rowsize; j++)
		    {
		      int col = a.GetRowIndices(i)[j];
		      if (inner->Test(col))
			rowelems += entrysize;
		    }
	      }
            else if (cluster)
	      {
		for (int j = 0; j < rowsize; j++ )
		  {
		    int col = a.GetRowIndices(i)[j];
		    if ( (*cluster)[i] == (*cluster)[col] && (*cluster)[i])
		      rowelems += entrysize;
		    else if ( i == col )
		      rowelems++;
		  }
	      }
	    else 
	      rowelems = rowsize * entrysize;

	    if (ci != -1)
	      for (int k = 0; k < entrysize; k++ )
		rowstart[ci*entrysize+k] = counter+rowelems*k+1;
	    

            if (inner)
              {
                for (int j=0; j<rowsize; j++ )
                  {
                    int col = a.GetRowIndices(i)[j];
                    int ccol = icompress[col];

                    if (inner->Test(i) && inner->Test(col))
                      {
                        TM entry = a(i,col);
                        for (int k=0; k<entrysize; k++ )
                          for (int l=0; l<entrysize; l++ )
                            {
                              indices[counter+rowelems*k+l] = ccol*entrysize+l+1;
                              matrix[counter+rowelems*k+l] = Access(entry,k,l);
                            }
                        counter+=entrysize;
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
		  const TM & entry = a.GetRowValues(i)(j);
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
	rowstart[compressed_height] = counter+1;
	
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

    nze = rowstart[compressed_height];
    
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
    cout << IM(3) << "call pardiso ..." << flush;


    


    // retvalue = 
    F77_FUNC(pardiso) ( pt, &maxfct, &mnum, &matrixtype, &phase, &compressed_height, 
			reinterpret_cast<double *>(&matrix[0]),
			&rowstart[0], &indices[0], NULL, &nrhs, params, &msglevel,
			NULL, NULL, &error );
    
    cout << IM(3) << " done" << endl;

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
    static Timer timer("Pardiso Solve");
    RegionTimer reg (timer);

    FlatVector<TVX> fx = x.FV<TVX> ();
    FlatVector<TVX> fy = y.FV<TVX> ();

    integer maxfct = 1, mnum = 1, phase = 33, msglevel = 0, error;
    integer nrhs = 1;

    integer * params = const_cast <integer*> (&hparams[0]);

    if (compressed)
      {
	Vector<TVX> hx(compress.Size());
	Vector<TVX> hy(compress.Size());

	/*
	for (int i = 0; i < compress.Size(); i++)
	  hx(i) = fx(compress[i]);
	*/
	hx = fx(compress);

	F77_FUNC(pardiso) ( const_cast<integer*>(pt), &maxfct, &mnum, 
			    const_cast<integer*>(&matrixtype),
			    &phase, const_cast<integer*>(&compressed_height), 
			    reinterpret_cast<double *>(&matrix[0]),
			    &rowstart[0], &indices[0],
			    NULL, &nrhs, params, &msglevel,
			    static_cast<double *>(hx.Data()), 
			    static_cast<double *>(hy.Data()), &error );

	fy = 0;
	fy(compress) = hy;
	/*
	for (int i = 0; i < compress.Size(); i++)
	  fy(compress[i]) = hy(i);
	*/
      }
    else
      {
	F77_FUNC(pardiso) ( const_cast<integer *>(pt), &maxfct, &mnum, 
			    const_cast<integer *>(&matrixtype),
			    &phase, const_cast<integer *>(&compressed_height), 
			    reinterpret_cast<double *>(&matrix[0]),
			    &rowstart[0], &indices[0],
			    NULL, &nrhs, params, &msglevel,
			    static_cast<double *>(fx.Data()), 
			    static_cast<double *>(fy.Data()), &error );

	for (int i = 0; i < compressed_height/entrysize; i++)
	  if (!used.Test(i))
	    for (int j=0; j<entrysize; j++ ) fy(i*entrysize+j) = 0.0;
      }

    if ( error != 0 )
      cout << "Apply Inverse: PARDISO returned error " << error << "!" << endl;

  }
  
  




  template<>
  void PardisoInverse<double,Complex,Complex> :: 
  Mult (const BaseVector & x, BaseVector & y) const
  {
    static Timer timer ("Pardiso Solve");
    RegionTimer reg (timer);

    // FlatVector<TVX> fx = x.FV<TVX> ();
    FlatVector<TVX> fy = y.FV<TVX> ();


    FlatMatrixFixWidth<2> hx(x.Size(), (double*)x.Memory());
    Matrix<> tx(2, x.Size());
    tx = Trans (hx);

    FlatMatrixFixWidth<2> hy(y.Size(), (double*)y.Memory());
    Matrix<> ty(2, y.Size());

      
    integer maxfct = 1, mnum = 1, phase = 33, msglevel = 0, error;
    integer nrhs = 2;

    integer * params = const_cast <integer*> (&hparams[0]);

    F77_FUNC(pardiso) ( const_cast<integer *>(pt), 
			&maxfct, &mnum, const_cast<integer *>(&matrixtype),
			&phase, const_cast<integer *>(&compressed_height), 
			reinterpret_cast<double *>(&matrix[0]),
			&rowstart[0], &indices[0],
			NULL, &nrhs, params, &msglevel, &tx(0,0), &ty(0,0),
			&error );

    if ( error != 0 )
      cout << "Apply Inverse: PARDISO returned error " << error << "!" << endl;

    hy = Trans (ty);

    for (int i=0; i<compressed_height/entrysize; i++)
      if (!used.Test(i))
	for (int j=0; j<entrysize; j++ ) fy(i*entrysize+j) = 0.0;
  }
  
  






  template<class TM, class TV_ROW, class TV_COL>
  ostream & PardisoInverse<TM,TV_ROW,TV_COL> :: Print (ostream & ost) const
  {
    cout << "PardisoInverse::Print not implemented!" << endl;
    return ost; 
  }


  template<class TM, class TV_ROW, class TV_COL>
  PardisoInverse<TM,TV_ROW,TV_COL> :: ~PardisoInverse()
  {
    integer maxfct = 1, mnum = 1, phase = -1, nrhs = 1, msglevel = 1, error;
    integer * params = const_cast <integer*> (&hparams[0]);

    //    cout << "call pardiso (clean up) ..." << endl;
    F77_FUNC(pardiso) ( pt, &maxfct, &mnum, &matrixtype, &phase, &compressed_height, NULL,
			&rowstart[0], &indices[0], NULL, &nrhs, params, &msglevel,
			NULL, NULL, &error );
    if ( error != 0 )
      cout << "Clean Up: PARDISO returned error " << error << "!" << endl;

    // delete [] rowstart;
    // delete [] indices;
    // delete [] matrix;
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

