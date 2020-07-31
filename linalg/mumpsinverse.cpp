/* *************************************************************************/
/* File:   mumpsinverse.cpp                                                */
/* Author: Joachim Schoeberl                                               */
/* Date:   May. 2009                                                       */
/* *************************************************************************/

// #define DEBUG
#ifdef USE_MUMPS

#include <la.hpp>
#include "mumpsinverse.hpp"

#include <comp.hpp>
// #include <parallelngs.hpp>

namespace ngla
{
  using namespace ngcomp;

#define JOB_INIT -1
#define JOB_END -2

#define JOB_ANALYSIS 1
#define JOB_FACTOR 2
#define JOB_SOLVE 3

#define USE_COMM_WORLD -987654

  
  template <class TM, class TV_ROW, class TV_COL>
  MumpsInverse<TM,TV_ROW,TV_COL> :: 
  MumpsInverse (const SparseMatrix<TM,TV_ROW,TV_COL> & a, 
                shared_ptr<BitArray> ainner,
                shared_ptr<const Array<int>> acluster,
                bool asymmetric)
    : comm(MPI_COMM_NULL, false)
  { 
    static Timer timer ("Mumps Inverse");
    static Timer timer_analysis ("Mumps Inverse - analysis");
    static Timer timer_factor ("Mumps Inverse - factor");
    RegionTimer reg (timer);


    symmetric = asymmetric;
    inner = ainner;
    cluster = acluster;

    auto pds = a.GetParallelDofs();
    if ( (pds != nullptr) && // if we are on an "only-me comm", take it
	 (pds->GetCommunicator().Size() == 1) )
      { comm = pds->GetCommunicator(); }
    else { // otherwise make an "only-me comm" from world-communicator
      NgMPI_Comm wcomm(MPI_COMM_WORLD);
      Array<int> procs = { wcomm.Rank() };
      comm = wcomm.SubCommunicator(procs);
    }

    int ntasks = comm.Size();
    int id = comm.Rank();
    
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
		FlatArray<int> rowind = a.GetRowIndices(i);

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




    for (int i = 0; i < 40; i++)
      mumps_id.icntl[i] = 0;


    mumps_id.job =JOB_INIT; 
    mumps_id.par = (ntasks == 1) ? 1 : 0;
    mumps_id.sym = symmetric ? 1 : 0;
    // mumps_id.comm_fortran=USE_COMM_WORLD;
    mumps_id.comm_fortran = MPI_Comm_c2f (comm);
    mumps_trait<TSCAL>::MumpsFunction (&mumps_id);

    // cout << IM(1) << "MUMPS version number is " << mumps_id.version_number << endl;


    /* Define the problem on the host */
    mumps_id.n   = height; 
    mumps_id.nz  = nze;
    mumps_id.irn = row_indices;
    mumps_id.jcn = col_indices;

    /*
      if (id == 0)
      {
      cout << "Mumps predefined values: ";
      for (int j = 0; j < 40; j++)
      cout << "ICNTL(" << j+1 << ") = " << mumps_id.icntl[j] << endl;
      }
    */

    mumps_id.icntl[0]=-1; 
    mumps_id.icntl[1]=-1; 
    mumps_id.icntl[2]=-1; 
    mumps_id.icntl[3]=0;
    mumps_id.icntl[6]=7;   // 0..min deg, 3..scotch 5..metis, 7..default
    mumps_id.icntl[12]=1;  // not using scalapck for root schur complement
    mumps_id.icntl[13]=60; // memory increase (in %) due to error -9
    mumps_id.icntl[27]=0;  // 0..default,  2..parallel analysis
    mumps_id.icntl[28]=0;  // 0..auto, 1..ptscotch 2..parmetis

    // mumps_id.comm_fortran=USE_COMM_WORLD;
    mumps_id.comm_fortran = MPI_Comm_c2f (comm);
    mumps_id.job = JOB_ANALYSIS;


    if (id == 0)
      cout << "analysis ... " << flush;

    timer_analysis.Start();
    mumps_trait<TSCAL>::MumpsFunction (&mumps_id);
    timer_analysis.Stop();

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

    MPI_Barrier (comm);

    timer_factor.Start();
    mumps_trait<TSCAL>::MumpsFunction (&mumps_id);
    timer_factor.Stop();

    if (mumps_id.infog[0] != 0)
      {
	cout << " factorization done" << endl;
	cout << "error-code = " << mumps_id.infog[0] << endl;
	cout << "info(1) = " << mumps_id.info[0] << endl;
	cout << "info(2) = " << mumps_id.info[1] << endl;
      }


    /*
      if ( error != 0 )
      {
      cout << "Setup and Factorization: Mumps returned error " << error << "!" << endl;
      throw Exception("MumpsInverse: Setup and Factorization failed.");
      }
    */



    
    if (id == 0)
      cout << " done " << endl;
    delete [] colstart;
    delete [] counter;
    // delete [] rhs;    


    delete [] col_indices;
    delete [] row_indices;
    delete [] matrix;
  }
  
  

  template <class TM, class TV_ROW, class TV_COL>
  void MumpsInverse<TM,TV_ROW,TV_COL> :: 
  Mult (const BaseVector & x, BaseVector & y) const
  {
    int id = comm.Rank();

    static int timer = NgProfiler::CreateTimer ("Mumps mult inverse");
    NgProfiler::RegionTimer reg (timer);


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
  }
  

  template <class TM, class TV_ROW, class TV_COL>
  MumpsInverse<TM,TV_ROW,TV_COL> :: ~MumpsInverse()
  {
    mumps_id.job=JOB_END; 
    mumps_trait<TSCAL>::MumpsFunction (&mumps_id);
  }



















  

  template <class TM, class TV>
  ParallelMumpsInverse<TM,TV> :: 
  ParallelMumpsInverse (const BaseSparseMatrix & ba, 
			shared_ptr<BitArray> inner,
			shared_ptr<const Array<int>> cluster,
			shared_ptr<ParallelDofs> pardofs,
			bool asymmetric)
    : BaseMatrix(pardofs)
  { 
    static Timer timer ("Mumps Inverse");
    static Timer timer_analysis ("Mumps Inverse - analysis");
    static Timer timer_factor ("Mumps Inverse - factor");
    RegionTimer reg (timer);

    const SparseMatrixTM<TM> & a = dynamic_cast<const SparseMatrixTM<TM> &> (ba);

    
    symmetric = asymmetric;
    // symmetric = true;
    // inner = ainner;
    // cluster = acluster;

    cout << IM(1) << "Mumps Parallel inverse, symmetric = " << symmetric << endl;

    NgMPI_Comm comm = pardofs->GetCommunicator();
    int ntasks = comm.Size();
    int id = comm.Rank();

    if (id == 0)
      {
	if ( ( inner && inner->Size() < a.Height() ) ||
	     ( cluster && cluster->Size() < a.Height() ) )
	  {
	    cout << "Mumps: Size of inner/cluster does not match matrix size!" << endl;
	    throw Exception("Invalid parameters inner/cluster. Thrown by ParallelMumpsInverse.");
	  }
	
	if ( int( mat_traits<TM>::WIDTH) != int(mat_traits<TM>::HEIGHT) )
	  {
	    cout << "Mumps: Each entry in the square matrix has to be a square matrix!" << endl;
	    throw Exception("No Square Matrix. Thrown by ParallelMumpsInverse.");
	  }
      }

    // find global dofs

    num_globdofs = 0;   // valid on id=0






   // consistent enumeration (new version)
    
    Array<int> global_nums;
    int num_glob_dofs;
    pardofs -> EnumerateGlobally (inner, global_nums, num_glob_dofs);

    int ndof = pardofs->GetNDofLocal();


    /*
    Array<int> global_nums(ndof);
    global_nums = -1;
    int num_master_dofs = 0;
    for (int i = 0; i < ndof; i++)
      if (pardofs -> IsMasterDof (i) && (!inner || (inner && inner->Test(i))))
	global_nums[i] = num_master_dofs++;
    
    Array<int> first_master_dof(ntasks);
    MPI_Allgather (&num_master_dofs, 1, MPI_INT, 
		   &first_master_dof[0], 1, MPI_INT, 
		   pardofs -> GetCommunicator());
    
    int num_glob_dofs = 0;
    for (int i = 0; i < ntasks; i++)
      {
	int cur = first_master_dof[i];
	first_master_dof[i] = num_glob_dofs;
	num_glob_dofs += cur;
      }
    
    for (int i = 0; i < ndof; i++)
      if (global_nums[i] != -1)
	global_nums[i] += first_master_dof[id];
    
    ScatterDofData (global_nums, *pardofs);
    */
    

    // copy to old variables ...
    num_globdofs = num_glob_dofs;
    loc2glob = global_nums;
    for (int row = 0; row < ndof; row++)
      if (!inner || inner->Test(row))
	select.Append (row);
    
    
    /*
      // test dofs, only for ParallelMeshDofs
    const ParallelMeshDofs & pmdofs = dynamic_cast<const ParallelMeshDofs&> (*pardofs);
    for (int row = 0; row < ndof; row++)
      if (!inner || inner->Test(row))
	if (global_nums[row] < 0)
	  {
	    Node node = pmdofs.GetDofNodes()[row];
	    cout << "illegal gobal num, id = " << id << ", localnum = " << row << 
	      ", on " << node << endl;
	    
	    Array<int> procs;
	    pmdofs.GetMeshAccess().GetDistantProcs (node, procs);
	    cout << "procs = " << procs << endl;
	  }
    */


    entrysize = mat_traits<TM>::HEIGHT; 
    iscomplex = mat_traits<TM>::IS_COMPLEX;

    Array<int> col_indices;
    Array<int> row_indices;
    Array<TSCAL> matrix;

    if (id != 0)
      {
	height = a.Height() * entrysize;
	col_indices.SetSize (a.NZE() * entrysize * entrysize);
	row_indices.SetSize (a.NZE() * entrysize * entrysize);
	matrix.SetSize (a.NZE() * entrysize * entrysize);


	if ( symmetric )
	  {
	    int ii = 0;
	    for (int i = 0; i < a.Height(); i++ )
	      {
		FlatArray<int> rowind = a.GetRowIndices(i);
		FlatVector<TM> values = a.GetRowValues(i);
		
		for (int j = 0; j < rowind.Size(); j++ )
		  {
		    int col = rowind[j];
		    
		    if (  (!inner && !cluster) ||
			  (inner && (inner->Test(i) && inner->Test(col) ) ) ||
			  (!inner && cluster &&
			   ((*cluster)[i] == (*cluster)[col] 
			    && (*cluster)[i] ))  )
		      {
			for (int l = 0; l < entrysize; l++ )
			  for (int k = 0; k < entrysize; k++)
			    {
			      if (i == col && k > l) continue;

			      int rowi = loc2glob[i]*entrysize+l+1;
			      int coli = loc2glob[col]*entrysize+k+1;

			      if (rowi >= coli)
				{
				  col_indices[ii] = coli;
				  row_indices[ii] = rowi;
				  matrix[ii] = Access(values[j],l,k);
				  ii++;
				}
			      else
				{
				  col_indices[ii] = rowi;
				  row_indices[ii] = coli;
				  matrix[ii] = Access(values[j],l,k);
				  ii++;
				}
			    }
		      }
		  }
	      }
	    nze = ii;
	  }
	else
	  {
	    int ii = 0;
	    for (int i = 0; i < a.Height(); i++ )
	      {
		FlatArray<int> rowind = a.GetRowIndices(i);
		FlatVector<TM> values = a.GetRowValues(i);
		
		for (int j = 0; j < rowind.Size(); j++ )
		  {
		    int col = rowind[j];
		    
		    /*
		    if (  (!inner && !cluster) ||
			  (inner && (inner->Test(i) && inner->Test(col) ) ) ||
			  (!inner && cluster &&
			   ((*cluster)[i] == (*cluster)[col] 
			    && (*cluster)[i] ))  )
		    */
		    if (loc2glob[i] != -1 && loc2glob[col] != -1)
		      {
			for (int l = 0; l < entrysize; l++ )
			  for (int k = 0; k < entrysize; k++)
			    {
			      int rowi = loc2glob[i]*entrysize+l+1;
			      int coli = loc2glob[col]*entrysize+k+1;

			      col_indices[ii] = coli;
			      row_indices[ii] = rowi;
			      matrix[ii] = Access(values[j],l,k);
			      ii++;
			    }
		      }
		  }
	      }
	    nze = ii;


#ifdef OLDNONSYM	   



	    // --- transform matrix to compressed column storage format ---

	    int * colstart = new int[height+1];
	    int * counter = new int[height];
	    for ( int i = 0; i < height; i++ ) 
	      counter[i] = 0;
	    for ( int i = 0; i < height; i++ ) 
	      colstart[i+1] = 0;
	    


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
	    /*
	    int * col_indices = new int[a.NZE() * entrysize * entrysize ];
	    int * row_indices = new int[a.NZE() * entrysize * entrysize ];
	    TSCAL * matrix = new TSCAL[a.NZE() * entrysize * entrysize ];      
	    */
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
#endif
	  }
      }
    else
      {
	height = 0;
	nze = 0;
      }

    /*
    for (int i = 0; i < 40; i++)
      mumps_id.icntl[i] = 0;
    */


    /*
    *testout << "mumps matrix: n = " << num_globdofs << ", nz = " << nze << endl;
    *testout << "loc2glob = " << loc2glob << endl;
    for (int i = 0; i < nze; i++)
      *testout << "a(" << row_indices[i] << "," << col_indices[i] << ") = " << matrix[i] << endl;
      */
    
    mumps_id.job = JOB_INIT; 
    mumps_id.par = (ntasks == 1) ? 1 : 0;
    mumps_id.sym = symmetric ? 2 : 0;    // 1 .. spd, 2 .. general symmetric
    //mumps_id.comm_fortran=USE_COMM_WORLD;
    mumps_id.comm_fortran = MPI_Comm_c2f (comm);


    MumpsFunction (mumps_id);

    // cout << IM(0) << "MUMPS version number is " << mumps_id.version_number << endl;

    
    /* distributed matrix definition */
    mumps_id.n   = num_globdofs * entrysize;  // only on host
    mumps_id.nz_loc  = nze;
    mumps_id.irn_loc = row_indices.Addr(0);
    mumps_id.jcn_loc = col_indices.Addr(0);


    mumps_id.icntl[0]=-1; 
    mumps_id.icntl[1]=-1; 
    mumps_id.icntl[2]=-1; 
    mumps_id.icntl[3]=1;

    if (getenv ("MUMPSMSG"))
      {
	int level = atoi (getenv("MUMPSMSG"));
	mumps_id.icntl[0]= 6; 
	mumps_id.icntl[1]= 6; 
	mumps_id.icntl[2]= 6; 
	mumps_id.icntl[3]= level;
      }


    /*
    // mumps_id.icntl[7]=7;   // BUG (??) 0..min deg, 3..scotch 5..metis, 7..default
    mumps_id.icntl[6]=7;   // 0..min deg, 3..scotch 5..metis, 7..default
    mumps_id.icntl[12]=1;  // 0..do use, 1..not using scalapck for root schur complement
    mumps_id.icntl[13]=200; // memory increase (in %) due to error -9
    mumps_id.icntl[17]=3;  // parallel input
    mumps_id.icntl[27]=0;  // 0..default, 1..seq, 2..parallel analysis
    mumps_id.icntl[28]=0;  // 0..auto, 2..parmetis
    */

    mumps_id.icntl[6]=7;   // 0..min deg, 3..scotch 5..metis, 7..default
    mumps_id.icntl[12]=0;  // 0..do use, 1..not using scalapck for root schur complement
    mumps_id.icntl[13]=50; // memory increase (in %) due to error -9
    mumps_id.icntl[17]=3;  // parallel input
    mumps_id.icntl[27]=0;  // 0..default, 1..seq, 2..parallel analysis
    mumps_id.icntl[28]=2;  // 0..auto, 2..parmetis


    // mumps_id.comm_fortran=USE_COMM_WORLD;
    mumps_id.comm_fortran = MPI_Comm_c2f (comm);
    mumps_id.job = JOB_ANALYSIS;

    cout << IM(1) << "analysis ... " << flush;

    timer_analysis.Start();
    MumpsFunction (mumps_id);
    timer_analysis.Stop();

    // cout << "num floating-point ops = " << mumps_id.rinfog[0] << endl;
    if (mumps_id.infog[0])
      {
	cout << "analysis done" << endl;
	cout << "error-code = " << mumps_id.infog[0] << flush;
      }

    /*
    cout << "use ordering = " << mumps_id.infog[6] << endl;
    cout << "parallel ordering = " << mumps_id.infog[31] << endl;
    */

    mumps_id.a_loc = (typename mumps_trait<TSCAL>::MUMPS_TSCAL*) (TSCAL*)matrix.Addr(0);
    mumps_id.job = JOB_FACTOR;
    
    cout << IM(1) << "factor ... " << flush;

    timer_factor.Start();
    MumpsFunction (mumps_id);
    timer_factor.Stop();

    if (mumps_id.infog[0] != 0)
      {
	cout << " factorization done" << endl;
	cout << "error-code = " << mumps_id.infog[0] << endl;
	cout << "info(1) = " << mumps_id.info[0] << endl;
	cout << "info(2) = " << mumps_id.info[1] << endl;
      }

    cout << IM(1) << " done " << endl;
    
  }
  
  

  template <class TM, class TV>
  void ParallelMumpsInverse<TM,TV> :: 
  Mult (const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("Parallelmumps mult inverse");
    RegionTimer reg (timer);

    x.Distribute();
    y.SetParallelStatus (CUMULATED);


    NgMPI_Comm comm = paralleldofs->GetCommunicator();
    int ntasks = comm.Size();
    int id = comm.Rank();


    if (id != 0)
      {
	FlatVector<TV> fx = x.FV<TV>();
	FlatVector<TV> fy = y.FV<TV>();

	Array<TV> hx(select.Size());
	for (int i = 0; i < select.Size(); i++)
	  hx[i] = fx(select[i]);
	Array<int> select_loc2glob(select.Size());
	for (int i = 0; i < select.Size(); i++)
	  select_loc2glob[i] = loc2glob[select[i]];
	
	comm.Send (select_loc2glob, 0, MPI_TAG_SOLVE);
	comm.Send (hx, 0, MPI_TAG_SOLVE);

	MUMPS_STRUC_C & ncid = const_cast<MUMPS_STRUC_C&> (mumps_id);
	ncid.job = JOB_SOLVE;
	mumps_trait<TSCAL>::MumpsFunction (&ncid);

	comm.Send (select_loc2glob, 0, MPI_TAG_SOLVE);
	comm.Recv (hx, 0, MPI_TAG_SOLVE);

	y = 0;
	for (int i = 0; i < select.Size(); i++)
	  fy(select[i]) = hx[i];
      }
    else
      {
	Vector<TV> rhs(num_globdofs);
	rhs = 0.0;
	
	for (int src = 1; src < ntasks; src++)
	  {
	    Array<int> loc2glob;
	    Array<TV> hx;
	    comm.Recv (loc2glob, src, MPI_TAG_SOLVE);
	    comm.Recv (hx, src, MPI_TAG_SOLVE);
	    for (int j = 0; j < loc2glob.Size(); j++)
	      rhs(loc2glob[j]) += hx[j];
	  } 
	
	MUMPS_STRUC_C & ncid = const_cast<MUMPS_STRUC_C&> (mumps_id);

	// ncid.rhs = (typename mumps_trait<TSCAL>::MUMPS_TSCAL*) &rhs(0);
	ncid.rhs = (typename mumps_trait<TSCAL>::MUMPS_TSCAL*) rhs.Data();
	ncid.job = JOB_SOLVE;

	mumps_trait<TSCAL>::MumpsFunction (&ncid);

	for (int src = 1; src < ntasks; src++)
	  {
	    Array<int> loc2glob;
	    comm.Recv (loc2glob, src, MPI_TAG_SOLVE);
	    Array<TV> hx(loc2glob.Size());

	    for (int j = 0; j < loc2glob.Size(); j++)
	      hx[j] = rhs(loc2glob[j]);
	    comm.Send (hx, src, MPI_TAG_SOLVE);
	  }
      }

  }
  

  template <class TM, class TV>
  ParallelMumpsInverse<TM,TV> :: ~ParallelMumpsInverse()
  {
    mumps_id.job=JOB_END; 
    mumps_trait<TSCAL>::MumpsFunction (&mumps_id);
  }


  template <class TM, class TV>
  AutoVector ParallelMumpsInverse<TM,TV> :: CreateVector () const
  {
    return make_unique<ParallelVVector<TV>> (paralleldofs->GetNDofLocal(), paralleldofs);
  }

  template <class TM, class TV>
  AutoVector ParallelMumpsInverse<TM,TV> :: CreateRowVector () const
  {
    return make_unique<ParallelVVector<TV>> (paralleldofs->GetNDofLocal(), paralleldofs);
  }

  template <class TM, class TV>
  AutoVector ParallelMumpsInverse<TM,TV> :: CreateColVector () const
  {
    return make_unique<ParallelVVector<TV>> (paralleldofs->GetNDofLocal(), paralleldofs);
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






  


  template class ParallelMumpsInverse<double>;
  template class ParallelMumpsInverse<Complex>;
  template class ParallelMumpsInverse<double,Complex>;
#if MAX_SYS_DIM >= 1
  template class ParallelMumpsInverse<Mat<1,1,double> >;
  template class ParallelMumpsInverse<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class ParallelMumpsInverse<Mat<2,2,double> >;
  template class ParallelMumpsInverse<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class ParallelMumpsInverse<Mat<3,3,double> >;
  template class ParallelMumpsInverse<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class ParallelMumpsInverse<Mat<4,4,double> >;
  template class ParallelMumpsInverse<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class ParallelMumpsInverse<Mat<5,5,double> >;
  template class ParallelMumpsInverse<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class ParallelMumpsInverse<Mat<6,6,double> >;
  template class ParallelMumpsInverse<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class ParallelMumpsInverse<Mat<7,7,double> >;
  template class ParallelMumpsInverse<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class ParallelMumpsInverse<Mat<8,8,double> >;
  template class ParallelMumpsInverse<Mat<8,8,Complex> >;
#endif

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



