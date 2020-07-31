/*********************************************************************/
/* File:   parallel_matrices.cpp                                     */
/* Author: Joachim Schoeberl                                         */
/* Date:   June 2011                                                 */
/*********************************************************************/

 
#include <la.hpp>
#include "../linalg/mumpsinverse.hpp"
#include <parallelngs.hpp>


namespace ngla
{

#ifdef PARALLEL
  
  template <typename TM> AutoVector MasterInverse<TM> :: CreateRowVector () const
  { return make_unique<ParallelVVector<double>> (paralleldofs->GetNDofLocal(), paralleldofs); }
  template <typename TM> AutoVector MasterInverse<TM> :: CreateColVector () const
  { return make_unique<ParallelVVector<double>> (paralleldofs->GetNDofLocal(), paralleldofs); }
  
  template <typename TM>
  MasterInverse<TM> :: MasterInverse (const SparseMatrixTM<TM> & mat, 
				      shared_ptr<BitArray> subset, 
				      shared_ptr<ParallelDofs> hpardofs)
    
    : BaseMatrix(hpardofs), loc2glob(hpardofs -> GetCommunicator().Size())
  {
    inv = nullptr;
    
    auto & comm = paralleldofs->GetCommunicator();
    int id = comm.Rank();
    int ntasks = comm.Size();

    // consistent enumeration
    
    int ndof = paralleldofs->GetNDofLocal();

    Array<int> global_nums(ndof);
    global_nums = -1;
    int num_master_dofs = 0;
    for (int i = 0; i < ndof; i++)
      if (paralleldofs -> IsMasterDof (i) && (!subset || (subset && subset->Test(i))))
	global_nums[i] = num_master_dofs++;
    

    Array<int> first_master_dof(ntasks);
    MPI_Allgather (&num_master_dofs, 1, MPI_INT, 
		   &first_master_dof[0], 1, MPI_INT, 
		   comm);

    
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

    paralleldofs -> ScatterDofData (global_nums);


    /*
    cout << "TESTING" << endl;
    for (int i = 0; i < ndof; i++)
      *testout << "dof" << i << " = " << pardofs->GetDofNodes()[i] << " free = " << subset->Test(i) << endl;
    *testout << "global_nums = " << endl << global_nums << endl;


    const MeshAccess & ma = pardofs->GetMeshAccess();

    for (int dest = 0; dest < ntasks; dest++)
      {
	*testout << "shared dofs with proc " << dest << endl;
	*testout << "exdofs " << endl << pardofs->GetExchangeDofs (dest) << endl;
	Array<int> procs;
	for (int i = 0; i < ma.GetNFaces(); i++)
	  {
	    ma.GetDistantProcs (Node(NT_FACE, i), procs);
	    if (procs.Contains (dest))
	      {
		Array<int> vnums;
		ma.GetFacePNums (i, vnums);
		int v1 = ma.GetGlobalNodeNum (Node(NT_VERTEX, vnums[0]));
		int v2 = ma.GetGlobalNodeNum (Node(NT_VERTEX, vnums[1]));
		int v3 = ma.GetGlobalNodeNum (Node(NT_VERTEX, vnums[2]));
		*testout << "F" << i << ": " << v1 << "-" << v2 << "-" << v3 << endl;
	      }
	  }
      }
    */


    /*
    for (int i = 0; i < ma.GetNEdges(); i++)
      {
	*testout << "local edge " << i << endl;
	int v1, v2;
	ma.GetEdgePNums (i, v1, v2);
	*testout << "loc pnts = " << v1 << "-" << v2 << ", glob pnts = " 
		 << ma.GetGlobalNodeNum (Node(NT_VERTEX, v1)) << "-" 
		 << ma.GetGlobalNodeNum (Node(NT_VERTEX, v2)) << endl;
	Array<int> procs;
	ma.GetDistantProcs (Node(NT_EDGE, i), procs);
	*testout << "dist procs = " << procs << endl;
      }

    for (int i = 0; i < ndof; i++)
      if ( (global_nums[i] == -1) && (!subset || (subset && subset->Test(i))))
	{
	  Node node = pardofs->GetDofNodes()[i];
	  cerr << "global enumeration problem, dof = " << i << ", node = " << node << endl;
	  *testout << "global enumeration problem, dof = " << i << ", node = " << node << endl;
	}
    */

    
    if (id != 0)
      {
	// const MeshAccess & ma = nodaldofs -> GetMeshAccess();

	int ndof = paralleldofs->GetNDofLocal();

	Array<int> rows, cols;
	Array<TM> vals;

	for (int row = 0; row < ndof; row++)
	  if (!subset || subset->Test(row))
	    select.Append (row);
	
	for (int row = 0; row < mat.Height(); row++)
	  if (!subset || subset->Test(row))
	    {
	      FlatArray<int> rcols = mat.GetRowIndices(row);
	      FlatVector<TM> rvals = mat.GetRowValues(row);

	      for (int j = 0; j < rcols.Size(); j++)
		if (!subset || subset->Test(rcols[j]))
		  {
                    /*
		    *testout << "send (" << row << "," << rcols[j] 
			     << "), global = (" <<  global_nums[row] << "," << global_nums[rcols[j]] << ")" << endl;
                    */
		    rows.Append (global_nums[row]);
		    cols.Append (global_nums[rcols[j]]);
		    vals.Append (rvals[j]);
		  }
	    }

	comm.Send (rows, 0, MPI_TAG_SOLVE);
	comm.Send (cols, 0, MPI_TAG_SOLVE);
	comm.Send (vals, 0, MPI_TAG_SOLVE);
	comm.Send (global_nums, 0, MPI_TAG_SOLVE);

#ifdef USE_MUMPS
	if (mat.GetInverseType() == MUMPS)
	  {
	    SparseMatrixSymmetric<TM> dummy_matrix (1,1);
	    inv = make_shared<MumpsInverse<TM>> (dummy_matrix, nullptr, nullptr, true);
	  }
	else
#endif 
	  inv = nullptr;
      }

    else
      {
	cout << IM(3) << "create masterinverse" << endl;

	bool symmetric = (dynamic_cast<const SparseMatrixSymmetric<TM>*>(&mat) != NULL);
	cout << IM(4) << "symmetric? " << symmetric << endl;

	Array<int> rows, cols;
	Array<TM> vals;
	HashTable<INT<1>, int> ht_globdofs(100000);
	// int num_globdofs = 0; 

	cout << IM(4) << "collect data" << flush;

	for (int src = 1; src < ntasks; src++)
	  {
	    Array<int> hrows, hcols;
	    Array<TM> hvals;
	    Array<int> hglobid;

	    comm.Recv (hrows, src, MPI_TAG_SOLVE);
	    comm.Recv (hcols, src, MPI_TAG_SOLVE);
	    comm.Recv (hvals, src, MPI_TAG_SOLVE);
	    comm.Recv (hglobid, src, MPI_TAG_SOLVE);

            /*
	    *testout << "got from P" << src << ":" << endl
		     << "rows " << endl << hrows << endl
		     << "cols " << endl << hcols << endl;
            */
            
	    for (int i = 0; i < hglobid.Size(); i ++)
	      {
		if (hglobid[i] == -1) continue;
		loc2glob.Add (src, hglobid[i]);
	      }

	    for (int i = 0; i < hrows.Size(); i++)
	      {
		if (hrows[i] < 0 || hcols[i] < 0) cerr << "Illegal (row/col)" << endl;
		rows.Append (hrows[i]);
		cols.Append (hcols[i]);
		vals.Append (hvals[i]);
	      }

	    cout << IM(5) << "\rmaster: got data from " << src << flush;
	  }
	cout << IM(5) << endl;


	cout << IM(4) << "now build graph" << endl;

	// build matrix
	DynamicTable<int> graph(num_glob_dofs);
	cout << IM(3) << "n = " << num_glob_dofs << endl;
	for (int i = 0; i < rows.Size(); i++)
	  {
	    int r = rows[i], c = cols[i];
	    if (symmetric && (r < c)) swap (r, c);
	    graph.AddUnique (r, c);
	  }

	// *testout << "graphi = " << endl << graph << endl;

	Array<int> els_per_row(num_glob_dofs);
	for (int i = 0; i < num_glob_dofs; i++)
	  els_per_row[i] = graph[i].Size();

	cout << IM(4) << "now build matrix" << endl;

	// SparseMatrixSymmetric<TM> matrix(els_per_row);
        auto matrix = symmetric ? make_shared<SparseMatrixSymmetric<TM>> (els_per_row)
	  : make_shared<SparseMatrix<TM>> (els_per_row);

	for (int i = 0; i < rows.Size(); i++)
	  {
	    int r = rows[i], c = cols[i];
	    if (symmetric && (r < c) ) swap (r, c);
	    matrix->CreatePosition(r, c);
	  }
	matrix->AsVector() = 0.0;

	for (int i = 0; i < rows.Size(); i++)
	  {
	    int r = rows[i], c = cols[i];
	    if (symmetric && (r < c)) swap (r, c);
	    (*matrix)(r,c) += vals[i];
	  }

	cout << IM(4) << "have matrix, now invert" << endl;

	matrix->SetInverseType (mat.GetInverseType());
	inv = matrix->InverseMatrix ();
      }

    // MPI_Barrier (ngs_comm);
  }

  template <typename TM>
  MasterInverse<TM> :: ~MasterInverse ()
  {
    ; // delete inv;
  }

  template <typename TM>
  void MasterInverse<TM> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    typedef typename mat_traits<TM>::TV_ROW TV;
    
    NgMPI_Comm comm = paralleldofs->GetCommunicator();
    int id = comm.Rank();
    int ntasks = comm.Size();

    bool is_x_cum = (dynamic_cast_ParallelBaseVector(x) . Status() == CUMULATED);
    // x.Distribute();
    y.Cumulate();

    if (id > 0)
      {
	FlatVector<TV> fx = x.FV<TV> ();
	FlatVector<TV> fy = y.FV<TV> ();

	Array<TV> lx (select.Size());
	for (int i = 0; i < select.Size(); i++)
	  lx[i] = fx(select[i]);
	
	MPI_Request request = comm.ISend (lx, 0, MPI_TAG_SOLVE);
	MPI_Request_free (&request);

	// only for MUMPS:
	if (inv)
	  y = (*inv) * x;
	
	request = comm.IRecv (lx, 0, MPI_TAG_SOLVE);
	MPI_Wait (&request, MPI_STATUS_IGNORE);


	for (int i = 0; i < select.Size(); i++)
	  fy(select[i]) += s * lx[i];
      }

    else
      {
	VVector<TV> hx(inv->Height());
	VVector<TV> hy(inv->Height());
	hx = 0.0;

	Array<int> sizes(ntasks);
	for (int i = 0; i < ntasks; i++)
	  sizes[i] = loc2glob[i].Size();

	Table<TV> exdata(sizes);


	for (int src = 1; src < ntasks; src++)
	  {
	    FlatArray<int> selecti = loc2glob[src];

	    Array<TV> lx(selecti.Size());
	    comm.Recv (lx, src, MPI_TAG_SOLVE);

	    if(is_x_cum) {
	      for (int i = 0; i < selecti.Size(); i++)
		hx(selecti[i]) = lx[i];
	    }
	    else {
	      for (int i = 0; i < selecti.Size(); i++)
		hx(selecti[i]) += lx[i];
	    }
	  }

	hy = (*inv) * hx;

	Array<MPI_Request> requ;
	for (int src = 1; src < ntasks; src++)
	  {
	    FlatArray<int> selecti = loc2glob[src];
	    for (int j = 0; j < selecti.Size(); j++)
	      exdata[src][j] = hy(selecti[j]);
	    requ.Append (comm.ISend (exdata[src], src, MPI_TAG_SOLVE));
	  }
	MyMPI_WaitAll (requ);
      }

    // if (is_x_cum)
    // dynamic_cast_ParallelBaseVector(x) . Cumulate(); // AllReduce(&hoprocs);

  }




  template class MasterInverse<double>;
  template class MasterInverse<Complex>;

#if MAX_SYS_DIM >= 1
  template class MasterInverse<Mat<1,1,double> >;
  template class MasterInverse<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class MasterInverse<Mat<2,2,double> >;
  template class MasterInverse<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class MasterInverse<Mat<3,3,double> >;
  template class MasterInverse<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class MasterInverse<Mat<4,4,double> >;
  template class MasterInverse<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class MasterInverse<Mat<5,5,double> >;
  template class MasterInverse<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class MasterInverse<Mat<6,6,double> >;
  template class MasterInverse<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class MasterInverse<Mat<7,7,double> >;
  template class MasterInverse<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class MasterInverse<Mat<8,8,double> >;
  template class MasterInverse<Mat<8,8,Complex> >;
#endif





  
#endif

  
  ParallelMatrix :: ParallelMatrix (shared_ptr<BaseMatrix> amat,
				    shared_ptr<ParallelDofs> arpardofs,
				    shared_ptr<ParallelDofs> acpardofs,
				    PARALLEL_OP aop)
    : BaseMatrix((arpardofs==acpardofs) ? arpardofs : nullptr), mat(amat),
      row_paralleldofs(arpardofs), col_paralleldofs(acpardofs)
  {
    op = aop;
    if(row_paralleldofs==col_paralleldofs)
      mat->SetParallelDofs (arpardofs);
    if (auto spmat = dynamic_pointer_cast<BaseSparseMatrix>(mat)) {
#ifdef USE_MUMPS
      spmat->SetInverseType(MUMPS);
#else
      spmat->SetInverseType(MASTERINVERSE);
#endif
    }
  }

  ParallelMatrix :: ParallelMatrix (shared_ptr<BaseMatrix> amat,
				    shared_ptr<ParallelDofs> apardofs,
				    PARALLEL_OP aop)
    : ParallelMatrix(amat, apardofs, apardofs)
  { }

  
  AutoVector ParallelMatrix :: CreateRowVector () const
  {
    if (IsComplex()) {
      if (row_paralleldofs == nullptr)
	return make_unique<S_ParallelBaseVectorPtr<Complex>>
	  (mat->Width(), paralleldofs->GetEntrySize(), paralleldofs, DISTRIBUTED);
      else
	return make_unique<S_ParallelBaseVectorPtr<Complex>>
	  (mat->Width(), row_paralleldofs->GetEntrySize(), row_paralleldofs, DISTRIBUTED);
    }
    else {
      if (row_paralleldofs == nullptr)
	return make_unique<S_ParallelBaseVectorPtr<double>>
	  (mat->Width(), paralleldofs->GetEntrySize(), paralleldofs, DISTRIBUTED);
      else
	return make_unique<S_ParallelBaseVectorPtr<double>>
	  (mat->Width(), row_paralleldofs->GetEntrySize(), row_paralleldofs, DISTRIBUTED);
    }
  }
  
  AutoVector ParallelMatrix :: CreateColVector () const
  {
    if (IsComplex()) {
      if (col_paralleldofs==nullptr)
	return make_unique<S_ParallelBaseVectorPtr<Complex>>
	  (mat->Height(), paralleldofs->GetEntrySize(), paralleldofs, DISTRIBUTED);
      else
	return make_unique<S_ParallelBaseVectorPtr<Complex>>
	  (mat->Height(), col_paralleldofs->GetEntrySize(), col_paralleldofs, DISTRIBUTED);
    }
    else {
      if (col_paralleldofs==nullptr)
	return make_unique<S_ParallelBaseVectorPtr<double>>
	  (mat->Height(), paralleldofs->GetEntrySize(), paralleldofs, DISTRIBUTED);
      else
	return make_unique<S_ParallelBaseVectorPtr<double>>
	  (mat->Height(), col_paralleldofs->GetEntrySize(), col_paralleldofs, DISTRIBUTED);
    }
  }


  ParallelMatrix :: ~ParallelMatrix ()
  {
    ; // delete &mat;
  }

  void ParallelMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    const auto & xpar = dynamic_cast_ParallelBaseVector(x);
    auto & ypar = dynamic_cast_ParallelBaseVector(y);
    if (op & char(2))
      x.Cumulate();
    else
      x.Distribute();
    if (op & char(1))
      y.Cumulate();
    else
      y.Distribute();
    mat->MultAdd (s, *xpar.GetLocalVector(), *ypar.GetLocalVector());
  }

  void ParallelMatrix :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    const auto & xpar = dynamic_cast_ParallelBaseVector(x);
    auto & ypar = dynamic_cast_ParallelBaseVector(y);
    if (op & char(1))
      x.Distribute();
    else
      x.Cumulate();
    if (op & char(2))
      y.Distribute();
    else
      y.Cumulate();
    mat->MultTransAdd (s, *xpar.GetLocalVector(), *ypar.GetLocalVector());
  }

  shared_ptr<BaseMatrix> ParallelMatrix :: CreateMatrix () const
  {
    return make_shared<ParallelMatrix> (mat->CreateMatrix(), paralleldofs);
  }

  AutoVector ParallelMatrix :: CreateVector () const
  {
    if(row_paralleldofs != col_paralleldofs) {
      throw Exception("ParallelMatrix::CreateVector called for nonsymmetric case (use CreateRowVector of CreateColVector)!");
    }
    if (IsComplex()) {
      return make_unique<S_ParallelBaseVectorPtr<Complex>>
	(mat->Width(), paralleldofs->GetEntrySize(), paralleldofs, DISTRIBUTED);
    }
    else {
      return make_unique<S_ParallelBaseVectorPtr<double>>
	(mat->Width(), paralleldofs->GetEntrySize(), paralleldofs, DISTRIBUTED);
    }
  }

  ostream & ParallelMatrix :: Print (ostream & ost) const
  {
    ost << *mat;
    return ost;
  }

  int ParallelMatrix :: VHeight() const
  {
    return mat->VHeight();
  }

  int ParallelMatrix :: VWidth() const
  {
    return mat->VWidth();
  }


  shared_ptr<BaseMatrix> ParallelMatrix::InverseMatrix (shared_ptr<BitArray> subset) const
  {
    shared_ptr<BaseMatrix> inv;

    Iterate<MAX_SYS_DIM> ( [&](auto i) -> void {
	constexpr int N = 1+i.value;
	if (inv) return;
	if constexpr(N == 1) {
	    inv = InverseMatrixTM<double> (subset); if (inv) return;
	    inv = InverseMatrixTM<Complex> (subset); if (inv) return;
	  }
	else {
	  inv = InverseMatrixTM<Mat<N> > (subset); if (inv) return;
	  inv = InverseMatrixTM<Mat<N,N,Complex> > (subset); if (inv) return;
	}
      });

    if (inv) return inv;

    throw Exception ("ParallelMatrix::Inverse(BitArray) not available, typeid(mat) = " 
		     + ToString (typeid(mat).name()));
  }


  
  template <typename TM>
  shared_ptr<BaseMatrix> ParallelMatrix::InverseMatrixTM (shared_ptr<BitArray> subset) const
  {
    const SparseMatrixTM<TM> * dmat = dynamic_cast<const SparseMatrixTM<TM>*> (mat.get());
    if (!dmat) return nullptr;

#ifdef USE_MUMPS
    bool symmetric = dynamic_cast<const SparseMatrixSymmetric<TM>*> (mat.get()) != NULL;
    if (mat->GetInverseType() == MUMPS)
      return make_shared<ParallelMumpsInverse<TM>> (*dmat, subset, nullptr, paralleldofs, symmetric);
    else 
#endif

#ifdef PARALLEL
      return make_shared<MasterInverse<TM>> (*dmat, subset, paralleldofs);
#endif
    throw Exception ("ParallelMatrix: don't know how to invert");
  }



  shared_ptr<BaseMatrix> ParallelMatrix::InverseMatrix (shared_ptr<const Array<int>> clusters) const
  {
    throw Exception ("ParallelMatrix::Inverse(cluster) not available");
  }

  INVERSETYPE ParallelMatrix::SetInverseType (INVERSETYPE ainversetype) const
  {
    return mat->SetInverseType (ainversetype);
  }

  INVERSETYPE ParallelMatrix::SetInverseType (string ainversetype) const
  {
    return mat->SetInverseType (ainversetype);
  }
  
  INVERSETYPE ParallelMatrix::GetInverseType () const
  {
    return mat->GetInverseType ();
  }



#ifdef PARALLEL





  FETI_Jump_Matrix :: FETI_Jump_Matrix (shared_ptr<ParallelDofs> apardofs, shared_ptr<ParallelDofs> au_paralleldofs)
    : BaseMatrix(apardofs), u_paralleldofs(au_paralleldofs)
  {
    
    size_t njs = 0;
    for (auto p:paralleldofs->GetDistantProcs())
      njs += paralleldofs->GetExchangeDofs(p).Size();

    Array<size_t> ones(njs);
    ones = 1;
    Table<int> dps(ones);
    njs = 0;
    for (auto p : paralleldofs->GetDistantProcs()) 
      for ([[maybe_unused]] auto d : paralleldofs->GetExchangeDofs(p)) 
	dps[njs++][0] = p;

    this->jump_paralleldofs = make_shared<ParallelDofs>(paralleldofs->GetCommunicator(), move(dps));
  }


  void FETI_Jump_Matrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    y.Distribute();
    size_t count = 0;
    /*
      for (auto p:paralleldofs->GetDistantProcs()) {
      auto exdofs = paralleldofs->GetExchangeDofs(p);
      if (p<MyMPI_GetId(paralleldofs->GetCommunicator())) {
      for (auto k:Range(exdofs.Size())) {
      y.FVDouble()[count++] -= s*x.FVDouble()[exdofs[k]];
      }
      }
      else {
      for (auto k:Range(exdofs.Size())) {
      y.FVDouble()[count++] += s*x.FVDouble()[exdofs[k]];
      }
      }
      }
    */
    auto me = paralleldofs->GetCommunicator().Rank();
    auto fx = x.FVDouble();
    auto fy = y.FVDouble();
    for (auto p : paralleldofs->GetDistantProcs())
      {
        double signed_s = (p < me) ? -s : s;
        for (auto exdof : paralleldofs->GetExchangeDofs(p))
          fy[count++] += signed_s * fx[exdof];
      }
  }

  void FETI_Jump_Matrix :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    x.Cumulate();
    size_t count = 0;
    for (auto p:paralleldofs->GetDistantProcs()) {
      auto exdofs = paralleldofs->GetExchangeDofs(p);
      if (p<paralleldofs->GetCommunicator().Rank()) {
	for (auto k:Range(exdofs.Size())) {
	  y.FVDouble()[exdofs[k]] -= s*x.FVDouble()[count++];
	}
      }
      else {
	for (auto k:Range(exdofs.Size())) {
	  y.FVDouble()[exdofs[k]] += s*x.FVDouble()[count++];
	}
      }
    }
  }

  
  AutoVector FETI_Jump_Matrix :: CreateRowVector () const
  {
    // throw Exception("Called FETI_Jump_Matrix :: CreateRowVector, this is not well defined");
    if (u_paralleldofs==nullptr) {
      return make_unique<VVector<double>>(VHeight());
    }
    return make_unique<ParallelVVector<double>> (u_paralleldofs->GetNDofLocal(),
						 u_paralleldofs);
  }
  
  AutoVector FETI_Jump_Matrix :: CreateColVector () const
  {
    return make_unique<ParallelVVector<double>> (jump_paralleldofs->GetNDofLocal(),
						 jump_paralleldofs);
  }
  
#endif
  
}

