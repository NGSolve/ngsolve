#ifdef HYPRE
/*********************************************************************/
/* File:   hypre_ams_precond.cpp                                     */
/* Author: Lukas Kogler                                              */
/* Date:   May 2017                                                  */
/*********************************************************************/

#include <solve.hpp>

namespace ngcomp
{

  HypreAMSPreconditioner :: HypreAMSPreconditioner (const PDE & pde, const Flags & aflags, const string & aname)
    : Preconditioner(&pde, aflags)
  {
    throw Exception ("HYPRE-AMS pde-file-interface not implemented; please use direct python-interface");
  }
  
  HypreAMSPreconditioner :: HypreAMSPreconditioner (const BaseMatrix & matrix, const shared_ptr<BitArray> afreedofs)
    : Preconditioner(shared_ptr<BilinearForm>(nullptr), Flags("not_register_for_auto_update"))
  {
    throw Exception ("HYPRE-AMS basematrix-interface not implemented; please use direct python-interface");
  }
  
  HypreAMSPreconditioner :: HypreAMSPreconditioner (shared_ptr<BilinearForm> abfa, const Flags & aflags,
						    const string aname)
    : Preconditioner(abfa, aflags, aname)
  {
    throw Exception ("HYPRE-AMS blf-interface not implemented; please use direct python-interface");
  }

  HypreAMSPreconditioner :: ~HypreAMSPreconditioner ()
  {
    ;
  }
  
  void HypreAMSPreconditioner :: Update()
  {
    throw Exception ("HYPRE-AMS Update not implemented; please use direct python interface");
    return;
  }

  void HypreAMSPreconditioner :: FinalizeLevel (const BaseMatrix * mat)
  {
    throw Exception ("HYPRE-AMS Finalielevel not implemented; please use direct python interface");
    return;
  }

  void HypreAMSPreconditioner :: Setup (const BaseMatrix & matrix)
  {
    throw Exception ("HYPRE-AMS Setup not implemented; please use direct python interface");
    return;
  }
  
  HypreAMSPreconditioner :: HypreAMSPreconditioner (shared_ptr<FESpace> ahcurlfes, shared_ptr<BilinearForm> abfa,
						    shared_ptr<FESpace> ah1fes, shared_ptr<BaseMatrix> agrad_mat,
						    shared_ptr<BilinearForm> abf_alpha, shared_ptr<BilinearForm> abf_beta, 
						    shared_ptr<BaseVector> ozz, shared_ptr<BaseVector> zoz, 
						    shared_ptr<BaseVector> zzo, Flags & flags)
    : Preconditioner( abfa, flags, "hypreams"), bfa(abfa), ngs_grad_mat(agrad_mat), bf_alpha(abf_alpha), bf_beta(abf_beta),
      hcurlfes(ahcurlfes), h1fes(ah1fes)
  {
    cout << IM(1) << "Setup AMS-Hypre preconditioner" << endl;
    static Timer t("AMS-hypre setup");
    RegionTimer reg(t);

    HYPRE_Int err;	

    this->rank = MyMPI_GetId();
    this->np = MyMPI_GetNTasks();

    this->parallel = (np>1);
    
    auto matrix = bfa->GetMatrixPtr();
    const ParallelMatrix* pmat = dynamic_cast<const ParallelMatrix*> (matrix.get());
    // const SparseMatrix<double>* mat_spm = dynamic_cast<SparseMatrix<double>*>(&pmat->GetMatrix());
    // const SparseMatrixTM<double>* mat_tm = dynamic_cast<SparseMatrix<double>*>(&pmat->GetMatrix());
    const SparseMatrix<double> * some_ptr = dynamic_cast<const SparseMatrix<double>*>(matrix.get());
    const SparseMatrix<double> & mat = (pmat!=nullptr)?dynamic_cast< const SparseMatrix<double> &>(pmat->GetMatrix()):dynamic_cast< const SparseMatrix<double> &>(*matrix.get());
    if (dynamic_cast< const SparseMatrixSymmetric<double> *> (&mat))
      throw Exception ("Please use fully stored sparse matrix for hypre (bf -nonsymmetric)");

    cout << "hi" << endl;

    auto compute_global_nums = [this](auto ndof, auto & fes, auto & pardofs, auto & global_nums,
				      auto & ilower, auto & iupper) 
      {
	global_nums.SetSize(ndof);
	if(ndof) global_nums = -1;
	int num_master_dofs = 0;
	for (int i = 0; i < ndof; i++)
	  if (pardofs -> IsMasterDof(i))
	    global_nums[i] = num_master_dofs++;
	Array<int> first_master_dof(this->np);
	MPI_Allgather (&num_master_dofs, 1, MPI_INT, 
		       &first_master_dof[0], 1, MPI_INT, 
		       pardofs -> GetCommunicator());    
	int num_glob_dofs = 0;
	for (int i = 0; i < this->np; i++)
	  {
	    int cur = first_master_dof[i];
	    first_master_dof[i] = num_glob_dofs;
	    num_glob_dofs += cur;
	  }
	first_master_dof.Append(num_glob_dofs);
	for (int i = 0; i < ndof; i++)
	  if (global_nums[i] != -1)
	    global_nums[i] += first_master_dof[this->rank];    
	ScatterDofData (global_nums, *pardofs);
	ilower = first_master_dof[this->rank];
	iupper = first_master_dof[this->rank+1]-1;
      };
    
    cout << "hi" << endl;

    this->hc_freedofs = this->hcurlfes->GetFreeDofs();
    this->hc_ndof = this->hc_freedofs->Size();
    this->h1_freedofs = this->h1fes->GetFreeDofs();
    this->h1_ndof = this->h1_freedofs->Size();
    
    if(parallel) {
      this->hc_pardofs = shared_ptr<ParallelDofs>(const_cast<ParallelDofs*>(abfa->GetMatrix().GetParallelDofs()), NOOP_Deleter);
      (void) compute_global_nums(this->hc_ndof, this->hcurlfes, this->hc_pardofs, this->hc_global_nums, this->hc_ilower, this->hc_iupper);
      this->h1_pardofs = shared_ptr<ParallelDofs>(const_cast<ParallelDofs*>(abf_alpha->GetMatrix().GetParallelDofs()), NOOP_Deleter);
      (void) compute_global_nums(this->h1_ndof, this->h1fes, this->h1_pardofs, this->h1_global_nums, this->h1_ilower, this->h1_iupper);
      cout << "global hcurl h1 dofs: " << this->hc_pardofs->GetNDofGlobal() << " " << this->h1_pardofs->GetNDofGlobal() << endl;
      this->hc_masterdofs.SetSize(this->hc_iupper-this->hc_ilower+1);
      int s = 0;
      for(auto k:Range(this->hc_ndof))
	if(this->hc_pardofs->IsMasterDof(k))
	  this->hc_masterdofs[s++] = k;
      this->buf_hc.SetSize(this->hc_masterdofs.Size());
      this->hc_intrange.SetSize(this->hc_iupper-this->hc_ilower+1);
      for(auto k:Range(this->hc_intrange.Size()))
	this->hc_intrange[k] = this->hc_ilower+k;
    }
    else {
      this->hc_global_nums.SetSize(hc_ndof);
      for(auto k:Range(this->hc_global_nums.Size()))
	this->hc_global_nums[k] = k;
      this->hc_ilower = 0;
      this->hc_iupper = this->hc_ndof-1;
      this->h1_global_nums.SetSize(h1_ndof);
      for(auto k:Range(this->h1_global_nums.Size()))
	this->h1_global_nums[k] = k;
      this->h1_ilower = 0;
      this->h1_iupper = this->h1_ndof-1;
      this->hc_masterdofs.SetSize(this->hc_ndof);
      for(auto k:Range(hc_masterdofs.Size()))
	hc_masterdofs[k] = k;
      this->buf_hc.SetSize(this->hc_ndof);
      this->hc_intrange.SetSize(this->hc_ndof);
      for(auto k:Range(this->hc_intrange.Size()))
	this->hc_intrange[k] = this->hc_ilower+k;
    }

    cout << "hi" << endl;

    cout << "h1 fd: " << h1_freedofs->NumSet() << "/" << h1_freedofs->Size() << endl << "hc fd: " << hc_freedofs->NumSet() << "/" << hc_freedofs->Size() << endl;
    
    cout << "nes: " << hcurlfes->GetMeshAccess()->GetNEdges() << endl;
    
    // create AMS solver
    HYPRE_AMSCreate(&this->precond);
    HYPRE_AMSSetPrintLevel(this->precond, 1);

    // set space dimension
    this->dimension = flags.GetNumFlag("dimension", 3);

    cout << "dimension: " << this->dimension << endl;

    HYPRE_AMSSetDimension(this->precond, this->dimension);
    /*
    auto check_IJ_mat_par = [](auto & ijmat, auto & pijmat, auto & ngsmat, auto & globnums, auto & globnums_c,
			       auto & pd) {      
      Array<int> rows(ngsmat->Height());
      Array<int> nc(ngsmat->Height());
      Array<int> grows(ngsmat->Height());
      for(auto k:Range(rows.Size())) {
	rows[k] = k;
	grows[k] = globnums[k];
      }
      HYPRE_IJMatrixGetRowCounts(ijmat, rows.Size(), &grows[0], &nc[0]);
      for(auto k:Range(rows.Size())) {
	if(!pd->IsMasterDof(k)) continue;
	auto cols = ngsmat->GetRowIndices(k);
	auto vals = ngsmat->GetRowValues(k);
	// Array<double> hvals(nc[k]);
	// HYPRE_IJMatrixGetValues(ijmat, 1, &nc[k], &rows[k], &cols[0], &hvals[0]);
	double* hvals = NULL;
	int* hcols = NULL;
	int hnc = nc[k];
	int err = 0;
	if( (err = HYPRE_ParCSRMatrixGetRow(pijmat, globnums[k], &hnc, &hcols, &hvals)) != 0)
	  cerr << "hypre returned with error from HYPRE_ParCSRMatrixGetRow!!" << endl;
		bool same = (hnc==cols.Size())?true:false;
	if(same)
	  for(auto j:Range(cols.Size()))
	    if(cols[j]!=hcols[j])
	      same = false;	
	if(!same) {
	  cout << endl;
	  cout << "row " << k << " -> " << globnums[k] << ", nc " << cols.Size() << ", hnc " << hnc << endl;
	  cout << " cols:  ";
	  for(auto c:cols)
	    cout << c << " ";
	  cout << endl;
	  cout << "vals:  ";
	  for(auto v:vals)
	    cout << v << " ";
	  cout << endl;
	  cout << "hcols: ";
	  for(auto l:Range(hnc))
	    cout << globnums_c.Pos(hcols[l]) << " ";
	  cout << endl;
	  cout << "hvals: ";
	  for(auto l:Range(hnc))
	    cout << hvals[l] << " ";
	  cout << endl;
	}
      }
    };      
    */    
    auto check_IJ_mat = [](auto & ijmat, auto & ngsmat) {      
      return; //below: sequential/full check
      cout << "check IJmat!" << endl;
      Array<int> rows(ngsmat->Height());
      Array<int> nc(ngsmat->Height());
      for(auto k:Range(rows.Size()))
	rows[k] = k;
      HYPRE_IJMatrixGetRowCounts(ijmat, rows.Size(), &rows[0], &nc[0]);
      for(auto k:Range(rows.Size()))
	if(nc[k]!=ngsmat->GetRowIndices(k).Size())
	  cout << "diff @ row " << k << ": " << nc[k] << " " << ngsmat->GetRowIndices(k).Size() << endl;
      for(auto k:Range(rows.Size())) {
	auto cols = ngsmat->GetRowIndices(k);
	auto vals = ngsmat->GetRowValues(k);
	Array<double> hvals(cols.Size());
	HYPRE_IJMatrixGetValues(ijmat, 1, &nc[k], &rows[k], &cols[0], &hvals[0]);
	for(auto j:Range(hvals.Size()))
	  hvals[j] -= vals[j];
	bool have_err = false;
	for(auto j:Range(hvals.Size()))
	  if(hvals[j]!=0)
	    have_err = true;
	if(have_err) {
	  cout << "row " << k << ": " << endl;
	  for(auto j:Range(cols.Size()))
	    cout << cols[j] << " ";
	  cout << endl;
	  for(auto j:Range(vals.Size()))
	    cout << vals[j] << " ";
	  cout << endl;
	  for(auto j:Range(vals.Size()))
	    cout << hvals[j] << " ";
	  cout << endl;
	}
      }      
      cout << "IJMAT checked!!" << endl;
      return;
    };

    //     create_IJ_from_blf(&this->A, &this->parcsr_A, this->bfa,
    // 		       this->hc_global_nums, this->hc_ilower, this->hc_iupper);
    // construct hypre-IJ-matrix from ngsolve matrix; dirichlet-rows have a "1" in the diagonal
    auto create_IJ_from_spmat = [this, &err, &check_IJ_mat] 
      (HYPRE_IJMatrix* ijmat, HYPRE_ParCSRMatrix* pijmat, 
       const SparseMatrix<double> & ngsmat, shared_ptr<ParallelDofs> pardofs,
       shared_ptr<BitArray> freedofs, Array<int> & global_nums,
       int &ilower, int & iupper)
      {
	int ndof = ngsmat.Height();
	
	cout << "ij from spmat, have glob nums?" << global_nums.Size() << endl;
	
	HYPRE_IJMatrixCreate(ngs_comm, ilower, iupper, ilower, iupper, ijmat);
	HYPRE_IJMatrixSetPrintLevel(*ijmat, 1);
	HYPRE_IJMatrixSetObjectType(*ijmat, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(*ijmat);

	int maxperow = 0;
	for(auto k:Range(ndof))
	  maxperow = max2((int)maxperow, (int)ngsmat.GetRowIndices(k).Size());
	Array<int> cols_global(maxperow);
	Array<double> vals_global(maxperow);
	for(int k:Range(ndof)) {
	  int row = global_nums[k];
	  if(row==-1) cout << "WARNING, globnum==-1" << endl;
	  int s = 0;
	  if(freedofs->Test(k)) {
	    auto cols = ngsmat.GetRowIndices(k);
	    auto vals = ngsmat.GetRowValues(k);
	    for(auto j:Range(cols.Size())) {
	    if(global_nums[cols[j]]==-1) cout << "Warning; GCOLNUM==-1" << endl;
	    cols_global[s] = global_nums[cols[j]];
	    vals_global[s++] = vals[j];
	    }
	  }
	  else {
	    s = 1;
	    cols_global[0] = global_nums[k];
	    vals_global[0] = 1.0;
	  }
	  // cout << "add to row " << row << " (" << s << " vals): " << endl;
	  // for(auto k:Range(s))
	  //   cout << k << " " << cols_global[k] << " " << vals_global[k] << endl;
	  // cout << endl;
	  HYPRE_IJMatrixAddToValues(*ijmat, 1, &s, &row, &cols_global[0], &vals_global[0]);
	}

	if( (err = HYPRE_IJMatrixAssemble(*ijmat)) != 0 )
	  cerr << "HYPRE_IJMatrixAssemble returned with error: " << err << endl;

	if( (err = HYPRE_IJMatrixGetObject(*ijmat, (void**) pijmat)) !=0)
	  cerr << "HYPRE_IJMatrixGetObject returned with error: " << err << endl;      
	
	cout << "mat done, check it.." << endl;
	const auto ptr = &ngsmat;
	(void) check_IJ_mat(*ijmat, ptr);	
	cout << "mat check done!" << endl;
	return;
      };

    
    auto create_IJ_from_blf = [&create_IJ_from_spmat](HYPRE_IJMatrix* ijmat, HYPRE_ParCSRMatrix* pijmat, 
						      shared_ptr<BilinearForm> & bf, Array<int> & global_nums,
						      int &ilower, int & iupper)
      {
	cout << "ij from blf" << endl;
	auto fes = bf->GetFESpace();
	auto fd = fes->GetFreeDofs();
	auto matrix = bf->GetMatrixPtr();
	auto pmat = dynamic_cast<const ParallelMatrix*>(matrix.get());
	shared_ptr<ParallelDofs> pardofs = (pmat!=nullptr)?shared_ptr<ParallelDofs>(const_cast<ParallelDofs*>(pmat->GetParallelDofs()), NOOP_Deleter):nullptr;
	auto spmat = dynamic_cast<const SparseMatrix<double>*>( ( (pmat!=nullptr)?&pmat->GetMatrix():matrix.get() ) );
	cout << "fd: " << fd->NumSet() << " out of " << fd->Size() << endl;
	if (dynamic_cast< const SparseMatrixSymmetric<double> *> (spmat))
	  throw Exception ("Please use fully stored sparse matrix for hypre (bf -nonsymmetric)");
	(void) create_IJ_from_spmat(ijmat, pijmat, *spmat, pardofs, fd, global_nums, ilower, iupper);
	return;
      };

    // cout << "alpha mat" << endl;
    create_IJ_from_blf(&this->alpha_mat, &this->parcsr_alpha_mat, bf_alpha, 
    		       this->h1_global_nums, this->h1_ilower, this->h1_iupper);    
    // cout << "alpha done" << endl;

    // set gradient matrix
    {
      cout << "h1: " << this->h1fes->GetNDof() << endl;
      cout << "hcurl: " << this->hcurlfes->GetNDof() << endl;
      cout << "local grad mat: " << ngs_grad_mat->Height() << " x " << ngs_grad_mat->Width() << endl;
      cout << "now grad mat.." << endl;
      HYPRE_IJMatrixCreate(ngs_comm, this->hc_ilower, this->hc_iupper, this->h1_ilower, this->h1_iupper, &this->grad_mat);
      HYPRE_IJMatrixSetPrintLevel(this->grad_mat, 1);
      HYPRE_IJMatrixSetObjectType(this->grad_mat, HYPRE_PARCSR);
      HYPRE_IJMatrixInitialize(this->grad_mat);
      auto spgrad = dynamic_cast<SparseMatrix<double>*>(this->ngs_grad_mat.get());
      int rmax = 0;
      for(int k:Range(this->ngs_grad_mat->Height()))
	rmax = max2(rmax, (int)spgrad->GetRowIndices(k).Size());
      Array<int> cols_global(rmax);
      Array<double> vals_global(rmax);
      int s = 0;
      for(auto k:Range(this->ngs_grad_mat->Height())) {
	auto cols = spgrad->GetRowIndices(k);
	auto vals = spgrad->GetRowValues(k);
	s = 0;
	for(auto j:Range(cols.Size())) {
	  if(parallel && (!this->hc_pardofs->IsMasterDof(k)) && (!this->h1_pardofs->IsMasterDof(cols[j]))) continue;
	  cols_global[s] = this->h1_global_nums[cols[j]];
	  vals_global[s++] = vals[j];
	}
	int row = this->hc_global_nums[k];
	if(s) HYPRE_IJMatrixAddToValues(this->grad_mat, 1, &s, &row, &cols_global[0], &vals[0]);
      }
      HYPRE_IJMatrixAssemble(this->grad_mat);
      HYPRE_IJMatrixGetObject(this->grad_mat, (void**) &this->parcsr_grad_mat);
      if( (err = HYPRE_AMSSetDiscreteGradient(this->precond, this->parcsr_grad_mat)) != 0)
	cerr << "HYPRE_AMSSetDiscreteGradient returned with err " << err << endl;
      cout << "grad mat done!" << endl;
      cout << "check gradmat!" << endl;
      const SparseMatrix<double>* Gmat = dynamic_cast<SparseMatrix<double>*>(ngs_grad_mat.get());
      (void)check_IJ_mat(this->grad_mat, Gmat);
      // if(this->parallel) (void) check_IJ_mat_par(this->grad_mat, this->parcsr_grad_mat, Gmat, this->hc_global_nums, this->h1_global_nums, thjis->hc_pardofs);
      cout << "grad mat done!" << endl;
      }
      
     // set coords for constant vecs
    auto create_hijvec  = [this](HYPRE_IJVector & v, HYPRE_ParVector & pv, shared_ptr<BaseVector> ngs_vec, 
				 Array<int> & global_nums, int ilower, int iupper)
      {
	if(parallel)
	  ngs_vec->Distribute();
	HYPRE_IJVectorCreate(ngs_comm, ilower, iupper,&v);
	HYPRE_IJVectorSetPrintLevel(v, 1);
	HYPRE_IJVectorSetObjectType(v, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(v);
	Array<double> zeros(global_nums.Size());
	if(zeros.Size()) {
	  zeros = 0.0;
	  HYPRE_IJVectorSetValues(v, zeros.Size(), &global_nums[0], &zeros[0]);
	  HYPRE_IJVectorAddToValues(v, zeros.Size(), &global_nums[0], &ngs_vec->FVDouble()[0]);
	  // HYPRE_IJVectorSetValues(v, global_nums.Size(), &global_nums[0], &ngs_vec->FVDouble()[0]);
	}
	HYPRE_IJVectorAssemble(v);
	HYPRE_IJVectorGetObject(v, (void **) &pv);
	return;
      };

    // x, y and z vectors
    VVector<double> x(h1_ndof);
    BaseVector & bx(x);
    VVector<double> y(h1_ndof);
    BaseVector & by(y);
    VVector<double> z(h1_ndof);
    BaseVector & bz(z);
    x.FVDouble() = 0.0;
    y.FVDouble() = 0.0;
    z.FVDouble() = 0.0;
    auto ma = hcurlfes->GetMeshAccess();
    for(auto k:Range(h1_ndof)) {
      if( parallel && (!this->h1_pardofs->IsMasterDof(k)) ) continue;
      Vec<3> co;
      ma->GetPoint(k, co);
      x.FVDouble()[k] = co[0];
      y.FVDouble()[k] = co[1];
      z.FVDouble()[k] = co[2];
    }
    HYPRE_IJVector hx;
    HYPRE_ParVector phx;
    (void) create_hijvec(hx, phx, shared_ptr<BaseVector>(&bx, NOOP_Deleter), this->h1_global_nums, this->h1_ilower, this->h1_iupper);
    HYPRE_IJVector hy;
    HYPRE_ParVector phy;
    (void) create_hijvec(hy, phy, shared_ptr<BaseVector>(&by, NOOP_Deleter), this->h1_global_nums, this->h1_ilower, this->h1_iupper);
    HYPRE_IJVector hz;
    HYPRE_ParVector phz;
    (void) create_hijvec(hz, phz, shared_ptr<BaseVector>(&bz, NOOP_Deleter), this->h1_global_nums, this->h1_ilower, this->h1_iupper);
    if( (err = HYPRE_AMSSetCoordinateVectors(this->precond, phx, phy, phz)) != 0)
      cerr << "HYPRE_AMSSetCoordinateVectors returned with error " << err << endl;

    // (1,0,0), (0,1,0) and (0,0,1)
    // HYPRE_IJVector hozz;
    // HYPRE_ParVector phozz;
    // (void) create_hijvec(hozz, phozz, ozz, this->hc_global_nums, this->hc_ilower, this->hc_iupper);
    // HYPRE_IJVector hzoz;
    // HYPRE_ParVector phzoz;
    // (void) create_hijvec(hzoz, phzoz, zoz, this->hc_global_nums, this->hc_ilower, this->hc_iupper);
    // HYPRE_IJVector hzzo;
    // HYPRE_ParVector phzzo;
    // (void) create_hijvec(hzzo, phzzo, zzo, this->hc_global_nums, this->hc_ilower, this->hc_iupper);
    // if( (err = HYPRE_AMSSetEdgeConstantVectors(this->precond, phozz, phzoz, phzzo)) != 0 )
    // cerr << "HYPRE_AMSSetEdgeConstantVectors returned with error " << err << endl;

    // set beta poisson mat (or tell AMS if beta is zero)
    this->beta_is_zero = flags.GetDefineFlag("beta_is_zero");
    cout << "beta = zero? " << this->beta_is_zero << endl;
    // if(this->beta_is_zero)
    //   HYPRE_AMSSetBetaPoissonMatrix(this->precond, NULL);
    // else {
    //   cout << "beta" << endl;
    create_IJ_from_blf(&this->beta_mat, &this->parcsr_beta_mat, bf_beta,
    		       this->h1_global_nums, this->h1_ilower, this->h1_iupper);
    //   cout << "done, set beta" << endl;
    // HYPRE_AMSSetBetaPoissonMatrix(this->precond, this->parcsr_beta_mat);
    //   cout << "done!" << endl;
    // }

    // set alpha poisson mat
    cout << "set alpha..." << endl;
    // HYPRE_AMSSetAlphaPoissonMatrix(this->precond, this->parcsr_alpha_mat);
    cout << "done!" << endl;
    

    // cout << "grad mat: " << endl << *ngs_grad_mat << endl;
    // cout << "a mat: " << endl << bfa->GetMatrix() << endl;
    
    // set up solver
    HYPRE_IJVector hb;
    HYPRE_ParVector par_b;
    HYPRE_IJVector hv;
    HYPRE_ParVector par_v;
    (void) create_hijvec(hb, par_b, zzo, this->hc_global_nums, this->hc_ilower, this->hc_iupper);
    (void) create_hijvec(hv, par_v, zzo, this->hc_global_nums, this->hc_ilower, this->hc_iupper);

    cout << "A-mat" << endl;
    create_IJ_from_blf(&this->A, &this->parcsr_A, this->bfa,
		       this->hc_global_nums, this->hc_ilower, this->hc_iupper);
    cout << "check amat!" << endl;
    const SparseMatrix<double>* Aspm = dynamic_cast<const SparseMatrix<double>*>(&(this->bfa->GetMatrix()));
    (void)check_IJ_mat(this->A, Aspm);
    cout << "check done!" << endl;
    cout << "setup ams solver..." << endl;
    HYPRE_AMSSetup(this->precond, this->parcsr_A, par_b, par_v);
    cout << "done!" << endl;

    // create working vecs...
    // (void) create_hijvec(this->b, this->par_b, zzo, this->hc_global_nums, this->hc_ilower, this->hc_iupper);
    // (void) create_hijvec(this->x, this->par_x, zzo, this->hc_global_nums, this->hc_ilower, this->hc_iupper);

    // ams options
    HYPRE_AMSSetTol(this->precond, 0);
    HYPRE_AMSSetMaxIter(this->precond,1);
    HYPRE_AMSSetCycleType(this->precond,1);

    HYPRE_AMSSetSmoothingOptions(this->precond, 4, 1, 1.0, 1.0);      //truncated l1-smoother
    HYPRE_AMSSetBetaAMGOptions(this->precond, 10, 10, 6, 0.5, 0, 0);  //symmetric l1-smoother
    HYPRE_AMSSetAlphaAMGOptions(this->precond, 10, 10, 6, 0.5, 0, 0); //symmetric l1-smoother


    cout << endl;
    cout << endl;
    cout << "hc: " << hc_ndof << " " << hc_ilower << " " << hc_iupper << " " << hc_iupper-hc_ilower+1 << endl;
    cout << endl;
    cout << "globnums: " << endl << this->hc_global_nums << endl;
    cout << endl;
    cout << "hc_intrange: " << endl << hc_intrange << endl;
    cout << endl;
    cout << "hc_masterdofs: " << endl << hc_masterdofs << endl;
    cout << endl;

    cout << endl;
    cout << endl;
    cout << "h1: " << h1_ndof << " " << h1_ilower << " " << h1_iupper << " " << h1_iupper-h1_ilower+1 << endl;
    cout << endl;
    cout << "globnums: " << endl << this->h1_global_nums << endl;
    cout << endl;

    HYPRE_IJMatrixPrint(A, "IJ.out.A");
    HYPRE_IJMatrixPrint(grad_mat, "IJ.out.grad_mat");
    // HYPRE_IJMatrixPrint(alpha_mat, "IJ.out.alpha");
    // HYPRE_IJMatrixPrint(beta_mat, "IJ.out.beta");

    
  }

  void HypreAMSPreconditioner :: Mult (const BaseVector & f, BaseVector & u) const
  {
    static Timer t("AMS-hypre mult");
    RegionTimer reg(t);
    
    if(parallel) {
      f.Cumulate();
      u.SetParallelStatus(DISTRIBUTED);
    }
    
    /**
       note: apparently we need to re-assemble vectors in every iteration...
       no idea why it does not work when reusing old ones...
    **/
    HYPRE_IJVector b;
    HYPRE_ParVector par_b;
    HYPRE_IJVector x;
    HYPRE_ParVector par_x;

    HYPRE_IJVectorCreate(ngs_comm, this->hc_ilower, this->hc_iupper,&b);
    HYPRE_IJVectorSetPrintLevel(b, 1);
    HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(b);

    HYPRE_IJVectorCreate(ngs_comm, this->hc_ilower, this->hc_iupper,&x);
    HYPRE_IJVectorSetPrintLevel(x, 1);
    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(x);

    Array<double> zeros(this->hc_intrange.Size());
    zeros = 0.0;
    
    if(this->hc_ndof) {
      for(auto k:Range(this->buf_hc.Size())) {
	buf_hc[k] = f.FVDouble()[this->hc_masterdofs[k]];
      }
      HYPRE_IJVectorSetValues(b, this->hc_intrange.Size(), &this->hc_intrange[0], &this->buf_hc[0]);
    }
    //HYPRE_IJVectorSetValues(b, this->hc_global_nums.Size(), &this->hc_global_nums[0], &f.FVDouble()[0]);
    HYPRE_IJVectorAssemble(b);
    HYPRE_IJVectorGetObject(b, (void**) &par_b);
    
    HYPRE_IJVectorSetValues(x, this->hc_intrange.Size(), &this->hc_intrange[0], &zeros[0]);
    HYPRE_IJVectorAssemble(x);
    HYPRE_IJVectorGetObject(x, (void**) &par_x);

    HYPRE_AMSSolve(this->precond, this->parcsr_A, par_b, par_x);
    
    cout << "get AMS sol!" << endl;
    HYPRE_IJVectorGetValues(x, this->hc_intrange.Size(), &this->hc_intrange[0], &this->buf_hc[0]);
    
    if(this->hc_ndof) u.FVDouble() = 0.0;    
    for(auto k:Range(this->buf_hc.Size()))
      u.FVDouble()[this->hc_masterdofs[k]] = this->buf_hc[k];

    cout << "have it!" << endl;
    // HYPRE_IJVectorGetValues(x, this->hc_ndof, &this->hc_global_nums[0], &zeros[0]);
    // for(auto k:Range(hc_ndof))
    //  u.FVDouble()[k] = zeros[k];

    if(parallel)
      u.Cumulate();

    return;
  }


  // dont need to register, those constructors are not implemented anyways...
  // static RegisterPreconditioner<HypreAMSPreconditioner> init_hyprepre ("AMShypre");
}
#endif
