
#ifdef HYPRE
/*********************************************************************/
/* File:   hypre_ams_precond.cpp                                     */
/* Author: Lukas Kogler                                              */
/* Date:   May 2017                                                  */
/*********************************************************************/

#include <solve.hpp>
#include "hypre_ams_precond.hpp"

namespace ngcomp
{

  /**
     Create hypre-ij-mat from ngsolve sparse-matrix.
     If matrix_cumulated==true, assumes that each subprocess has "full" values for a_ij
       (ex: discrete gradient matrix)
     Else, if matrix_cumulated==false, assumes that all have "partial" values
       (ex: locally assembled FEM matrix)
  **/
  void Create_IJMat_from_SPMat (HYPRE_IJMatrix* ijmat, HYPRE_ParCSRMatrix* pijmat, 
				const SparseMatrix<double> & ngsmat, 
				const shared_ptr<ParallelDofs> & row_pardofs, const shared_ptr<ParallelDofs> & col_pardofs,
				const shared_ptr<BitArray> & row_freedofs, const shared_ptr<BitArray> & col_freedofs,
				const int & row_ilower, const int & row_iupper, const int & col_ilower, const int & col_iupper,
				FlatArray<int> row_gnums, FlatArray<int> col_gnums,
				bool matrix_cumulated = false)
  {
    NgMPI_Comm comm = row_pardofs->GetCommunicator();
    HYPRE_IJMatrixCreate(comm, row_ilower, row_iupper, col_ilower, col_iupper, ijmat);
    HYPRE_IJMatrixSetPrintLevel(*ijmat, 1);
    HYPRE_IJMatrixSetObjectType(*ijmat, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(*ijmat);

    size_t nr = ngsmat.Height();
    size_t nc = ngsmat.Width();
    
    if( ( (row_pardofs) && (nr != row_pardofs->GetNDofLocal()) ) || (nr != row_gnums.Size()) ) {
      cerr << "# of rows in Create_IJ_from_spmat not consistent!" 
	   << nr << " " << row_pardofs->GetNDofLocal() << " " << row_gnums.Size() << endl;
      return;
    }
    if( ( (col_pardofs) && (nc != col_pardofs->GetNDofLocal()) ) || (nc != col_gnums.Size()) ) {
      cerr << "# of cols in Create_IJ_from_spmat not consistent!" << endl
	   << nc << " " << col_pardofs->GetNDofLocal() << " " << col_gnums.Size() << endl;
      return;
    }
    
    size_t maxperow = 0;
    for(auto k:Range(nr))
      maxperow = max2(maxperow, ngsmat.GetRowIndices(k).Size());
    
    Array<int> cg(maxperow);
    Array<double> vg(maxperow);
    int s = 0;
    for(auto k:Range(nr)) {
      auto cols = ngsmat.GetRowIndices(k);
      auto vals = ngsmat.GetRowValues(k);
      int row = row_gnums[k];
      s = 0;
      if( (!row_freedofs) || row_freedofs->Test(k) ) {
	for(auto j:Range(cols.Size())) {
	  if( (row_pardofs) && (matrix_cumulated) && (!row_pardofs->IsMasterDof(k)) && (!col_pardofs->IsMasterDof(cols[j])) )
	    continue; //if someone else is master of both dofs, they have to have the fill val!
	  if( (col_freedofs) && (!col_freedofs->Test(cols[j])) )
	    continue;
	  cg[s] = col_gnums[cols[j]];
	  vg[s++] = vals[j];
	}
      }
      else { //dirichlet DOF - only a 1 @ diag
	if( (row_pardofs) && (!row_pardofs->IsMasterDof(k)) )
	  continue;
	s = 1;
	cg[0] = row;
	vg[0] = 1.0;
      }
    
      if(!s) continue; //no entries to add to matrix!
      HYPRE_IJMatrixAddToValues(*ijmat, 1, &s, &row, cg.Data(), vg.Data());
    }    
    HYPRE_IJMatrixAssemble(*ijmat);
    HYPRE_IJMatrixGetObject(*ijmat, (void**) pijmat);

    return;
  };

  /**
     Creates hypre-IJvector from ngsolve-vector and copies values.
   **/
  void Create_IJVec_from_BVec(NgMPI_Comm & comm, HYPRE_IJVector & v, HYPRE_ParVector & pv, double* vals,
			      Array<int> & global_nums, int ilower, int iupper, bool full_vals) {    
    HYPRE_IJVectorCreate(comm, ilower, iupper,&v);
    HYPRE_IJVectorSetPrintLevel(v, 1);
    HYPRE_IJVectorSetObjectType(v, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(v);
    Array<double> zeros(global_nums.Size());
    if(zeros.Size()) {
      zeros = 0.0;
      if(full_vals) //cumul. or not parallel
	HYPRE_IJVectorSetValues(v, global_nums.Size(), global_nums.Data(), vals);
      else
	HYPRE_IJVectorAddToValues(v, zeros.Size(), global_nums.Data(), vals);
    }
    HYPRE_IJVectorAssemble(v);
    HYPRE_IJVectorGetObject(v, (void **) &pv);
    return;
  }
    


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
    : Preconditioner(abfa, aflags, aname), bfa(abfa), bf_beta(nullptr), bf_alpha(nullptr),
      x(NULL), par_x(NULL), b(NULL), par_b(NULL), dimension(aflags.GetNumFlag("dimension", 3)),
      beta_is_zero(flags.GetDefineFlag("beta_is_zero"))
  {
    /** Get HC-FESpace **/
    shared_ptr<HCurlHighOrderFESpace> hcfes = dynamic_pointer_cast<HCurlHighOrderFESpace>(abfa->GetFESpace());
    if(hcfes==nullptr)
      throw Exception ("HYPRE-AMS Setup needs HCurlHighOrder-FESpace");

    /** Get H1-FESpace + discrete gradient matrix **/
    this->hcurlfes = hcfes;
    auto h1f1 = hcfes->CreateGradientSpace();
    this->ngs_grad_mat = hcfes->CreateGradient(*h1f1);
    this->h1fes = h1f1;
    
    return;
  }

  HypreAMSPreconditioner :: ~HypreAMSPreconditioner ()
  {
    ;
  }
  
  void HypreAMSPreconditioner :: Update()
  {
    /** Call setup! **/
    this->Setup();
    return;
  }

  void HypreAMSPreconditioner :: FinalizeLevel (const BaseMatrix * mat)
  {
    /** Call setup! **/
    this->Setup();
    return;
  }

  HypreAMSPreconditioner :: HypreAMSPreconditioner (shared_ptr<FESpace> ahcurlfes, shared_ptr<BilinearForm> abfa,
						    shared_ptr<FESpace> ah1fes, shared_ptr<BaseMatrix> agrad_mat,
						    shared_ptr<BilinearForm> abf_alpha, shared_ptr<BilinearForm> abf_beta, 
						    Flags & flags)
    : Preconditioner( abfa, flags, "hypreams"), bfa(abfa), ngs_grad_mat(agrad_mat), bf_alpha(abf_alpha), bf_beta(abf_beta),
      hcurlfes(ahcurlfes), h1fes(ah1fes), dimension(flags.GetNumFlag("dimension", 3)), beta_is_zero(flags.GetDefineFlag("beta_is_zero"))
  {
    this-> Setup();
  }


  void HypreAMSPreconditioner :: Setup ()
  {

    cout << IM(1) << "Setup Hypre-AMS preconditioner" << endl;

    static Timer t("Hypre-AMS setup");
    RegionTimer reg(t);

    HYPRE_Int err;	

    this->parallel = this->hcurlfes->GetParallelDofs()!=nullptr;
    NgMPI_Comm comm = this->parallel ? this->hcurlfes->GetParallelDofs()->GetCommunicator() : NgMPI_Comm(MPI_COMM_WORLD);
    this->rank = comm.Rank();
    this->np = comm.Size();
    this->parallel = (np>1);

    cout << IM(2) << "NP: " << np << endl;

    /** Some initial book-keeping **/
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
		       first_master_dof.Data(), 1, MPI_INT,
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
	ScatterDofData (global_nums, pardofs);
	ilower = first_master_dof[this->rank];
	iupper = first_master_dof[this->rank+1]-1;
      };    
    this->hc_freedofs = this->hcurlfes->GetFreeDofs();
    this->hc_ndof = this->hc_freedofs->Size();
    this->h1_freedofs = this->h1fes->GetFreeDofs();
    this->h1_ndof = this->h1_freedofs->Size();
    if(parallel) {
      this->hc_pardofs = this->hcurlfes->GetParallelDofs();
      (void) compute_global_nums(this->hc_ndof, this->hcurlfes, this->hc_pardofs, this->hc_global_nums, this->hc_ilower, this->hc_iupper);
      this->h1_pardofs = this->h1fes->GetParallelDofs();
      (void) compute_global_nums(this->h1_ndof, this->h1fes, this->h1_pardofs, this->h1_global_nums, this->h1_ilower, this->h1_iupper);
      this->hc_masterdofs.SetSize(this->hc_iupper-this->hc_ilower+1);
      int s = 0;
      for(auto k:Range(this->hc_ndof))
	if(this->hc_pardofs->IsMasterDof(k))
	  this->hc_masterdofs[s++] = k;
      this->buf_hc.SetSize(this->hc_masterdofs.Size());
      this->buf_z.SetSize(this->hc_masterdofs.Size());
      this->buf_z = 0.0;
      this->hc_intrange.SetSize(this->hc_iupper-this->hc_ilower+1);
      for(auto k:Range(this->hc_intrange.Size()))
	this->hc_intrange[k] = this->hc_ilower+k;
      ParallelVVector<double> v1(this->hc_pardofs);
      BaseVector & bv1(v1);
      v1.FVDouble() = 0.0;
      for(auto k:Range(this->hc_ndof))
	if(this->hc_freedofs->Test(k))
	  v1.FVDouble()[k] = 1;
      bv1.SetParallelStatus(CUMULATED);
      int fdhc = (int) InnerProduct(bv1,bv1);
      ParallelVVector<double> v2(this->h1_pardofs);
      BaseVector & bv2(v2);
      v2.FVDouble() = 0.0;
      for(auto k:Range(this->h1_ndof))
	if(this->h1_freedofs->Test(k))
	  v2.FVDouble()[k] = 1;
      bv2.SetParallelStatus(CUMULATED);
      int fdh1 = (int) InnerProduct(bv2,bv2);      
      cout << IM(2) << "hcurl freedofs: " << fdhc << " out of " << this->hc_pardofs->GetNDofGlobal() << endl;
      cout << IM(2) << "h1    freedofs: " << fdh1 << " out of " << this->h1_pardofs->GetNDofGlobal() << endl;
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
      this->buf_z.SetSize(this->hc_ndof);
      this->buf_z = 0.0;
      this->hc_intrange.SetSize(this->hc_ndof);
      for(auto k:Range(this->hc_intrange.Size()))
	this->hc_intrange[k] = this->hc_ilower+k;
      cout << IM(2) << "hcurl freedofs: " << this->hc_freedofs->NumSet() << " out of " << this->hc_freedofs->Size() << endl;
      cout << IM(2) << "h1    freedofs: " << this->h1_freedofs->NumSet() << " out of " << this->h1_freedofs->Size() << endl;
    }


    /** Setup AMS-solver **/

    HYPRE_AMSCreate(&this->precond); // create AMS-solver
    HYPRE_AMSSetPrintLevel(this->precond, 1); // print-level, 0=none
    HYPRE_AMSSetDimension(this->precond, this->dimension); // set dimension for solver

    /** Gradient Matrix **/
    const SparseMatrix<double>* spgrad = dynamic_cast<SparseMatrix<double>*>(this->ngs_grad_mat.get());
    (void) Create_IJMat_from_SPMat(&this->grad_mat, &this->parcsr_grad_mat, *spgrad,
				   this->hc_pardofs, this->h1_pardofs,
				   nullptr, nullptr,
				   this->hc_ilower, this->hc_iupper,
				   this->h1_ilower, this->h1_iupper,
				   this->hc_global_nums, this->h1_global_nums, true);
    if( (err = HYPRE_AMSSetDiscreteGradient(this->precond, this->parcsr_grad_mat)) != 0)
      cerr << "id = " << rank << ", HYPRE_AMSSetDiscreteGradient returned with err " << err << endl;
    
    /** alpha-poisson Matrix **/
    if(bf_alpha!=nullptr) {
      // HYPRE_AMSSetAlphaPoissonMatrix(this->precond, this->parcsr_alpha_mat);
    }
    
    /** beta-poisson Matrix**/
    if(this->beta_is_zero) {
      HYPRE_AMSSetBetaPoissonMatrix(this->precond, NULL);
    }
    else if(this->bf_beta!=nullptr) {
      // HYPRE_AMSSetBetaPoissonMatrix(this->precond, this->parcsr_beta_mat);
    }

    /** x, y and z vectors **/
    {
      HYPRE_IJVector h[3];
      HYPRE_ParVector ph[3];
      Array<int> rs(3);
      rs = this->h1_ndof;
      Table<double> c(rs);
      double* cp[3];
      if(h1_ndof)
	for(auto k:Range(3)) cp[k] = c[k].Data();
      else
	for(auto k:Range(3)) cp[k] = NULL;
      for(auto k:Range(3)) c[k] = 0.0;
      auto ma = hcurlfes->GetMeshAccess();
      for(auto k:Range(h1_ndof)) {
	if( parallel && (!this->h1_pardofs->IsMasterDof(k)) ) continue;
	Vec<3> co;
	ma->GetPoint(k, co);
	for(auto j:Range(3)) 
	  c[j][k] = co[j];
      }
      for(auto k:Range(3))
	(void) Create_IJVec_from_BVec (comm, h[k], ph[k], cp[k], this->h1_global_nums, this->h1_ilower, this->h1_iupper, true);
      if( (err = HYPRE_AMSSetCoordinateVectors(this->precond, ph[0], ph[1], ph[2])) != 0)
	cerr << "HYPRE_AMSSetCoordinateVectors returned with error " << err << endl;      
    }

    /** working vectors **/
    Array<double> zeros(this->hc_ndof+1); //+1 for rank 0...
    zeros = 0.0;
    (void) Create_IJVec_from_BVec (comm, this->b, this->par_b,  zeros.Data(), this->hc_global_nums, this->hc_ilower, this->hc_iupper, true);
    (void) Create_IJVec_from_BVec (comm, this->x, this->par_x,  zeros.Data(), this->hc_global_nums, this->hc_ilower, this->hc_iupper, true);
    
    /** main system matrix **/
    const auto & matrix = this->bfa->GetMatrix();
    const ParallelMatrix* pma = dynamic_cast<const ParallelMatrix*>(&matrix);
    const SparseMatrix<double>* spma = dynamic_cast<const SparseMatrix<double>*>( (pma==nullptr) ? (&matrix) : (pma->GetMatrix().get()) );
    if (dynamic_cast< const SparseMatrixSymmetric<double> *> (spma))
      throw Exception ("Please use fully stored sparse matrix for hypre (bf -nonsymmetric)");
    (void) Create_IJMat_from_SPMat(&this->A, &this->parcsr_A, *spma,
				   this->hc_pardofs, this->hc_pardofs,
				   this->hc_freedofs, this->hc_freedofs,
				   this->hc_ilower, this->hc_iupper,
				   this->hc_ilower, this->hc_iupper,
				   this->hc_global_nums, this->hc_global_nums);
    HYPRE_AMSSetup(this->precond, this->parcsr_A, this->par_b, this->par_x); //call setup
    
    
    /** Set a couple of AMS options **/

    // ams options
    HYPRE_AMSSetTol(this->precond, 0);
    HYPRE_AMSSetMaxIter(this->precond,1);

    // HYPRE_AMSSetCycleType(this->precond,1);
    // HYPRE_AMSSetSmoothingOptions(this->precond, 4, 1, 1.0, 1.0);      //truncated l1-smoother
    // HYPRE_AMSSetBetaAMGOptions(this->precond, 10, 10, 6, 0.5, 0, 0);  //symmetric l1-smoother
    // HYPRE_AMSSetAlphaAMGOptions(this->precond, 10, 10, 6, 0.5, 0, 0); //symmetric l1-smoother

    // HYPRE_IJMatrixPrint(A, "IJ.out.2A");
    // HYPRE_IJMatrixPrint(grad_mat, "IJ.out.2grad_mat");
    // HYPRE_IJMatrixPrint(alpha_mat, "IJ.out.alpha");
    // HYPRE_IJMatrixPrint(beta_mat, "IJ.out.beta");

    cout << IM(1) << "Hypre-AMS setup done!" << endl;
    return;
  }

  void HypreAMSPreconditioner :: Mult (const BaseVector & f, BaseVector & u) const
  {
    static Timer t("AMS-hypre mult");
    RegionTimer reg(t);
    
    if(parallel) {
      f.Cumulate();
      u.SetParallelStatus(DISTRIBUTED);
    }
    
    /** re-initialize vectors **/
    HYPRE_IJVectorInitialize(this->b);
    HYPRE_IJVectorInitialize(this->x);
    if(this->hc_ndof) {
      HYPRE_IJVectorSetValues(this->x, this->hc_intrange.Size(), this->hc_intrange.Data(), this->buf_z.Data());
      for(auto k:Range(this->buf_hc.Size()))
	if(this->hc_freedofs->Test(this->hc_masterdofs[k]))
	  buf_hc[k] = f.FVDouble()[this->hc_masterdofs[k]];
	else
	  buf_hc[k] = 0.0;
      HYPRE_IJVectorSetValues(this->b, this->hc_intrange.Size(), this->hc_intrange.Data(), this->buf_hc.Data());
    }
    HYPRE_IJVectorAssemble(this->b);    
    HYPRE_IJVectorGetObject(this->b, (void **) &this->par_b);
    HYPRE_IJVectorAssemble(this->x);
    HYPRE_IJVectorGetObject(this->x, (void **) &this->par_x);

    /** call AMS-solver **/
    HYPRE_AMSSolve(this->precond, this->parcsr_A, this->par_b, this->par_x);
    
    /** get sol **/
    HYPRE_IJVectorGetValues(x, this->hc_intrange.Size(), this->hc_intrange.Data(), this->buf_hc.Data());    
    if(this->hc_ndof) u.FVDouble() = 0.0;    
    for(auto k:Range(this->buf_hc.Size()))
      u.FVDouble()[this->hc_masterdofs[k]] = this->buf_hc[k];

    /*
    // check for positiveness of PC
      auto ip = InnerProduct(u,f);    
    if(ip<0)
      cout << IM(1) << "                IP: " << ip << "  !!!!!!!!" << endl;
    else
      cout << IM(1) << "                IP: " << ip << endl;
    */
    
    return;
  }

  // don't need to register, those constructors are not implemented anyways...
  static RegisterPreconditioner<HypreAMSPreconditioner> init_hyprepre ("hypre_ams");
}
#endif
