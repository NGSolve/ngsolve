#ifdef HYPRE

#ifndef FILE_AMS_HYPRE_PRECOND
#define FILE_AMS_HYPRE_PRECOND


/*********************************************************************/
/* File:   hypre_precond.hpp                                         */
/* Author: Lukas Kolger                                              */
/* Date:   May 2017                                                  */
/*********************************************************************/


#include "HYPRE_utilities.h"

#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"


namespace ngcomp
{

  
class HypreAMSPreconditioner : public Preconditioner
{

  // HYPRE-side

  HYPRE_Solver precond; // the AMS solver
  
  // discrete gradient matrix
  HYPRE_IJMatrix grad_mat;
  HYPRE_ParCSRMatrix parcsr_grad_mat;

  // stiffness matrix alpha<curl,curl>+beta<,>
  HYPRE_IJMatrix A; 
  HYPRE_ParCSRMatrix parcsr_A;

  // alpha <grad,grad>
  HYPRE_IJMatrix alpha_mat;
  HYPRE_ParCSRMatrix parcsr_alpha_mat;

  // beta <grad, grad>
  bool beta_is_zero = false;
  HYPRE_IJMatrix beta_mat;
  HYPRE_ParCSRMatrix parcsr_beta_mat;

  // working vecs
  HYPRE_IJVector b;
  HYPRE_ParVector par_b;
  HYPRE_IJVector x;
  HYPRE_ParVector par_x;



  // NGSOLVE-side

  int rank, np;
  
  shared_ptr<BilinearForm> bfa; // BLF for the system
  shared_ptr<BilinearForm> bf_alpha; // alpha-mat - currently UNUSED!
  shared_ptr<BilinearForm> bf_beta; // beta-mat - currently UNUSED!
  shared_ptr<BaseMatrix> ngs_grad_mat; //discrete grad-mat

  int dimension = 3; //spatial dimension

  shared_ptr<FESpace> hcurlfes;
  shared_ptr<BitArray> hc_freedofs;
  shared_ptr<ParallelDofs> hc_pardofs;
  int hc_ndof;
  Array<int> hc_global_nums;
  int hc_ilower, hc_iupper;
  Array<int> hc_masterdofs;
  Array<double> buf_hc;
  Array<double> buf_z;
  Array<int> hc_intrange;

  shared_ptr<FESpace> h1fes;
  shared_ptr<BitArray> h1_freedofs;
  shared_ptr<ParallelDofs> h1_pardofs;
  int h1_ndof;
  Array<int> h1_global_nums;
  int h1_ilower, h1_iupper;

  bool parallel;


public:

  // dummy constructors - not implemented
  HypreAMSPreconditioner (const BaseMatrix & matrix, const shared_ptr<BitArray> afreedofs); 

  // actually implemented!
  HypreAMSPreconditioner (shared_ptr<BilinearForm> bfa, const Flags & aflags,
			  const string aname = "precond");
  HypreAMSPreconditioner (shared_ptr<FESpace> ahcurlfes, shared_ptr<BilinearForm> abfa,
			  shared_ptr<FESpace> ah1fes, shared_ptr<BaseMatrix> agrad_mat,
			  shared_ptr<BilinearForm> abf_alpha, shared_ptr<BilinearForm> abf_beta, 
			  Flags & flags);

  ~HypreAMSPreconditioner ();
	
  virtual void FinalizeLevel (const ngla::BaseMatrix * mat = NULL) override;
  virtual void Update() override;
  virtual void Mult (const BaseVector & f, BaseVector & u) const override;
  virtual int VHeight() const override { return hc_ndof; }
  virtual int VWidth() const override { return hc_ndof; }
  virtual const BaseMatrix & GetAMatrix() const override { return bfa->GetMatrix(); }
  virtual const BaseMatrix & GetMatrix() const override { return *this; }

  virtual const char * ClassName() const override
  { return "HYPRE AMG Preconditioner"; }

private:
  void Setup ();
};

}



#endif
#endif
