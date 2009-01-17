#ifndef FILE_SMOOTHER
#define FILE_SMOOTHER

/*********************************************************************/
/* File:   smoother.hh                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   20. Apr. 2000                                             */
/*********************************************************************/

/** 
    Smoothing iteration for multigrid method.
    Pure virtual base class.
*/
class Smoother
{
protected:
  /// additive or multiplicative smooting
  int additive; 
  /// should coarse levels be updated, too ?
  int updateall;
  /// 
  Flags flags;
public: 
  /// Constructor
  Smoother();
  /// Constructor
  Smoother(const Flags & aflags); 
  /// Destructor
  virtual ~Smoother();
  
  /// Update smoother (fine level or all levels)
  virtual void Update (bool force_update = 0) = 0;

  /// Do steps iterations of pre-smoothing
  virtual void PreSmooth (int level, ngla::BaseVector & u, 
			  const ngla::BaseVector & f, int steps) const = 0;


  /// Do steps iterations of pre-smoothing
  virtual void PreSmoothResiduum (int level, ngla::BaseVector & u, 
				  const ngla::BaseVector & f, 
				  ngla::BaseVector & res, 
				  int steps) const 
  {
    PreSmooth (level, u, f, steps);
    Residuum (level, u, f, res);
  }


  /// Do steps iterations of post-smoothing
  virtual void PostSmooth (int level, ngla::BaseVector & u, 
		     const ngla::BaseVector & f, int steps) const = 0;
  /// Apply the preconditioning action (additive or multiplicative)
  virtual void Precond (int level, const ngla::BaseVector & f, ngla::BaseVector & u) const;

  /// 
  virtual void Residuum (int level, ngla::BaseVector & u, 
			 const ngla::BaseVector & f, 
			 ngla::BaseVector & d) const = 0;

  ///
  void SetAdditive () { additive = 1; }
  ///
  void SetUpdateAll (int ua) { updateall = ua; }
  ///
  void SetMultiplicative () { additive = 0; }
  ///
  int Additive () const { return additive; }

  ///
  virtual ngla::BaseVector * CreateVector(int level) const = 0;

  virtual void MemoryUsage (ARRAY<MemoryUsageStruct*> & mu) const
  { ; }
};


/**
   Gauss-Seidel smoother.
   Common relaxation of unknwons in node.
 */
class GSSmoother : public Smoother
{
  ///
  const MeshAccess & ma;
  ///
  const BilinearForm & biform;
  ///
  ARRAY<ngla::BaseJacobiPrecond*> jac;
  
public:
  ///
  GSSmoother (const MeshAccess & ama,
	      const BilinearForm & abiform);
  ///
  virtual ~GSSmoother();
  
  ///
  virtual void Update (bool force_update = 0);
  ///
  virtual void PreSmooth (int level, ngla::BaseVector & u, 
			  const ngla::BaseVector & f, int steps) const;
  ///
  virtual void PostSmooth (int level, ngla::BaseVector & u, 
			   const ngla::BaseVector & f, int steps) const;

  virtual void PreSmoothResiduum (int level, ngla::BaseVector & u, 
				  const ngla::BaseVector & f, 
				  ngla::BaseVector & res, 
				  int steps) const;

  ///
  virtual void Residuum (int level, ngla::BaseVector & u, 
			 const ngla::BaseVector & f, ngla::BaseVector & d) const;
  ///
  virtual ngla::BaseVector * CreateVector(int level) const;
};


/**
   Anisotropic smoother.
   Common relaxation of vertically aligned nodes.
 */
class AnisotropicSmoother : public Smoother
{
  ///
  const MeshAccess & ma;
  ///
  const BilinearForm & biform;
  ///
  ARRAY<BaseBlockJacobiPrecond*> jac;
  
public:
  ///
  AnisotropicSmoother (const MeshAccess & ama,
		       const BilinearForm & abiform);
  ///
  virtual ~AnisotropicSmoother();
  
  ///
  virtual void Update (bool forace_update = 0);
  ///
  virtual void PreSmooth (int level, ngla::BaseVector & u, 
			  const ngla::BaseVector & f, int steps) const;
  ///
  virtual void PostSmooth (int level, ngla::BaseVector & u, 
			   const ngla::BaseVector & f, int steps) const;
  ///
  virtual void Residuum (int level, ngla::BaseVector & u, 
			 const ngla::BaseVector & f, ngla::BaseVector & d) const;
  ///
  virtual ngla::BaseVector * CreateVector(int level) const;
};



/**
   Element-by-element smoothing.
   Useful ?
 */
class EBESmoother : public Smoother
{
  ///
  const MeshAccess & ma;
  ///
  const BilinearForm & biform;
  ///
  ARRAY<BaseBlockJacobiPrecond*> jac;
  
public:
  ///
  EBESmoother (const MeshAccess & ama,
		       const BilinearForm & abiform);
  ///
  virtual ~EBESmoother();
  
  ///
  virtual void Update (bool forace_update = 0);
  ///
  virtual void PreSmooth (int level, ngla::BaseVector & u, 
			  const ngla::BaseVector & f, int steps) const;
  ///
  virtual void PostSmooth (int level, ngla::BaseVector & u, 
			   const ngla::BaseVector & f, int steps) const;
  ///
  virtual void Residuum (int level, ngla::BaseVector & u, 
			 const ngla::BaseVector & f, ngla::BaseVector & d) const;
  ///
  virtual ngla::BaseVector * CreateVector(int level) const;
};



/*
///   Arnold-Falk-Winther smoother for 3D H(curl) problems.
class EdgeSmoother : public Smoother
{
  ///
  const MeshAccess & ma;
  ///
  const NedelecFESpace & space;
  ///
  const BilinearForm & biform;
  ///
  ARRAY<BaseBlockJacobiPrecond*> jac;
  
public:
  ///
  EdgeSmoother (const MeshAccess & ama,
		const NedelecFESpace & aspace,
		const BilinearForm & abiform);
  ///
  virtual ~EdgeSmoother();
  
  ///
  virtual void Update ();
  ///
  virtual void PreSmooth (int level, ngla::BaseVector & u, 
			  const ngla::BaseVector & f, int steps) const;
  ///
  virtual void PostSmooth (int level, ngla::BaseVector & u, 
			   const ngla::BaseVector & f, int steps) const;
  ///
  virtual void Precond (int level, const ngla::BaseVector & f, ngla::BaseVector & u) const;
  ///
  virtual void Residuum (int level, ngla::BaseVector & u, 
			 const ngla::BaseVector & f, ngla::BaseVector & d) const;
  ///
  virtual ngla::BaseVector * CreateVector(int level) const;

  ///
  friend class EVCoarseGrid;
};





///   R. Hiptmair's smoother for 3D H(curl) problems
class HiptmairSmoother : public Smoother
{
protected:
  ///
  const MeshAccess & ma;
  ///
  const NedelecFESpace & space;
  ///
  const BilinearForm & biform;
  ///
  const BilinearForm & biformpot;
  ///
  ARRAY<BaseBlockJacobiPrecond*> jac;
  ///
  ARRAY<BaseBlockJacobiPrecond*> jacpot;
  ///
  ngla::BaseMatrix * potcoarseinv;
  
public:
  ///
  HiptmairSmoother (const MeshAccess & ama,
		    const NedelecFESpace & aspace,
		    const BilinearForm & abiform,
		    const BilinearForm & abiformpot);
  ///
  virtual ~HiptmairSmoother();
  
  ///
  virtual void Update ();
  ///
  virtual void PreSmooth (int level, ngla::BaseVector & u, 
			  const ngla::BaseVector & f, int steps) const;
  ///
  virtual void PostSmooth (int level, ngla::BaseVector & u, 
			   const ngla::BaseVector & f, int steps) const;
  ///
  virtual void Precond (int level, const ngla::BaseVector & f, ngla::BaseVector & u) const;
  ///
  virtual void Residuum (int level, ngla::BaseVector & u, 
			 const ngla::BaseVector & f, ngla::BaseVector & d) const;
  ///
  virtual ngla::BaseVector * CreateVector(int level) const;

  ///
  void AddGradient (double fac, int level, const ngla::BaseVector & pot, ngla::BaseVector & grad) const;
  ///
  void ApplyGradientT (int level, const ngla::BaseVector & gradt, ngla::BaseVector & pott) const;
};

*/







/**
   Block-Gauss-Seidel smoother.
   Blocks are defined by underlying FESpace.
 */
class BlockSmoother : public Smoother
{
  ///
  const MeshAccess & ma;
  ///
  const BilinearForm & biform;
  ///
  const LinearForm * constraint;
  ///
  ARRAY<ngla::BaseBlockJacobiPrecond*> jac;
  ///
  ARRAY<BaseMatrix*> inv;
  ///
  ARRAY<int> * direct;
  
public:
  ///
  BlockSmoother (const MeshAccess & ama,
		 const BilinearForm & abiform, const Flags & aflags);
  ///
  BlockSmoother (const MeshAccess & ama,
		 const BilinearForm & abiform,
		 const LinearForm & aconstraint, const Flags & aflags);
  ///
  virtual ~BlockSmoother();
  
  ///
  virtual void Update (bool forace_update = 0);
  ///
  virtual void PreSmooth (int level, BaseVector & u, 
			  const BaseVector & f, int steps) const;
  ///
  virtual void PreSmoothResiduum (int level, ngla::BaseVector & u, 
				  const ngla::BaseVector & f, 
				  ngla::BaseVector & res, 
				  int steps) const;

  ///
  virtual void PostSmooth (int level, ngla::BaseVector & u, 
			   const ngla::BaseVector & f, int steps) const;
  ///
  virtual void Precond (int level, const ngla::BaseVector & f, ngla::BaseVector & u) const;
  ///
  virtual void Residuum (int level, ngla::BaseVector & u, 
			 const ngla::BaseVector & f, ngla::BaseVector & d) const;
  ///
  virtual ngla::BaseVector * CreateVector(int level) const;

  virtual void MemoryUsage (ARRAY<MemoryUsageStruct*> & mu) const;
};



/**
   Hiptmair's smoother for H(curl)
 */
class PotentialSmoother : public Smoother
{
  ///
  const MeshAccess & ma;
  ///
  const BilinearForm & biform;
  ///
  ARRAY<BaseBlockJacobiPrecond*> jac;
  ///
  ARRAY<SparseMatrix<double>*> gradient;
  ///
  ARRAY<BaseSparseMatrix*> potmat;
  ///
  ARRAY<BaseBlockJacobiPrecond*> potjac;
  ///
  ARRAY<BaseMatrix*> inv;
public:
  ///
  PotentialSmoother (const MeshAccess & ama,
		     const BilinearForm & abiform);
  ///
  virtual ~PotentialSmoother();
  
  ///
  virtual void Update (bool forace_update = 0);
  ///
  virtual void PreSmooth (int level, ngla::BaseVector & u, 
			  const ngla::BaseVector & f, int steps) const;
  ///
  virtual void PostSmooth (int level, ngla::BaseVector & u, 
			   const ngla::BaseVector & f, int steps) const;
  ///
  virtual void Precond (int level, const ngla::BaseVector & f, ngla::BaseVector & u) const;
  ///
  virtual void Residuum (int level, ngla::BaseVector & u, 
			 const ngla::BaseVector & f, ngla::BaseVector & d) const;
  ///
  virtual BaseVector * CreateVector(int level) const;

  virtual void MemoryUsage (ARRAY<MemoryUsageStruct*> & mu) const;
};






/*
///   Experimental smoother
class StabEdgeSmoother : public Smoother
{
  ///
  const MeshAccess & ma;
  ///
  const NedelecFESpace & space;
  ///
  const NodalFESpace & nodalspace;
  ///
  const BilinearForm & biform;
  ///
  const BilinearForm & biforml2;
  ///
  ARRAY<ngla::BaseJacobiPrecond*> jacl2;
  ///
  ARRAY<ngla::BaseJacobiPrecond*> jacedge;
  ///
  ARRAY<double> lami;
public:
  ///
  StabEdgeSmoother (const MeshAccess & ama,
		    const NedelecFESpace & aspace,
		    const BilinearForm & abiform,
		    const NodalFESpace & anodalspace,
		    const BilinearForm & abiforml2);
  ///
  virtual ~StabEdgeSmoother();
  
  ///
  virtual void Update ();
  ///
  virtual void PreSmooth (int level, ngla::BaseVector & u, 
			  const ngla::BaseVector & f, int steps) const;
  ///
  virtual void PostSmooth (int level, ngla::BaseVector & u, 
			   const ngla::BaseVector & f, int steps) const;
  ///
  virtual void Residuum (int level, ngla::BaseVector & u, 
			 const ngla::BaseVector & f, ngla::BaseVector & d) const;
  ///
  virtual ngla::BaseVector * CreateVector(int level) const;
};
*/






/**
   Local, symmetric Uzawa iteration.
 */
class MixedSmoother : public Smoother
{
  ///
  const MeshAccess & ma;
  ///
  const BilinearForm & biforma;
  ///
  const BilinearForm & biformb;
  ///
  ARRAY<ngla::BaseJacobiPrecond*> jac;
  
  ///
  ARRAY<double> tau1;
  ///
  ARRAY<double> tau2;
public:

  ///
  MixedSmoother (const MeshAccess & ama,
		 const BilinearForm & abiforma,
		 const BilinearForm & abiformb);
  ///
  virtual ~MixedSmoother();
  
  ///
  virtual void Update (bool forace_update = 0);
  ///
  virtual void PreSmooth (int level, ngla::BaseVector & u, 
			  const ngla::BaseVector & f, int steps) const;
  ///
  virtual void PostSmooth (int level, ngla::BaseVector & u, 
			   const ngla::BaseVector & f, int steps) const;
  ///
  virtual void Residuum (int level, ngla::BaseVector & u, 
			 const ngla::BaseVector & f, ngla::BaseVector & d) const;
  ///
  virtual ngla::BaseVector * CreateVector(int level) const;
};



/*
///   Vanka smoother.
class VankaSmoother : public Smoother
{
  ///
  const MeshAccess & ma;
  ///
  const BilinearForm & biforma;
  ///
  const BilinearForm & biformb;

  ///
  ARRAY<SaddlePointJacobiPrecond<SysMatrix2d,SysVector2d,SysVector1d> *> jac;
public:

  ///
  VankaSmoother (const MeshAccess & ama,
		 const BilinearForm & abiforma,
		 const BilinearForm & abiformb);
  ///
  virtual ~VankaSmoother();
  
  ///
  virtual void Update ();
  ///
  virtual void PreSmooth (int level, ngla::BaseVector & u, 
			  const ngla::BaseVector & f, int steps) const;
  ///
  virtual void PostSmooth (int level, ngla::BaseVector & u, 
			   const ngla::BaseVector & f, int steps) const;
  ///
  virtual void Residuum (int level, ngla::BaseVector & u, 
			 const ngla::BaseVector & f, ngla::BaseVector & d) const;
  ///
  virtual ngla::BaseVector * CreateVector(int level) const;
};
*/



/**
   Matrix - vector multiplication by smoothing step.
 */
class SmoothingPreconditioner : public ngla::BaseMatrix
{
  ///
  const Smoother & smoother;
  ///
  int level;
public:
  ///
  SmoothingPreconditioner (const Smoother & asmoother,
			   int alevel = 0);
  ///
  virtual void Mult (const ngla::BaseVector & f, ngla::BaseVector & u) const;
  ///
  virtual ngla::BaseVector * CreateVector () const;
};


#endif
