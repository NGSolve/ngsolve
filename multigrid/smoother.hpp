#ifndef FILE_SMOOTHER
#define FILE_SMOOTHER

/*********************************************************************/
/* File:   smoother.hh                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   20. Apr. 2000                                             */
/*********************************************************************/

namespace ngmg
{

  /** 
      Smoothing iteration for multigrid method.
      Pure virtual base class.
  */
  class NGS_DLL_HEADER Smoother
  {
  protected:
    /// additive or multiplicative smoothing
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
    virtual AutoVector CreateVector(int level) const = 0;

    virtual Array<MemoryUsage> GetMemoryUsage () const { return Array<MemoryUsage>(); }
  private:
  MemoryTracer mt = { "Smoother" };
  public:
  const MemoryTracer& GetMemoryTracer() const { return mt; }
  };


  /**
     Gauss-Seidel smoother.
     Common relaxation of unknowns in node.
  */
  class GSSmoother : public Smoother
  {
    ///
    // const MeshAccess & ma;
    ///
    const BilinearForm & biform;
    ///
    Array<shared_ptr<BaseJacobiPrecond>> jac;
  
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
    virtual AutoVector CreateVector(int level) const;
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
    Array<shared_ptr<BaseBlockJacobiPrecond>> jac;
  
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
    virtual AutoVector CreateVector(int level) const;
  };





  /**
     Block-Gauss-Seidel smoother.
     Blocks are defined by underlying FESpace.
  */
  class BlockSmoother : public Smoother
  {
    ///
    // const MeshAccess & ma;
    ///
    const BilinearForm & biform;
    ///
    const LinearForm * constraint;
    ///
    Array<shared_ptr<BaseBlockJacobiPrecond>> jac;
    ///
    Array<shared_ptr<BaseMatrix>> inv;
    ///
    shared_ptr<Array<int>> direct;
    shared_ptr<Array<int>> userdefined_direct;

    Array<shared_ptr<Table<int>>> smoothing_blocks;

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
    virtual AutoVector CreateVector(int level) const;

    virtual Array<MemoryUsage> GetMemoryUsage () const;

    void SetDirectSolverCluster(shared_ptr<Array<int>> cluster)
    {  userdefined_direct = cluster; }
  };





#ifdef XXX_OBSOLETE
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
    virtual AutoVector CreateVector () const;
  };
#endif

}

#endif
