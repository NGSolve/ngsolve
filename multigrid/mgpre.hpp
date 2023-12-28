
#ifndef FILE_MGPRE
#define FILE_MGPRE

/*********************************************************************/
/* File:   mgpre.hh                                                  */
/* Author: Joachim Schoeberl                                         */
/* Date:   20. Apr. 2000                                             */
/*********************************************************************/

#include <bilinearform.hpp>
#include <la.hpp>

namespace ngmg
{
  using namespace ngla;
  using namespace ngcomp;
  using ngla::BaseMatrix;
  
  /** 
      Multigrid preconditioner
  */

  class Smoother;
  ///
  class Prolongation;
  

  ///
  class NGS_DLL_HEADER MultigridPreconditioner : public BaseMatrix
  {

  public:
    enum COARSETYPE { EXACT_COARSE, CG_COARSE, SMOOTHING_COARSE, USER_COARSE };

  private:
    shared_ptr<BilinearForm> biform;
    shared_ptr<MeshAccess> ma;
  
    shared_ptr<Smoother> smoother;
    shared_ptr<Prolongation> prolongation;
    shared_ptr<BaseMatrix> coarsegridpre;
    
    ///
    double checksumcgpre;
    ///
    COARSETYPE coarsetype;
    ///
    int cycle, incsmooth, smoothingsteps;
    ///
    int coarsesmoothingsteps;
    ///
    int updateall;
    /// creates a new smoother for each update
    bool update_always; 
    /// for robust prolongation
    bool harmonic_extension_prolongation = false;
    Array<shared_ptr<BaseMatrix>> he_prolongation;
  public:
    ///
    MultigridPreconditioner (shared_ptr<BilinearForm> abiform,
			     shared_ptr<Smoother> asmoother,
			     shared_ptr<Prolongation> aprolongation);
    ///
    ~MultigridPreconditioner ();
    ///
    virtual bool IsComplex() const override
    { return biform->GetFESpace()->IsComplex(); }

    ///
    void SetSmoothingSteps (int sstep);
    ///
    void SetCycle (int c);
    ///
    void SetIncreaseSmoothingSteps (int incsm);
    ///
    void SetCoarseType (COARSETYPE ctyp);
    ///
    void SetCoarseGridPreconditioner (shared_ptr<BaseMatrix> acoarsegridpre);
    ///
    void SetCoarseSmoothingSteps (int cstep);

    void SetUpdateAll (int ua = 1);
    ///
    void SetUpdateAlways (bool ua = 1) { update_always = ua; }
    ///
    void SetHarmonicExtensionProlongation (bool he = true)
    { harmonic_extension_prolongation = he; }
    
    ///
    virtual void Update () override;

    ///
    virtual void Mult (const BaseVector & x, BaseVector & y) const override;

    ///
    void MGM (int level, BaseVector & u, 
	      const BaseVector & f, int incsm = 1) const;
    ///
    AutoVector CreateRowVector () const override
    { return biform->GetMatrix().CreateColVector(); }
    AutoVector CreateColVector () const override
    { return biform->GetMatrix().CreateRowVector(); }
  
    ///
    const Smoother & GetSmoother() const
    { return *smoother; }
    ///
    Smoother & GetSmoother()
    { return *smoother; }
    ///
    const Prolongation & GetProlongation() const
    { return *prolongation; }


    virtual int VHeight() const override
    {
      return biform->GetMatrix().Height();
    }

    virtual int VWidth() const override
    {
      return biform->GetMatrix().VWidth();
    }

    virtual Array<MemoryUsage> GetMemoryUsage () const override;

  private:
    MemoryTracer mt = { "MultiGridPreconditioner" };
  public:
    const MemoryTracer& GetMemoryTracer() const { return mt; }
  };







  ///
  class NGS_DLL_HEADER TwoLevelMatrix : public BaseMatrix
  {
    ///
    const BaseMatrix * mat;
    ///
    const BaseMatrix * cpre;
    ///
    shared_ptr<Smoother> smoother;
    shared_ptr<BaseMatrix> embedding;
    ///
    int level;
    ///
    int smoothingsteps;
  public:
    ///
    TwoLevelMatrix (const BaseMatrix * amat, 
		    const BaseMatrix * acpre, 
		    shared_ptr<Smoother> asmoother, int alevel);
    ///
    ~TwoLevelMatrix ();

    void SetEmbedding (shared_ptr<BaseMatrix> aembedding)
    {
      embedding = aembedding;
    }
    
    virtual bool IsComplex() const override { return mat->IsComplex(); }

    ///

    virtual void Mult (const BaseVector & x, BaseVector & y) const override;
    ///
    AutoVector CreateRowVector () const override { return mat->CreateColVector(); }
    AutoVector CreateColVector () const override { return mat->CreateRowVector(); }
    ///
    virtual ostream & Print (ostream & s) const override;
    ///
    void SetSmoothingSteps(int ass) { smoothingsteps = ass; }
    ///
    const Smoother & GetSmoother() const
    { return *smoother; }
    ///
    Smoother & GetSmoother()
    { return *smoother; }
    ///
    virtual void Update() override;

    virtual int VHeight() const override
    {
      return mat->Height();
    }

    virtual int VWidth() const override
    {
      return mat->VWidth();
    }
  

    virtual Array<MemoryUsage> GetMemoryUsage () const override;
  };

}

#endif
