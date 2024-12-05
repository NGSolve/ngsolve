#ifndef FILE_CG_NEW
#define FILE_CG_NEW


/**************************************************************************/
/* File:   cg.hh                                                          */
/* Author: Joachim Schoeberl                                              */
/* Date:   5. Jul. 96                                                     */
/**************************************************************************/


#include "basematrix.hpp"

namespace ngla
{


  /**
     Krylov Space Solver
  */ 
  class KrylovSpaceSolver : public BaseMatrix
  {
  protected:
    ///
    shared_ptr<BaseMatrix> a, c;
    ///
    double prec;
    ///
    int maxsteps;
    ///
    int steps;
    ///
    int initialize;
    ///
    bool stop_absolute;
    ///
    int printrates;
    ///
    int absoluteRes;
    ///
    bool useseed;

  public:
    ///
    NGS_DLL_HEADER KrylovSpaceSolver();
    ///
    NGS_DLL_HEADER KrylovSpaceSolver(shared_ptr<BaseMatrix> aa);
    ///
    NGS_DLL_HEADER KrylovSpaceSolver(shared_ptr<BaseMatrix> aa, shared_ptr<BaseMatrix> ac);
    ///
    void SetMatrix (shared_ptr<BaseMatrix> aa)
    { a = aa; }
    ///
    void SetPrecond (shared_ptr<BaseMatrix> ac)
    { c = ac; }
    ///

    bool IsComplex() const override { return a->IsComplex(); }
    /// 
    void SetMaxSteps (int amaxsteps)
    { maxsteps = amaxsteps; }

    ///
    void SetPrecision (double aprec)
    { prec = aprec; stop_absolute = 0; }
    ///
    void SetAbsolutePrecision (double aprec)
    {  prec = aprec; stop_absolute = 1; }
    ///
    void SetRelativePrecision (double aprec)
    {  prec = aprec; stop_absolute = 0; }

    double GetPrecision() const { return prec; }
    int GetMaxSteps() const { return maxsteps; }


    void SetPrintRates (int pr = 1)
    { printrates = pr; }
    ///
    void SetInitialize (int ai)
    { initialize = ai; }
    ///

    void UseSeed(const bool useit = true)
    { useseed = useit; }

    ///
    int GetSteps () const
    { return steps; }

    ///
    NGS_DLL_HEADER AutoVector CreateRowVector() const override { return a->CreateColVector(); }
    NGS_DLL_HEADER AutoVector CreateColVector() const override { return a->CreateRowVector(); }


    int VHeight() const override {return a->VWidth();}
    int VWidth() const override {return a->VHeight();}
  };
  

  /// The conjugate gradient solver
  template <class IPTYPE>
  class CGSolver : public KrylovSpaceSolver
  {
  protected:
    ///
    void MultiMult (const BaseVector & f, BaseVector & u, const int dim) const;
    ///
    void MultiMultSeed (const BaseVector & f, BaseVector & u, const int dim) const;
  public:
    typedef typename SCAL_TRAIT<IPTYPE>::SCAL SCAL;
    ///
    CGSolver () 
      : KrylovSpaceSolver () { ; }
    ///
    CGSolver (shared_ptr<BaseMatrix> aa)
      : KrylovSpaceSolver (aa) { ; }

    ///
    CGSolver (shared_ptr<BaseMatrix> aa, shared_ptr<BaseMatrix> ac)
      : KrylovSpaceSolver (aa, ac) { ; }

    ///
    NGS_DLL_HEADER virtual void Mult (const BaseVector & v, BaseVector & prod) const;
  };


  /// The BiCGStab solver
  template <class IPTYPE>
  class NGS_DLL_HEADER BiCGStabSolver : public KrylovSpaceSolver
  {
  public:
    typedef typename SCAL_TRAIT<IPTYPE>::SCAL SCAL;
    ///
    BiCGStabSolver () 
      : KrylovSpaceSolver () { ; }
    ///
    BiCGStabSolver (shared_ptr<BaseMatrix> aa)
      : KrylovSpaceSolver (aa) { ; }

    ///
    BiCGStabSolver (shared_ptr<BaseMatrix> aa, shared_ptr<BaseMatrix> ac)
      : KrylovSpaceSolver (aa, ac) { ; }

    ///
    virtual void Mult (const BaseVector & v, BaseVector & prod) const;
  };
  




  //   Simple iteration solver
  template <class IPTYPE>
  class NGS_DLL_HEADER SimpleIterationSolver : public KrylovSpaceSolver
  {
  public:
    typedef typename SCAL_TRAIT<IPTYPE>::SCAL SCAL;
  private:
    SCAL tau;
  public:
    ///
    SimpleIterationSolver ()
      : KrylovSpaceSolver() { tau = 1; }
    ///
    SimpleIterationSolver (shared_ptr<BaseMatrix> aa)
      : KrylovSpaceSolver (aa) { tau = 1; }
    ///
    SimpleIterationSolver (shared_ptr<BaseMatrix> aa, shared_ptr<BaseMatrix> ac)
      : KrylovSpaceSolver (aa, ac) { tau = 1; }
    ///
    virtual void Mult (const BaseVector & v, BaseVector & prod) const;

    void SetTau (SCAL atau) { tau = atau; }
  };



  /// The conjugate gradient solver
  template <class IPTYPE>
  class NGS_DLL_HEADER GMRESSolver : public KrylovSpaceSolver
  {
  public:
    typedef typename SCAL_TRAIT<IPTYPE>::SCAL SCAL;
    ///
    GMRESSolver () 
      : KrylovSpaceSolver () { ; }
    ///
    GMRESSolver (shared_ptr<BaseMatrix> aa)
      : KrylovSpaceSolver (aa) { ; }

    ///
    GMRESSolver (shared_ptr<BaseMatrix> aa, shared_ptr<BaseMatrix> ac)
      : KrylovSpaceSolver (aa, ac) { ; }

    ///
    virtual void Mult (const BaseVector & v, BaseVector & prod) const;
  };
  




  /// The quasi-minimal residual (QMR) solver
  template <class IPTYPE>
  class NGS_DLL_HEADER QMRSolver : public KrylovSpaceSolver
  {
    int status;
    const BaseMatrix * c2;
  public:
    typedef typename SCAL_TRAIT<IPTYPE>::SCAL SCAL;
    ///
    QMRSolver () 
      : KrylovSpaceSolver (), c2(0) { ; }
    ///
    QMRSolver (shared_ptr<BaseMatrix> aa)
      : KrylovSpaceSolver (aa), c2(0) { ; }

    ///
    QMRSolver (shared_ptr<BaseMatrix> aa, shared_ptr<BaseMatrix> ac)
      : KrylovSpaceSolver (aa, ac), c2(0) { ; }

    ///
    virtual void Mult (const BaseVector & v, BaseVector & prod) const;
  };
  




  /*

  // Conjugate residual solver for symmetric, indefinite matrices


  class CRSolver : public CGSolver
  {
  public:
  ///
  CRSolver ();
  ///
  CRSolver (const BaseMatrix & aa);
  ///
  CRSolver (const BaseMatrix & aa, const BaseMatrix & ac);
  ///
  ~CRSolver ();
  ///
  virtual void Mult (const BaseVector & v, BaseVector & prod) const;
  };



  //  CG Solver for normal problem 


  class NormalCGSolver : public CGSolver
  {
  public:
  ///
  NormalCGSolver ();
  ///
  NormalCGSolver (const BaseMatrix & aa);
  ///
  NormalCGSolver (const BaseMatrix & aa, const BaseMatrix & ac);
  ///
  ~NormalCGSolver ();
  ///
  virtual void Mult (const BaseVector & v, BaseVector & prod) const;
  };



  /// gradient method for non-symmetric, indefinite matrices
  class MinResGradientSolver : public CGSolver
  {
  public:
  ///
  MinResGradientSolver ();
  ///
  MinResGradientSolver (const BaseMatrix & aa);
  ///
  MinResGradientSolver (const BaseMatrix & aa, const BaseMatrix & ac);
  ///
  ~MinResGradientSolver ();
  ///
  virtual void Mult (const BaseVector & v, BaseVector & prod) const;
  };



  //  Bi-CG Stab solver for non-symmetric matrices

  class BiCGStabSolver : public CGSolver
  {
  public:
  ///
  BiCGStabSolver ();
  ///
  BiCGStabSolver (const BaseMatrix & aa);
  ///
  BiCGStabSolver (const BaseMatrix & aa, const BaseMatrix & ac);
  ///
  ~BiCGStabSolver ();
  ///
  virtual void Mult (const BaseVector & v, BaseVector & prod) const;
  };









  //  QMR solver for non-symmetric matrices

  class NGS_DLL_HEADER QMRSolver : public CGSolver
  {
  int status;
  ///
  const BaseMatrix * c2;
  public:
  ///
  QMRSolver ();
  ///
  QMRSolver (const BaseMatrix & aa);
  ///
  QMRSolver (const BaseMatrix & aa, const BaseMatrix & ac);
  ///
  ~QMRSolver ();
  ///
  virtual void Mult (const BaseVector & v, BaseVector & prod) const;
  };






  //  complex QMR solver for non-symmetric matrices

  class ComplexQMRSolver : public CGSolver
  {
  public:
  ///
  int status;
  ///
  const BaseMatrix * c2;
  public:
  ///
  ComplexQMRSolver ();
  ///
  ComplexQMRSolver (const BaseMatrix & aa);
  ///
  ComplexQMRSolver (const BaseMatrix & aa, const BaseMatrix & ac);
  ///
  ~ComplexQMRSolver ();
  ///
  virtual void Mult (const BaseVector & v, BaseVector & prod) const;
  };

  */


}

#endif
