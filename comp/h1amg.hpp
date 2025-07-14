#ifndef H1AMGxx_HPP_
#define H1AMGxx_HPP_

// #include <la.hpp>
#include <basematrix.hpp>
#include <sparsematrix.hpp>
#include <preconditioner.hpp>

namespace ngcomp
{
  using namespace ngla;



  struct H1AMG_Parameters
  {
    int verbose = 0;
    int smoothing_steps = 3;
    int coarsenings_per_level = 3; // number of binary coarsenings per level
    bool block_smoother = true;    // block or point smoother ?
    bool use_smoothed_prolongation = true;
    int max_coarse = 10;
    int max_level = 20;
  };

  
  template <class SCAL>
  class NGS_DLL_HEADER H1AMG_Matrix : public ngla::BaseMatrix
  {
    size_t size;
    std::shared_ptr<ngla::SparseMatrixTM<SCAL>> mat;
    std::shared_ptr<ngla::BaseBlockJacobiPrecond> smoother;
    std::shared_ptr<ngla::SparseMatrixTM<double>> prolongation, restriction;
    std::shared_ptr<ngla::BaseMatrix> coarse_precond;
    int smoothing_steps = 1;

  public:
    H1AMG_Matrix (std::shared_ptr<ngla::SparseMatrixTM<SCAL>> amat,
                  std::shared_ptr<ngcore::BitArray> freedofs,
                  ngcore::FlatArray<ngcore::IVec<2>> e2v,
                  ngcore::FlatArray<double> edge_weights,
                  ngcore::FlatArray<double> vertex_weights,
                  const H1AMG_Parameters & param,                  
                  size_t level);

    virtual int VHeight() const override { return size; }
    virtual int VWidth() const override { return size; }
    virtual bool IsComplex() const override { return is_same<SCAL,Complex>(); }
    
    virtual AutoVector CreateRowVector () const override { return mat->CreateColVector(); }
    virtual AutoVector CreateColVector () const override { return mat->CreateRowVector(); }

    virtual void Mult (const ngla::BaseVector & b, ngla::BaseVector & x) const override;
  };





  template <class SCAL>
  class H1AMG_Preconditioner : public Preconditioner
  {
    shared_ptr<BitArray> freedofs;
    shared_ptr<H1AMG_Matrix<SCAL>> mat;

    ParallelHashTable<IVec<2>,double> edge_weights_ht;
    ParallelHashTable<IVec<1>,double> vertex_weights_ht;
    
    H1AMG_Parameters param;
    
  public:

    static shared_ptr<Preconditioner> CreateBF (shared_ptr<BilinearForm> bfa, const Flags & flags, const string & name)
    {
      if (bfa->GetFESpace()->IsComplex())
        return make_shared<H1AMG_Preconditioner<Complex>> (bfa, flags, name);
      else
        return make_shared<H1AMG_Preconditioner<double>> (bfa, flags, name);
    }

    static DocInfo GetDocu ();    

    H1AMG_Preconditioner (shared_ptr<BilinearForm> abfa, const Flags & aflags,
                          const string aname = "H1AMG_cprecond")
      : Preconditioner (abfa, aflags, aname)
    {
      if (is_same<SCAL,double>::value)
        cout << IM(3) << "Create H1AMG" << endl;
      else
        cout << IM(3) << "Create H1AMG, complex" << endl;
      
      param.verbose = int(flags.GetNumFlag("verbose", 0));
      param.smoothing_steps = int(flags.GetNumFlag("smoothingsteps", 1));    
      param.use_smoothed_prolongation = flags.GetDefineFlagX("smoothedprolongation").IsMaybeTrue();
      param.max_coarse = int(flags.GetNumFlag("maxcoarse", 10));
      param.max_level = int(flags.GetNumFlag("maxlevel", 20));
    }

    virtual void InitLevel (shared_ptr<BitArray> _freedofs) override
    {
      freedofs = _freedofs;
    }

    virtual void FinalizeLevel (const BaseMatrix * matrix) override;

    virtual void AddElementMatrix (FlatArray<int> dnums,
                                   FlatMatrix<SCAL> elmat,
                                   ElementId id,
                                   LocalHeap & lh) override;

    virtual void Update () override { ; }

    virtual const BaseMatrix & GetMatrix() const override 
    {
      return *mat;
    }

  };


  
}

#endif // H1AMG_HPP_
