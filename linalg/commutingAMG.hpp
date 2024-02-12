#ifndef FILE_COMMUTING_AMG
#define FILE_COMMUTING_AMG

/* *************************************************************************/
/* File:   commuting_amg.hoo                                               */
/* Author: Joachim Schoeberl                                               */
/* Date:   15. Aug. 2002                                                   */
/* *************************************************************************/

namespace ngla
{

  class CommutingAMG : public BaseMatrix
  {
  protected:
    const BaseSparseMatrix * pmat;
    CommutingAMG * recAMG;

    SparseMatrixTM<double> * prol;
  
    shared_ptr<BaseSparseMatrix> coarsemat;
    shared_ptr<BaseJacobiPrecond> jacobi;
    shared_ptr<BaseBlockJacobiPrecond> bjacobi;
    shared_ptr<BaseMatrix> inv;

  public:
    CommutingAMG ()
    { ; }
    virtual ~CommutingAMG () 
    { ; }

    virtual bool IsComplex() const override { return false; }

    virtual void ComputeMatrices (const BaseSparseMatrix & mat) = 0;
    virtual void Mult (const BaseVector & x, BaseVector & y) const override = 0;


    virtual int VHeight() const override { return pmat->Height(); }
    virtual int VWidth() const override { return pmat->Width(); }

    virtual size_t NZE() const override = 0;
    virtual AutoVector CreateVector () const override
    {
      return pmat->CreateColVector();
    }
    virtual AutoVector CreateColVector () const override
    {
      return pmat->CreateRowVector();
    }
    virtual AutoVector CreateRowVector () const override 
    {
      return pmat->CreateColVector();
    }
  };


  class AMG_H1 : public CommutingAMG
  {
  public:

    AMG_H1 (const BaseMatrix & sysmat,
	    Array<ngstd::IVec<2> > & e2v,
	    Array<double> & weighte,
	    int levels);

    virtual ~AMG_H1 ();

    virtual void ComputeMatrices (const BaseSparseMatrix & mat);
    virtual size_t NZE() const;

    virtual void Mult (const BaseVector & x, BaseVector & y) const;
  };




  class AMG_HCurl : public CommutingAMG
  {
    SparseMatrixTM<double> * grad;

    shared_ptr<BaseSparseMatrix> h1mat;
    AMG_H1 * h1AMG;

  public:

    AMG_HCurl (const BaseMatrix & sysmat,
	       const Array<Vec<3> > & vertices,
	       Array<ngstd::IVec<2> > & e2v,
	       Array<ngstd::IVec<4> > & f2v,
	       Array<double> & weighte,
	       Array<double> & weightf,
	       int levels);

    virtual ~AMG_HCurl ();

    virtual void ComputeMatrices (const BaseSparseMatrix & mat);
    virtual size_t NZE() const;

    virtual void Mult (const BaseVector & x, BaseVector & y) const;
  };


}


#endif
