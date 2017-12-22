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

    virtual bool IsComplex() const { return false; }

    virtual void ComputeMatrices (const BaseSparseMatrix & mat) = 0;
    virtual void Mult (const BaseVector & x, BaseVector & y) const = 0;


    virtual int VHeight() const { return pmat->Height(); }
    virtual int VWidth() const { return pmat->Width(); }

    virtual int NZE() const = 0;
    virtual AutoVector CreateVector () const
    {
      return pmat->CreateVector();
    }
  };


  class AMG_H1 : public CommutingAMG
  {
  public:

    AMG_H1 (const BaseMatrix & sysmat,
	    Array<ngstd::INT<2> > & e2v,
	    Array<double> & weighte,
	    int levels);

    virtual ~AMG_H1 ();

    virtual void ComputeMatrices (const BaseSparseMatrix & mat);
    virtual int NZE() const;

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
	       Array<ngstd::INT<2> > & e2v,
	       Array<ngstd::INT<4> > & f2v,
	       Array<double> & weighte,
	       Array<double> & weightf,
	       int levels);

    virtual ~AMG_HCurl ();

    virtual void ComputeMatrices (const BaseSparseMatrix & mat);
    virtual int NZE() const;

    virtual void Mult (const BaseVector & x, BaseVector & y) const;
  };


}


#endif
