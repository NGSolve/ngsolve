#ifdef CUDA


#include <cusparse.h>

namespace ngla
{

  class UnifiedVector : public S_BaseVector<double>
  {
    // using int size;
    double * host_data;
    double * dev_data;
    mutable bool host_uptodate;
    mutable bool dev_uptodate;
    
  public:
    UnifiedVector (int asize);
    
    BaseVector & operator= (double d);
    BaseVector & operator= (BaseVector & v2);

    template <typename T2>
    UnifiedVector & operator= (const VVecExpr<T2> & v)
    {
      BaseVector::operator= (v);
      return *this;
    }
    

    virtual BaseVector & Scale (double scal);
    virtual BaseVector & SetScalar (double scal);
    virtual BaseVector & Set (double scal, const BaseVector & v);
    virtual BaseVector & Add (double scal, const BaseVector & v);

    virtual double InnerProduct (const BaseVector & v2) const;


    void UpdateHost () const;
    void UpdateDevice () const;


    virtual ostream & Print (ostream & ost) const;    
    virtual AutoVector CreateVector () const;

    virtual FlatVector<double> FVDouble () const;
    virtual FlatVector<Complex> FVComplex () const;
    virtual void * Memory() const throw ();


    virtual void GetIndirect (const FlatArray<int> & ind, 
			      const FlatVector<double> & v) const;
    virtual void GetIndirect (const FlatArray<int> & ind, 
			      const FlatVector<Complex> & v) const;

    
    friend class DevSparseMatrix;
    friend class DevJacobiPreconditioner;
  };

  class DevSparseMatrix : public BaseMatrix
  {
    cusparseMatDescr_t * descr;
    int * dev_ind;
    int * dev_col;
    double * dev_val;
    int height, width, nze;
  public:
    DevSparseMatrix (const SparseMatrix<double> & mat);
    virtual void Mult (const BaseVector & x, BaseVector & y) const;
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
  };


  class DevJacobiPreconditioner : public BaseMatrix
  {
    // should be like this:
    // double * dev_diag;
    // int size;

    // stored as sparse matrix ...
    cusparseMatDescr_t * descr;
    int * dev_ind;
    int * dev_col;
    double * dev_val;
    int height, width, nze;
    

  public:
    DevJacobiPreconditioner (const SparseMatrix<double> & mat, const BitArray & freedofs);
    virtual void Mult (const BaseVector & x, BaseVector & y) const;
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
  };

}

#endif
