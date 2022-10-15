#include <cusparse.h>
#include <la.hpp>

namespace ngla
{
  void InitCuLinalg();

  class DevSparseMatrix;

  class UnifiedVector : public S_BaseVector<double>
  {
    /* using int size; */
    double * host_data;
    double * dev_data;
    cusparseDnVecDescr_t descr;
  
    mutable bool host_uptodate;
    mutable bool dev_uptodate;
    
  public:
    UnifiedVector (int asize);
    UnifiedVector (const BaseVector& vec);

    virtual ~UnifiedVector();

    BaseVector & operator= (double d);
    BaseVector & operator= (const BaseVector & v2);

    template <typename T2>
    UnifiedVector & operator= (const VVecExpr<T2> & v)
    {
      BaseVector::operator= (v);
      return *this;
    }

    const double & operator [] (const int ind) const;
    double & operator [] (const int ind);

    const cusparseDnVecDescr_t& GetDescr() const;

    cusparseDnVecDescr_t& GetDescr();

    virtual BaseVector & Scale (double scal);
    virtual BaseVector & SetScalar (double scal);
    virtual BaseVector & Set (double scal, const BaseVector & v);
    virtual BaseVector & Add (double scal, const BaseVector & v);

    double InnerProduct (const BaseVector & v2, bool conjugate=false) const;

    void UpdateHost () const;
    void UpdateDevice () const;

    virtual ostream & Print (ostream & ost) const;    
    virtual ostream & PrintStatus (ostream & ost) const;
    /* virtual void PrintDevice () const; */
    virtual AutoVector CreateVector () const;

    virtual FlatVector<double> FVDouble () const;
    virtual FlatVector<Complex> FVComplex () const;
    virtual void * Memory() const throw ();

    virtual void GetIndirect (const FlatArray<int> & ind, 
            const FlatVector<double> & v) const;
    virtual void GetIndirect (const FlatArray<int> & ind, 
            const FlatVector<Complex> & v) const;

    virtual double* DevData() const
    { 
      return dev_data; 
    }
    
    friend class DevDMatrix;
    friend class DevSparseMatrix;
    friend class DevJacobiPrecond;
  };



  /* AutoVector CreateUnifiedVector(size_t size); */

  class DevMatrix : public BaseMatrix
  {
  public:
    DevMatrix() { }

    virtual AutoVector CreateRowVector() const { return make_unique<UnifiedVector>(Width()); }
    virtual AutoVector CreateColVector() const { return make_unique<UnifiedVector>(Height()); }
  };

  shared_ptr<BaseMatrix> CreateDevMatrix (BaseMatrix &mat);
  shared_ptr<BaseMatrix> CreateDevMatrix (Matrix<> &mat);


  class DevSparseMatrix : public DevMatrix
  {
  protected:
    //cusparseMatDescr_t * descr;
    cusparseSpMatDescr_t descr;
    int * dev_ind;
    int * dev_col;
    double * dev_val;
    int height, width, nze;
  public:
    DevSparseMatrix () { }
    DevSparseMatrix (const SparseMatrix<double> & mat);
    virtual ~DevSparseMatrix ();

    // TODO: implement?
    /* virtual shared_ptr<BaseMatrix> InverseMatrix */ 
    /* virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override; */
    /* // y += s L * x */
    /* virtual void MultAdd1 (const BaseVector & x, BaseVector & y, */ 
    /*                        const BitArray * ainner = NULL, */
    /*                        const Array<int> * acluster = NULL) const override; */

    /* // y += s (D + L^T) * x) */
    /* virtual void MultAdd2 (double s, const BaseVector & x, BaseVector & y, */
    /*                        const BitArray * ainner = NULL, */
    /*                        const Array<int> * acluster = NULL) const override; */


    virtual void Mult (const BaseVector & x, BaseVector & y) const;
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
    /* virtual void Scale (double d); */

    virtual AutoVector CreateRowVector () const
    {
      return UnifiedVector(width).CreateVector();
    }

    virtual AutoVector CreateColVector () const
    {
      return UnifiedVector(height).CreateVector();
    }

    virtual int VHeight() const { return height; }
    virtual int VWidth() const { return width; }
  };

  shared_ptr<DevSparseMatrix> MatMult (const DevSparseMatrix& mata, const DevSparseMatrix& matb);

  // dense device matrix
  class DevDMatrix : public DevMatrix
  {
  private:
    size_t height, width;
    double* dev_data;
    /* cusparseDnMatDescr_t descr; */

  public:
    DevDMatrix ();
    DevDMatrix (size_t height, size_t width);
    DevDMatrix (const Matrix<>& mat);
    DevDMatrix (const DevDMatrix& mat);
    ~DevDMatrix ();

    /* const DevDMatrix & operator= (const Expr<TBxx> & mat) const; */
    const DevDMatrix & operator= (double d) const;
    const DevDMatrix & operator= (const DevDMatrix & mat) const;

    int VHeight() const { return height; }
    int VWidth() const { return width; }

    virtual AutoVector CreateRowVector () const;
    virtual AutoVector CreateColVector () const;

    virtual void Add (const BaseMatrix& b);
    virtual void Scale (double d);
    virtual void Mult (const BaseVector& x, BaseVector& y) const;
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;

    void SetZero ();

    double* DevData () const;

    virtual ostream & Print (ostream & ost) const;    
  };

  shared_ptr<DevDMatrix> MatMult (const DevDMatrix& mata, const DevDMatrix& matb);

  class DevEBEMatrix : public DevMatrix
  {
  private:
    size_t height, width;
    DevDMatrix devmat;
    Table<int> col_dnums, row_dnums;
    bool disjointrows, disjointcols;
  public:
    DevEBEMatrix (const ConstantElementByElementMatrix& mat);
    ~DevEBEMatrix ();

    virtual AutoVector CreateRowVector () const;
    virtual AutoVector CreateColVector () const;

    void MultAdd (double s, const UnifiedVector & x, UnifiedVector & y);
    /* void Scale (double d); */
  };

  // currently uses Mult and MultAdd of DevSparseMatrix
  class DevJacobiPrecond : public DevSparseMatrix
  {
  private:
    /* cusparseSpMatDescr_t descr; */
    shared_ptr<BitArray> inner;
    double* dev_invdiag;
  public:
    DevJacobiPrecond (const SparseMatrix<double> & amat, 
      shared_ptr<BitArray> ainner=nullptr, bool use_par=true);

    /* DevJacobiPrecond (const JacobiPrecond<double> & amat); */

    virtual ~DevJacobiPrecond ();

    /* void Mult (const BaseVector & x, BaseVector & y) const; */
    /* void MultAdd (double s, const BaseVector & x, BaseVector & y) const; */

    /* int VHeight() const override { return height; } */
    /* int VWidth() const override { return height; } */

  };


  // old version
  /* class DevJacobiPreconditioner : public BaseMatrix */
  /* { */
  /*   // should be like this: */
  /*   // double * dev_diag; */
  /*   // int size; */

  /*   // stored as sparse matrix ... */
  /*   /1* cusparseMatDescr_t * descr; *1/ */

    /* // used to be cusparseMatSpDescr_t... not sure about the difference */
    /* cusparseSpMatDescr_t descr; */
  /*   int * dev_ind; */
  /*   int * dev_col; */
  /*   double * dev_val; */
  /*   int height, width, nze; */
    

  /* public: */
  /*   DevJacobiPreconditioner (const SparseMatrix<double> & mat, const BitArray & freedofs); */
  /*   virtual void Mult (const BaseVector & x, BaseVector & y) const; */
  /*   virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const; */
  /* }; */

}
