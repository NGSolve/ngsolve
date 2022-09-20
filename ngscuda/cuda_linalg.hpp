#include <cusparse.h>
#include <la.hpp>

namespace ngla
{
  class DevSparseMatrix;

  class UnifiedVector : public S_BaseVector<double>
  {
    // using int size;
    double * host_data;
    double * dev_data;
    cusparseDnVecDescr_t descr;
  
    mutable bool host_uptodate;
    mutable bool dev_uptodate;
    
  public:
    UnifiedVector (int asize);
    UnifiedVector (const BaseVector& vec);
    UnifiedVector (const UnifiedVector& vec);
    virtual ~UnifiedVector();

    void initialize_unified (size_t size);
    
    BaseVector & operator= (double d);
    BaseVector & operator= (const BaseVector & v2);
    UnifiedVector & operator= (const UnifiedVector & v2);

    template <typename T2>
    UnifiedVector & operator= (const VVecExpr<T2> & v)
    {
      BaseVector::operator= (v);
      return *this;
    }

    size_t Size () const throw();

    const double & operator [] (const int ind) const;
    double & operator [] (const int ind);

    const cusparseDnVecDescr_t& GetDescr() const;

    cusparseDnVecDescr_t& GetDescr();

    virtual BaseVector & Scale (double scal);
    virtual BaseVector & SetScalar (double scal);
    virtual BaseVector & Set (double scal, const BaseVector & v);
    virtual BaseVector & Add (double scal, const BaseVector & v);

    virtual BaseVector & operator- () const;

    double InnerProduct (const BaseVector & v2) const;

    void UpdateHost () const;
    void UpdateDevice () const;

    virtual ostream & Print (ostream & ost) const;    
    /* virtual void PrintDevice () const; */
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

/*   BaseVector& operator+ (const UnifiedVector&, const BaseVector&); */
/*   BaseVector& operator+ (const BaseVector&, const UnifiedVector&); */
/*   BaseVector& operator+ (const UnifiedVector&, const UnifiedVector&); */

/*   BaseVector& operator- (const UnifiedVector&, const BaseVector&); */
/*   BaseVector& operator- (const BaseVector&, const UnifiedVector&); */
/*   BaseVector& operator- (const UnifiedVector&, const UnifiedVector&); */

/*   BaseVector& operator* (double, const UnifiedVector&); */
/*   BaseVector& operator* (double, const UnifiedVector&); */
/*   BaseVector& operator* (const UnifiedVector& v, const DevSparseMatrix& mat); */

  AutoVector CreateUnifiedVector(size_t size);

  class DevMatrix : public BaseMatrix
  {
  public:
    DevMatrix() { }

    virtual AutoVector CreateRowVector() const { return UnifiedVector(Width()).CreateVector(); }
    virtual AutoVector CreateColVector() const { return UnifiedVector(Height()).CreateVector(); }

  };

  shared_ptr<DevMatrix> CreateDevMatrix (BaseMatrix &mat);
  shared_ptr<DevMatrix> CreateDevMatrix (Matrix<> &mat);


  class DevSparseMatrix : public DevMatrix
  {
    //cusparseMatDescr_t * descr;
    cusparseSpMatDescr_t descr;
    int * dev_ind;
    int * dev_col;
    double * dev_val;
    int height, width, nze;
  public:
    DevSparseMatrix (const SparseMatrix<double> & mat);
    ~DevSparseMatrix ();

    virtual void Mult (const BaseVector & x, BaseVector & y) const;
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;

    /* DevSparseMatrix operator- () const; */
    /* void Scale (double d); */

    virtual AutoVector CreateRowVector () const
    {
      return UnifiedVector(width).CreateVector();
    }

    virtual AutoVector CreateColVector () const
    {
      return UnifiedVector(height).CreateVector();
    }

    void Scale (double d);
    void Mult (const DevSparseMatrix& a);
  };


  /* DevSparseMatrix & operator+ (const DevSparseMatrix& a, const DevSparseMatrix& b); */
  /* DevSparseMatrix & operator- (const DevSparseMatrix& a, const DevSparseMatrix& b); */
  /* DevSparseMatrix & operator* (const DevSparseMatrix& a, const DevSparseMatrix& b); */
  /* DevSparseMatrix & operator* (const DevSparseMatrix& a, double d); */

  // dense device matrix
  class DevDMatrix : public DevMatrix
  {
  private:
    size_t height, width;

    cusparseDnMatDescr_t descr;
    double* host_data;
    double* dev_data;
  public:
    DevDMatrix ();
    DevDMatrix (const Matrix<>& mat);
    ~DevDMatrix ();

    int VHeight() const { return height; }
    int VWidth() const { return width; }

    virtual AutoVector CreateRowVector () const;
    virtual AutoVector CreateColVector () const;

    void Add (const DevDMatrix& b);
    void Mult (const DevDMatrix& b, DevDMatrix& c);
    void Scale (double d);

    void SetZero ();
    double* DevData () const;

    virtual ostream & Print (ostream & ost) const;    
  };

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
    void Scale (double d);
  };

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
