#include <cusparse.h>
#include <la.hpp>

namespace ngla
{

	/* TODO:
	 *  -) always create cusparseDnVecDescr_t?
	 *   .) alternative: static, similiar to handle
	 *   .) alternative: create derived class "CuspUnifiedVector"?
	 *  -) keep cusparseSpMV_bufferSize? (similiar to handle)
	 *   .) alternative: add to devMatrix or create MultMatVec class
	 *  -) why was there no destructor? (cudaFree!)
	 *  */
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
		UnifiedVector (BaseVector& vec);
		virtual ~UnifiedVector();
    
    BaseVector & operator= (double d);
    BaseVector & operator= (BaseVector & v2);

    template <typename T2>
    UnifiedVector & operator= (const VVecExpr<T2> & v)
    {
      BaseVector::operator= (v);
      return *this;
    }
    
		/*
		 *	TODO:
		 *		avoid cpy between device and host
		 *		maybe leave the update to the user?
		 * */
		const double & operator [] (const int ind) const
		{
			UpdateHost(); 
			return host_data[ind];
		}

		double & operator [] (const int ind)
		{
			UpdateHost(); 
			dev_uptodate = false; // TODO: not sure, that this is the best approach
			return host_data[ind];
		}

		int Size () const
		{
			return size;
		}

		const cusparseDnVecDescr_t& GetDescr() const
		{
			return descr;
		}

		cusparseDnVecDescr_t& GetDescr()
		{
			return descr;
		}

    virtual BaseVector & Scale (double scal);
    virtual BaseVector & SetScalar (double scal);
    virtual BaseVector & Set (double scal, const BaseVector & v);
    virtual BaseVector & Add (double scal, const BaseVector & v);

    virtual double InnerProduct (const BaseVector & v2) const;


    void UpdateHost () const;
    void UpdateDevice () const;


    virtual ostream & Print (ostream & ost) const;    
		virtual void PrintDevice () const;
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

		// TODO: check this. just placeholders 
		virtual AutoVector CreateRowVector () const
		{
			return UnifiedVector(width).CreateVector();
		}

		virtual AutoVector CreateColVector () const
		{
			return UnifiedVector(height).CreateVector();
		}
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
