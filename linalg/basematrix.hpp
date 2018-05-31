#ifndef FILE_NGS_BASEMATRIX
#define FILE_NGS_BASEMATRIX


/*********************************************************************/
/* File:   basematrix.hpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

namespace ngla
{


  // sets the solver which is used for InverseMatrix
  enum INVERSETYPE { PARDISO, PARDISOSPD, SPARSECHOLESKY, SUPERLU, SUPERLU_DIST, MUMPS, MASTERINVERSE, UMFPACK };
  extern string GetInverseName (INVERSETYPE type);

  /**
     The base for all matrices in the linalg.
  */
  class NGS_DLL_HEADER BaseMatrix : public enable_shared_from_this_virtual<BaseMatrix>
  {
  protected:
    shared_ptr<ParallelDofs> paralleldofs;

  protected:
    /// 
    BaseMatrix ();
    /// 
    // BaseMatrix (const BaseMatrix & amat);
    //
    BaseMatrix (shared_ptr<ParallelDofs> aparalleldofs); 

  public:
    /// 
    virtual ~BaseMatrix ();
    /// virtual function must be overloaded
    virtual int VHeight() const;

    /// virtual function must be overloaded
    virtual int VWidth() const;

    /// inline function VHeight
    int Height() const
    {
      return VHeight();
    }
  
    /// inline function VWidth
    int Width() const
    {
      return VWidth();
    }

    /// is matrix complex ?
    virtual bool IsComplex() const { return false; }
    
    /// scalar assignment
    BaseMatrix & operator= (double s)
    {
      AsVector().SetScalar(s);
      return *this;
    }

    /// linear access of matrix memory
    virtual BaseVector & AsVector();
    /// linear access of matrix memory
    virtual const BaseVector & AsVector() const;
    ///
    virtual void SetZero();

    virtual ostream & Print (ostream & ost) const;
    virtual void MemoryUsage (Array<MemoryUsageStruct*> & mu) const;
    virtual size_t NZE () const;
    // virtual const void * Data() const;
    // virtual void * Data();
    
    template <typename T>
      shared_ptr<T> SharedFromThis()
    { return dynamic_pointer_cast<T> (shared_from_this()); }
    /// whatever it means ... e.g. refactor sparse factorization
    virtual void Update() { ; } 
    /// creates matrix of same type
    virtual shared_ptr<BaseMatrix> CreateMatrix () const;
    /// creates matrix of same type
    // virtual BaseMatrix * CreateMatrix (const Array<int> & elsperrow) const;
    /// creates a matching vector, size = width
    virtual AutoVector CreateRowVector () const;
    /// creates a matching vector, size = height
    virtual AutoVector CreateColVector () const;
    /// creates a matching vector (for square matrices)
    virtual AutoVector CreateVector () const;

    /// y = matrix * x. Multadd should be implemented, instead
    virtual void Mult (const BaseVector & x, BaseVector & y) const;
    /// y += s matrix * x
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
    /// y += s matrix * x
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const;
  
    /// y += s Trans(matrix) * x
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const;
    /// y += s Trans(matrix) * x
    virtual void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const;




    /**
       to split mat x vec for symmetric matrices
       only rows with inner or cluster true need by added (but more can be ...)
    */
    virtual void MultAdd1 (double s, const BaseVector & x, BaseVector & y,
			   const BitArray * ainner = NULL,
			   const Array<int> * acluster = NULL) const;

    /// only cols with inner or cluster true need by added (but more can be ...)
    virtual void MultAdd2 (double s, const BaseVector & x, BaseVector & y,
			   const BitArray * ainner = NULL,
			   const Array<int> * acluster = NULL) const;


    void SetParallelDofs (shared_ptr<ParallelDofs> pardofs) { paralleldofs = pardofs; }
    shared_ptr<ParallelDofs> GetParallelDofs () const { return paralleldofs; }

    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = nullptr) const;
    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<const Array<int>> clusters) const;
    virtual INVERSETYPE SetInverseType ( INVERSETYPE ainversetype ) const;
    virtual INVERSETYPE SetInverseType ( string ainversetype ) const;
    virtual INVERSETYPE  GetInverseType () const;

    virtual void DoArchive (Archive & ar);
    
  private:
    BaseMatrix & operator= (const BaseMatrix & m2) { return *this; }
  };






  /// specifies the scalar type.
  template <typename SCAL>
  class NGS_DLL_HEADER S_BaseMatrix : virtual public BaseMatrix
  {
  public:
    ///
    S_BaseMatrix ();
    ///
    virtual ~S_BaseMatrix ();

    virtual bool IsComplex() const { return false; }
  };

  // specifies the scalar type Complex.
  template <>
  class S_BaseMatrix<Complex> : virtual public BaseMatrix
  {
  public:
    ///
    S_BaseMatrix ();
    ///
    virtual ~S_BaseMatrix ();
    virtual bool IsComplex() const { return true; }
    
    /// calls MultAdd (Complex s);
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
    /// must be overloaded
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const;
  
    /// calls MultTransAdd (Complex s);
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const;
    /// should be overloaded
    virtual void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const;
  };







  /* *************************** Matrix * Vector ******************** */


  /// 
  class VMatVecExpr
  {
    const BaseMatrix & m;
    const BaseVector & x;
  
  public:
    VMatVecExpr (const BaseMatrix & am, const BaseVector & ax) : m(am), x(ax) { ; }

    template <class TS>
    void AssignTo (TS s, BaseVector & v) const
    { 
      CheckSize (v);
      /*
      if (m.Height() != v.Size() || m.Width() != x.Size())
	throw Exception (ToString ("matrix-vector: size does not fit\n") +
                         "matrix-type = " + typeid(m).name() +
			 "Matrix:     " + ToString(m.Height()) + " x " + ToString(m.Width()) + "\n"
			 "Vector in : " + ToString(x.Size()) + "\n"
			 "Vector res: " + ToString(v.Size()));
      */
      m.Mult (x, v);
      v *= s;
    }

    template <class TS>
    void AddTo (TS s, BaseVector & v) const
    { 
      CheckSize (v);
      /*
      if (m.Height() != v.Size() || m.Width() != x.Size())
	throw Exception ("matrix-vector MultAdd: size does not fit");
      */
      m.MultAdd (s, x, v);
    }

    NGS_DLL_HEADER void CheckSize (BaseVector & dest_vec) const;
  };


  /// BaseMatrix times Vector - expression template
  inline VVecExpr<VMatVecExpr>
  operator* (const BaseMatrix & a, const BaseVector & b)
  {
    return VMatVecExpr (a, b);
  }


  /* ************************** Transpose ************************* */

  /**
     The Transpose of a BaseMatrix.
  */
  class Transpose : public BaseMatrix
  {
    const BaseMatrix & bm;
    shared_ptr<BaseMatrix> spbm;
  public:
    ///
    Transpose (const BaseMatrix & abm) : bm(abm) { ; }
    Transpose (shared_ptr<BaseMatrix> aspbm) : bm(*aspbm), spbm(aspbm) { ; }
    ///
    virtual bool IsComplex() const override { return bm.IsComplex(); }

    virtual AutoVector CreateRowVector () const override { return bm.CreateColVector(); }
    virtual AutoVector CreateColVector () const override { return bm.CreateRowVector(); }
    
    ///
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      bm.MultTransAdd (s, x, y);
    }
    ///
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const override
    {
      bm.MultTransAdd (s, x, y);
    }
    ///
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      bm.MultAdd (s, x, y);
    }
    ///
    virtual void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const override
    {
      bm.MultAdd (s, x, y);
    }  

    virtual int VHeight() const override { return bm.VWidth(); }
    virtual int VWidth() const override { return bm.VHeight(); }


    virtual ostream & Print (ostream & ost) const override
    {
      ost << "Transpose of " << endl;
      bm.Print(ost);
      return ost;
    }
  };



  /* ************************** Product ************************* */

  /// action of product of two matrices 
  class ProductMatrix : public BaseMatrix
  {
    const BaseMatrix & bma;
    const BaseMatrix & bmb;
    shared_ptr<BaseMatrix> spbma;
    shared_ptr<BaseMatrix> spbmb;
    mutable AutoVector tempvec;
  public:
    ///
    ProductMatrix (const BaseMatrix & abma, const BaseMatrix & abmb)
      : bma(abma), bmb(abmb), tempvec(abmb.CreateColVector())
    { ; }
    ProductMatrix (shared_ptr<BaseMatrix> aspbma, shared_ptr<BaseMatrix> aspbmb)
      : bma(*aspbma), bmb(*aspbmb), spbma(aspbma), spbmb(aspbmb),
        tempvec(aspbmb->CreateColVector())
    { ; }
    ///
    virtual bool IsComplex() const override { return bma.IsComplex() || bmb.IsComplex(); }

    virtual AutoVector CreateRowVector () const override { return bmb.CreateRowVector(); }
    virtual AutoVector CreateColVector () const override { return bma.CreateColVector(); }
    
    ///
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      bmb.Mult (x, tempvec);
      bma.MultAdd (s, tempvec, y);
    }
    ///
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const override
    {
      bmb.Mult (x, tempvec);
      bma.MultAdd (s, tempvec, y);
    }
    ///
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      tempvec = 0.0;
      bma.MultTransAdd (1, x, tempvec);
      bmb.MultTransAdd (s, tempvec, y);
    }
    ///
    virtual void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const override
    {
      tempvec = 0.0;      
      bma.MultTransAdd (1, x, tempvec);
      bmb.MultTransAdd (s, tempvec, y);
    }  

    virtual int VHeight() const override { return bma.VHeight(); }
    virtual int VWidth() const override { return bmb.VWidth(); }

    virtual ostream & Print (ostream & ost) const override
    {
      ost << "Product of" << endl;
      bma.Print(ost);
      bmb.Print(ost);
      return ost;
    }
  };




  
  class VScaleMatrix : public BaseMatrix
  {
    const BaseMatrix & bm;
    double scale;
  public:
    ///
    VScaleMatrix (const BaseMatrix & abm, double ascale) : bm(abm), scale(ascale) { ; }
    virtual bool IsComplex() const { return bm.IsComplex(); } 
    ///
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const
    {
      bm.MultAdd (s*scale, x, y);
    }
    ///
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const 
    {
      bm.MultAdd (s*scale, x, y);
    }
    ///
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
    {
      bm.MultTransAdd (s*scale, x, y);
    }
    ///
    virtual void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const
    {
      bm.MultTransAdd (s*scale, x, y);
    }  

    virtual int VHeight() const { return bm.VHeight(); }
    virtual int VWidth() const { return bm.VWidth(); }
  };

  inline VScaleMatrix operator* (double d, const BaseMatrix & m)
  {
    return VScaleMatrix (m, d);
  }

  /* *********************** operator<< ********************** */

  /// output operator for matrices
  inline ostream & operator<< (ostream & ost, const BaseMatrix & m)
  {
    return m.Print(ost);
  }

}

#endif
