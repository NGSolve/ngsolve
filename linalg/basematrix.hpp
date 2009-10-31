#ifndef FILE_NGS_BASEMATRIX
#define FILE_NGS_BASEMATRIX


/*********************************************************************/
/* File:   basematrix.hpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

namespace ngla
{


  /**
     The base for all matrices in the linalg.
  */
  class BaseMatrix
  {
  protected:

  public:
    /// constructor
    BaseMatrix ();
    /// copy-constructor
    BaseMatrix ( const BaseMatrix & amat );
    //
    //  BaseMatrix ( const ngparallel::ParallelDofs * aparalleldofs ); 
    /// destructor
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
  
    virtual ostream & Print (ostream & ost) const;
    virtual void MemoryUsage (Array<MemoryUsageStruct*> & mu) const;

    // virtual const void * Data() const;
    // virtual void * Data();

    /// creates matrix of same type
    virtual BaseMatrix * CreateMatrix () const;
    /// creates matrix of same type
    virtual BaseMatrix * CreateMatrix (const Array<int> & elsperrow) const;
    /// creates a matching vector, size = width
    virtual BaseVector * CreateRowVector () const;
    /// creates a matching vector, size = height
    virtual BaseVector * CreateColVector () const;
    /// creates a matching vector (for square matrices)
    virtual BaseVector * CreateVector () const;

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




    // to split mat x vec for symmetric matrices
    virtual void MultAdd1 (double s, const BaseVector & x, BaseVector & y) const
    {
      MultAdd (s, x, y);
    }

    virtual void MultAdd2 (double s, const BaseVector & x, BaseVector & y) const
    {
      ;
    }

    /*
    // parallel methods --> ParallelBaseMatrix

    virtual BaseMatrix * ConsistentMat () { cerr << "ERROR -- BaseMatrix::ConsistentMat() called" << endl; return 0; }
    virtual const BaseMatrix * ConsistentMat () const  { cerr << "ERROR -- BaseMatrix::ConsistentMat() called" << endl; return 0; }

    virtual void SetConsistentMat ( BaseMatrix * aconsistentmat )
    { cerr << "ERROR -- BaseMatrix::SetConsistentMat called" << endl; }

    virtual void AllocateConsistentMat ()
    { cerr << "ERROR -- BaseMatrix::AllocateConsistentMat called" << endl; }

    virtual void  AllocateConsistentMat ( const ngla::MatrixGraph & graph )
    { cerr << "ERROR -- BaseMatrix::AllocateConsistentMat(const MatrixGraph&) called" << endl; }

    virtual void CalcConsistentMat () 
    { cerr << "ERROR -- BaseMatrix::CalcConsistentMat called" << endl; }

    virtual void SetParallelDofs ( ngparallel::ParallelDofs * aparalleldofs );

    virtual const ngparallel::ParallelDofs * GetParallelDofs ( ) const { return paralleldofs; }

    virtual bool IsParallelMatrix() const
    {
    (*testout) << "PARALLELDOFS " <<  paralleldofs << endl;
    if ( paralleldofs ) return true;
    else return false;
    }
    */
  };

  class ParallelBaseMatrix : virtual public BaseMatrix
  {
  protected:
    const ngparallel::ParallelDofs * paralleldofs;

  public:

    /// constructor
    ParallelBaseMatrix ();
    /// copy-constructor
    ParallelBaseMatrix ( const ParallelBaseMatrix & amat );
    ///
    ParallelBaseMatrix ( const ngparallel::ParallelDofs * aparalleldofs );
    /// destructor
    virtual ~ParallelBaseMatrix ();

    virtual BaseMatrix * ConsistentMat () { cerr << "ERROR -- ParallelBaseMatrix::ConsistentMat() called" << endl; return 0; }
    virtual const BaseMatrix * ConsistentMat () const  
    { cerr << "ERROR -- ParallelBaseMatrix::ConsistentMat() called" << endl; return 0; }

    virtual void SetConsistentMat ( BaseMatrix * aconsistentmat )
    { cerr << "ERROR -- ParallelBaseMatrix::SetConsistentMat called" << endl; }

    virtual void AllocateConsistentMat ()
    { cerr << "ERROR -- ParallelBaseMatrix::AllocateConsistentMat called" << endl; }

    virtual void  AllocateConsistentMat ( const class MatrixGraph & graph )
    { cerr << "ERROR -- ParallelBaseMatrix::AllocateConsistentMat(const MatrixGraph&) called" << endl; }

    virtual void CalcConsistentMat (LocalHeap & lh) 
    { cerr << "ERROR -- ParallelBaseMatrix::CalcConsistentMat called" << endl; }

    virtual void SetParallelDofs ( ngparallel::ParallelDofs * aparalleldofs );

    virtual const ngparallel::ParallelDofs * GetParallelDofs ( ) const { return paralleldofs; }

    virtual bool IsParallelMatrix() const
    {
      (*testout) << "PARALLELDOFS " <<  paralleldofs << endl;
      if ( paralleldofs ) return true;
      else return false;
    }

  };


  /// specifies the scalar type.
  template <typename SCAL>
  class S_BaseMatrix :  virtual public BaseMatrix
  {
  public:
    ///
    S_BaseMatrix ();
    //
    //S_BaseMatrix (const ngparallel::ParallelDofs * aparallelsdof);
    ///
    virtual ~S_BaseMatrix ();
  };


  // specifies the scalar type Complex.
  template <>
  class S_BaseMatrix<Complex> : virtual public BaseMatrix
  {
  public:
    ///
    S_BaseMatrix ();
    // S_BaseMatrix (const ngparallel::ParallelDofs * aparallelsdof);
    ///
    virtual ~S_BaseMatrix ();

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
#ifdef DEBUG
      if (m.Height() != v.Size() || m.Width() != x.Size())
	throw Exception ("matrix-vector: size does not fit");
#endif
      m.Mult (x, v);
      v *= s;
    }

    template <class TS>
    void AddTo (TS s, BaseVector & v) const
    { 
#ifdef DEBUG
      if (m.Height() != v.Size() || m.Width() != x.Size())
	throw Exception ("matrix-vector MultAdd: size does not fit");
#endif
      m.MultAdd (s, x, v);
    }
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
  public:
    ///
    Transpose (const BaseMatrix & abm) : bm(abm) { ; }
    ///
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const
    {
      bm.MultTransAdd (s, x, y);
    }
    ///
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const 
    {
      bm.MultTransAdd (s, x, y);
    }
    ///
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
    {
      bm.MultAdd (s, x, y);
    }
    ///
    virtual void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const
    {
      bm.MultAdd (s, x, y);
    }  
  };

  /* *********************** operator<< ********************** */

  /// output operator for matrices
  inline ostream & operator<< (ostream & ost, const BaseMatrix & m)
  {
    return m.Print(ost);
  }

}

#endif
