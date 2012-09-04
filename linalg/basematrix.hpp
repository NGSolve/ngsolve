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
  enum INVERSETYPE { PARDISO, PARDISOSPD, SPARSECHOLESKY, SUPERLU, SUPERLU_DIST, MUMPS, MASTERINVERSE };


  /**
     The base for all matrices in the linalg.
  */
  class NGS_DLL_HEADER BaseMatrix
  {
  protected:
    const ParallelDofs * paralleldofs;

  public:
    /// 
    BaseMatrix ();
    /// 
    // BaseMatrix (const BaseMatrix & amat);
    //
    BaseMatrix (ParallelDofs * aparalleldofs); 
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


    void SetParallelDofs (const ParallelDofs * pardofs) { paralleldofs = pardofs; }
    const ParallelDofs * GetParallelDofs () const { return paralleldofs; }

    virtual BaseMatrix * InverseMatrix (const BitArray * subset = 0) const;
    virtual BaseMatrix * InverseMatrix (const Array<int> * clusters) const;
    virtual INVERSETYPE SetInverseType ( INVERSETYPE ainversetype ) const;
    virtual INVERSETYPE SetInverseType ( string ainversetype ) const;
    virtual INVERSETYPE  GetInverseType () const;
    
    
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



  class VScaleMatrix : public BaseMatrix
  {
    const BaseMatrix & bm;
    double scale;
  public:
    ///
    VScaleMatrix (const BaseMatrix & abm, double ascale) : bm(abm), scale(ascale) { ; }
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
