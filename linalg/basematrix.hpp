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
    mutable char safety_check = 0;
    bool is_complex = false;
    
  protected:
    /// 
    BaseMatrix ();
    /// 
    BaseMatrix (shared_ptr<ParallelDofs> aparalleldofs); 

  public:
    /// 
    virtual ~BaseMatrix ();
    /// virtual function must be overloaded
    virtual int VHeight() const;

    /// virtual function must be overloaded
    virtual int VWidth() const;

    /// inline function VHeight
    size_t Height() const
    {
      return VHeight();
    }
  
    /// inline function VWidth
    size_t Width() const
    {
      return VWidth();
    }

    virtual tuple<size_t, size_t> Shape() const { return { Height(), Width() }; }

    virtual xbool IsSymmetric() const { return maybe; }

    /// is matrix complex ?
    virtual bool IsComplex() const { return is_complex; }
    
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
    virtual Array<MemoryUsage> GetMemoryUsage () const;
    virtual size_t NZE () const;

    template <typename T>
      shared_ptr<T> SharedFromThis()
    { return dynamic_pointer_cast<T> (shared_from_this()); }
    /// whatever it means ... e.g. refactor sparse factorization
    virtual void Update();
    /// creates matrix of same type
    virtual shared_ptr<BaseMatrix> CreateMatrix () const;
    /// creates a matching vector, size = width
    virtual AutoVector CreateRowVector () const = 0;
    /// creates a matching vector, size = height
    virtual AutoVector CreateColVector () const = 0;
    /// creates a matching vector (for square matrices)
    [[deprecated("use CreateRowVector or CreateColVector instead")]]
    virtual AutoVector CreateVector () const;

    virtual AutoVector Evaluate(BaseVector & v) const 
    {
      auto res = CreateColVector();
      Mult (v, res);
      return res;
    }
    
    /// y = matrix * x. 
    virtual void Mult (const BaseVector & x, BaseVector & y) const;
    ///
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const;
    /// y += s matrix * x
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
    /// y += s matrix * x
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const;
  
    /// y += s Trans(matrix) * x
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const;
    /// y += s Trans(matrix) * x
    virtual void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const;
   /// y += s Trans(matrix) * x
    virtual void MultConjTransAdd (Complex s, const BaseVector & x, BaseVector & y) const;

    /// y += alpha M x
    virtual void MultAdd (FlatVector<double> alpha, const MultiVector & x, MultiVector & y) const;
    
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
    virtual optional<NgMPI_Comm> GetCommunicator() const { return nullopt; }
    
    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = nullptr) const;
    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<const Array<int>> clusters) const;
    virtual INVERSETYPE SetInverseType ( INVERSETYPE ainversetype ) const;
    virtual INVERSETYPE SetInverseType ( string ainversetype ) const;
    virtual INVERSETYPE  GetInverseType () const;
    virtual void SetInverseFlags (const Flags & flags) { ; }
    virtual shared_ptr<BaseMatrix> DeleteZeroElements(double tol) const
    {
      throw Exception (string("DeleteZeroElements not overloaded, type =")+typeid(*this).name());
    }
    
    
    virtual void DoArchive (Archive & ar);

    
    template <typename TSCAL>
      Matrix<TSCAL> ToDense() const;

    // time per run
    double Timing (int runs = 10) const;
    
    class OperatorInfo
    {
    public:
      string name = "undef";
      size_t height = 0, width = 0;
      Array<const BaseMatrix*> childs;
      OperatorInfo() = default;
      OperatorInfo(string aname, size_t ah, size_t aw)
        : name(aname), height(ah), width(aw) { } 
    };
    
    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const;
    void PrintOperatorInfo (ostream & ost, int level = 0) const;

    // base class checks for sizes, derived BlockMatrix and ParallelMatrix check more 
    virtual xbool SameShape (BaseMatrix & other) const;
    // *this * other
    virtual xbool CanCompose (BaseMatrix & other) const;
    
    virtual shared_ptr<BaseMatrix> CreateDeviceMatrix() const;
    static std::map<type_index, function<shared_ptr<BaseMatrix>(const BaseMatrix&)>> devmatcreator;
    static void RegisterDeviceMatrixCreator (type_index type,
                                             function<shared_ptr<BaseMatrix>(const BaseMatrix&)> creator)
    {
      devmatcreator[type] = creator;
    }
    
  private:
    BaseMatrix & operator= (const BaseMatrix & m2) { return *this; }

    MemoryTracer mt = { "BaseMatrix" };
  public:
    const MemoryTracer& GetMemoryTracer() const { return mt; }
  };






  /// specifies the scalar type.
  template <typename SCAL>
  class NGS_DLL_HEADER S_BaseMatrix : virtual public BaseMatrix
  {
  public:
    S_BaseMatrix () = default;
    virtual ~S_BaseMatrix () = default;
    // virtual bool IsComplex() const { return false; }
  };

  // specifies the scalar type Complex.
  template <>
  class S_BaseMatrix<Complex> : virtual public BaseMatrix
  {
  public:
    ///
    S_BaseMatrix () { is_complex = true; }
    virtual ~S_BaseMatrix () = default;
    
    // virtual bool IsComplex() const { return true; }

    /*
    /// calls MultAdd (Complex s);
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
    /// must be overloaded
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const;
  
    /// calls MultTransAdd (Complex s);
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const;
    /// should be overloaded
    virtual void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const;
    */
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



  class NGS_DLL_HEADER DynamicMatVecExpression : public DynamicBaseExpression
  {
    shared_ptr<BaseMatrix> m;
    shared_ptr<BaseVector> v;

    AutoVector CreateVector() const override
    { return m->CreateColVector(); }    

    AutoVector Evaluate() const override
    {
      return m->Evaluate(*v);
    }
    
    void AssignTo (double s, BaseVector & v2) const override
    {
      m->Mult(*v, v2);
      v2 *= s;
    }
    void AddTo (double s, BaseVector & v2) const override
    {
      m->MultAdd (s, *v, v2);
    }
    void AssignTo (Complex s, BaseVector & v2) const override
    {
      m->Mult(*v, v2);
      v2 *= s;
    }
    void AddTo (Complex s, BaseVector & v2) const override
    {
      m->MultAdd (s, *v, v2);
    }
  public:
    DynamicMatVecExpression (shared_ptr<BaseMatrix> am, shared_ptr<BaseVector> av)
      : m(am), v(av) { } 
  };

  

  /* ************************** Transpose ************************* */

  /**
     The Transpose of a BaseMatrix.
  */
  class NGS_DLL_HEADER Transpose : public BaseMatrix
  {
    const BaseMatrix & bm;
    shared_ptr<BaseMatrix> spbm;
  public:
    ///
    Transpose (const BaseMatrix & abm) : bm(abm) { ; }
    Transpose (shared_ptr<BaseMatrix> aspbm) : bm(*aspbm), spbm(aspbm) { ; }
    ///
    virtual bool IsComplex() const override { return bm.IsComplex(); }
    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override;
    
    virtual AutoVector CreateRowVector () const override { return bm.CreateColVector(); }
    virtual AutoVector CreateColVector () const override { return bm.CreateRowVector(); }

    virtual void Mult (const BaseVector & x, BaseVector & y) const override
    {
      bm.MultTrans (x, y);
    }
    
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override
    {
      bm.Mult (x, y);
    }
    
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

    auto SPtrMat() const { return spbm; }

    virtual ostream & Print (ostream & ost) const override
    {
      ost << "Transpose of " << endl;
      bm.Print(ost);
      return ost;
    }

    virtual shared_ptr<BaseMatrix> CreateDeviceMatrix() const override
    {
      return make_shared<Transpose>(bm.CreateDeviceMatrix());
    }

  };



  /* ************************** ConjTrans ************************* */

  /**
     The conjugate transpose of a BaseMatrix.
  */
  class NGS_DLL_HEADER ConjTrans : public BaseMatrix
  {
    shared_ptr<BaseMatrix> spbm;
  public:
    ConjTrans (shared_ptr<BaseMatrix> aspbm) : spbm(aspbm) { ; }
    ///
    virtual bool IsComplex() const override { return spbm->IsComplex(); }

    virtual AutoVector CreateRowVector () const override { return spbm->CreateColVector(); }
    virtual AutoVector CreateColVector () const override { return spbm->CreateRowVector(); }

    virtual void Mult (const BaseVector & x, BaseVector & y) const override
    {
      y = 0.0;
      spbm->MultConjTransAdd (1, x, y);
    }
    
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override
    {
      throw Exception("Trans of ConjTrans not available");
    }
    
    ///
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      spbm->MultConjTransAdd (s, x, y);
    }
    ///
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const override
    {
      spbm->MultConjTransAdd (s, x, y);
    }
    ///
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      throw Exception("Trans of ConjTrans not available");
    }
    ///
    virtual void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const override
    {
      throw Exception("Trans of ConjTrans not available");      
    }  

    virtual int VHeight() const override { return spbm->VWidth(); }
    virtual int VWidth() const override { return spbm->VHeight(); }


    virtual ostream & Print (ostream & ost) const override
    {
      ost << "ConjTrans of " << endl;
      spbm->Print(ost);
      return ost;
    }
  };


  

  /* ************************** Product ************************* */

  /// action of product of two matrices 
  class NGS_DLL_HEADER ProductMatrix : public BaseMatrix
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
      : bma(*aspbma), bmb(*aspbmb), spbma(aspbma), spbmb(aspbmb)
        // tempvec(aspbmb->CreateColVector())
    {
      try
        {
          tempvec.AssignPointer(bmb.CreateColVector());
        }
      catch (Exception & e)
        {
          tempvec.AssignPointer(bma.CreateRowVector());          
        }
    }
    ///
    virtual bool IsComplex() const override { return bma.IsComplex() || bmb.IsComplex(); }
    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override;

    virtual AutoVector CreateRowVector () const override { return bmb.CreateRowVector(); }
    virtual AutoVector CreateColVector () const override { return bma.CreateColVector(); }

    auto SPtrA() const { return spbma; }
    auto SPtrB() const { return spbmb; }

    
    ///
    virtual void Mult (const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("ProductMatrix::Mult"); RegionTimer reg(t);      
      bmb.Mult (x, tempvec);
      bma.Mult (tempvec, y);
    }

    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("ProductMatrix::Mult"); RegionTimer reg(t);      
      bma.MultTrans (x, tempvec);
      bmb.MultTrans (tempvec, y);
    }

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("ProductMatrix::MultAdd"); RegionTimer reg(t);      
      bmb.Mult (x, tempvec);
      bma.MultAdd (s, tempvec, y);
    }
    ///
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("ProductMatrix::MultAdd complex"); RegionTimer reg(t);            
      bmb.Mult (x, tempvec);
      bma.MultAdd (s, tempvec, y);
    }
    ///
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("ProductMatrix::MultTransAdd"); RegionTimer reg(t);            
      bma.MultTrans (x, tempvec);
      bmb.MultTransAdd (s, tempvec, y);
    }
    ///
    virtual void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("ProductMatrix::MultTransAdd complex"); RegionTimer reg(t);
      bma.MultTrans (x, tempvec);
      bmb.MultTransAdd (s, tempvec, y);
    }  

    virtual void MultAdd (FlatVector<double> alpha, const MultiVector & x, MultiVector & y) const override
    {
      static Timer t("ProductMatrix::MultAdd(mv)"); RegionTimer reg(t);
      auto tempvec = shared_ptr<BaseVector>(bmb.CreateColVector())->CreateMultiVector(x.Size());
      *tempvec = 0;
      Vector ones(x.Size());
      ones = 1.0;
      bmb.MultAdd (ones, x, *tempvec);
      bma.MultAdd (alpha, *tempvec, y);
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

    virtual shared_ptr<BaseMatrix> CreateDeviceMatrix() const override
    {
      return make_shared<ProductMatrix>(bma.CreateDeviceMatrix(), bmb.CreateDeviceMatrix());
    }

  };


  /* ************************** Sum ************************* */

  /// action of product of two matrices 
  class NGS_DLL_HEADER SumMatrix : public BaseMatrix
  {
    const BaseMatrix & bma;
    const BaseMatrix & bmb;
    shared_ptr<BaseMatrix> spbma;
    shared_ptr<BaseMatrix> spbmb;
    double a, b;
  public:
    ///
    SumMatrix (const BaseMatrix & abma, const BaseMatrix & abmb,
               double aa = 1, double ab = 1)
      : bma(abma), bmb(abmb), a(aa), b(ab)
    { ; }
    SumMatrix (shared_ptr<BaseMatrix> aspbma, shared_ptr<BaseMatrix> aspbmb,
               double aa = 1, double ab = 1);

    ///
    virtual bool IsComplex() const override { return bma.IsComplex() || bmb.IsComplex(); }

    auto SPtrA() const { return spbma; }
    auto SPtrB() const { return spbmb; }
    
    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override;

    virtual AutoVector CreateRowVector () const override
    {
      try
        {
          return bma.CreateRowVector();
        }
      catch (Exception & e)
        {
          return bmb.CreateRowVector();          
        }
    }
    virtual AutoVector CreateColVector () const override
    {
      try
        {
          return bma.CreateColVector();
        }
      catch (Exception & e)
        {
          return bmb.CreateColVector();          
        }
    }

    virtual void Mult (const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("SumMatrix::Mult"); RegionTimer reg(t);
      if (a == 1)
        bma.Mult (x, y);
      else
        {
          y = 0.0;
          bma.MultAdd (a, x, y);
        }
      bmb.MultAdd (b, x, y);
    }

    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("SumMatrix::MultTrans"); RegionTimer reg(t);
      if (a == 1)
        bma.MultTrans (x, y);
      else
        {
          y = 0.0;
          bma.MultTransAdd (a, x, y);
        }
      bmb.MultTransAdd (b, x, y);
    }
    
    ///
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("SumMatrix::MultAdd"); RegionTimer reg(t);
      bma.MultAdd (a*s, x, y);
      bmb.MultAdd (b*s, x, y);
    }
    ///
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("SumMatrix::MultAdd complex"); RegionTimer reg(t);      
      bma.MultAdd (a*s, x, y);
      bmb.MultAdd (b*s, x, y);
    }
    ///
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("SumMatrix::MultTransAdd"); RegionTimer reg(t);      
      bma.MultTransAdd (a*s, x, y);
      bmb.MultTransAdd (b*s, x, y);
    }
    ///
    virtual void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("SumMatrix::MultAdd complex"); RegionTimer reg(t);      
      bma.MultTransAdd (a*s, x, y);
      bmb.MultTransAdd (b*s, x, y);
    }  

    /// y += alpha M x
    virtual void MultAdd (FlatVector<double> alpha, const MultiVector & x, MultiVector & y) const override
    {
      static Timer t("SumMatrix::MultAdd(mv)"); RegionTimer reg(t);
      bma.MultAdd (Vector(a*alpha), x, y);
      bmb.MultAdd (Vector(b*alpha), x, y);
    }


    
    virtual int VHeight() const override
    {
      try
        {
          return bma.VHeight();
        }
      catch (Exception &)
        {
          return bmb.VHeight();
        }
    }
    
    virtual int VWidth() const override
    {
      try
        {
          return bma.VWidth();
        }
      catch (Exception &)
        {
          return bmb.VWidth();
        }
    }

    virtual ostream & Print (ostream & ost) const override
    {
      ost << "Sum of" << endl;
      ost << "Scale a = " << a << endl;
      bma.Print(ost);
      ost << "Scale b = " << b << endl;
      bmb.Print(ost);
      return ost;
    }

    virtual shared_ptr<BaseMatrix> CreateDeviceMatrix() const override
    {
      return make_shared<SumMatrix>(bma.CreateDeviceMatrix(), bmb.CreateDeviceMatrix(), a, b);
    }

  };


  /* ************************** Scale ************************* */

  template <typename TSCAL>
  class VScaleMatrix : public BaseMatrix
  {
    const BaseMatrix & bm;
    shared_ptr<BaseMatrix> spbm;
    TSCAL scale;
  public:
    ///
    VScaleMatrix (const BaseMatrix & abm, TSCAL ascale) : bm(abm), scale(ascale) { ; }
    VScaleMatrix (shared_ptr<BaseMatrix> aspbm, TSCAL ascale)
      : bm(*aspbm), spbm(aspbm), scale(ascale) { ; }
    virtual bool IsComplex() const override
    { return bm.IsComplex() || typeid(TSCAL)==typeid(Complex); } 
    ///
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("ScaleMatrix::MultAdd"); RegionTimer reg(t);
      bm.MultAdd (s*scale, x, y);
    }
    ///
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("ScaleMatrix::MultAdd complex"); RegionTimer reg(t);      
      bm.MultAdd (s*scale, x, y);
    }
    ///
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("ScaleMatrix::MultTransAdd"); RegionTimer reg(t);      
      bm.MultTransAdd (s*scale, x, y);
    }
    ///
    virtual void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("ScaleMatrix::MultTransAdd complex"); RegionTimer reg(t);      
      bm.MultTransAdd (s*scale, x, y);
    }  

    virtual void MultAdd (FlatVector<double> alpha, const MultiVector & x, MultiVector & y) const override
    {
      static Timer t("ScaleMatrix::MultAdd(mv)"); RegionTimer reg(t);
      if constexpr (is_same<TSCAL, double>())
                     {
                       bm.MultAdd (Vector(scale*alpha), x, y);
                     }
      else
        BaseMatrix::MultAdd(alpha, x, y);
    }



    
    virtual int VHeight() const override { return bm.VHeight(); }
    virtual int VWidth() const override { return bm.VWidth(); }
    virtual AutoVector CreateRowVector () const override { return bm.CreateRowVector(); }
    virtual AutoVector CreateColVector () const override { return bm.CreateColVector(); }
    virtual ostream & Print (ostream & ost) const override
    {
      ost << "Scale with " << scale << ":" << endl;
      bm.Print(ost);
      return ost;
    }

    auto SPtrMat() const { return spbm; }
    
    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override
    {
      OperatorInfo info;
      info.name = "ScaleMatrix, scale = "+ToString(scale);
      info.height = Height();
      info.width = Width();
      info.childs += &bm;
      return info;
    }
    
    virtual shared_ptr<BaseMatrix> CreateDeviceMatrix() const override
    {
      return make_shared<VScaleMatrix<TSCAL>>(bm.CreateDeviceMatrix(), scale);
    }
  };
  
  inline VScaleMatrix<double> operator* (double d, const BaseMatrix & m)
  {
    return VScaleMatrix<double> (m, d);
  }


  /* ************************** Identity ************************* */
  
  class NGS_DLL_HEADER IdentityMatrix : public BaseMatrix
  {
    bool has_format;
    size_t size;
    // bool is_complex;
  public:
    ///
    IdentityMatrix ()
      : has_format(false) { ; }
    IdentityMatrix (size_t asize, bool ais_complex)
      : has_format(true), size(asize) { is_complex=ais_complex; }
    
    // virtual bool IsComplex() const override { return is_complex; }
    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override;
    
    ///
    virtual void Mult (const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("IdentityMatrix::Mult"); RegionTimer reg(t);
      y = x;
    }
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("IdentityMatrix::MultTrans"); RegionTimer reg(t);
      y = x;
    }
    ///
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("IdentityMatrix::MultAdd"); RegionTimer reg(t);
      y += s*x;
    }
    ///
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("IdentityMatrix::MultAdd Complex"); RegionTimer reg(t);
      y += s*x;
    }
    ///
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("IdentityMatrix::MultTransAdd"); RegionTimer reg(t);
      y += s*x;
    }
    ///
    virtual void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("IdentityMatrix::MultTransAdd Complex"); RegionTimer reg(t);
      y += s*x;      
    }  

    virtual int VHeight() const override
    {
      if (has_format) return size;
      throw Exception("Identity: no Height");
    }
    virtual int VWidth() const override
    {
      if (has_format) return size;
      throw Exception("Identity: no Width");
    }
    virtual AutoVector CreateRowVector () const override
    {
      if (has_format)
        return CreateBaseVector(size, is_complex, 1);
      throw Exception("Identity: no RowVector");
    }
    virtual AutoVector CreateColVector () const override
    {
      if (has_format)
        return CreateBaseVector(size, is_complex, 1);
      throw Exception("Identity: no ColVector");
    }

    virtual ostream & Print (ostream & ost) const override
    {
      ost << "Identity" << endl;
      return ost;
    }
    
  };

  
  /* *********************** operator<< ********************** */

  // default is ProductMatrix, but optimizations for
  // ParallelMatrices
  // Embedding Matrices
  // ....
  shared_ptr<BaseMatrix> ComposeOperators (shared_ptr<BaseMatrix> a,
                                           shared_ptr<BaseMatrix> b);
  shared_ptr<BaseMatrix> AddOperators (shared_ptr<BaseMatrix> a,
                                       shared_ptr<BaseMatrix> b,
                                       double faca, double facb);

  inline shared_ptr<BaseMatrix> operator* (shared_ptr<BaseMatrix> a,
                                           shared_ptr<BaseMatrix> b)
  {
    return ComposeOperators(a,b);
  }

  inline shared_ptr<BaseMatrix> operator+ (shared_ptr<BaseMatrix> a,
                                           shared_ptr<BaseMatrix> b)
  {
    return AddOperators(a,b,1,1);
  }

  
  shared_ptr<BaseMatrix> TransposeOperator (shared_ptr<BaseMatrix> mat);
  
  /// output operator for matrices
  inline ostream & operator<< (ostream & ost, const BaseMatrix & m)
  {
    return m.Print(ost);
  }

}

#endif
