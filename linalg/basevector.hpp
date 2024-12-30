#ifndef FILE_BASEVECTOR
#define FILE_BASEVECTOR

/*********************************************************************/
/* File:   basevector.hpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   7. Feb. 2003                                              */
/*********************************************************************/

#include <bla.hpp>
#include <core/mpi_wrapper.hpp>
// #include "paralleldofs.hpp"

namespace ngla
{
  using namespace ngbla;
  
  class BaseVector;
  class AutoVector;
  class MultiVector;

  class ParallelDofs;
  
  class DofRange : public T_Range<size_t>
  {
    shared_ptr<ParallelDofs> pardofs;
  public:
    DofRange () { }
    DofRange (T_Range<size_t> range, shared_ptr<ParallelDofs> apardofs)
      : T_Range<size_t>(range), pardofs(apardofs) { ; }
    DofRange (size_t end, shared_ptr<ParallelDofs> apardofs)
      : T_Range<size_t>(0, end), pardofs(apardofs) { ; }
    shared_ptr<ParallelDofs> GetParallelDofs() const { return pardofs; }
  };
  

  
  
  template <class SCAL> class S_BaseVector;

  class NGS_DLL_HEADER ComplexConjugate;
  class NGS_DLL_HEADER ComplexConjugate2;

  template<class IPTYPE>
  class SCAL_TRAIT
  {
  public:
    typedef double SCAL;
  };

  template<> class SCAL_TRAIT<Complex>
  {
  public:
    typedef Complex SCAL;
  };

  template<> class SCAL_TRAIT<ComplexConjugate>
  {
  public:
    typedef Complex SCAL;
  };

  template<> class SCAL_TRAIT<ComplexConjugate2>
  {
  public:
    typedef Complex SCAL;
  };



  /**
     Base class to linalg expression templates
  */
  template <class T>
  class VVecExpr 
  {
    const T data;
  public:
    ///
    VVecExpr (const T & d) : data(d) { ; }

    /// assign s * vector-expression data to v
    template <class TS>
    void AssignTo (TS s , BaseVector & v) const { data.AssignTo(s, v); }

    /// add s * vector-expression data to v
    template <class TS>
    void AddTo (TS s, BaseVector & v) const { data.AddTo(s, v); }
  };


  enum PARALLEL_STATUS { DISTRIBUTED, CUMULATED, NOT_PARALLEL };
  inline ostream & operator<< (ostream & ost, PARALLEL_STATUS stat)
  {
    switch (stat)
      {
      case DISTRIBUTED: ost << "distributed"; break;
      case CUMULATED: ost << "cumulated"; break;
      default: ost << "sequential";
      }
    return ost;
  }
  

  /**
     Base vector for linalg
  */
  class NGS_DLL_HEADER BaseVector : public enable_shared_from_this_virtual<BaseVector>
  {
  protected:
    /// size of vector
    size_t size;
    /// number of doubles per entry
    int entrysize = 1;
    ///
    BaseVector () { ; }
    
  public:
    ///
    virtual ~BaseVector () { ; }
    
    ///
    template <typename T> 
      BaseVector & operator= (const VVecExpr<T> & v)
    {
      v.AssignTo (1.0, *this);
      return *this;
    }
    
    ///
    BaseVector & operator= (const BaseVector & v)
    {
      Set (1.0, v);
      return *this;
    }
    ///
    BaseVector & operator= (double s)
    {
      SetScalar (s);
      return *this;
    }
    ///
    BaseVector & operator= (Complex s)
    {
      SetScalar (s);
      return *this;
    }

    virtual void SetZero ()
    {
      SetScalar(0);
    }
    
    ///
    template <typename T>
    BaseVector & operator+= (const VVecExpr<T> & v)
    {
      v.AddTo (1.0, *this);
      return *this;
    }

    ///
    BaseVector & operator+= (const BaseVector & v)
    {
      Add (1.0, v);
      return *this;
    }

    ///
    template <typename T>
    BaseVector & operator-= (const VVecExpr<T> & v)
    {
      v.AddTo (-1.0, *this);
      return *this;
    }

    ///
    BaseVector & operator-= (const BaseVector & v)
    {
      Add (-1.0, v);
      return *this;
    }

    ///
    BaseVector & operator*= (double s)
    {
      return Scale (s);
    }

    ///
    BaseVector & operator*= (Complex s)
    {
      return Scale (s);
    }

    ///
    BaseVector & operator/= (double s)
    {
      if (s == 0)
	throw Exception ("BaseVector::operator/=: division by zero");
      return Scale (1/s);
    }

    ///
    BaseVector & operator/= (Complex s)
    {
      if (s == 0.0)
	throw Exception ("BaseVector::operator/=: division by zero");
      return Scale (1.0/s);
    }

    template <class SCAL>
    S_BaseVector<SCAL> & Spec()
    {
      return dynamic_cast<S_BaseVector<SCAL>&> (*this);
    }
  
    template <class SCAL>
    const S_BaseVector<SCAL> & Spec() const
    {
      return dynamic_cast<const S_BaseVector<SCAL>&> (*this);
    }

    size_t Size() const throw () { return size; }
    T_Range<size_t> Range() const { return T_Range<size_t> (0, size); }
    // one entry has the size of that many doubles
    int EntrySize() const throw () { return entrysize; }
    // one entry has the size of that many scalars (double or complex)
    virtual int EntrySizeScal() const throw () = 0;
    virtual void * Memory () const = 0;
    virtual FlatVector<double> FVDouble () const = 0;
    virtual FlatVector<Complex> FVComplex () const = 0;

    template <typename SCAL = double>
    FlatSysVector<SCAL> SV () const
    {
      return FlatSysVector<SCAL> (Size(), EntrySize() * sizeof(double)/sizeof(SCAL), (SCAL*)Memory());
    }

    template <typename T>
      FlatVector<T> FV () const;
    
    /*
    template <class TSCAL>
    TSCAL InnerProduct (const BaseVector & v2) const 
    {
      return dynamic_cast<const S_BaseVector<TSCAL>&> (*this) . 
	InnerProduct (v2);
    }
    */

    virtual double InnerProductD (const BaseVector & v2) const;
    virtual Complex InnerProductC (const BaseVector & v2, bool conjuagte = false) const;
    
    virtual double L2Norm () const;
    virtual bool IsComplex() const { return false; }

    virtual BaseVector & Scale (double scal);
    virtual BaseVector & Scale (Complex scal);

    virtual BaseVector & SetScalar (double scal);
    virtual BaseVector & SetScalar (Complex scal);

    virtual BaseVector & Set (double scal, const BaseVector & v);
    virtual BaseVector & Set (Complex scal, const BaseVector & v);

    virtual BaseVector & Add (double scal, const BaseVector & v);
    virtual BaseVector & Add (Complex scal, const BaseVector & v);

    virtual ostream & Print (ostream & ost) const;
    virtual void Save(ostream & ost) const;
    virtual void Load(istream & ist);
    virtual void SaveText(ostream & ost) const;
    virtual void LoadText(istream & ist);

    virtual Array<MemoryUsage> GetMemoryUsage () const;
    virtual size_t CheckSum () const;
    // 
    // virtual shared_ptr<BaseVector> CreateVector () const = 0;
    virtual AutoVector CreateVector () const = 0;
    virtual unique_ptr<MultiVector> CreateMultiVector (size_t cnt) const;

    virtual void SetRandom ();

    inline AutoVector Range (size_t begin, size_t end) const;
    inline AutoVector Range (size_t end) const;
    // { return Range(T_Range(begin, end)); }
    virtual AutoVector Range (T_Range<size_t> range) const;
    virtual AutoVector Range (DofRange range) const;
      // { return Range(T_Range<size_t>(range)); }

    static bool IsRegularIndex (int index) { return index >= 0; }
    virtual void GetIndirect (FlatArray<int> ind, 
                              FlatVector<double> v) const = 0;
    virtual void GetIndirect (FlatArray<int> ind, 
                              FlatVector<Complex> v) const = 0;
    void SetIndirect (FlatArray<int> ind, FlatVector<double> v);
    void SetIndirect (FlatArray<int> ind, FlatVector<Complex> v);
    void AddIndirect (FlatArray<int> ind, FlatVector<double> v, bool use_atomic = false);
    void AddIndirect (FlatArray<int> ind, FlatVector<Complex> v, bool use_atomic = false);

    virtual shared_ptr<BaseVector> GetLocalVector () const 
    { return const_cast<BaseVector*>(this)->shared_from_this(); }
    
    virtual void Cumulate () const;
    virtual void Distribute() const;
    virtual PARALLEL_STATUS GetParallelStatus () const;
    virtual void SetParallelStatus (PARALLEL_STATUS stat) const;
    virtual optional<NgMPI_Comm> GetCommunicator() const { return nullopt; }


    virtual shared_ptr<BaseVector> CreateDeviceVector(bool unified) const;
    static std::map<type_index, function<shared_ptr<BaseVector>(const BaseVector&,bool)>> devveccreator;
    static void RegisterDeviceVectorCreator (type_index type,
                                             function<shared_ptr<BaseVector>(const BaseVector&,bool)> creator)
    {
      devveccreator[type] = creator;
    }

    
    const MemoryTracer& GetMemoryTracer() const { return mt; }
  private:
    MemoryTracer mt = { "BaseVector" };
  };
  

  AutoVector CreateBaseVector(size_t size, bool is_complex = false, int es = 1);
  
  
  class NGS_DLL_HEADER AutoVector 
  {
    shared_ptr<BaseVector> vec;
  public:
    AutoVector () = default;
    AutoVector (AutoVector && av2) = default; //  : vec(std::move(av2.vec)) { } 
    
    AutoVector (shared_ptr<BaseVector> hvec) : vec(std::move(hvec)) { }
    
    AutoVector (unique_ptr<BaseVector> hvec) : vec(std::move(hvec)) { } 

    template<typename U>
    AutoVector (unique_ptr<U> hvec) : vec(std::move(hvec)) { } 

    ~AutoVector();

    auto Size() const { return vec->Size(); }
    
    template <typename T> 
    BaseVector & operator= (const VVecExpr<T> & v)
    {
      v.AssignTo (1.0, *vec);
      return *this;
    }

    ///
    BaseVector & operator= (const BaseVector & v)
    {
      vec->Set (1.0, v);
      return *this;
    }
    ///
    BaseVector & operator= (const AutoVector & v)
    {
      vec->Set (1.0, *v);
      return *this;
    }

    template <typename T>
    auto & operator+= (const VVecExpr<T> & v)
    {
      (*vec) += v;
      return *this;
    }

    auto & operator+= (const BaseVector & v)
    {
      (*vec) += v;
      return *this;
    }
    
    template <typename T>
    auto & operator-= (const VVecExpr<T> & v)
    {
      (*vec) -= v;
      return *this;
    }

    auto & operator-= (const BaseVector & v)
    {
      (*vec) -= v;
      return *this;
    }

    auto & operator*= (double s)
    {
      (*vec) *= s;
      return *this;
    }

    ///
    auto & operator*= (Complex s)
    {
      (*vec) *= s;
      return *this;
    }

    auto & operator/= (double s)
    {
      (*vec) /= s;
      return *this;
    }

    ///
    auto & operator/= (Complex s)
    {
      (*vec) /= s;
      return *this;
    }

    auto & SetRandom ()
    {
      vec->SetRandom();
      return *this;
    }

    
    ///
    BaseVector & AssignPointer (AutoVector && v)
    {
      vec = std::move(v.vec);
      // size = v.size;
      // entrysize = v.entrysize;
      return *this;
    }
    ///
    BaseVector & operator= (double s)
    {
      vec->SetScalar (s);
      return *this;
    }
    ///
    BaseVector & operator= (Complex s)
    {
      vec->SetScalar (s);
      return *this;
    }
    
    operator shared_ptr<BaseVector> () && { return std::move(vec); }
    operator shared_ptr<BaseVector> () & { return vec; }
    BaseVector & operator* () { return *vec; }
    const BaseVector & operator* () const { return *vec; }
    operator BaseVector & () { return *vec; }
    operator const BaseVector & () const { return *vec; }

    AutoVector Range (size_t begin, size_t end) const { return vec->Range(begin,end); }
    AutoVector Range (size_t end) const { return vec->Range(0,end); }
    AutoVector Range (T_Range<size_t> range) const { return vec->Range(range); }
    
    template <typename T>
    auto FV () const { return vec->FV<T>(); }
    

    void * Memory () const throw () 
    {
      return vec->Memory();
    }

    FlatVector<double> FVDouble () const 
    {
      return vec->FVDouble();
    }
    
    FlatVector<Complex> FVComplex () const
    {
      return vec->FVComplex();
    }

    AutoVector CreateVector () const
    {
      return vec->CreateVector();
    }

    double InnerProductD (const BaseVector & v2) const
    {
      return vec->InnerProductD (v2);
    }

    Complex InnerProductC (const BaseVector & v2, bool conjugate) const
    {
      return vec->InnerProductC (v2, conjugate);
    }

    double L2Norm () const
    {
      return vec->L2Norm();
    }

    bool IsComplex() const 
    {
      return vec->IsComplex();
    }

    BaseVector & Scale (double scal)
    {
      return vec->Scale(scal);
    }

    BaseVector & Scale (Complex scal)
    {
      return vec->Scale(scal);
    }

    BaseVector & SetScalar (double scal)
    {
      return vec->SetScalar(scal);
    }
    BaseVector & SetScalar (Complex scal)
    {
      return vec->SetScalar(scal);
    }

    BaseVector & Set (double scal, const BaseVector & v)
    {
      return vec->Set (scal,v);
    }
    BaseVector & Set (Complex scal, const BaseVector & v)
    {
      return vec->Set (scal,v);
    }

    BaseVector & Add (double scal, const BaseVector & v)
    {
      return vec->Add (scal,v);
    }
    BaseVector & Add (Complex scal, const BaseVector & v)
    {
      return vec->Add (scal,v);
    }

    ostream & Print (ostream & ost) const
    {
      return vec->Print (ost);
    }


    void GetIndirect (FlatArray<int> ind, 
                      FlatVector<double> v) const
    {
      vec -> GetIndirect (ind, v);
    }
    void GetIndirect (FlatArray<int> ind, 
                      FlatVector<Complex> v) const
    {
      vec -> GetIndirect (ind, v);
    }

    void SetIndirect (FlatArray<int> ind, FlatVector<double> v)
    {
      vec->SetIndirect (ind,v);
    }
    void SetIndirect (FlatArray<int> ind, FlatVector<Complex> v)
    {
      vec->SetIndirect (ind,v);      
    }
    
    void AddIndirect (FlatArray<int> ind, FlatVector<double> v, bool use_atomic = false)
    {
      vec->AddIndirect (ind, v, use_atomic);
    }
    void AddIndirect (FlatArray<int> ind, FlatVector<Complex> v, bool use_atomic = false)
    {
      vec->AddIndirect (ind, v, use_atomic);
    }
    

    void Cumulate () const 
    { vec -> Cumulate(); }

    void Distribute() const
    { vec -> Distribute(); }

    PARALLEL_STATUS GetParallelStatus () const
    { return vec -> GetParallelStatus(); }

    void SetParallelStatus (PARALLEL_STATUS stat) const
    { vec -> SetParallelStatus(stat); }
  };

  AutoVector BaseVector::Range (size_t begin, size_t end) const
  {
    return Range(T_Range(begin, end));
  }

  AutoVector BaseVector::Range (size_t end) const
  {
    return Range (0, end);
  }
  

  template <>
  inline FlatVector<double> BaseVector::FV<double> () const
  {
    return FVDouble();
  }

  template <>
  inline FlatVector<Complex> BaseVector::FV<Complex> () const
  {
    return FVComplex();
  }

  template <typename T>
  inline FlatVector<T> BaseVector::FV () const
  {
    typedef typename mat_traits<T>::TSCAL TSCAL;
    return FlatVector<T> (Size(), static_cast<T*> (static_cast<void*>(FV<TSCAL>().Addr(0))));
  }






  /**
     Decision between double or Complex
  */



  template <class SCAL>
  class NGS_DLL_HEADER S_BaseVector : virtual public BaseVector
  {
  public:
    S_BaseVector () throw () { ; }
    virtual ~S_BaseVector() { ; }

    S_BaseVector & operator= (double s);
    virtual BaseVector & SetScalar (double scal) override;

    virtual bool IsComplex() const override
    { return typeid(SCAL) == typeid(Complex); }
    
    virtual int EntrySizeScal() const throw () override
    { return EntrySize() * sizeof(double)/sizeof(SCAL); }
    
    virtual SCAL InnerProduct (const BaseVector & v2, bool conjugate = false) const;

    virtual double InnerProductD (const BaseVector & v2) const override;
    virtual Complex InnerProductC (const BaseVector & v2, bool conjugate = false) const override;


    virtual FlatVector<double> FVDouble () const override;
    virtual FlatVector<Complex> FVComplex () const override;

    virtual FlatVector<SCAL> FVScal () const 
    {
      return FlatVector<SCAL> (size * entrysize * sizeof(double)/sizeof(SCAL), 
                               (SCAL*)Memory());
    }


    virtual void GetIndirect (FlatArray<int> ind, 
                              FlatVector<double> v) const override;
    virtual void GetIndirect (FlatArray<int> ind, 
                              FlatVector<Complex> v) const override;

  };


  template <>
  double S_BaseVector<double> :: InnerProduct (const BaseVector & v2, bool conjugate) const;


#if !defined(FILE_BASEVECTOR_CPP)
  extern template class S_BaseVector<double>;
  extern template class S_BaseVector<Complex>;
#endif




  class BlockVector;
  extern NGS_DLL_HEADER BlockVector & dynamic_cast_BlockVector (BaseVector & x);
  extern NGS_DLL_HEADER const BlockVector & dynamic_cast_BlockVector (const BaseVector & x);

  class BlockVector : public BaseVector
  {
    Array<shared_ptr<BaseVector>> vecs;
    BitArray ispar;
    NgMPI_Comm comm;
  public:
    BlockVector (const Array<shared_ptr<BaseVector>> & avecs);

    size_t NBlocks() const throw () { return vecs.Size(); }
    virtual int EntrySizeScal() const throw () override { return vecs[0]->EntrySizeScal(); }
    shared_ptr<BaseVector> & operator[] (size_t i) const { return vecs[i]; }

    void * Memory () const override;
    FlatVector<double> FVDouble () const override;
    FlatVector<Complex> FVComplex () const override;
    void GetIndirect (FlatArray<int> ind,
                      FlatVector<double> v) const override;
    void GetIndirect (FlatArray<int> ind,
                      FlatVector<Complex> v) const override;

    bool IsComplex() const override;

    AutoVector CreateVector () const override;

    double InnerProductD (const BaseVector & v2) const override;
    Complex InnerProductC (const BaseVector & v2,
                           bool conjugate = false) const override;
    double L2Norm () const override;
    
    BaseVector & Scale (double scal) override;
    BaseVector & Scale (Complex scal) override;
    BaseVector & SetScalar (double scal) override;
    BaseVector & SetScalar (Complex scal) override;

    ostream & Print (ostream & ost) const override;

    BaseVector & Set (double scal, const BaseVector & v) override;
    BaseVector & Add (double scal, const BaseVector & v) override;
    BaseVector & Set (Complex scal, const BaseVector & v) override;
    BaseVector & Add (Complex scal, const BaseVector & v) override;
  };

  







  

  /* ********************* Expression templates ***************** */



  template <> class VVecExpr<BaseVector>
  {
    const BaseVector & v;
  public:
    VVecExpr (const BaseVector & av) : v(av) { ; }

    template <class TS>
    void AssignTo (TS s, BaseVector & v2) const { v2.Set (s, v); }
    template <class TS>
    void AddTo (TS s, BaseVector & v2) const { v2.Add (s,  v); }
  };



  /* ***************************** VSumExpr ************************** */

  ///
  template <class TA, class TB>
  class VSumExpr
  {
    const TA a;
    const TB b;
  
  public:
    VSumExpr (const TA & aa, const TB & ab) : a(aa), b(ab) { ; }

    template <class TS>
    void AssignTo (TS s, BaseVector & v) const
    { 
      a.AssignTo (s, v);
      b.AddTo (s, v);
    }
    template <class TS>
    void AddTo (TS s, BaseVector & v) const
    { 
      a.AddTo (s, v);
      b.AddTo (s, v);
    }
  };


  /*
  inline VVecExpr<VSumExpr<VVecExpr<BaseVector>, VVecExpr<BaseVector> > >
  operator+ (const BaseVector & a, const BaseVector & b)
  {
    typedef VSumExpr<VVecExpr<BaseVector>, VVecExpr<BaseVector> > TRES;
    return TRES (a, b);
  }
  */
  inline auto operator+ (const BaseVector & a, const BaseVector & b)
  {
    return VVecExpr<VSumExpr<VVecExpr<BaseVector>,VVecExpr<BaseVector>>>{{a,b}};
  }

  /*
  template <class TA>
  inline VVecExpr<VSumExpr<VVecExpr<TA>, VVecExpr<BaseVector> > >
  operator+ (const VVecExpr<TA> & a, const BaseVector & b)
  {
    typedef VSumExpr<VVecExpr<TA>, VVecExpr<BaseVector> > TRES;
    return TRES (a, b);
  }
  */
  template <class TA>
  inline auto operator+ (const VVecExpr<TA> & a, const BaseVector & b)
  {
    return VSumExpr<VVecExpr<TA>,VVecExpr<BaseVector>> (a,b);
  }
  

  template <class TB>
  inline VVecExpr<VSumExpr<VVecExpr<BaseVector>, VVecExpr<TB> > >
  operator+ (const BaseVector & a, const VVecExpr<TB> & b)
  {
    typedef VSumExpr<VVecExpr<BaseVector>, VVecExpr<TB> > TRES;
    return TRES (a, b);
  }

  template <class TA, class TB>
  inline VVecExpr<VSumExpr<VVecExpr<TA>, VVecExpr<TB> > >
  operator+ (const VVecExpr<TA> & a, const VVecExpr<TB> & b)
  {
    typedef VSumExpr<VVecExpr<TA>, VVecExpr<TB> > TRES;
    return TRES (a, b);
  }








  /* ***************************** VSubExpr ************************** */

  ///
  template <class TA, class TB>
  class VSubExpr
  {
    const TA a;
    const TB b;
  
  public:
    VSubExpr (const TA & aa, const TB & ab) : a(aa), b(ab) { ; }


    template <class TS>
    void AssignTo (TS s, BaseVector & v) const
    { 
      a.AssignTo (s, v);
      b.AddTo (-s, v);
    }
    template <class TS>
    void AddTo (TS s, BaseVector & v) const
    { 
      a.AddTo (s, v);
      b.AddTo (-s, v);
    }
  };



  inline VVecExpr<VSubExpr<VVecExpr<BaseVector>, VVecExpr<BaseVector> > >
  operator- (const BaseVector & a, const BaseVector & b)
  {
    typedef VSubExpr<VVecExpr<BaseVector>, VVecExpr<BaseVector> > TRES;
    return TRES (a, b);
  }

  template <class TA>
  inline VVecExpr<VSubExpr<VVecExpr<TA>, VVecExpr<BaseVector> > >
  operator- (const VVecExpr<TA> & a, const BaseVector & b)
  {
    typedef VSubExpr<VVecExpr<TA>, VVecExpr<BaseVector> > TRES;
    return TRES (a, b);
  }

  template <class TB>
  inline VVecExpr<VSubExpr<VVecExpr<BaseVector>, VVecExpr<TB> > >
  operator- (const BaseVector & a, const VVecExpr<TB> & b)
  {
    typedef VSubExpr<VVecExpr<BaseVector>, VVecExpr<TB> > TRES;
    return TRES (a, b);
  }

  template <class TA, class TB>
  inline VVecExpr<VSubExpr<VVecExpr<TA>, VVecExpr<TB> > >
  operator- (const VVecExpr<TA> & a, const VVecExpr<TB> & b)
  {
    typedef VSubExpr<VVecExpr<TA>, VVecExpr<TB> > TRES;
    return TRES (a, b);
  }





  /* ************************* Scal * Vec ******************** */


  ///
  template <class TA, class TSCAL>
  class VScaleExpr
  {
    const TA a;
    const TSCAL scal;
  
  public:
    VScaleExpr (const TA & aa, const TSCAL & as) : a(aa), scal(as) { ; }


    template <class TS>
    void AssignTo (TS s, BaseVector & v) const
    { 
      a.AssignTo (scal * s, v);
    }
    template <class TS>
    void AddTo (TS s, BaseVector & v) const
    { 
      a.AddTo (scal * s, v);
    }
  };


  inline VVecExpr<VScaleExpr<VVecExpr<BaseVector>, double> >
  operator* (const BaseVector & a, const double & b)
  {
    typedef VScaleExpr<VVecExpr<BaseVector>, double> TRES;
    return TRES (a, b);
  }

  template <class TA>
  inline VVecExpr<VScaleExpr<VVecExpr<TA>, double> >
  operator* (const VVecExpr<TA> & a, const double & b)
  {
    typedef VScaleExpr<VVecExpr<TA>, double> TRES;
    return TRES (a, b);
  }



  inline VVecExpr<VScaleExpr<VVecExpr<BaseVector>, Complex> >
  operator* (const BaseVector & a, const Complex & b)
  {
    typedef VScaleExpr<VVecExpr<BaseVector>, Complex> TRES;
    return TRES (a, b);
  }

  template <class TA>
  inline VVecExpr<VScaleExpr<VVecExpr<TA>, Complex> >
  operator* (const VVecExpr<TA> & a, const Complex & b)
  {
    typedef VScaleExpr<VVecExpr<TA>, Complex> TRES;
    return TRES (a, b);
  }





  inline VVecExpr<VScaleExpr<VVecExpr<BaseVector>, double> >
  operator* (const double & b, const BaseVector & a)
  {
    typedef VScaleExpr<VVecExpr<BaseVector>, double> TRES;
    return TRES (a, b);
  }

  template <class TA>
  inline VVecExpr<VScaleExpr<VVecExpr<TA>, double> >
  operator* (const double & b, const VVecExpr<TA> & a)
  {
    typedef VScaleExpr<VVecExpr<TA>, double> TRES;
    return TRES (a, b);
  }



  inline VVecExpr<VScaleExpr<VVecExpr<BaseVector>, Complex> >
  operator* (const Complex & b, const BaseVector & a)
  {
    typedef VScaleExpr<VVecExpr<BaseVector>, Complex> TRES;
    return TRES (a, b);
  }

  template <class TA>
  inline VVecExpr<VScaleExpr<VVecExpr<TA>, Complex> >
  operator* (const Complex & b, const VVecExpr<TA> & a)
  {
    typedef VScaleExpr<VVecExpr<TA>, Complex> TRES;
    return TRES (a, b);
  }


  template <class TA>
  inline VVecExpr<VScaleExpr<VVecExpr<TA>,double > >
  operator- (const VVecExpr<TA> & a)
  {
    typedef VScaleExpr<VVecExpr<TA>, double> TRES;
    return TRES (a, -1);
  }




  /* *********************** operator<< ********************** */

  ///
  inline ostream & operator<< (ostream & ost, const BaseVector & v)
  {
    return v.Print(ost);
  }

  ///
  template <typename T = double>
  inline T InnerProduct (const BaseVector & v1, const BaseVector & v2, bool conjugate = false)
  {
    // return dynamic_cast<const S_BaseVector<double>&>(v1).InnerProduct(v2);
    if constexpr (is_same<T,double>::value)
                   return v1.InnerProductD(v2);
    else
      return v1.InnerProductC(v2, conjugate);
  }

  ///
  template <class IPTYPE>
  inline typename SCAL_TRAIT<IPTYPE>::SCAL S_InnerProduct (const BaseVector & v1, const BaseVector & v2)
  {
    return dynamic_cast<const S_BaseVector<typename SCAL_TRAIT<IPTYPE>::SCAL>&>(v1).InnerProduct(v2); 
  }

  template <> inline double 
  S_InnerProduct<double> (const BaseVector & v1, const BaseVector & v2)
  {
    return v1.InnerProductD (v2);
  }

  template <> inline Complex 
  S_InnerProduct<Complex> (const BaseVector & v1, const BaseVector & v2)
  {
    return v1.InnerProductC (v2);
  }


  template <> inline Complex 
  S_InnerProduct<ComplexConjugate> (const BaseVector & v1, const BaseVector & v2)
  {
    return v1.InnerProductC(v2, true);
    // return InnerProduct( v1.FVComplex(), Conj(v2.FVComplex()) );
  }

  template <>
  inline Complex S_InnerProduct<ComplexConjugate2> (const BaseVector & v1, const BaseVector & v2)
  {
    return v2.InnerProductC(v1, true);
    // return InnerProduct( v2.FVComplex(), Conj(v1.FVComplex()) );
  }

  ///
  inline double L2Norm (const BaseVector & v)
  {
    return v.L2Norm();
  }










  class DynamicBaseExpression 
  {
  protected:
  public:
    DynamicBaseExpression () { } 
    virtual ~DynamicBaseExpression() { }
    virtual AutoVector CreateVector() const = 0;
    virtual AutoVector Evaluate() const
    {
      auto vec = CreateVector();
      AssignTo (1, vec);
      return vec;
    }
    virtual void AssignTo (double s, BaseVector & v2) const = 0;
    virtual void AddTo (double s, BaseVector & v2) const = 0;
    virtual void AssignTo (Complex s, BaseVector & v2) const = 0;
    virtual void AddTo (Complex s, BaseVector & v2) const = 0;
  };


  class DynamicVecExpression : public DynamicBaseExpression
  {
  protected:
    shared_ptr<BaseVector> a;
  public:
    DynamicVecExpression (shared_ptr<BaseVector> aa) : a(aa) { ; }
    AutoVector CreateVector() const override
    { return a->CreateVector(); }
    void AssignTo (double s, BaseVector & v2) const override
    { v2.Set (s, *a); }
    void AddTo (double s, BaseVector & v2) const override
    { v2.Add (s, *a); }
    void AssignTo (Complex s, BaseVector & v2) const override
    { v2.Set (s, *a); }
    void AddTo (Complex s, BaseVector & v2) const override
    { v2.Add (s, *a); }
  };

  class DynamicSumExpression : public DynamicBaseExpression
  {
    shared_ptr<DynamicBaseExpression> a,b;
    
    AutoVector CreateVector() const override
    { return a->CreateVector(); }    
    void AssignTo (double s, BaseVector & v2) const override
    {
      a->AssignTo(s, v2);
      b->AddTo(s, v2);
    }
    void AddTo (double s, BaseVector & v2) const override
    {
      a->AddTo(s, v2);
      b->AddTo(s, v2);
    }
    void AssignTo (Complex s, BaseVector & v2) const override
    {
      a->AssignTo(s, v2);
      b->AddTo(s, v2);
    }
    void AddTo (Complex s, BaseVector & v2) const override
    {
      a->AddTo(s, v2);
      b->AddTo(s, v2);
    }
  public:
    DynamicSumExpression (shared_ptr<DynamicBaseExpression> aa,
                          shared_ptr<DynamicBaseExpression> ab)
      : a(aa), b(ab) { ; } 
  };

  class DynamicSubExpression : public DynamicBaseExpression
  {
    shared_ptr<DynamicBaseExpression> a,b;
    
    AutoVector CreateVector() const override
    { return a->CreateVector(); }    

    void AssignTo (double s, BaseVector & v2) const override
    {
      a->AssignTo(s, v2);
      b->AddTo(-s, v2);
    }
    void AddTo (double s, BaseVector & v2) const override
    {
      a->AddTo(s, v2);
      b->AddTo(-s, v2);
    }
    void AssignTo (Complex s, BaseVector & v2) const override
    {
      a->AssignTo(s, v2);
      b->AddTo(-s, v2);
    }
    void AddTo (Complex s, BaseVector & v2) const override
    {
      a->AddTo(s, v2);
      b->AddTo(-s, v2);
    }
  public:
    DynamicSubExpression (shared_ptr<DynamicBaseExpression> aa,
                          shared_ptr<DynamicBaseExpression> ab)
      : a(aa), b(ab) { ; } 
  };

  template <typename T>
  class DynamicScaleExpression : public DynamicBaseExpression
  {
    T scale;
    shared_ptr<DynamicBaseExpression> a;

    AutoVector CreateVector() const override
    { return a->CreateVector(); }    
    
    void AssignTo (double s, BaseVector & v2) const override
    {
      a->AssignTo(s*scale, v2);
    }
    void AddTo (double s, BaseVector & v2) const override
    {
      a->AddTo(s*scale, v2);
    }
    void AssignTo (Complex s, BaseVector & v2) const override
    {
      a->AssignTo(s*scale, v2);
    }
    void AddTo (Complex s, BaseVector & v2) const override
    {
      a->AddTo(s*scale, v2);
    }
  public:
    DynamicScaleExpression (T ascale, shared_ptr<DynamicBaseExpression> aa)
      : scale(ascale), a(aa) { ; } 
  };



  
  class DynamicVectorExpression 
  {
    shared_ptr<DynamicBaseExpression> ve;
  public:
    DynamicVectorExpression() { } 
    DynamicVectorExpression (shared_ptr<DynamicBaseExpression> ave) : ve(ave) { }
    DynamicVectorExpression (shared_ptr<BaseVector> v)
      : ve(make_shared<DynamicVecExpression>(v)) { }

    AutoVector Evaluate() const
    {
      return ve->Evaluate();
      /*
      auto vec = ve->CreateVector();
      ve->AssignTo (1, vec);
      return vec;
      */
    }

    AutoVector CreateVector() const
    {
      return ve->CreateVector();
    }
    
    void AssignTo (double s, BaseVector & v2) const
    { ve->AssignTo(s,v2); }
    void AddTo (double s, BaseVector & v2) const
    { ve->AddTo(s,v2); }
    void AssignTo (Complex s, BaseVector & v2) const
    { ve->AssignTo(s,v2); }
    void AddTo (Complex s, BaseVector & v2) const
    { ve->AddTo(s,v2); }
    auto Ptr() const { return ve; }
  };

  inline auto operator+ (DynamicVectorExpression a, DynamicVectorExpression b)
  {
    return DynamicVectorExpression(make_shared<DynamicSumExpression>(a.Ptr(),b.Ptr()));
  }

  inline auto operator- (DynamicVectorExpression a, DynamicVectorExpression b)
  {
    return DynamicVectorExpression(make_shared<DynamicSubExpression>(a.Ptr(),b.Ptr()));
  }

  template <typename T>
  inline auto operator* (T s, DynamicVectorExpression v)
  {
    return DynamicVectorExpression(make_shared<DynamicScaleExpression<T>>(s, v.Ptr()));
  }





  
}

#endif
