namespace ngla
{

class BaseVector;
typedef auto_ptr<BaseVector> TempVector;
template <class SCAL> class S_BaseVector;


class ComplexConjugate;
class ComplexConjugate2;

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


enum PARALLEL_STATUS { DISTRIBUTED, CUMULATED, NOT_PARALLEL };


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


/**
   Base vector for linalg
 */
class NGS_DLL_HEADER BaseVector
{
protected:
  /// size of vector
  int size;
  /// number of doubles per entry
  int entrysize;

public:
  ///
  BaseVector () throw ()  { ; }
  ///
  virtual ~BaseVector () throw ();

  ///
  template <typename T> 
  BaseVector & operator= (const VVecExpr<T> & v)
  {
    v.AssignTo (1.0, *this);
    return *this;
  }

  ///
  BaseVector & operator= (const BaseVector & v);
  ///
  BaseVector & operator= (double s);
  ///
  BaseVector & operator= (Complex s);

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

  int Size() const throw () { return size; }
  int EntrySize() const throw () { return entrysize; }
  virtual void * Memory () const throw () = 0;

  virtual FlatVector<double> FVDouble () const throw() = 0;
  virtual FlatVector<Complex> FVComplex () const throw() = 0;

  template <class TSCAL>
  TSCAL InnerProduct (const BaseVector & v2) const 
  {
    return dynamic_cast<const S_BaseVector<TSCAL>&> (*this) . 
      InnerProduct (v2);
  }

  virtual double L2Norm () const
  {
    return ngbla::L2Norm (FVDouble());
  }

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

  virtual void MemoryUsage (Array<MemoryUsageStruct*> & mu) const;

  // create vector, procs is the set of processors on which the vector exists
  // default 0 pointer means all procs
  virtual BaseVector * CreateVector ( const Array<int> * procs = 0) const;


  virtual void SetRandom ();

  // virtual BaseVector * Range (int begin, int end);
  virtual BaseVector * Range (int begin, int end) const;


  void GetIndirect (const FlatArray<int> & ind, 
		    const FlatVector<double> & v) const;
  void GetIndirect (const FlatArray<int> & ind, 
		    const FlatVector<Complex> & v) const;
  void SetIndirect (const FlatArray<int> & ind, 
		    const FlatVector<double> & v);
  void SetIndirect (const FlatArray<int> & ind, 
		    const FlatVector<Complex> & v);
  void AddIndirect (const FlatArray<int> & ind, 
		    const FlatVector<double> & v);
  void AddIndirect (const FlatArray<int> & ind, 
		    const FlatVector<Complex> & v);

  template<int S>
  void GetIndirect (const Array<int> & ind, 
		    FlatVector< Vec<S,double> > & v) const;
  template<int S>
  void GetIndirect (const Array<int> & ind, 
		    FlatVector< Vec<S,Complex> > & v) const;
  template<int S>
  void AddIndirect (const Array<int> & ind, 
		    const FlatVector< Vec<S,double> > & v);
  template<int S>
  void AddIndirect (const Array<int> & ind, 
		    const FlatVector< Vec<S,Complex> > & v);
  
  


  virtual const PARALLEL_STATUS Status () const 
  { return NOT_PARALLEL; }

  virtual void SetStatus ( PARALLEL_STATUS astatus )
  { 
    if ( astatus == CUMULATED || astatus == DISTRIBUTED )
      cerr << "ERROR -- BaseVector::SetStatus(status) called for parallel status" << endl;
  }

  virtual class ngparallel::ParallelDofs * GetParallelDofs () const
  { 
    cerr << "WARNING -- GetParallelDofs called for BaseVector, is not parallel" << endl; 
    return 0;
  }

  virtual void PrintParallelDofs() const
  { cerr << "ERROR -- PrintParallelDofs called for BaseVector, is not parallel" << endl; }

  virtual bool IsParallelVector () const
  { 
    if ( this->Status() == NOT_PARALLEL ) return false;
    return true;
  }

  virtual void PrintStatus ( ostream & ost ) const
  { cerr << "ERROR -- PrintStatus called for BaseVector, is not parallel" << endl; }

  virtual void AllReduce ( Array<int> * reduceprocs, Array<int> * sendtoprocs=0 ) const
  { cerr << "ERROR -- AllReduce called for BaseVector, is not parallel" << endl; }

  virtual void Distribute() const
  { cerr << "ERROR -- Distribute called for BaseVector, is not parallel" << endl; }

  virtual void ISend ( const int dest, int & request )
  { cerr << "ERROR -- ISend called for BaseVector, is not parallel" << endl; }
  virtual void Send ( const int dest )
  { cerr << "ERROR -- Send called for BaseVector, is not parallel" << endl; }

  virtual void IRecvVec ( const int dest, int & request )
  { cerr << "ERROR -- IRecvVec called for BaseVector, is not parallel" << endl; }

  virtual void SetParallelDofs ( ngparallel::ParallelDofs * aparalleldofs, const Array<int> * procs =0 )
  { 
    if ( aparalleldofs == 0 ) return;
    cerr << "ERROR -- SetParallelDofs called for BaseVector, is not parallel" << endl; 
  }


};



/**
   Decision between double or Complex
 */
template <class SCAL>
class NGS_DLL_HEADER S_BaseVector : public BaseVector
{
public:
  S_BaseVector () throw () { ; }
  virtual ~S_BaseVector() throw() { ; }

  S_BaseVector & operator= (double s);

#ifdef PARALLEL
  SCAL InnerProduct (const BaseVector & v2) const;
#else
  inline SCAL InnerProduct (const BaseVector & v2) const
  {
    return ngbla::InnerProduct (FVScal(), 
				dynamic_cast<const S_BaseVector&>(v2).FVScal());
  }
#endif


  virtual FlatVector<double> FVDouble () const throw();
  virtual FlatVector<Complex> FVComplex () const throw();

  virtual FlatVector<SCAL> FVScal () const throw() 
  {
    return FlatVector<SCAL> (size * entrysize, Memory());
  }

#ifdef PARALLEL
  /*
  virtual void SetParallelDofs ( ngparallel::ParallelDofs & aparalleldofs, const Array<int> * procs  ) = 0;

  /// values from reduceprocs are added up,
  /// vectors in sendtoprocs are set to the cumulated values
  virtual void AllReduce ( Array<int> * reduceprocs, Array<int> * sendtoprocs = 0 ) const = 0;

  virtual void Distribute() const = 0;

  virtual void ISend ( const int dest, MPI_Request & request ) = 0;
  virtual void Send ( const int dest ) = 0;

//   template <class T>
//   virtual void IRecvVec (  Array<T> & s, const int dest, MPI_Request & request ) = 0;
*/
 #endif
};




template <>
class NGS_DLL_HEADER S_BaseVector<Complex> : public BaseVector
{
public:
  S_BaseVector () throw() { ; }
  ~S_BaseVector () throw() { ; }


#ifdef PARALLEL
  Complex InnerProduct (const BaseVector & v2) const;
#else
  Complex InnerProduct (const BaseVector & v2) const
  {
    return ngbla::InnerProduct (FVScal(), 
				dynamic_cast<const S_BaseVector&>(v2).FVScal());
  }
#endif

  virtual FlatVector<double> FVDouble () const throw();
  virtual FlatVector<Complex> FVComplex () const throw();
  virtual FlatVector<Complex> FVScal () const throw() 
  {
    return FlatVector<Complex> (size * entrysize/2, Memory());
  }

#ifdef PARALLEL
  /*
  virtual void SetParallelDofs ( ngparallel::ParallelDofs & aparalleldofs, const Array<int> * procs  ) = 0;

  /// values from reduceprocs are added up,
  /// vectors in sendtoprocs are set to the cumulated values
  virtual void AllReduce ( Array<int> * reduceprocs, Array<int> * sendtoprocs = 0) const = 0;

  virtual void Distribute() const = 0;

  virtual void ISend ( const int dest, MPI_Request & request ) = 0;
  virtual void Send ( const int dest ) = 0;

  */
#endif
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



inline VVecExpr<VSumExpr<VVecExpr<BaseVector>, VVecExpr<BaseVector> > >
operator+ (const BaseVector & a, const BaseVector & b)
{
  typedef VSumExpr<VVecExpr<BaseVector>, VVecExpr<BaseVector> > TRES;
  return TRES (a, b);
}

template <class TA>
inline VVecExpr<VSumExpr<VVecExpr<TA>, VVecExpr<BaseVector> > >
operator+ (const VVecExpr<TA> & a, const BaseVector & b)
{
  typedef VSumExpr<VVecExpr<TA>, VVecExpr<BaseVector> > TRES;
  return TRES (a, b);
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
inline double InnerProduct (const BaseVector & v1, const BaseVector & v2)
{
  return dynamic_cast<const S_BaseVector<double>&>(v1).InnerProduct(v2); 
}


///
template <class IPTYPE>
inline typename SCAL_TRAIT<IPTYPE>::SCAL S_InnerProduct (const BaseVector & v1, const BaseVector & v2)
{
  return dynamic_cast<const S_BaseVector<typename SCAL_TRAIT<IPTYPE>::SCAL>&>(v1).InnerProduct(v2); 
}

template <>
inline Complex 
S_InnerProduct<ComplexConjugate> (const BaseVector & v1, const BaseVector & v2)
{
  return InnerProduct( v1.FVComplex(), Conj(v2.FVComplex()) );
}

template <>
inline Complex S_InnerProduct<ComplexConjugate2> (const BaseVector & v1, const BaseVector & v2)
{
  return InnerProduct( v2.FVComplex(), Conj(v1.FVComplex()) );
}

///
inline double L2Norm (const BaseVector & v)
{
  return v.L2Norm();
}

}
