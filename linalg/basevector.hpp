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

enum PARALLEL_STATUS { DISTRIBUTED, CUMULATED, NOT_PARALLEL };



#ifdef PARALLEL

class ParallelDofs
  {
  protected:
    int ndof;

    /// proc 2 dofs
    Table<int> * exchangedofs;

    /// dof 2 procs
    Table<int> * dist_procs;

    /// mpi-datatype to send exchange dofs
    Array<MPI_Datatype> mpi_t;

    /// is this the master process ?
    BitArray ismasterdof;

    MPI_Comm comm;

  public:
    /*
    ParallelDofs (const MeshAccess & ama, const Array<Node> & adofnodes, 
		  int dim = 1, bool iscomplex = false);
    */

    virtual ~ParallelDofs()
    {
      ;
    }


    int GetNTasks() const { return exchangedofs->Size(); }

    const FlatArray<int>  GetExchangeDofs (int proc) const
    { return (*exchangedofs)[proc]; }

    const FlatArray<int>  GetDistantProcs (int dof) const
    { return (*dist_procs)[dof]; }

    bool IsMasterDof ( int localdof ) const
    { return ismasterdof.Test(localdof); }

    int GetNDof () const { return ndof; }

    int GetNDofGlobal () const
    {
      if (GetNTasks() == 1) return ndof;
      int nlocal = 0;
      for (int i = 0; i < ndof; i++)
	if (ismasterdof.Test(i)) nlocal++;
      return ngparallel::MyMPI_AllReduce (nlocal);
    }

    bool IsExchangeProc ( int proc ) const
    { return (*exchangedofs)[proc].Size() != 0; }

    MPI_Datatype MyGetMPI_Type ( int dest ) const
    { return mpi_t[dest]; }

    MPI_Comm GetCommunicator () const { return comm; }
  };

#else
class ParallelDofs 
{
protected:
  int ndof;
  
public:

    int GetNDofGlobal () const { return ndof; }
  
};

#endif








class NGS_DLL_HEADER BaseVector
{
protected:
  /// size of vector
  int size;
  /// number of doubles per entry
  int entrysize;
  ///
  const ParallelDofs * paralleldofs;
public:
  ///
  BaseVector () throw ()  
    : paralleldofs (NULL) 
  { ; }
  ///
  virtual ~BaseVector ();

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

  virtual FlatVector<double> FVDouble () const = 0;
  virtual FlatVector<Complex> FVComplex () const = 0;

  template <typename T>
  FlatVector<T> FV () const;

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

  virtual BaseVector * Range (int begin, int end) const;
  virtual BaseVector * Range (IntRange range) const;

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
		    FlatVector< Vec<S,double> > & v) const
  { 
    FlatVector<double> fv = FVDouble();
    // int es = EntrySize();
    for (int i = 0, ii = 0; i < ind.Size(); i++)
      if (ind[i] != -1)
	{
	  int base = S * ind[i];
	  for (int j = 0; j < S; j++)
	    v[ii++] = fv[base++];
	}
      else
	{
	  for (int j = 0; j < S; j++)
	    v[ii++] = 0;
	}
  }

    
  template<int S>
  void GetIndirect (const Array<int> & ind, 
		    FlatVector< Vec<S,Complex> > & v) const
  { 
    FlatVector<Complex> fv = FVComplex();
    // int es = EntrySize() / 2;
    for (int i = 0, ii = 0; i < ind.Size(); i++)
      if (ind[i] != -1)
	{
	  int base = S * ind[i];
	  for (int j = 0; j < S; j++)
	    v[ii++] = fv[base++];
	}
      else
	{
	  for (int j = 0; j < S; j++)
	    v[ii++] = 0.0;
	}
  }


  template<int S>
  void AddIndirect (const Array<int> & ind, 
		    const FlatVector< Vec<S,double> > & v)
  { 
    FlatVector<double> fv = FVDouble();
    // int es = EntrySize();
    
    for (int i = 0; i < ind.Size(); i++)
      if (ind[i] != -1)
	{
	  int base = S * ind[i];
	  for (int j = 0; j < S; j++)
	    fv[base++] += v[i](j);
	}
   }

  template<int S>
  void AddIndirect (const Array<int> & ind, 
		    const FlatVector< Vec<S,Complex> > & v)
  { 
    FlatVector<Complex> fv = FVComplex();
    // if(EntrySize() != 2*S)
    //       throw Exception("BaseVector::AddIndirect() wrong dimensions");

    for (int i = 0; i < ind.Size(); i++)
      if (ind[i] != -1)
	{
	  int base = S * ind[i];
	  for (int j = 0; j < S; j++)
	    fv[base++] += v[i](j);
	}
  }
  
  
  virtual void Cumulate () const { ; }
  virtual void Distribute() const { ; }
  virtual PARALLEL_STATUS GetParallelStatus () const
  { return NOT_PARALLEL; }
  virtual void SetParallelStatus (PARALLEL_STATUS stat) const { ; }
};



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
  return FlatVector<T> (Size(), FV<TSCAL>().Addr(0));
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

  virtual SCAL InnerProduct (const BaseVector & v2) const
  {
    return ngbla::InnerProduct (FVScal(), 
				dynamic_cast<const S_BaseVector&>(v2).FVScal());
  }

  virtual FlatVector<double> FVDouble () const;
  virtual FlatVector<Complex> FVComplex () const;

  virtual FlatVector<SCAL> FVScal () const throw() 
  {
    return FlatVector<SCAL> (size * entrysize, Memory());
  }
};




template <>
class NGS_DLL_HEADER S_BaseVector<Complex> : virtual public BaseVector
{
public:
  S_BaseVector () throw() { ; }
  ~S_BaseVector () { ; }

  virtual Complex InnerProduct (const BaseVector & v2) const
  {
    return ngbla::InnerProduct (FVScal(), 
				dynamic_cast<const S_BaseVector&>(v2).FVScal());
  }

  virtual FlatVector<double> FVDouble () const throw();
  virtual FlatVector<Complex> FVComplex () const throw();
  virtual FlatVector<Complex> FVScal () const throw() 
  {
    return FlatVector<Complex> (size * entrysize/2, Memory());
  }
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
