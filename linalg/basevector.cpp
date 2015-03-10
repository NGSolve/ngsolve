/*********************************************************************/
/* File:   basevector.cpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   7. Feb. 2003                                              */
/*********************************************************************/

/* 
   base class in vector hierarchy
*/

#define FILE_BASEVECTOR_CPP

#include <la.hpp>

#include "../ngstd/taskhandler.hpp"

namespace ngla
{

  
  TaskHandler * task_handler = nullptr;
  void SetTaskHandler (TaskHandler * handler)
  {
    task_handler = handler;
  }




  /*
  BaseVector :: BaseVector ()
    : paralleldofs (NULL) 
  {
    ;
  }
  
  BaseVector :: ~BaseVector () 
  { 
    ;
  }

  BaseVector & BaseVector :: operator= (const BaseVector & v)
  {
    Set (1.0, v);
    return *this;
  }

  BaseVector & BaseVector :: operator= (double s)
  {
    SetScalar (s);
    return *this;
  }

  BaseVector & BaseVector :: operator= (Complex s)
  {
    SetScalar (s);
    return *this;
  }
  */

  double BaseVector :: L2Norm () const
  {
    return ngbla::L2Norm (FVDouble());
  }

  BaseVector & BaseVector :: Scale (double scal)
  {
    if (scal == 1) return *this;

    static Timer t("BaseVector::Scale");
    RegionTimer reg(t);

    // FVDouble() *= scal;

    FlatVector<double> me = FVDouble();

    switch (omp_status)
      {
      case OMP_STATUS::NONE:
        me *= scal;
        break;

      case OMP_STATUS::USE:
      case OMP_STATUS::TASKS:
#pragma omp for
        for (int i = 0; i < me.Size(); i++)
          me(i) *= scal;
        break;

      case OMP_STATUS::CREATE:
#pragma omp parallel for
        for (int i = 0; i < me.Size(); i++)
          me(i) *= scal;
        break;

      case OMP_STATUS::CREATE_IF_NOT_PARALLEL:
        if (omp_in_parallel())
#pragma omp for
          for (int i = 0; i < me.Size(); i++)
            me(i) *= scal;
        else
#pragma omp parallel for
          for (int i = 0; i < me.Size(); i++)
            me(i) *= scal;
        break;
      }


    return *this;
  }
  BaseVector & BaseVector :: Scale (Complex scal)
  {
    FVComplex() *= scal;
    return *this;
  }

  BaseVector & BaseVector :: SetScalar (double scal)
  {
    static Timer t("BaseVector::SetScalar");
    RegionTimer reg(t);

    if (task_handler)
      {
	FlatVector<> fv = FVDouble();

	task_handler -> CreateTask 
	  ( [fv,scal] (int thd)
	    {
	      int n = fv.Size();
	      int num_thds = omp_get_num_threads();
	      IntRange r (thd*n/num_thds, (thd+1)*n/num_thds);		  
	      fv.Range(r) = scal;
	    } );

	return *this; 
      }


    FVDouble() = scal;
    return *this;
  }

  BaseVector & BaseVector :: SetScalar (Complex scal)
  {
    FVComplex() = scal;
    return *this;
  }

  BaseVector & BaseVector :: Set (double scal, const BaseVector & v)
  {
    static Timer t("BaseVector::Set");
    RegionTimer reg(t);

    if(Size() != v.Size())
      throw Exception (string ("BaseVector::Set: size of me = ") + ToString(Size()) + " != size of other = " + ToString(v.Size()));


    FlatVector<double> me = FVDouble();
    FlatVector<double> you = v.FVDouble();
    

    switch (omp_status)
      {
      case OMP_STATUS::NONE:
        me = scal * you;
        break;

      case OMP_STATUS::USE:
      case OMP_STATUS::TASKS:
#pragma omp for
        for (int i = 0; i < me.Size(); i++)
          me(i) = scal * you(i);
        break;

      case OMP_STATUS::CREATE:
#pragma omp parallel for
        for (int i = 0; i < me.Size(); i++)
          me(i) = scal * you(i);
        break;

      case OMP_STATUS::CREATE_IF_NOT_PARALLEL:
        if (omp_in_parallel())
#pragma omp for
          for (int i = 0; i < me.Size(); i++)
            me(i) = scal * you(i);
        else
#pragma omp parallel for
          for (int i = 0; i < me.Size(); i++)
            me(i) = scal * you(i);
        break;
      }

    return *this;
  }

  BaseVector & BaseVector :: Set (Complex scal, const BaseVector & v)
  {
    if(Size() != v.Size())
        throw Exception (string ("BaseVector::Set: size of me = ") + ToString(Size() + " != size of other = " + ToString(v.Size())));
    FVComplex() = scal * v.FVComplex();
    return *this;
  }
    


  inline void VectorAddRec (FlatVector<> a, FlatVector<> b, double scal)
  {
    int n = a.Size();
    if (n <= 32000)
      {
	a += scal * b;
	return;
      }

    int n1 = n/2;
#pragma omp task
    {
      VectorAddRec (a.Range(0,n1), b.Range(0,n1), scal);
    }
    //#pragma omp task
    {
      VectorAddRec (a.Range(n1,n), b.Range(n1,n), scal);
    }
    //#pragma omp taskwait
  }


  BaseVector & BaseVector :: Add (double scal, const BaseVector & v)
  {


    if (task_handler)
      {
	static Timer t("BaseVector::Add (taskhandler)");
	RegionTimer reg(t);
	
	FlatVector<double> me = FVDouble();
	FlatVector<double> you = v.FVDouble();
	
	task_handler -> CreateTask 
	  ( [me,you] (int thd)
	    {
	      int n = me.Size();
	      int num_thds = omp_get_num_threads();
	      IntRange r (thd*n/num_thds, (thd+1)*n/num_thds);		  
	      me.Range(r) += you.Range(r);
	    } );

	return *this;
      }
    


    static Timer t("BaseVector::Add");
    RegionTimer reg(t);

    if(Size() != v.Size())
        throw Exception (string ("BaseVector::Add: size of me = ") + ToString(Size() + " != size of other = " + ToString(v.Size())));
    // FVDouble() += scal * v.FVDouble();
    

    FlatVector<double> me = FVDouble();
    FlatVector<double> you = v.FVDouble();

    switch (omp_status)
      {
      case OMP_STATUS::NONE:
        me = scal * you;
        break;

      case OMP_STATUS::USE:
      case OMP_STATUS::TASKS:
#pragma omp for
        for (int i = 0; i < me.Size(); i++)
          me(i) += scal * you(i);
        break;

      case OMP_STATUS::CREATE:
#pragma omp parallel for
        for (int i = 0; i < me.Size(); i++)
          me(i) += scal * you(i);
        break;

      case OMP_STATUS::CREATE_IF_NOT_PARALLEL:
        if (omp_in_parallel())
	  {

	    /*
	      #pragma omp for
	      for (int i = 0; i < me.Size(); i++)
	      me(i) += scal * you(i);
	    */
	    auto r = OmpSplit (ngstd::Range (me.Size()));
	    me.Range(r) += scal * you.Range(r);
#pragma omp barrier

	  }
        else
#pragma omp parallel for
          for (int i = 0; i < me.Size(); i++)
            me(i) += scal * you(i);
        break;
      }

    return *this;
  }

  BaseVector & BaseVector :: Add (Complex scal, const BaseVector & v)
  {
    if(Size() != v.Size())
        throw Exception (string ("BaseVector::Add: size of me = ") + ToString(Size() + " != size of other = " + ToString(v.Size())));
    FVComplex() += scal * v.FVComplex();
    return *this;
  }


  double BaseVector :: InnerProductD (const BaseVector & v2) const
  {
    return dynamic_cast<const S_BaseVector<double>&> (*this) . 
      InnerProduct (v2);
  }
  
  Complex BaseVector :: InnerProductC (const BaseVector & v2) const
  {
    return dynamic_cast<const S_BaseVector<Complex>&> (*this) . 
      InnerProduct (v2);
  }



  AutoVector BaseVector ::Range (int begin, int end) const
  {
    throw Exception ("BaseVector::Range const called");
  }

  AutoVector BaseVector :: Range (IntRange range) const
  {
    throw Exception ("BaseVector::Range (IntRange) const called");
  }


  ostream & BaseVector :: Print (ostream & ost) const
  {
    throw Exception ("BaseVector::Print called");
  }
  
  void BaseVector :: Save(ostream & ost) const
  {
    FlatVector<double> fv = FVDouble();
    for (int i = 0; i < fv.Size(); i++)
      SaveBin (ost, fv(i));
  }

  void BaseVector :: Load(istream & ist) 
  {
    FlatVector<double> fv = FVDouble();
    for (int i = 0; i < fv.Size(); i++)
      LoadBin (ist, fv(i));
  }

  void BaseVector :: SaveText(ostream & ost) const
  {
    FlatVector<double> fv = FVDouble();
    for (int i = 0; i < fv.Size(); i++)
      ost << fv(i) << " ";
  }

  void BaseVector :: LoadText(istream & ist) 
  {
    FlatVector<double> fv = FVDouble();
    for (int i = 0; i < fv.Size(); i++)
      ist >> fv(i);
  }


  void BaseVector :: MemoryUsage (Array<MemoryUsageStruct*> & mu) const
  { 
    ;
  }

  /*
    BaseVector * BaseVector :: CreateVector ( const Array<int> * procs ) const
    {
    cout << "Create vec called for base class" << endl;
    return 0;
    }
  */

  void BaseVector :: SetRandom () 
  {
    FlatVector<double> fv = FVDouble();
    for (int i = 0; i < fv.Size(); i++)
      fv(i) = double (rand()) / RAND_MAX;
  }
  
  /*  
      template<int S>
      void BaseVector :: GetIndirect (const Array<int> & ind, 
      FlatVector< Vec<S,double> > & v) const 
      { 
      FlatVector<double> fv = FVDouble();
      if(EntrySize() != S)
      throw Exception("BaseVector::GetIndirect() wrong dimensions");

      //int ii = 0;
      for (int i = 0; i < ind.Size(); i++)
      if (ind[i] != -1)
      {
      int base = S * ind[i];
      for (int j = 0; j < S; j++)
      v[i](j) = fv[base++];
      }
      else
      {
      for (int j = 0; j < S; j++)
      v[i](j) = 0;
      }
      }

      template<int S>
      void BaseVector :: GetIndirect (const Array<int> & ind, 
      FlatVector< Vec<S,Complex> > & v) const 
      { 
      FlatVector<Complex> fv = FVComplex();
      if(EntrySize() != 2*S)
      throw Exception("BaseVector::GetIndirect() wrong dimensions");

      //int ii = 0;
      for (int i = 0; i < ind.Size(); i++)
      if (ind[i] != -1)
      {
      int base = S * ind[i];
      for (int j = 0; j < S; j++)
      v[i](j) = fv[base++];
      }
      else
      {
      for (int j = 0; j < S; j++)
      v[i](j) = 0;
      }
      }
  */




  template<>
  FlatVector<Complex> S_BaseVector<Complex> :: FVComplex () const
  {
    FlatVector<Complex> fv = FVScal();
    return FlatVector<Complex> (fv.Size() * sizeof(Complex)/sizeof(Complex),
                                reinterpret_cast<Complex*> (&fv(0)));
  }


  template<> double S_BaseVector<double> :: InnerProductD (const BaseVector & v2) const
  {
    return InnerProduct(v2);
  }
  template<> Complex S_BaseVector<double> :: InnerProductC (const BaseVector & v2) const
  {
    throw Exception ("InnerProductC for real vector");
  }

  template<> double S_BaseVector<Complex> :: InnerProductD (const BaseVector & v2) const
  {
    throw Exception ("InnerProductD for complex vector");
  }
  template<> Complex S_BaseVector<Complex> :: InnerProductC (const BaseVector & v2) const
  {
    return InnerProduct(v2);
  }




  template<>
  void S_BaseVector<double> :: GetIndirect (const FlatArray<int> & ind, 
					    const FlatVector<double> & v) const 
  { 
    FlatSysVector<double> lsv(Size(), EntrySize(), &FVDouble()(0));
    FlatSysVector<double> sv(ind.Size(), EntrySize(), &v(0));

    for (int i = 0; i < ind.Size(); i++)
      if (ind[i] != -1)
	sv(i) = lsv(ind[i]);
      else
	sv(i) = -1.0;
  }
  
  template<>
  void S_BaseVector<double> :: GetIndirect (const FlatArray<int> & ind, 
					    const FlatVector<Complex> & v) const 
  { 
    FlatSysVector<double> lsv(Size(), EntrySize(), &FVDouble()(0));
    FlatSysVector<Complex> sv(ind.Size(), EntrySize(), &v(0));

    for (int i = 0; i < ind.Size(); i++)
      if (ind[i] != -1)
	sv(i) = lsv(ind[i]);
      else
	sv(i) = -1.0;
    /*
    FlatVector<Complex> fv = FVComplex();
    int es = EntrySize() / 2;
    int ii = 0;
    for (int i = 0; i < ind.Size(); i++)
      if (ind[i] != -1)
	{
	  int base = es * ind[i];
	  for (int j = 0; j < es; j++)
	    v[ii++] = fv[base++];
	}
      else
	{
	  for (int j = 0; j < es; j++)
	    v[ii++] = 0;
	}
    */
  }
  

  template<>
  void S_BaseVector<Complex> :: GetIndirect (const FlatArray<int> & ind, 
					     const FlatVector<double> & v) const 
  { 
    throw Exception ("BaseVector<Complex>::GetIndirect<double> called");
  }
  
  template<>
  void S_BaseVector<Complex> :: GetIndirect (const FlatArray<int> & ind, 
					     const FlatVector<Complex> & v) const 
  { 
    FlatVector<Complex> fv = FVComplex();
    int es = EntrySize() / 2;
    int ii = 0;
    for (int i = 0; i < ind.Size(); i++)
      if (ind[i] != -1)
	{
	  int base = es * ind[i];
	  for (int j = 0; j < es; j++)
	    v[ii++] = fv[base++];
	}
      else
	{
	  for (int j = 0; j < es; j++)
	    v[ii++] = 0;
	}
  }
  







  void BaseVector :: SetIndirect (FlatArray<int> ind, 
				  FlatVector<double> v) 
  { 
    FlatSysVector<double> lsv(Size(), EntrySize(), &FVDouble()(0));
    FlatSysVector<double> sv(ind.Size(), EntrySize(), &v(0));

    for (int i = 0; i < ind.Size(); i++)
      if (ind[i] != -1)
	lsv(ind[i]) = sv(i);

    /*
      FlatVector<double> fv = FVDouble();
      int es = EntrySize();
      int ii = 0;
      for (int i = 0; i < ind.Size(); i++)
      if (ind[i] != -1)
      {
      int base = es * ind[i];
      for (int j = 0; j < es; j++)
      fv[base++] = v[ii++];
      }
      else
      ii += es;
    */
  }

  void BaseVector :: SetIndirect (FlatArray<int> ind, 
				  FlatVector<Complex> v) 
  { 
    FlatVector<Complex> fv = FVComplex();
    int es = EntrySize() / 2;
    int ii = 0;
    for (int i = 0; i < ind.Size(); i++)
      if (ind[i] != -1)
	{
	  int base = es * ind[i];
	  for (int j = 0; j < es; j++)
	    fv[base++] = v[ii++];
	}
      else
	ii += es;
  }

  /*
    template<int S>
    void BaseVector :: AddIndirect (const Array<int> & ind, 
    const FlatVector< Vec<S,double> > & v) 
    { 
    FlatVector<double> fv = FVDouble();
    int es = EntrySize();
    
    for (int i = 0; i < ind.Size(); i++)
    if (ind[i] != -1)
    {
    int base = es * ind[i];
    for (int j = 0; j < es; j++)
    fv[base++] += v[i](j);
    }
    }

    template<int S>
    void BaseVector :: AddIndirect (const Array<int> & ind, 
    const FlatVector< Vec<S,Complex> > & v)
    { 
    FlatVector<Complex> fv = FVComplex();
    if(EntrySize() != 2*S)
    throw Exception("BaseVector::AddIndirect() wrong dimensions");

    for (int i = 0; i < ind.Size(); i++)
    if (ind[i] != -1)
    {
    int base = S * ind[i];
    for (int j = 0; j < S; j++)
    fv[base++] += v[i](j);
    }
    }
  */  

  void BaseVector :: AddIndirect (FlatArray<int> ind, 
				  FlatVector<double> v) 
  { 
    FlatSysVector<double> lsv(Size(), EntrySize(), &FVDouble()(0));
    FlatSysVector<double> sv(ind.Size(), EntrySize(), &v(0));

    for (int i = 0; i < ind.Size(); i++)
      if (ind[i] != -1)
	lsv(ind[i]) += sv(i);
  }

  void BaseVector :: AddIndirect (FlatArray<int> ind, 
				  FlatVector<Complex> v)
  { 
    FlatVector<Complex> fv = FVComplex();
    int es = EntrySize() / 2;
    int ii = 0;
    for (int i = 0; i < ind.Size(); i++)
      if (ind[i] != -1)
	{
	  int base = es * ind[i];
	  for (int j = 0; j < es; j++)
	    fv[base++] += v[ii++];
	}
      else
	ii += es;
  }


  void BaseVector :: Cumulate () const { ; }
  void BaseVector :: Distribute() const { ; }
  PARALLEL_STATUS BaseVector :: GetParallelStatus () const { return NOT_PARALLEL; }
  void BaseVector :: SetParallelStatus (PARALLEL_STATUS stat) const { ; }
  



  /**
     Decision between double or Complex
  */



  template <class SCAL>
  S_BaseVector<SCAL> & S_BaseVector<SCAL> :: operator= (double s)
  {
    SetScalar (s);
    return *this;
  }


  template <class SCAL>
  SCAL S_BaseVector<SCAL> :: InnerProduct (const BaseVector & v2) const
  {
    static Timer t("S_BaseVector::InnerProduct");
    RegionTimer reg(t);

    return ngbla::InnerProduct (FVScal(), v2.FV<SCAL>());
    // dynamic_cast<const S_BaseVector&>(v2).FVScal());
  }


  double InnerProductRec (FlatVector<> a, FlatVector<> b)
  {
    int n = a.Size();
    if (n < 4096)
      {
	return InnerProduct (a, b);
      }

    int n1 = n/2;
    double sum1, sum2;
#pragma omp task shared(sum1) 
    {
      sum1 = InnerProductRec (a.Range(0,n1), b.Range(0,n1));
    }
    sum2 = InnerProductRec (a.Range(n1,n), b.Range(n1,n));

#pragma omp taskwait
    return sum1+sum2;
  }



  template <>
  double S_BaseVector<double> :: InnerProduct (const BaseVector & v2) const
  {
    if (omp_status == OMP_STATUS::TASKS)
      {
        static Timer t("S_BaseVector<double>::InnerProduct (TASKS)");
        RegionTimer reg(t);
        
        FlatVector<> me = FVDouble();
        FlatVector<> you = FVDouble();
        return InnerProductRec (me, you);
      }

    
    static Timer t("S_BaseVector<double>::InnerProduct");
    RegionTimer reg(t);
    return ngbla::InnerProduct (FVScal(), v2.FV<double>());
  }







  //template <> 
  /*
  Complex S_BaseVector<Complex> :: InnerProduct (const BaseVector & v2) const
  {
    return ngbla::InnerProduct (FVScal(), 
                                dynamic_cast<const S_BaseVector&>(v2).FVScal());
  }
  */


  template <class SCAL>
  FlatVector<double> S_BaseVector<SCAL> :: FVDouble () const 
  {
    return FlatVector<double> (size * entrysize, Memory());
    /*
      FlatVector<SCAL> fv = FVScal();
      return FlatVector<SCAL> (fv.Size() * sizeof(SCAL)/sizeof(double),
      reinterpret_cast<double*> (&fv(0)));
    */
  }

  template <class SCAL>
  FlatVector<Complex> S_BaseVector<SCAL> :: FVComplex () const
  {
    throw Exception ("FVComplex called for real vector");
  }



  template<>
  FlatVector<double> S_BaseVector<Complex> :: FVDouble () const
  {
    FlatVector<Complex> fv = FVScal();
    return FlatVector<double> (fv.Size() * sizeof(Complex)/sizeof(double),
                               reinterpret_cast<double*> (&fv(0)));
  }





  template <typename TSCAL>
  // shared_ptr<BaseVector> S_BaseVectorPtr<TSCAL> :: CreateVector () const
  AutoVector S_BaseVectorPtr<TSCAL> :: CreateVector () const
  {
    switch (es)
      {
      case 1: return make_shared<VVector<TSCAL>> (this->size);
      case 2: return make_shared<VVector<Vec<2,TSCAL>>> (this->size);
      case 3: return make_shared<VVector<Vec<3,TSCAL>>> (this->size);
      }
    return make_shared<S_BaseVectorPtr<TSCAL>> (this->size, es);
  }

  template <typename TSCAL>
  AutoVector S_BaseVectorPtr<TSCAL> :: Range (int begin, int end) const
  {
    return make_shared<S_BaseVectorPtr<TSCAL>> (end-begin, es, pdata+begin*es);
  }
  
  template <typename TSCAL>
  AutoVector S_BaseVectorPtr<TSCAL> :: Range (IntRange range) const
  {
    return make_shared<S_BaseVectorPtr<TSCAL>> (range.Size(), es, pdata+range.First()*es);
  }
  

  template class S_BaseVector<double>;
  // template class S_BaseVector<Complex>;
  
  template class VFlatVector<double>;
  
  template class S_BaseVectorPtr<double>;
  template class S_BaseVectorPtr<Complex>;


  template class VVector<double>;
  template class VVector<Complex>;

}
