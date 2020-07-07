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

#ifdef PARALLEL
#include "../parallel/parallelvector.hpp"
#endif


namespace ngla
{

  using namespace ngbla;
  
#ifdef PARALLEL
  class ParallelBaseVector;
#endif


  unique_ptr<MultiVector> BaseVector :: CreateMultiVector (size_t cnt) const
  {
    return make_unique<MultiVector> (CreateVector(), cnt);
  }

  
  double BaseVector :: L2Norm () const
  {
    static Timer t("BaseVector::L2Norm");
    RegionTimer reg(t);

    auto me = FVDouble();
    t.AddFlops (me.Size());

    /*
    atomic<double> sum(0.0);
    ParallelForRange ( me.Range(),
                       [&] (IntRange r) 
                       {
                         double mysum = ngbla::L2Norm2 (me.Range(r));
                         sum += mysum;
                       });
    */
    double parts[16];
    ParallelJob ([me,&parts] (TaskInfo ti)
                 {
                   auto r = ngstd::Range(me).Split (ti.task_nr, ti.ntasks);
                   parts[ti.task_nr] = ngbla::L2Norm2 (me.Range(r));
                 }, 16);
    double sum = 0;
    for (double part : parts) sum += part;
    return sqrt(double(sum));
  }

  BaseVector & BaseVector :: Scale (double scal)
  {
    if (scal == 1) return *this;

    auto me = FVDouble();

    static Timer t("BaseVector::Scale");
    RegionTimer reg(t);
    t.AddFlops (me.Size());

    ParallelFor ( me.Range(),
                  [me,scal] (size_t i) { me(i) *= scal; });

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
    
    auto fv = FVDouble();
    t.AddFlops (fv.Size());

    ParallelFor ( fv.Range(),
                  [fv,scal] (size_t i) { fv(i) = scal; },
                  TasksPerThread(1), TotalCosts(fv.Size())
                  );
    
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


    auto me = FVDouble();
    auto you = v.FVDouble();

    if (&me(0) == &you(0) && scal==1.0) return *this;
    
    t.AddFlops (me.Size());

    ParallelFor ( me.Range(),
                  [me,you,scal] (size_t i) { me(i) = scal * you(i); });
    
    return *this;
  }

  BaseVector & BaseVector :: Set (Complex scal, const BaseVector & v)
  {
    if (Size() != v.Size())
      throw Exception (string ("BaseVector::Set: size of me = ") +
                       ToString(Size()) + " != size of other = " + ToString(v.Size()));

    if (v.IsComplex())
      FVComplex() = scal * v.FVComplex();
    else
      FVComplex() = scal * v.FVDouble();      
    return *this;
  }
    


  BaseVector & BaseVector :: Add (double scal, const BaseVector & v)
  {
    static Timer t("BaseVector::Add");
    RegionTimer reg(t);
    
    auto me = FVDouble();
    auto you = v.FVDouble();

    if (me.Size() != you.Size())
      throw Exception (string ("BaseVector::Add: size of me = ") +
                       ToString(Size()) + " != size of other = " + ToString(v.Size()));
    
    t.AddFlops (me.Size());

    ParallelFor (me.Range(),
                 [me,you,scal] (size_t i) { me(i) += scal * you(i); });
    
    return *this;
  }

  BaseVector & BaseVector :: Add (Complex scal, const BaseVector & v)
  {
    if(Size() != v.Size())
      throw Exception (string ("BaseVector::Add: size of me = ") +
                       ToString(Size()) + " != size of other = " + ToString(v.Size()));

    if (v.IsComplex())
      FVComplex() += scal * v.FVComplex();
    else
      FVComplex() += scal * v.FVDouble();      
    return *this;
  }


  double BaseVector :: InnerProductD (const BaseVector & v2) const
  {
    return dynamic_cast<const S_BaseVector<double>&> (*this) . 
      InnerProduct (v2);
  }
  
  Complex BaseVector :: InnerProductC (const BaseVector & v2, bool conjugate) const
  {
    return dynamic_cast<const S_BaseVector<Complex>&> (*this) . 
      InnerProduct (v2, conjugate);
  }


  AutoVector BaseVector ::Range (size_t begin, size_t end) const
  {
    throw Exception ("BaseVector::Range const called");
  }

  AutoVector BaseVector :: Range (T_Range<size_t> range) const
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
    for (size_t i = 0; i < fv.Size(); i++)
      SaveBin (ost, fv(i));
  }

  void BaseVector :: Load(istream & ist) 
  {
    FlatVector<double> fv = FVDouble();
    for (size_t i = 0; i < fv.Size(); i++)
      LoadBin (ist, fv(i));
  }

  void BaseVector :: SaveText(ostream & ost) const
  {
    FlatVector<double> fv = FVDouble();
    for (size_t i = 0; i < fv.Size(); i++)
      ost << fv(i) << " ";
  }

  void BaseVector :: LoadText(istream & ist) 
  {
    FlatVector<double> fv = FVDouble();
    for (size_t i = 0; i < fv.Size(); i++)
      ist >> fv(i);
  }

  size_t BaseVector :: CheckSum () const
  {
    size_t sum = 0;
    auto fv = FVDouble();
    for (auto i : ngstd::Range(fv))
      {
        double val = fv(i);
        sum += *reinterpret_cast<size_t*> (&val);
      }
    return sum;
  }

  Array<MemoryUsage> BaseVector :: GetMemoryUsage () const
  { 
    return Array<MemoryUsage>();
  }

  void BaseVector :: SetRandom () 
  {
    FlatVector<double> fv = FVDouble();
    for (size_t i = 0; i < fv.Size(); i++)
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

  AutoVector CreateBaseVector(size_t size, bool is_complex, int es)
  {
    shared_ptr<BaseVector> res;
    if(es > 1)
      {
        if(is_complex)
          res = make_shared<S_BaseVectorPtr<Complex>> (size, es);
        else
          res = make_shared<S_BaseVectorPtr<double>> (size, es);
        return res;
      }
    
    if (is_complex)
      res = make_shared<VVector<Complex>> (size);
    else
      res = make_shared<VVector<double>> (size);
    return res;
  }


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
  template<> Complex S_BaseVector<double> :: InnerProductC (const BaseVector & v2, bool conjugate) const
  {
    throw Exception ("InnerProductC for real vector");
  }

  template<> double S_BaseVector<Complex> :: InnerProductD (const BaseVector & v2) const
  {
    throw Exception ("InnerProductD for complex vector");
  }
  template<> Complex S_BaseVector<Complex> :: InnerProductC (const BaseVector & v2, bool conjugate) const
  {
    return InnerProduct(v2, conjugate);
  }




  template<>
  void S_BaseVector<double> :: GetIndirect (FlatArray<int> ind, 
                                            FlatVector<double> v) const 
  {
    if (EntrySize() == 1)
      {
        FlatVector<double> lsv(Size(), &FVDouble()(0));
        for (auto i : ind.Range())
          {
            int index = ind[i];
            v(i) = IsRegularIndex(index) ? lsv(index) : 0;
          }
        /*
        int i = 0;
        double temp[8];
        for ( ; i + 7 < ind.Size(); i+=8)
          {
            for (int i2 = 0; i2 < 8; i2++)
              {
                int index = ind[i+i2];
                temp[i2] = (index != -1) ? lsv(index) : 0;
              }
            for (int i2 = 0; i2 < 8; i2++)
              v(i+i2) = temp[i2];
          }
        for ( ; i < ind.Size(); i++)
          {
            int index = ind[i];
            v(i) = (index != -1) ? lsv(index) : 0;
          }
        */
      }
    else
      {
        FlatSysVector<double> lsv(Size(), EntrySize(), &FVDouble()(0));
        FlatSysVector<double> sv(ind.Size(), EntrySize(), &v(0));
        
        for (size_t i = 0; i < ind.Size(); i++)
          if (IsRegularIndex(ind[i]))
            sv(i) = lsv(ind[i]);
          else
            sv(i) = -1.0;
      }
  }
  
  template<>
  void S_BaseVector<double> :: GetIndirect (FlatArray<int> ind, 
					    FlatVector<Complex> v) const 
  { 
    FlatSysVector<double> lsv(Size(), EntrySize(), &FVDouble()(0));
    FlatSysVector<Complex> sv(ind.Size(), EntrySize(), &v(0));

    for (size_t i = 0; i < ind.Size(); i++)
      if (IsRegularIndex(ind[i]))
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
  void S_BaseVector<Complex> :: GetIndirect (FlatArray<int> ind, 
                                             FlatVector<double> v) const 
  { 
    throw Exception ("BaseVector<Complex>::GetIndirect<double> called");
  }
  
  template<>
  void S_BaseVector<Complex> :: GetIndirect (FlatArray<int> ind, 
                                             FlatVector<Complex> v) const 
  { 
    FlatVector<Complex> fv = FVComplex();
    int es = EntrySize() / 2;
    int ii = 0;
    for (size_t i = 0; i < ind.Size(); i++)
      if (IsRegularIndex(ind[i]))
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

    for (size_t i = 0; i < ind.Size(); i++)
      if (IsRegularIndex(ind[i]))
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
    for (size_t i = 0; i < ind.Size(); i++)
      if (IsRegularIndex(ind[i]))
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
				  FlatVector<double> v, bool use_atomic) 
  {
    if (EntrySize() == 1)
      {
        FlatVector<double> lsv(Size(), &FVDouble()(0));

        if (!use_atomic)
          {
            for (size_t i = 0; i < ind.Size(); i++)
              if (IsRegularIndex(ind[i]))
                lsv(ind[i]) += v(i);
          }
        else
          {
            for (size_t i = 0; i < ind.Size(); i++)
              if (IsRegularIndex(ind[i]))
                AtomicAdd (lsv(ind[i]), v(i));
            // lsv(ind[i]) += v(i);
          }
      }
    else
      {
        FlatSysVector<double> lsv(Size(), EntrySize(), &FVDouble()(0));
        FlatSysVector<double> sv(ind.Size(), EntrySize(), &v(0));
        
        for (size_t i = 0; i < ind.Size(); i++)
          if (IsRegularIndex(ind[i]))
            lsv(ind[i]) += sv(i);
      }
  }

  void BaseVector :: AddIndirect (FlatArray<int> ind, 
				  FlatVector<Complex> v, bool use_atomic)
  { 
    FlatVector<Complex> fv = FVComplex();
    int es = EntrySize() / 2;

    if (es == 1)
      {
        if (!use_atomic)
          {
            for (size_t i = 0; i < ind.Size(); i++)
              if (IsRegularIndex(ind[i]))
                fv(ind[i]) += v(i);
          }
        else
          {
            for (size_t i = 0; i < ind.Size(); i++)
              if (IsRegularIndex(ind[i]))
                AtomicAdd (fv(ind[i]), v(i));
          }
      }
    else
      {
    
        int ii = 0;
        for (size_t i = 0; i < ind.Size(); i++)
          if (IsRegularIndex(ind[i]))
            {
              int base = es * ind[i];
              for (int j = 0; j < es; j++)
                fv[base++] += v[ii++];
            }
          else
            ii += es;
      }
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
  BaseVector & S_BaseVector<SCAL> :: SetScalar (double scal)
  {
    static Timer t("S_BaseVector::SetScalar");
    RegionTimer reg(t);
    
    auto me = FVScal();
    ParallelForRange (me.Size(),
                      [me, scal] (IntRange r) { me.Range(r) = scal; });
    
    // FVScal() = scal;
    return *this;
  }
  
  template <class SCAL>
  SCAL S_BaseVector<SCAL> :: InnerProduct (const BaseVector & v2, bool conjugate) const
  {
    static Timer t("S_BaseVector::InnerProduct");
    RegionTimer reg(t);

    if (conjugate)
      return ngbla::InnerProduct (Conj(FVScal()), v2.FV<SCAL>());
    else
      return ngbla::InnerProduct (FVScal(), v2.FV<SCAL>());
    // dynamic_cast<const S_BaseVector&>(v2).FVScal());
  }



  template <>
  double S_BaseVector<double> :: InnerProduct (const BaseVector & v2, bool conjugate) const
  {
    static Timer t("BaseVector::InnerProduct (taskhandler)");
    RegionTimer reg(t);

    auto me = FVDouble();
    auto you = v2.FVDouble();
	
    t.AddFlops (me.Size());
    /*
    atomic<double> scal(0);

    ParallelForRange (ngstd::Range(me.Size()),
                      [me,you,&scal] (IntRange r)
                      {
                        double myscal = ngbla::InnerProduct (me.Range(r), you.Range(r));
                        scal += myscal;
                      } );
    */
    double parts[16];
    ParallelJob ([me,you,&parts] (TaskInfo ti)
                 {
                   auto r = ngstd::Range(me).Split (ti.task_nr, ti.ntasks);
                   parts[ti.task_nr] =  ngbla::InnerProduct (me.Range(r), you.Range(r));
                 }, 16);
    double scal = 0;
    for (double part : parts) scal += part;
    return scal;
  }




  template <class SCAL>
  FlatVector<double> S_BaseVector<SCAL> :: FVDouble () const 
  {
    return FlatVector<double> (size * entrysize, Memory());
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



  BlockVector & dynamic_cast_BlockVector (BaseVector & x)
  {
    AutoVector * ax = dynamic_cast<AutoVector*> (&x);
    if (ax) return dynamic_cast<BlockVector&> (**ax);
    return dynamic_cast<BlockVector&> (x);
  }
  
  const BlockVector & dynamic_cast_BlockVector (const BaseVector & x)
  {
    const AutoVector * ax = dynamic_cast<const AutoVector*> (&x);
    if (ax) return dynamic_cast<const BlockVector&> (**ax);
    return dynamic_cast<const BlockVector&> (x);
  }


  BlockVector :: BlockVector (const Array<shared_ptr<BaseVector>> & avecs)
    : vecs(avecs), ispar(vecs.Size())
  {
    size = 0;
    for (auto & vec:vecs)
      size += vec->Size();
#ifdef PARALLEL
    ispar.Clear();
    for (size_t k = 0; k<vecs.Size(); k++) {
      auto stat = vecs[k]->GetParallelStatus();
      if ( stat==NOT_PARALLEL ) continue;
      ispar.SetBit(k);
      auto * pv = dynamic_cast_ParallelBaseVector(vecs[k].get());
      comm = pv->GetParallelDofs()->GetCommunicator();
      /*
      auto vcomm = pv->GetParallelDofs()->GetCommunicator();
      if (comm==MPI_COMM_NULL)
	comm = vcomm;
      else if (comm != vcomm)
	throw Exception("Tried to construct a BlockVector with components in different MPI-Communicators!!");
      */
    }
#endif
  }

  void * BlockVector :: Memory () const
  { throw Exception("BlockVector::Memory is not useful"); }
  FlatVector<double> BlockVector :: FVDouble () const 
  { throw Exception("BlockVector::FVDouble is not useful"); }
  FlatVector<Complex> BlockVector :: FVComplex () const
  { throw Exception("BlockVector::FVComplex is not useful"); }
  void BlockVector :: GetIndirect (FlatArray<int> ind, 
                                   FlatVector<double> v) const
  { throw Exception("BlockVector::GetIndirect is not useful"); }      
  void BlockVector :: GetIndirect (FlatArray<int> ind, 
                                   FlatVector<Complex> v) const 
  { throw Exception("BlockVector::GetIndirect is not useful"); }      
  
  // not yet implemented properly for complex components!!
  bool BlockVector :: IsComplex() const
  { return vecs[0]->IsComplex(); }
  
  AutoVector BlockVector :: CreateVector () const
  {
    Array<shared_ptr<BaseVector>> v2;
    for (auto & v : vecs)
      v2 += v->CreateVector();
    return make_shared<BlockVector> (v2);
  }
  
  double BlockVector :: InnerProductD (const BaseVector & v2) const
  {
    double pp = 0;
    double ps = 0;
    const auto & v2b = dynamic_cast_BlockVector(v2);
    for (size_t k = 0; k<vecs.Size(); k++) {
      auto p = vecs[k]->InnerProductD(*v2b[k]);
      if (ispar.Test(k)) pp += p;
      else ps += p;
    }
    // if all vectors are sequential, do not reduce reduce
    // if (MPI_Comm(comm) == MPI_COMM_NULL) return ps;
    return pp + comm.AllReduce(ps, MPI_SUM);
  }

  Complex BlockVector :: InnerProductC (const BaseVector & v2,
                                        bool conjugate) const
  {
    Complex pp = 0;
    Complex ps = 0;
    const auto & v2b = dynamic_cast_BlockVector(v2);
    for (size_t k = 0; k<vecs.Size(); k++) {
      auto p = vecs[k]->InnerProductC(*v2b[k], conjugate);
      if (ispar.Test(k)) pp += p;
      else ps += p;
    }
    // if all vectors are sequential, do not reduce reduce
    // if (MPI_Comm(comm) == MPI_COMM_NULL) return ps;
    return pp + comm.AllReduce(ps, MPI_SUM);
  }


  double BlockVector :: L2Norm () const
  {
    double sum = 0;
    for (size_t k = 0; k < vecs.Size(); k++)
      sum += sqr(vecs[k]->L2Norm());
    return sqrt(sum);
  }

  
  BaseVector & BlockVector :: Scale (double scal)
  {
    for (auto i : ngstd::Range(vecs))
      *vecs[i] *= scal;
    return *this;
  }

  BaseVector & BlockVector :: Scale (Complex scal)
  {
    for(auto i : ngstd::Range(vecs))
      vecs[i]->Scale(scal);
    return *this;
  }
  
  BaseVector & BlockVector :: SetScalar (double scal)
  {
    for (auto i : ngstd::Range(vecs))
      vecs[i]->SetScalar(scal);
    return *this;
  }

  BaseVector & BlockVector :: SetScalar (Complex scal)
  {
    for (auto i : ngstd::Range(vecs))
      vecs[i]->SetScalar(scal);
    return *this;
  }

  ostream & BlockVector :: Print (ostream & ost) const
  {
    for (auto i : ngstd::Range(vecs))
      vecs[i] -> Print(ost);
    return ost;
  }
  
  BaseVector & BlockVector :: Set (double scal, const BaseVector & v)
  {
    auto & bv = dynamic_cast_BlockVector(v);
    for (size_t i : ngstd::Range(vecs))
      vecs[i] -> Set(scal, *bv[i]);
    return *this;
  }

  BaseVector & BlockVector :: Set(Complex scal, const BaseVector & v)
  {
    auto & bv = dynamic_cast_BlockVector(v);
    for (size_t i : ngstd::Range(vecs))
      vecs[i] -> Set(scal, *bv[i]);
    return *this;
  }

  BaseVector & BlockVector :: Add (double scal, const BaseVector & v)
  {
    auto & bv = dynamic_cast_BlockVector(v);
    for (size_t i : ngstd::Range(vecs))
      vecs[i] -> Add(scal, *bv[i]);
    return *this;
  }

  BaseVector & BlockVector :: Add (Complex scal, const BaseVector & v)
  {
    auto & bv = dynamic_cast_BlockVector(v);
    for (size_t i : ngstd::Range(vecs))
      vecs[i] -> Add(scal, *bv[i]);
    return *this;
  }
  
  template <typename TSCAL>
  S_BaseVectorPtr<TSCAL> :: ~S_BaseVectorPtr ()
  {
    if (ownmem) delete [] pdata;
  }

  template <typename TSCAL>
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
  ostream & S_BaseVectorPtr<TSCAL> :: Print (ostream & ost) const 
  {
    if (es == 1)
      ost << FlatVector<TSCAL> (this->size, pdata) << endl;
    else
      ost << FlatSysVector<TSCAL> (this->size, es, pdata);
    return ost;
  }
  
  
  template <typename TSCAL>
  Array<MemoryUsage> S_BaseVectorPtr<TSCAL> :: GetMemoryUsage () const
  {
    if (ownmem)
      return { { "Vector", sizeof(TSCAL)*es*this->size, 1 } };
    else
      return Array<MemoryUsage>();
  }
  
  

  template <typename TSCAL>
  AutoVector S_BaseVectorPtr<TSCAL> :: Range (size_t begin, size_t end) const
  {
    return make_shared<S_BaseVectorPtr<TSCAL>> (end-begin, es, pdata+begin*es);
  }
  
  template <typename TSCAL>
  AutoVector S_BaseVectorPtr<TSCAL> :: Range (T_Range<size_t> range) const
  {
    return make_shared<S_BaseVectorPtr<TSCAL>> (range.Size(), es, pdata+range.First()*es);
  }

  template class S_BaseVector<double>;
  template class S_BaseVector<Complex>;
  
  template class VFlatVector<double>;
  
  template class S_BaseVectorPtr<double>;
  template class S_BaseVectorPtr<Complex>;

  template class VVector<double>;
  template class VVector<Complex>;
}
