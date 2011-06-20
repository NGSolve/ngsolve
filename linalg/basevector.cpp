/*********************************************************************/
/* File:   basevector.cpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   7. Feb. 2003                                              */
/*********************************************************************/

/* 
   base class in vector hierarchy
*/


#include <la.hpp>


namespace ngla
{
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


  BaseVector & BaseVector :: Scale (double scal)
  {
    FVDouble() *= scal;
    return *this;
  }
  BaseVector & BaseVector :: Scale (Complex scal)
  {
    FVComplex() *= scal;
    return *this;
  }

  BaseVector & BaseVector :: SetScalar (double scal)
  {
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
    FVDouble() = scal * v.FVDouble();
    return *this;
  }

  BaseVector & BaseVector :: Set (Complex scal, const BaseVector & v)
  {
    FVComplex() = scal * v.FVComplex();
    return *this;
  }
    
  BaseVector & BaseVector :: Add (double scal, const BaseVector & v)
  {
    FVDouble() += scal * v.FVDouble();
    return *this;
  }

  BaseVector & BaseVector :: Add (Complex scal, const BaseVector & v)
  {
    FVComplex() += scal * v.FVComplex();
    return *this;
  }


  BaseVector * BaseVector ::Range (int begin, int end) const
  {
    throw Exception ("BaseVector::Range const called");
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
      {
	ost << fv(i) << " ";
      }
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

  BaseVector * BaseVector :: CreateVector ( const Array<int> * procs ) const
  {
    cout << "Create vec called for base class" << endl;
    return 0;
  }


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

  void BaseVector :: GetIndirect (const FlatArray<int> & ind, 
				  const FlatVector<double> & v) const 
  { 
    FlatVector<double> fv = FVDouble();
    int es = EntrySize();
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
  
  void BaseVector :: GetIndirect (const FlatArray<int> & ind, 
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
  
  void BaseVector :: SetIndirect (const FlatArray<int> & ind, 
				  const FlatVector<double> & v) 
  { 
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
  }

  void BaseVector :: SetIndirect (const FlatArray<int> & ind, 
				  const FlatVector<Complex> & v) 
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

  void BaseVector :: AddIndirect (const FlatArray<int> & ind, 
				  const FlatVector<double> & v) 
  { 
    FlatVector<double> fv = FVDouble();
    int es = EntrySize();
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

  void BaseVector :: AddIndirect (const FlatArray<int> & ind, 
				  const FlatVector<Complex> & v)
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
  

  /*
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<2,double> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<3,double> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<4,double> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<5,double> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<6,double> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<7,double> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<8,double> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<9,double> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<10,double> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<11,double> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<12,double> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<13,double> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<14,double> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<15,double> > & v) const;

  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<2,Complex> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<3,Complex> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<4,Complex> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<5,Complex> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<6,Complex> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<7,Complex> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<8,Complex> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<9,Complex> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<10,Complex> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<11,Complex> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<12,Complex> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<13,Complex> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<14,Complex> > & v) const;
  template void BaseVector::GetIndirect(const Array<int> & ind, 
					   FlatVector< Vec<15,Complex> > & v) const;




  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<2,double> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<3,double> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<4,double> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<5,double> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<6,double> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<7,double> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<8,double> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<9,double> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<10,double> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<11,double> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<12,double> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<13,double> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<14,double> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<15,double> > & v);

  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<2,Complex> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<3,Complex> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<4,Complex> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<5,Complex> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<6,Complex> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<7,Complex> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<8,Complex> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<9,Complex> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<10,Complex> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<11,Complex> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<12,Complex> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<13,Complex> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<14,Complex> > & v);
  template void BaseVector::AddIndirect(const Array<int> & ind, 
					   const FlatVector< Vec<15,Complex> > & v);
  */






/**
   Decision between double or Complex
 */

/*
template <class SCAL>
S_BaseVector<SCAL> :: S_BaseVector () throw()
{ 
  ;
}

template <class SCAL>
S_BaseVector<SCAL> :: ~S_BaseVector () throw()
{ 
  ;
}
*/


template <class SCAL>
S_BaseVector<SCAL> & S_BaseVector<SCAL> :: operator= (double s)
{
  SetScalar (s);
  return *this;
}

/*
template <class SCAL>
SCAL S_BaseVector<SCAL> :: InnerProduct (const BaseVector & v2) const
{
  throw Exception ("Inner Product called for S_BaseVector");
}
*/

template <class SCAL>
FlatVector<double> S_BaseVector<SCAL> :: FVDouble () const 
{
  return FlatVector<SCAL> (size * entrysize, Memory());
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




FlatVector<double> S_BaseVector<Complex> :: FVDouble () const throw()
{
  FlatVector<Complex> fv = FVScal();
  return FlatVector<double> (fv.Size() * sizeof(Complex)/sizeof(double),
			     reinterpret_cast<double*> (&fv(0)));
}


FlatVector<Complex> S_BaseVector<Complex> :: FVComplex () const throw()
{
  FlatVector<Complex> fv = FVScal();
  return FlatVector<Complex> (fv.Size() * sizeof(Complex)/sizeof(Complex),
			      reinterpret_cast<Complex*> (&fv(0)));
}




template class S_BaseVector<double>;
//SZ	template class S_BaseVector<Complex>;


}
