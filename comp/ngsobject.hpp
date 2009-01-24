#ifndef FILE_NGSOBJECT
#define FILE_NGSOBJECT

/*********************************************************************/
/* File:   ngsobject.hh                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   2. Aug. 2000                                              */
/*********************************************************************/

/** 
    NGSolve base class
*/

class NGS_Object
{
protected:
  string name;
  Flags flaglist; // define flags
  const MeshAccess & ma;
  double time;

  /// profiling
  int timer;
public:
  /// 
  NGS_Object (const MeshAccess & ama, const string & aname = "noname", bool parseflags=false)
  : name(aname), ma(ama), time(0), skipCleanUp(0)
  { 
    timer = NgProfiler::CreateTimer (aname);
  }
  
  NGS_Object (const NGS_Object& obj)
  : name(obj.name), ma(obj.ma), time(obj.time), timer(obj.timer), 
    skipCleanUp(obj.skipCleanUp), flaglist(obj.flaglist)  
  { ; }
  
  virtual ~NGS_Object () { ; }
  ///
  void SetName (const string & aname)
  { 
    name = aname; 
    NgProfiler::SetName (timer, name); 
  }

  ///
  const string & GetName () const
  { return name; }

  ///
  const MeshAccess & GetMeshAccess() const
  { return ma; }

  virtual string GetClassName () const
  {
    return typeid(*this).name();
  }

  virtual void PrintReport (ostream & ost)
  {
    ost << typeid(*this).name();
  }

  virtual void MemoryUsage (ARRAY<MemoryUsageStruct*> & mu) const
  {
    cout << "MemoryUsage not implemented for class " << GetClassName() << endl;
  }

  double GetTime () const { return time; }
  void SetTime (double t) { time = t; }

  int GetTimer () const { return timer; }

  bool skipCleanUp; 
  
  void DefineStringFlag(const char* s, const char* val="");
  void DefineNumFlag(const char* s, double val=0);
  void DefineDefineFlag(const char* s);
  void DefineStringListFlag(const char* s);
  void DefineNumListFlag(const char* s);
  int ParseFlags(const Flags& flags); 
};



#ifdef OLD

#define CreateVecObject1S(dest, Object, dim, SCAL, arg) \
switch (dim) {  \
 case 1: dest = new Object<SCAL>(arg); break; \
 case 2: dest = new Object<Vec<2,SCAL> >(arg); break; \
 case 3: dest = new Object<Vec<3,SCAL> >(arg); break; \
 case 4: dest = new Object<Vec<4,SCAL> >(arg); break; \
 case 5: dest = new Object<Vec<5,SCAL> >(arg); break; \
 case 6: dest = new Object<Vec<6,SCAL> >(arg); break; \
 case 7: dest = new Object<Vec<7,SCAL> >(arg); break; \
 case 8: dest = new Object<Vec<8,SCAL> >(arg); break; \
 case 9: dest = new Object<Vec<9,SCAL> >(arg); break; \
 case 12: dest = new Object<Vec<12,SCAL> >(arg); break; \
 case 18: dest = new Object<Vec<18,SCAL> >(arg); break; \
 case 24: dest = new Object<Vec<24,SCAL> >(arg); break; \
}

#define CreateVecObject1(dest, Object, dim, iscomplex, arg) \
if (iscomplex) \
  CreateVecObject1S (dest, Object, dim, Complex, arg) \
else \
 CreateVecObject1S (dest, Object, dim, double, arg);




#define CreateVecObject2S(dest, Object, dim, SCAL, arg, arg2) \
switch (dim) {  \
 case 1: dest = new Object<SCAL>(arg, arg2); break; \
 case 2: dest = new Object<Vec<2,SCAL> >(arg, arg2); break; \
 case 3: dest = new Object<Vec<3,SCAL> >(arg, arg2); break; \
 case 4: dest = new Object<Vec<4,SCAL> >(arg, arg2); break; \
 case 5: dest = new Object<Vec<5,SCAL> >(arg, arg2); break; \
 case 6: dest = new Object<Vec<6,SCAL> >(arg, arg2); break; \
 case 7: dest = new Object<Vec<7,SCAL> >(arg, arg2); break; \
 case 8: dest = new Object<Vec<8,SCAL> >(arg, arg2); break; \
 case 9: dest = new Object<Vec<9,SCAL> >(arg, arg2); break; \
 case 12: dest = new Object<Vec<12,SCAL> >(arg, arg2); break; \
 case 18: dest = new Object<Vec<18,SCAL> >(arg, arg2); break; \
 case 24: dest = new Object<Vec<24,SCAL> >(arg, arg2); break; \
}

#define CreateVecObject2(dest, Object, dim, iscomplex, arg, arg2) \
if (iscomplex) \
  CreateVecObject2S (dest, Object, dim, Complex, arg, arg2) \
else \
  CreateVecObject2S (dest, Object, dim, double, arg, arg2);




#ifdef OLD

#define CreateVecObject3S(dest, Object, dim, SCAL, arg, arg2, arg3) \
switch (dim) {  \
 case 1: dest = new Object<SCAL>(arg, arg2, arg3); break; \
 case 2: dest = new Object<Vec<2,SCAL> >(arg, arg2, arg3); break; \
 case 3: dest = new Object<Vec<3,SCAL> >(arg, arg2, arg3); break; \
 case 4: dest = new Object<Vec<4,SCAL> >(arg, arg2, arg3); break; \
 case 5: dest = new Object<Vec<5,SCAL> >(arg, arg2, arg3); break; \
 case 6: dest = new Object<Vec<6,SCAL> >(arg, arg2, arg3); break; \
 case 7: dest = new Object<Vec<7,SCAL> >(arg, arg2, arg3); break; \
 case 8: dest = new Object<Vec<8,SCAL> >(arg, arg2, arg3); break; \
 case 9: dest = new Object<Vec<9,SCAL> >(arg, arg2, arg3); break; \
 case 12: dest = new Object<Vec<12,SCAL> >(arg, arg2, arg3); break; \
 case 18: dest = new Object<Vec<18,SCAL> >(arg, arg2, arg3); break; \
 case 24: dest = new Object<Vec<24,SCAL> >(arg, arg2, arg3); break; \
}

#define CreateVecObject3(dest, Object, dim, iscomplex, arg, arg2, arg3) \
if (iscomplex) \
  CreateVecObject3S (dest, Object, dim, Complex, arg, arg2, arg3) \
else \
  CreateVecObject3S (dest, Object, dim, double, arg, arg2, arg3);

#endif

#endif








template <template <class T> class Object, class Base, class SCAL, class ARG, class ARG2, class ARG3, int ACTDIM>
class TCreateVecObject3S { 
public:
  static Base * Create (int dim, ARG & arg, ARG2 & arg2, ARG3 & arg3)
  {
    if (dim == ACTDIM) return new Object<ngbla::Vec<ACTDIM,SCAL> > (arg, arg2, arg3);
    else return TCreateVecObject3S<Object, Base, SCAL, ARG, ARG2, ARG3, ACTDIM-1>::Create(dim, arg, arg2, arg3);
  }
};

template <template <class T> class Object, class Base, class SCAL, class ARG, class ARG2, class ARG3>
class TCreateVecObject3S<Object, Base, SCAL, ARG, ARG2, ARG3, 1> {
public:
  static Base * Create (int dim, ARG & arg, ARG2 & arg2, ARG3 & arg3)
  { 
    if (dim == 1) return new Object<SCAL> (arg, arg2, arg3);

    stringstream err;
    err << "illegal CreateVecObject3, dim = " << dim << endl;
    throw Exception (err.str());
  }
};

template <template <class T> class Object, class Base, class ARG, class ARG2, class ARG3>
Base * CreateVecObject (int dim, bool iscomplex, ARG & arg, ARG2 & arg2, ARG3 & arg3)
{
  if (!iscomplex)
    return TCreateVecObject3S<Object, Base, double, ARG, ARG2, ARG3, 24>::Create (dim, arg, arg2, arg3);
  else
    return TCreateVecObject3S<Object, Base, Complex, ARG, ARG2, ARG3, 24>::Create (dim, arg, arg2, arg3);
}










template <template <class T> class Object, class Base, class SCAL, class ARG, int ACTDIM>
class TCreateVecObjectS {
  Base * Create (int dim, ARG & arg)
  {
    if (dim == ACTDIM) return new Object<ngbla::Vec<ACTDIM,SCAL> > (arg);
    else return TCreateVecObjectS<Object, Base, SCAL, ARG, ACTDIM-1>::Create(dim, arg);
  }
};

template <template <class T> class Object, class Base, class SCAL, class ARG>
class TCreateVecObjectS<Object, Base, SCAL, ARG, 1> {
  Base * Create (int dim, ARG & arg)
  { 
    if (dim == 1) return new Object<SCAL> (arg);

    stringstream err;
    err << "illegal CreateVecObject, dim = " << dim << endl;
    throw Exception (err.str());
  }
};
/*
template <template <class T> class Object, class Base, class SCAL, class ARG>
Base * CreateVecObjectS (int dim, ARG & arg)
{
  switch (dim)
    {
    case 1: return new Object<SCAL> (arg);
    case 2: return new Object<ngbla::Vec<2,SCAL> > (arg);
    case 3: return new Object<ngbla::Vec<3,SCAL> > (arg);
    case 4: return new Object<ngbla::Vec<4,SCAL> > (arg);
    case 5: return new Object<ngbla::Vec<5,SCAL> > (arg);
    case 6: return new Object<ngbla::Vec<6,SCAL> > (arg);
    case 7: return new Object<ngbla::Vec<7,SCAL> > (arg);
    case 8: return new Object<ngbla::Vec<8,SCAL> > (arg);
    case 9: return new Object<ngbla::Vec<9,SCAL> > (arg);
    case 12: return new Object<ngbla::Vec<12,SCAL> > (arg);
    case 18: return new Object<ngbla::Vec<18,SCAL> > (arg);
    case 24: return new Object<ngbla::Vec<24,SCAL> > (arg);
    }

  stringstream err;
  err << "illegal CreateVecObject, dim = " << dim << endl;
  throw Exception (err.str());
}
*/

template <template <class T> class Object, class Base, class ARG>
Base * CreateVecObject (int dim, bool iscomplex, ARG & arg)
{
  if (!iscomplex)
    return TCreateVecObjectS<Object, Base, double, ARG, 12>::Create (dim, arg);
  else
    return TCreateVecObjectS<Object, Base, Complex, ARG, 12>::Create (dim, arg);
}

/*    
template <template <class T> class Object, class Base, class ARG>
Base * CreateVecObject (int dim, bool iscomplex, ARG & arg)
{
  if (!iscomplex)
    {
      switch (dim)
	{
	case 1: return new Object<double> (arg);
	case 2: return new Object<ngbla::Vec<2> > (arg);
	case 3: return new Object<ngbla::Vec<3> > (arg);
	case 4: return new Object<ngbla::Vec<4> > (arg);
	case 5: return new Object<ngbla::Vec<5> > (arg);
	case 6: return new Object<ngbla::Vec<6> > (arg);
	case 7: return new Object<ngbla::Vec<7> > (arg);
	case 8: return new Object<ngbla::Vec<8> > (arg);
	case 9: return new Object<ngbla::Vec<9> > (arg);
	case 12: return new Object<ngbla::Vec<12> > (arg);
	case 18: return new Object<ngbla::Vec<18> > (arg);
	case 24: return new Object<ngbla::Vec<24> > (arg);
	}
    }
  else
    {
      switch (dim)
	{
	case 1: return new Object<Complex> (arg);
	case 2: return new Object<ngbla::Vec<2,Complex> > (arg);
	case 3: return new Object<ngbla::Vec<3,Complex> > (arg);
	case 4: return new Object<ngbla::Vec<4,Complex> > (arg);
	case 5: return new Object<ngbla::Vec<5,Complex> > (arg);
	case 6: return new Object<ngbla::Vec<6,Complex> > (arg);
	case 7: return new Object<ngbla::Vec<7,Complex> > (arg);
	case 8: return new Object<ngbla::Vec<8,Complex> > (arg);
	case 9: return new Object<ngbla::Vec<9,Complex> > (arg);
	case 10: return new Object<ngbla::Vec<10,Complex> > (arg);
	case 11: return new Object<ngbla::Vec<11,Complex> > (arg);
	case 12: return new Object<ngbla::Vec<12,Complex> > (arg);
	case 18: return new Object<ngbla::Vec<18,Complex> > (arg);
	case 24: return new Object<ngbla::Vec<24,Complex> > (arg);
	}
    }

  stringstream err;
  err << "illegal CreateVecObject, dim = " << dim << endl;
  throw Exception (err.str());
}
*/



template <template <class T> class Object, class Base, class ARG, class ARG2>
Base * CreateVecObject (int dim, bool iscomplex, ARG & arg, ARG2 & arg2)
{
  if (!iscomplex)
    {
      switch (dim)
	{
	case 1: return new Object<double> (arg,arg2);
	case 2: return new Object<ngbla::Vec<2> > (arg,arg2);
	case 3: return new Object<ngbla::Vec<3> > (arg,arg2);
	case 4: return new Object<ngbla::Vec<4> > (arg,arg2);
	case 5: return new Object<ngbla::Vec<5> > (arg,arg2);
	case 6: return new Object<ngbla::Vec<6> > (arg,arg2);
	case 7: return new Object<ngbla::Vec<7> > (arg,arg2);
	case 8: return new Object<ngbla::Vec<8> > (arg,arg2);
	case 9: return new Object<ngbla::Vec<9> > (arg,arg2);
	case 10: return new Object<ngbla::Vec<10> > (arg,arg2);
	case 11: return new Object<ngbla::Vec<11> > (arg,arg2);
	case 12: return new Object<ngbla::Vec<12> > (arg,arg2);
	case 18: return new Object<ngbla::Vec<18> > (arg,arg2);
	case 24: return new Object<ngbla::Vec<24> > (arg,arg2);
	}
    }
  else
    {
      switch (dim)
	{
	case 1: return new Object<Complex> (arg,arg2);
	case 2: return new Object<ngbla::Vec<2,Complex> > (arg,arg2);
	case 3: return new Object<ngbla::Vec<3,Complex> > (arg,arg2);
	case 4: return new Object<ngbla::Vec<4,Complex> > (arg,arg2);
	case 5: return new Object<ngbla::Vec<5,Complex> > (arg,arg2);
	case 6: return new Object<ngbla::Vec<6,Complex> > (arg,arg2);
	case 7: return new Object<ngbla::Vec<7,Complex> > (arg,arg2);
	case 8: return new Object<ngbla::Vec<8,Complex> > (arg,arg2);
	case 9: return new Object<ngbla::Vec<9,Complex> > (arg,arg2);
	case 10: return new Object<ngbla::Vec<10,Complex> > (arg,arg2);
	case 11: return new Object<ngbla::Vec<11,Complex> > (arg,arg2);
	case 12: return new Object<ngbla::Vec<12,Complex> > (arg,arg2);
	case 18: return new Object<ngbla::Vec<18,Complex> > (arg,arg2);
	case 24: return new Object<ngbla::Vec<24,Complex> > (arg,arg2);
	}
    }

  stringstream err;
  err << "illegal CreateVecObject, dim = " << dim << endl;
  throw Exception (err.str());
}



/*

template <template <class T> class Object, class Base, class SCAL, class ARG, class ARG2, class ARG3>
Base * CreateVecObjectS (int dim, ARG & arg, ARG2 & arg2, ARG3 & arg3)
{
  switch (dim)
    {
    case 1: return new Object<SCAL> (arg,arg2,arg3);
    case 2: return new Object<ngbla::Vec<2,SCAL> > (arg,arg2,arg3);
    case 3: return new Object<ngbla::Vec<3,SCAL> > (arg,arg2,arg3);
    case 4: return new Object<ngbla::Vec<4,SCAL> > (arg,arg2,arg3);
    case 5: return new Object<ngbla::Vec<5,SCAL> > (arg,arg2,arg3);
    case 6: return new Object<ngbla::Vec<6,SCAL> > (arg,arg2,arg3);
    case 7: return new Object<ngbla::Vec<7,SCAL> > (arg,arg2,arg3);
    case 8: return new Object<ngbla::Vec<8,SCAL> > (arg,arg2,arg3);
    case 9: return new Object<ngbla::Vec<9,SCAL> > (arg,arg2,arg3);
    case 10: return new Object<ngbla::Vec<10,SCAL> > (arg,arg2,arg3);
    case 11: return new Object<ngbla::Vec<11,SCAL> > (arg,arg2,arg3);
    case 12: return new Object<ngbla::Vec<12,SCAL> > (arg,arg2,arg3);
    case 18: return new Object<ngbla::Vec<18,SCAL> > (arg,arg2,arg3);
    case 24: return new Object<ngbla::Vec<24,SCAL> > (arg,arg2,arg3);
    }

  stringstream err;
  err << "illegal CreateVecObject, dim = " << dim << endl;
  throw Exception (err.str());
}


template <template <class T> class Object, class Base, class ARG, class ARG2, class ARG3>
Base * CreateVecObject (int dim, bool iscomplex, ARG & arg, ARG2 & arg2, ARG3 & arg3)
{
  if (!iscomplex)
    return CreateVecObjectS<Object, Base, double, ARG, ARG2, ARG3> (dim, arg, arg2, arg3);
  else
    return CreateVecObjectS<Object, Base, Complex, ARG, ARG2, ARG3> (dim, arg, arg2, arg3);
}
*/





#if MAX_SYS_DIM <= 1

#define CreateMatObject2S(dest, Object, dim1, dim2, SCAL, arg, arg2) \
switch (dim1) {  \
 case 1: dest = new Object<SCAL>(arg, arg2); break; \
} 

#endif

#if MAX_SYS_DIM == 2

#define CreateMatObject2S(dest, Object, dim1, dim2, SCAL, arg, arg2) \
switch (dim1) {  \
 case 1: dest = new Object<SCAL>(arg, arg2); break; \
 case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break; \
} 

#endif


#if MAX_SYS_DIM == 3

#define CreateMatObject2S(dest, Object, dim1, dim2, SCAL, arg, arg2) \
switch (dim1) {  \
 case 1: dest = new Object<SCAL>(arg, arg2); break; \
 case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break; \
 case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2); break;  \
} 

#endif



#if MAX_SYS_DIM == 4

#define CreateMatObject2S(dest, Object, dim1, dim2, SCAL, arg, arg2) \
switch (dim1) {  \
 case 1: dest = new Object<SCAL>(arg, arg2); break; \
 case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break; \
 case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2); break;  \
 case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2); break;  \
} 

#endif



#if MAX_SYS_DIM >= 5

#define CreateMatObject2S(dest, Object, dim1, dim2, SCAL, arg, arg2) \
switch (dim1) {  \
 case 1: dest = new Object<SCAL>(arg, arg2); break; \
 case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break; \
 case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2); break;  \
 case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2); break;  \
 case 5: dest = new Object<Mat<5,5,SCAL> >(arg, arg2); break;  \
 case 6: dest = new Object<Mat<6,6,SCAL> >(arg, arg2); break;  \
 case 7: dest = new Object<Mat<7,7,SCAL> >(arg, arg2); break;  \
 case 8: dest = new Object<Mat<8,8,SCAL> >(arg, arg2); break;  \
} 

#endif


#define CreateMatObject2(dest, Object, dim1, dim2, iscomplex, arg, arg2) \
if (iscomplex) \
  CreateMatObject2S (dest, Object, dim1, dim2, Complex, arg, arg2) \
else \
 CreateMatObject2S (dest, Object, dim1, dim2, double, arg, arg2);



#if MAX_SYS_DIM <= 1

#define CreateSymMatObject2S(dest, Object, dim1, SCAL, arg, arg2) \
switch (dim1) {  \
 case 1: dest = new Object<SCAL>(arg, arg2); break; \
} 

#endif

#if MAX_SYS_DIM == 2

#define CreateSymMatObject2S(dest, Object, dim1, SCAL, arg, arg2) \
switch (dim1) {  \
 case 1: dest = new Object<SCAL>(arg, arg2); break; \
 case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break; \
} 

#endif

#if MAX_SYS_DIM == 3

#define CreateSymMatObject2S(dest, Object, dim1, SCAL, arg, arg2) \
switch (dim1) {  \
 case 1: dest = new Object<SCAL>(arg, arg2); break; \
 case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break; \
 case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2); break;  \
} 

#endif

#if MAX_SYS_DIM == 4

#define CreateSymMatObject2S(dest, Object, dim1, SCAL, arg, arg2) \
switch (dim1) {  \
 case 1: dest = new Object<SCAL>(arg, arg2); break; \
 case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break; \
 case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2); break;  \
 case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2); break;  \
} 

#endif


#if MAX_SYS_DIM >= 5

#define CreateSymMatObject2S(dest, Object, dim1, SCAL, arg, arg2) \
switch (dim1) {  \
 case 1: dest = new Object<SCAL>(arg, arg2); break; \
 case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break; \
 case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2); break;  \
 case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2); break;  \
 case 5: dest = new Object<Mat<5,5,SCAL> >(arg, arg2); break;  \
 case 6: dest = new Object<Mat<6,6,SCAL> >(arg, arg2); break;  \
 case 7: dest = new Object<Mat<7,7,SCAL> >(arg, arg2); break;  \
 case 8: dest = new Object<Mat<8,8,SCAL> >(arg, arg2); break;  \
} 

#endif



#define CreateSymMatObject2(dest, Object, dim, iscomplex, arg, arg2) \
if (iscomplex) \
  CreateSymMatObject2S (dest, Object, dim, Complex, arg, arg2) \
else \
 CreateSymMatObject2S (dest, Object, dim, double, arg, arg2);





















#if MAX_SYS_DIM <= 1

#define CreateMatObject3S(dest, Object, dim1, dim2, SCAL, arg, arg2, arg3) \
switch (dim1) {  \
 case 1: dest = new Object<SCAL>(arg, arg2, arg3); break; \
} 

#endif

#if MAX_SYS_DIM == 2

#define CreateMatObject3S(dest, Object, dim1, dim2, SCAL, arg, arg2, arg3) \
switch (dim1) {  \
 case 1: dest = new Object<SCAL>(arg, arg2, arg3); break; \
 case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break; \
} 

#endif


#if MAX_SYS_DIM == 3

#define CreateMatObject3S(dest, Object, dim1, dim2, SCAL, arg, arg2, arg3) \
switch (dim1) {  \
 case 1: dest = new Object<SCAL>(arg, arg2, arg3); break; \
 case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break; \
 case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2, arg3); break;  \
} 

#endif



#if MAX_SYS_DIM == 4

#define CreateMatObject3S(dest, Object, dim1, dim2, SCAL, arg, arg2, arg3) \
switch (dim1) {  \
 case 1: dest = new Object<SCAL>(arg, arg2, arg3); break; \
 case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break; \
 case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2, arg3); break;  \
 case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2, arg3); break;  \
} 

#endif



#if MAX_SYS_DIM >= 5

#define CreateMatObject3S(dest, Object, dim1, dim2, SCAL, arg, arg2, arg3) \
switch (dim1) {  \
 case 1: dest = new Object<SCAL>(arg, arg2, arg3); break; \
 case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break; \
 case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2, arg3); break;  \
 case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2, arg3); break;  \
 case 5: dest = new Object<Mat<5,5,SCAL> >(arg, arg2, arg3); break;  \
 case 6: dest = new Object<Mat<6,6,SCAL> >(arg, arg2, arg3); break;  \
 case 7: dest = new Object<Mat<7,7,SCAL> >(arg, arg2, arg3); break;  \
 case 8: dest = new Object<Mat<8,8,SCAL> >(arg, arg2, arg3); break;  \
} 

#endif


#define CreateMatObject3(dest, Object, dim1, dim2, iscomplex, arg, arg2, arg3) \
if (iscomplex) \
  CreateMatObject3S (dest, Object, dim1, dim2, Complex, arg, arg2, arg3) \
else \
 CreateMatObject3S (dest, Object, dim1, dim2, double, arg, arg2, arg3);



#if MAX_SYS_DIM <= 1

#define CreateSymMatObject3S(dest, Object, dim1, SCAL, arg, arg2, arg3) \
switch (dim1) {  \
 case 1: dest = new Object<SCAL>(arg, arg2, arg3); break; \
} 

#endif

#if MAX_SYS_DIM == 2

#define CreateSymMatObject3S(dest, Object, dim1, SCAL, arg, arg2, arg3) \
switch (dim1) {  \
 case 1: dest = new Object<SCAL>(arg, arg2, arg3); break; \
 case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break; \
} 

#endif

#if MAX_SYS_DIM == 3

#define CreateSymMatObject3S(dest, Object, dim1, SCAL, arg, arg2, arg3) \
switch (dim1) {  \
 case 1: dest = new Object<SCAL>(arg, arg2, arg3); break; \
 case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break; \
 case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2, arg3); break;  \
} 

#endif

#if MAX_SYS_DIM == 4

#define CreateSymMatObject3S(dest, Object, dim1, SCAL, arg, arg2, arg3) \
switch (dim1) {  \
 case 1: dest = new Object<SCAL>(arg, arg2, arg3); break; \
 case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break; \
 case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2, arg3); break;  \
 case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2, arg3); break;  \
} 

#endif


#if MAX_SYS_DIM >= 5

#define CreateSymMatObject3S(dest, Object, dim1, SCAL, arg, arg2, arg3) \
switch (dim1) {  \
 case 1: dest = new Object<SCAL>(arg, arg2, arg3); break; \
 case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break; \
 case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2, arg3); break;  \
 case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2, arg3); break;  \
 case 5: dest = new Object<Mat<5,5,SCAL> >(arg, arg2, arg3); break;  \
 case 6: dest = new Object<Mat<6,6,SCAL> >(arg, arg2, arg3); break;  \
 case 7: dest = new Object<Mat<7,7,SCAL> >(arg, arg2, arg3); break;  \
 case 8: dest = new Object<Mat<8,8,SCAL> >(arg, arg2, arg3); break;  \
} 

#endif




#if MAX_CACHEBLOCKS < 2

#define CreateSymMatObject4S(dest, Object, dim1, blocksize, SCAL, arg, arg2, arg3) \
switch (dim1) { \
 case 1: \
   switch (blocksize) { \
   case 1: \
     dest = new Object<SCAL>(arg, arg2, arg3); break; \
   } \
   break; \
}


#else

#if MAX_CACHEBLOCKS < 3

#define CreateSymMatObject4S(dest, Object, dim1, blocksize, SCAL, arg, arg2, arg3) \
switch (dim1) { \
 case 1: \
   switch (blocksize) { \
   case 1: \
     dest = new Object<SCAL>(arg, arg2, arg3); break; \
   case 2: \
     dest = new Object<SCAL, Vec<2,SCAL> >(arg, arg2, arg3); break; \
   } \
   break; \
}


#else
#if MAX_CACHEBLOCKS < 5

#define CreateSymMatObject4S(dest, Object, dim1, blocksize, SCAL, arg, arg2, arg3) \
switch (dim1) { \
 case 1: \
   switch (blocksize) { \
   case 1: \
     dest = new Object<SCAL>(arg, arg2, arg3); break; \
   case 2: \
     dest = new Object<SCAL, Vec<2,SCAL> >(arg, arg2, arg3); break; \
   case 3: \
     dest = new Object<SCAL, Vec<3,SCAL> >(arg, arg2, arg3); break; \
   case 4: \
     dest = new Object<SCAL, Vec<4,SCAL> >(arg, arg2, arg3); break; \
   } \
   break; \
}


#else

#define CreateSymMatObject4S(dest, Object, dim1, blocksize, SCAL, arg, arg2, arg3) \
switch (dim1) { \
 case 1: \
   switch (blocksize) { \
   case 1: \
     dest = new Object<SCAL>(arg, arg2, arg3); break; \
   case 2: \
     dest = new Object<SCAL, Vec<2,SCAL> >(arg, arg2, arg3); break; \
   case 3: \
     dest = new Object<SCAL, Vec<3,SCAL> >(arg, arg2, arg3); break; \
   case 4: \
     dest = new Object<SCAL, Vec<4,SCAL> >(arg, arg2, arg3); break; \
   case 5: \
     dest = new Object<SCAL, Vec<5,SCAL> >(arg, arg2, arg3); break; \
   case 6: \
     dest = new Object<SCAL, Vec<6,SCAL> >(arg, arg2, arg3); break; \
   case 7: \
     dest = new Object<SCAL, Vec<7,SCAL> >(arg, arg2, arg3); break; \
   case 8: \
     dest = new Object<SCAL, Vec<8,SCAL> >(arg, arg2, arg3); break; \
   case 9: \
     dest = new Object<SCAL, Vec<9,SCAL> >(arg, arg2, arg3); break; \
   case 10: \
     dest = new Object<SCAL, Vec<10,SCAL> >(arg, arg2, arg3); break; \
   case 11: \
     dest = new Object<SCAL, Vec<11,SCAL> >(arg, arg2, arg3); break; \
   case 12: \
     dest = new Object<SCAL, Vec<12,SCAL> >(arg, arg2, arg3); break; \
   case 13: \
     dest = new Object<SCAL, Vec<13,SCAL> >(arg, arg2, arg3); break; \
   case 14: \
     dest = new Object<SCAL, Vec<14,SCAL> >(arg, arg2, arg3); break; \
   case 15: \
     dest = new Object<SCAL, Vec<15,SCAL> >(arg, arg2, arg3); break; \
   } \
   break; \
}

#endif
#endif
#endif





#define CreateSymMatObject3(dest, Object, dim, iscomplex, arg, arg2, arg3) \
if (iscomplex) \
  CreateSymMatObject3S (dest, Object, dim, Complex, arg, arg2, arg3) \
else \
 CreateSymMatObject3S (dest, Object, dim, double, arg, arg2, arg3);



#define CreateSymMatObject4(dest, Object, dim, blocksize, iscomplex, arg, arg2, arg3) \
if (iscomplex)							\
  CreateSymMatObject4S (dest, Object, dim, blocksize, Complex, arg, arg2, arg3) \
else \
  CreateSymMatObject4S (dest, Object, dim, blocksize, double, arg, arg2, arg3);













/*

template <template <class T> class Object, class Base, class SCAL, class ARG, class ARG2>
Base * S_CreateMatObject (int dim1, int dim2, ARG & arg, ARG2 & arg2)
{
  switch (dim1)
    {
    case 1:
      return new Object<SCAL> (arg,arg2);
    case 2:
      return new Object<ngbla::Mat<2,2,SCAL> > (arg,arg2);
    case 3:
      return new Object<ngbla::Mat<3,3,SCAL> > (arg,arg2);
    case 4:
      return new Object<ngbla::Mat<4,4,SCAL> > (arg,arg2);
    }

  stringstream err;
  err << "illegal CreateVecObject, dim = " << dim << endl;
  throw Exception (err.str());
}


template <template <class T> class Object, class Base, class ARG, class ARG2>
Base * CreateMatObject (int dim1, int dim2, bool iscomplex, ARG & arg, ARG2 & arg2)
{
  if (!iscomplex)
    S_CreateMatObject<Object, Base, double, ARG,ARG2> (dim1, dim2, arg, arg2);
  else
    S_CreateMatObject<Object, Base, Complex, ARG,ARG2> (dim1, dim2, arg, arg2);
}





template <template <class T> class Object, class Base, class SCAL, class ARG, class ARG2>
Base * S_CreateMatObjectSym (int dim, ARG & arg, ARG2 & arg2)
{
  switch (dim)
    {
    case 1:
      return new Object<SCAL> (arg,arg2);
    case 2:
      return new Object<ngbla::Mat<2,2,SCAL> > (arg,arg2);
    case 3:
      return new Object<ngbla::Mat<3,3,SCAL> > (arg,arg2);
    case 4:
      return new Object<ngbla::Mat<4,4,SCAL> > (arg,arg2);
    }

  stringstream err;
  err << "illegal CreateVecObject, dim = " << dim << endl;
  throw Exception (err.str());
}


template <template <class T> class Object, class Base, class ARG, class ARG2>
Base * CreateMatObjectSym (int dim, bool iscomplex, ARG & arg, ARG2 & arg2)
{
  if (!iscomplex)
    S_CreateMatObjectSym<Object, Base, double, ARG,ARG2> (dim, arg, arg2);
  else
    S_CreateMatObjectSym<Object, Base, Complex, ARG,ARG2> (dim, arg, arg2);
}
*/







#endif
