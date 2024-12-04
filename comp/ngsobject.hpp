#ifndef FILE_NGSOBJECT
#define FILE_NGSOBJECT

/*********************************************************************/
/* File:   ngsobject.hpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   2. Aug. 2000                                              */
/*********************************************************************/


#include "meshaccess.hpp"


namespace ngcomp
{
  class DocInfo
  {
  public:
    string short_docu;
    string long_docu;
    std::vector<tuple<string,string>> arguments;

    string & Arg(string name)
    {
      arguments.push_back ( tuple<string,string> (name, "none") );
      return get<1> (arguments.back());
    }

    NGS_DLL_HEADER string GetPythonDocString() const;
  };


  
  /** 
      NGSolve base class
  */
  class NGS_DLL_HEADER NGS_Object : public enable_shared_from_this_virtual<NGS_Object>
  {
  protected:
    /// name of the instance
    string name;

    /// the valid flags for this class
    Flags flaglist;

    // the flags of this class 
    Flags flags;

    /// access to the mesh
    shared_ptr<MeshAccess> ma;

    /// profiling
    Timer<> timer;

    /// keep the object
    bool skipCleanUp; 

    /// timestamp of ngs-objects
    static size_t global_timestamp;
    size_t timestamp = 0;
  public:

    /// 
    NGS_Object (shared_ptr<MeshAccess> ama, Flags aflags, const string & aname = "noname", 
                bool checkflags = false)
      : name(aname), flags(aflags), ma(ama), timer(aname), skipCleanUp(0)
    { 
      ;
    }
  
    ///
    virtual ~NGS_Object ();

    ///
    void SetName (const string & aname)
    { 
      name = aname; 
      NgProfiler::SetName (timer, name); 
    }

    ///
    const string & GetName () const
    {
      return name; 
    }

    const Flags & GetFlags() const { return flags; }

    auto GetTimeStamp() const { return timestamp; } 
    static auto GetNextTimeStamp() { return ++global_timestamp; }
    ///
    const shared_ptr<MeshAccess> & GetMeshAccess() const
    { 
      return ma; 
    }

    virtual string GetClassName () const;
    virtual void PrintReport (ostream & ost) const;
    virtual Array<MemoryUsage> GetMemoryUsage () const;
    
    Timer<> & GetTimer () { return timer; }
    const Timer<> & GetTimer () const { return timer; }

    bool SkipCleanUp () { return skipCleanUp; }
    void SetSkipCleanUp (bool skip = true) { skipCleanUp=skip; }

  protected:
    void DefineStringFlag (const char* s); // , const char* val="");
    void DefineNumFlag (const char* s); // , double val=0);
    void DefineDefineFlag (const char* s);
    void DefineStringListFlag (const char* s);
    void DefineNumListFlag (const char* s);
    int CheckFlags (const Flags& flags); 

    MemoryTracer mem_tracer = {name};

  public:
    const MemoryTracer & GetMemoryTracer() const
    {
      return mem_tracer;
    }

  };


  inline ostream & operator<< (ostream & ost, const NGS_Object & obj)
  {
    obj.PrintReport (ost);
    return ost;
  }
  


  
  
  template <template <class T> class Object, class Base, int ACTDIM, typename ... ARG>
  class TCreateVecObjectS {
  public:
    static Base * Create (int dim, bool iscomplex, ARG & ... arg)
    {
      if (dim == ACTDIM) 
        {
          if (iscomplex)
            return new Object<Vec<ACTDIM,Complex> > (arg...);
          else
            return new Object<Vec<ACTDIM,double> > (arg...);
        }
      else return TCreateVecObjectS<Object, Base, ACTDIM-1, ARG...>::Create(dim, iscomplex, arg...);
    }
  };
  
  template <template <class T> class Object, class Base, typename ... ARG>
  class TCreateVecObjectS<Object, Base, 1, ARG...> {
  public:
    static Base * Create (int dim, bool iscomplex, ARG & ... arg)
    { 
      if (dim == 1) 
        {
          if (iscomplex)
            return new Object<Complex> (arg...);
          else
            return new Object<double> (arg...);
        }
      throw Exception ("illegal CreateVecObject, dim = "
                       + ToString(dim) + '\n');
    }
  };
  
  template <template <class T> class Object, class Base, typename ... ARG>
  Base * CreateVecObject (int dim, bool iscomplex, ARG && ... arg)
  {
    return TCreateVecObjectS<Object, Base, 12, ARG...>::Create (dim, iscomplex, arg ...);
  }
  




  
  template <template <class T> class Object, class Base, int ACTDIM, typename ... ARG>
  class TCreateSharedVecObjectS {
  public:
    static shared_ptr<Base> Create (int dim, bool iscomplex, ARG & ... arg)
    {
      if (dim == ACTDIM) 
        {
          if (iscomplex)
            return make_shared<Object<Vec<ACTDIM,Complex>>> (arg...);
          else
            return make_shared<Object<Vec<ACTDIM,double>>> (arg...);
        }
      else return TCreateSharedVecObjectS<Object, Base, ACTDIM-1, ARG...>::Create(dim, iscomplex, arg...);
    }
  };
  
  template <template <class T> class Object, class Base, typename ... ARG>
  class TCreateSharedVecObjectS<Object, Base, 1, ARG...> {
  public:
    static shared_ptr<Base> Create (int dim, bool iscomplex, ARG & ... arg)
    { 
      if (dim == 1) 
        {
          if (iscomplex)
            return make_shared<Object<Complex>> (arg...);
          else
            return make_shared<Object<double>> (arg...);
        }
      throw Exception ("illegal CreateVecObject, dim = "
                       + ToString(dim) + '\n');
    }
  };
  
  template <template <class T> class Object, class Base, typename ... ARG>
  shared_ptr<Base> CreateSharedVecObject (int dim, bool iscomplex, ARG && ... arg)
  {
    return TCreateSharedVecObjectS<Object, Base, 12, ARG...>::Create (dim, iscomplex, arg ...);
  }
  



  

  /*
    template <template <class T> class Object, class Base, class SCAL, class ARG, int ACTDIM>
    class TCreateVecObjectS {
    Base * Create (int dim, ARG & arg)
    {
    if (dim == ACTDIM) return new Object<Vec<ACTDIM,SCAL> > (arg);
    else return TCreateVecObjectS<Object, Base, SCAL, ARG, ACTDIM-1>::Create(dim, arg);
    }
    };

    template <template <class T> class Object, class Base, class SCAL, class ARG>
    class TCreateVecObjectS<Object, Base, SCAL, ARG, 1> {
    Base * Create (int dim, ARG & arg)
    { 
    if (dim == 1) return new Object<SCAL> (arg);
    throw Exception ("illegal CreateVecObject, dim = "
    + ToString(dim) + '\n');
    }
    };

    template <template <class T> class Object, class Base, class ARG>
    Base * CreateVecObject (int dim, bool iscomplex, ARG & arg)
    {
    if (!iscomplex)
    return TCreateVecObjectS<Object, Base, double, ARG, 12>::Create (dim, arg);
    else
    return TCreateVecObjectS<Object, Base, Complex, ARG, 12>::Create (dim, arg);
    }






    template <template <class T> class Object, class Base, class SCAL, 
    class ARG, class ARG2, class ARG3, int ACTDIM>
    class TCreateVecObject3S { 
    public:
    static Base * Create (int dim, ARG & arg, ARG2 & arg2, ARG3 & arg3)
    {
    if (dim == ACTDIM) return new Object<Vec<ACTDIM,SCAL> > (arg, arg2, arg3);
    else return TCreateVecObject3S<Object, Base, SCAL, ARG, ARG2, ARG3, ACTDIM-1>::Create(dim, arg, arg2, arg3);
    }
    };

    template <template <class T> class Object, class Base, class SCAL, 
    class ARG, class ARG2, class ARG3>
    class TCreateVecObject3S<Object, Base, SCAL, ARG, ARG2, ARG3, 1> {
    public:
    static Base * Create (int dim, ARG & arg, ARG2 & arg2, ARG3 & arg3)
    { 
    if (dim == 1) return new Object<SCAL> (arg, arg2, arg3);
    throw Exception ("illegal CreateVecObject3, dim = "
    + ToString(dim) + '\n');
    }
    };

    template <template <class T> class Object, class Base, class ARG, class ARG2, class ARG3>
    Base * CreateVecObject (int dim, bool iscomplex, ARG & arg, ARG2 & arg2, ARG3 & arg3)
    {
    if (!iscomplex)
    return TCreateVecObject3S<Object, Base, double, ARG, ARG2, ARG3, 8>::Create (dim, arg, arg2, arg3);
    else
    return TCreateVecObject3S<Object, Base, Complex, ARG, ARG2, ARG3, 8>::Create (dim, arg, arg2, arg3);
    }

  */



  /*
    template <template <class T, class TV> class Object, class Base, class SCAL, 
    class ARG, class ARG2, class ARG3, int ACTDIM>
    class TCreateSymMatObject3S { 
    public:
    static Base * Create (int dim, ARG & arg, ARG2 & arg2, ARG3 & arg3)
    {
    if (dim == ACTDIM) return new Object<Mat<ACTDIM,ACTDIM,SCAL>,Vec<ACTDIM,SCAL> > (arg, arg2, arg3);
    else return TCreateSymMatObject3S<Object, Base, SCAL, ARG, ARG2, ARG3, ACTDIM-1>::Create(dim, arg, arg2, arg3);
    }
    };

    template <template <class T, class TV> class Object, class Base, class SCAL, 
    class ARG, class ARG2, class ARG3>
    class TCreateSymMatObject3S<Object, Base, SCAL, ARG, ARG2, ARG3, 1> {
    public:
    static Base * Create (int dim, ARG & arg, ARG2 & arg2, ARG3 & arg3)
    { 
    if (dim == 1) return new Object<SCAL,SCAL> (arg, arg2, arg3);
    throw Exception ("illegal CreateMatObject3, dim = "
    + ToString(dim) + '\n');
    }
    };

    template <template <class T, class TV> class Object, class Base, class ARG, class ARG2, class ARG3>
    Base * CreateSymMatObject (int dim, bool iscomplex, ARG & arg, ARG2 & arg2, ARG3 & arg3)
    {
    if (!iscomplex)
    return TCreateSymMatObject3S<Object, Base, double, ARG, ARG2, ARG3, MAX_SYS_DIM>::Create (dim, arg, arg2, arg3);
    else
    return TCreateSymMatObject3S<Object, Base, Complex, ARG, ARG2, ARG3, MAX_SYS_DIM>::Create (dim, arg, arg2, arg3);
    }
  */



  template <template <typename...T> class Object, class Base, int ACTDIM, typename ... ARG>
  class TCreateSymMatObjectS {
  public:
    static Base * Create (int dim, bool iscomplex, ARG & ... arg)
    {
      if (dim == ACTDIM) 
        {
          if (iscomplex)
            return new Object<Mat<ACTDIM,ACTDIM,Complex> > (arg...);
          else
            return new Object<Mat<ACTDIM,ACTDIM,double> > (arg...);
        }
      else return TCreateSymMatObjectS<Object, Base, ACTDIM-1, ARG...>::Create(dim, iscomplex, arg...);
    }
  };
  
  template <template <typename...T> class Object, class Base, typename ... ARG>
  class TCreateSymMatObjectS<Object, Base, 1, ARG...> {
  public:
    static Base * Create (int dim, bool iscomplex, ARG & ... arg)
    { 
      if (dim == 1) 
        {
          if (iscomplex)
            return new Object<Complex> (arg...);
          else
            return new Object<double> (arg...);
        }
      throw Exception ("illegal CreateSymMatObject, dim = "
                       + ToString(dim) + '\n');
    }
  };
  
  template <template <typename...T> class Object, class Base, typename ... ARG>
  Base * CreateSymMatObject (int dim, bool iscomplex, ARG && ... arg)
  {
    return TCreateSymMatObjectS<Object, Base, MAX_SYS_DIM, ARG...>::Create (dim, iscomplex, arg ...);
  }






  /*


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

  throw Exception ("illegal CreateVecObject, dim = " + ToString(dim) + "\n");
  }
  */





  

#if MAX_SYS_DIM <= 1

#define CreateMatObject2S(dest, Object, dim1, dim2, SCAL, arg, arg2)    \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2); break;                    \
  } 

#endif

#if MAX_SYS_DIM == 2

#define CreateMatObject2S(dest, Object, dim1, dim2, SCAL, arg, arg2)    \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2); break;                    \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break;          \
  } 

#endif


#if MAX_SYS_DIM == 3

#define CreateMatObject2S(dest, Object, dim1, dim2, SCAL, arg, arg2)    \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2); break;                    \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break;          \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2); break;          \
  } 

#endif



#if MAX_SYS_DIM == 4

#define CreateMatObject2S(dest, Object, dim1, dim2, SCAL, arg, arg2)    \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2); break;                    \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break;          \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2); break;          \
  case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2); break;          \
  } 

#endif



#if MAX_SYS_DIM == 5

#define CreateMatObject2S(dest, Object, dim1, dim2, SCAL, arg, arg2)    \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2); break;                    \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break;          \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2); break;          \
  case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2); break;          \
  case 5: dest = new Object<Mat<5,5,SCAL> >(arg, arg2); break;          \
  } 

#endif

#if MAX_SYS_DIM == 6

#define CreateMatObject2S(dest, Object, dim1, dim2, SCAL, arg, arg2)    \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2); break;                    \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break;          \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2); break;          \
  case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2); break;          \
  case 5: dest = new Object<Mat<5,5,SCAL> >(arg, arg2); break;          \
  case 6: dest = new Object<Mat<6,6,SCAL> >(arg, arg2); break;          \
  } 

#endif

#if MAX_SYS_DIM == 7

#define CreateMatObject2S(dest, Object, dim1, dim2, SCAL, arg, arg2)    \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2); break;                    \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break;          \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2); break;          \
  case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2); break;          \
  case 5: dest = new Object<Mat<5,5,SCAL> >(arg, arg2); break;          \
  case 6: dest = new Object<Mat<6,6,SCAL> >(arg, arg2); break;          \
  case 7: dest = new Object<Mat<7,7,SCAL> >(arg, arg2); break;          \
  } 

#endif

#if MAX_SYS_DIM >= 8

#define CreateMatObject2S(dest, Object, dim1, dim2, SCAL, arg, arg2)    \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2); break;                    \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break;          \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2); break;          \
  case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2); break;          \
  case 5: dest = new Object<Mat<5,5,SCAL> >(arg, arg2); break;          \
  case 6: dest = new Object<Mat<6,6,SCAL> >(arg, arg2); break;          \
  case 7: dest = new Object<Mat<7,7,SCAL> >(arg, arg2); break;          \
  case 8: dest = new Object<Mat<8,8,SCAL> >(arg, arg2); break;          \
  } 

#endif


#define CreateMatObject2(dest, Object, dim1, dim2, iscomplex, arg, arg2) \
  if (iscomplex)                                                        \
    CreateMatObject2S (dest, Object, dim1, dim2, Complex, arg, arg2)    \
    else                                                                \
      CreateMatObject2S (dest, Object, dim1, dim2, double, arg, arg2);



#if MAX_SYS_DIM <= 1

#define CreateSymMatObject2S(dest, Object, dim1, SCAL, arg, arg2)       \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2); break;                    \
  } 

#endif

#if MAX_SYS_DIM == 2

#define CreateSymMatObject2S(dest, Object, dim1, SCAL, arg, arg2)       \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2); break;                    \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break;          \
  } 

#endif

#if MAX_SYS_DIM == 3

#define CreateSymMatObject2S(dest, Object, dim1, SCAL, arg, arg2)       \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2); break;                    \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break;          \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2); break;          \
  } 

#endif

#if MAX_SYS_DIM == 4

#define CreateSymMatObject2S(dest, Object, dim1, SCAL, arg, arg2)       \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2); break;                    \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break;          \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2); break;          \
  case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2); break;          \
  } 

#endif


#if MAX_SYS_DIM == 5

#define CreateSymMatObject2S(dest, Object, dim1, SCAL, arg, arg2)       \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2); break;                    \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break;          \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2); break;          \
  case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2); break;          \
  case 5: dest = new Object<Mat<5,5,SCAL> >(arg, arg2); break;          \
  } 

#endif

#if MAX_SYS_DIM == 6

#define CreateSymMatObject2S(dest, Object, dim1, SCAL, arg, arg2)       \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2); break;                    \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break;          \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2); break;          \
  case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2); break;          \
  case 5: dest = new Object<Mat<5,5,SCAL> >(arg, arg2); break;          \
  case 6: dest = new Object<Mat<6,6,SCAL> >(arg, arg2); break;          \
  } 

#endif

#if MAX_SYS_DIM == 7

#define CreateSymMatObject2S(dest, Object, dim1, SCAL, arg, arg2)       \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2); break;                    \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break;          \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2); break;          \
  case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2); break;          \
  case 5: dest = new Object<Mat<5,5,SCAL> >(arg, arg2); break;          \
  case 6: dest = new Object<Mat<6,6,SCAL> >(arg, arg2); break;          \
  case 7: dest = new Object<Mat<7,7,SCAL> >(arg, arg2); break;          \
  } 

#endif

#if MAX_SYS_DIM >= 8

#define CreateSymMatObject2S(dest, Object, dim1, SCAL, arg, arg2)       \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2); break;                    \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2); break;          \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2); break;          \
  case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2); break;          \
  case 5: dest = new Object<Mat<5,5,SCAL> >(arg, arg2); break;          \
  case 6: dest = new Object<Mat<6,6,SCAL> >(arg, arg2); break;          \
  case 7: dest = new Object<Mat<7,7,SCAL> >(arg, arg2); break;          \
  case 8: dest = new Object<Mat<8,8,SCAL> >(arg, arg2); break;          \
  } 

#endif



#define CreateSymMatObject2(dest, Object, dim, iscomplex, arg, arg2)    \
  if (iscomplex)                                                        \
    CreateSymMatObject2S (dest, Object, dim, Complex, arg, arg2)        \
    else                                                                \
      CreateSymMatObject2S (dest, Object, dim, double, arg, arg2);





















#if MAX_SYS_DIM <= 1

#define CreateMatObject3S(dest, Object, dim1, dim2, SCAL, arg, arg2, arg3) \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2, arg3); break;              \
  } 

#endif

#if MAX_SYS_DIM == 2

#define CreateMatObject3S(dest, Object, dim1, dim2, SCAL, arg, arg2, arg3) \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2, arg3); break;              \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break;    \
  } 

#endif


#if MAX_SYS_DIM == 3

#define CreateMatObject3S(dest, Object, dim1, dim2, SCAL, arg, arg2, arg3) \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2, arg3); break;              \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break;    \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2, arg3); break;    \
  } 

#endif



#if MAX_SYS_DIM == 4

#define CreateMatObject3S(dest, Object, dim1, dim2, SCAL, arg, arg2, arg3) \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2, arg3); break;              \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break;    \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2, arg3); break;    \
  case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2, arg3); break;    \
  } 

#endif



#if MAX_SYS_DIM == 5

#define CreateMatObject3S(dest, Object, dim1, dim2, SCAL, arg, arg2, arg3) \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2, arg3); break;              \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break;    \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2, arg3); break;    \
  case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2, arg3); break;    \
  case 5: dest = new Object<Mat<5,5,SCAL> >(arg, arg2, arg3); break;    \
  } 

#endif

#if MAX_SYS_DIM == 6

#define CreateMatObject3S(dest, Object, dim1, dim2, SCAL, arg, arg2, arg3) \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2, arg3); break;              \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break;    \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2, arg3); break;    \
  case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2, arg3); break;    \
  case 5: dest = new Object<Mat<5,5,SCAL> >(arg, arg2, arg3); break;    \
  case 6: dest = new Object<Mat<6,6,SCAL> >(arg, arg2, arg3); break;    \
  } 

#endif

#if MAX_SYS_DIM == 7

#define CreateMatObject3S(dest, Object, dim1, dim2, SCAL, arg, arg2, arg3) \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2, arg3); break;              \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break;    \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2, arg3); break;    \
  case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2, arg3); break;    \
  case 5: dest = new Object<Mat<5,5,SCAL> >(arg, arg2, arg3); break;    \
  case 6: dest = new Object<Mat<6,6,SCAL> >(arg, arg2, arg3); break;    \
  case 7: dest = new Object<Mat<7,7,SCAL> >(arg, arg2, arg3); break;    \
  } 

#endif

#if MAX_SYS_DIM >= 8

#define CreateMatObject3S(dest, Object, dim1, dim2, SCAL, arg, arg2, arg3) \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2, arg3); break;              \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break;    \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2, arg3); break;    \
  case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2, arg3); break;    \
  case 5: dest = new Object<Mat<5,5,SCAL> >(arg, arg2, arg3); break;    \
  case 6: dest = new Object<Mat<6,6,SCAL> >(arg, arg2, arg3); break;    \
  case 7: dest = new Object<Mat<7,7,SCAL> >(arg, arg2, arg3); break;    \
  case 8: dest = new Object<Mat<8,8,SCAL> >(arg, arg2, arg3); break;    \
  } 

#endif


#define CreateMatObject3(dest, Object, dim1, dim2, iscomplex, arg, arg2, arg3) \
  if (iscomplex)                                                        \
    CreateMatObject3S (dest, Object, dim1, dim2, Complex, arg, arg2, arg3) \
    else                                                                \
      CreateMatObject3S (dest, Object, dim1, dim2, double, arg, arg2, arg3);



#if MAX_SYS_DIM <= 1

#define CreateSymMatObject3S(dest, Object, dim1, SCAL, arg, arg2, arg3) \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2, arg3); break;              \
  } 

#endif

#if MAX_SYS_DIM == 2

#define CreateSymMatObject3S(dest, Object, dim1, SCAL, arg, arg2, arg3) \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2, arg3); break;              \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break;    \
  } 

#endif

#if MAX_SYS_DIM == 3

#define CreateSymMatObject3S(dest, Object, dim1, SCAL, arg, arg2, arg3) \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2, arg3); break;              \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break;    \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2, arg3); break;    \
  } 

#endif

#if MAX_SYS_DIM == 4

#define CreateSymMatObject3S(dest, Object, dim1, SCAL, arg, arg2, arg3) \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2, arg3); break;              \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break;    \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2, arg3); break;    \
  case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2, arg3); break;    \
  } 

#endif

#if MAX_SYS_DIM == 5

#define CreateSymMatObject3S(dest, Object, dim1, SCAL, arg, arg2, arg3) \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2, arg3); break;              \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break;    \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2, arg3); break;    \
  case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2, arg3); break;    \
  case 5: dest = new Object<Mat<5,5,SCAL> >(arg, arg2, arg3); break;    \
  } 

#endif

#if MAX_SYS_DIM == 6

#define CreateSymMatObject3S(dest, Object, dim1, SCAL, arg, arg2, arg3) \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2, arg3); break;              \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break;    \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2, arg3); break;    \
  case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2, arg3); break;    \
  case 5: dest = new Object<Mat<5,5,SCAL> >(arg, arg2, arg3); break;    \
  case 6: dest = new Object<Mat<6,6,SCAL> >(arg, arg2, arg3); break;    \
  } 

#endif

#if MAX_SYS_DIM == 7

#define CreateSymMatObject3S(dest, Object, dim1, SCAL, arg, arg2, arg3) \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2, arg3); break;              \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break;    \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2, arg3); break;    \
  case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2, arg3); break;    \
  case 5: dest = new Object<Mat<5,5,SCAL> >(arg, arg2, arg3); break;    \
  case 6: dest = new Object<Mat<6,6,SCAL> >(arg, arg2, arg3); break;    \
  case 7: dest = new Object<Mat<7,7,SCAL> >(arg, arg2, arg3); break;    \
  } 

#endif

#if MAX_SYS_DIM >= 8

#define CreateSymMatObject3S(dest, Object, dim1, SCAL, arg, arg2, arg3) \
  switch (dim1) {                                                       \
  case 1: dest = new Object<SCAL>(arg, arg2, arg3); break;              \
  case 2: dest = new Object<Mat<2,2,SCAL> >(arg, arg2, arg3); break;    \
  case 3: dest = new Object<Mat<3,3,SCAL> >(arg, arg2, arg3); break;    \
  case 4: dest = new Object<Mat<4,4,SCAL> >(arg, arg2, arg3); break;    \
  case 5: dest = new Object<Mat<5,5,SCAL> >(arg, arg2, arg3); break;    \
  case 6: dest = new Object<Mat<6,6,SCAL> >(arg, arg2, arg3); break;    \
  case 7: dest = new Object<Mat<7,7,SCAL> >(arg, arg2, arg3); break;    \
  case 8: dest = new Object<Mat<8,8,SCAL> >(arg, arg2, arg3); break;    \
  } 

#endif

  










#if MAX_CACHEBLOCKS < 2

#define CreateSymMatObject4S(dest, Object, dim1, blocksize, SCAL, arg, arg2, arg3) \
  switch (dim1) {                                                       \
  case 1:                                                               \
    switch (blocksize) {                                                \
    case 1:                                                             \
      dest = new Object<SCAL>(arg, arg2, arg3); break;                  \
    }                                                                   \
    break;                                                              \
  }


#else

#if MAX_CACHEBLOCKS < 3

#define CreateSymMatObject4S(dest, Object, dim1, blocksize, SCAL, arg, arg2, arg3) \
  switch (dim1) {                                                       \
  case 1:                                                               \
    switch (blocksize) {                                                \
    case 1:                                                             \
      dest = new Object<SCAL>(arg, arg2, arg3); break;                  \
    case 2:                                                             \
      dest = new Object<SCAL, Vec<2,SCAL> >(arg, arg2, arg3); break;    \
    }                                                                   \
    break;                                                              \
  }


#else
#if MAX_CACHEBLOCKS < 5

#define CreateSymMatObject4S(dest, Object, dim1, blocksize, SCAL, arg, arg2, arg3) \
  switch (dim1) {                                                       \
  case 1:                                                               \
    switch (blocksize) {                                                \
    case 1:                                                             \
      dest = new Object<SCAL>(arg, arg2, arg3); break;                  \
    case 2:                                                             \
      dest = new Object<SCAL, Vec<2,SCAL> >(arg, arg2, arg3); break;    \
    case 3:                                                             \
      dest = new Object<SCAL, Vec<3,SCAL> >(arg, arg2, arg3); break;    \
    case 4:                                                             \
      dest = new Object<SCAL, Vec<4,SCAL> >(arg, arg2, arg3); break;    \
    }                                                                   \
    break;                                                              \
  }


#else

#define CreateSymMatObject4S(dest, Object, dim1, blocksize, SCAL, arg, arg2, arg3) \
  switch (dim1) {                                                       \
  case 1:                                                               \
    switch (blocksize) {                                                \
    case 1:                                                             \
      dest = new Object<SCAL>(arg, arg2, arg3); break;                  \
    case 2:                                                             \
      dest = new Object<SCAL, Vec<2,SCAL> >(arg, arg2, arg3); break;    \
    case 3:                                                             \
      dest = new Object<SCAL, Vec<3,SCAL> >(arg, arg2, arg3); break;    \
    case 4:                                                             \
      dest = new Object<SCAL, Vec<4,SCAL> >(arg, arg2, arg3); break;    \
    case 5:                                                             \
      dest = new Object<SCAL, Vec<5,SCAL> >(arg, arg2, arg3); break;    \
    case 6:                                                             \
      dest = new Object<SCAL, Vec<6,SCAL> >(arg, arg2, arg3); break;    \
    case 7:                                                             \
      dest = new Object<SCAL, Vec<7,SCAL> >(arg, arg2, arg3); break;    \
    case 8:                                                             \
      dest = new Object<SCAL, Vec<8,SCAL> >(arg, arg2, arg3); break;    \
    case 9:                                                             \
      dest = new Object<SCAL, Vec<9,SCAL> >(arg, arg2, arg3); break;    \
    case 10:                                                            \
      dest = new Object<SCAL, Vec<10,SCAL> >(arg, arg2, arg3); break;   \
    case 11:                                                            \
      dest = new Object<SCAL, Vec<11,SCAL> >(arg, arg2, arg3); break;   \
    case 12:                                                            \
      dest = new Object<SCAL, Vec<12,SCAL> >(arg, arg2, arg3); break;   \
    case 13:                                                            \
      dest = new Object<SCAL, Vec<13,SCAL> >(arg, arg2, arg3); break;   \
    case 14:                                                            \
      dest = new Object<SCAL, Vec<14,SCAL> >(arg, arg2, arg3); break;   \
    case 15:                                                            \
      dest = new Object<SCAL, Vec<15,SCAL> >(arg, arg2, arg3); break;   \
    }                                                                   \
    break;                                                              \
  }

#endif
#endif
#endif





#define CreateSymMatObject3(dest, Object, dim, iscomplex, arg, arg2, arg3) \
  if (iscomplex)                                                        \
    CreateSymMatObject3S (dest, Object, dim, Complex, arg, arg2, arg3)  \
    else                                                                \
      CreateSymMatObject3S (dest, Object, dim, double, arg, arg2, arg3);



#define CreateSymMatObject4(dest, Object, dim, blocksize, iscomplex, arg, arg2, arg3) \
  if (iscomplex)							\
    CreateSymMatObject4S (dest, Object, dim, blocksize, Complex, arg, arg2, arg3) \
    else                                                                \
      CreateSymMatObject4S (dest, Object, dim, blocksize, double, arg, arg2, arg3);





}



#endif
