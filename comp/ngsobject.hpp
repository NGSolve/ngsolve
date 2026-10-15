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
  };


  inline ostream & operator<< (ostream & ost, const NGS_Object & obj)
  {
    obj.PrintReport (ost);
    return ost;
  }



  /*
    Object<TM> for the matrix entry type of a system with dim components:
    TM = SCAL for dim 1, else Mat<dim,dim,SCAL>.
    nullptr if dim exceeds MAX_SYS_DIM.
  */
  template <typename Base, template <typename...> class Object, typename SCAL, typename ... ARGS>
  shared_ptr<Base> CreateSysMatObject (int dim, const ARGS & ... args)
  {
    shared_ptr<Base> res;
    if (dim < 1 || dim > MAX_SYS_DIM) return res;
    Switch<MAX_SYS_DIM> (dim-1, [&] (auto DM)
    {
      constexpr int D = DM+1;
      using TM = std::conditional_t<D==1, SCAL, Mat<D,D,SCAL>>;
      res = make_shared<Object<TM>> (args...);
    });
    return res;
  }

  template <typename Base, template <typename...> class Object, typename ... ARGS>
  shared_ptr<Base> CreateSysMatObject (int dim, bool iscomplex, const ARGS & ... args)
  {
    if (iscomplex)
      return CreateSysMatObject<Base,Object,Complex> (dim, args...);
    return CreateSysMatObject<Base,Object,double> (dim, args...);
  }

}



#endif
