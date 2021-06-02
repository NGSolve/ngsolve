#ifndef FILE_NGS_UTILS
#define FILE_NGS_UTILS

/**************************************************************************/
/* File:   ngs_utils.hpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   Oct. 14                                                        */
/**************************************************************************/

#ifdef WIN32
#include <windows.h>
#else // WIN32
#include <dlfcn.h>
#endif //WIN32


namespace ngstd
{

  
  template <typename T>
  INLINE T RemoveConst (const T & x)
  {
    return x;
  }
  

  // ////////////////    integral constants

  template <int N>
  using IC = integral_constant<int,N>;

  // Class to handle/load shared libraries
  class SharedLibrary
  {
    string lib_name;

#ifdef WIN32
    HINSTANCE lib;
#else // WIN32
    void *lib;
#endif // WIN32

  public:
    SharedLibrary() : lib(nullptr) {}
    SharedLibrary(string lib_name_) : lib(nullptr)
    {
      Load(lib_name_);
    }

    SharedLibrary(const SharedLibrary &) = delete;
    SharedLibrary & operator =(const SharedLibrary &) = delete;

    ~SharedLibrary()
    {
      Unload();
    }

    template <typename TFunc>
    TFunc GetFunction( string func_name )
    {
      return reinterpret_cast<TFunc>(GetRawFunction(func_name));
    }

    void Load( string alib_name )
    {
      Unload();
      lib_name = alib_name;
#ifdef WIN32
      lib = LoadLibrary(lib_name.c_str());
      if (!lib) throw std::runtime_error(string("Could not load library ") + lib_name);
#else // WIN32
      lib = dlopen(lib_name.c_str(), RTLD_NOW);
      if(lib == nullptr) throw std::runtime_error(dlerror());
#endif // WIN32
    }

    void Unload() {
      if(lib)
      {
#ifdef WIN32
        FreeLibrary(lib);
#else // WIN32
        int rc = dlclose(lib);
        if(rc != 0) cerr << "Failed to close library " << lib_name << endl;
#endif // WIN32
      }
    }

    void* GetRawFunction( string func_name )
    {
#ifdef WIN32
      void* func = GetProcAddress(lib, func_name.c_str());
      if(func == nullptr)
        throw std::runtime_error(string("Could not find function ") + func_name + " in library " + lib_name);
#else // WIN32
      void* func = dlsym(lib, func_name.c_str());
      if(func == nullptr)
          throw std::runtime_error(dlerror());
#endif // WIN32

      return func;
    }
  };
}
namespace std
{
  template <int I1, int I2>
  constexpr INLINE integral_constant<int,I1+I2> 
  operator+ (integral_constant<int,I1> /* i1 */,
             integral_constant<int,I2> /* i2 */)
  {
    return integral_constant<int,I1+I2>();
  }


  // may be used as an index e.g. for FlatVector
  template <int N>
  struct is_integral<ngstd::IC<N>>
  {
    enum { value = 1 };
  };



  

#define aligned_alloca(size,align)  (( (size_t)alloca(size+align-1)+align-1) & -align)
#ifdef VLA
#define STACK_ARRAY(TYPE,VAR,SIZE) TYPE VAR[SIZE]
#else
#define STACK_ARRAY(TYPE,VAR,SIZE) TYPE * VAR = (TYPE*)aligned_alloca((SIZE)*sizeof(TYPE), alignof(TYPE))
#endif




  // technique from 
  // http://stackoverflow.com/questions/14939190/boost-shared-from-this-and-multiple-inheritance

  template<typename T>
  struct enable_shared_from_this_virtual;
  
  class enable_shared_from_this_virtual_base : public std::enable_shared_from_this<enable_shared_from_this_virtual_base>
  {
    typedef std::enable_shared_from_this<enable_shared_from_this_virtual_base> base_type;
    template<typename T>
    friend struct enable_shared_from_this_virtual;

    std::shared_ptr<enable_shared_from_this_virtual_base> shared_from_this()
    {
      return base_type::shared_from_this();
    }
    std::shared_ptr<enable_shared_from_this_virtual_base const> shared_from_this() const
    {
      return base_type::shared_from_this();
    }
  public:
    virtual ~enable_shared_from_this_virtual_base() { }
  };
  
  template<typename T>
  struct enable_shared_from_this_virtual: virtual enable_shared_from_this_virtual_base
  {
    typedef enable_shared_from_this_virtual_base base_type;
    
  public:
    std::shared_ptr<T> shared_from_this()
    {
      std::shared_ptr<T> result(base_type::shared_from_this(), static_cast<T*>(this));
      return result;
    }
    
    std::shared_ptr<T const> shared_from_this() const
    {
      std::shared_ptr<T const> result(base_type::shared_from_this(), static_cast<T const*>(this));
      return result;
    }
  };



  template <typename T>
  struct has_call_operator
  {
      template <typename C> static std::true_type check( decltype( sizeof(&C::operator() )) ) { return std::true_type(); }
      template <typename> static std::false_type check(...) { return std::false_type(); }
      typedef decltype( check<T>(sizeof(char)) ) type;
      static constexpr type value = type();
  };

  
}


#endif
