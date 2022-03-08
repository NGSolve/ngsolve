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

#include <filesystem>
#include <string>

namespace ngstd
{
  using namespace std;
  
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
    filesystem::path lib_name;
    optional<filesystem::path> directory_to_delete = nullopt;

#ifdef WIN32
    HINSTANCE lib = nullptr;
#else // WIN32
    void *lib = nullptr;
#endif // WIN32

  public:
    SharedLibrary() = default;
    SharedLibrary(const filesystem::path & lib_name_, optional<filesystem::path> directory_to_delete_ = nullopt );

    SharedLibrary(const SharedLibrary &) = delete;
    SharedLibrary & operator =(const SharedLibrary &) = delete;

    ~SharedLibrary();

    template <typename TFunc>
    TFunc GetFunction( string func_name )
    {
      return reinterpret_cast<TFunc>(GetRawFunction(func_name));
    }

    void Load( const filesystem::path & lib_name_ );
    void Unload();
    void* GetRawFunction( string func_name );
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
