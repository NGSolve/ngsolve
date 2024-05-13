#ifndef FILE_NGS_UTILS
#define FILE_NGS_UTILS

/**************************************************************************/
/* File:   ngs_utils.hpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   Oct. 14                                                        */
/**************************************************************************/


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

  // template <int N>
  // using IC = integral_constant<int,N>;

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
