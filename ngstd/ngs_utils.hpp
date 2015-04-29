#ifndef FILE_NGS_UTILS
#define FILE_NGS_UTILS

/**************************************************************************/
/* File:   ngs_utils.hpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   Oct. 14                                                        */
/**************************************************************************/


namespace ngstd
{

  
  template <typename T>
  INLINE T RemoveConst (const T & x)
  {
    return x;
  }
  

  //////////////////////////////////////////////////////////////////////
  // Lambda to function pointer conversion,  M. Hochsteger
  
  template <typename Function>
  struct function_traits
    : public function_traits<decltype(&Function::operator())> 
  { };
  
  template <typename ClassType, typename ReturnType, typename... Args>
  struct function_traits<ReturnType(ClassType::*)(Args...) const> 
  {
    typedef ReturnType (*pointer)(Args...);
    typedef ReturnType return_type;
  };
  
  template <typename Function>
  typename function_traits<Function>::pointer
  FunctionPointer (const Function& lambda) 
  {
    return static_cast<typename function_traits<Function>::pointer>(lambda);
  }

  template <typename Function>
  typename function_traits<Function>::return_type
  GetReturnValue (Function func) 
  {
    typename function_traits<Function>::return_type *dummy;
    return *dummy;
  }




  



  // ////////////////    integral constants
}



namespace std
{

  template <int I1, int I2>
  constexpr INLINE integral_constant<int,I1+I2> 
  operator+ (integral_constant<int,I1> i1,
             integral_constant<int,I2> i2)
  {
    return integral_constant<int,I1+I2>();
  }


}


#endif
