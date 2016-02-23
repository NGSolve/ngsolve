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

  template <int N>
  using IC = integral_constant<int,N>;

  /*
  template <int I1, int I2>
  constexpr INLINE IC<I1+I2> operator+ (IC<I1> i1, IC<I2> i2) { return IC<I1+I2>(); }
  */
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



 /*
  INLINE void MyAtomicAdd(atomic<double> & sum, double val) {
    auto current = sum.load();
    while (!sum.compare_exchange_weak(current, current + val))
      ;
  }
  */
  INLINE atomic<double> & operator+= (atomic<double> & sum, double val) {
    auto current = sum.load();
    while (!sum.compare_exchange_weak(current, current + val))
      ;
    return sum;
  }

  INLINE atomic<double> & operator-= (atomic<double> & sum, double val) {
    auto current = sum.load();
    while (!sum.compare_exchange_weak(current, current - val))
      ;
    return sum;
  }

  template<typename T>
  INLINE atomic<T> & AsAtomic (T & d)
  {
    return reinterpret_cast<atomic<T>&> (d);
  }
  
  INLINE void MyAtomicAdd (double & sum, double val)
  {
    AsAtomic(sum) += val;
    /*
    auto & asum = reinterpret_cast<atomic<double>&>(sum);
    auto current = asum.load();
    while (!asum.compare_exchange_weak(current, current + val))
      ;
    */
  }
  
}


#endif
