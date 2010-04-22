#ifndef FILE_DYNAMICMEM
#define FILE_DYNAMICMEM

/**************************************************************************/
/* File:   dynamicmem.hpp                                                 */
/* Author: Joachim Schoeberl                                              */
/* Date:   12. Feb. 2003                                                  */
/**************************************************************************/


namespace netgen
{

  class BaseDynamicMem;

  class BaseDynamicMem
  {
  private:
    static BaseDynamicMem *first;
    static BaseDynamicMem *last;

    BaseDynamicMem *prev, *next;
    size_t size;
    char * ptr;
    char * name;

  protected:
    NGS_DLL_HEADER BaseDynamicMem ();
    NGS_DLL_HEADER ~BaseDynamicMem ();
    NGS_DLL_HEADER void Alloc (size_t s);
    NGS_DLL_HEADER void ReAlloc (size_t s);
    NGS_DLL_HEADER void Free ();
    char * Ptr() const { return ptr; }
//    const char * Ptr() const { return ptr; }
    NGS_DLL_HEADER void Swap (BaseDynamicMem & m2);
  public:
    NGS_DLL_HEADER void SetName (const char * aname);
    NGS_DLL_HEADER static void Print ();
    NGS_DLL_HEADER static void GetUsed (int nr, char * ch);
  };


  template <typename T>
  class DynamicMem : public BaseDynamicMem
  {
  public:
    DynamicMem ()
      : BaseDynamicMem () 
    {
      ;
    }
    DynamicMem (size_t s)
      : BaseDynamicMem () 
    {
      Alloc (s);
    }
    void Alloc (size_t s)
    {
      BaseDynamicMem::Alloc (sizeof(T) * s);
    }
    void ReAlloc (size_t s)
    {
      BaseDynamicMem::ReAlloc (sizeof(T) * s);
    }
    void Free ()
    {
      BaseDynamicMem::Free ();
    }

/*
	const T * Ptr() const
    {
      return reinterpret_cast<const T*> (BaseDynamicMem::Ptr());
    }
*/
    T * Ptr() const
    {
      return reinterpret_cast<T*> (BaseDynamicMem::Ptr());
    }

	/*
    operator const T* () const
    {
      return reinterpret_cast<const T*> (BaseDynamicMem::Ptr());
    }
*/
    operator T* () const
    {
      return reinterpret_cast<T*> (BaseDynamicMem::Ptr());
    }

    void Swap (DynamicMem<T> & m2)
    {
      BaseDynamicMem::Swap (m2);
    }
  protected:
    DynamicMem (const DynamicMem & m);
    DynamicMem & operator= (const DynamicMem & m);
  };
}



#endif
