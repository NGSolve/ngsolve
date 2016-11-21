#ifndef FILE_NGS_Array
#define FILE_NGS_Array

/**************************************************************************/
/* File:   array.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

#ifdef DEBUG
#define CHECK_RANGE
#endif

#include "tuple.hpp"

namespace ngstd
{

  /**
     Exception thrown by array range check.
     Only thrown when compiled with RANGE_CHECK
  */
  class ArrayRangeException : public Exception
  {
  public:
    ArrayRangeException () : Exception("ArrayRangeException\n") { ; }
  };



  /*
    Some class which can be treated as array
   */
  template <typename T> // , typename TA = T>
  class BaseArrayObject
  {
  public:
    INLINE BaseArrayObject() { ; }

    INLINE const T & Spec() const { return static_cast<const T&> (*this); }
    INLINE size_t Size() const { return Spec().Size(); }
    template <typename T2>
    INLINE bool Contains(const T2 & el) const
    {
      for (size_t i = 0; i < Size(); i++)
        if (Spec()[i] == el)
          return true;
      return false;
    }
    template <typename T2>
    INLINE int Pos(const T2 & el) const
    {
      for (int i = 0; i < Size(); i++)
        if (Spec()[i] == el)
          return i;
      return -1;
    }
    // INLINE auto operator[] (int i) -> decltype (Spec()[i]) { return Spec()[i]; }
    // INLINE auto operator[] (int i) const -> decltype (T::operator[](i)) { return Spec()[i]; }
  };


  template <typename T>
  inline ostream & operator<< (ostream & ost, const BaseArrayObject<T> & array)
  {
    for (auto i : array.Range())
      ost << i << ":" << array.Spec()[i] << endl;
    return ost;
  }



  template <typename AO>
  class AOWrapperIterator
  {
    const AO & ao;
    size_t ind;
  public:
    INLINE AOWrapperIterator (const AO &  aao, size_t ai) 
      : ao(aao), ind(ai) { ; }
    INLINE AOWrapperIterator operator++ (int) 
    { return AOWrapperIterator(ao, ind++); }
    INLINE AOWrapperIterator operator++ ()
    { return AOWrapperIterator(ao, ++ind); }
    INLINE auto operator*() const -> decltype(ao[ind]) { return ao[ind]; }
    INLINE auto operator*() -> decltype(ao[ind]) { return ao[ind]; }
    INLINE bool operator != (AOWrapperIterator d2) { return ind != d2.ind; }
  };


  
  template <typename T>
  class AOWrapper : public BaseArrayObject<AOWrapper<T>>
  {
    T ar;
  public:
    INLINE AOWrapper (T aar) : ar(aar) { ; }
    INLINE operator T () const { return ar; }
    INLINE size_t Size() const { return ar.Size(); }
    INLINE auto operator[] (size_t i) { return ar[i]; }
    INLINE auto operator[] (size_t i) const { return ar[i]; }
    INLINE AOWrapperIterator<AOWrapper> begin () { return AOWrapperIterator<AOWrapper> (*this, 0); }
    INLINE AOWrapperIterator<AOWrapper> end () { return AOWrapperIterator<AOWrapper> (*this, Size()); }
  };

  template <typename T>
  INLINE AOWrapper<const T&> ArrayObject (const T & ar)
  {
    return AOWrapper<const T&> (ar);
  }

  template <typename T>
  INLINE AOWrapper<T> ArrayObject (T && ar)
  {
    return AOWrapper<T> (ar);
  }







  /**
     nothing more but a new type for a C array.
     return value for Addr - operator of array 
  */
  template <class T>
  class CArray
  {
  protected:
    /// the data
    T * data;
  public:

    /// initialize array 
    INLINE CArray () { data = 0; }

    /// provide size and memory
    INLINE CArray (T * adata) 
      : data(adata) { ; }

    /// Access array
    INLINE T & operator[] (size_t i) const
    {
      return data[i]; 
    }

    INLINE operator T* () const { return data; }
  };

  template <class T> class FlatArray;


  template <typename TELEM, typename TSIZE>
  class ArrayIterator
  {
    FlatArray<TELEM> ar;
    TSIZE ind;
  public:
    INLINE ArrayIterator (FlatArray<TELEM> aar, TSIZE ai) 
      : ar(aar), ind(ai) { ; }
    INLINE ArrayIterator operator++ (int) 
    { return ArrayIterator(ar, ind++); }
    INLINE ArrayIterator operator++ ()
    { return ArrayIterator(ar, ++ind); }
    INLINE TELEM operator*() const { return ar[ind]; }
    INLINE TELEM & operator*() { return ar[ind]; }
    INLINE bool operator != (ArrayIterator d2) { return ind != d2.ind; }
  };



  template <typename TSIZE>
  class ArrayRangeIterator
  {
    TSIZE ind;
  public:
    INLINE ArrayRangeIterator (TSIZE ai) : ind(ai) { ; }
    INLINE ArrayRangeIterator operator++ (int) { return ind++; }
    INLINE ArrayRangeIterator operator++ () { return ++ind; }
    INLINE TSIZE operator*() const { return ind; }
    INLINE TSIZE Index() { return ind; }
    INLINE operator TSIZE () const { return ind; }
    INLINE bool operator != (ArrayRangeIterator d2) { return ind != d2.ind; }
  };

  /// a range of intergers
  template <typename T>
  class T_Range : public BaseArrayObject <T_Range<T>>
  {
    T first, next;
  public: 
    INLINE T_Range () { ; }
    INLINE T_Range (T n) : first(0), next(n) {;}
    INLINE T_Range (T f, T n) : first(f), next(n) {;}
    template <typename T2>
      INLINE T_Range(T_Range<T2> r2) : first(r2.First()), next(r2.Next()) { ; }
    INLINE T First() const { return first; }
    INLINE T Next() const { return next; }
    INLINE T Size() const { return next-first; }
    INLINE T operator[] (T i) const { return first+i; }
    INLINE bool Contains (T i) const { return ((i >= first) && (i < next)); }

    INLINE ArrayRangeIterator<T> begin() const { return first; }
    INLINE ArrayRangeIterator<T> end() const { return next; }

    T_Range Split (size_t nr, int tot) const
    {
      T diff = next-first;
      return T_Range (first + nr * diff / tot,
                      first + (nr+1) * diff / tot);
    }
    // INLINE operator IntRange () const { return IntRange(first,next); }
  };

  using IntRange = T_Range<size_t>;

  template <typename T>
  INLINE T_Range<T> Range (T a, T b)
  {
    return T_Range<T>(a,b);
  }

  template <typename T>
  INLINE T_Range<T> Range_impl (T n, std::true_type)
  {
    return T_Range<T> (0, n);
  }

  template <typename T>
  INLINE auto Range_impl (const T & ao, std::false_type)
    -> T_Range<decltype(ao.Size())>
  {
    return T_Range<decltype(ao.Size())> (0, ao.Size());
  }

  template <typename T>
  auto Range(const T & x) -> decltype(Range_impl(x, std::is_integral<T>())) {
    return Range_impl(x, std::is_integral<T>());
  }


  INLINE IntRange operator+ (const IntRange & range, int shift)
  {
    return IntRange (range.First()+shift, range.Next()+shift);
  }

  INLINE IntRange operator+ (int shift, const IntRange & range)
  {
    return IntRange (range.First()+shift, range.Next()+shift);
  }

  INLINE IntRange operator* (int scale, const IntRange & range)
  {
    return IntRange (scale*range.First(), scale*range.Next());
  }

  INLINE IntRange operator* (const IntRange & range, int scale)
  {
    return IntRange (scale*range.First(), scale*range.Next());
  }

  template <typename TI>
  inline ostream & operator<< (ostream & s, T_Range<TI> ir)
  {
    s << "[" << ir.First() << "," << ir.Next() << ")";
    return s;
  }

  template <typename ... ARGS>
  ostream & operator<< (ostream & ost, Tuple<IntRange, ARGS...> tup)
  {
    ost << tup.Head() << ", " << tup.Tail();
    return ost;
  }



  template <typename T, typename INDEX_ARRAY>
  class IndirectArray : public BaseArrayObject<IndirectArray<T, INDEX_ARRAY> >
  {
    FlatArray<T> ba;
    const INDEX_ARRAY & ia;

  public:
    INLINE IndirectArray (FlatArray<T> aba,
                          const INDEX_ARRAY & aia)
      : ba(aba), ia(aia) { ; }
    
    INLINE size_t Size() const { return ia.Size(); }
    INLINE T operator[] (size_t i) const { return ba[ia.Spec()[i]]; }
    INLINE T & operator[] (size_t i) { return ba[ia.Spec()[i]]; }

    INLINE IndirectArray operator= (const T & val) 
    {
      for (auto i : Range(Size()))
        (*this)[i] = val;
      return IndirectArray (ba, ia);
    }

    template <typename T2>
    INLINE IndirectArray operator= (const BaseArrayObject<T2> & a2) 
    {
      for (auto i : Range(Size()))
	(*this)[i] = a2.Spec()[i];
      return IndirectArray (ba, ia);
    }

  };


  /**
     A simple array container.
     Array represented by size and data-pointer.
     No memory allocation and deallocation, must be provided by user.
     Helper functions for printing. 
     Optional range check by macro CHECK_RANGE
  */
  template <class T>
  class FlatArray : public BaseArrayObject<FlatArray<T> >
  {
  protected:
    /// the size
    size_t size;
    /// the data
    T * __restrict data;
  public:

    /// initialize array 
    INLINE FlatArray () { ; } // size = 0; data = 0; }

    /// copy constructor allows size-type conversion 
    INLINE FlatArray (const FlatArray<T> & a2)  
      : size(a2.Size()), data(&a2[0]) { ; } 

    /// provide size and memory
    INLINE FlatArray (size_t asize, T * adata) 
      : size(asize), data(adata) { ; }
    
    /// memory from local heap
    INLINE FlatArray(size_t asize, Allocator & lh)
      : size(asize), data(new (lh) T[asize])
    { ; }

    INLINE FlatArray(size_t asize, LocalHeap & lh)
      : size(asize), data (lh.Alloc<T> (asize))
    { ; }

    /// the size
    INLINE size_t Size() const { return size; }

    /// Fill array with value val
    INLINE const FlatArray & operator= (const T & val) const
    {
      for (size_t i = 0; i < size; i++)
        data[i] = val;
      return *this;
    }

    /// copies array
    INLINE const FlatArray & operator= (const FlatArray & a2) const
    {
      for (size_t i = 0; i < size; i++) data[i] = a2[i];
      return *this;
    }

    template <typename T2>
    INLINE const FlatArray & operator= (const BaseArrayObject<T2> & a2) const
    {
      T * hdata = data;
      for (size_t i = 0; i < size; i++) hdata[i] = a2.Spec()[i];
      return *this;
    }

    INLINE const FlatArray & operator= (const std::function<T(int)> & func) const
    {
      for (size_t i = 0; i < size; i++)
        data[i] = func(i);
      return *this;
    }

    /// copies pointers
    INLINE const FlatArray & Assign (const FlatArray & a2)
    {
      size = a2.size;
      data = a2.data;
      return *this;
    }

    /// assigns memory from local heap
    INLINE const FlatArray & Assign (size_t asize, LocalHeap & lh)
    {
      size = asize;
      data = lh.Alloc<T> (asize);
      return *this;
    }

    /// Access array. range check by macro CHECK_RANGE
    INLINE T & operator[] (size_t i) const
    {
#ifdef CHECK_RANGE
      if (i < 0 || i >= size)
        throw RangeException ("FlatArray::operator[]", i, 0, size-1);
#endif
      return data[i]; 
    }
  
    INLINE T_Range<size_t> Range () const
    {
      return T_Range<size_t> (0, Size());
    }
    
    INLINE const CArray<T> Addr (int pos) const
    { return CArray<T> (data+pos); }

    // const CArray<T> operator+ (int pos)
    // { return CArray<T> (data+pos); }
    INLINE T * operator+ (size_t pos) const { return data+pos; }

    /// access last element. check by macro CHECK_RANGE
    T & Last () const
    {
#ifdef CHECK_RANGE
      if (!size)
        throw Exception ("Array should not be empty");
#endif

      return data[size-1];
    }

    /// takes sub-array starting from position pos
    INLINE const FlatArray<T> Part (size_t pos)
    {
      return FlatArray<T> (size-pos, data+pos);
    }

    /// takes subsize elements starting from position pos
    INLINE const FlatArray<T> Part (size_t pos, size_t subsize)
    {
      return FlatArray<T> (subsize, data+pos);
    }

    /// takes range starting from position start of end-start elements
    INLINE FlatArray<T> Range (size_t start, size_t end) const
    {
      return FlatArray<T> (end-start, data+start);
    }

    /// takes range starting from position start of end-start elements
    INLINE FlatArray<T> Range (T_Range<size_t> range) const
    {
      return FlatArray<T> (range.Size(), data+range.First());
    }

    /// takes range starting from position start of end-start elements
    INLINE const FlatArray<T> operator[] (IntRange range) const
    {
      return FlatArray<T> (range.Size(), data+range.First());
    }
    
    template <typename TI1>
    IndirectArray<T, BaseArrayObject<TI1> > 
    operator[] (const BaseArrayObject<TI1> & ind_array) const
    {
      return IndirectArray<T, BaseArrayObject<TI1> > (*this, ind_array);
    }

    /// first position of element elem, returns -1 if element not contained in array 
    INLINE int Pos(const T & elem) const
    {
      int pos = -1;
      for(int i=0; pos==-1 && i < size; i++)
        if(elem == (*this)[i]) pos = i;
      return pos;
    }

    /// does the array contain element elem ?
    INLINE bool Contains(const T & elem) const
    {
      return ( Pos(elem) >= 0 );
    }
    
    ArrayIterator<T, size_t> begin() const
    { return ArrayIterator<T,size_t> (*this, 0); }
    ArrayIterator<T, size_t> end() const
    { return ArrayIterator<T,size_t> (*this, size); }

    /*
    FlatArray<T,TSIZE> OmpSplit() const
    {
      return Range(ngstd::OmpSplit(Range()));
    }
    */
  };



  /// print array
  template <class T>
  inline ostream & operator<< (ostream & s, const FlatArray<T> & a)
  {
    for (auto i : a.Range())
      s << i << ": " << a[i] << "\n";
    return s;
  }

  /// have arrays the same contents ?
  template <class T1, class T2>
  inline bool operator== (const FlatArray<T1> & a1,
                          const FlatArray<T2> & a2)
  {
    if (a1.Size () != a2.Size()) return 0;
    for (size_t i = 0; i < a1.Size(); i++)
      if (a1[i] != a2[i]) return false;
    return true;
  }
		 


  /** 
      Dynamic array container.
   
      Array<T> is an automatically increasing array container.
      The allocated memory doubles on overflow. 
      Either the container takes care of memory allocation and deallocation,
      or the user provides one block of data.
  */
  template <class T>
  class Array : public FlatArray<T>
  {
  protected:
    /// physical size of array
    size_t allocsize;
    /// that's the data we have to delete, nullptr for not owning the memory
    T * mem_to_delete;

    using FlatArray<T>::size;
    using FlatArray<T>::data;

  public:
    /// Generate array of logical and physical size asize
    INLINE explicit Array()
      : FlatArray<T> (0, nullptr)
    {
      allocsize = 0; 
      mem_to_delete = nullptr;
    }

    INLINE explicit Array(size_t asize)
      : FlatArray<T> (asize, new T[asize])
    {
      allocsize = asize; 
      mem_to_delete = data;
    }


    /// Generate array in user data
    INLINE Array(size_t asize, T* adata)
      : FlatArray<T> (asize, adata)
    {
      allocsize = asize; 
      mem_to_delete = nullptr;
    }

    /// Generate array in user data
    template <typename ALLOCATOR>
    INLINE Array(size_t asize, ALLOCATOR & lh)
      : FlatArray<T> (asize, lh)
    {
      allocsize = asize; 
      mem_to_delete = nullptr;
    }

    INLINE Array (Array && a2) 
    {
      size = a2.size; 
      data = a2.data;
      allocsize = a2.allocsize;
      mem_to_delete = a2.mem_to_delete;
      a2.size = 0;
      a2.allocsize = 0;
      a2.data = nullptr;
      a2.mem_to_delete = nullptr;
    }

    /// array copy 
    INLINE explicit Array (const Array & a2)
      : FlatArray<T> (a2.Size(), a2.Size() ? new T[a2.Size()] : nullptr)
    {
      allocsize = size;
      mem_to_delete = data;
      for (size_t i = 0; i < size; i++)
        (*this)[i] = a2[i];
    }

    
    template <typename TA>
    explicit Array (const BaseArrayObject<TA> & a2)
      : FlatArray<T> (a2.Size(), 
                      a2.Size() ? new T[a2.Size()] : nullptr)
    {
      allocsize = size;
      mem_to_delete = data;
      for (size_t i = 0; i < size; i++)
        (*this)[i] = a2.Spec()[i];
    }

    Array (std::initializer_list<T> list) 
      : FlatArray<T> (list.size(), 
                      list.size() ? new T[list.size()] : NULL)
    {
      allocsize = size;
      mem_to_delete = data;
      size_t cnt = 0;
      for (auto val : list)
        data[cnt++] = val;
      /*
      for (auto i = list.begin(); i < list.end(); i++, cnt++)
        data[cnt] = *i;
      */
    }

    /// array merge-copy
    explicit Array (const Array<T> & a2, const Array<T> & a3)
      : FlatArray<T> (a2.Size()+a3.Size(), 
                      a2.Size()+a3.Size() ? new T[a2.Size()+a3.Size()] : 0)
    {
      allocsize = size;
      mem_to_delete = data;
      for(size_t i = 0; i <  a2.Size(); i++)
        (*this)[i] = a2[i];
      for (size_t i = a2.Size(), j=0; i < size; i++,j++)
        (*this)[i] = a3[j];
    }

    /// if responsible, deletes memory
    INLINE ~Array()
    {
      delete [] mem_to_delete;
    }

    /// we tell the compiler that there is no need for deleting the array ..
    INLINE void NothingToDelete () 
    { 
      mem_to_delete = nullptr;
    }

    /// Change logical size. If necessary, do reallocation. Keeps contents.
    INLINE void SetSize(size_t nsize)
    {
      if (nsize > allocsize) ReSize (nsize);
      size = nsize; 
    }

    ///
    INLINE void SetSize0()
    {
      size = 0; 
    }

    /// Change physical size. Keeps logical size. Keeps contents.
    INLINE void SetAllocSize (size_t nallocsize)
    {
      if (nallocsize > allocsize)
        ReSize (nallocsize);
    }

    /// Change physical size. Keeps logical size. Keeps contents.
    INLINE size_t AllocSize () const
    {
      return allocsize;
    }


    /// assigns memory from local heap
    INLINE const Array & Assign (size_t asize, LocalHeap & lh)
    {
      delete [] mem_to_delete;
      size = allocsize = asize;
      data = lh.Alloc<T> (asize);
      mem_to_delete = nullptr;
      return *this;
    }

    /// Add element at end of array. reallocation if necessary.
    INLINE size_t Append (const T & el)
    {
      if (size == allocsize) 
        ReSize (size+1);
      data[size] = el;
      size++;
      return size;
    }

    /// Add element at end of array. reallocation if necessary.
    INLINE size_t Append (T && el)
    {
      if (size == allocsize) 
        ReSize (size+1);
      data[size] = move(el);
      size++;
      return size;
    }


    INLINE Array<T> & operator += (const T & el)
    {
      Append (el);
      return *this;
    }


    /// Append array at end of array. reallocation if necessary.
    // int Append (const Array<T> & source)
    INLINE size_t Append (FlatArray<T> source)
    {
      if(size + source.Size() >= allocsize)
        ReSize (size + source.Size() + 1);

      size_t i,j;
      for(i = size, j=0; j<source.Size(); i++, j++)
        data[i] = source[j];

      size += source.Size();
      return size;
    }



    /// Delete element i. Move last element to position i.
    INLINE void DeleteElement (size_t i)
    {
#ifdef CHECK_RANGE
      // RangeCheck (i);
#endif

      data[i] = data[size-1];
      size--;
    }


    /// Delete element i. Move all remaining elements forward
    INLINE void RemoveElement (size_t i)
    {
      for(size_t j = i; j < this->size-1; j++)
	this->data[j] = this->data[j+1];
      this->size--;
    }


    /// Delete last element. 
    INLINE void DeleteLast ()
    {
#ifdef CHECK_RANGE
      //    CheckNonEmpty();
#endif
      size--;
    }

    /// Deallocate memory
    INLINE void DeleteAll ()
    {
      delete [] mem_to_delete;
      mem_to_delete = NULL;
      data = 0;
      size = allocsize = 0;
    }

    /// Fill array with val
    INLINE Array & operator= (const T & val)
    {
      FlatArray<T>::operator= (val);
      return *this;
    }

    /// array copy
    INLINE Array & operator= (const Array & a2)
    {
      SetSize (a2.Size());
      for (size_t i = 0; i < size; i++)
        (*this)[i] = a2[i];
      return *this;
    }

    /// steal array 
    INLINE Array & operator= (Array && a2)
    {
      ngstd::Swap (size, a2.size);
      ngstd::Swap (data, a2.data);
      ngstd::Swap (allocsize, a2.allocsize);
      ngstd::Swap (mem_to_delete, a2.mem_to_delete);
      return *this;
    }


    /// array copy
    INLINE Array & operator= (const FlatArray<T> & a2)
    {
      SetSize (a2.Size());
      for (size_t i = 0; i < size; i++)
        (*this)[i] = a2[i];
      return *this;
    }

    /*
    /// fill array with first, first+1, ... 
    Array & operator= (const IntRange & range)
    {
      SetSize (range.Size());
      for (int i = 0; i < size; i++)
        (*this)[i] = range.First()+i;
      return *this;
    }
    */
    template <typename T2>
    Array & operator= (const BaseArrayObject<T2> & a2)
    {
      SetSize (a2.Spec().Size());
      for (size_t i = 0; i < size; i++)
        (*this)[i] = a2.Spec()[i];
      return *this;
    }

    template <typename ...ARGS>
    Array & operator= (Tuple<ARGS...> tup)
    {
      SetSize (ArraySize (tup));
      StoreToArray (*this, tup);
      return *this;
    }

    Array & operator= (std::initializer_list<T> list)
    {
      *this = Array<T> (list); 
      return *this;
    }
      

    INLINE void Swap (Array & b)
    {
      ngstd::Swap (size, b.size);
      ngstd::Swap (data, b.data);
      ngstd::Swap (allocsize, b.allocsize);
      // ngstd::Swap (ownmem, b.ownmem);
      ngstd::Swap (mem_to_delete, b.mem_to_delete);
    }

  private:

    /// resize array, at least to size minsize. copy contents
    INLINE void ReSize (size_t minsize);
    INLINE void ResetSize (size_t new_size);
     /*
    {
      TSIZE nsize = 2 * allocsize;
      if (nsize < minsize) nsize = minsize;

      if (data)
        {
          T * p = new T[nsize];
	
          TSIZE mins = (nsize < size) ? nsize : size; 
          memcpy (p, data, mins * sizeof(T));

          if (ownmem) delete [] data;
          ownmem = 1;
          data = p;
        }
      else
        {
          data = new T[nsize];
          ownmem = 1;
        }
    
      allocsize = nsize;
    }
    */
  };

  template <class T> 
  INLINE void Array<T> :: ResetSize (size_t new_size)
  {
    if(new_size<=allocsize)
      {
	size = new_size;
	return;
      }

    delete [] mem_to_delete;
    data = new T[new_size];
    allocsize = new_size;
    size = new_size;
    mem_to_delete = data;
  }
    
  /// resize array, at least to size minsize. copy contents
  template <class T> 
  INLINE void Array<T> :: ReSize (size_t minsize)
  {
    size_t nsize = 2 * allocsize;
    if (nsize < minsize) nsize = minsize;
    
    T * hdata = data;
    data = new T[nsize];

    if (hdata)
      {
        size_t mins = (nsize < size) ? nsize : size;
#if defined(__GNUG__) && __GNUC__ < 5 && !defined(__clang__)
        for (size_t i = 0; i < mins; i++) data[i] = move(hdata[i]);
#else
        if (std::is_trivially_copyable<T>::value)
          memcpy (data, hdata, sizeof(T)*mins);
        else
          for (size_t i = 0; i < mins; i++) data[i] = move(hdata[i]);
#endif
        delete [] mem_to_delete;
      }

    mem_to_delete = data;
    allocsize = nsize;
  }

  //extern template class Array<int,int>;
  

  /**
     Array with static and dynamic memory management.
     Declares a static array which size is given by the template parameter.
     If the dynamic size fits into the static size, use static memory, 
     otherwise perform dynamic allocation
  */
  template <class T, int S> 
  class ArrayMem : public Array<T>
  {
    T mem[S];    

    using Array<T>::size;
    using Array<T>::allocsize;
    using Array<T>::data;
    using Array<T>::mem_to_delete;
    // using Array<T>::ownmem;

  public:
    /// Generate array of logical and physical size asize
    explicit ArrayMem(size_t asize = 0)    
      : Array<T> (S, mem)
    {
      size = asize;
      if (asize > S)
        {
          data = new T[asize];
          allocsize = size;
          // ownmem = 1;
          mem_to_delete = data;
        }
    }

    /// copies from Array a2
    explicit ArrayMem(const Array<T> & a2)
      : Array<T> (S, (T*)mem)
    {
      Array<T>::operator= (a2);
    }

    /// copies from ArrayMem a2
    explicit ArrayMem(const ArrayMem & a2)
      : Array<T> (S, (T*)mem)
    {
      Array<T>::operator= (a2);
    }
  

    ArrayMem & operator= (const T & val)
    {
      FlatArray<T>::operator= (val);
      return *this;
    }

    /// array copy
    ArrayMem & operator= (const FlatArray<T> & a2)
    {
      this->SetSize (a2.Size());
      for (int i = 0; i < size; i++)
        (*this)[i] = a2[i];
      return *this;
    }


    template <typename T2>
    ArrayMem & operator= (const BaseArrayObject<T2> & a2)
    {
      this->SetSize (a2.Spec().Size());
      for (int i = 0; i < size; i++)
        (*this)[i] = a2.Spec()[i];
      return *this;
    }

  };





  template <typename ... ARGS>
  int ArraySize (Tuple<ARGS...> tup)  
  { return 0;}
  
  template <typename ... ARGS>
  int ArraySize (Tuple<int,ARGS...> tup) 
  { return 1+ArraySize(tup.Tail()); }
  
  template <typename ... ARGS>
  int ArraySize (Tuple<IntRange,ARGS...> tup) 
  { return tup.Head().Size()+ArraySize(tup.Tail()); }

  
  template <typename ... ARGS>
  void StoreToArray (FlatArray<int> a, Tuple<ARGS...> tup) { ; }
  
  template <typename ... ARGS>
  void StoreToArray (FlatArray<int> a, Tuple<int,ARGS...> tup)
  {
    a[0] = tup.Head();
    StoreToArray (a.Range(1, a.Size()), tup.Tail());
  }
  
  template <typename ... ARGS>
  void StoreToArray (FlatArray<int> a, Tuple<IntRange,ARGS...> tup)
  {
    IntRange r = tup.Head();
    a.Range(0,r.Size()) = r;
    StoreToArray (a.Range(r.Size(), a.Size()), tup.Tail());
  }

  /*
  template <typename T> template <typename ...ARGS>
  INLINE Array<T> & Array<T> :: operator= (Tuple<ARGS...> tup)
  {
    SetSize (ArraySize (tup));
    StoreToArray (*this, tup);
  }
  */

  /*
  /// append integers to array
  inline Array<int> & operator+= (Array<int> & array, const IntRange & range)
  {
    int oldsize = array.Size();
    int s = range.Next() - range.First();
    
    array.SetSize (oldsize+s);

    for (int i = 0; i < s; i++)
      array[oldsize+i] = range.First()+i;

    return array;
  }
  */
  

  template <typename T2>
  inline Array<int> & operator+= (Array<int> & array, const BaseArrayObject<T2> & a2)
  {
    int oldsize = array.Size();
    int s = a2.Spec().Size();
    
    array.SetSize (oldsize+s);

    for (int i = 0; i < s; i++)
      array[oldsize+i] = a2.Spec()[i];

    return array;
  }




  /// bubble sort array
  template <class T>
  inline void BubbleSort (const FlatArray<T> & data)
  {
    T hv;
    for (int i = 0; i < data.Size(); i++)
      for (int j = i+1; j < data.Size(); j++)
        if (data[i] > data[j])
          {
            hv = data[i];
            data[i] = data[j];
            data[j] = hv;
          }
  }

  /// bubble sort array
  template <class T, class S>
  inline void BubbleSort (FlatArray<T> & data, FlatArray<S> & slave)
  {
    for (int i = 0; i < data.Size(); i++)
      for (int j = i+1; j < data.Size(); j++)
	if (data[i] > data[j])
	  {
	    T hv = data[i];
	    data[i] = data[j];
	    data[j] = hv;

	    S hvs = slave[i];
	    slave[i] = slave[j];
	    slave[j] = hvs;
	  }
  }




  template <class T, typename TLESS>
  void QuickSort (FlatArray<T> data, TLESS less)
  {
    if (data.Size() <= 1) return;

    int i = 0;
    int j = data.Size()-1;

    T midval = data[ (i+j)/2 ];
  
    do
      {
        while (less (data[i], midval)) i++;
        while (less (midval, data[j])) j--;

        if (i <= j)
          {
	    Swap (data[i], data[j]);
            i++; j--;
          }
      }
    while (i <= j);

    QuickSort (data.Range (0, j+1), less);
    QuickSort (data.Range (i, data.Size()), less);
  }

  template <typename T>
  INLINE bool DefaultLess (const T & a, const T & b)
  {
    return a < b;
  }

  template <typename T>
  class DefaultLessCl
  {
  public:
    bool operator() (const T & a, const T & b) const
    {
      return a < b;
    }
  };



  template <class T>
  INLINE void QuickSort (FlatArray<T> data)
  {
    QuickSort (data, DefaultLessCl<T>());
  }



  template <class T, typename TLESS>
  void QuickSortI (FlatArray<T> data, FlatArray<int> index, TLESS less)
  {
    if (index.Size() <= 1) return;

    int i = 0;
    int j = index.Size()-1;

    int midval = index[ (i+j)/2 ];
  
    do
      {
        while (less (data[index[i]],data[midval])  ) i++;
        while (less (data[midval],  data[index[j]])) j--;

        if (i <= j)
          {
	    Swap (index[i], index[j]);
            i++; j--;
          }
      }
    while (i <= j);

    QuickSortI (data, index.Range (0, j+1), less);
    QuickSortI (data, index.Range (i, index.Size()), less);
  }


  template <class T>
  INLINE void QuickSortI (FlatArray<T> data, FlatArray<int> index)
  {
    QuickSortI (data, index, DefaultLessCl<T>());
  }





  template <typename T>
  INLINE T xxxRemoveRef (const T & x)
  {
    return x;
  }

  template <class TA1, class TA2> 
  class SumArray : public BaseArrayObject<SumArray<TA1,TA2>>
  {
    const TA1 & a1;
    const TA2 & a2;
  public:
    SumArray (const TA1 & aa1, const TA2 & aa2) : a1(aa1), a2(aa2) { ; }
    int Size() const { return a1.Size()+a2.Size(); }
    auto operator[] (int i) const -> decltype (xxxRemoveRef (a1[0])) 
    {
      return (i < a1.Size()) ? a1[i] : a2[i-a1.Size()];
    }
  };

  template <class TA1, class TA2> 
  SumArray<TA1,TA2> operator+ (const BaseArrayObject<TA1> & a1,
                               const BaseArrayObject<TA2> & a2)
  {
    return SumArray<TA1,TA2> (a1.Spec(), a2.Spec());
  }
                               

  // head-tail array
  template <size_t S, typename T>
  class HTArray
  {
    HTArray<S-1,T> tail;
    T head;
  public:
    HTArray () = default;
    HTArray (const HTArray &) = default;
    HTArray & operator= (const HTArray &) = default;

    T * Ptr () { return tail.Ptr(); }
    T & operator[] (size_t i) { return Ptr()[i]; }

    const T * Ptr () const { return tail.Ptr(); }
    const T & operator[] (size_t i) const { return Ptr()[i]; }
  };

  template <typename T>
  class HTArray<1,T>
  {
    T head;
  public:
    HTArray () = default;
    HTArray (const HTArray &) = default;
    HTArray & operator= (const HTArray &) = default;

    T * Ptr () { return &head; }
    T & operator[] (size_t i) { return Ptr()[i]; }

    const T * Ptr () const { return &head; }
    const T & operator[] (size_t i) const { return Ptr()[i]; }
  };

  template <typename T>
  class HTArray<0,T>
  {
    T head; // dummy variable
  public:
    HTArray () = default;
    HTArray (const HTArray &) = default;
    HTArray & operator= (const HTArray &) = default;

    T * Ptr () { return &head; }
    T & operator[] (size_t i) { return Ptr()[i]; }

    const T * Ptr () const { return &head; }
    const T & operator[] (size_t i) const { return Ptr()[i]; }
  };

  template<size_t S, typename T>
  const T * operator+ (const HTArray<S,T> & ar, size_t i)
  {
    return ar.Ptr()+i;
  }
  template<size_t S, typename T>
  T * operator+ (HTArray<S,T> & ar, size_t i)
  {
    return ar.Ptr()+i;
  }
  

  template <typename T> 
  Archive & operator & (Archive & archive, Array<T> & a)
  {
    if (archive.Output())
      archive << a.Size();
    else
      {
        size_t size;
        archive & size;
        a.SetSize (size);
      }

    /*
    for (auto & ai : a)
      archive & ai;
    */
    archive.Do (&a[0], a.Size());
    return archive;
  }
}


#endif

