#ifndef FILE_NGS_Array
#define FILE_NGS_Array

/**************************************************************************/
/* File:   array.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/


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
    CArray () { data = 0; }

    /// provide size and memory
    CArray (T * adata) 
      : data(adata) { ; }

    /// Access array
    T & operator[] (int i) const
    {
      return data[i]; 
    }

    operator T* () const { return data; }
  };







  /**
     A simple array container.
     Array represented by size and data-pointer.
     No memory allocation and deallocation, must be provided by user.
     Helper functions for printing. 
     Optional range check by macro RANGE_CHECK
  */
  template <class T>
  class FlatArray
  {
  protected:
    /// the size
    int size;
    /// the data
    T * data;
  public:

    /// initialize array 
    FlatArray () { ; } // size = 0; data = 0; }

    /// provide size and memory
    FlatArray (int asize, T * adata) 
      : size(asize), data(adata) { ; }

    /// memory from local heap
    FlatArray(int asize, LocalHeap & lh)
      : size(asize),
        // data(static_cast<T*> (lh.Alloc(asize*sizeof(T))))
        data (lh.Alloc<T> (asize))
    { ; }

    /// the size
    int Size() const { return size; }


    /// Fill array with value val
    const FlatArray & operator= (const T & val) const
    {
      for (int i = 0; i < size; i++)
        data[i] = val;
      return *this;
    }

    /// copies pointers
    const FlatArray & operator= (const FlatArray & a2)
    {
      size = a2.size;
      data = a2.data;
      return *this;
    }


    /*
   /// access array. range check by macro CHECK_RANGE
   T & operator[] (int i) 
   { 

   #ifdef CHECK_RANGE
   if (i < 0 || i >= size)
   throw RangeException ("FlatArray::operator[]", i, 0, size-1);
   #endif

   return data[i]; 
   }
    */

    /// Access array. range check by macro CHECK_RANGE
    T & operator[] (int i) const
    {

#ifdef CHECK_RANGE
      if (i < 0 || i >= size)
        throw RangeException ("FlatArray::operator[]", i, 0, size-1);
#endif

      return data[i]; 
    }
  

    const CArray<T> Addr (int pos)
    {
      return CArray<T> (data+pos);
    }

    /*
    // icc does not like this version
    T * const Addr (int i) const
    {
    return data+i;
    }
    */


    /*
   /// access last element. check by macro CHECK_RANGE
   T & Last ()
   {
   #ifdef CHECK_RANGE
   if (!size)
   throw Exception ("Array should not be empty");
   #endif

   return data[size-1];
   }
    */

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
    const FlatArray<T> Part (int pos)
    {
      return FlatArray<T> (size-pos, data+pos);
    }

    /// takes subsize elements starting from position pos
    const FlatArray<T> Part (int pos, int subsize)
    {
      return FlatArray<T> (subsize, data+pos);
    }

    /// takes range starting from position start of end-start elements
    const FlatArray<T> Range (int start, int end) const
    {
      return FlatArray<T> (end-start, data+start);
    }

    /// first position of element elem, returns -1 if element not contained in array 
    int Pos(const T & elem) const
    {
      int pos = -1;
      for(int i=0; pos==-1 && i < size; i++)
        if(elem == (*this)[i]) pos = i;
      return pos;
    }

    /// does the array contain element elem ?
    bool Contains(const T & elem) const
    {
      return ( Pos(elem) >= 0 );
    }

  };



  /// print array
  template <class T>
  ostream & operator<< (ostream & s, const FlatArray<T> & a)
  {
    for (int i = 0; i < a.Size(); i++)
      s << i << ": " << a[i] << "\n";
    return s;
  }

  /// have arrays the same contents ?
  template <class T1, class T2>
  inline bool operator== (const FlatArray<T1> & a1,
                          const FlatArray<T2> & a2)
  {
    if (a1.Size () != a2.Size()) return 0;
    for (int i = 0; i < a1.Size(); i++)
      if (a1[i] != a2[i]) return 0;
    return 1;
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
    int allocsize;
    /// memory is responsibility of container
    bool ownmem;

    using FlatArray<T>::size;
    using FlatArray<T>::data;
  public:
    /// Generate array of logical and physical size asize
    explicit Array(int asize = 0)
      : FlatArray<T> (asize, asize ? new T[asize] : 0)
    {
      allocsize = asize; 
      ownmem = 1;
    }

    /// Generate array in user data
    Array(int asize, T* adata)
      : FlatArray<T> (asize, adata)
    {
      allocsize = asize; 
      ownmem = 0;
    }

    /// array copy 
    explicit Array (const Array<T> & a2)
      : FlatArray<T> (a2.Size(), a2.Size() ? new T[a2.Size()] : 0)
    {
      allocsize = size;
      ownmem = 1;
      for (int i = 0; i < size; i++)
        (*this)[i] = a2[i];
    }

    /// array merge-copy
    explicit Array (const Array<T> & a2, const Array<T> & a3)
      : FlatArray<T> (a2.Size()+a3.Size(), 
                      a2.Size()+a3.Size() ? new T[a2.Size()+a3.Size()] : 0)
    {
      allocsize = size;
      ownmem = 1;
      for(int i = 0; i <  a2.Size(); i++)
        (*this)[i] = a2[i];
      for (int i = a2.Size(), j=0; i < size; i++,j++)
        (*this)[i] = a3[j];
    }

    /// if responsible, deletes memory
    ~Array()
    {
      if (ownmem) delete [] data;
    }

    /// Change logical size. If necessary, do reallocation. Keeps contents.
    void SetSize(int nsize)
    {
      if (nsize > allocsize) ReSize (nsize);
      size = nsize; 
    }

    /// Change physical size. Keeps logical size. Keeps contents.
    void SetAllocSize (int nallocsize)
    {
      if (nallocsize > allocsize)
        ReSize (nallocsize);
    }

    /// Add element at end of array. reallocation if necessary.
    int Append (const T & el)
    {
      if (size == allocsize) 
        ReSize (size+1);
      data[size] = el;
      size++;
      return size;
    }

    Array<T> & operator += (const T & el)
    {
      Append (el);
      return *this;
    }


    /// Append array at end of array. reallocation if necessary.
    int Append (const Array<T> & source)
    {
      if(size + source.Size() >= allocsize)
        ReSize (size + source.Size() + 1);

      int i,j;
      for(i = size, j=0; j<source.Size(); i++, j++)
        data[i] = source[j];

      size += source.Size();
      return size;
    }



    /// Delete element i. Move last element to position i.
    void DeleteElement (int i)
    {
#ifdef CHECK_RANGE
      // RangeCheck (i);
#endif

      data[i] = data[size-1];
      size--;
    }

    /// Delete last element. 
    void DeleteLast ()
    {
#ifdef CHECK_RANGE
      //    CheckNonEmpty();
#endif

      size--;
    }

    /// Deallocate memory
    void DeleteAll ()
    {
      if (ownmem) delete [] data;
      data = 0;
      size = allocsize = 0;
    }

    /// Fill array with val
    Array & operator= (const T & val)
    {
      FlatArray<T>::operator= (val);
      return *this;
    }

    /// array copy
    Array & operator= (const Array & a2)
    {
      SetSize (a2.Size());
      for (int i = 0; i < size; i++)
        (*this)[i] = a2[i];
      return *this;
    }

    /// array copy
    Array & operator= (const FlatArray<T> & a2)
    {
      SetSize (a2.Size());
      for (int i = 0; i < size; i++)
        (*this)[i] = a2[i];
      return *this;
    }



  private:

    /// resize array, at least to size minsize. copy contents
    void ReSize (int minsize)
    {
      int nsize = 2 * allocsize;
      if (nsize < minsize) nsize = minsize;

      if (data)
        {
          T * p = new T[nsize];
	
          int mins = (nsize < size) ? nsize : size; 
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
  };


  /**
     Array with static and dynamic memory management.
     Declares a static array which size is given by the template parameter.
     If the dynamic size fits into the static size, use static memory, 
     otherwise perform dynamic allocation
  */
  template <class T, int S> 
  class ArrayMem : public Array<T>
  {
    // T mem[S];                     // should be best, but calls trivial default constructor 
    // char mem[S*sizeof(T)];     // avoids calling the array default-constructor (icc)
    // static memory
    double mem[(S*sizeof(T)+7) / 8];   // alignment (on ia64 machines)

    using Array<T>::size;
    using Array<T>::data;
    using Array<T>::ownmem;

  public:
    /// Generate array of logical and physical size asize
    explicit ArrayMem(int asize = 0)    
      : Array<T> (S, static_cast<T*> (static_cast<void*>(&mem[0])))
    {
      size = asize;
      if (asize > S)
        {
          data = new T[asize];
          ownmem = 1;
        }
      // this->SetSize (asize);
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
    ArrayMem & operator= (const Array<T> & a2)
    {
      SetSize (a2.Size());
      for (int i = 0; i < size; i++)
        (*this)[i] = a2[i];
      return *this;
    }

  };







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





  template <class T>
  void QuickSort (const FlatArray<T> & data)
  {
    if (data.Size() <= 1) return;

    int i = 0;
    int j = data.Size()-1;

    T midval = data[ (i+j)/2 ];
  
    do
      {
        while (data[i] < midval) i++;
        while (midval < data[j]) j--;
      
        if (i <= j)
          {
            T hv = data[i];
            data[i] = data[j];
            data[j] = hv;
            i++; j--;
          }
      }
    while (i <= j);

    QuickSort (data.Range (0, j+1));
    QuickSort (data.Range (i, data.Size()));

    for (int i = 0; i < data.Size()-1; i++)
      if (data[i] > data[i+1]) cerr << "heap sort is wrong !!" << endl;
  }


}


#endif

