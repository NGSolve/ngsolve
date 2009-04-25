#ifndef FILE_SYMMETRICMATRIX
#define FILE_SYMMETRICMATRIX

/****************************************************************************/
/* File:   symmetricmatrix.hpp                                              */
/* Author: Joachim Schoeberl                                                */
/* Date:   12. Dec. 2004                                                    */
/****************************************************************************/


namespace ngbla
{



  /// A symmetric band-matrix
  template <class T>
  class FlatSymmetricMatrix
  {
  protected:
    ///
    int n;
    ///
    T *data;
  public:
    typedef typename mat_traits<T>::TV_COL TV;

    ///
    FlatSymmetricMatrix () { ; }

    ///
    FlatSymmetricMatrix (int an, T * adata)
      : n(an), data(adata)
    { ; }

    /// set size, and assign mem
    void AssignMemory (int an, T * mem) throw()
    {
      n = an;
      data = mem;
    }
  


    ///
    void Mult (const FlatVector<TV> & x, FlatVector<TV> & y) const
    {
      for (int i = 0; i < n; i++)
        {
          TV xi = x(i);
          TV sum = (*this)(i,i) * xi;
          for (int j = 0; j < i; j++)
            {
              sum += (*this)(i,j) * x(j);
              y(j) += Trans((*this)(i,j)) * xi;
            }
          y(i) = sum;
        }
    }

    ///
    ostream & Print (ostream & ost) const
    {
      for (int i = 0; i < n; i++)
        {
          for (int j = 0; j < n; j++)
            if (j <= i)
              ost << setw(8) << (*this)(i,j) << " ";
            else 
              ost << setw(8) << "sym" << " ";
          ost << endl;
        }
      return *this;
    }

    ///
    int Height() const { return n; }

  
    ///
    const T & operator() (int i, int j) const
    { return data[ i * (i+1) / 2 + j ]; }

    ///
    T & operator() (int i, int j) 
    { return data[ i * (i+1) / 2 + j ]; }

  
    ///
    FlatSymmetricMatrix & operator= (const T & val)
    {
      int nel = n * (n+1) / 2;
      for (int i = 0; i < nel; i++)
        data[i] = val;
      return *this;
    }
  };


  template<typename T>
  inline std::ostream & operator<< (std::ostream & s, const FlatSymmetricMatrix<T> & m)
  {
    m.Print (s);
    return s;
  }



  ///  A symmetric band-matrix with memory management
  template <class T>
  class SymmetricMatrix : public FlatSymmetricMatrix<T>
  {
  public:
    typedef typename mat_traits<T>::TV_COL TV;

    ///
    SymmetricMatrix (int an)
      : FlatSymmetricMatrix<T> (an, new T[an*(an+1)/2])
    { ; }

    ///
    ~SymmetricMatrix ()
    { delete [] this->data; }

    ///
    SymmetricMatrix & operator= (const T & val)
    {
      int nel = this->n * (this->n+1) / 2;
      for (int i = 0; i < nel; i++)
        this->data[i] = val;
      return *this;
    }
  };







}


#endif
