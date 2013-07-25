#ifndef FILE_CHOLESKY
#define FILE_CHOLESKY

/****************************************************************************/
/* File:   cholesky.hpp                                                     */
/* Author: Joachim Schoeberl                                                */
/* Date:   25. Mar. 2000, 16. June 2002                                     */
/****************************************************************************/

namespace ngbla
{

  /**
     The Cholesky-factorization of a symmetric dense matrix.
     A = L D L^T
  */
  template <class T>
  class FlatCholeskyFactors
  {
  protected:
    /// matrix size
    int n;
    /// left factor
    T * lfact;
    /// inverse diagonal
    T * diag;
  public:
    // typedef typename mat_traits<T>::TV_COL TV;
    /// Factor the matrix A
    FlatCholeskyFactors (const FlatMatrix<T> & a, T * data)
    {
      diag = data;
      Factor (a);
    }

    /// Factor the matrix A
    FlatCholeskyFactors (const FlatMatrix<T> & a, LocalHeap & lh)
    {
      diag = (T*)lh.Alloc(sizeof(T)*RequiredMem(a.Height()));
      Factor (a);
    }

    ///
    NGS_DLL_HEADER void Factor (const FlatMatrix<T> & a);
    /// Multiply with the inverse of A 
    template <typename TV1, typename TV2>
    // NGS_DLL_HEADER void Mult (SliceVector<TV> x, SliceVector<TV> y) const
    NGS_DLL_HEADER void Mult (TV1 && x, TV2 && y) const
    {   
      // TV sum, val;
      // decltype (y(0)) sum, val;


      const T *pj;

      for (int i = 0; i < n; i++)
        y(i) = x(i);
      
      for (int i = 0; i < n; i++)
        {
          auto sum = y(i);
          
          pj = PRow(i);
          for (int j = 0; j < i; ++j)
            sum -= pj[j] * y(j);
          
          y(i) = sum;
        }
      
      for (int i = 0; i < n; i++)
        {
          auto sum = diag[i] * y(i);
          y(i) = sum;
      }
      
      for (int i = n-1; i >= 0; i--)
        {
          pj = PRow(i);
          auto val = y(i);
          for (int j = 0; j < i; ++j)
            y(j) -= pj[j] * val;
        }
    }
    
    /// Print factorization
    NGS_DLL_HEADER ostream & Print (ostream & ost) const;


    /// computes required memory
    static int RequiredMem (int n)
    { return n*(n+1)/2; }

  private:
    /// first element in row
    T * PRow (int i) const { return lfact + (i*(i-1)) / 2; }
  };


  ///  output operator.
  template<typename T>
  inline std::ostream & operator<< (std::ostream & s, const FlatCholeskyFactors<T> & m)
  {
    m.Print (s);
    return s;
  }



  template <class T>
  class CholeskyFactors : public FlatCholeskyFactors<T>
  {
  public:
    /// Factor the matrix A
    CholeskyFactors (const FlatMatrix<T> & a)
      : FlatCholeskyFactors<T> (a, new T[this->RequiredMem(a.Height())])
    { ; }
    /// Delete memory
    ~CholeskyFactors ()
    {
      delete [] this->diag;
    }
  };

}

#endif
