#ifndef FILE_BANDMATRIX
#define FILE_BANDMATRIX

/****************************************************************************/
/* File:   bandmatrix.hpp                                                   */
/* Author: Joachim Schoeberl                                                */
/* Date:   14. Aug. 2002                                                    */
/****************************************************************************/

namespace ngbla
{

  /**
     A symmetric band-matrix.
  */
  template <class T = double>
  class FlatSymBandMatrix
  {
  protected:
    /// the matrix size
    int n;
    /// number of bands in the triangular matrix
    int bw;
    /// the matrix entries
    T *data;
  public:
    /// the according vector type
    typedef typename mat_traits<T>::TV_COL TV;

    /// Construction of FlatSymBandMatirx
    FlatSymBandMatrix (int an, int abw, T * adata)
      : n(an), bw(abw), data(adata)
    { ; }

    /// Matrix vector multiplication
    void Mult (const FlatVector<TV> & x, FlatVector<TV> & y) const
    {
      for (int i = 0; i < n; i++)
        y(i) = (*this)(i,i) * x(i);
      for (int i = 0; i < n; i++)
        for (int j = max2(i-bw+1, 0); j < i; j++)
          {
            y(i) += (*this)(i,j) * x(j);
            y(j) += Trans((*this)(i,j)) * x(i);
          }
    }

    /// Print matrix to stream
    ostream & Print (ostream & ost) const
    {
      for (int i = 0; i < n; i++)
        {
          for (int j = 0; j < n; j++)
            if (Used (i, j))
              ost << setw(8) << (*this)(i,j) << " ";
            else if (Used (j,i))
              ost << setw(8) << "sym" << " ";
            else
              ost << setw(8) << 0;
          ost << endl;
        }
      return ost;
    }

    /// matrix height
    int Height() const { return n; }

    /// matrix width
    int BandWidth() const { return bw; }
  
    /// access operator. Assumes that $i \geq j$ referes to an element in the band
    const T & operator() (int i, int j) const
    { return data[i * bw + j - i + bw-1]; }

    /// access operator. Assumes that $i \geq j$ referes to an element in the band
    T & operator() (int i, int j) 
    { return data[i * bw + j - i + bw-1]; }

    /// check whether i, j refers to a valid element
    bool Used (int i, int j) const
    {
      return (n > i && i >= j && j >= 0 && i-j < bw);
    }
  
    /// assigns a constant value
    FlatSymBandMatrix & operator= (const T & val)
    {
      for (int i = 0; i < bw * n; i++)
        data[i] = val;
      return *this;
    }


    /// computes required memory
    static int RequiredMem (int n, int bw)
    { return n*bw; }
  };



  /// output operator
  template<typename T>
  inline std::ostream & operator<< (std::ostream & s, const FlatSymBandMatrix<T> & m)
  {
    m.Print (s);
    return s;
  }




  /**
     A symmetric band-matrix with memory management.
  */
  template <class T = double>
  class SymBandMatrix : public FlatSymBandMatrix<T>
  {
  public:
    typedef typename mat_traits<T>::TV_COL TV;

    /// Generates a symmetric band matrix
    SymBandMatrix (int an, int abw)
      : FlatSymBandMatrix<T> (an, abw, new T[an*abw])
    { ; }

    /// Deletes matrix
    ~SymBandMatrix ()
    { delete [] this->data; }

    /// assigns a constant value
    SymBandMatrix & operator= (const T & val)
    {
      for (int i = 0; i < this->bw * this->n; i++)
        this->data[i] = val;
      return *this;
    }
  };













  /**
     Cholesky factors of a band matrix.
     This class does not provide memory management.


     storage:

     lfact (bw = 3)

     \begin{verbatim}
     d0                
     0   d1            
     1    2   d2        
     3    4   d3   
     \end{verbatim}
  */
  template <class T = double>
  class FlatBandCholeskyFactors
  {
  protected:
    /// matrix dimension
    int n;
    /// number of bands in the triangular matrix
    int bw;
    /// matrix matrix data, first diags, than lfact
    T * mem; 
  public:
    // typedef typename mat_traits<T>::TV_COL TV;
    typedef typename mat_traits<T>::TSCAL TSCAL;

    /// assign dimension, bandwidth and memory
    FlatBandCholeskyFactors (int an, int abw, T * amem)
    { n = an, bw = abw, mem = amem; }

    /// default constructor
    FlatBandCholeskyFactors ()
    { n = bw = 0; mem = 0; }

    /// factor bandmatrix a
    NGS_DLL_HEADER void Factor (const FlatSymBandMatrix<T> & a);

    /// solve with factored matrices
    template <class TVX, class TVY>
    void Mult (const FlatVector<TVX> & x, FlatVector<TVY> & y) const
    {
      // const TVX * hx = x.Addr(0);
      // TVY * hy = y.Addr(0);
      FlatVector<TVX> hx = x;
      FlatVector<TVY> hy = y;
      const T * hm = &mem[0];

      for (int i = 0; i < n; i++)
        hy[i] = hx[i];

      int i, jj = n;
      for (i = 0; i < bw-1; i++)
        {
          typedef typename mat_traits<TVY>::TSCAL TTSCAL;
          TVY sum = TTSCAL(0.0);  

          for (int j = 0; j < i; j++, jj++)
            sum += hm[jj] * hy[j];

          hy[i] -= sum;
        }

      for (  ; i < n; i++)
        {
          typedef typename mat_traits<TVY>::TSCAL TTSCAL;
          TVY sum = TTSCAL(0.0);

          for (int j = i-bw+1; j < i; j++, jj++)
            sum += hm[jj] * hy[j];

          hy[i] -= sum;
        }

      for (int i = 0; i < n; i++)
        {
          TVY sum = mem[i] * hy[i];
          hy[i] = sum;
        }

      // jj = n + (n-1) * (bw-1) - bw*(bw-1)/2;   
      for (i = n-1; i >= bw-1; i--)
        {
          jj -= bw-1;
          TVY val = hy[i];

          int firstj = i-bw+1;
          for (int j = 0; j < bw-1; j++)
            hy[firstj+j] -= Trans (mem[jj+j]) * val;
        }
    
      for (  ; i >= 0; i--) 
        {
          jj -= i;
          TVY val = hy[i];

          for (int j = 0; j < i; j++)
            hy[j] -= Trans (mem[jj+j]) * val;
        }
    }



    /// print matrix factors
    ostream & Print (ostream & ost) const;

    /// compute linear position of matrix element (i,j)
    int Index (int i, int j) const
    {
      if (i < bw)
        return n + (i * (i-1)) / 2 + j;
      else
        return n + i * (bw-2) + j - ((bw-1)*(bw-2))/2;   
    }
  
    /// matrix element (i,j), (i,j) must be a valid position
    const T & operator() (int i, int j) const
    { 
      if (i < bw)
        return mem[n + (i * (i-1)) / 2 + j];
      else
        return mem[n + i * (bw-2) + j - ((bw-1)*(bw-2))/2];   
    }

    /// matrix element (i,j), (i,j) must be a valid position
    T & operator() (int i, int j) 
    {
      if (i < bw)
        return mem[n + (i * (i-1)) / 2 + j];
      else
        return mem[n + i * (bw-2) + j - ((bw-1)*(bw-2))/2];   
    }

    /// matrix size
    int Size() const { return n; }
    /// band-width of triangular matrix
    int BandWidth() const { return bw; }
    /// computes required memory
    static int RequiredMem (int n, int bw)
    { return n*bw - (bw * (bw-1)) / 2 + n; }
  };
  

  ///  output operator.
  template<typename T>
  inline std::ostream & operator<< (std::ostream & s, const FlatBandCholeskyFactors<T> & m)
  {
    m.Print (s);
    return s;
  }





  /**
     Cholesky factors of a band matrix.
  */
  template <class T = double>
  class BandCholeskyFactors : public FlatBandCholeskyFactors<T>
  {
  public:
    /// allocate memory and factor the matrix a
    BandCholeskyFactors (const SymBandMatrix<T> & a)
      : FlatBandCholeskyFactors<T> (a.Height(), 
                                    a.BandWidth(),
                                    new T[FlatBandCholeskyFactors<T>::RequiredMem (a.Height(), a.BandWidth())])
    {
      this->Factor (a);
    }

    /// delete memory
    ~BandCholeskyFactors ()
    {
      delete [] this->mem;
    }
  };

}

#endif
