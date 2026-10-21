#ifndef FILE_NGS_SPARSEMATRIX
#define FILE_NGS_SPARSEMATRIX

/**************************************************************************/
/* File:   sparsematrix.hpp                                               */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Oct. 94, 15 Jan. 02                                        */
/**************************************************************************/


#include "vvector.hpp"
#include "basematrix.hpp"

namespace ngla
{

  // the vector types are not template arguments, their scalar type is dispatched at run time
  template<class TM>
  class SparseMatrix;

  // the former storage base class is SparseMatrix itself
  template<class TM>
  using SparseMatrixTM = SparseMatrix<TM>;

  template<class TM>
  class SparseMatrixSymmetric;



  class BaseJacobiPrecond;

  template<class TM>
  class JacobiPrecond;

  template<class TM>
  class JacobiPrecondSymmetric;


  class BaseBlockJacobiPrecond;

  template<class TM>
  class BlockJacobiPrecond;

  template<class TM>
  class BlockJacobiPrecondSymmetric;





  


  /** 
      The graph of a sparse matrix.
  */
  class NGS_DLL_HEADER MatrixGraph : public BaseMatrix
  {
  public:
    typedef int ColIdx;
  protected:
    /// number of rows
    size_t size;
    /// with of matrix
    size_t width;
    /// non-zero elements
    size_t nze; 

    /// column numbers
    // Array<int, size_t> colnr;
    NumaDistributedArray<ColIdx> colnr;

    /// pointer to first in row
    Array<size_t> firsti;
  
    /// row has same non-zero elements as previous row
    Array<size_t> same_nze;
    
    /// balancing for multi-threading
    Partitioning balance;

    /// owner of arrays ?
    bool owner;

  public:
    /// arbitrary number of els/row
    MatrixGraph (FlatArray<int> elsperrow, size_t awidth);
    /// matrix of height as, uniform number of els/row
    MatrixGraph (size_t as, int max_elsperrow);    
    /// shadow matrix graph
    MatrixGraph (const MatrixGraph & graph); // , bool stealgraph = false);
    /// move-constuctor
    MatrixGraph (MatrixGraph && graph);
    /// 
    MatrixGraph (size_t size, size_t width,
                 FlatTable<int> rowelements, FlatTable<int> colelements, bool symmetric);
    /// 
    // MatrixGraph (const Table<int> & dof2dof, bool symmetric);
    virtual ~MatrixGraph ();

    /// eliminate unused columne indices (was never implemented)
    void Compress();
  
    /// returns position of Element (i, j), exception for unused
    size_t GetPosition (size_t i, size_t j) const;
    
    /// returns position of Element (i, j), -1 for unused
    size_t GetPositionTest (size_t i, size_t j) const;

    /// find positions of n sorted elements, overwrite pos, exception for unused
    void GetPositionsSorted (size_t row, size_t n, int * pos) const;

    /// returns position of new element
    size_t CreatePosition (size_t i, size_t j);

    size_t Size() const { return size; }

    size_t NZE() const override { return nze; }

    // full col-index array
    FlatArray<ColIdx> GetColIndices() const { return colnr; }
    
    // col-indices of the i-th row
    FlatArray<ColIdx> GetRowIndices(size_t i) const
    { return FlatArray<ColIdx> (firsti[i+1]-firsti[i], colnr+firsti[i]); }      

    size_t First (int i) const { return firsti[i]; }
    FlatArray<size_t> GetFirstArray () const  { return firsti; } 

    void FindSameNZE();
    void CalcBalancing ();
    const Partitioning & GetBalancing() const { return balance; } 

    void EmbedHeight (size_t starti, size_t newheight);
    void EmbedWidth (size_t starti, size_t newwidth);
    
    ostream & Print (ostream & ost) const override;

    virtual Array<MemoryUsage> GetMemoryUsage () const override;    

    // the graph knows the shape only
    VecFormat RowFormat () const override { return VecFormat (width); }
    VecFormat ColFormat () const override { return VecFormat (size); }
  };


  





  /// base class for all sparse matrices
  class NGS_DLL_HEADER BaseSparseMatrix : public MatrixGraph
  {
  protected:
    /// sparse direct solver
    Flags inverseflags;
    bool spd = false;
    /// scalar of the vectors the matrix creates, unset: the matrix scalar
    optional<Scalar> vecscal;
    
  public:
    using MatrixGraph::MatrixGraph;
    using MatrixGraph::ColIdx;
    /*
    BaseSparseMatrix (int as, int max_elsperrow)
      : MatrixGraph (as, max_elsperrow)  
    { ; }
    
    BaseSparseMatrix (const Array<int> & elsperrow, int awidth)
      : MatrixGraph (elsperrow, awidth) 
    { ; }

    BaseSparseMatrix (int size, int width, const Table<int> & rowelements, 
		      const Table<int> & colelements, bool symmetric)
      : MatrixGraph (size, width, rowelements, colelements, symmetric)
    { ; }
    */

    BaseSparseMatrix (const MatrixGraph & agraph)
      : MatrixGraph (agraph)
    { ; }

    BaseSparseMatrix (MatrixGraph && agraph)
      : MatrixGraph (std::move(agraph))
    { ; }
    
    BaseSparseMatrix (const BaseSparseMatrix & amat)
      : MatrixGraph (amat), vecscal(amat.vecscal)
    { ; }   

    BaseSparseMatrix (BaseSparseMatrix && amat)
      : MatrixGraph (std::move(amat)), vecscal(amat.vecscal)
    { ; }

    virtual ~BaseSparseMatrix ();

    BaseSparseMatrix & operator= (double s)
    {
      AsVector() = s;
      return *this;
    }

    BaseSparseMatrix & Add (double s, const BaseSparseMatrix & m2)
    {
      AsVector() += s * m2.AsVector();
      return *this;
    }

    virtual shared_ptr<BaseSparseMatrix> CreateSparseMatrix() const override 
    {
      return dynamic_pointer_cast<BaseSparseMatrix>(const_cast<BaseSparseMatrix*>(this)->shared_from_this());
    }

    virtual shared_ptr<BaseJacobiPrecond> CreateJacobiPrecond (shared_ptr<BitArray> inner = nullptr) const 
    {
      throw Exception ("BaseSparseMatrix::CreateJacobiPrecond");
    }

    virtual shared_ptr<BaseBlockJacobiPrecond>
      CreateBlockJacobiPrecond (shared_ptr<Table<int>> blocks,
                                const BaseVector * constraint = 0,
                                bool parallel  = 1,
                                shared_ptr<BitArray> freedofs = NULL) const
    { 
      throw Exception ("BaseSparseMatrix::CreateBlockJacobiPrecond");
    }

    virtual shared_ptr<BaseSparseMatrix> CreateTranspose(bool sorted = true) const
    {
      throw Exception ("BaseSparseMatrix::CreateTranspose");      
    }
      
    virtual shared_ptr<BaseMatrix>
      InverseMatrix (shared_ptr<BitArray> subset = nullptr) const override
    { 
      throw Exception ("BaseSparseMatrix::CreateInverse called");
    }

    virtual shared_ptr<BaseMatrix>
      InverseMatrix (shared_ptr<const Array<int>> clusters) const override
    { 
      throw Exception ("BaseSparseMatrix::CreateInverse called");
    }

    virtual shared_ptr<BaseSparseMatrix> Restrict (const SparseMatrix<double> & prol,
                                                   shared_ptr<BaseSparseMatrix> cmat = nullptr ) const
    {
      throw Exception ("BaseSparseMatrix::Restrict");
    }

    virtual shared_ptr<BaseSparseMatrix> Reorder (const Array<size_t> & reorder) const
    {
      throw Exception ("BaseSparseMatrix::Reorder");
    }
    
    virtual void SetInverseFlags (const Flags & flags) override
    { inverseflags = flags; }
    virtual const Flags & GetInverseFlags () const { return inverseflags; } 
    
    // a real matrix acting on complex vectors creates complex vectors
    void SetVectorScalar (Scalar ascal) { vecscal = ascal; }
    void CopyVectorScalar (const BaseSparseMatrix & other) { vecscal = other.vecscal; }

    void SetSPD (bool aspd = true) { spd = aspd; }
    bool IsSPD () const { return spd; }
    virtual size_t NZE () const override { return nze; }
    virtual tuple<int,int> EntrySizes() const = 0;

    virtual shared_ptr<BaseMatrix> DeleteZeroElements(double tol) const override
    {
      throw Exception ("DeleteZeroElements not overloaded");
    }
  };



  /*
    Vector types of the sparse matrix family: the entry shape is given by the
    matrix entry, the scalar type is taken from the vectors at run time.
  */

  // scalar types allowed for y += s * A x, with matrix scalar TSCAL
  template <typename TSCAL, typename TS, typename TX, typename TY>
  constexpr bool ValidVectorScalars ()
  {
    constexpr bool valid = requires (TY & ey, TS es, TSCAL em, TX ex) { ey += es * (em * ex); };
    if constexpr (is_same_v<TSCAL,float> || is_same_v<TX,float> || is_same_v<TY,float>)
      return valid && is_same_v<TSCAL,float> && is_same_v<TX,float> && is_same_v<TY,float>;
    else
      return valid;
  }

  // the entry size of the matrix, or a flat vector of the same total size
  NGS_DLL_HEADER void CheckEntrySize (const BaseVector & v, int es, size_t n, const char * name);
  NGS_DLL_HEADER void ThrowIllegalScalars (const char * name, Scalar mat, const BaseVector & x, const BaseVector & y);
  NGS_DLL_HEADER void CheckSameScalar (std::initializer_list<const BaseVector*> vecs, const char * name);

  // the factor in the scalar type T of the result, a complex factor needs complex vectors
  template <typename T>
  T ScalarAs (Scalar s, const char * name)
  {
    return std::visit ([&] (auto v) -> T
    {
      if constexpr (is_same_v<decltype(v),Complex> && !is_same_v<T,Complex>)
        throw Exception (string(name) + ": complex factor for real vectors");
      else
        return T(v);
    }, s);
  }

  /*
    calls f.template operator()<TVX,TVY>(sy) with the run-time vector types of x and y,
    and the factor s in the scalar type of y
  */
  template <typename TSCAL, typename TVX0, typename TVY0, typename FUNC>
  void DispatchVectorTypes (Scalar s, const BaseVector & x, size_t nx, const BaseVector & y, size_t ny,
                            const char * name, FUNC && f)
  {
    std::visit ([&] (auto sx, auto sy)
    {
      typedef decltype(sx) TX;
      typedef decltype(sy) TY;
      if constexpr (ValidVectorScalars<TSCAL,TY,TX,TY>())
        {
          typedef WithScalar<TVX0,TX> TVX;
          typedef WithScalar<TVY0,TY> TVY;
          CheckEntrySize (x, sizeof(TVX)/sizeof(TX), nx, name);
          CheckEntrySize (y, sizeof(TVY)/sizeof(TY), ny, name);
          f.template operator()<TVX,TVY> (ScalarAs<TY> (s, name));
        }
      else
        ThrowIllegalScalars (name, TSCAL(0), x, y);
    }, x.GetScalarType(), y.GetScalarType());
  }

  // calls f.template operator()<TV>() for several vectors of the same run-time type
  template <typename TSCAL, typename TV0, typename FUNC>
  void DispatchVectorType (std::initializer_list<const BaseVector*> vecs, size_t n,
                           const char * name, FUNC && f)
  {
    CheckSameScalar (vecs, name);
    const BaseVector & v0 = **vecs.begin();
    std::visit ([&] (auto sv)
    {
      typedef decltype(sv) TS;
      if constexpr (ValidVectorScalars<TSCAL,double,TS,TS>())
        {
          typedef WithScalar<TV0,TS> TV;
          for (auto v : vecs)
            CheckEntrySize (*v, sizeof(TV)/sizeof(TS), n, name);
          f.template operator()<TV> ();
        }
      else
        ThrowIllegalScalars (name, TSCAL(0), v0, v0);
    }, v0.GetScalarType());
  }
  
  template <typename TSCAL>
  class NGS_DLL_HEADER S_BaseSparseMatrix : public BaseSparseMatrix
  {
  protected:
    int entry_height, entry_width;
    // in general entry_size=entry_height*entry_width, but matrix entry type could
    // also be a diagonal matrix
    int entry_size;
    VFlatVector<TSCAL> asvec;
    using MatrixGraph::ColIdx;    
  public:
    using BaseSparseMatrix::BaseSparseMatrix;
    /*
    void SetEntrySize (int eh, int ew, int es)
    {
      entry_height = eh;
      entry_width = ew;
      entry_size = es;
    }
    */
    
    virtual BaseVector & AsVector() override
    {
      return asvec; 
    }

    virtual const BaseVector & AsVector() const override
    { 
      return asvec; 
    }

    // tuple<int,int> EntryShape() const { return { entry_height, entry_width }; }
    
    FlatVector<TSCAL> GetRowValue (int row, int j)
    {
      TSCAL * p = asvec(entry_size * (firsti[row] + j)).Data(0);
      return FlatVector<TSCAL> (entry_size, p);
    }
    
    FlatMatrix<TSCAL> GetRowValueMat (int row, int j) const
    {
      // TSCAL * p = asvec(entry_size * (firsti[row] + j)).Data(0);
      TSCAL * p = &asvec(entry_size * (firsti[row] + j));
      return FlatMatrix<TSCAL> (entry_height, entry_width, p);
    }
  };

  
  /// A general, sparse matrix
  template<class TM>
  class NGS_DLL_HEADER SparseMatrix : public S_BaseSparseMatrix<typename mat_traits<TM>::TSCAL>
  {
  protected:
    typedef S_BaseSparseMatrix<typename mat_traits<TM>::TSCAL> BASE;
    using BASE::owner;
    using BASE::size;
    using BASE::width;
    using BASE::nze;
    using BASE::asvec;
  public:
    using BASE::firsti;
    using BASE::colnr;
    using BASE::balance;
    using BASE::CreatePosition;
    using BASE::GetPositionTest;
    using BASE::FindSameNZE;
    using BASE::AsVector;
    using typename BASE::ColIdx;

    typedef TM TENTRY;
    typedef typename mat_traits<TM>::TSCAL TSCAL;
    // entry shapes of the vectors, the vector scalar is chosen at run time
    typedef typename mat_traits<TM>::TV_ROW TV_ROW;
    typedef typename mat_traits<TM>::TV_COL TV_COL;

    // Array<TM, size_t> data;
    NumaDistributedArray<TM> data;
  protected:
    TM nul;
    bool hermitian = false;
    
  public:
    void SetEntrySize()
    {
      BASE::entry_height = ngbla::Height<TM>();
      BASE::entry_width = ngbla::Width<TM>();
      BASE::entry_size = sizeof(TM)/sizeof(TSCAL);
      BASE::is_complex = ngbla::IsComplex<TM>();
      asvec.AssignMemory (nze*sizeof(TM)/sizeof(TSCAL), (void*)data.Addr(0));
    }

    ///
    SparseMatrix (int as, int max_elsperrow)
      : BASE (as, max_elsperrow), data(nze), nul(TSCAL(0))
    { SetEntrySize (); }

    SparseMatrix (const Array<int> & elsperrow)
      : SparseMatrix (elsperrow, elsperrow.Size()) { ; }

    SparseMatrix (const Array<int> & elsperrow, int awidth)
      : BASE (elsperrow, awidth), data(nze), nul(TSCAL(0))
    { SetEntrySize (); }

    SparseMatrix (int height, int width, FlatTable<int> rowelements, 
		  FlatTable<int> colelements, bool symmetric)
      : BASE (height, width, rowelements, colelements, symmetric), data(nze), nul(TSCAL(0))
    { SetEntrySize (); }

    SparseMatrix (const MatrixGraph & agraph)
      : BASE (agraph), data(nze), nul(TSCAL(0))
    {
      SetEntrySize ();
      FindSameNZE();
    }

    SparseMatrix (MatrixGraph && agraph)
      : BASE (std::move(agraph)), data(nze), nul(TSCAL(0))
    {
      SetEntrySize ();
      FindSameNZE();
    }

    SparseMatrix (const SparseMatrix & amat)
      : BASE (amat), data(nze), nul(TSCAL(0))
    {
      SetEntrySize ();
      AsVector() = amat.AsVector();
    }

    SparseMatrix (SparseMatrix && amat)
      : BASE (std::move(amat)), nul(TSCAL(0))
    {
      data.Swap(amat.data);
      SetEntrySize ();
    }

    static shared_ptr<SparseMatrix> CreateFromCOO (FlatArray<int> i, FlatArray<int> j,
                                                   FlatArray<TM> val, size_t h, size_t w);
      
    virtual ~SparseMatrix ();

    size_t Height() const { return size; }
    size_t Width() const { return width; }
    virtual int VHeight() const override { return size; }
    virtual int VWidth() const override { return width; }
    virtual tuple<size_t, size_t> Shape() const override { return { size, width }; }
    
    // the scalar of the vectors the matrix creates
    Scalar GetVectorScalar () const { return this->vecscal ? *this->vecscal : Scalar(TSCAL(0)); }

    bool IsHermitian() const { return hermitian; }
    void SetHermitian(bool herm) { hermitian = herm; }
    
    TM & operator[] (size_t i)  { return data[i]; }
    const TM & operator[] (size_t i) const { return data[i]; }

    TM & operator() (size_t row, size_t col)
    {
      return data[CreatePosition(row, col)];
    }

    const TM & operator() (size_t row, size_t col) const
    {
      size_t pos = GetPositionTest (row,col);
      if (pos != numeric_limits<size_t>::max())
	return data[pos];
      else
	return nul;
    }

    void PrefetchRow (size_t rownr) const;

    // full value array
    FlatVector<TM> GetValues() const { return FlatVector<TM> (data.Size(), data.Addr(0)); }

    FlatVector<TM> GetRowValues(int i) const
    { return FlatVector<TM> (firsti[i+1]-firsti[i], data+firsti[i]); }

    static bool IsRegularIndex (int index) { return index >= 0; }
    virtual void AddElementMatrix(FlatArray<int> dnums1, 
                                  FlatArray<int> dnums2, 
                                  BareSliceMatrix<TSCAL> elmat,
                                  bool use_atomic = false);

    virtual void AddElementMatrixSymmetric(FlatArray<int> dnums,
                                           BareSliceMatrix<TSCAL> elmat,
                                           bool use_atomic = false);

    virtual void SetZero() override;

    ///
    virtual ostream & Print (ostream & ost) const override;

    ///
    virtual Array<MemoryUsage> GetMemoryUsage () const override;    

    virtual tuple<int,int> EntrySizes() const override { return { ngbla::Height<TM>(), ngbla::Width<TM>() }; }

    virtual shared_ptr<BaseMatrix> CreateMatrix () const override;
    // a DeviceSparseMatrix for scalar TM if a gpu backend is registered
    virtual shared_ptr<BaseMatrix> CreateDeviceMatrix () const override;
    ///
    virtual AutoVector CreateVector () const override;
    VecFormat RowFormat () const override;
    VecFormat ColFormat () const override;


    BaseMatrix::OperatorInfo GetOperatorInfo () const override
    { return { string("SparseMatrix")+typeid(TM).name()+" (nze="+ToString(this->NZE())+")", this->Height(), this->Width() }; }
    
    virtual shared_ptr<BaseJacobiPrecond>
      CreateJacobiPrecond (shared_ptr<BitArray> inner) const override;
    
    virtual shared_ptr<BaseBlockJacobiPrecond>
      CreateBlockJacobiPrecond (shared_ptr<Table<int>> blocks,
                                const BaseVector * constraint = 0,
                                bool parallel = 1,
                                shared_ptr<BitArray> freedofs = NULL) const override;

    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = nullptr) const override;
    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<const Array<int>> clusters) const override;

    virtual shared_ptr<BaseSparseMatrix> Restrict (const SparseMatrix<double> & prol,
					 shared_ptr<BaseSparseMatrix> cmat = nullptr) const override;
    
    virtual shared_ptr<BaseSparseMatrix> Reorder (const Array<size_t> & reorder) const override;
    
    virtual shared_ptr<BaseSparseMatrix> CreateTranspose(bool sorted = true) const override;

    virtual shared_ptr<BaseMatrix> DeleteZeroElements(double tol) const override;
    
    /// row times vector, the scalar of the result is the product of matrix and vector scalars
    template <typename TVX>
    INLINE auto RowTimesVector (int row, const FlatVector<TVX> vec) const
    {
      typedef decltype(TSCAL()*typename mat_traits<TVX>::TSCAL()) TPROD;
      WithScalar<TV_COL,TPROD> sum = TPROD(0);
      for (size_t j = firsti[row]; j < firsti[row+1]; j++)
	sum += data[j] * vec(colnr[j]);
      return sum;
    }

    ///
    template <typename TEL, typename TVX>
    void AddRowTransToVector (int row, const TEL & el, FlatVector<TVX> vec) const
    {
      WithScalar<TV_COL,typename mat_traits<TVX>::TSCAL> hel = el;
      size_t first = firsti[row];
      size_t last = firsti[row+1];

      const ColIdx * colpi = colnr.Addr(0);
      const TM * datap = data.Addr(0);

      for (size_t j = first; j < last; j++)
        vec[colpi[j]] += Trans(datap[j]) * hel; 
    }

    template <typename TEL, typename TVX>
    void AddRowTransToVectorAtomic (int row, const TEL & el, FlatVector<TVX> vec) const
    {
      WithScalar<TV_COL,typename mat_traits<TVX>::TSCAL> hel = el;
      size_t first = firsti[row];
      size_t last = firsti[row+1];

      const ColIdx * colpi = colnr.Addr(0);
      const TM * datap = data.Addr(0);

      for (size_t j = first; j < last; j++)
        // vec[colpi[j]] += Trans(datap[j]) * el;
        AtomicAdd (vec[colpi[j]], Trans(datap[j]) * hel);
    }

    
    ///
    template <typename TEL, typename TVX>
    void AddRowConjTransToVector (int row, const TEL & el, FlatVector<TVX> vec) const
    {
      WithScalar<TV_COL,typename mat_traits<TVX>::TSCAL> hel = el;
      size_t first = firsti[row];
      size_t last = firsti[row+1];

      const ColIdx * colpi = colnr.Addr(0);
      const TM * datap = data.Addr(0);

      for (size_t j = first; j < last; j++)
        vec[colpi[j]] += Conj(Trans(datap[j])) * hel; 
    }


    using BaseMatrix::MultAdd;
    using BaseMatrix::MultTransAdd;
    virtual void MultAdd (Scalar s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (Scalar s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultConjTransAdd (Complex s, const BaseVector & x, BaseVector & y) const override;

    virtual void MultAdd (FlatVector<double> alpha, const MultiVector & x, MultiVector & y) const override;

    
    virtual void MultAdd1 (double s, const BaseVector & x, BaseVector & y,
			   const BitArray * ainner = NULL,
			   const Array<int> * acluster = NULL) const override;
    
    virtual void DoArchive (Archive & ar) override;

  protected:
    void T_MultTransAdd (Scalar s, const BaseVector & x, BaseVector & y, bool conjugate) const;
  };

  

  /// A symmetric sparse matrix
  template<class TM>
  class NGS_DLL_HEADER SparseMatrixSymmetric : public SparseMatrix<TM>
  {

  public:
    using SparseMatrix<TM>::firsti;
    using SparseMatrix<TM>::colnr;
    using SparseMatrix<TM>::data;
    using BaseSparseMatrix::ColIdx;    

    typedef typename mat_traits<TM>::TSCAL TSCAL;
    typedef typename mat_traits<TM>::TV_ROW TV;
    typedef TV TV_COL;
    typedef TV TV_ROW;

    SparseMatrixSymmetric (int as, int max_elsperrow)
      : SparseMatrix<TM> (as, max_elsperrow)
    { ; }
  
    SparseMatrixSymmetric (const Array<int> & elsperrow)
      : SparseMatrix<TM> (elsperrow, elsperrow.Size())
    { ; }

    SparseMatrixSymmetric (int size, FlatTable<int> rowelements)
      : SparseMatrix<TM> (size, size, rowelements, rowelements, true)
    { ; }

    SparseMatrixSymmetric (const MatrixGraph & agraph);
    SparseMatrixSymmetric (MatrixGraph && agraph);    

    SparseMatrixSymmetric (const SparseMatrixSymmetric & amat)
      : SparseMatrix<TM> (amat)
    { 
      this->AsVector() = amat.AsVector(); 
    }

    SparseMatrixSymmetric (const SparseMatrix<TM> & amat)
      : SparseMatrix<TM> (amat)
      { 
        this->AsVector() = amat.AsVector(); 
      }
    
    ///
    virtual ~SparseMatrixSymmetric ();
    virtual shared_ptr<BaseMatrix> CreateDeviceMatrix () const override;

    SparseMatrixSymmetric & operator= (double s)
    {
      this->AsVector() = s;
      return *this;
    }

    virtual shared_ptr<BaseMatrix> CreateMatrix () const override
    {
      return make_shared<SparseMatrixSymmetric> (*this);
    }

    virtual shared_ptr<BaseSparseMatrix> CreateTranspose (bool sorted = true) const override
    {
      return make_shared<SparseMatrixSymmetric> (*this);
    }

    virtual void AddElementMatrix(FlatArray<int> dnums1, 
				  FlatArray<int> dnums2, 
				  BareSliceMatrix<TSCAL> elmat,
                                  bool use_atomic = false) override
    {
      this->AddElementMatrixSymmetric (dnums1, elmat, use_atomic);
    }
    
    virtual shared_ptr<BaseJacobiPrecond> CreateJacobiPrecond (shared_ptr<BitArray> inner) const override;
    virtual shared_ptr<BaseBlockJacobiPrecond>
      CreateBlockJacobiPrecond (shared_ptr<Table<int>> blocks,
                                const BaseVector * constraint = 0,
                                bool parallel  = 1,
                                shared_ptr<BitArray> freedofs = NULL) const override;

    virtual shared_ptr<BaseSparseMatrix> Restrict (const SparseMatrix<double> & prol,
					 shared_ptr<BaseSparseMatrix> cmat = nullptr) const override;

    ///
    using SparseMatrix<TM>::MultAdd;
    using SparseMatrix<TM>::MultTransAdd;
    virtual void MultAdd (Scalar s, const BaseVector & x, BaseVector & y) const override;

    virtual void MultTransAdd (Scalar s, const BaseVector & x, BaseVector & y) const override
    {
      MultAdd (s, x, y);
    }
    // via MultTransAdd with conjugated vectors, the stored triangle alone is not enough
    virtual void MultConjTransAdd (Complex s, const BaseVector & x, BaseVector & y) const override
    {
      BaseMatrix::MultConjTransAdd (s, x, y);
    }


    /*
      y += s L * x
    */
    virtual void MultAdd1 (double s, const BaseVector & x, BaseVector & y,
			   const BitArray * ainner = NULL,
			   const Array<int> * acluster = NULL) const override;


    /*
      y += s (D + L^T) * x
    */
    virtual void MultAdd2 (double s, const BaseVector & x, BaseVector & y,
			   const BitArray * ainner = NULL,
			   const Array<int> * acluster = NULL) const override;
    




    using SparseMatrix<TM>::RowTimesVector;
    using SparseMatrix<TM>::AddRowTransToVector;


    template <typename TVX>
    auto RowTimesVectorNoDiag (int row, const FlatVector<TVX> vec) const
    {
      typedef decltype(TSCAL()*typename mat_traits<TVX>::TSCAL()) TPROD;
      WithScalar<TV,TPROD> sum = TPROD(0);

      size_t last = firsti[row+1];
      size_t first = firsti[row];
      if (last == first) return sum;
      if (colnr[last-1] == row) last--;

      for (size_t j = first; j < last; j++)
	sum += data[j] * vec(colnr[j]);
      return sum;
    }

    template <typename TEL, typename TVX>
    void AddRowTransToVectorNoDiag (int row, const TEL & el, FlatVector<TVX> vec) const
    {
      WithScalar<TV,typename mat_traits<TVX>::TSCAL> hel = el;
      size_t first = firsti[row];
      size_t last = firsti[row+1];

      if (first == last) return;
      if (this->colnr[last-1] == row) last--;

      for (size_t j = first; j < last; j++)
        vec[colnr[j]] += Trans(data[j]) * hel;
    }
  
    BaseSparseMatrix & AddMerge (double s, const SparseMatrixSymmetric  & m2);

    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = nullptr) const override;
    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<const Array<int>> clusters) const override;
  };

  // specialized in sparsematrix.cpp, declared before anything instantiates them
  template <> shared_ptr<BaseSparseMatrix>
  SparseMatrix<double> :: Restrict (const SparseMatrix<double> & prol, shared_ptr<BaseSparseMatrix> cmat) const;
  template <> shared_ptr<BaseSparseMatrix>
  SparseMatrix<Complex> :: Restrict (const SparseMatrix<double> & prol, shared_ptr<BaseSparseMatrix> cmat) const;
  template <> shared_ptr<BaseSparseMatrix>
  SparseMatrixSymmetric<double> :: Restrict (const SparseMatrix<double> & prol, shared_ptr<BaseSparseMatrix> cmat) const;
  template <> void
  SparseMatrix<double> :: MultAdd (FlatVector<double> alpha, const MultiVector & x, MultiVector & y) const;

  [[deprecated("Use sparsematrix->CreateTranspose() instead!")]]            
  NGS_DLL_HEADER shared_ptr<SparseMatrix<double>> TransposeMatrix (const SparseMatrix<double> & mat);

  NGS_DLL_HEADER shared_ptr<SparseMatrix<double>>
  MatAdd (double sa, const SparseMatrix<double> & mata,
          double sb, const SparseMatrix<double> & matb);

  NGS_DLL_HEADER shared_ptr<SparseMatrix<Complex>>
  MatAdd (Complex sa, const SparseMatrix<Complex> & mata,
          Complex sb, const SparseMatrix<Complex> & matb);

  
  NGS_DLL_HEADER shared_ptr<SparseMatrix<double>>
  MatMult (const SparseMatrix<double> & mata, const SparseMatrix<double> & matb, bool sort_output = true);
  NGS_DLL_HEADER shared_ptr<SparseMatrix<float>>
  MatMult (const SparseMatrix<float> & mata, const SparseMatrix<float> & matb, bool sort_output = true);
  NGS_DLL_HEADER shared_ptr<SparseMatrix<Complex>>
  MatMult (const SparseMatrix<Complex> & mata, const SparseMatrix<Complex> & matb, bool sort_output = true);

#ifdef GOLD
#include <sparsematrix_spec.hpp>
#endif

shared_ptr<BaseMatrix> CreateSparseMatrixInverse(shared_ptr<const BaseSparseMatrix> A,
                                                 shared_ptr<BitArray> subset,
                                                 shared_ptr<const Array<int>> clusters);

#ifdef FILE_SPARSEMATRIX_CPP
#define SPARSEMATRIX_EXTERN
#else
#define SPARSEMATRIX_EXTERN extern
  

  SPARSEMATRIX_EXTERN template class SparseMatrix<double>;
  SPARSEMATRIX_EXTERN template class SparseMatrix<Complex>;
  
  SPARSEMATRIX_EXTERN template class SparseMatrixSymmetric<double>;
  SPARSEMATRIX_EXTERN template class SparseMatrixSymmetric<Complex>;

#define INST_SPMS(N) \
  SPARSEMATRIX_EXTERN template class SparseMatrix<Mat<N, N, double>>; \
  SPARSEMATRIX_EXTERN template class SparseMatrix<Mat<1, N, double>>; \
  SPARSEMATRIX_EXTERN template class SparseMatrix<Mat<N, 1, double>>; \
  SPARSEMATRIX_EXTERN template class SparseMatrix<Mat<N, N, Complex>>; \
  SPARSEMATRIX_EXTERN template class SparseMatrix<Mat<1, N, Complex>>; \
  SPARSEMATRIX_EXTERN template class SparseMatrix<Mat<N, 1, Complex>>; \
  SPARSEMATRIX_EXTERN template class SparseMatrixSymmetric<Mat<N, N, double>>; \
  SPARSEMATRIX_EXTERN template class SparseMatrixSymmetric<Mat<N, N, Complex>>;


  NGS_FOR_SYS_DIMS(INST_SPMS)

#undef INST_SPMS
#undef SPARSEMATRIX_EXTERN


#endif



  template <typename TSCAL>
  class NGS_DLL_HEADER SparseBlockMatrix : public S_BaseSparseMatrix<TSCAL>
  {
    size_t bheight, bwidth;
    NumaDistributedArray<TSCAL> data;
    
    typedef S_BaseSparseMatrix<TSCAL> BASE;
    using BASE::firsti;
    using BASE::colnr;
    using BASE::owner;
    using BASE::size;
    using BASE::width;
    using BASE::nze;
    using BASE::balance;
    using BASE::asvec;
    
  public:
    using BASE::CreatePosition;
    using BASE::GetPositionTest;
    using BASE::FindSameNZE;
    // using BASE::SetEntrySize;
    using BASE::AsVector;
    
    SparseBlockMatrix (const MatrixGraph & agraph, size_t abheight, size_t abwidth)
      : BASE (agraph), bheight(abheight), bwidth(abwidth),
      data(nze*bheight*bwidth)
        { 
          // SetEntrySize (bheight, bwidth, bheight*bwidth);
          BASE::entry_height = bheight;
          BASE::entry_width = bwidth;
          BASE::entry_size = bheight*bwidth;
          BASE::is_complex = ngbla::IsComplex<TSCAL>();
          
          asvec.AssignMemory (nze*bheight*bwidth, (void*)data.Addr(0));
          // FindSameNZE();
        }

    SparseBlockMatrix (MatrixGraph && agraph, size_t abheight, size_t abwidth)
      : BASE (std::move(agraph)), bheight(abheight), bwidth(abwidth),
      data(nze*bheight*bwidth)
        { 
          // SetEntrySize (bheight, bwidth, bheight*bwidth);
          BASE::entry_height = bheight;
          BASE::entry_width = bwidth;
          BASE::entry_size = bheight*bwidth;
          BASE::is_complex = ngbla::IsComplex<TSCAL>();
          
          asvec.AssignMemory (nze*bheight*bwidth, (void*)data.Addr(0));
          // FindSameNZE();
        }
    

    
    tuple<int,int> EntrySizes() const override { return { bheight, bwidth }; }
    
    
    VecFormat RowFormat () const override { return VecFormat (this->width, TSCAL(0), int(bwidth)); }
    VecFormat ColFormat () const override { return VecFormat (this->size, TSCAL(0), int(bheight)); }

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;

    
    void SetZero() override
    {
      data = TSCAL(0);
    }
    
    virtual void AddElementMatrix(FlatArray<int> dnums1, 
                                  FlatArray<int> dnums2, 
                                  BareSliceMatrix<TSCAL> elmat,
                                  bool use_atomic = false);

    ostream & Print (ostream & ost) const override;
  };


  

}

#endif
