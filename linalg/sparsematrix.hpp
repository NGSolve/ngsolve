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

  template<class TM>
  class SparseMatrixTM ;

  template<class TM, 
	   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class SparseMatrix;


  template<class TM, 
	   class TV = typename mat_traits<TM>::TV_ROW>
  class SparseMatrixSymmetric;



  class BaseJacobiPrecond;

  template<class TM, 
	   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class JacobiPrecond;

  template<class TM, 
	   class TV = typename mat_traits<TM>::TV_ROW>
  class JacobiPrecondSymmetric;


  class BaseBlockJacobiPrecond;

  template<class TM, 
	   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class BlockJacobiPrecond;

  template<class TM, 
	   class TV = typename mat_traits<TM>::TV_ROW>
  class BlockJacobiPrecondSymmetric;





  
#ifdef USE_PARDISO
  const INVERSETYPE default_inversetype = PARDISO;
#else
#ifdef USE_MUMPS
  const INVERSETYPE default_inversetype = MUMPS;
#else
#ifdef USE_UMFPACK
  const INVERSETYPE default_inversetype = UMFPACK;
#else
  const INVERSETYPE default_inversetype = SPARSECHOLESKY;
#endif
#endif
#endif


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

    size_t NZE() const { return nze; }

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
    
    ostream & Print (ostream & ost) const;

    virtual Array<MemoryUsage> GetMemoryUsage () const;    

    const MemoryTracer & GetMemoryTracer() const
    {
      return mem_tracer;
    }

    virtual AutoVector CreateRowVector () const
    { throw Exception("MatrixGraph::CreateRowVector called"); }
    
    virtual AutoVector CreateColVector () const
    { throw Exception("MatrixGraph::CreateRowVector called"); }

  private:
    
    MemoryTracer mem_tracer = {"MatrixGraph",
      colnr, "colnr",
      firsti, "firsti",
      same_nze, "same_nze"
    };
  };


  





  /// base class for all sparse matrices
  class NGS_DLL_HEADER BaseSparseMatrix : public MatrixGraph
  {
  protected:
    /// sparse direct solver
    mutable INVERSETYPE inversetype = default_inversetype;   
    Flags inverseflags;
    bool spd = false;
    
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
    
    BaseSparseMatrix (const BaseSparseMatrix & amat)
      : MatrixGraph (amat)
    { ; }   

    BaseSparseMatrix (BaseSparseMatrix && amat)
      : MatrixGraph (std::move(amat))
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

    virtual shared_ptr<BaseSparseMatrix> CreateTranspose() const
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

    virtual shared_ptr<BaseSparseMatrix> Restrict (const SparseMatrixTM<double> & prol,
                                                   shared_ptr<BaseSparseMatrix> cmat = nullptr ) const
    {
      throw Exception ("BaseSparseMatrix::Restrict");
    }

    virtual shared_ptr<BaseSparseMatrix> Reorder (const Array<size_t> & reorder) const
    {
      throw Exception ("BaseSparseMatrix::Reorder");
    }
    
    virtual INVERSETYPE SetInverseType ( INVERSETYPE ainversetype ) const override
    {

      INVERSETYPE old_invtype = inversetype;
      inversetype = ainversetype; 
      return old_invtype;
    }

    virtual INVERSETYPE SetInverseType ( string ainversetype ) const override;

    virtual INVERSETYPE  GetInverseType () const override
    { return inversetype; }
    virtual void SetInverseFlags (const Flags & flags) override
    { inverseflags = flags; }
    virtual const Flags & GetInverseFlags () const { return inverseflags; } 
    
    void SetSPD (bool aspd = true) { spd = aspd; }
    bool IsSPD () const { return spd; }
    virtual size_t NZE () const override { return nze; }
    virtual tuple<int,int> EntrySizes() const = 0;

    virtual shared_ptr<BaseMatrix> DeleteZeroElements(double tol) const override
    {
      throw Exception ("DeleteZeroElements not overloaded");
    }
  };

  
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
  class NGS_DLL_HEADER SparseMatrixTM : public S_BaseSparseMatrix<typename mat_traits<TM>::TSCAL>
  {
  protected:
    // Array<TM, size_t> data;
    NumaDistributedArray<TM> data;
    TM nul;
    bool hermitian = false;
    
    typedef S_BaseSparseMatrix<typename mat_traits<TM>::TSCAL> BASE;
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
    using typename BASE::ColIdx;
    
    void SetEntrySize()
    {
      // BASE::SetEntrySize (ngbla::Height<TM>(), ngbla::Width<TM>(), sizeof(TM)/sizeof(TSCAL));
      BASE::entry_height = ngbla::Height<TM>();
      BASE::entry_width = ngbla::Width<TM>();
      BASE::entry_size = sizeof(TM)/sizeof(TSCAL);
      BASE::is_complex = ngbla::IsComplex<TM>();
    }
    
  public:
    typedef TM TENTRY;
    typedef typename mat_traits<TM>::TSCAL TSCAL;

    SparseMatrixTM (int as, int max_elsperrow)
      : BASE (as, max_elsperrow),
	data(nze), nul(TSCAL(0))
    {
      // SetEntrySize (mat_traits<TM>::HEIGHT, mat_traits<TM>::WIDTH, sizeof(TM)/sizeof(TSCAL));
      // SetEntrySize (Height<TM>(), Width<TM>(), sizeof(TM)/sizeof(TSCAL));
      SetEntrySize ();
      asvec.AssignMemory (nze*sizeof(TM)/sizeof(TSCAL), (void*)data.Addr(0));
      GetMemoryTracer().Track(*static_cast<MatrixGraph*>(this), "MatrixGraph",
                              data, "data");
      GetMemoryTracer().SetName("SparseMatrix");
    }

    SparseMatrixTM (const Array<int> & elsperrow, int awidth)
      : BASE (elsperrow, awidth), 
	data(nze), nul(TSCAL(0))
    {
      // SetEntrySize (mat_traits<TM>::HEIGHT, mat_traits<TM>::WIDTH, sizeof(TM)/sizeof(TSCAL));
      SetEntrySize ();      
      asvec.AssignMemory (nze*sizeof(TM)/sizeof(TSCAL), (void*)data.Addr(0));
      GetMemoryTracer().Track(*static_cast<MatrixGraph*>(this), "MatrixGraph",
                              data, "data");
      GetMemoryTracer().SetName("SparseMatrix");

    }

    SparseMatrixTM (int size, int width, const Table<int> & rowelements, 
		    const Table<int> & colelements, bool symmetric)
      : BASE (size, width, rowelements, colelements, symmetric), 
	data(nze), nul(TSCAL(0))
    { 
      // SetEntrySize (mat_traits<TM>::HEIGHT, mat_traits<TM>::WIDTH, sizeof(TM)/sizeof(TSCAL));
      SetEntrySize ();            
      asvec.AssignMemory (nze*sizeof(TM)/sizeof(TSCAL), (void*)data.Addr(0));
      GetMemoryTracer().Track(*static_cast<MatrixGraph*>(this), "MatrixGraph",
                              data, "data");
      GetMemoryTracer().SetName("SparseMatrix");

    }

    SparseMatrixTM (const MatrixGraph & agraph)
      : BASE (agraph), 
	data(nze), nul(TSCAL(0))
    { 
      // SetEntrySize (mat_traits<TM>::HEIGHT, mat_traits<TM>::WIDTH, sizeof(TM)/sizeof(TSCAL));
      SetEntrySize ();      
      asvec.AssignMemory (nze*sizeof(TM)/sizeof(TSCAL), (void*)data.Addr(0));
      FindSameNZE();
      GetMemoryTracer().Track(*static_cast<MatrixGraph*>(this), "MatrixGraph",
                              data, "data");
      GetMemoryTracer().SetName("SparseMatrix");
    }

    SparseMatrixTM (MatrixGraph && agraph)
      : BASE (std::move(agraph)), 
	data(nze), nul(TSCAL(0))
    { 
      // SetEntrySize (mat_traits<TM>::HEIGHT, mat_traits<TM>::WIDTH, sizeof(TM)/sizeof(TSCAL));
      SetEntrySize ();      
      asvec.AssignMemory (nze*sizeof(TM)/sizeof(TSCAL), (void*)data.Addr(0));
      FindSameNZE();
      GetMemoryTracer().Track(*static_cast<MatrixGraph*>(this), "MatrixGraph",
                              data, "data");
      GetMemoryTracer().SetName("SparseMatrix");
    }
    
    SparseMatrixTM (const SparseMatrixTM & amat)
      : BASE (amat), 
      data(nze), nul(TSCAL(0))
    {
      // SetEntrySize (mat_traits<TM>::HEIGHT, mat_traits<TM>::WIDTH, sizeof(TM)/sizeof(TSCAL));
      SetEntrySize ();            
      asvec.AssignMemory (nze*sizeof(TM)/sizeof(TSCAL), (void*)data.Addr(0));      
      AsVector() = amat.AsVector();
      GetMemoryTracer().Track(*static_cast<MatrixGraph*>(this), "MatrixGraph",
                              data, "data");
      GetMemoryTracer().SetName("SparseMatrix");
    }

    SparseMatrixTM (SparseMatrixTM && amat)
      : BASE (std::move(amat)), nul(TSCAL(0))
    {
      // SetEntrySize (mat_traits<TM>::HEIGHT, mat_traits<TM>::WIDTH, sizeof(TM)/sizeof(TSCAL));
      SetEntrySize ();
      GetMemoryTracer().Track(*static_cast<MatrixGraph*>(this), "MatrixGraph",
                              data, "data");
      GetMemoryTracer().SetName("SparseMatrix");
      data.Swap(amat.data);
      asvec.AssignMemory (nze*sizeof(TM)/sizeof(TSCAL), (void*)data.Addr(0));            
    }

    static shared_ptr<SparseMatrixTM> CreateFromCOO (FlatArray<int> i, FlatArray<int> j,
                                                     FlatArray<TM> val, size_t h, size_t w);
      
    virtual ~SparseMatrixTM ();

    size_t Height() const { return size; }
    size_t Width() const { return width; }
    virtual int VHeight() const override { return size; }
    virtual int VWidth() const override { return width; }
    virtual tuple<size_t, size_t> Shape() const override { return { size, width }; }
    
    
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

    virtual AutoVector CreateVector () const override;
    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;

    // virtual tuple<int,int> EntrySizes() const override { return { mat_traits<TM>::HEIGHT, mat_traits<TM>::WIDTH }; }
    virtual tuple<int,int> EntrySizes() const override { return { ngbla::Height<TM>(), ngbla::Width<TM>() }; }
    
    shared_ptr<BaseSparseMatrix>
      CreateTransposeTM (const function<shared_ptr<SparseMatrixTM<decltype(ngbla::Trans(TM()))>>(const Array<int>&, int)> & creator) const;

  public:
    using BaseMatrix::GetMemoryTracer;
  };
  



  template<class TM, class TV_ROW, class TV_COL>
  class NGS_DLL_HEADER SparseMatrix :  public SparseMatrixTM<TM>
  {
  public:
    using SparseMatrixTM<TM>::firsti;
    using SparseMatrixTM<TM>::colnr;
    using SparseMatrixTM<TM>::data;
    using SparseMatrixTM<TM>::balance;
    using typename SparseMatrixTM<TM>::ColIdx;

    typedef typename mat_traits<TM>::TSCAL TSCAL;
    typedef TV_ROW TVX;
    typedef TV_COL TVY;

    ///
    SparseMatrix (int as, int max_elsperrow)
      : SparseMatrixTM<TM> (as, max_elsperrow) { ; }

    SparseMatrix (const Array<int> & aelsperrow)
      : SparseMatrixTM<TM> (aelsperrow, aelsperrow.Size()) { ; }

    SparseMatrix (const Array<int> & aelsperrow, int awidth)
      : SparseMatrixTM<TM> (aelsperrow, awidth) { ; }

    SparseMatrix (int height, int width, const Table<int> & rowelements, 
		  const Table<int> & colelements, bool symmetric)
      : SparseMatrixTM<TM> (height, width, rowelements, colelements, symmetric) { ; }

    SparseMatrix (const MatrixGraph & agraph);
    SparseMatrix (MatrixGraph && agraph);

    SparseMatrix (const SparseMatrix & amat)
      : SparseMatrixTM<TM> (amat) { ; }

    SparseMatrix (const SparseMatrixTM<TM> & amat)
      : SparseMatrixTM<TM> (amat) { ; }

    SparseMatrix (SparseMatrixTM<TM> && amat)
      : SparseMatrixTM<TM> (std::move(amat)) { ; }

    virtual shared_ptr<BaseMatrix> CreateMatrix () const override;
    // virtual BaseMatrix * CreateMatrix (const Array<int> & elsperrow) const;
    ///
    virtual AutoVector CreateVector () const override;
    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;


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

    virtual shared_ptr<BaseSparseMatrix> Restrict (const SparseMatrixTM<double> & prol,
					 shared_ptr<BaseSparseMatrix> cmat = nullptr) const override;
    
    virtual shared_ptr<BaseSparseMatrix> Reorder (const Array<size_t> & reorder) const override;
    
    virtual shared_ptr<BaseSparseMatrix> CreateTranspose() const override
    {
      return this->CreateTransposeTM
        ( [](const Array<int> & elsperrow, int width) -> shared_ptr<SparseMatrixTM<decltype(Trans(TM()))>>
          { return make_shared<SparseMatrix<decltype(Trans(TM())), TV_COL, TV_ROW>> (elsperrow, width); } );
    }

    virtual shared_ptr<BaseMatrix> DeleteZeroElements(double tol) const override;
    
    ///
    inline TVY RowTimesVector (int row, const FlatVector<TVX> vec) const
    {
      typedef typename mat_traits<TVY>::TSCAL TTSCAL;
      TVY sum = TTSCAL(0);
      for (size_t j = firsti[row]; j < firsti[row+1]; j++)
	sum += data[j] * vec(colnr[j]);
      return sum;
    }

    ///
    void AddRowTransToVector (int row, TVY el, FlatVector<TVX> vec) const
    {
      size_t first = firsti[row];
      size_t last = firsti[row+1];

      const ColIdx * colpi = colnr.Addr(0);
      const TM * datap = data.Addr(0);

      for (size_t j = first; j < last; j++)
        vec[colpi[j]] += Trans(datap[j]) * el; 
    }
    
    ///
    void AddRowConjTransToVector (int row, TVY el, FlatVector<TVX> vec) const
    {
      size_t first = firsti[row];
      size_t last = firsti[row+1];

      const ColIdx * colpi = colnr.Addr(0);
      const TM * datap = data.Addr(0);

      for (size_t j = first; j < last; j++)
        vec[colpi[j]] += Conj(Trans(datap[j])) * el; 
    }


    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultConjTransAdd (Complex s, const BaseVector & x, BaseVector & y) const override;

    virtual void MultAdd (FlatVector<double> alpha, const MultiVector & x, MultiVector & y) const override;

    
    virtual void MultAdd1 (double s, const BaseVector & x, BaseVector & y,
			   const BitArray * ainner = NULL,
			   const Array<int> * acluster = NULL) const override;
    
    virtual void DoArchive (Archive & ar) override;
  };

  

  /// A symmetric sparse matrix
  template<class TM, class TV>
  class NGS_DLL_HEADER SparseMatrixSymmetric : public SparseMatrix<TM,TV,TV>
  {

  public:
    using SparseMatrixTM<TM>::firsti;
    using SparseMatrixTM<TM>::colnr;
    using SparseMatrixTM<TM>::data;
    using BaseSparseMatrix::ColIdx;    

    typedef typename mat_traits<TM>::TSCAL TSCAL;
    typedef TV TV_COL;
    typedef TV TV_ROW;
    typedef TV TVY;
    typedef TV TVX;

    SparseMatrixSymmetric (int as, int max_elsperrow)
      : SparseMatrix<TM,TV,TV> (as, max_elsperrow)
    { ; }
  
    SparseMatrixSymmetric (const Array<int> & elsperrow)
      : SparseMatrix<TM,TV,TV> (elsperrow, elsperrow.Size())
    { ; }

    SparseMatrixSymmetric (int size, const Table<int> & rowelements)
      : SparseMatrix<TM,TV,TV> (size, size, rowelements, rowelements, true)
    { ; }

    SparseMatrixSymmetric (const MatrixGraph & agraph);
    SparseMatrixSymmetric (MatrixGraph && agraph);    

    SparseMatrixSymmetric (const SparseMatrixSymmetric & amat)
      : SparseMatrix<TM,TV,TV> (amat)
    { 
      this->AsVector() = amat.AsVector(); 
    }

    SparseMatrixSymmetric (const SparseMatrixTM<TM> & amat)
      : SparseMatrix<TM,TV,TV> (amat)
      { 
        this->AsVector() = amat.AsVector(); 
      }
    
    ///
    virtual ~SparseMatrixSymmetric ();

    SparseMatrixSymmetric & operator= (double s)
    {
      this->AsVector() = s;
      return *this;
    }

    virtual shared_ptr<BaseMatrix> CreateMatrix () const override
    {
      return make_shared<SparseMatrixSymmetric> (*this);
    }

    /*
    virtual BaseMatrix * CreateMatrix (const Array<int> & elsperrow) const
    {
      return new SparseMatrix<TM,TV,TV>(elsperrow);
    }
    */

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

    virtual shared_ptr<BaseSparseMatrix> Restrict (const SparseMatrixTM<double> & prol,
					 shared_ptr<BaseSparseMatrix> cmat = nullptr) const override;

    ///
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;

    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      MultAdd (s, x, y);
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
    




    using SparseMatrix<TM,TV,TV>::RowTimesVector;
    using SparseMatrix<TM,TV,TV>::AddRowTransToVector;


    TV_COL RowTimesVectorNoDiag (int row, const FlatVector<TVX> vec) const
    {
      size_t last = firsti[row+1];
      size_t first = firsti[row];
      if (last == first) return TVY(0);
      if (colnr[last-1] == row) last--;

      typedef typename mat_traits<TVY>::TSCAL TTSCAL;
      TVY sum = TTSCAL(0);

      for (size_t j = first; j < last; j++)
	sum += data[j] * vec(colnr[j]);
      return sum;
    }

    void AddRowTransToVectorNoDiag (int row, TVY el, FlatVector<TVX> vec) const
    {
      size_t first = firsti[row];
      size_t last = firsti[row+1];

      if (first == last) return;
      if (this->colnr[last-1] == row) last--;

      for (size_t j = first; j < last; j++)
        vec[colnr[j]] += Trans(data[j]) * el;
    }
  
    BaseSparseMatrix & AddMerge (double s, const SparseMatrixSymmetric  & m2);

    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = nullptr) const override;
    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<const Array<int>> clusters) const override;
  };

  [[deprecated("Use sparsematrix->CreateTranspose() instead!")]]            
  NGS_DLL_HEADER shared_ptr<SparseMatrixTM<double>> TransposeMatrix (const SparseMatrixTM<double> & mat);

  NGS_DLL_HEADER shared_ptr<SparseMatrixTM<double>>
  MatAdd (double sa, const SparseMatrixTM<double> & mata,
          double sb, const SparseMatrixTM<double> & matb);
  
  NGS_DLL_HEADER shared_ptr<SparseMatrixTM<double>>
  MatMult (const SparseMatrixTM<double> & mata, const SparseMatrixTM<double> & matb, bool sort_output = true);
  NGS_DLL_HEADER shared_ptr<SparseMatrixTM<Complex>>
  MatMult (const SparseMatrixTM<Complex> & mata, const SparseMatrixTM<Complex> & matb, bool sort_output = true);

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
  SPARSEMATRIX_EXTERN template class SparseMatrix<double, Complex, Complex>;
  
  SPARSEMATRIX_EXTERN template class SparseMatrixSymmetric<double>;
  SPARSEMATRIX_EXTERN template class SparseMatrixSymmetric<Complex>;
  SPARSEMATRIX_EXTERN template class SparseMatrixSymmetric<double, Complex>;

#define INST_SPMS(N) \
  SPARSEMATRIX_EXTERN template class SparseMatrix<Mat<N, N, double>>; \
  SPARSEMATRIX_EXTERN template class SparseMatrix<Mat<1, N, double>>; \
  SPARSEMATRIX_EXTERN template class SparseMatrix<Mat<N, 1, double>>; \
  SPARSEMATRIX_EXTERN template class SparseMatrix<Mat<N, N, Complex>>; \
  SPARSEMATRIX_EXTERN template class SparseMatrix<Mat<1, N, Complex>>; \
  SPARSEMATRIX_EXTERN template class SparseMatrix<Mat<N, 1, Complex>>; \
  SPARSEMATRIX_EXTERN template class SparseMatrixSymmetric<Mat<N, N, double>>; \
  SPARSEMATRIX_EXTERN template class SparseMatrixSymmetric<Mat<N, N, Complex>>;


#if MAX_SYS_DIM >= 1
  INST_SPMS(1);
#endif
#if MAX_SYS_DIM >= 2
  INST_SPMS(2);
#endif
#if MAX_SYS_DIM >= 3
  INST_SPMS(3);
#endif
#if MAX_SYS_DIM >= 4
  INST_SPMS(4);
#endif
#if MAX_SYS_DIM >= 5
  INST_SPMS(5);
#endif
#if MAX_SYS_DIM >= 6
  INST_SPMS(6);
#endif
#if MAX_SYS_DIM >= 7
  INST_SPMS(7);
#endif
#if MAX_SYS_DIM >= 8
  INST_SPMS(8);
#endif

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
          GetMemoryTracer().Track(*static_cast<MatrixGraph*>(this), "MatrixGraph",
                                  data, "data");
          GetMemoryTracer().SetName("SparseMatrix");
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
          GetMemoryTracer().Track(*static_cast<MatrixGraph*>(this), "MatrixGraph",
                                  data, "data");
          GetMemoryTracer().SetName("SparseMatrix");
        }
    

    
    tuple<int,int> EntrySizes() const override { return { bheight, bwidth }; }
    
    AutoVector CreateRowVector () const override
    {
      return AutoVector(make_shared<S_BaseVectorPtr<TSCAL>> (this->width, this->bwidth));
    }
    
    AutoVector CreateColVector () const override
    {
      return AutoVector(make_shared<S_BaseVectorPtr<TSCAL>> (this->size, this->bheight));
    }

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

    
    const MemoryTracer & GetMemoryTracer() const
    {
      return mem_tracer;
    }
    
  private:
    MemoryTracer mem_tracer =
      {"MatrixGraph",
       colnr, "colnr",
       firsti, "firsti"
       // same_nze, "same_nze"
      };
    
  };


  

}

#endif
