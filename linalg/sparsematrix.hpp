#ifndef FILE_NGS_SPARSEMATRIX
#define FILE_NGS_SPARSEMATRIX

/* ************************************************************************/
/* File:   sparsematrix.hpp                                               */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Oct. 94, 15 Jan. 02                                        */
/* ************************************************************************/

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

#ifdef PARALLEL
  class ParallelBaseMatrix;

  template<class TM, 
	   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class ParallelSparseMatrix;

  template<class TM, 
	   class TV = typename mat_traits<TM>::TV_ROW>
  class ParallelSparseMatrixSymmetric;
#endif


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
 
#ifdef PARALLEL
  template<class TM, 
	   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class ParallelBlockJacobiPrecond;

  template<class TM, 
	   class TV = typename mat_traits<TM>::TV_ROW>
  class ParallelBlockJacobiPrecondSymmetric;
#endif



  template<class TM, 
	   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class SparseCholesky;

  template<class TM, 
	   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class PardisoInverse;

  template<class TM, 
	   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class MumpsInverse;

  template<class TM, 
	   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class SuperLUInverse;


  template<class TM, 
	   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class SuperLU_DIST_Inverse;
  //template <class TM> class SuperLUInverse;

  // sets the solver which is used for InverseMatrix
  enum INVERSETYPE { PARDISO, PARDISOSPD, SPARSECHOLESKY, SUPERLU, SUPERLU_DIST, MUMPS };






  /** 
      The graph of a sparse matrix.
  */
  class MatrixGraph
  {
  protected:
    /// number of rows
    int size;

    /// non-zero elements
    int nze; 

    /// column numbers
    // int * colnr;
    DynamicMem<int> colnr;

    /// pointer to first in row
    // int * firsti;
    DynamicMem<int> firsti;
  
    /// position of diagonal element
    //   int * diagi;
    DynamicMem<int> diagi;

    /// row has same non-zero elements as previous row
    Array<int> same_nze;

    /// owner of arrays ?
    bool owner;

    /// sparse direct solver
    INVERSETYPE inversetype;

  public:
    /// matrix of hight as, uniform number of els/row
    MatrixGraph (int as, int max_elsperrow);
    /// arbitrary number of els/row
    MatrixGraph (const Array<int> & elsperrow);
    /// shadow matrix graph
    MatrixGraph (const MatrixGraph & graph, bool stealgraph);
    /// 
    MatrixGraph (int size, const Table<int> & elements, bool symmetric);

    virtual ~MatrixGraph ();

    /// eliminate unused columne indices
    void Compress();
  
    /// returns position of Element (i, j), exception for unused
    int GetPosition (int i, int j) const;

    /// returns position of Element (i, j), -1 for unused
    int GetPositionTest (int i, int j) const;

    /// find positions of n sorted elements, overwrite pos, exception for unused
    void GetPositionsSorted (int row, int n, int * pos) const;

    /// returns position of new element
    int CreatePosition (int i, int j);

    int Size() const { return size; }

    int NZE() const { return nze; }

    FlatArray<const int> GetRowIndices(int i) const
    { return FlatArray<const int> (firsti[i+1]-firsti[i], colnr+firsti[i]); }

    FlatArray<int> GetRowIndices(int i) 
    { return FlatArray<int> (firsti[i+1]-firsti[i], colnr+firsti[i]); }


    int * GetRowIndicesPointer (int i) 
    { return colnr+firsti[i]; }

    int First (int i) const { return firsti[i]; }

    void FindSameNZE();

    ostream & Print (ostream & ost) const;

    virtual void MemoryUsage (Array<MemoryUsageStruct*> & mu) const;

    virtual void SetInverseType ( INVERSETYPE ainversetype )
    { inversetype = ainversetype; }

    virtual void SetInverseType ( string ainversetype ) const;

    virtual INVERSETYPE  GetInverseType () const
    { return inversetype; }

    virtual int GetNRowEntries ( const int & row )
    { return (firsti[row+1] - firsti[row]); }

#ifdef ASTRID
    MatrixGraph * CreatePartialGraph ( BitArray & rows, BitArray & cols );
#endif
  };








  /// A virtual base class for all sparse matrices
  class BaseSparseMatrix : virtual public BaseMatrix, 
			   public MatrixGraph
  {

  public:
    BaseSparseMatrix (int as, int max_elsperrow)
      : MatrixGraph (as, max_elsperrow)
    { ; }

    BaseSparseMatrix (const Array<int> & elsperrow)
      : MatrixGraph (elsperrow)
    { ; }

    BaseSparseMatrix (const MatrixGraph & agraph, bool stealgraph)
      : MatrixGraph (agraph, stealgraph)
    { ; }   


    BaseSparseMatrix (const BaseSparseMatrix & amat)
      : BaseMatrix(amat), MatrixGraph (amat, 0)
    { ; }   



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


    virtual BaseJacobiPrecond * CreateJacobiPrecond (const BitArray * inner = 0) const 
    {
      throw Exception ("BaseSparseMatrix::CreateJacobiPrecond");
    }

    virtual BaseBlockJacobiPrecond * 
    CreateBlockJacobiPrecond (Table<int> & blocks,
			      const BaseVector * constraint = 0,
			      const ngcomp::Preconditioner * acoarsegridprecond = 0, bool parallel  = 1 ) const
    { 
      throw Exception ("BaseSparseMatrix::CreateBlockJacobiPrecond");
    }


    virtual BaseMatrix * 
    InverseMatrix (BitArray * subset = 0) const
    { 
      throw Exception ("BaseSparseMatrix::CreateInverse called");
    }

    virtual BaseMatrix * 
    InverseMatrix (Array<int> * clusters) const
    { 
      throw Exception ("BaseSparseMatrix::CreateInverse called");
    }

    virtual BaseSparseMatrix * Restrict (const SparseMatrixTM<double> & prol,
					 BaseSparseMatrix* cmat = NULL ) const
    {
      throw Exception ("BaseSparseMatrix::Restrict");
    }

    //virtual void AllReduceHO () const = 0;


  };




  /// A general, sparse matrix
  template<class TM>
  class SparseMatrixTM : public BaseSparseMatrix, 
			 public S_BaseMatrix<typename mat_traits<TM>::TSCAL>
  {
  protected:
    DynamicMem<TM> data;
    VFlatVector<typename mat_traits<TM>::TSCAL> asvec;
    TM nul;

  public:
    typedef typename mat_traits<TM>::TSCAL TSCAL;

    SparseMatrixTM (int as, int max_elsperrow);
    SparseMatrixTM (const Array<int> & elsperrow);
    SparseMatrixTM (const MatrixGraph & agraph, bool stealgraph);
    SparseMatrixTM (const SparseMatrixTM & amat);
    virtual ~SparseMatrixTM ();

    int Height() const { return size; }
    int Width() const { return size; }
    virtual int VHeight() const { return size; }
    virtual int VWidth() const { return size; }

    TM & operator[] (int i)  { return data[i]; }
    const TM & operator[] (int i) const { return data[i]; }


    FlatArray<const TM> GetRow(int i) const
    { 
      return FlatArray<const TM> (firsti[i+1]-firsti[i], data+firsti[i]); 
    }

    TM & operator() (int row, int col)
    {
      return data[CreatePosition(row, col)];
    }

    const TM & operator() (int row, int col) const
    {
      int pos = GetPositionTest (row,col);
      if (pos != -1)
	return data[pos];
      else
	return nul;
    }

    FlatVector<const TM> GetRowValues(int i) const
    { return FlatVector<const TM> (firsti[i+1]-firsti[i], &data[firsti[i]]); }

    FlatVector<TM> GetRowValues(int i) 
    { return FlatVector<TM> (firsti[i+1]-firsti[i], &data[firsti[i]]); }

    const TM & GetRowValue(int i, int j) const
    { return data[firsti[i]+j]; }

    TM & GetRowValue(int i, int j) 
    { return data[firsti[i]+j]; }


    void AddElementMatrix(const Array<int> & dnums1, const Array<int> & dnums2, const FlatMatrix<TSCAL> & elmat);

    virtual BaseVector & AsVector() 
    {
      asvec.AssignMemory (nze*sizeof(TM)/sizeof(TSCAL), (void*)&data[0]);
      return asvec; 
    }

    virtual const BaseVector & AsVector() const
    { 
      const_cast<VFlatVector<TSCAL>&> (asvec).
	AssignMemory (nze*sizeof(TM)/sizeof(TSCAL), (void*)&data[0]);
      return asvec; 
    }

    ///
    virtual ostream & Print (ostream & ost) const;

    ///
    virtual void MemoryUsage (Array<MemoryUsageStruct*> & mu) const;


  };



  template<class TM, class TV_ROW, class TV_COL>
  class SparseMatrix : virtual public SparseMatrixTM<TM>
  {

  public:
    typedef typename mat_traits<TM>::TSCAL TSCAL;
    typedef TV_ROW TVX;
    typedef TV_COL TVY;
    typedef typename mat_scale_type<TVX,Complex>::TMAT TVXC;
    typedef typename mat_scale_type<TVY,Complex>::TMAT TVYC;

    ///
    SparseMatrix (int as, int max_elsperrow)
      : SparseMatrixTM<TM> (as, max_elsperrow) { ; }

    SparseMatrix (const Array<int> & elsperrow)
      : SparseMatrixTM<TM> (elsperrow) { ; }

    SparseMatrix (const MatrixGraph & agraph, bool stealgraph)
      : SparseMatrixTM<TM> (agraph, stealgraph) { ; }

    SparseMatrix (const SparseMatrix & amat)
      : SparseMatrixTM<TM> (amat) { ; }

    SparseMatrix (const SparseMatrixTM<TM> & amat)
      : SparseMatrixTM<TM> (amat) { ; }

    virtual BaseMatrix * CreateMatrix () const;
    virtual BaseMatrix * CreateMatrix (const Array<int> & elsperrow) const;
    ///
    virtual BaseVector * CreateVector () const;

    virtual BaseJacobiPrecond * 
    CreateJacobiPrecond (const BitArray * inner) const
    { 
      return new JacobiPrecond<TM,TV_ROW,TV_COL> (*this, inner);
    }

    virtual BaseBlockJacobiPrecond * 
    CreateBlockJacobiPrecond (Table<int> & blocks,
			      const BaseVector * constraint = 0,
			      const ngcomp::Preconditioner * acoarsegridprecond = 0, bool parallel  = 1 ) const
    { 
      return new BlockJacobiPrecond<TM,TV_ROW,TV_COL> (*this, blocks );
    }

    virtual BaseMatrix * InverseMatrix (BitArray * subset = 0) const;

    virtual BaseMatrix * InverseMatrix (Array<int> * clusters) const;

    virtual BaseSparseMatrix * Restrict (const SparseMatrixTM<double> & prol,
					 BaseSparseMatrix* cmat = NULL ) const
    { return 0; }


  
    ///
    // template<class TEL>
    inline TVY RowTimesVector (int row, const FlatVector<TVX> & vec) const
    {
      /*
	TVY sum = TSCAL(0);
	for (int j = firsti[row]; j < firsti[row+1]; j++)
	sum += data[j] * vec(colnr[j]);
	return sum;
      */

      int nj = (this->firsti)[row+1] - (this->firsti)[row];
      const int * colpi = (this->colnr)+(this->firsti)[row];
      const TM * datap = (this->data)+(this->firsti)[row];
      const TVX * vecp = &vec(0);
    
      typedef typename mat_traits<TVY>::TSCAL TTSCAL;
      TVY sum = TTSCAL(0);
      for (int j = 0; j < nj; j++, colpi++, datap++)
	sum += *datap * vecp[*colpi];
      return sum;
    }

    ///
    void AddRowTransToVector (int row, TVY el, FlatVector<TVX> & vec) const
    {
      /*
      // cannot be optimized by the compiler, since hvec could alias
      // the matrix --> keyword 'restrict' will hopefully solve this problem
      
      for (int j = firsti[row]; j < firsti[row+1]; j++)
      hvec(colnr[j]) += Trans(data[j]) * el;
      */

      int nj = (this->firsti)[row+1] - (this->firsti)[row];
      const int * colpi = (this->colnr)+(this->firsti)[row];
      const TM * datap = (this->data)+(this->firsti)[row];
      TVX * vecp = &vec(0);
      for (int j = 0; j < nj; j++, colpi++, datap++)
	vecp[*colpi] += Trans(*datap) * el;
    }


    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const;
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const;
    virtual void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const;


  };




  /// A symmetric sparse matrix
  template<class TM>
  class SparseMatrixSymmetricTM : virtual public SparseMatrixTM<TM>
  {

  protected:

    SparseMatrixSymmetricTM (int as, int max_elsperrow)
      : SparseMatrixTM<TM> (as, max_elsperrow) { ; }

    SparseMatrixSymmetricTM (const Array<int> & elsperrow)
      : SparseMatrixTM<TM> (elsperrow) { ; }

    SparseMatrixSymmetricTM (const MatrixGraph & agraph, bool stealgraph)
      : SparseMatrixTM<TM> (agraph, stealgraph) { ; }

    SparseMatrixSymmetricTM (const SparseMatrixSymmetricTM & amat)
      : SparseMatrixTM<TM> (amat) { ; }


  public:
    typedef typename mat_traits<TM>::TSCAL TSCAL;
    virtual void AddElementMatrix(const Array<int> & dnums, const FlatMatrix<TSCAL> & elmat);
  };


  /// A symmetric sparse matrix
  template<class TM, class TV>
  class SparseMatrixSymmetric : virtual public SparseMatrixSymmetricTM<TM>, 
				virtual public SparseMatrix<TM,TV,TV>
  {

  public:
    typedef typename mat_traits<TM>::TSCAL TSCAL;
    typedef TV TV_COL;
    typedef TV TV_ROW;
    typedef TV TVY;
    typedef TV TVX;

    /*
   ///
   SparseMatrixSymmetric (int as, int max_elsperrow);
   ///
   SparseMatrixSymmetric (const Array<int> & elsperrow);
   ///
   SparseMatrixSymmetric (const MatrixGraph & agraph, bool stealgraph);
   ///
   SparseMatrixSymmetric (const SparseMatrixSymmetric & agraph);
    */

    SparseMatrixSymmetric (int as, int max_elsperrow)
      : SparseMatrixTM<TM> (as, max_elsperrow) , 
	SparseMatrixSymmetricTM<TM> (as, max_elsperrow),
	SparseMatrix<TM,TV,TV> (as, max_elsperrow)
    { ; }
  
    SparseMatrixSymmetric (const Array<int> & elsperrow)
      : SparseMatrixTM<TM> (elsperrow), 
	SparseMatrixSymmetricTM<TM> (elsperrow),
	SparseMatrix<TM,TV,TV> (elsperrow)
    { ; }

    SparseMatrixSymmetric (const MatrixGraph & agraph, bool stealgraph)
      : SparseMatrixTM<TM> (agraph, stealgraph), 
	SparseMatrixSymmetricTM<TM> (agraph, stealgraph),
	SparseMatrix<TM,TV,TV> (agraph, stealgraph)
    { ; }


    SparseMatrixSymmetric (const SparseMatrixSymmetric & amat)
      : SparseMatrixTM<TM> (amat), 
	SparseMatrixSymmetricTM<TM> (amat),
	SparseMatrix<TM,TV,TV> (amat)
    { 
      this->AsVector() = amat.AsVector(); 
    }

    SparseMatrixSymmetric (const SparseMatrixSymmetricTM<TM> & amat)
      : SparseMatrixTM<TM> (amat), 
	SparseMatrixSymmetricTM<TM> (amat),
	SparseMatrix<TM,TV,TV> (amat)
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

    virtual BaseMatrix * CreateMatrix () const
    {
      return new SparseMatrixSymmetric (*this);
    }

    virtual BaseMatrix * CreateMatrix (const Array<int> & elsperrow) const
    {
      return new SparseMatrix<TM,TV,TV>(elsperrow);
    }

    virtual BaseJacobiPrecond * CreateJacobiPrecond (const BitArray * inner) const
    { 
      return new JacobiPrecondSymmetric<TM,TV> (*this, inner);
    }

    virtual BaseBlockJacobiPrecond * 
    CreateBlockJacobiPrecond (Table<int> & blocks,
			      const BaseVector * constraint = 0,
			      const ngcomp::Preconditioner * acoarsegridprecond = 0, bool parallel  = 1 ) const
    { 
      if (!constraint)
	return new BlockJacobiPrecondSymmetric<TM,TV> (*this, blocks);
      else
	return new BlockJacobiPrecondSymmetric<TM,TV> 
	  (*this, dynamic_cast<const T_BaseVector<TVX>&>(*constraint).FV(),
	   blocks);
    }


    virtual BaseSparseMatrix * Restrict (const SparseMatrixTM<double> & prol,
					 BaseSparseMatrix* cmat = NULL ) const;

    ///
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;

    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
    {
      MultAdd (s, x, y);
    }


    /*
      y += s L * x
    */
    virtual void MultAdd1 (double s, const BaseVector & x, BaseVector & y) const;


    /*
      y += s (D + L^T) * x
    */
    virtual void MultAdd2 (double s, const BaseVector & x, BaseVector & y) const;





    TV_COL RowTimesVectorNoDiag (int row, const FlatVector<TVX> & vec) const
    {
      int last = this->firsti[row+1];
      if (this->colnr[last-1] == row) last--;

      /*
	TV_COL sum = TSCAL(0);
	for (int j = this->firsti[row]; j < last; j++)
	sum += this->data[j] * vec(this->colnr[j]);
	return sum;
      */

      const int * colpi = this->colnr+this->firsti[row];
      const TM * datap = this->data+this->firsti[row];
      const TVX * vecp = &vec(0);

      typedef typename mat_traits<TVY>::TSCAL TTSCAL;
      TVY sum = TTSCAL(0);
      for (int j = this->firsti[row]; j < last; j++, colpi++, datap++)
	sum += *datap * vecp[*colpi];
      return sum;
    }

    void AddRowTransToVectorNoDiag (int row, TVY el, FlatVector<TVX> & vec) const
    {
      /*
	int last = firsti[row+1];
	if (colnr[last-1] == row) last--;
    
	for (int j = firsti[row]; j < last; j++)
	vec(colnr[j]) += Trans(data[j]) * el;
      */
      int last = this->firsti[row+1];
      if (this->colnr[last-1] == row) last--;

      int nj = last - this->firsti[row];
      const int * colpi = this->colnr+this->firsti[row];
      const TM * datap = this->data+this->firsti[row];
      TVX * vecp = &vec(0);
      for (int j = 0; j < nj; j++, colpi++, datap++)
	vecp[*colpi] += Trans(*datap) * el;
    }
  
    BaseSparseMatrix & AddMerge (double s, const SparseMatrixSymmetric  & m2);

    virtual BaseMatrix * InverseMatrix (BitArray * subset = 0) const;

    virtual BaseMatrix * InverseMatrix (Array<int> * clusters) const;

  };











  template<class TM>
  class VarBlockSparseMatrix : public BaseSparseMatrix, 
			       public S_BaseMatrix<typename mat_traits<TM>::TSCAL>
  {
    Array<int> block2linear;
    Array<int> data_index;
    DynamicMem<TM> data;

  public:

    typedef typename mat_traits<TM>::TSCAL TSCAL;
    typedef typename mat_traits<TM>::TV_COL TV_COL;
    typedef typename mat_traits<TM>::TV_ROW TV_ROW;
    typedef typename mat_traits<TM>::TV_COL TVY;
    typedef typename mat_traits<TM>::TV_ROW TVX;



    VarBlockSparseMatrix (Array<int> & elsperrow, 
			  Array<int> & ablock2linear, 
			  Array<int> & linear2block, 
			  const SparseMatrix<TM> & sm);
    static VarBlockSparseMatrix * Create (const SparseMatrix<TM> & sm);
    virtual ~VarBlockSparseMatrix ();
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;


  };




  // ***************************************************************

  //               PARALLEL MATRICES

#ifdef PARALLEL

  template<class TM, class TV_ROW, class TV_COL>
  class ParallelSparseMatrix : virtual public SparseMatrix<TM,TV_ROW,TV_COL>, 
			       public ParallelBaseMatrix
  {
    SparseMatrix<TM,TV_ROW,TV_COL> * consistentmat;

  public:
    typedef typename mat_traits<TM>::TSCAL TSCAL;
    typedef TV_ROW TVX;
    typedef TV_COL TVY;
    typedef typename mat_scale_type<TVX,Complex>::TMAT TVXC;
    typedef typename mat_scale_type<TVY,Complex>::TMAT TVYC;


    ParallelSparseMatrix (int as, int max_elsperrow)
      :  SparseMatrixTM<TM> (as, max_elsperrow) ,SparseMatrix<TM,TV_ROW,TV_COL> (as, max_elsperrow)
    { consistentmat = 0; }

    ParallelSparseMatrix (const Array<int> & elsperrow)
      :  SparseMatrixTM<TM> (elsperrow) ,SparseMatrix<TM,TV_ROW,TV_COL> (elsperrow) 
    { consistentmat = 0; }

    ParallelSparseMatrix (const MatrixGraph & agraph, bool stealgraph)
      :  SparseMatrixTM<TM> (agraph, stealgraph) ,SparseMatrix<TM,TV_ROW,TV_COL> (agraph, stealgraph) 
    { consistentmat = 0; }

    ParallelSparseMatrix (const MatrixGraph & agraph, bool stealgraph,
			  const ngparallel::ParallelDofs * aparalleldofs)
      :  SparseMatrixTM<TM> (agraph, stealgraph) ,
	 SparseMatrix<TM,TV_ROW,TV_COL> (agraph, stealgraph),
	 ParallelBaseMatrix (aparalleldofs)
    { consistentmat = 0; }

    ParallelSparseMatrix (const ParallelSparseMatrix & amat);

    ParallelSparseMatrix (int as, int max_elsperrow, SparseMatrix<TM,TV_ROW,TV_COL> * aconsistentmat)
      : SparseMatrixTM<TM> (as, max_elsperrow) ,SparseMatrix<TM,TV_ROW,TV_COL>(as, max_elsperrow) 
    { consistentmat = aconsistentmat; }

    ParallelSparseMatrix (const Array<int> & elsperrow, SparseMatrix<TM,TV_ROW,TV_COL> * aconsistentmat)
      : SparseMatrixTM<TM> (elsperrow),SparseMatrix<TM,TV_ROW,TV_COL>(elsperrow)  
    { consistentmat = aconsistentmat; }

    ParallelSparseMatrix (const MatrixGraph & agraph, bool stealgraph, 
			  SparseMatrix<TM,TV_ROW,TV_COL> * aconsistentmat)
      : SparseMatrixTM<TM> (agraph, stealgraph),SparseMatrix<TM,TV_ROW,TV_COL>(agraph,stealgraph) 
    { consistentmat = aconsistentmat; }

    ParallelSparseMatrix (const MatrixGraph & agraph, bool stealgraph, 
			  const ngparallel::ParallelDofs * aparalleldofs,
			  SparseMatrix<TM,TV_ROW,TV_COL> * aconsistentmat)
      : SparseMatrixTM<TM> (agraph, stealgraph),
	SparseMatrix<TM,TV_ROW,TV_COL>(agraph,stealgraph),
	ParallelBaseMatrix (aparalleldofs) 
    { consistentmat = aconsistentmat; }

    ~ParallelSparseMatrix ()
    { if ( consistentmat ) delete consistentmat; }

    virtual BaseMatrix * ConsistentMat () { return consistentmat; }
    virtual const BaseMatrix * ConsistentMat () const  { return consistentmat; }

    virtual void SetConsistentMat ( BaseMatrix * aconsistentmat )
    { consistentmat = dynamic_cast<SparseMatrix<TM,TV_ROW,TV_COL> *> (aconsistentmat); }

    virtual void AllocateConsistentMat ();
    virtual void  AllocateConsistentMat ( const MatrixGraph & graph );

    virtual void CalcConsistentMat (LocalHeap & lh);

    virtual BaseBlockJacobiPrecond * 
    CreateBlockJacobiPrecond ( Table<int> & blocks,
			       const BaseVector * constraint = 0, 
			       const ngcomp::Preconditioner * acoarsegridprecond = 0, bool parallel = 1) const;
    //   { 
    //     if ( parallel )
    //       return new ParallelBlockJacobiPrecond<TM,TV_ROW,TV_COL> 
    // 	(*this, blocks, &acoarsegridprecond->GetMatrix());
    //     else
    //       return new BlockJacobiPrecond<TM,TV_ROW,TV_COL>
    // 	(*this, blocks);
    //   }

    virtual BaseMatrix * CreateMatrix () const;
    virtual BaseMatrix * CreateMatrix (const Array<int> & elsperrow) const;
    ///
    virtual BaseVector * CreateVector () const;

    virtual BaseMatrix * InverseMatrix (BitArray * subset = 0) const;

    virtual BaseMatrix * InverseMatrix (Array<int> * clusters) const;

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const;
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const;
    virtual void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const;

  };



  template<class TM, class TV>
  class ParallelSparseMatrixSymmetric : public SparseMatrixSymmetric<TM,TV>,
					public ParallelSparseMatrix<TM,TV,TV>
  {
    SparseMatrixSymmetric<TM,TV> * consistentmat;
  public:
    typedef typename mat_traits<TM>::TSCAL TSCAL;
    typedef TV TV_COL;
    typedef TV TV_ROW;
    typedef TV TVY;
    typedef TV TVX;


    ParallelSparseMatrixSymmetric (int as, int max_elsperrow)
      : SparseMatrixTM<TM> (as, max_elsperrow) , 
	SparseMatrix<TM,TV,TV> (as, max_elsperrow),
	SparseMatrixSymmetricTM<TM> (as, max_elsperrow),
	ParallelSparseMatrix<TM,TV,TV> (as, max_elsperrow),
	SparseMatrixSymmetric<TM,TV> (as, max_elsperrow) 
    {     consistentmat = 0; }
  
    ParallelSparseMatrixSymmetric (const Array<int> & elsperrow)
      : SparseMatrixTM<TM> (elsperrow) , 
	SparseMatrix<TM,TV,TV> (elsperrow),
	SparseMatrixSymmetricTM<TM> (elsperrow),
	SparseMatrixSymmetric<TM,TV> (elsperrow),
	ParallelSparseMatrix<TM,TV,TV> (elsperrow)
    {     consistentmat = 0; }

    ParallelSparseMatrixSymmetric (const MatrixGraph & agraph, bool stealgraph)
      : SparseMatrixTM<TM> (agraph, stealgraph) , 
	SparseMatrixSymmetricTM<TM> (agraph, stealgraph),
	SparseMatrix<TM,TV,TV> (agraph, stealgraph),
	SparseMatrixSymmetric<TM,TV> (agraph, stealgraph),
	ParallelSparseMatrix<TM,TV,TV> (agraph, stealgraph)
    {     consistentmat = 0; }

    ParallelSparseMatrixSymmetric (const MatrixGraph & agraph, bool stealgraph,
				   const ngparallel:: ParallelDofs * aparalleldofs)
      : SparseMatrixTM<TM> (agraph, stealgraph) , 
	SparseMatrixSymmetricTM<TM> (agraph, stealgraph),
	SparseMatrix<TM,TV,TV> (agraph, stealgraph),
	SparseMatrixSymmetric<TM,TV> (agraph, stealgraph),
	ParallelSparseMatrix<TM,TV,TV> (agraph, stealgraph,aparalleldofs)
    {     consistentmat = 0; }


    ParallelSparseMatrixSymmetric (const ParallelSparseMatrixSymmetric & amat);

    ParallelSparseMatrixSymmetric (int as, int max_elsperrow, SparseMatrixSymmetric<TM,TV> * aconsistentmat )
      : SparseMatrixTM<TM> (as, max_elsperrow) , 
	SparseMatrix<TM,TV,TV> (as, max_elsperrow),
	SparseMatrixSymmetricTM<TM> (as, max_elsperrow),
	SparseMatrixSymmetric<TM,TV> (as, max_elsperrow) ,
	ParallelSparseMatrix<TM,TV,TV> (as, max_elsperrow)
    {     consistentmat = aconsistentmat; }
  
    ParallelSparseMatrixSymmetric (const Array<int> & elsperrow, SparseMatrixSymmetric<TM,TV> * aconsistentmat )
      : SparseMatrixTM<TM> (elsperrow) , 
	SparseMatrix<TM,TV,TV> (elsperrow),
	SparseMatrixSymmetricTM<TM> (elsperrow),
	SparseMatrixSymmetric<TM,TV> (elsperrow),
	ParallelSparseMatrix<TM,TV,TV> (elsperrow)
    {     consistentmat = aconsistentmat; }

    ParallelSparseMatrixSymmetric (const MatrixGraph & agraph, bool stealgraph, SparseMatrixSymmetric<TM,TV> * aconsistentmat )
      : SparseMatrixTM<TM> (agraph, stealgraph) , 
	SparseMatrix<TM,TV,TV> (agraph, stealgraph),
	SparseMatrixSymmetricTM<TM> (agraph, stealgraph),
	ParallelSparseMatrix<TM,TV,TV> (agraph, stealgraph),
	SparseMatrixSymmetric<TM,TV> (agraph, stealgraph)
    {     consistentmat = aconsistentmat; }

    ParallelSparseMatrixSymmetric (const MatrixGraph & agraph, bool stealgraph, 
				   const ngparallel::ParallelDofs * aparalleldofs,
				   SparseMatrixSymmetric<TM,TV> * aconsistentmat )
      : SparseMatrixTM<TM> (agraph, stealgraph) , 
	SparseMatrix<TM,TV,TV> (agraph, stealgraph),
	SparseMatrixSymmetricTM<TM> (agraph, stealgraph),
	ParallelSparseMatrix<TM,TV,TV> (agraph, stealgraph,aparalleldofs),
	SparseMatrixSymmetric<TM,TV> (agraph, stealgraph)
    {     consistentmat = aconsistentmat; }








    ///
    virtual ~ParallelSparseMatrixSymmetric ()
    {
      if ( consistentmat )
	delete consistentmat;
    }

    virtual BaseBlockJacobiPrecond * 
    CreateBlockJacobiPrecond ( Table<int> & blocks,
			       const BaseVector * constraint = 0, 
			       const ngcomp::Preconditioner * acoarsegridprecond = 0, bool parallel  = 1) const;
    //   { 
    //     if ( parallel )
    //       {
    // 	if (!constraint)
    // 	  return new ParallelBlockJacobiPrecondSymmetric<TM,TV> 
    // 	    (*this, blocks, &acoarsegridprecond->GetMatrix());
    // 	else
    // 	  return new ParallelBlockJacobiPrecondSymmetric<TM,TV> 
    // 	    (*this,
    // 	     dynamic_cast<const T_BaseVector<TVX>&>(*constraint).FV(),
    // 	     blocks, &acoarsegridprecond->GetMatrix());
    //       }
    //     else
    //       {
    // 	if (!constraint)
    // 	  return new BlockJacobiPrecondSymmetric<TM,TV> 
    // 	    (*this, blocks);
    // 	else
    // 	  return new BlockJacobiPrecondSymmetric<TM,TV> 
    // 	    (*this,
    // 	     dynamic_cast<const T_BaseVector<TVX>&>(*constraint).FV(),
    // 	     blocks);
    //       }

    //   }

    virtual BaseMatrix * ConsistentMat () { return consistentmat; }
    virtual const BaseMatrix * ConsistentMat () const  { return consistentmat; }

    virtual void SetConsistentMat ( BaseMatrix * aconsistentmat )
    { consistentmat = dynamic_cast<SparseMatrixSymmetric<TM,TV> *> (aconsistentmat); }

    virtual void AllocateConsistentMat ();
    virtual void  AllocateConsistentMat ( const MatrixGraph & graph );

    virtual void CalcConsistentMat (LocalHeap & lh);

    virtual BaseMatrix * CreateMatrix () const;
    virtual BaseMatrix * CreateMatrix (const Array<int> & elsperrow) const;
    ///
    virtual BaseVector * CreateVector () const;

    virtual BaseMatrix * InverseMatrix (BitArray * subset = 0) const;

    virtual BaseMatrix * InverseMatrix (Array<int> * clusters) const;

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
    { MultAdd ( s, x, y ); }

    /*
      y += s L * x
    */
    virtual void MultAdd1 (double s, const BaseVector & x, BaseVector & y) const;

    /*
      y += s (D + L^T) * x
    */
    virtual void MultAdd2 (double s, const BaseVector & x, BaseVector & y) const;



  };

#endif // parallel 

}

#endif
