#ifndef FILE_NGS_PARALLEL_SPARSEMATRIX
#define FILE_NGS_PARALLEL_SPARSEMATRIX

/* ************************************************************************/
/* File:   parallelsparsematrix.hpp                                       */
/* Author: Astrid Sinwel, Joachim Schoeberl                               */
/* Date:   2007,2011                                                      */
/* ************************************************************************/

namespace ngla
{

#ifdef PARALLEL

  class ParallelBaseMatrix : virtual public BaseMatrix
  {
  protected:
    const ngparallel::ParallelDofs * paralleldofs;

  public:
    ///
    ParallelBaseMatrix () : paralleldofs(NULL) { ; }

    ///
    ParallelBaseMatrix ( const ParallelBaseMatrix & amat )
      : paralleldofs ( amat.GetParallelDofs() ) { ; }

    ///
    ParallelBaseMatrix ( const ngparallel::ParallelDofs * aparalleldofs )
      : paralleldofs ( aparalleldofs ) { ; }

    /// destructor
    virtual ~ParallelBaseMatrix ();


    virtual BaseMatrix * ConsistentMat () 
    { cerr << "ERROR -- ParallelBaseMatrix::ConsistentMat() called" << endl; return 0; }

    virtual const BaseMatrix * ConsistentMat () const  
    { cerr << "ERROR -- ParallelBaseMatrix::ConsistentMat() called" << endl; return 0; }

    virtual void SetConsistentMat ( BaseMatrix * aconsistentmat )
    { cerr << "ERROR -- ParallelBaseMatrix::SetConsistentMat called" << endl; }

    virtual void AllocateConsistentMat ()
    { cerr << "ERROR -- ParallelBaseMatrix::AllocateConsistentMat called" << endl; }

    virtual void  AllocateConsistentMat ( const class MatrixGraph & graph )
    { cerr << "ERROR -- ParallelBaseMatrix::AllocateConsistentMat(const MatrixGraph&) called" << endl; }

    virtual void CalcConsistentMat (LocalHeap & lh) 
    { cerr << "ERROR -- ParallelBaseMatrix::CalcConsistentMat called" << endl; }

    virtual void SetParallelDofs ( ngparallel::ParallelDofs * aparalleldofs )
    { paralleldofs = aparalleldofs; }


    virtual const ngparallel::ParallelDofs * GetParallelDofs ( ) const { return paralleldofs; }

    virtual bool IsParallelMatrix() const
    {
      (*testout) << "PARALLELDOFS " <<  paralleldofs << endl;
      if ( paralleldofs ) return true;
      else return false;
    }

  };














  template<class TM, 
	   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class ParallelSparseMatrix;

  template<class TM, 
	   class TV = typename mat_traits<TM>::TV_ROW>
  class ParallelSparseMatrixSymmetric;




  template<class TM, 
	   class TV_ROW = typename mat_traits<TM>::TV_ROW, 
	   class TV_COL = typename mat_traits<TM>::TV_COL>
  class ParallelBlockJacobiPrecond;

  template<class TM, 
	   class TV = typename mat_traits<TM>::TV_ROW>
  class ParallelBlockJacobiPrecondSymmetric;




  // ***************************************************************

  //               PARALLEL MATRICES

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
			       const ngcomp::Preconditioner * acoarsegridprecond = 0, bool parallel = 1,
			       const BitArray * freedofs = NULL) const;
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

    virtual BaseMatrix * InverseMatrix (const BitArray * subset = 0) const;

    virtual BaseMatrix * InverseMatrix (const Array<int> * clusters) const;

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
	SparseMatrixSymmetricTM<TM> (as, max_elsperrow),
	SparseMatrix<TM,TV,TV> (as, max_elsperrow),
	SparseMatrixSymmetric<TM,TV> (as, max_elsperrow), 
	ParallelSparseMatrix<TM,TV,TV> (as, max_elsperrow)
    {     consistentmat = 0; }
  
    ParallelSparseMatrixSymmetric (const Array<int> & elsperrow)
      : SparseMatrixTM<TM> (elsperrow) , 
	SparseMatrixSymmetricTM<TM> (elsperrow),
	SparseMatrix<TM,TV,TV> (elsperrow),
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
	SparseMatrixSymmetricTM<TM> (as, max_elsperrow),
	SparseMatrix<TM,TV,TV> (as, max_elsperrow),
	SparseMatrixSymmetric<TM,TV> (as, max_elsperrow) ,
	ParallelSparseMatrix<TM,TV,TV> (as, max_elsperrow)
    {     consistentmat = aconsistentmat; }
  
    ParallelSparseMatrixSymmetric (const Array<int> & elsperrow, SparseMatrixSymmetric<TM,TV> * aconsistentmat )
      : SparseMatrixTM<TM> (elsperrow) , 
	SparseMatrixSymmetricTM<TM> (elsperrow),
	SparseMatrix<TM,TV,TV> (elsperrow),
	SparseMatrixSymmetric<TM,TV> (elsperrow),
	ParallelSparseMatrix<TM,TV,TV> (elsperrow)
    {     consistentmat = aconsistentmat; }

    ParallelSparseMatrixSymmetric (const MatrixGraph & agraph, bool stealgraph, SparseMatrixSymmetric<TM,TV> * aconsistentmat )
      : SparseMatrixTM<TM> (agraph, stealgraph) , 
	SparseMatrixSymmetricTM<TM> (agraph, stealgraph),
	SparseMatrix<TM,TV,TV> (agraph, stealgraph),
	SparseMatrixSymmetric<TM,TV> (agraph, stealgraph),
	ParallelSparseMatrix<TM,TV,TV> (agraph, stealgraph)
    {     consistentmat = aconsistentmat; }

    ParallelSparseMatrixSymmetric (const MatrixGraph & agraph, bool stealgraph, 
				   const ngparallel::ParallelDofs * aparalleldofs,
				   SparseMatrixSymmetric<TM,TV> * aconsistentmat )
      : SparseMatrixTM<TM> (agraph, stealgraph) , 
	SparseMatrixSymmetricTM<TM> (agraph, stealgraph),
	SparseMatrix<TM,TV,TV> (agraph, stealgraph),
	SparseMatrixSymmetric<TM,TV> (agraph, stealgraph),
	ParallelSparseMatrix<TM,TV,TV> (agraph, stealgraph,aparalleldofs)
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
			       const ngcomp::Preconditioner * acoarsegridprecond = 0, bool parallel  = 1,
			       const BitArray * freedofs = NULL) const;
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

    virtual BaseMatrix * InverseMatrix (const BitArray * subset = 0) const;

    virtual BaseMatrix * InverseMatrix (const Array<int> * clusters) const;

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











  class ParallelBaseBlockJacobiPrecond : virtual public BaseBlockJacobiPrecond,
					 virtual public ParallelBaseMatrix
  {
  protected:
    bool usecoarsegrid;
    const BaseMatrix * coarsegridprecond;
  
    Array<int> color;
    int ncolors;

  public:
    /// 
    ParallelBaseBlockJacobiPrecond (Table<int> & ablocktable, 
				    const ngparallel::ParallelDofs * aparalleldofs,
				    const ngcomp::Preconditioner * acoarsegridprecond = 0,
				    const BitArray * freedofs  = NULL);

    /// deletes the inverse matrix coarsegridprecond
    virtual ~ParallelBaseBlockJacobiPrecond ();

    virtual void ColorBlocks ();
  
    virtual void  MarkLockedDofs ( const int block, BitArray & lockeddofs ); 

    virtual bool IsFreeBlock ( const int block, const BitArray & lockeddofs ) const ;

    virtual bool IsLocalBlock ( const int block ) const;

    virtual void MyMPI_SendBitarrays ( const Array< BitArray* > & dofcolors, const int dest ) const ;

    virtual void MyMPI_RecvBitarrays ( Array< BitArray* > & dofcolors, const int dest ) const;

  };


#endif





#ifdef PARALLEL
  /**
     A parallel block-Jacobi preconditioner.
  */
  template <class TM, class TV_ROW, class TV_COL>
  class ParallelBlockJacobiPrecond : virtual public ParallelBaseBlockJacobiPrecond,
				     virtual public S_BaseMatrix<typename mat_traits<TM>::TSCAL>
  {
  protected:
    /// a reference to the matrix
    const ParallelSparseMatrix<TM,TV_ROW,TV_COL> & mat;
    /// inverses of the small blocks
    Array<Matrix<TM>*> invdiag;

  public:
    // typedef typename mat_traits<TM>::TV_ROW TVX;
    typedef TV_ROW TVX;
    typedef typename mat_traits<TM>::TSCAL TSCAL;

    ///
    ParallelBlockJacobiPrecond (const ParallelSparseMatrix<TM,TV_ROW,TV_COL> & amat,
				Table<int> & ablocktable, const ngcomp::Preconditioner * acoarsegridprecond,
				const BitArray * freedofs = NULL);
    ///
    virtual ~ParallelBlockJacobiPrecond ();

    int Height() const { return mat.Height(); }
    virtual int VHeight() const { return mat.Height(); }
    int Width() const { return mat.Width(); }
    virtual int VWidth() const { return mat.Width(); }

    ///
    virtual BaseVector * CreateVector () const 
    {
      return mat.CreateVector();
    }


    ///
    virtual void MultAdd (TSCAL s, const BaseVector & x, BaseVector & y) const ;


    ///
    virtual void MultTransAdd (TSCAL s, const BaseVector & x, BaseVector & y) const; 


    ///
    virtual void GSSmooth (BaseVector & x, const BaseVector & b,
			   int steps = 1) const ;

  
    virtual void GSSmoothResiduum (BaseVector & x, const BaseVector & b,
				   BaseVector & res, int steps = 1) const; 
  
    ///
    virtual void GSSmoothBack (BaseVector & x, const BaseVector & b,
			       int steps = 1) const ;
 
    ///
    virtual void GSSmoothNumbering (BaseVector & x, const BaseVector & b,
				    const Array<int> & numbering, 
				    int forward = 1) const
    {
      ;
    }

    virtual void MemoryUsage (Array<MemoryUsageStruct*> & mu) const
    {
      int nels = 0;
      for (int i = 0; i < blocktable.Size(); i++)
	{
	  int bs = blocktable[i].Size();
	  nels += bs*bs;
	}
      mu.Append (new MemoryUsageStruct ("BlockJac", nels*sizeof(TM), blocktable.Size()));
    }

    virtual void InitCoarseType ( string act, const BitArray * freedofs);
  };

#endif




#ifdef SYMCHOLESKY


#ifdef PARALLEL
  ///
  template <class TM, class TV>
  class ParallelBlockJacobiPrecondSymmetric : 
    virtual public ParallelBaseBlockJacobiPrecond,
    virtual public S_BaseMatrix<typename mat_traits<TM>::TSCAL>
    // virtual public ParallelBaseMatrix
  {
  protected:
    const ParallelSparseMatrixSymmetric<TM,TV> & mat;
    Array<CholeskyFactors<TM>*> invdiag;

    const SparseMatrixSymmetric<TM,TV> * consistentmat;

  public:

    typedef TV TVX;
    typedef typename mat_traits<TM>::TSCAL TSCAL;

    ParallelBlockJacobiPrecondSymmetric (const ParallelSparseMatrixSymmetric<TM,TV> & amat,
					 Table<int> & ablocktable, const ngcomp::Preconditioner * acoarsegridprecond = 0,
					 const BitArray * freedofs);

    ParallelBlockJacobiPrecondSymmetric (const ParallelSparseMatrixSymmetric<TM,TV> & amat, 
					 const FlatVector<TVX> & constraint,
					 Table<int> & ablocktable, const ngcomp::Preconditioner * acoarsegridprecond = 0,
					 const BitArray * freedofs);
    ///
    virtual ~ParallelBlockJacobiPrecondSymmetric ();

    virtual void MultAdd (TSCAL s, const BaseVector & x, BaseVector & y) const;

    ///
    virtual void MultTransAdd (TSCAL s, const BaseVector & x, 
			       BaseVector & y) const 
    {
      MultAdd (s, x, y);
    }

    ///
    virtual BaseVector * CreateVector () const 
    {
      return mat.CreateVector();
    }


    int Height() const { return mat.Height(); }
    virtual int VHeight() const { return mat.Height(); }
    int Width() const { return mat.Width(); }
    virtual int VWidth() const { return mat.Width(); }


    virtual void GSSmooth (BaseVector & x, const BaseVector & b,
			   int steps = 1) const ;

  
    virtual void GSSmoothPartial (BaseVector & x, const BaseVector & b,
				  BaseVector & y) const ;

    ///
    virtual void GSSmoothResiduum (BaseVector & x, const BaseVector & b,
				   BaseVector & res, int steps = 1) const ;
  
    ///
    virtual void GSSmoothBack (BaseVector & x, const BaseVector & b,
			       int steps = 1) const ;



    virtual void GSSmoothBackPartial (BaseVector & x, const BaseVector & b,
				      BaseVector & y) const ;



    void SmoothBlock (int i, 
		      FlatVector<TVX> & x,
		      const FlatVector<TVX> & b,
		      FlatVector<TVX> & y) const
    {
      FlatArray<int> row = blocktable[i];

      int bs = row.Size();
      if (bs == 0) return;
    
      VectorMem<500,TVX> di (bs);
      VectorMem<500,TVX> wi (bs, &wimem[0]);

    
      // di = P_i (y - L x)
      for (int j = 0; j < bs; j++)
	di(j) = y(row[j]) - mat.RowTimesVectorNoDiag (row[j], x);

      invdiag[i]->Mult (di, wi);

      // x += P_i w
      // y -= (D L^t) P_i w
      for (int j = 0; j < bs; j++)
	{
	  x(row[j]) += wi(j);
	  mat.AddRowTransToVector (row[j], -wi(j), y);
	}
    }

  
    ///
    virtual void GSSmoothNumbering (BaseVector & x, const BaseVector & b,
				    const Array<int> & numbering, 
				    int forward = 1) const
    {
      ;
    }
  
  
    virtual void MemoryUsage (Array<MemoryUsageStruct*> & mu) const
    {
      int nels = 0;
      for (int i = 0; i < blocktable.Size(); i++)
	{
	  int bs = blocktable[i].Size();
	  nels += (bs * (bs+1))/2;
	}
      mu.Append (new MemoryUsageStruct ("BlockJacSym", nels*sizeof(TM), blocktable.Size()));
    }

    virtual void InitCoarseType ( string act, const BitArray * freedofs) ;

  };

#endif

#else


#ifdef PARALLEL

  ///
  template <class TM, class TV>
  class ParallelBlockJacobiPrecondSymmetric : 
    virtual public ParallelBaseBlockJacobiPrecond,
    virtual public S_BaseMatrix<typename mat_traits<TM>::TSCAL>
    // virtual public ParallelBaseMatrix
  {
  protected:
    const ParallelSparseMatrixSymmetric<TM,TV> & mat;

    // 
    enum { NBLOCKS = 20 };
    DynamicMem<int> blockstart, blocksize, blockbw;
    DynamicMem<TM> data[NBLOCKS];
    bool lowmem;

    //   bool directsolver;
    //   BaseMatrix * coarsegridprecond;
  public:
  
    typedef TV TVX;
    typedef typename mat_traits<TM>::TSCAL TSCAL;
  
    ///
    ParallelBlockJacobiPrecondSymmetric (const ParallelSparseMatrixSymmetric<TM,TV> & amat,
					 Table<int> & ablocktable, const ngcomp::Preconditioner * acoarsegridprecond = 0,
					 const BitArray * freedofs = NULL);
    ///
    ParallelBlockJacobiPrecondSymmetric (const ParallelSparseMatrixSymmetric<TM,TV> & amat, 
					 const FlatVector<TVX> & constraint,
					 Table<int> & ablocktable, const ngcomp::Preconditioner * acoarsegridprecond = 0,
					 const BitArray * freedofs = NULL);
    ///
    virtual ~ParallelBlockJacobiPrecondSymmetric ();
  
    FlatBandCholeskyFactors<TM> InvDiag (int i) const
    {
      return FlatBandCholeskyFactors<TM> (blocksize[i], 
					  blockbw[i], 
					  const_cast<TM*>(data[i%NBLOCKS]+blockstart[i]));
    }

    void ComputeBlockFactor (FlatArray<int> block, int bw, FlatBandCholeskyFactors<TM> & inv) const;
  
    void ComputeBlockFactorParallel(FlatArray<int> block, int bw, FlatBandCholeskyFactors<TM> & inv) const;

    virtual void MultAdd (TSCAL s, const BaseVector & x, BaseVector & y) const;

    ///
    virtual void MultTransAdd (TSCAL s, const BaseVector & x, BaseVector & y) const 
    {
      MultAdd (s, x, y);
    }


    ///
    virtual BaseVector * CreateVector () const 
    {
      return mat.CreateVector();
    }


    int Height() const { return mat.Height(); }
    virtual int VHeight() const { return mat.Height(); }
    int Width() const { return mat.Width(); }
    virtual int VWidth() const { return mat.Width(); }


    virtual void GSSmooth (BaseVector & x, const BaseVector & b,
			   int steps = 1) const ;

  
    virtual void GSSmoothPartial (BaseVector & x, const BaseVector & b,
				  BaseVector & y) const ;

    ///
    virtual void GSSmoothResiduum (BaseVector & x, const BaseVector & b,
				   BaseVector & res, int steps = 1) const ;
  
    ///
    virtual void GSSmoothBack (BaseVector & x, const BaseVector & b,
			       int steps = 1) const ;



    virtual void GSSmoothBackPartial (BaseVector & x, const BaseVector & b,
				      BaseVector & y) const ;


    /// x ... solution, y ... partial residual
    void SmoothBlock (int i, 
		      FlatVector<TVX> & x,
		      // const FlatVector<TVX> & b,
		      FlatVector<TVX> & y) const;

    ///
    virtual void GSSmoothNumbering (BaseVector & x, const BaseVector & b,
				    const Array<int> & numbering, 
				    int forward = 1) const
    {
      ;
    }

    virtual void InitCoarseType ( string act, const BitArray * freedofs) ;

  
  };
#endif // parallel

#endif










#ifdef PARALLEL

  template <typename TM>
  class MasterInverse : public ParallelBaseMatrix
  {
    BaseMatrix * inv;
    const BitArray * subset;
    DynamicTable<int> loc2glob;
    Array<int> select;
  public:
    MasterInverse (const SparseMatrixTM<TM> & mat, const BitArray * asubset, const ParallelDofs * pardofs);
    virtual ~MasterInverse ();
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
  };




  class ParallelMatrix : public BaseMatrix
  {
    const BaseMatrix & mat;
    const ParallelDofs & pardofs;
  public:
    ParallelMatrix (const BaseMatrix & amat, const ParallelDofs & apardofs)
      : mat(amat), pardofs(apardofs) { ; }

    virtual ~ParallelMatrix ();
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const;
  };






#endif






#endif




}
