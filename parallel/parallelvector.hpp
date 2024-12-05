#ifndef FILE_NGS_PARALLELVECTOR
#define FILE_NGS_PARALLELVECTOR

/* ************************************************************************/
/* File:   parallelvector.hpp                                             */
/* Author: Astrid Sinwel, Joachim Schoeberl                               */
/* Date:   2007,2011                                                      */
/* ************************************************************************/

#include <vvector.hpp>
#include <multivector.hpp>
#include <paralleldofs.hpp>

namespace ngla
{
  
  class NGS_DLL_HEADER ParallelBaseVector : virtual public BaseVector
  {
  protected:
    mutable PARALLEL_STATUS status;
    shared_ptr<ParallelDofs> paralleldofs;    
    shared_ptr<BaseVector> local_vec;
    
    mutable NgMPI_Requests sreqs;
    mutable NgMPI_Requests rreqs;

  public:
    ParallelBaseVector ()
    { ; }

    template <typename T> 
    BaseVector & operator= (const VVecExpr<T> & v)
    {
      v.AssignTo (1.0, *this);
      return *this;
    }

    virtual PARALLEL_STATUS Status () const { return status; }

    virtual void SetStatus ( PARALLEL_STATUS astatus ) const  
    {
      status = astatus;
    }

    virtual PARALLEL_STATUS GetParallelStatus () const override { return Status(); }
    virtual void SetParallelStatus (PARALLEL_STATUS stat) const override { SetStatus (stat); }
    virtual optional<NgMPI_Comm> GetCommunicator() const override
    {
      if (paralleldofs)
        return paralleldofs->GetCommunicator();
      else
        return nullopt;
    }


    virtual shared_ptr<ParallelDofs> GetParallelDofs () const
    {
      return paralleldofs; 
    }
    
    virtual bool IsParallelVector () const
    {
      return (this->Status() != NOT_PARALLEL);
    }
    
    virtual BaseVector & SetScalar (double scal) override;
    virtual BaseVector & SetScalar (Complex scal) override;
    virtual void SetZero () override;
    
    virtual BaseVector & Set (double scal, const BaseVector & v) override;
    virtual BaseVector & Set (Complex scal, const BaseVector & v) override;

    virtual BaseVector & Add (double scal, const BaseVector & v) override;
    virtual BaseVector & Add (Complex scal, const BaseVector & v) override;

    void PrintStatus ( ostream & ost ) const;

    virtual shared_ptr<BaseVector> GetLocalVector () const override
    { return local_vec; }
    
    virtual void Cumulate () const override; 
    
    virtual void Distribute() const override = 0;
    // { cerr << "ERROR -- Distribute called for BaseVector, is not parallel" << endl; }
    
    virtual NgMPI_Request ISend ( int dest ) const;
    // virtual void Send ( int dest ) const;
    
    virtual NgMPI_Request IRecvVec ( int dest ) = 0;
    // { cerr << "ERROR -- IRecvVec called for BaseVector, is not parallel" << endl; }

    // virtual void RecvVec ( int dest )
    // { cerr << "ERROR -- IRecvVec called for BaseVector, is not parallel" << endl; }
    
    virtual void AddRecvValues( int sender ) = 0;
    // { cerr << "ERROR -- AddRecvValues called for BaseVector, is not parallel" << endl; }

    virtual void SetParallelDofs (shared_ptr<ParallelDofs> aparalleldofs) = 0;
    // const Array<int> * procs = 0) = 0;
  };




  inline ParallelBaseVector * dynamic_cast_ParallelBaseVector (BaseVector * x)
  {
    // cout << "my dynamic * cast" << endl;
    AutoVector * ax = dynamic_cast<AutoVector*> (x);
    if (ax)
      return dynamic_cast<ParallelBaseVector*> (&**ax);
    return dynamic_cast<ParallelBaseVector*> (x);
  }

  inline const ParallelBaseVector * dynamic_cast_ParallelBaseVector (const BaseVector * x)
  {
    // cout << "my dynamic const * cast" << endl;
    const AutoVector * ax = dynamic_cast<const AutoVector*> (x);
    if (ax)
      { 
        // cout << "is an autovector" << endl; 
        return dynamic_cast<const ParallelBaseVector*> (&**ax);
      }
    return dynamic_cast<const ParallelBaseVector*> (x);
  }

  inline ParallelBaseVector & dynamic_cast_ParallelBaseVector (BaseVector & x)
  {
    // cout << "my dynamic cast" << endl;
    AutoVector * ax = dynamic_cast<AutoVector*> (&x);
    if (ax)
      return dynamic_cast<ParallelBaseVector&> (**ax);
    return dynamic_cast<ParallelBaseVector&> (x);
  }

  inline const ParallelBaseVector & dynamic_cast_ParallelBaseVector (const BaseVector & x)
  {
    // cout << "my dynamic cast" << endl;
    const AutoVector * ax = dynamic_cast<const AutoVector*> (&x);
    if (ax)
      return dynamic_cast<const ParallelBaseVector&> (**ax);
    return dynamic_cast<const ParallelBaseVector&> (x);
  }






  
  template <class SCAL>
  class NGS_DLL_HEADER S_ParallelBaseVector 
    : virtual public S_BaseVector<SCAL>, 
      virtual public ParallelBaseVector
  {
  protected:
    virtual SCAL InnerProduct (const BaseVector & v2, bool conjugate = false) const;
    virtual BaseVector & SetScalar (double scal)
    { return ParallelBaseVector::SetScalar(scal); }
  };


  template <class SCAL>
  class NGS_DLL_HEADER S_ParallelBaseVectorPtr
    : virtual public S_BaseVectorPtr<SCAL>, 
      virtual public S_ParallelBaseVector<SCAL>
  {
  protected:
    typedef SCAL TSCAL;
    using ParallelBaseVector :: status;
    using ParallelBaseVector :: paralleldofs;

    using ParallelBaseVector :: sreqs;
    using ParallelBaseVector :: rreqs;

    Table<SCAL> recvvalues;

    using S_BaseVectorPtr<TSCAL> :: pdata;
    using ParallelBaseVector :: local_vec;

  public:
    S_ParallelBaseVectorPtr (int as, int aes, void * adata, shared_ptr<ParallelDofs> apd, PARALLEL_STATUS stat) throw();
    S_ParallelBaseVectorPtr (int as, int aes, shared_ptr<ParallelDofs> apd, PARALLEL_STATUS stat) throw();

    virtual ~S_ParallelBaseVectorPtr ();
    virtual void SetParallelDofs (shared_ptr<ParallelDofs> aparalleldofs) override; 

    virtual void Distribute() const override;
    virtual ostream & Print (ostream & ost) const override;
    virtual AutoVector Range (T_Range<size_t> range) const override;
    virtual AutoVector Range (DofRange range) const override;
    
    virtual NgMPI_Request IRecvVec ( int dest ) override;
    // virtual void  RecvVec ( int dest );
    virtual void AddRecvValues( int sender ) override;
    virtual AutoVector CreateVector () const override;
    virtual unique_ptr<MultiVector> CreateMultiVector (size_t cnt) const override;
    
    virtual double L2Norm () const override;
  };
 



  template <typename T = double>
  class ParallelVVector : // public VVector<T>, 
			  public S_ParallelBaseVectorPtr<typename mat_traits<T>::TSCAL>
  {
    typedef typename mat_traits<T>::TSCAL TSCAL;
    enum { ES = sizeof(T) / sizeof(TSCAL) };

    using S_BaseVectorPtr<TSCAL> :: pdata;
    using ParallelBaseVector :: local_vec;

  public:
    [[deprecated("too much info, use ParallelVVector(pardofs) instead")]]
    explicit ParallelVVector (int as, shared_ptr<ParallelDofs> aparalleldofs,
			      PARALLEL_STATUS astatus = CUMULATED)
    : S_BaseVectorPtr<TSCAL> (as, ES), // VVector<T> (as), 
	S_ParallelBaseVectorPtr<TSCAL> (as, ES, aparalleldofs, astatus)
    { local_vec = make_shared<VFlatVector<T>>(as, (T*)pdata); }

    explicit ParallelVVector (shared_ptr<ParallelDofs> aparalleldofs,
			      PARALLEL_STATUS astatus = CUMULATED)
      : S_BaseVectorPtr<TSCAL> (aparalleldofs->GetNDofLocal(), ES), 
      // VVector<T> (aparalleldofs->GetNDofLocal()), 
	S_ParallelBaseVectorPtr<TSCAL> (aparalleldofs->GetNDofLocal(), ES, aparalleldofs, astatus)
    { local_vec = make_shared<VFlatVector<T>>(aparalleldofs->GetNDofLocal(), (T*)pdata); }


    virtual ~ParallelVVector() throw()
    { ; }
  };


  template <typename T = double>
  class ParallelVFlatVector : public VFlatVector<T>,
			      public S_ParallelBaseVectorPtr<typename mat_traits<T>::TSCAL>
  {
    typedef typename mat_traits<T>::TSCAL TSCAL;
    enum { ES = sizeof(T) / sizeof(TSCAL) };

    using S_BaseVectorPtr<TSCAL> :: pdata;
    using ParallelBaseVector :: local_vec;

  public:
    explicit ParallelVFlatVector (int as, T * adata, 
				  shared_ptr<ParallelDofs> aparalleldofs, 
				  PARALLEL_STATUS astatus)
    : S_BaseVectorPtr<TSCAL> (as, ES, adata),
      VFlatVector<T> (as, adata),
      S_ParallelBaseVectorPtr<TSCAL> (as, ES, aparalleldofs, astatus)
    { local_vec = make_shared<VFlatVector<T>>(aparalleldofs->GetNDofLocal(), (T*)pdata); }

    explicit ParallelVFlatVector ()
    : S_BaseVectorPtr<TSCAL> (0, ES, NULL),
      S_ParallelBaseVectorPtr<TSCAL> (0, ES, NULL)
    { local_vec = make_shared<VFlatVector<T>>(0, NULL); }
      
    virtual ~ParallelVFlatVector() throw()
    { ; }
  };



  extern AutoVector CreateParallelVector (shared_ptr<ParallelDofs> pardofs, PARALLEL_STATUS status);
  
}

#endif
