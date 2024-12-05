#ifndef FILE_NGS_PARALLEL_MATRICES
#define FILE_NGS_PARALLEL_MATRICES

/* ************************************************************************/
/* File:   parallelmatrices.hpp                                           */
/* Author: Astrid Sinwel, Joachim Schoeberl                               */
/* Date:   2007,2011                                                      */
/* ************************************************************************/


#include <basematrix.hpp>
#include <sparsematrix.hpp>
#include <paralleldofs.hpp>

namespace ngla
{

  // enum with char type issue in pybind11:
  // https://github.com/pybind/pybind11/issues/1820
  enum PARALLEL_OP : uint8_t { D2D = 0,   // 00
                               D2C = 1,   // 01
                               C2D = 2,   // 10
                               C2C = 3 }; // 11

  inline ostream & operator<< (ostream & ost, PARALLEL_OP op)
  {
    switch (op)
      {
      case D2D: ost << "D2D"; break;
      case D2C: ost << "D2C"; break;
      case C2D: ost << "C2D"; break;
      case C2C: ost << "C2C"; break;
      default: ost << "undefined parallelop";
      }
    return ost;
  }

  
  inline PARALLEL_STATUS RowType (PARALLEL_OP op)
  {
    if (op == C2D || op == C2C)
      return CUMULATED;
    else
      return DISTRIBUTED;
  }

  inline PARALLEL_STATUS ColType (PARALLEL_OP op)
  {
    if (op == D2C || op == C2C)
      return CUMULATED;
    else
      return DISTRIBUTED;
  }

  inline PARALLEL_STATUS InvertType (PARALLEL_STATUS stat)
  {
    if (stat == NOT_PARALLEL) return NOT_PARALLEL;
    return (stat == CUMULATED) ? DISTRIBUTED : CUMULATED;
  }
  
  inline PARALLEL_OP ParallelOp (PARALLEL_STATUS stat_row, PARALLEL_STATUS stat_col)
  {
    if (stat_row == CUMULATED)
      return (stat_col == CUMULATED) ? C2C : C2D;
    else
      return (stat_col == CUMULATED) ? D2C : D2D;
  }

  
  class ParallelMatrix : public BaseMatrix
  {
    shared_ptr<BaseMatrix> mat;
    shared_ptr<ParallelDofs> row_paralleldofs, col_paralleldofs;
    PARALLEL_OP op;
    
  public:
    ParallelMatrix (shared_ptr<BaseMatrix> amat, shared_ptr<ParallelDofs> apardofs,
		    PARALLEL_OP op = C2D);
    // : mat(*amat), pardofs(*apardofs) 
    // {const_cast<BaseMatrix&>(mat).SetParallelDofs (apardofs);}

    ParallelMatrix (shared_ptr<BaseMatrix> amat, shared_ptr<ParallelDofs> arpardofs,
		    shared_ptr<ParallelDofs> acpardofs, PARALLEL_OP op = C2D);

    virtual ~ParallelMatrix () override;
    virtual bool IsComplex() const override { return mat->IsComplex(); }
    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override;
    
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override ;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    virtual BaseVector & AsVector() override { return mat->AsVector(); }
    virtual const BaseVector & AsVector() const override { return mat->AsVector(); }

    shared_ptr<BaseMatrix> GetMatrix() const { return mat; }
    virtual shared_ptr<BaseMatrix> CreateMatrix () const override;
    virtual AutoVector CreateVector () const override;
    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;

    virtual ostream & Print (ostream & ost) const override;

    virtual int VHeight() const override;
    virtual int VWidth() const override;

    // virtual const ParallelDofs * GetParallelDofs () const {return &pardofs;}
    shared_ptr<ParallelDofs> GetRowParallelDofs () const { return row_paralleldofs; }
    shared_ptr<ParallelDofs> GetColParallelDofs () const { return col_paralleldofs; }

    PARALLEL_OP GetOpType () const { return op; }
    virtual optional<NgMPI_Comm> GetCommunicator() const override
    {
      if (row_paralleldofs)
        return row_paralleldofs->GetCommunicator();
      else
        return nullopt;
    }

    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = 0) const override;
    template <typename TM>
    shared_ptr<BaseMatrix> InverseMatrixTM (shared_ptr<BitArray> subset = 0) const;

    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<const Array<int>> clusters) const override;
    virtual INVERSETYPE SetInverseType ( INVERSETYPE ainversetype ) const override;
    virtual INVERSETYPE SetInverseType ( string ainversetype ) const override;
    virtual INVERSETYPE GetInverseType () const override;

    virtual shared_ptr<BaseMatrix> DeleteZeroElements(double tol) const override
    {
      return make_shared<ParallelMatrix> (mat->DeleteZeroElements(tol),
                                          row_paralleldofs, col_paralleldofs, op);
    }
    
  };



  class CumulationOperator : public BaseMatrix
  {
    shared_ptr<ParallelDofs> pardofs;
    
  public:
    CumulationOperator (shared_ptr<ParallelDofs> apardofs)
      : pardofs(apardofs) { } 
    ~CumulationOperator() override;
    bool IsComplex() const override { return pardofs->IsComplex(); } 
    void Mult (const BaseVector & x, BaseVector & y) const override ;
    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override ;
    void MultTrans (const BaseVector & x, BaseVector & y) const override;
    void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    AutoVector CreateRowVector () const override;
    AutoVector CreateColVector () const override;

    int VHeight() const override;
    int VWidth() const override;

    ostream & Print (ostream & ost) const override;
    PARALLEL_OP GetOpType () const { return PARALLEL_OP::D2C; }
  };


  

  
#ifdef PARALLEL


  template <typename TM>
  class MasterInverse : public BaseMatrix
  {
    shared_ptr<BaseMatrix> inv;
    shared_ptr<BitArray> subset;
    DynamicTable<int> loc2glob;
    Array<int> select;
    string invtype;
    //shared_ptr<ParallelDofs> pardofs;
  public:
    MasterInverse (const SparseMatrixTM<TM> & mat, shared_ptr<BitArray> asubset, 
		   shared_ptr<ParallelDofs> apardofs);
    virtual ~MasterInverse () override;
    virtual bool IsComplex() const override { return inv->IsComplex(); } 
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;

    virtual int VHeight() const override { return paralleldofs->GetNDofLocal(); }
    virtual int VWidth() const override { return paralleldofs->GetNDofLocal(); }

    AutoVector CreateRowVector() const override;
    AutoVector CreateColVector() const override;
  };


  
  class FETI_Jump_Matrix : public BaseMatrix
  {
  public:
    FETI_Jump_Matrix (shared_ptr<ParallelDofs> pardofs, shared_ptr<ParallelDofs> au_paralleldofs = nullptr);
    
    virtual bool IsComplex() const override { return false; }
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    
    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;

    shared_ptr<ParallelDofs> GetRowParallelDofs () const { return u_paralleldofs; }
    shared_ptr<ParallelDofs> GetColParallelDofs () const { return jump_paralleldofs; }
    
    virtual int VHeight() const override { return paralleldofs->GetNDofLocal(); }
    virtual int VWidth()  const override { return jump_paralleldofs->GetNDofLocal(); }

  protected:

    shared_ptr<ParallelDofs> jump_paralleldofs;
    shared_ptr<ParallelDofs> u_paralleldofs;

  };

#endif
}

#endif
