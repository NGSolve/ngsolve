#ifndef FILE_NGS_PARALLEL_MATRICES
#define FILE_NGS_PARALLEL_MATRICES

/* ************************************************************************/
/* File:   parallelmatrices.hpp                                           */
/* Author: Astrid Sinwel, Joachim Schoeberl                               */
/* Date:   2007,2011                                                      */
/* ************************************************************************/

namespace ngla
{

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

    virtual AutoVector CreateVector () const override;
  };


  class ParallelMatrix : public BaseMatrix
  {
    shared_ptr<BaseMatrix> mat;
    // const ParallelDofs & pardofs;

    shared_ptr<ParallelDofs> row_paralleldofs, col_paralleldofs;

  public:
    ParallelMatrix (shared_ptr<BaseMatrix> amat, shared_ptr<ParallelDofs> apardofs);
    // : mat(*amat), pardofs(*apardofs) 
    // {const_cast<BaseMatrix&>(mat).SetParallelDofs (apardofs);}

    ParallelMatrix (shared_ptr<BaseMatrix> amat, shared_ptr<ParallelDofs> arpardofs,
		    shared_ptr<ParallelDofs> acpardofs);
    
    virtual ~ParallelMatrix () override;
    virtual bool IsComplex() const override { return mat->IsComplex(); } 
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

    shared_ptr<ParallelDofs> GetRowParallelDofs () const { return row_paralleldofs; }
    shared_ptr<ParallelDofs> GetColParallelDofs () const { return col_paralleldofs; }

    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<BitArray> subset = 0) const override;
    template <typename TM>
    shared_ptr<BaseMatrix> InverseMatrixTM (shared_ptr<BitArray> subset = 0) const;

    virtual shared_ptr<BaseMatrix> InverseMatrix (shared_ptr<const Array<int>> clusters) const override;
    virtual INVERSETYPE SetInverseType ( INVERSETYPE ainversetype ) const override;
    virtual INVERSETYPE SetInverseType ( string ainversetype ) const override;
    virtual INVERSETYPE GetInverseType () const override;
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


  class FETIDP_Constraint_Matrix : public BaseMatrix
  {
  public:
    FETIDP_Constraint_Matrix (const Table<size_t> & dofs, const Table<int> & dps,
			      const Table<double> & vals, shared_ptr<ParallelDofs> pardofs);
			      
    virtual bool IsComplex() const override { return false; }
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    
    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;

    shared_ptr<ParallelDofs> GetRowParallelDofs () const { return GetParallelDofs(); }
    shared_ptr<ParallelDofs> GetColParallelDofs () const { return mu_paralleldofs; }

    virtual ostream & Print (ostream & ost) const override
    {
      ost << "FETI-DP constraints; matrix: " << endl;
      ost << *mat << endl;
      return ost;
    }

  protected:
    shared_ptr<SparseMatrix<double> > mat;
    shared_ptr<ParallelDofs> mu_paralleldofs;
  };
  
#endif
}

#endif
