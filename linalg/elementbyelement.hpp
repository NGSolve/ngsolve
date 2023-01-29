#ifndef FILE_NGS_ELEMENTBYELEMENT
#define FILE_NGS_ELEMENTBYELEMENT

/* ************************************************************************/
/* File:   elementbyelement.hpp                                           */
/* Author: Joachim Schoeberl                                              */
/* Date:   June 2010							  */
/* ************************************************************************/

/* 
   Element by element matrix
*/

namespace ngla
{


  template <class SCAL>
  class NGS_DLL_HEADER ElementByElementMatrix : public BaseMatrix
  {
    Array<FlatMatrix<SCAL> > elmats;
    Array<FlatArray<int> > rowdnums;
    Array<FlatArray<int> > coldnums;
    int height;
    int width;
    int ne;
    bool symmetric;
    bool disjointrows;
    bool disjointcols;
    BitArray clone;
    int max_row_size = 0;
    int max_col_size = 0;

    Array<int> allrow, allcol;
    Array<SCAL> allvalues;
  public:
    ElementByElementMatrix (int h, int ane, bool isymmetric=false);
    ElementByElementMatrix (int h, int w, int ane, bool isymmetric=false);
    ElementByElementMatrix (int h, int w, int ane, bool isymmetric, bool adisjointrows, bool adisjointcols);
    ElementByElementMatrix (int h, int ane, bool isymmetric, bool adisjointrows, bool adisjointcols)
      : ElementByElementMatrix(h, h, ane, isymmetric, adisjointrows, adisjointcols) {};

    // allocate all memory at once
    ElementByElementMatrix (size_t h, size_t w,
                            FlatArray<int> nrowi, FlatArray<int> ncoli,
                            bool isymmetric, bool adisjointrows, bool adisjointcols);
    
    ~ElementByElementMatrix();
    bool IsComplex() const override { return typeid(SCAL)==typeid(Complex); }

    void SetDisjointRows(bool newval){disjointrows=newval;}
    void SetDisjointCols(bool newval){disjointcols=newval;}
    int VHeight() const override { return height; }
    int VWidth() const override { return width; }

    AutoVector CreateRowVector () const override { return make_unique<VVector<double>> (width); } 
    AutoVector CreateColVector () const override { return make_unique<VVector<double>> (height); }

    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    void AddElementMatrix (int elnr,
                           FlatArray<int> dnums1,
                           FlatArray<int> dnums2,
                           BareSliceMatrix<SCAL> elmat);
			   
    void AddCloneElementMatrix(int elnr,
                           const FlatArray<int> & dnums1,
			   const FlatArray<int> & dnums2,
			   int refelnr);

    BaseVector & AsVector() override
    {
      throw Exception ("Cannot access ebe-matrix AsVector");
      // return *new VVector<double> (1);
    }

    const FlatMatrix<SCAL> GetElementMatrix( int elnum ) const
    {
      return elmats[elnum];
    }

    const FlatArray<int> GetElementRowDNums ( int elnum ) const
    {
      return rowdnums[elnum]; 
    }

    const FlatArray<int> GetElementColumnDNums ( int elnum ) const
    {
      if (symmetric)
	return rowdnums[elnum]; 
      else
	return coldnums[elnum]; 
    }

    ostream & Print (ostream & ost) const override;

    virtual BaseBlockJacobiPrecond * 
    CreateBlockJacobiPrecond (Table<int> & blocks,
			      const BaseVector * constraint = 0, int * paralleloptions = 0) const;
//     { 
//       return new BlockJacobiPrecond_ElByEl<SCAL> (*this, blocks );
//     }

    using BaseMatrix::InverseMatrix;
    virtual shared_ptr<BaseMatrix> 
    InverseMatrix (BitArray * subset = 0) const;

    size_t GetNZE () const
    {
      size_t nze = 0;
      for (size_t i = 0; i < elmats.Size(); i++)
	if (!clone.Test(i))
	  nze += elmats[i].Height()*elmats[i].Width();
      return nze;
    }
    
    size_t NZE () const override { return GetNZE(); }

  private:
    void InitMemoryTracing() const;
  };  


  class NGS_DLL_HEADER ConstantElementByElementMatrix : public BaseMatrix
  {
    size_t h, w;
    Matrix<> matrix;
    Table<int> col_dnums;
    Table<int> row_dnums;
    bool disjoint_rows, disjoint_cols;
    Table<int> row_coloring, col_coloring;
  public:
    ConstantElementByElementMatrix (size_t ah, size_t aw, Matrix<> amatrix,
                                    Table<int> acol_dnums, Table<int> arow_dnums);

    virtual int VHeight() const override { return h; }
    virtual int VWidth() const override { return w; }

    virtual BaseMatrix::OperatorInfo GetOperatorInfo () const override;
    
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    
    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;
    
    FlatMatrix<> GetMatrix() const { return matrix; }
    FlatTable<int> GetRowDNums() const { return row_dnums; }
    FlatTable<int> GetColDNums() const { return col_dnums; }

    FlatTable<int> GetRowColoring() const { return row_coloring; }    
    FlatTable<int> GetColColoring() const { return col_coloring; }    
  };

  class NGS_DLL_HEADER StructuredElementByElementMatrix : public BaseMatrix
  {
    size_t num;
    Matrix<> matrix;
  public:
    StructuredElementByElementMatrix (size_t anum, Matrix<> amatrix)
      : num(anum), matrix(amatrix) { ; } 
    
    virtual int VHeight() const override { return num*matrix.Height(); }
    virtual int VWidth() const override { return num*matrix.Width(); }
    
    virtual void Mult (const BaseVector & x, BaseVector & y) const override;
    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override;
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    
    virtual AutoVector CreateRowVector () const override
    {
      return make_unique<VVector<>> (num*matrix.Width());
    }
      
    virtual AutoVector CreateColVector () const override
    {
      return make_unique<VVector<>> (num*matrix.Height());
    }

    const Matrix<> & GetMatrix() const { return matrix; }
  };
  

}

#endif
