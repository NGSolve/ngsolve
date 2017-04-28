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
  public:
    ElementByElementMatrix (int h, int ane, bool isymmetric=false);
    ElementByElementMatrix (int h, int w, int ane, bool isymmetric=false);
    ElementByElementMatrix (int h, int w, int ane, bool isymmetric, bool adisjointrows, bool adisjointcols);
    ElementByElementMatrix (int h, int ane, bool isymmetric, bool adisjointrows, bool adisjointcols)
      : ElementByElementMatrix(h, h, ane, isymmetric, adisjointrows, adisjointcols) {};
    ~ElementByElementMatrix();
    virtual bool IsComplex() const { return typeid(SCAL)==typeid(Complex); }

    void SetDisjointRows(bool newval){disjointrows=newval;}
    void SetDisjointCols(bool newval){disjointcols=newval;}
    virtual int VHeight() const { return height; }
    virtual int VWidth() const { return width; }

    virtual AutoVector CreateVector () const
    {
      return make_shared<VVector<double>> (height);
    }

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const;

    void AddElementMatrix (int elnr,
                           FlatArray<int> dnums1,
                           FlatArray<int> dnums2,
                           BareSliceMatrix<SCAL> elmat);
			   
    void AddCloneElementMatrix(int elnr,
                           const FlatArray<int> & dnums1,
			   const FlatArray<int> & dnums2,
			   int refelnr);

    virtual BaseVector & AsVector() 
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

    virtual ostream & Print (ostream & ost) const;

    virtual BaseBlockJacobiPrecond * 
    CreateBlockJacobiPrecond (Table<int> & blocks,
			      const BaseVector * constraint = 0, int * paralleloptions = 0) const;
//     { 
//       return new BlockJacobiPrecond_ElByEl<SCAL> (*this, blocks );
//     }

    using BaseMatrix::InverseMatrix;
    virtual shared_ptr<BaseMatrix> 
    InverseMatrix (BitArray * subset = 0) const;

    int GetNZE () const
    {
      int nze = 0;
      for (int i = 0; i < elmats.Size(); i++)
	if (!clone.Test(i))
	  nze += elmats[i].Height()*elmats[i].Width();
      return nze;
    }

  };  

}

#endif
