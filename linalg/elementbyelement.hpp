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
  class ElementByElementMatrix : public BaseMatrix
  {
    Array<FlatMatrix<SCAL> > elmats;
    Array<FlatArray<int> > dnums;
    int height;
    int ne;

  public:
    ElementByElementMatrix (int h, int ane);

    virtual int VHeight() const { return height; }
    virtual int VWidth() const { return height; }

    virtual BaseVector * CreateVector () const
    {
      return new VVector<double> (height);
    }

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const;

    void AddElementMatrix (int elnr,
                           const Array<int> & dnums1,
			   const Array<int> & dnums2,
			   const FlatMatrix<SCAL> & elmat);
    

    virtual BaseVector & AsVector() 
    {
      return *new VVector<double> (1);
    }

    const FlatMatrix<SCAL> GetElementMatrix( int elnum ) const
    {
      return elmats[elnum];
    }

    const FlatArray<int> GetElementDNums ( int elnum ) const
    {
      return dnums[elnum]; 
    }

    virtual ostream & Print (ostream & ost) const
    {
      ost << "Element-by-Element Matrix:" << endl;
      return ost;
    }

    virtual BaseBlockJacobiPrecond * 
    CreateBlockJacobiPrecond (Table<int> & blocks,
			      const BaseVector * constraint = 0, int * paralleloptions = 0) const;
//     { 
//       return new BlockJacobiPrecond_ElByEl<SCAL> (*this, blocks );
//     }


    virtual BaseMatrix * 
    InverseMatrix (BitArray * subset = 0) const;

  };  

}

#endif
