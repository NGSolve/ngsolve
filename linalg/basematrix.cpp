/*********************************************************************/
/* File:   basematrix.cpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

/* 
   base class in matrix hierarchy
*/

#include <la.hpp>

namespace ngla
{
  BaseMatrix :: BaseMatrix()
    : paralleldofs (NULL)
  {
    ;
  }

  BaseMatrix :: BaseMatrix (const ParallelDofs * aparalleldofs)
    : paralleldofs ( aparalleldofs )
  {     
    ;
  }
  
  BaseMatrix :: ~BaseMatrix ()
  {
    ;
  }

  int BaseMatrix :: VHeight() const
  {
    stringstream str;
    str << "Height not available, type = " << typeid(*this).name() << endl;    
    throw Exception (str.str());
  }
  
  int BaseMatrix :: VWidth() const
  {
    stringstream str;
    str << "Width not available, type = " << typeid(*this).name() << endl;    
    throw Exception (str.str());
  }

  BaseVector & BaseMatrix :: AsVector()
  {
    Exception e("AsVector called for basematrix, type = ");
    e.Append (typeid(*this).name());
    throw e;
  }

  const BaseVector & BaseMatrix :: AsVector() const
  {
    Exception e("AsVector called for basematrix, type = ");
    e.Append (typeid(*this).name());
    throw e;
  }
  
  void BaseMatrix :: SetZero()
  {
    AsVector() = 0;
  }

  ostream & BaseMatrix :: Print (ostream & ost) const
  {
    return (ost << "Print base-matrix" << endl);
  }

  void BaseMatrix :: MemoryUsage (Array<MemoryUsageStruct*> & mu) const
  { ; }

  /*
  const void * BaseMatrix :: Data() const
  {
    cout << "Data not available, type = " << typeid(*this).name() << endl;
    return 0;
  }

  void * BaseMatrix :: Data() 
  {
    cout << "Data not available, type = " << typeid(*this).name() << endl;
    return 0;
  }
  */

  shared_ptr<BaseMatrix> BaseMatrix :: CreateMatrix () const
  {
    throw Exception ("BaseMatrix::CraeteMatrix called");
  }
  /*
  BaseMatrix * BaseMatrix :: CreateMatrix (const Array<int> & elsperrow) const
  {
    throw Exception ("BaseMatrix::CraeteMatrix called");
  }
  */
  AutoVector BaseMatrix :: CreateVector () const
  {
    cout << "BaseMatrix::CreateVector called" << endl;
    return shared_ptr<BaseVector>();
  }

  AutoVector BaseMatrix :: CreateRowVector () const
  {
    return CreateVector();
  }

  AutoVector BaseMatrix :: CreateColVector () const
  {
    return CreateVector();
  }




  
  void BaseMatrix :: Mult (const BaseVector & x, BaseVector & y) const
  {
    y = 0;
    MultAdd (1, x, y);
  }

  void BaseMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    //    cout << "Warning: BaseMatrix::MultAdd(double), this = " << typeid(*this).name() << endl;
    auto temp = y.CreateVector();
    Mult (x, *temp);
    y += s * *temp;
  }

  void BaseMatrix :: MultAdd (Complex s, const BaseVector & x, BaseVector & y) const 
  {
    stringstream err;
    err << "BaseMatrix::MultAdd (Complex) called, type = " 
	<< typeid(*this).name();
    throw Exception (err.str());
  }
  
  void BaseMatrix :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    cout << "warning: BaseMatrix::MultTransAdd(double) calls MultAdd, ";
    cout << "type = " << typeid(*this).name() << endl;
    MultAdd (s, x, y);
    return;

    stringstream err;
    err << "BaseMatrix::MultTransAdd (double) called, type = " 
	<< typeid(*this).name();
    throw Exception (err.str());
  }

  void BaseMatrix :: MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const
  {
    //    cout << "warning: BaseMatrix::MultTransAdd(complex) calls MultAdd" << endl;
    MultAdd (s, x, y);
    return;

    stringstream err;
    err << "BaseMatrix::MultTransAdd (Complex) called, type = " 
	<< typeid(*this).name();
    throw Exception (err.str());
  }

   // to split mat x vec for symmetric matrices
  void BaseMatrix :: MultAdd1 (double s, const BaseVector & x, BaseVector & y,
			       const BitArray * ainner,
			       const Array<int> * acluster) const
  {
    MultAdd (s, x, y);
  }

  void BaseMatrix :: MultAdd2 (double s, const BaseVector & x, BaseVector & y,
			       const BitArray * ainner,
			       const Array<int> * acluster) const
  {
    ;
  }
  



  shared_ptr<BaseMatrix> BaseMatrix :: InverseMatrix (const BitArray * subset) const 
  {
    cerr << "BaseMatrix::InverseMatrix not available" << endl;
    return NULL;
  }
  
  shared_ptr<BaseMatrix> BaseMatrix :: InverseMatrix (const Array<int> * clusters) const
  {
    cerr << "BaseMatrix::InverseMatrix not available" << endl;
    return NULL;
  }
  
  INVERSETYPE BaseMatrix :: SetInverseType ( INVERSETYPE ainversetype ) const
  {
    cerr << "BaseMatrix::SetInverseType not available" << endl;
    return SPARSECHOLESKY;
  }
  
  INVERSETYPE BaseMatrix :: SetInverseType ( string ainversetype ) const
  {
    cerr << "BaseMatrix::SetInverseType not available" << endl;
    return SPARSECHOLESKY;
  }
  
  INVERSETYPE BaseMatrix :: GetInverseType () const
  {
    cerr << "BaseMatrix::GetInverseType not available" << endl;
    return SPARSECHOLESKY;
  }

  void BaseMatrix :: DoArchive (Archive & ar)
  {
    ;
  }


  template<>
  S_BaseMatrix<double> :: S_BaseMatrix () 
  { ; }


  template<>
  S_BaseMatrix<double> :: ~S_BaseMatrix () 
  { ; }

  S_BaseMatrix<Complex> :: S_BaseMatrix () 
  { ; }

  S_BaseMatrix<Complex> :: ~S_BaseMatrix () 
  { ; }


  void S_BaseMatrix<Complex> :: 
  MultAdd (double s, const BaseVector & x, BaseVector & y) const 
  {
    MultAdd (Complex(s), x, y);
  }

  void S_BaseMatrix<Complex> :: 
  MultAdd (Complex s, const BaseVector & x, BaseVector & y) const 
  {
    stringstream err;
    err << "S_BaseMatrix<Complex>::MultAdd (Complex) called, type = " 
	<< typeid(*this).name();
    throw Exception (err.str());
  }
  
  void S_BaseMatrix<Complex> :: 
  MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    MultTransAdd (Complex(s), x, y);
  }

  void S_BaseMatrix<Complex> :: 
  MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const
  {
    stringstream err;
    err << "S_BaseMatrix<Complex>::MultTransAdd (Complex) called, type = " 
	<< typeid(*this).name();
    throw Exception (err.str());
  }


  void VMatVecExpr :: CheckSize (BaseVector & dest_vec) const
  {
    if (m.Height() != dest_vec.Size() || m.Width() != x.Size())
      throw Exception (ToString ("matrix-vector: size does not fit\n") +
                       "matrix-type = " + typeid(m).name() + "\n" +
                       "matrix:     " + ToString(m.Height()) + " x " + ToString(m.Width()) + "\n"
                       "vector in : " + ToString(x.Size()) + "\n"
                       "vector res: " + ToString(dest_vec.Size()));

  }



  string GetInverseName (INVERSETYPE type)
  {
    switch (type)
      {
      case PARDISO:         return "pardiso";
      case PARDISOSPD:      return "pardisospd";
      case SPARSECHOLESKY:  return "sparsecholesky";
      case SUPERLU:         return "superlu";
      case SUPERLU_DIST:    return "superlu_dist";
      case MUMPS:           return "mumps";
      case MASTERINVERSE:   return "masterinverse";
      case UMFPACK:         return "umfpack";
      }
    return "";
  }

}
