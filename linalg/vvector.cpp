/*********************************************************************/
/* File:   vvector.cpp                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Aug. 2002                                             */
/*********************************************************************/

/* 
   Implementation of VVector
*/
#include <la.hpp>

namespace ngla
{
  template <typename T>
  VFlatVector<T> :: VFlatVector () throw()
  { 
    this->entrysize = sizeof(T) / sizeof(double);
  }

  template <typename T>
  VFlatVector<T> :: VFlatVector (int as, T * adata) throw()
    : data(adata)
  { 
    this->size = as;
    this->entrysize = sizeof(T) / sizeof(double);

  }



  template <typename T>
  VFlatVector<T> :: ~VFlatVector() throw()
  { 
    ;
  }


  template <typename T>
  VVector<T> :: VVector (int as) 
    : data(as)
  { 
    this->size = as;
    this->entrysize = sizeof(T) / sizeof(double);
    // data = new T[s];
    data.SetName ("VVector");
    // mem_total_alloc_vector += as * sizeof(T);


  }


  template <typename T>
  VVector<T> :: ~VVector() throw()
  { 
    ;
  }

  template <typename T>
  BaseVector * T_BaseVector<T> :: CreateVector ( const Array<int> * procs ) const
  {
    VVector<T> * vec = new VVector<T> (this->size);
    return vec;
  }
// #else 
//   template <typename T>
//   BaseVector * T_BaseVector<T> :: CreateVector ( const Array<int> * procs ) const
//   {
//     VVector<T> * vec;
//     ParallelVVector<T> * parvec;
//     const ParallelVVector<T> * thisparvec = dynamic_cast<const ParallelVVector<T> *> (this);
 
//     if ( ! thisparvec )
//       {
// 	vec = new VVector<T> (this->size);
// 	return vec;
//       }
//     else
//       {
// 	parvec = new ParallelVVector<T> (this->size);
// 	parvec->SetStatus( thisparvec->Status() );

// 	*testout << "procs... " << procs << endl;
// 	if( thisparvec->Status() != NOT_PARALLEL )
// 	  {
// 	    ParallelDofs * aparalleldofs = thisparvec-> GetParallelDofs();
// 	    parvec->SetParallelDofs ( aparalleldofs, procs );
// 	  }
// 	return parvec;
//       }
//   }
// #endif


  template <typename T>
  void VVector<T> :: SetSize(int as)
  {
    if (this->size == as) return;
    this->size = as;
    data.Alloc (this->size);
  }








  /*
    template <typename T>
    VVector<FlatVector<T> > :: VVector (int as, int abs) 
    : data(as*abs)
    { 
    this->size = as;
    this->entrysize = sizeof(T) / sizeof(double);
    bs = abs;

    data.SetName ("VVector");
    }

    template <typename T>
    VVector<FlatVector<T> > :: ~VVector() throw()
    { ; }

    template <typename T>
    void VVector<FlatVector<T> > :: SetSize(int as, int abs)
    {
    if (this->size == as && bs == abs) return;
    this->size = as;
    bs = abs;
    data.Alloc (this->size * bs);
    }

  */


  template < class T >
  ostream & T_BaseVector<T> :: Print (ostream & ost) const
  {
    ost << "addr = " << &FV()(0) << endl;
// #ifdef PARALLEL
//     const PARALLEL_STATUS status = this->Status();
//     if ( status == NOT_PARALLEL )
//       ost << "NOT PARALLEL" << endl;
//     else if ( status == DISTRIBUTED )
//       ost << "DISTRIBUTED" <<endl;
//     else if ( status == CUMULATED )
//       ost << "CUMULATED" << endl;
// #endif
    return (ost << FV() << endl);
  }







template class T_BaseVector<double>;
template class T_BaseVector<Complex>;
template class T_BaseVector<Vec<1,double> >;
template class T_BaseVector<Vec<1,Complex> >;
template class T_BaseVector<Vec<2,double> >;
template class T_BaseVector<Vec<2,Complex> >;
template class T_BaseVector<Vec<3,double> >;
template class T_BaseVector<Vec<3,Complex> >;
template class T_BaseVector<Vec<4,double> >;
template class T_BaseVector<Vec<4,Complex> >;

template class T_BaseVector<Vec<5,double> >;
template class T_BaseVector<Vec<5,Complex> >;
template class T_BaseVector<Vec<6,double> >;
template class T_BaseVector<Vec<6,Complex> >;
template class T_BaseVector<Vec<7,double> >;
template class T_BaseVector<Vec<7,Complex> >;
template class T_BaseVector<Vec<8,double> >;
template class T_BaseVector<Vec<8,Complex> >;
template class T_BaseVector<Vec<9,double> >;
template class T_BaseVector<Vec<9,Complex> >;

template class T_BaseVector<Vec<12,double> >;
template class T_BaseVector<Vec<12,Complex> >;
template class T_BaseVector<Vec<18,double> >;
template class T_BaseVector<Vec<18,Complex> >;
template class T_BaseVector<Vec<24,double> >;
template class T_BaseVector<Vec<24,Complex> >;

template class VFlatVector<double>;
template class VFlatVector<Complex>;
template class VFlatVector<Vec<1,double> >;
template class VFlatVector<Vec<1,Complex> >;
template class VFlatVector<Vec<2,double> >;
template class VFlatVector<Vec<2,Complex> >;
template class VFlatVector<Vec<3,double> >;
template class VFlatVector<Vec<3,Complex> >;
template class VFlatVector<Vec<4,double> >;
template class VFlatVector<Vec<4,Complex> >;

template class VFlatVector<Vec<5,double> >;
template class VFlatVector<Vec<5,Complex> >;
template class VFlatVector<Vec<6,double> >;
template class VFlatVector<Vec<6,Complex> >;
template class VFlatVector<Vec<7,double> >;
template class VFlatVector<Vec<7,Complex> >;
template class VFlatVector<Vec<8,double> >;
template class VFlatVector<Vec<8,Complex> >;

  /*
template class VFlatVector<Vec<9,double> >;
template class VFlatVector<Vec<9,Complex> >;
template class VFlatVector<Vec<12,double> >;
template class VFlatVector<Vec<12,Complex> >;
  */
  
  /*
template class VFlatVector<Vec<18,double> >;
template class VFlatVector<Vec<18,Complex> >;
template class VFlatVector<Vec<24,double> >;
template class VFlatVector<Vec<24,Complex> >;
  */

template class VVector<double>;
template class VVector<Complex>;
template class VVector<Vec<1,double> >;
template class VVector<Vec<1,Complex> >;
template class VVector<Vec<2,double> >;
template class VVector<Vec<2,Complex> >;
template class VVector<Vec<3,double> >;
template class VVector<Vec<3,Complex> >;
template class VVector<Vec<4,double> >;
template class VVector<Vec<4,Complex> >;
template class VVector<Vec<5,double> >;
template class VVector<Vec<5,Complex> >;
template class VVector<Vec<6,double> >;
template class VVector<Vec<6,Complex> >;

template class VVector<Vec<7,double> >;
template class VVector<Vec<7,Complex> >;
template class VVector<Vec<8,double> >;
template class VVector<Vec<8,Complex> >;
  /*
template class VVector<Vec<9,double> >;
template class VVector<Vec<9,Complex> >;
template class VVector<Vec<10,double> >;
template class VVector<Vec<10,Complex> >;
template class VVector<Vec<11,double> >;
template class VVector<Vec<11,Complex> >;
template class VVector<Vec<12,double> >;
template class VVector<Vec<12,Complex> >;
  */

  /*
template class VVector<Vec<13,double> >;
template class VVector<Vec<13,Complex> >;
template class VVector<Vec<14,double> >;
template class VVector<Vec<14,Complex> >;
template class VVector<Vec<15,double> >;
template class VVector<Vec<15,Complex> >;
template class VVector<Vec<16,double> >;
template class VVector<Vec<16,Complex> >;
template class VVector<Vec<17,double> >;
template class VVector<Vec<17,Complex> >;
template class VVector<Vec<18,double> >;
template class VVector<Vec<18,Complex> >;

template class VVector<Vec<19,double> >;
template class VVector<Vec<19,Complex> >;
template class VVector<Vec<20,double> >;
template class VVector<Vec<20,Complex> >;
template class VVector<Vec<21,double> >;
template class VVector<Vec<21,Complex> >;
template class VVector<Vec<22,double> >;
template class VVector<Vec<22,Complex> >;
template class VVector<Vec<23,double> >;
template class VVector<Vec<23,Complex> >;
template class VVector<Vec<24,double> >;
template class VVector<Vec<24,Complex> >;
  */
  
  // template class VFlatVector<FlatVector<double> >;
  // template class VVector<FlatVector<double> >;





/*
  template class T_BaseVector<Mat<1,1,double> >;
  template class T_BaseVector<Mat<1,1,Complex> >;
  template class T_BaseVector<Mat<2,2,double> >;
  template class T_BaseVector<Mat<2,2,Complex> >;
  template class T_BaseVector<Mat<3,3,double> >;
  template class T_BaseVector<Mat<3,3,Complex> >;
  template class T_BaseVector<Mat<4,4,double> >;
  template class T_BaseVector<Mat<4,4,Complex> >;

  template class T_BaseVector<Mat<5,5,double> >;
  template class T_BaseVector<Mat<5,5,Complex> >;
  template class T_BaseVector<Mat<6,6,double> >;
  template class T_BaseVector<Mat<6,6,Complex> >;
  template class T_BaseVector<Mat<7,7,double> >;
  template class T_BaseVector<Mat<7,7,Complex> >;
  template class T_BaseVector<Mat<8,8,double> >;
  template class T_BaseVector<Mat<8,8,Complex> >;


  template class VFlatVector<Mat<1,1,double> >;
  template class VFlatVector<Mat<1,1,Complex> >;
  template class VFlatVector<Mat<2,2,double> >;
  template class VFlatVector<Mat<2,2,Complex> >;
  template class VFlatVector<Mat<3,3,double> >;
  template class VFlatVector<Mat<3,3,Complex> >;
  template class VFlatVector<Mat<4,4,double> >;
  template class VFlatVector<Mat<4,4,Complex> >;

  template class VFlatVector<Mat<5,5,double> >;
  template class VFlatVector<Mat<5,5,Complex> >;
  template class VFlatVector<Mat<6,6,double> >;
  template class VFlatVector<Mat<6,6,Complex> >;
  template class VFlatVector<Mat<7,7,double> >;
  template class VFlatVector<Mat<7,7,Complex> >;
  template class VFlatVector<Mat<8,8,double> >;
  template class VFlatVector<Mat<8,8,Complex> >;



  template class VVector<Mat<1,1,double> >;
  template class VVector<Mat<1,1,Complex> >;
  template class VVector<Mat<2,2,double> >;
  template class VVector<Mat<2,2,Complex> >;
  template class VVector<Mat<3,3,double> >;
  template class VVector<Mat<3,3,Complex> >;
  template class VVector<Mat<4,4,double> >;
  template class VVector<Mat<4,4,Complex> >;

  template class VVector<Mat<5,5,double> >;
  template class VVector<Mat<5,5,Complex> >;
  template class VVector<Mat<6,6,double> >;
  template class VVector<Mat<6,6,Complex> >;
  template class VVector<Mat<7,7,double> >;
  template class VVector<Mat<7,7,Complex> >;
  template class VVector<Mat<8,8,double> >;
  template class VVector<Mat<8,8,Complex> >;
*/
}
