/*********************************************************************/
/* File:   elementbyelement.cpp                                      */
/* Author: Joachim Schoeberl                                         */
/* Date:   June 2010                                                 */
/*********************************************************************/

/* 
   Element by element matrix
*/

#include <la.hpp>
namespace ngla
{
  using namespace ngla;

  
  template <class SCAL> class ElementByElementMatrix;

  
  template <class SCAL>
  ElementByElementMatrix<SCAL> :: ElementByElementMatrix (int h, int ane) 
  {
    height = h; 
    ne = ane; 
    elmats.SetSize(ne);
    dnums.SetSize(ne);
    for (int i = 0; i < ne; i++)
      {
        elmats[i].AssignMemory (0, 0, NULL);
        dnums[i] = FlatArray<int> (0, NULL);
      }
  }
  
  
  template <class SCAL>
  void ElementByElementMatrix<SCAL> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static int timer = NgProfiler::CreateTimer ("EBE-matrix::MultAdd");
    NgProfiler::RegionTimer reg (timer);

    int maxs = 0;
    for (int i = 0; i < dnums.Size(); i++)
      maxs = max2 (maxs, dnums[i].Size());

    ArrayMem<SCAL, 100> mem1(maxs), mem2(maxs);
      
    FlatVector<SCAL> vx = dynamic_cast<const S_BaseVector<SCAL> & >(x).FVScal();
    FlatVector<SCAL> vy = dynamic_cast<S_BaseVector<SCAL> & >(y).FVScal();

    for (int i = 0; i < dnums.Size(); i++)
      {
        FlatArray<int> di (dnums[i]);
        FlatVector<SCAL> hv1(di.Size(), &mem1[0]);
        FlatVector<SCAL> hv2(di.Size(), &mem2[0]);
	  
        for (int j = 0; j < di.Size(); j++)
          hv1(j) = vx (di[j]);

        hv2 = elmats[i] * hv1;
        hv2 *= s;

        for (int j = 0; j < dnums[i].Size(); j++)
          vy (di[j]) += hv2[j];
      }
  }

  template <class SCAL>
  BaseMatrix *  ElementByElementMatrix<SCAL> :: InverseMatrix ( BitArray * subset ) const
  {
    cout << "wird das tatsaechlich verwendet ???" << endl;
    ElementByElementMatrix<SCAL> * invmat = new ElementByElementMatrix<SCAL> (height, ne);

    int maxs = 0;
    for (int i = 0; i < dnums.Size(); i++)
      maxs = max2 (maxs, dnums[i].Size());

    LocalHeap lh (maxs*maxs*sizeof(SCAL)+100);

    for ( int i = 0; i < dnums.Size(); i++ )
      {
        int nd = dnums[i] . Size();
        FlatMatrix<SCAL> mat(nd, nd, lh);
        Array<int> dnumsarray(nd);
        for ( int j = 0; j < nd; j++ )
          dnumsarray[j] = dnums[i][j];
        mat = elmats[i];

        LapackInverse(mat);

        invmat -> AddElementMatrix(i, dnumsarray, dnumsarray, mat);
        lh.CleanUp();
      }
  
    return invmat;
  }
  
  template <class SCAL>
  void ElementByElementMatrix<SCAL> :: AddElementMatrix (int elnr,
                                                         const Array<int> & dnums1,
                                                         const Array<int> & dnums2,
                                                         const FlatMatrix<SCAL> & elmat)
  {
    ArrayMem<int,50> used;
    for (int i = 0; i < dnums1.Size(); i++)
      if (dnums1[i] >= 0) used.Append(i);

    int s = used.Size();

    FlatMatrix<SCAL> mat (s, new SCAL[s*s]);
    for (int i = 0; i < s; i++)
      for (int j = 0; j < s; j++)
        mat(i,j) = elmat(used[i], used[j]);

    FlatArray<int> dn(s, new int[s]);
    for (int i = 0; i < s; i++)
      dn[i] = dnums1[used[i]];

    if (elnr < elmats.Size())
      {
        dnums[elnr] = dn;
        elmats[elnr].AssignMemory (s, s, &mat(0,0));
      }
  }


  template <class SCAL>
  BaseBlockJacobiPrecond * ElementByElementMatrix<SCAL> :: 
  CreateBlockJacobiPrecond (Table<int> & blocks,
                            const BaseVector * constraint, int * paralleloptions) const
  { 
    return 0;//new helmholtz_exp_cpp::BlockJacobiPrecond_ElByEl<SCAL> (*this, blocks );
  }

  template class ElementByElementMatrix<double>;
  template class ElementByElementMatrix<Complex>;  
  
}
