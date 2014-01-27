/*********************************************************************/
/* File:   diffop.cpp                                                */
/* Author: Start                                                     */
/* Date:   24. Nov. 2009                                             */
/*********************************************************************/


 
#include <fem.hpp>
#include "diffop_impl.hpp"



namespace ngfem
{


  void DifferentialOperator ::
  CalcMatrix (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatMatrix<double,ColMajor> mat, 
              LocalHeap & lh) const 
  {
    cerr << "DifferentialOperator::CalcMatrix called for base class, type = " 
         << typeid(*this).name()
         << endl;
  }

  void DifferentialOperator ::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationPoint & mip,
         FlatVector<double> x, 
         FlatVector<double> flux,
         LocalHeap & lh) const
  {
    cout << "called base class apply, type = " << typeid(*this).name() << endl;
    FlatMatrix<double,ColMajor> mat(Dim(), x.Size(), lh);
    CalcMatrix (fel, mip, mat, lh);
    flux = mat * x;
  }
  
  void DifferentialOperator ::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationPoint & mip,
         FlatVector<Complex> x, 
         FlatVector<Complex> flux,
         LocalHeap & lh) const
  {
    cout << "called base class apply, complex" << endl;
    FlatMatrix<double,ColMajor> mat(Dim(), x.Size(), lh);
    CalcMatrix (fel, mip, mat, lh);
    flux = mat * x;
  }
  
  void DifferentialOperator ::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationRule & mir,
         FlatVector<double> x, 
         FlatMatrix<double> flux,
         LocalHeap & lh) const
  {
    for (int i = 0; i < mir.Size(); i++)
      Apply (fel, mir[i], x, flux.Row(i), lh);
  }
  
  void DifferentialOperator ::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationRule & mir,
         FlatVector<Complex> x, 
         FlatMatrix<Complex> flux,
         LocalHeap & lh) const
  {
    for (int i = 0; i < mir.Size(); i++)
      Apply (fel, mir[i], x, flux.Row(i), lh);
  }
  
  
  void DifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<double> flux,
              FlatVector<double> x, 
              LocalHeap & lh) const 
  {
    cout << "called base class apply trans" << endl;
    FlatMatrix<double,ColMajor> mat(Dim(), x.Size(), lh);
    CalcMatrix (fel, mip, mat, lh);
    x = Trans(mat) * flux;
  }
  
  
  void DifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<Complex> flux,
              FlatVector<Complex> x, 
              LocalHeap & lh) const 
  {
    cout << "called base class apply trans, complex" << endl;
    FlatMatrix<double,ColMajor> mat(Dim(), x.Size(), lh);
    CalcMatrix (fel, mip, mat, lh);
    x = Trans(mat) * flux;
  }
  
  
  void DifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationRule & mir,
              FlatMatrix<double> flux,
              FlatVector<double> x, 
              LocalHeap & lh) const 
  {
    FlatVector<double> hx(x.Size(), lh);
    x = 0.0;
    for (int i = 0; i < mir.Size(); i++)
      {
        ApplyTrans (fel, mir[i], flux.Row(i), hx, lh);
        x += hx;
      }
  }
  
  void DifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationRule & mir,
              FlatMatrix<Complex> flux,
              FlatVector<Complex> x, 
              LocalHeap & lh) const 
  {
    FlatVector<Complex> hx(x.Size(), lh);
    x = 0.0;
    for (int i = 0; i < mir.Size(); i++)
      {
        ApplyTrans (fel, mir[i], flux.Row(i), hx, lh);
        x += hx;
      }
  }
  




  BlockDifferentialOperator ::
  ~BlockDifferentialOperator ()
  {
    delete &diffop;
  }


  void BlockDifferentialOperator ::
  CalcMatrix (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatMatrix<double,ColMajor> mat, 
              LocalHeap & lh) const 
  {
    HeapReset hr(lh);
    FlatMatrix<double,ColMajor> mat1(diffop.Dim(), fel.GetNDof(), lh);
    diffop.CalcMatrix (fel, mip, mat1, lh);
    mat = 0;
    
    if (comp == -1)
      for (int i = 0; i < mat1.Height(); i++)
        for (int j = 0; j < mat1.Width(); j++)
          for (int k = 0; k < dim; k++)
            mat(dim*i+k, dim*j+k) = mat1(i,j);
    else
      for (int i = 0; i < mat1.Height(); i++)
        for (int j = 0; j < mat1.Width(); j++)
          mat(dim*i+comp, dim*j+comp) = mat1(i,j);
  }
  

  void BlockDifferentialOperator ::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationPoint & mip,
         FlatVector<double> x, 
         FlatVector<double> flux,
         LocalHeap & lh) const
  {
    FlatVector<> hx(fel.GetNDof(), lh);
    FlatVector<> hflux(diffop.Dim(), lh);
    
    if (comp == -1)
      {
        for (int k = 0; k < dim; k++)
          {
            hx = x.Slice(k, dim);
            diffop.Apply(fel, mip, hx, hflux, lh);
            flux.Slice(k,dim) = hflux;
          }
      }
    else
      {
        hx = x.Slice(comp, dim);
        diffop.Apply(fel, mip, hx, hflux, lh);
        flux.Slice(comp,dim) = hflux;
      }
  }



  /*
  template class T_DifferentialOperator<DiffOpId<1> >;
  template class T_DifferentialOperator<DiffOpId<2> >;
  template class T_DifferentialOperator<DiffOpId<3> >;

  template class T_DifferentialOperator<DiffOpIdBoundary<1> >;
  template class T_DifferentialOperator<DiffOpIdBoundary<2> >;
  template class T_DifferentialOperator<DiffOpIdBoundary<3> >;


  template class T_DifferentialOperator<DiffOpGradient<1> >;
  template class T_DifferentialOperator<DiffOpGradient<2> >;
  template class T_DifferentialOperator<DiffOpGradient<3> >;

  template class T_DifferentialOperator<DiffOpIdEdge<2> >;
  template class T_DifferentialOperator<DiffOpIdEdge<3> >;

  template class T_DifferentialOperator<DiffOpCurlEdge<2> >;
  template class T_DifferentialOperator<DiffOpCurlEdge<3> >;
  */


}
