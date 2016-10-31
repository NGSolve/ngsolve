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
              SliceMatrix<double,ColMajor> mat, 
              LocalHeap & lh) const 
  {
    throw Exception (string("Error: DifferentialOperator::CalcMatrix called for base class, type = ")
                     + typeid(*this).name());
  }

  void DifferentialOperator ::
  CalcMatrix (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              SliceMatrix<Complex,ColMajor> mat, 
              LocalHeap & lh) const 
  {
    throw Exception (string("Error: DifferentialOperator::CalcMatrix<Complex> called for base class, type = ")
                     + typeid(*this).name());
  }

  void DifferentialOperator ::
  CalcMatrix (const FiniteElement & fel,
              const BaseMappedIntegrationRule & mir,
              SliceMatrix<double,ColMajor> mat,   
              LocalHeap & lh) const
  {
    int dim = Dim(); 
    for (int i = 0; i < mir.Size(); i++)
      CalcMatrix (fel, mir[i], mat.Rows(i*dim, (i+1)*dim), lh);
  }

  void DifferentialOperator ::
  CalcMatrix (const FiniteElement & fel,
              const BaseMappedIntegrationRule & mir,
              SliceMatrix<Complex,ColMajor> mat,   
              LocalHeap & lh) const
  {
    int dim = Dim(); 
    for (int i = 0; i < mir.Size(); i++)
      CalcMatrix (fel, mir[i], mat.Rows(i*dim, (i+1)*dim), lh);
  }
  
  void DifferentialOperator ::
  CalcMatrix (const FiniteElement & fel,
              const SIMD_BaseMappedIntegrationRule & mir,
              ABareMatrix<double> mat) const
  {
    throw ExceptionNOSIMD(string("Error: DifferentialOperator::CalcMatrix does not support SIMD, type = ")
                          + typeid(*this).name());
    
  }
  
  void DifferentialOperator ::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationPoint & mip,
         FlatVector<double> x, 
         FlatVector<double> flux,
         LocalHeap & lh) const
  {
#ifndef FASTCOMPILE
    cout << "called base class apply, type = " << typeid(*this).name() << endl;
#endif
    HeapReset hr(lh);
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
#ifndef FASTCOMPILE
    cout << "called base class apply, complex" << endl;
#endif
    HeapReset hr(lh);
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
  Apply (const FiniteElement & bfel,
         const SIMD_BaseMappedIntegrationRule & bmir,
         BareSliceVector<double> x, 
         ABareSliceMatrix<double> flux) const
  // LocalHeap & lh) const
  {
    throw Exception (string("DifferentialOperator :: Apply ( ... SIMD ... ) not overloaded for class ")
                     + typeid(*this).name());
  }
  
  
  void DifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<double> flux,
              FlatVector<double> x, 
              LocalHeap & lh) const 
  {
#ifndef FASTCOMPILE
    cout << "called base class apply trans, type = " << typeid(*this).name() << endl;
#endif
    HeapReset hr(lh);
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
#ifndef FASTCOMPILE
    cout << "called base class apply trans complex, type = " << typeid(*this).name() << endl;
#endif
    HeapReset hr(lh);
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
    HeapReset hr(lh);
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
    HeapReset hr(lh);
    FlatVector<Complex> hx(x.Size(), lh);
    x = 0.0;
    for (int i = 0; i < mir.Size(); i++)
      {
        ApplyTrans (fel, mir[i], flux.Row(i), hx, lh);
        x += hx;
      }
  }

  void DifferentialOperator ::
  AddTrans (const FiniteElement & bfel,
            const SIMD_BaseMappedIntegrationRule & bmir,
            ABareMatrix<double> flux,
            BareSliceVector<double> x) const
  // LocalHeap & lh) const
  {
    throw Exception (string("DifferentialOperator :: AddTrans ( ... SIMD ... ) not overloaded") +
                     + typeid(*this).name());
  }

  

  BlockDifferentialOperator :: ~BlockDifferentialOperator ()  { ; }


  void BlockDifferentialOperator ::
  CalcMatrix (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              SliceMatrix<double,ColMajor> mat, 
              LocalHeap & lh) const 
  {
    HeapReset hr(lh);
    FlatMatrix<double,ColMajor> mat1(diffop->Dim(), fel.GetNDof(), lh);
    diffop->CalcMatrix (fel, mip, mat1, lh);
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
  CalcMatrix (const FiniteElement & fel,
              const SIMD_BaseMappedIntegrationRule & mir,
              ABareMatrix<double> mat) const
  {
    throw ExceptionNOSIMD("BlockDifferentialOperator::CalcMatrix does not support SIMD");
  }
  

  void BlockDifferentialOperator ::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationPoint & mip,
         FlatVector<double> x, 
         FlatVector<double> flux,
         LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatVector<> hx(fel.GetNDof(), lh);
    FlatVector<> hflux(diffop->Dim(), lh);
    
    if (comp == -1)
      {
        for (int k = 0; k < dim; k++)
          {
            hx = x.Slice(k, dim);
            diffop->Apply(fel, mip, hx, hflux, lh);
            flux.Slice(k,dim) = hflux;
          }
      }
    else
      {
        hx = x.Slice(comp, dim);
        diffop->Apply(fel, mip, hx, hflux, lh);
        flux.Slice(comp,dim) = hflux;
      }
  }


  void BlockDifferentialOperator ::
  Apply (const FiniteElement & fel,
         const SIMD_BaseMappedIntegrationRule & mir,
         BareSliceVector<double> x, 
         ABareSliceMatrix<double> flux) const
  // LocalHeap & lh) const
  {
    if (comp == -1)
      for (int k = 0; k < dim; k++)
        diffop->Apply(fel, mir, x.Slice(k, dim), flux.RowSlice(k,dim));
    else
      diffop->Apply(fel, mir, x.Slice(comp, dim), flux.RowSlice(comp,dim));
  }
  
  
  void BlockDifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<double> flux,
              FlatVector<double> x, 
              LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatVector<> hx(fel.GetNDof(), lh);
    FlatVector<> hflux(diffop->Dim(), lh);
    
    if (comp == -1)
      {
        for (int k = 0; k < dim; k++)
          {
            hflux = flux.Slice(k, dim);
            diffop->ApplyTrans(fel, mip, hflux, hx, lh);
            x.Slice(k,dim) = hx;
          }
      }
    else
      {
        hflux = flux.Slice(comp, dim);
        diffop->ApplyTrans(fel, mip, hflux, hx,lh);
        x = 0.0;
        x.Slice(comp,dim) = hx;
      }
  }

  void BlockDifferentialOperator ::
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationPoint & mip,
              FlatVector<Complex> flux,
              FlatVector<Complex> x, 
              LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatVector<Complex> hx(fel.GetNDof(), lh);
    FlatVector<Complex> hflux(diffop->Dim(), lh);
    
    if (comp == -1)
      {
        for (int k = 0; k < dim; k++)
          {
            hflux = flux.Slice(k, dim);
            diffop->ApplyTrans(fel, mip, hflux, hx, lh);
            x.Slice(k,dim) = hx;
          }
      }
    else
      {
        hflux = flux.Slice(comp, dim);
        diffop->ApplyTrans(fel, mip, hflux, hx,lh);
        x = 0.0;
        x.Slice(comp,dim) = hx;
      }
  }

    
  void BlockDifferentialOperator ::
  AddTrans (const FiniteElement & fel,
            const SIMD_BaseMappedIntegrationRule & mir,
            ABareMatrix<double> flux,
            BareSliceVector<double> x) const
  // LocalHeap & lh) const
  {
    if (comp == -1)
      for (int k = 0; k < dim; k++)
        diffop->AddTrans(fel, mir, flux.RowSlice(k,dim), x.Slice(k,dim));
    else
      diffop->AddTrans(fel, mir, flux.RowSlice(comp,dim), x.Slice(comp,dim));
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
