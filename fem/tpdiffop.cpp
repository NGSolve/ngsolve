#include <fem.hpp>
#include "tpintrule.hpp"
#include "tpdiffop.hpp"

    
namespace ngfem {
    void ProlongateCoefficientFunction :: Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const
    {
        const TPMappedIntegrationRule & tpmir = dynamic_cast<const TPMappedIntegrationRule &>(ir);
        auto & irs = tpmir.GetIRs();
        coef->Evaluate(*irs[1-prolongateto],values);
        if(prolongateto == 1)
            for(int i=tpmir.GetIRs()[0]->Size()-1;i>=0;i--)
                values.Rows(i*tpmir.GetIRs()[1]->Size(),(i+1)*tpmir.GetIRs()[1]->Size()) = values.Row(i)(0);
        if(prolongateto == 0)
            for(int i=1;i<tpmir.GetIRs()[0]->Size();i++)
                values.Rows(i*tpmir.GetIRs()[1]->Size(),(i+1)*tpmir.GetIRs()[1]->Size()) = values.Rows(0,tpmir.GetIRs()[1]->Size());
    }


    void TPDifferentialOperator :: 
    Apply (const FiniteElement & fel,
        const BaseMappedIntegrationRule & mir,
        FlatVector<double> x, 
        FlatMatrix<double> flux,
        LocalHeap & lh) const
    {
      const TPHighOrderFE & tpfel = dynamic_cast<const TPHighOrderFE &>(fel);
      const TPMappedIntegrationRule & tpmir = dynamic_cast<const TPMappedIntegrationRule &>(mir);
      auto & elements = tpfel.elements;
      auto & irs = tpmir.GetIRs();
      int ndof0 = elements[0]->GetNDof();
      int ndof1 = elements[1]->GetNDof();
      int dim0 = evaluators[0]->Dim();
      int dim1 = evaluators[1]->Dim();
      int nip0 = irs[0]->IR().Size();
      int nip1 = irs[1]->IR().Size();
      FlatMatrix<double, ColMajor> shape0( nip0*dim0, ndof0,lh ); 
      FlatMatrix<double, ColMajor> shape1( nip1*dim1, ndof1,lh ); 
      FlatMatrix<double> fcoefs( ndof0, ndof1, &x(0) ); //TODO slicematrix
      
      evaluators[0]->CalcMatrix( *elements[0], *irs[0], shape0, lh ); 
      evaluators[1]->CalcMatrix( *elements[1], *irs[1], shape1, lh ); 
      if(dim0 == 1)
      {
        FlatMatrix<double,ColMajor> fvals( nip1*dim1, nip0*dim0, &flux(0) );
        FlatMatrix<double,ColMajor> helper(ndof1,nip0*dim0,lh);
        helper = Trans(fcoefs) * Trans(shape0); // (a*b'*c')' = c*(a*b')' = c*b*a'
        fvals = shape1 * helper;
      }
      if(dim0 > 1)
      {
        FlatMatrix<double,RowMajor> fvals( nip1*dim1, nip0*dim0, lh );
        FlatMatrix<double,RowMajor> helper(ndof1,nip0*dim0,lh);
        helper = Trans(fcoefs) * Trans(shape0);
        fvals = shape1 * helper;
        for(int i=0;i<nip0;i++)
          flux.Rows(i*nip1,(i+1)*nip1) = fvals.Cols(dim0*i,dim0*(i+1));         
      }
    }
    

    void TPDifferentialOperator :: 
    ApplyTrans (const FiniteElement & fel,
        const BaseMappedIntegrationRule & mir,
        FlatMatrix<double> flux,
        FlatVector<double> x, 
        LocalHeap & lh) const
    {
      const TPHighOrderFE & tpfel = dynamic_cast<const TPHighOrderFE &>(fel);
      const TPMappedIntegrationRule & tpmir = dynamic_cast<const TPMappedIntegrationRule &>(mir);

      auto & elements = tpfel.elements;
      //auto & irs = tpmir.irs;
      auto & irs = tpmir.GetIRs();

      int ndof0 = elements[0]->GetNDof();
      int ndof1 = elements[1]->GetNDof();
      int dim0 = evaluators[0]->Dim();
      int dim1 = evaluators[1]->Dim();
      int nip0 = irs[0]->IR().Size();
      int nip1 = irs[1]->IR().Size();
      FlatMatrix<double, ColMajor> shape0( nip0*dim0, ndof0,lh ); 
      FlatMatrix<double, ColMajor> shape1( nip1*dim1, ndof1,lh ); 
      evaluators[0]->CalcMatrix( *elements[0], *irs[0], shape0, lh ); 
      evaluators[1]->CalcMatrix( *elements[1], *irs[1], shape1, lh ); 

      if(dim0 == 1)
      {
        FlatMatrix<double> fvals( nip0*dim0, nip1*dim1, &flux(0,0) );
        FlatMatrix<double> fcoefs( ndof0, ndof1, &x(0) ); //TODO slicematrix
        FlatMatrix<double> helper(nip0*dim0,ndof1,lh);
        helper = fvals*shape1;
        fcoefs = Trans(shape0) * helper;
      }
      else      
      {
        FlatMatrix<double> fvals( nip0, nip1*dim0, &flux(0,0) );
        FlatMatrix<double> fcoefs( ndof0, ndof1, &x(0) ); //TODO slicematrix
        FlatMatrix<double> helper(nip0*dim0,ndof1,lh);
        
        FlatMatrix<double> fvals1( nip0*dim0, nip1*dim1, lh );
        for(int i=0;i<nip1;i++)
          for(int j=0;j<nip0;j++)
            fvals1.Rows(dim0*j,dim0*(j+1)).Col(i) = (fvals.Cols(dim0*i,dim0*(i+1)).Row(j));
        
        helper = fvals1 * shape1;
        fcoefs = Trans(shape0) * helper;
      }
    }

    void TPDifferentialOperator ::
    ApplyX(const FiniteElement &fel,
        const BaseMappedIntegrationRule & mirx,
        FlatMatrix<double> flux,
        SliceMatrix<double> x,
        LocalHeap & lh) const
        {
          int dimx = evaluators[0]->Dim();
          int dimy = evaluators[1]->Dim();
          int nipx = mirx.IR().Size();
          FlatMatrix<double, ColMajor> bmatx( nipx*dimx, fel.GetNDof(),lh );
          evaluators[0]->CalcMatrix(fel,mirx,bmatx,lh);
          if(dimx == 1)
          {
            FlatMatrix<> resultmat(nipx*dimx,x.Width(), &flux(0,0));
            resultmat = bmatx*x;
          }
          else
          {
            FlatMatrix<> resultmat(nipx*dimx,x.Width(), lh);
            resultmat = bmatx*x;
            for(int k=0;k<x.Height();k+=2)
              flux.Rows(k*nipx,(k+1)*nipx) = Trans(resultmat.Rows(k,k+1));
          }
        }

    void TPDifferentialOperator ::
    ApplyXTrans(const FiniteElement &fel,
        const BaseMappedIntegrationRule & mirx,
        FlatMatrix<double> flux,
        SliceMatrix<double> x,
        LocalHeap & lh) const
        {
          int dimx = evaluators[0]->Dim();
          int dimy = evaluators[1]->Dim();
          int nipx = mirx.IR().Size();
          int nipy = flux.Height()/(dimy*dimx*nipx);
          FlatMatrix<double, ColMajor> bmatx( nipx*dimx, fel.GetNDof(),lh );
          evaluators[0]->CalcMatrix(fel,mirx,bmatx,lh);
          if(dimx == 1)
          {
            FlatMatrix<> proxyvaluesasmat(nipx, nipy*dimy, &flux(0,0));
            x = Trans(bmatx)*proxyvaluesasmat;
          }
          else
          {
            FlatMatrix<double> proxyvaluesasmat( nipx, nipy*dimx, &flux(0,0) );
            FlatMatrix<double> proxyvaluesasmat1( nipx*dimx, nipy*dimy, lh );
            for(int i=0;i<nipy;i++)
              for(int j=0;j<nipx;j++)
                proxyvaluesasmat1.Rows(dimx*j,dimx*(j+1)).Col(i) = (proxyvaluesasmat.Cols(dimx*i,dimx*(i+1)).Row(j));
            x = proxyvaluesasmat1*bmatx;
          }
        }

    void TPDifferentialOperator ::
    ApplyY(const FiniteElement &fel,
        const BaseMappedIntegrationRule & miry,
        FlatMatrix<double> flux,
        SliceMatrix<double> x,
        LocalHeap & lh) const
        {
          int dimx = evaluators[0]->Dim();
          int dimy = evaluators[1]->Dim();
          int nipy = miry.IR().Size();
          FlatMatrix<double, ColMajor> bmaty( nipy*dimy, fel.GetNDof(),lh );
          evaluators[1]->CalcMatrix(fel,miry,bmaty,lh);
          if(dimx == 1)
          {
            FlatMatrix<> resultmat(x.Height(),nipy*dimy, &flux(0,0));
            resultmat = x*Trans(bmaty);
          }
          else
          {
            FlatMatrix<> resultmat(x.Height(),nipy*dimy, lh);
            resultmat = x*Trans(bmaty);
            for(int k=0;k<x.Height();k+=2)
              flux.Rows(k*nipy,(k+1)*nipy) = Trans(resultmat.Rows(k,k+1));
          }
        }

    void TPDifferentialOperator :: 
    ApplyYTrans(const FiniteElement &fel,
        const BaseMappedIntegrationRule & miry,
        FlatMatrix<double> flux,
        SliceMatrix<double> x,
        LocalHeap & lh) const
        {
          int dimx = evaluators[0]->Dim();
          int dimy = evaluators[1]->Dim();
          int nipy = miry.IR().Size();
          int nipx = flux.Height()/(dimx*dimy*nipy);
          FlatMatrix<double, ColMajor> bmaty( nipy*dimy, fel.GetNDof(),lh );
          evaluators[1]->CalcMatrix(fel,miry,bmaty,lh);
          if(dimx == 1)
          {
            FlatMatrix<> proxyvaluesasmat(nipx, nipy*dimy, &flux(0,0));
            x = proxyvaluesasmat*bmaty;
          }
          else
          {
            FlatMatrix<double> proxyvaluesasmat( nipx, nipy*dimx, &flux(0,0) );
            FlatMatrix<double> proxyvaluesasmat1( nipx*dimx, nipy*dimy, lh );
            for(int i=0;i<nipy;i++)
              for(int j=0;j<nipx;j++)
                proxyvaluesasmat1.Rows(dimx*j,dimx*(j+1)).Col(i) = (proxyvaluesasmat.Cols(dimx*i,dimx*(i+1)).Row(j));
            x = proxyvaluesasmat1*bmaty;
          }
        }
}
