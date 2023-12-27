#include "tpintrule.hpp"
#include "tpdiffop.hpp"

    
namespace ngfem 
{

  double ProlongateCoefficientFunction :: Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    IntegrationPoint iphelp = ip.IP();
    // cout << iphelp << endl;
    const ElementTransformation * trafo = &(ip.GetTransformation());
    DimMappedIntegrationPoint<1> miphelp(iphelp(0), *trafo);
    if(prolongateto == 0)
      {
        if(dimx == 1)
          {
            ip.GetPoint()[0] = ip.GetPoint()[1];
            if(dimy == 2)
              {
                ip.GetPoint()[1] = ip.GetPoint()[2];
                ip.GetPoint()[2] = 0.0;
              }
          }
        else if(dimx == 2)
          ip.GetPoint()[0] = ip.GetPoint()[2];
      }
    return coef->Evaluate(ip);
  }

  void ProlongateCoefficientFunction :: Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<double> values) const
  {
    const TPMappedIntegrationRule * tpmir = static_cast<const TPMappedIntegrationRule *>(&ir);
    auto & irs = tpmir->GetIRs();
    auto tempvals = values.AddSize(irs[0]->Size()*irs[1]->Size(),Dimension());
    coef->Evaluate(*irs[1-prolongateto],tempvals.Rows(0,irs[1-prolongateto]->Size()) );        
    if(prolongateto == 1)
      for(int i=irs[0]->Size()-1;i>=0;i--)
        tempvals.Rows(i*irs[1]->Size(),(i+1)*irs[1]->Size()) = tempvals.Row(i)(0);
    if(prolongateto == 0)
      for(int i=1;i<irs[0]->Size();i++)
        tempvals.Rows(i*irs[1]->Size(),(i+1)*irs[1]->Size()) = tempvals.Rows(0,irs[1]->Size());
  }

  void ProlongateCoefficientFunction :: EvaluateStdRule (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const
  {
    for(int i=0;i< ir.Size();i++)
      values(i,0) = Evaluate(ir[i]);
    // return;
    // if(prolongateto == 0)
      // for(int i : Range(ir.Size()) )
      // {
        // cout << ir[i].GetPoint() << endl;
        // if(dimx == 1)
        // {
          // ir[i].GetPoint()[0] = ir[i].GetPoint()[1];
          // if(dimy == 2)
            // ir[i].GetPoint()[1] = ir[i].GetPoint()[2];
        // }
        // else if(dimx == 2)
          // ir[i].GetPoint()[0] = ir[i].GetPoint()[2];
        // cout << ir[i].GetPoint() << endl;
      // }
    // cout << "Evaluating"<<endl;
    // // coef->Evaluate(ir,values);
    // cout << "Evaluating"<<endl;
  }
  
  void TPDifferentialOperator :: Apply (
            const FiniteElement & fel,
            const BaseMappedIntegrationRule & mir,
            BareSliceVector<double> x, 
            BareSliceMatrix<double> flux,
            LocalHeap & lh) const
  {
    const TPHighOrderFE & tpfel = static_cast<const TPHighOrderFE &>(fel);
    const TPMappedIntegrationRule & tpmir = static_cast<const TPMappedIntegrationRule &>(mir);
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
    FlatMatrix<double> fcoefs( ndof0, ndof1, &x(0) );
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
        flux.Rows(i*nip1,(i+1)*nip1).AddSize(nip1,dim0) = fvals.Cols(dim0*i,dim0*(i+1));
    }
  }
    

  void TPDifferentialOperator :: ApplyTrans (
            const FiniteElement & fel,
            const BaseMappedIntegrationRule & mir,
            FlatMatrix<double> flux,
            BareSliceVector<double> x, 
            LocalHeap & lh) const
  {
    const TPHighOrderFE & tpfel = static_cast<const TPHighOrderFE &>(fel);
    const TPMappedIntegrationRule & tpmir = static_cast<const TPMappedIntegrationRule &>(mir);
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
    evaluators[0]->CalcMatrix( *elements[0], *irs[0], shape0, lh ); 
    evaluators[1]->CalcMatrix( *elements[1], *irs[1], shape1, lh );
    if(dim0 == 1)
    {
      FlatMatrix<double> fvals( nip0*dim0, nip1*dim1, &flux(0,0) );
      FlatMatrix<double> fcoefs( ndof0, ndof1, &x(0) );
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

  void TPDifferentialOperator :: ApplyX(
            const FiniteElement &fel,
            const BaseMappedIntegrationRule & mirx,
            FlatMatrix<double> flux,
            SliceMatrix<double> x,
            LocalHeap & lh) const
  {
    int dimx = evaluators[0]->Dim();
    int dimy = evaluators[1]->Dim();
    int nipx = mirx.IR().Size();
    int nipy = x.Width()/dimy;
    // for(int i:Range(mirx.IR().Size()) )
      // cout << mirx[i].GetPoint()<<endl;
    FlatMatrix<double, ColMajor> bmatx( nipx*dimx, fel.GetNDof(),lh );
    evaluators[0]->CalcMatrix(fel,mirx,bmatx,lh);
    if(dimx == 1)
    {
      FlatMatrix<> resultmat(nipx*dimx,x.Width(), &flux(0,0));
      resultmat = bmatx*x | Lapack;
    }
    else
    {
      FlatMatrix<> resultmat(nipx*dimx,x.Width(), lh);
      resultmat = bmatx*x | Lapack;
      // cout << "resultmat = "<<endl<<resultmat<<endl;
      for(int k=0;k<flux.Height()/nipy;k++)
        flux.Rows(k*nipy,(k+1)*nipy) = Trans(resultmat.Rows(dimx*k,dimx*(k+1)));
      // cout << "flux = "<<endl<<flux<<endl;
    }
  }

  void TPDifferentialOperator :: ApplyXTrans(
            const FiniteElement &fel,
            const BaseMappedIntegrationRule & mirx,
            FlatMatrix<double> flux,
            SliceMatrix<double> x,
            LocalHeap & lh) const
  {
    int dimx = evaluators[0]->Dim();
    int dimy = evaluators[1]->Dim();
    int nipx = mirx.IR().Size();
    int nipy = flux.Height()/nipx;
    FlatMatrix<double, ColMajor> bmatx( nipx*dimx, fel.GetNDof(),lh );
    evaluators[0]->CalcMatrix(fel,mirx,bmatx,lh);
    if(dimx == 1)
    {
      FlatMatrix<> proxyvaluesasmat(nipx, nipy*dimy, &flux(0,0));
      // cout << "Input as mat = "<<endl<<proxyvaluesasmat<<endl;
      x = Trans(bmatx)*proxyvaluesasmat | Lapack;
    }
    else
    {
      FlatMatrix<double> proxyvaluesasmat( nipx, nipy*dimx, &flux(0,0) );
      FlatMatrix<double> proxyvaluesasmat1( nipx*dimx, nipy*dimy, lh );
      for(int i=0;i<nipy;i++)
        for(int j=0;j<nipx;j++)
          proxyvaluesasmat1.Rows(dimx*j,dimx*(j+1)).Col(i) = (proxyvaluesasmat.Cols(dimx*i,dimx*(i+1)).Row(j));
      x = proxyvaluesasmat1*bmatx | Lapack;
    }
  }

  void TPDifferentialOperator :: ApplyY(
            const FiniteElement &fel,
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
      FlatMatrix<double, ColMajor> resultmat(x.Height(),nipy*dimy, lh);
      resultmat = x*Trans(bmaty);
      for(int k=0;k<x.Height()/dimx;k++)
        flux.Rows(k*nipy,(k+1)*nipy) = Trans(resultmat.Rows(dimx*k,dimx*(k+1)));
    }
  }

  void TPDifferentialOperator :: ApplyYTrans(
            const FiniteElement &fel,
            const BaseMappedIntegrationRule & miry,
            FlatMatrix<double> flux,
            SliceMatrix<double> x,
            LocalHeap & lh) const
  {
    int dimx = evaluators[0]->Dim();
    int dimy = evaluators[1]->Dim();
    int nipy = miry.IR().Size();
    int nipx = x.Height()/dimx;
    FlatMatrix<double, ColMajor> bmaty( nipy*dimy, fel.GetNDof(),lh );
    evaluators[1]->CalcMatrix(fel,miry,bmaty,lh);
    if(dimx == 1)
    {
      FlatMatrix<> proxyvaluesasmat(nipx, nipy*dimy, &flux(0,0));
      x = proxyvaluesasmat*bmaty | Lapack;
    }
    else
    {
      FlatMatrix<double> proxyvaluesasmat( nipx, nipy*dimx, &flux(0,0) );
      FlatMatrix<double> proxyvaluesasmat1( nipx*dimx, nipy*dimy, lh );
      for(int i=0;i<nipy;i++)
        for(int j=0;j<nipx;j++)
          proxyvaluesasmat1.Rows(dimx*j,dimx*(j+1)).Col(i) = (proxyvaluesasmat.Cols(dimx*i,dimx*(i+1)).Row(j));
      x = proxyvaluesasmat1*bmaty | Lapack;
    }
  }


  void TPBlockDifferentialOperator :: Apply (
            const FiniteElement & fel,
            const BaseMappedIntegrationRule & mir,
            BareSliceVector<double> x, 
            BareSliceMatrix<double> flux,
            LocalHeap & lh) const
  {
    const TPHighOrderFE & tpfel = static_cast<const TPHighOrderFE &>(fel);
    const TPMappedIntegrationRule & tpmir = static_cast<const TPMappedIntegrationRule &>(mir);
    auto & elements = tpfel.elements;
    auto & evalx = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(0);
    auto & evaly = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(1);
    auto & irs = tpmir.GetIRs();
    int ndof0 = elements[0]->GetNDof();
    int ndof1 = elements[1]->GetNDof();
    int dim0 = evalx->Dim();
    int dim1 = evaly->Dim();
    int nip0 = irs[0]->IR().Size();
    int nip1 = irs[1]->IR().Size();
    FlatMatrix<double, ColMajor> shape0( nip0*dim0, ndof0,lh );
    FlatMatrix<double, ColMajor> shape1( nip1*dim1, ndof1,lh );
    FlatMatrix<double> fcoefs( ndof0, ndof1*BlockDim(), &x(0) );
    evalx->CalcMatrix( *elements[0], *irs[0], shape0, lh );
    evaly->CalcMatrix( *elements[1], *irs[1], shape1, lh );
    if(dim0 == 1)
    {
      FlatMatrix<double,ColMajor> helper1(nip0*dim0,ndof1*BlockDim(),lh);
      FlatMatrix<> helper2(ndof1,nip0*dim0*BlockDim(),&helper1(0,0));
      FlatMatrix<double,ColMajor> result(nip1*dim1,nip0*dim0*BlockDim(),lh);
      // FlatMatrix<double,ColMajor> fluxCM(flux.Height(),flux.Width(),&result(0,0));
      FlatMatrix<double,ColMajor> fluxCM(nip0*nip1, dim1,&result(0,0));
      helper1 = shape0 * fcoefs;
      result =  shape1*helper2;
      flux.AddSize(nip0*nip1, dim1) = fluxCM;
    }
    if(dim0 > 1)
    {
      // FlatMatrix<double,RowMajor> fvals( nip1*dim1, nip0*dim0*BlockDim(), lh );
      // FlatMatrix<double,RowMajor> helper(ndof1,nip0*dim0,lh);
      // helper = Trans(fcoefs) * Trans(shape0);
      // fvals = shape1 * helper;
      // for(int i=0;i<nip0;i++)
        // flux.Rows(i*nip1,(i+1)*nip1) = fvals.Cols(dim0*i,dim0*(i+1));

      // FlatMatrix<Vec<3>,RowMajor> fvals( nip1*dim1, nip0*dim0, lh );
      // FlatMatrix<double,RowMajor> fvalsdb( nip1*dim1, nip0*dim0*BlockDim(), reinterpret_cast<double * >(&fvals(0,0)) );
      // FlatMatrix<Vec<3>,RowMajor> helper(ndof1,nip0*dim0,lh);
      // FlatMatrix<Vec<3>> fcoefs1(ndof0,ndof1,reinterpret_cast<Vec<3> *>(&fcoefs(0,0)));
      // helper = fcoefs1 * Trans(shape0);
      // fvals = shape1 * helper;
      // for(int i=0;i<nip0;i++)
        // flux.Rows(i*nip1,(i+1)*nip1) = fvalsdb.Cols(dim0*i*BlockDim(),dim0*(i+1)*BlockDim());
    }      
  }
    

  void TPBlockDifferentialOperator :: ApplyTrans (
            const FiniteElement & fel,
            const BaseMappedIntegrationRule & mir,
            FlatMatrix<double> flux,
            BareSliceVector<double> x, 
            LocalHeap & lh) const
  {
    const TPHighOrderFE & tpfel = static_cast<const TPHighOrderFE &>(fel);
    const TPMappedIntegrationRule & tpmir = static_cast<const TPMappedIntegrationRule &>(mir);
    auto & elements = tpfel.elements;
    auto & evalx = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(0);
    auto & evaly = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(1);
    auto & irs = tpmir.GetIRs();
    int ndof0 = elements[0]->GetNDof();
    int ndof1 = elements[1]->GetNDof();
    int dim0 = evalx->Dim();
    int dim1 = evaly->Dim();
    int nip0 = irs[0]->IR().Size();
    int nip1 = irs[1]->IR().Size();
    FlatMatrix<double, ColMajor> shape0( nip0*dim0, ndof0,lh ); 
    FlatMatrix<double, ColMajor> shape1( nip1*dim1, ndof1,lh ); 
    evalx->CalcMatrix( *elements[0], *irs[0], shape0, lh ); 
    evaly->CalcMatrix( *elements[1], *irs[1], shape1, lh ); 

    if(dim0 == 1)
    {
      // FlatMatrix<Vec<3> > fvalsvec( nip0*dim0, nip1*dim1, reinterpret_cast<Vec<3> * >(&flux(0,0)) );
      // FlatMatrix<Vec<3> > fcoefsvec( ndof0, ndof1, reinterpret_cast<Vec<3> * >(&x(0,0)) );
      // FlatMatrix<Vec<3> > helpervec(nip0*dim0,ndof1*BlockDim(),lh);
      // helpervec = fvalsvec*shape1;
      // fcoefsvec = Trans(shape0)*helpervec;
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

  void TPBlockDifferentialOperator :: ApplyX(
            const FiniteElement &fel,
            const BaseMappedIntegrationRule & mirx,
            FlatMatrix<double> flux,
            SliceMatrix<double> x,
            LocalHeap & lh) const
  {
    auto & evalx = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(0);
    int dimx = evalx->Dim();
    int nipx = mirx.IR().Size();
    FlatMatrix<double, ColMajor> bmatx( nipx*dimx, fel.GetNDof(),lh );
    evalx->CalcMatrix(fel,mirx,bmatx,lh);
    if(dimx == 1)
    {
      FlatMatrix<> x_reshape(1.0/BlockDim()*x.Height(),BlockDim()*x.Width(),&x(0,0));
      FlatMatrix<> resultmat(bmatx.Height(),x_reshape.Width(),&flux(0,0));
      resultmat = bmatx*x_reshape;
    }
    else
    {
      // FlatMatrix<> resultmat(nipx*dimx,x.Width(), lh);
      // resultmat = bmatx*x;
      FlatMatrix<> x_reshape(1.0/BlockDim()*x.Height(),BlockDim()*x.Width(),&x(0,0));
      FlatMatrix<> resultmat(bmatx.Height(),x_reshape.Width(),lh);
      resultmat = bmatx*x_reshape;
      // cout << "TPBlockDifferentialOperator::ApplyX"<<endl;
      // cout << "x = "<<endl<<x << endl;
      // cout << "x_reshpae = "<<endl<<x_reshape << endl;
      // cout << "Result mat = "<<endl<<resultmat1 << endl;
      // cout << "flux = "<<endl<<flux << endl;
      
      for(int k=0;k<x.Height();k+=2)
        flux.Rows(k*nipx,(k+1)*nipx) = Trans(resultmat.Rows(dimx*k,dimx*(k+1)));

    }
  }

  void TPBlockDifferentialOperator :: ApplyXTrans(
            const FiniteElement &fel,
            const BaseMappedIntegrationRule & mirx,
            FlatMatrix<double> flux,
            SliceMatrix<double> x,
            LocalHeap & lh) const
  {
    auto & evalx = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(0);
    auto & evaly = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(1);
    int dimx = evalx->Dim();
    int dimy = evaly->Dim();
    int nipx = mirx.IR().Size();
    int nipy = x.Width();
    FlatMatrix<double, ColMajor> bmatx( nipx*dimx, fel.GetNDof(),lh );
    evalx->CalcMatrix(fel,mirx,bmatx,lh);
    if(dimx == 1)
    {
      FlatMatrix<> resultmat(fel.GetNDof(),BlockDim()*nipy*dimy,&x(0,0));
      resultmat = Trans(bmatx)*flux;
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

  void TPBlockDifferentialOperator :: ApplyY(
            const FiniteElement &fel,
            const BaseMappedIntegrationRule & miry,
            FlatMatrix<double> flux,
            SliceMatrix<double> x,
            LocalHeap & lh) const
  {
    auto & evalx = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(0);
    auto & evaly = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(1);
    int dimx = evalx->Dim();
    int dimy = evaly->Dim();
    int nipy = miry.IR().Size();
    FlatMatrix<double, ColMajor> bmaty( nipy*dimy, fel.GetNDof(),lh );
    evaly->CalcMatrix(fel,miry,bmaty,lh);
    if(dimx == 1)
    {
      FlatMatrix<double> x_reshape(1.0/BlockDim()*x.Height(),BlockDim()*x.Width(),&x(0,0));
      FlatMatrix<double,ColMajor> fluxCM(flux.Height(),flux.Width(),lh);
      FlatMatrix<double,ColMajor> resultmat(bmaty.Height(),x_reshape.Width(),&fluxCM(0,0));
      resultmat = bmaty * x_reshape;
      flux = fluxCM;
    }
    else
    {
      // FlatMatrix<Vec<3>, ColMajor> resultmat(x.Height(),nipy*dimy, lh);
      // FlatMatrix<double, ColMajor> resultmat_dbl(x.Height(),BlockDim()*nipy*dimy, reinterpret_cast<double *>(&resultmat(0,0)));
      // SliceMatrix< Vec<3> > x_vec(x.Height(),0.5*x.Width(),0.5*x.Dist(),reinterpret_cast<Vec<3> * >(&x(0,0)));
      // resultmat = x_vec*Trans(bmaty);
      // for(int k=0;k<x.Height()/dimx;k++)
        // flux.Rows(k*nipy,(k+1)*nipy) = Trans(resultmat_dbl.Rows(dimx*k*BlockDim(),dimx*(k+1)*BlockDim()));
      // cout << "TPBlockDifferentialOperator::ApplyY" << endl;
      // cout << "x = "<< endl<<x << endl;
      // cout << "x_vec = "<<endl<<x_vec<<endl;
      // cout << "resultmat = "<<endl<<resultmat<<endl;
      // cout << "flux = "<<endl<<flux<<endl;
    }
  }

  void TPBlockDifferentialOperator :: ApplyYTrans(
            const FiniteElement &fel,
            const BaseMappedIntegrationRule & miry,
            FlatMatrix<double> flux,
            SliceMatrix<double> x,
            LocalHeap & lh) const
  {
    auto & evalx = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(0);
    auto & evaly = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(1);
    int dimx = evalx->Dim();
    int dimy = evaly->Dim();
    int nipy = miry.IR().Size();
    int nipx = x.Height()/dimx;
    FlatMatrix<double, ColMajor> bmaty( nipy*dimy, fel.GetNDof(),lh );
    evaly->CalcMatrix(fel,miry,bmaty,lh);
    if(dimx == 1)
    {
      FlatMatrix<double> helper(flux.Width(),flux.Height(),lh);
      helper = Trans(flux);
      FlatMatrix<double,RowMajor> helper2(BlockDim()*nipx,nipy*dimy,&helper(0,0));
      FlatMatrix<double,ColMajor> temp(BlockDim()*nipx,bmaty.Width(),lh);
      temp = helper2*bmaty;
      FlatMatrix<double,ColMajor> xt(x.Height(),x.Width(),&temp(0,0));
      x = xt;
    }
    else
    {
      // FlatMatrix<Vec<3> > proxyvaluesasmat( nipx, nipy*dimx, reinterpret_cast<Vec<3> *>(&flux(0,0)));
      // FlatMatrix<Vec<3> > proxyvaluesasmat1( nipx*dimx, nipy*dimy, lh );
      // for(int i=0;i<nipy;i++)
        // for(int j=0;j<nipx;j++)
      // proxyvaluesasmat1.Rows(dimx*j,dimx*(j+1)).Col(i) = (proxyvaluesasmat.Cols(dimx*i,dimx*(i+1)).Row(j));
      // SliceMatrix<Vec<3> > x_vec(x.Height(), 1.0/3.0*x.Width(),1.0/3.0*x.Dist(), reinterpret_cast<Vec<3> *>(&x(0,0)));
      // x_vec = proxyvaluesasmat1*bmaty;
      // cout << "TPBlockDifferentialOperator::ApplyYTrans" << endl;
      // cout << "flux (in) = "<<endl<<flux<<endl;
      // cout << "fluxasmat = "<<endl<<proxyvaluesasmat<<endl;
      // cout << "fluxasmat = "<<endl<<proxyvaluesasmat1<<endl;
      // cout << "x_vec = "<<endl<<x_vec<<endl;
      // cout << "x = "<< endl<<x << endl;
    }
  }













  void TPBlockDifferentialOperator2 :: Apply (
            const FiniteElement & fel,
            const BaseMappedIntegrationRule & mir,
            BareSliceVector<double> x, 
            BareSliceMatrix<double> flux,
            LocalHeap & lh) const
  {
    const TPHighOrderFE & tpfel = static_cast<const TPHighOrderFE &>(fel);
    const TPMappedIntegrationRule & tpmir = static_cast<const TPMappedIntegrationRule &>(mir);
    auto & elements = tpfel.elements;
    auto & evalx = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(0);
    auto & evaly = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(1);
    auto & irs = tpmir.GetIRs();
    int ndof0 = elements[0]->GetNDof();
    int ndof1 = elements[1]->GetNDof();
    int dim0 = evalx->Dim();
    int dim1 = evaly->Dim();
    int nip0 = irs[0]->IR().Size();
    int nip1 = irs[1]->IR().Size();
    FlatMatrix<double, ColMajor> shape0( nip0*dim0, ndof0,lh );
    FlatMatrix<double, ColMajor> shape1( nip1*dim1, ndof1,lh );
    FlatMatrix<double> fcoefs( ndof0, ndof1*BlockDim(), &x(0) );
    evalx->CalcMatrix( *elements[0], *irs[0], shape0, lh );
    evaly->CalcMatrix( *elements[1], *irs[1], shape1, lh );
    flux.AddSize(nip0*nip1, dim0*dim1) = 0.0;
    FlatMatrix<> dim0geq1(nip1*dim1,nip0*dim0*BlockDim(),lh);
    dim0geq1 = 0.0;
    if(dim0 == 1)
    {
      for(int i : Range(BlockDim()) )
      {
        DoubleSliceMatrix<> fcomp(ndof0,ndof1,ndof1*BlockDim(),BlockDim(), &x(i));
        FlatMatrix<> fcomp_calc(ndof0,ndof1,lh);
        FlatMatrix<> fvals_calc(nip0*dim0,nip1*dim1,lh);
        DoubleSliceMatrix<double> fvals( nip0*dim0, nip1*dim1,BlockDim()*nip1*dim1,BlockDim(), &flux(i) );
        FlatMatrix<double,ColMajor> helper(ndof1,nip0*dim0,lh);
        // helper = Trans(fcomp) * Trans(shape0);
        // fvals = Trans(shape1 * helper);
        fcomp_calc = fcomp;
        helper = Trans(fcomp_calc) * Trans(shape0);
        fvals_calc = Trans(shape1 * helper);
        fvals = fvals_calc;
      }
    }
    if(dim0 > 1)
    {
      for(int i: Range(BlockDim()) )
      {
        DoubleSliceMatrix<> fcomp(ndof0,ndof1,ndof1*BlockDim(),BlockDim(), &x(i));
        DoubleSliceMatrix<double> fvals( nip1*dim1, nip0*dim0,BlockDim()*nip0*dim0,BlockDim(), &dim0geq1(0,i) );
        FlatMatrix<> fcomp_calc(ndof0, ndof1, lh);
        FlatMatrix<double> fvals_calc( nip1*dim1, nip0*dim0, lh);
        FlatMatrix<double,RowMajor> helper(ndof1,nip0*dim0,lh);
        // helper = Trans(fcomp) * Trans(shape0);
        // fvals = (shape1 * helper);
        fcomp_calc = fcomp;
        helper = Trans(fcomp_calc) * Trans(shape0);
        fvals_calc = (shape1 * helper);
        fvals = fvals_calc;
      }
      
      for(int i=0;i<nip0;i++)
        flux.Rows(i*nip1,(i+1)*nip1).AddSize(nip1,dim0*BlockDim()) = dim0geq1.Cols(i*BlockDim()*dim0,(i+1)*BlockDim()*dim0);
    }
  }


  void TPBlockDifferentialOperator2 :: ApplyTrans (
            const FiniteElement & fel,
            const BaseMappedIntegrationRule & mir,
            FlatMatrix<double> flux,
            BareSliceVector<double> x, 
            LocalHeap & lh) const
  {
    const TPHighOrderFE & tpfel = static_cast<const TPHighOrderFE &>(fel);
    const TPMappedIntegrationRule & tpmir = static_cast<const TPMappedIntegrationRule &>(mir);
    auto & elements = tpfel.elements;
    auto & evalx = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(0);
    auto & evaly = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(1);
    auto & irs = tpmir.GetIRs();
    int ndof0 = elements[0]->GetNDof();
    int ndof1 = elements[1]->GetNDof();
    int dim0 = evalx->Dim();
    int dim1 = evaly->Dim();
    int nip0 = irs[0]->IR().Size();
    int nip1 = irs[1]->IR().Size();
    FlatMatrix<double, ColMajor> shape0( nip0*dim0, ndof0,lh ); 
    FlatMatrix<double, ColMajor> shape1( nip1*dim1, ndof1,lh ); 
    evalx->CalcMatrix( *elements[0], *irs[0], shape0, lh ); 
    evaly->CalcMatrix( *elements[1], *irs[1], shape1, lh );
    // x = 0.0;
    if(dim0 == 1)
    {
      for(int i: Range(BlockDim()) )
      {
        DoubleSliceMatrix<> fcomp(ndof0,ndof1,ndof1*BlockDim(),BlockDim(), &x(i));
        DoubleSliceMatrix<double> fvals( nip0*dim0, nip1*dim1,BlockDim()*nip1*dim1,BlockDim(), &flux(i) );
        FlatMatrix<> fcomp_calc(ndof0,ndof1, lh);
        FlatMatrix<double> fvals_calc( nip0*dim0, nip1*dim1, lh);
        FlatMatrix<double> helper(nip0*dim0,ndof1,lh);
        fvals_calc = fvals;
        helper = fvals_calc*shape1;
        fcomp_calc = Trans(shape0) * helper;
        fcomp = fcomp_calc;
      }
    }
    else
    {
      for(int i: Range(BlockDim()) )
      {
        DoubleSliceMatrix<> fcomp(ndof0,ndof1,ndof1*BlockDim(),BlockDim(), &x(i));
        DoubleSliceMatrix<double> fvals( nip0*dim0, nip1*dim0,BlockDim()*nip1*dim0,BlockDim(), &flux(i) );
        FlatMatrix<double> helper(nip0*dim0,ndof1,lh);
        FlatMatrix<double> fvals_calc( nip0*dim0, nip1*dim1, lh );
        FlatMatrix<> fcomp_calc(ndof0,ndof1,lh);
        for(int i=0;i<nip1;i++)
          for(int j=0;j<nip0;j++)
            fvals_calc.Rows(dim0*j,dim0*(j+1)).Col(i) = (fvals.Cols(dim0*i,dim0*(i+1)).Row(j));
        helper = fvals_calc * shape1;
        fcomp_calc = Trans(shape0) * helper;
        fcomp = fcomp_calc;
      }
    }
  }

  void TPBlockDifferentialOperator2 :: ApplyX(
            const FiniteElement &fel,
            const BaseMappedIntegrationRule & mirx,
            FlatMatrix<double> flux,
            SliceMatrix<double> x,
            LocalHeap & lh) const
  {
    auto & evalx = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(0);
    auto & evaly = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(1);
    int dimx = evalx->Dim();
    int dimy = evaly->Dim();
    int nipx = mirx.IR().Size();
    FlatMatrix<double, ColMajor> bmatx( nipx*dimx, fel.GetNDof(),lh );
    evalx->CalcMatrix(fel,mirx,bmatx,lh);
    int nipy = x.Width()/dimy;
    flux = 0.0;
    FlatMatrix<> dim0geq1(nipx,dimx*x.Width()*BlockDim(),lh);
    if(dimx == 1)
    {
      for(int i: Range(BlockDim()) )
      {
        SliceMatrix<> xcomp(fel.GetNDof(),x.Width(),BlockDim()*x.Width(),&x(i,0));
        DoubleSliceMatrix<> resultmat(bmatx.Height(),xcomp.Width(),flux.Width()*nipy,flux.Width()/dimy, &flux(0,i));
        FlatMatrix<> resultmat_calc(bmatx.Height(),xcomp.Width(),lh);
        FlatMatrix<> xcomp_calc(fel.GetNDof(),x.Width(),lh);        
        xcomp_calc = xcomp;
        resultmat_calc = bmatx*xcomp_calc;
        resultmat = resultmat_calc;
        // resultmat = bmatx*xcomp;
      }
    }
    else
    {
      for(int i : Range(BlockDim()) )
      {
        SliceMatrix<> xcomp(fel.GetNDof(),x.Width(),BlockDim()*x.Width(),&x(i,0));
        DoubleSliceMatrix<double> resultmat(bmatx.Height(),x.Width(),BlockDim()*x.Width(),BlockDim(),&dim0geq1(0,i));
        FlatMatrix<double> resultmat_comp(bmatx.Height(),x.Width(),lh);
        resultmat_comp = bmatx*xcomp;
        resultmat = resultmat_comp;
      }
      for(int k=0;k<flux.Height();k++)
        flux.Rows(k*nipy,(k+1)*nipy) = (dim0geq1.Cols(dimx*k*BlockDim(),dimx*(k+1)*BlockDim()));
    }
  }

  void TPBlockDifferentialOperator2 :: ApplyXTrans(
            const FiniteElement &fel,
            const BaseMappedIntegrationRule & mirx,
            FlatMatrix<double> flux,
            SliceMatrix<double> x,
            LocalHeap & lh) const
  {
    auto & evalx = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(0);
    auto & evaly = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(1);
    int dimx = evalx->Dim();
    int dimy = evaly->Dim();
    int nipx = mirx.IR().Size();
    int nipy = x.Width()/dimy;
    FlatMatrix<double, ColMajor> bmatx( nipx*dimx, fel.GetNDof(),lh );
    evalx->CalcMatrix(fel,mirx,bmatx,lh);
    if(dimx == 1)
    {
      for(int i: Range(BlockDim()) )
      {
        SliceMatrix<> resultmat(fel.GetNDof(),x.Width(),BlockDim()*x.Width(),&x(i,0));
        DoubleSliceMatrix<> fluxcomp(nipx,x.Width(),BlockDim()*dimy*nipy, BlockDim(), &flux(0,i));
        FlatMatrix<> fluxcomp_calc(nipx,x.Width(),lh);
        FlatMatrix<> resultmat_calc(fel.GetNDof(),x.Width(),lh);
        fluxcomp_calc = fluxcomp;
        resultmat_calc = Trans(bmatx)*fluxcomp_calc;
        resultmat = resultmat_calc;
        // resultmat = Trans(bmatx)*fluxcomp;
      }
    }
    else
    {
      cout << "Also in here!! "<<endl;
      cout << "This branch may yield undesirable results!!!!"<<endl;
      FlatMatrix<double> proxyvaluesasmat( nipx, nipy*dimx, &flux(0,0) );
      FlatMatrix<double> proxyvaluesasmat1( nipx*dimx, nipy*dimy, lh );
      for(int i=0;i<nipy;i++)
        for(int j=0;j<nipx;j++)
          proxyvaluesasmat1.Rows(dimx*j,dimx*(j+1)).Col(i) = (proxyvaluesasmat.Cols(dimx*i,dimx*(i+1)).Row(j));
      x = proxyvaluesasmat1*bmatx;
    }
  }

  void TPBlockDifferentialOperator2 :: ApplyY(
            const FiniteElement &fel,
            const BaseMappedIntegrationRule & miry,
            FlatMatrix<double> flux,
            SliceMatrix<double> x,
            LocalHeap & lh) const
  {
    auto & evalx = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(0);
    auto & evaly = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(1);
    int dimx = evalx->Dim();
    int dimy = evaly->Dim();
    int nipy = miry.IR().Size();
    int ndofy = fel.GetNDof();
    FlatMatrix<double, ColMajor> bmaty( nipy*dimy, fel.GetNDof(),lh );
    evaly->CalcMatrix(fel,miry,bmaty,lh);
    flux = 0.0;
    FlatMatrix<> dim0geq1(nipy*dimy,x.Width()*BlockDim(),lh);
    dim0geq1 = 0.0;
    if(dimx == 1)
    {
      for(int i: Range(BlockDim()) )
      {
        SliceMatrix<> fcomp(ndofy,x.Width(),x.Width()*BlockDim(), &x(i,0));
        DoubleSliceMatrix<double> resultmat(fcomp.Width(),bmaty.Height(),BlockDim()*nipy*dimy,BlockDim(), &flux(0,i));
        FlatMatrix<double> resultmat_calc(fcomp.Width(),bmaty.Height(),lh);
        FlatMatrix<> fcomp_calc(x.Width(),ndofy,lh);
        fcomp_calc = Trans(fcomp);
        resultmat_calc = fcomp_calc*Trans(bmaty);
        resultmat = resultmat_calc;
        // resultmat_calc = Trans(bmaty * fcomp);
        // resultmat = resultmat_calc;
      }
      // FlatMatrix<double, RowMajor> resultmat_calc(bmaty.Height(),x.Width()*BlockDim(),lh);
      // FlatMatrix<double, RowMajor> fcomp_calc(ndofy,x.Width()*BlockDim(),lh);
      // for(int i: Range(BlockDim()) )
      // {
        // SliceMatrix<> fcomp(ndofy,x.Width(),x.Width()*BlockDim(), &x(i,0));
        // fcomp_calc.Cols(i*x.Width(),(i+1)*x.Width() ) = fcomp;
      // }
      // resultmat_calc = bmaty*fcomp_calc;
      // for(int i : Range(BlockDim()) )
      // {
        // DoubleSliceMatrix<double> resultmat(x.Width()*BlockDim(),bmaty.Height(),BlockDim()*nipy*dimy,BlockDim(), &flux(0,i));
        // resultmat = Trans(resultmat_calc).Rows(i*resultmat.Height(),(i+1)*resultmat.Height());
      // }
    }
    else
    {
      for(int i: Range(BlockDim()) )
      {
        SliceMatrix<> fcomp(ndofy,x.Width(),x.Width()*BlockDim(), &x(i,0));
        DoubleSliceMatrix<double> resultmat(bmaty.Height(),x.Width(),BlockDim()*x.Width(),BlockDim(),&dim0geq1(0,i));
        FlatMatrix<double> resultmat_calc(bmaty.Height(),x.Width(),lh);
        resultmat_calc =  bmaty * fcomp;
        resultmat = resultmat_calc;
      }
      for(int k=0;k<x.Width()/dimx;k++)
        flux.Rows(k*nipy,(k+1)*nipy) = (dim0geq1.Cols(dimx*k*BlockDim(),dimx*(k+1)*BlockDim()));
    }
  }

  void TPBlockDifferentialOperator2 :: ApplyYTrans(
            const FiniteElement &fel,
            const BaseMappedIntegrationRule & miry,
            FlatMatrix<double> flux,
            SliceMatrix<double> x,
            LocalHeap & lh) const
  {
    auto & evalx = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(0);
    auto & evaly = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(1);
    int dimx = evalx->Dim();
    int dimy = evaly->Dim();
    int nipy = miry.IR().Size();
    int nipx = x.Height()/dimx;
    FlatMatrix<double, ColMajor> bmaty( nipy*dimy, fel.GetNDof(),lh );
    evaly->CalcMatrix(fel,miry,bmaty,lh);
    if(dimx == 1)
    {
      for(int i : Range(BlockDim()) )
      {
        DoubleSliceMatrix<double> fvals( x.Height(), nipy*dimy,BlockDim()*nipy*dimy,BlockDim(), &flux(i) );
        DoubleSliceMatrix<> fcomp(x.Height(), fel.GetNDof(),x.Dist(),BlockDim(), &x(0,i));
        FlatMatrix<double> fvals_calc( x.Height(), nipy*dimy,lh );
        FlatMatrix<> fcomp_calc(x.Height(), fel.GetNDof(),lh);
        fvals_calc = fvals;
        fcomp_calc = fvals_calc*bmaty;
        fcomp = fcomp_calc;
        // fcomp = fvals*bmaty;
      }
      // FlatMatrix<double> fvals_calc( x.Height()*BlockDim(), nipy*dimy,lh );
      // FlatMatrix<> fcomp_calc(x.Height()*BlockDim(), fel.GetNDof(),lh);      
      // for(int i : Range(BlockDim()) )
      // {
        // DoubleSliceMatrix<double> fvals( x.Height(), nipy*dimy,BlockDim()*nipy*dimy,BlockDim(), &flux(i) );
        // fvals_calc.Rows(i*fvals.Height(),(i+1)*fvals.Height()) = fvals;
      // }
      // fcomp_calc = fvals_calc*bmaty;
      // for(int i : Range(BlockDim()) )
      // {
        // DoubleSliceMatrix<> fcomp(x.Height(), fel.GetNDof(),x.Dist(),BlockDim(), &x(0,i));
        // fcomp = fcomp_calc.Rows(i*fcomp.Height(),(i+1)*fcomp.Height());
      // }
    }
    else
    {
      FlatMatrix<> fvalstotal(nipx*dimx*BlockDim(), nipy*dimy,lh);
      for(int i=0;i<nipy;i++)
        for(int j=0;j<nipx;j++)
          fvalstotal.Rows(dimx*j*BlockDim(),dimx*(j+1)*BlockDim()).Col(i) = (flux.Cols(dimx*i*BlockDim(),dimx*(i+1)*BlockDim()).Row(j*nipy));
      for(int i: Range(BlockDim()) )
      {
        SliceMatrix<> fvals1(nipx*dimx,nipy*dimy,BlockDim()*nipy*dimy,&fvalstotal(i,0));
        DoubleSliceMatrix<> xcomp(nipx*dimx,fel.GetNDof(),x.Dist(),BlockDim(), &x(0,i));
        FlatMatrix<> xcomp_calc(nipx*dimx,fel.GetNDof(),lh);
        xcomp_calc = fvals1*bmaty;
        xcomp = xcomp_calc;
      }
    }
  }


}




