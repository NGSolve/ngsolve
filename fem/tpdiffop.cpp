#include <fem.hpp>
#include "tpintrule.hpp"
#include "tpdiffop.hpp"

    
namespace ngfem 
{

  void ProlongateCoefficientFunction :: Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const
  {
    const TPMappedIntegrationRule * tpmir = static_cast<const TPMappedIntegrationRule *>(&ir);
    auto & irs = tpmir->GetIRs();
    coef->Evaluate(*irs[1-prolongateto],values);
    if(prolongateto == 1)
      for(int i=irs[0]->Size()-1;i>=0;i--)
        values.Rows(i*irs[1]->Size(),(i+1)*irs[1]->Size()) = values.Row(i)(0);
    if(prolongateto == 0)
      for(int i=1;i<irs[0]->Size();i++)
        values.Rows(i*irs[1]->Size(),(i+1)*irs[1]->Size()) = values.Rows(0,irs[1]->Size());
  }

  void TPDifferentialOperator :: Apply (
            const FiniteElement & fel,
            const BaseMappedIntegrationRule & mir,
            FlatVector<double> x, 
            FlatMatrix<double> flux,
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
    // cout << "In = "<< endl << fcoefs << endl;
    if(dim0 == 1)
    {
      FlatMatrix<double,ColMajor> fvals( nip1*dim1, nip0*dim0, &flux(0) );
      FlatMatrix<double,ColMajor> helper(ndof1,nip0*dim0,lh);
      helper = Trans(fcoefs) * Trans(shape0); // (a*b'*c')' = c*(a*b')' = c*b*a'
      fvals = shape1 * helper;
      // cout << flux << endl;
    }
    if(dim0 > 1)
    {
      FlatMatrix<double,RowMajor> fvals( nip1*dim1, nip0*dim0, lh );
      FlatMatrix<double,RowMajor> helper(ndof1,nip0*dim0,lh);
      helper = Trans(fcoefs) * Trans(shape0);
      fvals = shape1 * helper;
      for(int i=0;i<nip0;i++)
        flux.Rows(i*nip1,(i+1)*nip1) = fvals.Cols(dim0*i,dim0*(i+1));
      // cout << "Flux = "<<endl<<flux << endl;
    }
  }
    

  void TPDifferentialOperator :: ApplyTrans (
            const FiniteElement & fel,
            const BaseMappedIntegrationRule & mir,
            FlatMatrix<double> flux,
            FlatVector<double> x, 
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
    // cout << "TPDifferentialOperator::ApplyTrans"<<endl;
    // cout << "Dims = ("<<dim0<<", "<<dim1<<")"<<endl;
    // cout << "In = "<<endl<<flux << endl;
    if(dim0 == 1)
    {
      FlatMatrix<double> fvals( nip0*dim0, nip1*dim1, &flux(0,0) );
      // cout << "Flux as mat = "<<endl<<fvals << endl;
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
      // cout << "fcoef mat = "<<endl<<fcoefs<< endl;
    }
    // cout << "Out = "<<endl<<x << endl;
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
        flux.Rows(k*nipx,(k+1)*nipx) = Trans(resultmat.Rows(dimx*k,dimx*(k+1)));
      // cout << "TPDifferentialOperator::ApplyX"<< endl;
      // cout << "x = "   << endl << x    << endl;
      // cout << "flux = "<< endl << flux << endl;        
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
      cout << "resultmat = "<< endl<<resultmat << endl;
    }
    else
    {
      // cout << "TPDifferentialOperator::ApplyY" << endl;
      FlatMatrix<double, ColMajor> resultmat(x.Height(),nipy*dimy, lh);
      resultmat = x*Trans(bmaty);
      for(int k=0;k<x.Height()/dimx;k++)
        flux.Rows(k*nipy,(k+1)*nipy) = Trans(resultmat.Rows(dimx*k,dimx*(k+1)));
      // cout << "x = "<<endl<< x << endl;
      cout << "resultmat = "<< endl<<resultmat << endl;
      // cout << "flux = "<<endl<< flux << endl;
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
      x = proxyvaluesasmat*bmaty;
      // cout << "Input as Matrix = "<<endl<<proxyvaluesasmat<<endl;
    }
    else
    {
      FlatMatrix<double> proxyvaluesasmat( nipx, nipy*dimx, &flux(0,0) );
      FlatMatrix<double> proxyvaluesasmat1( nipx*dimx, nipy*dimy, lh );
      for(int i=0;i<nipy;i++)
        for(int j=0;j<nipx;j++)
          proxyvaluesasmat1.Rows(dimx*j,dimx*(j+1)).Col(i) = (proxyvaluesasmat.Cols(dimx*i,dimx*(i+1)).Row(j));
      x = proxyvaluesasmat1*bmaty;
      // cout << "nipx = "<<nipx<<" nipy = "<<nipy << " dimx = "<<dimx << " dimy = "<<dimy << endl;
      // cout << "TPBlockDifferentialOperator::ApplyYTrans" << endl;
      // cout << "flux (in) = "<<endl<<flux<<endl;
      // cout << "fluxasmat = "<<endl<<proxyvaluesasmat<<endl;
      // cout << "fluxasmat = "<<endl<<proxyvaluesasmat1<<endl;
      // cout << "x = "<< endl<<x << endl;      
    }
  }


  void TPBlockDifferentialOperator :: Apply (
            const FiniteElement & fel,
            const BaseMappedIntegrationRule & mir,
            FlatVector<double> x, 
            FlatMatrix<double> flux,
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
      FlatMatrix<double,ColMajor> fluxCM(flux.Height(),flux.Width(),&result(0,0));
      helper1 = shape0 * fcoefs;
      result =  shape1*helper2;
      flux = fluxCM;
    }
    if(dim0 > 1)
    {
      // FlatMatrix<double,RowMajor> fvals( nip1*dim1, nip0*dim0*BlockDim(), lh );
      // FlatMatrix<double,RowMajor> helper(ndof1,nip0*dim0,lh);
      // helper = Trans(fcoefs) * Trans(shape0);
      // fvals = shape1 * helper;
      // for(int i=0;i<nip0;i++)
        // flux.Rows(i*nip1,(i+1)*nip1) = fvals.Cols(dim0*i,dim0*(i+1));

      FlatMatrix<Vec<3>,RowMajor> fvals( nip1*dim1, nip0*dim0, lh );
      FlatMatrix<double,RowMajor> fvalsdb( nip1*dim1, nip0*dim0*BlockDim(), reinterpret_cast<double * >(&fvals(0,0)) );
      FlatMatrix<Vec<3>,RowMajor> helper(ndof1,nip0*dim0,lh);
      FlatMatrix<Vec<3>> fcoefs1(ndof0,ndof1,reinterpret_cast<Vec<3> *>(&fcoefs(0,0)));
      helper = fcoefs1 * Trans(shape0);
      fvals = shape1 * helper;
      for(int i=0;i<nip0;i++)
        flux.Rows(i*nip1,(i+1)*nip1) = fvalsdb.Cols(dim0*i*BlockDim(),dim0*(i+1)*BlockDim());
    }      
  }
    

  void TPBlockDifferentialOperator :: ApplyTrans (
            const FiniteElement & fel,
            const BaseMappedIntegrationRule & mir,
            FlatMatrix<double> flux,
            FlatVector<double> x, 
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
      FlatMatrix<Vec<3> > fvalsvec( nip0*dim0, nip1*dim1, reinterpret_cast<Vec<3> * >(&flux(0,0)) );
      FlatMatrix<Vec<3> > fcoefsvec( ndof0, ndof1, reinterpret_cast<Vec<3> * >(&x(0,0)) );
      FlatMatrix<Vec<3> > helpervec(nip0*dim0,ndof1*BlockDim(),lh);
      helpervec = fvalsvec*shape1;
      fcoefsvec = Trans(shape0)*helpervec;
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
    auto & evaly = static_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(1);
    int dimx = evalx->Dim();
    int dimy = evaly->Dim();
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
      FlatMatrix<Vec<3>, ColMajor> resultmat(x.Height(),nipy*dimy, lh);
      FlatMatrix<double, ColMajor> resultmat_dbl(x.Height(),BlockDim()*nipy*dimy, reinterpret_cast<double *>(&resultmat(0,0)));
      SliceMatrix< Vec<3> > x_vec(x.Height(),0.5*x.Width(),0.5*x.Dist(),reinterpret_cast<Vec<3> * >(&x(0,0)));
      resultmat = x_vec*Trans(bmaty);
      for(int k=0;k<x.Height()/dimx;k++)
        flux.Rows(k*nipy,(k+1)*nipy) = Trans(resultmat_dbl.Rows(dimx*k*BlockDim(),dimx*(k+1)*BlockDim()));
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
      FlatMatrix<Vec<3> > proxyvaluesasmat( nipx, nipy*dimx, reinterpret_cast<Vec<3> *>(&flux(0,0)));
      FlatMatrix<Vec<3> > proxyvaluesasmat1( nipx*dimx, nipy*dimy, lh );
      for(int i=0;i<nipy;i++)
        for(int j=0;j<nipx;j++)
      proxyvaluesasmat1.Rows(dimx*j,dimx*(j+1)).Col(i) = (proxyvaluesasmat.Cols(dimx*i,dimx*(i+1)).Row(j));
      SliceMatrix<Vec<3> > x_vec(x.Height(), 1.0/3.0*x.Width(),1.0/3.0*x.Dist(), reinterpret_cast<Vec<3> *>(&x(0,0)));
      x_vec = proxyvaluesasmat1*bmaty;
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
            FlatVector<double> x, 
            FlatMatrix<double> flux,
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
    flux = 0.0;
    FlatMatrix<> dim0geq1(nip1*dim1,nip0*dim0*BlockDim(),lh);
    dim0geq1 = 0.0;
    if(dim0 == 1)
    {
      for(int i : Range(BlockDim()) )
      {
        DoubleSliceMatrix<> fcomp(ndof0,ndof1,ndof1*BlockDim(),BlockDim(), &x(i));
        DoubleSliceMatrix<double> fvals( nip0*dim0, nip1*dim1,BlockDim()*nip1*dim1,BlockDim(), &flux(i) );
        FlatMatrix<double,ColMajor> helper(ndof1,nip0*dim0,lh);
        helper = Trans(fcomp) * Trans(shape0);
        fvals = Trans(shape1 * helper);
      }
    }
    if(dim0 > 1)
    {
      for(int i: Range(BlockDim()) )
      {
        DoubleSliceMatrix<> fcomp(ndof0,ndof1,ndof1*BlockDim(),BlockDim(), &x(i));
        DoubleSliceMatrix<double> fvals( nip1*dim1, nip0*dim0,BlockDim()*nip0*dim0,BlockDim(), &dim0geq1(0,i) );
        FlatMatrix<double,RowMajor> helper(ndof1,nip0*dim0,lh);
        helper = Trans(fcomp) * Trans(shape0);
        fvals = (shape1 * helper);
      }
      for(int i=0;i<nip0;i++)
        flux.Rows(i*nip1,(i+1)*nip1) = dim0geq1.Cols(i*BlockDim()*dim0,(i+1)*BlockDim()*dim0);
    }
  }


  void TPBlockDifferentialOperator2 :: ApplyTrans (
            const FiniteElement & fel,
            const BaseMappedIntegrationRule & mir,
            FlatMatrix<double> flux,
            FlatVector<double> x, 
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
    // cout << "TPBlockDifferentialOperator2::ApplyTrans"<<endl;
    // cout << "Dims = ("<<dim0<<", "<<dim1<<")"<<endl;
    // cout << "In = "<<endl<< flux << endl;
    x = 0.0;
    for(int i: Range(BlockDim()) )
    {
      if(dim0 == 1)
      {
        DoubleSliceMatrix<> fcomp(ndof0,ndof1,ndof1*BlockDim(),BlockDim(), &x(i));
        DoubleSliceMatrix<double> fvals( nip0*dim0, nip1*dim1,BlockDim()*nip1*dim1,BlockDim(), &flux(i) );
        FlatMatrix<double> helper(nip0*dim0,ndof1,lh);
        helper = fvals*shape1;
        fcomp = Trans(shape0) * helper;
        // FlatMatrix<Vec<3> > fvalsvec( nip0*dim0, nip1*dim1, reinterpret_cast<Vec<3> * >(&flux(0,0)) );
        // FlatMatrix<Vec<3> > fcoefsvec( ndof0, ndof1, reinterpret_cast<Vec<3> * >(&x(0,0)) );
        // FlatMatrix<Vec<3> > helpervec(nip0*dim0,ndof1*BlockDim(),lh);
        // helpervec = fvalsvec*shape1;
        // fcoefsvec = Trans(shape0)*helpervec;
      }
      else
      {
        DoubleSliceMatrix<> fcomp(ndof0,ndof1,ndof1*BlockDim(),BlockDim(), &x(i));
        DoubleSliceMatrix<double> fvals( nip0*dim0, nip1*dim0,BlockDim()*nip1*dim0,BlockDim(), &flux(i) );
        FlatMatrix<double> helper(nip0*dim0,ndof1,lh);
        FlatMatrix<double> fvals1( nip0*dim0, nip1*dim1, lh );
        for(int i=0;i<nip1;i++)
          for(int j=0;j<nip0;j++)
            fvals1.Rows(dim0*j,dim0*(j+1)).Col(i) = (fvals.Cols(dim0*i,dim0*(i+1)).Row(j));
        
        helper = fvals1 * shape1;
        fcomp = Trans(shape0) * helper;
      }
      // cout << "fcoefs = "<<endl<<x << endl;
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
    FlatMatrix<> dim0geq1(nipy*BlockDim(),flux.Width(),lh);
    dim0geq1 = 0.0;
    
    if(dimx == 1)
    {
      for(int i: Range(BlockDim()) )
      {
        SliceMatrix<> fcomp(ndofy,x.Width(),x.Width()*BlockDim(), &x(i,0));
        DoubleSliceMatrix<double> resultmat(fcomp.Width(),bmaty.Height(),BlockDim()*BlockDim()*dimy,BlockDim(), &flux(0,i));
        resultmat = Trans(bmaty * fcomp);
      }
    }
    else
    {
      for(int i: Range(BlockDim()) )
      {
        SliceMatrix<> fcomp(ndofy,x.Width(),x.Width()*BlockDim(), &x(i,0));
        DoubleSliceMatrix<double> resultmat(bmaty.Height(),fcomp.Width(),BlockDim()*dimx*3,BlockDim(),&dim0geq1(0,i));
        FlatMatrix<double> resultmat1(bmaty.Height(),fcomp.Width(),lh);
        resultmat =  (bmaty * fcomp);
        resultmat1 =  (bmaty * fcomp);
        cout << "resultmat (DoubleSliceMat) = "<<endl<<resultmat << endl;
        cout << "resultmat (FlatMat) = "<<endl<<resultmat1 << endl;
        cout << "dim0geq1 = "<<endl<<dim0geq1 << endl;
        // FlatMatrix<Vec<3>, ColMajor> resultmat(x.Height(),nipy*dimy, lh);
        // FlatMatrix<double, ColMajor> resultmat_dbl(x.Height(),BlockDim()*nipy*dimy, reinterpret_cast<double *>(&resultmat(0,0)));
        // SliceMatrix< Vec<3> > x_vec(x.Height(),0.5*x.Width(),0.5*x.Dist(),reinterpret_cast<Vec<3> * >(&x(0,0)));
        // resultmat = x_vec*Trans(bmaty);
      }
      for(int k=0;k<x.Height()/dimx;k++)
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
      FlatMatrix<Vec<3> > proxyvaluesasmat( nipx, nipy*dimx, reinterpret_cast<Vec<3> *>(&flux(0,0)));
      FlatMatrix<Vec<3> > proxyvaluesasmat1( nipx*dimx, nipy*dimy, lh );
      for(int i=0;i<nipy;i++)
        for(int j=0;j<nipx;j++)
      proxyvaluesasmat1.Rows(dimx*j,dimx*(j+1)).Col(i) = (proxyvaluesasmat.Cols(dimx*i,dimx*(i+1)).Row(j));
      SliceMatrix<Vec<3> > x_vec(x.Height(), 1.0/3.0*x.Width(),1.0/3.0*x.Dist(), reinterpret_cast<Vec<3> *>(&x(0,0)));
      x_vec = proxyvaluesasmat1*bmaty;
    }
  }

 
}




