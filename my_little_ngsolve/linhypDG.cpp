
/*
  
Solver for the linear hyperbolic equation

du/dt  +  div (b u) = 0

by an explicit time-stepping method

*/



#include <solve.hpp>
// #include <compatibility.hpp>
using namespace ngsolve;
using ngfem::ELEMENT_TYPE;



template <int D>
class NumProcLinearHyperbolic : public NumProc
{
protected:
  shared_ptr<CoefficientFunction> cfflow;
  shared_ptr<GridFunction> gfu;

  double dt;
  double tend;

  Timer timer_element, timer_facet, timer_mass;

  class FacetData
  {
  public:
    int elnr[2];
    int facetnr[2];
    Vector<> flown;

    FacetData (int nip)
      : flown(nip) { ; }
  };

  class ElementData
  {
  public:
    MatrixFixWidth<D> flowip;
    Matrix<> invmass;

    ElementData (int ndof, int nip)
      : flowip(nip), invmass(ndof) { ; }
  };

  Array<FacetData*> facetdata;
  Array<ElementData*> elementdata;

    
public:
    
  NumProcLinearHyperbolic (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde), 
      timer_element("convection - time element"), 
      timer_facet("convection - time facet"),
      timer_mass("convection - time mass")
  {
    gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", "u"));
    cfflow = apde->GetCoefficientFunction (flags.GetStringFlag ("flow", "flow"));

    dt = flags.GetNumFlag ("dt", 0.001);
    tend = flags.GetNumFlag ("tend", 1);
  }



  virtual void Do(LocalHeap & lh)
  {
    cout << "solve conservation equation" << endl;



    // prepare ...

    const L2HighOrderFESpace & fes = 
      dynamic_cast<const L2HighOrderFESpace&> (*gfu->GetFESpace());


    int ne = ma->GetNE();
    int nf = ma->GetNFacets();

    elementdata.SetSize (ne);
    facetdata.SetSize (nf);

    ConstantCoefficientFunction one(1);
    MassIntegrator<D> bfi(&one);

    for (int i = 0; i < ne; i++)
      {
	HeapReset hr(lh);
	
	const DGFiniteElement<D> & fel = dynamic_cast<const DGFiniteElement<D>&> (fes.GetFE (i, lh));
	const IntegrationRule ir(fel.ElementType(), 2*fel.Order());

        const_cast<DGFiniteElement<D>&> (fel).PrecomputeShapes (ir);
        const_cast<DGFiniteElement<D>&> (fel).PrecomputeTrace ();


	MappedIntegrationRule<D,D> mir(ir, ma->GetTrafo (i, 0, lh), lh);
	
	elementdata[i] = new ElementData (fel.GetNDof(), ir.Size());
	ElementData & edi = *elementdata[i];

	cfflow -> Evaluate (mir, FlatMatrix<> (edi.flowip));
			    
	for (int j = 0; j < ir.Size(); j++)
	  {
	    Vec<D> flow = mir[j].GetJacobianInverse() * edi.flowip.Row(j);
	    flow *= ir[j].Weight() * mir[j].GetMeasure();		
	    edi.flowip.Row(j) = flow;
	  }


	FlatMatrix<> mass(fel.GetNDof(), lh);
	bfi.CalcElementMatrix (fel, ma->GetTrafo(i, 0, lh), mass, lh);
	CalcInverse (mass, edi.invmass);
      }



    Array<int> elnums, fnums, vnums;
    
    for (int i = 0; i < nf; i++)
      {
	HeapReset hr(lh);
	
	const DGFiniteElement<D-1> & felfacet = 
	  dynamic_cast<const DGFiniteElement<D-1>&> (fes.GetFacetFE (i, lh));
	IntegrationRule ir (felfacet.ElementType(), 2*felfacet.Order());
        const_cast<DGFiniteElement<D-1>&> (felfacet).PrecomputeShapes (ir);
	

	facetdata[i] = new FacetData (ir.Size());
	FacetData & fai = *facetdata[i];

	ma->GetFacetElements (i, elnums);
	
	fai.elnr[1] = -1;
	for (int j = 0; j < elnums.Size(); j++)
	  {
	    fai.elnr[j] = elnums[j];
	    ma->GetElFacets (elnums[j], fnums);
	    for (int k = 0; k < fnums.Size(); k++)
	      if (fnums[k] == i) fai.facetnr[j] = k;
	  }

	
	ELEMENT_TYPE eltype = ma->GetElType(elnums[0]);

	ma->GetElVertices (elnums[0], vnums);
	Facet2ElementTrafo transform(eltype, vnums); 
	FlatVec<D> normal_ref = ElementTopology::GetNormals(eltype) [fai.facetnr[0]];
	
	int nip = ir.Size();
	
	// transform facet coordinates to element coordinates
	IntegrationRule & irt = transform(fai.facetnr[0], ir, lh);  
	MappedIntegrationRule<D,D> mir(irt, ma->GetTrafo(elnums[0], 0, lh), lh);
	
	FlatMatrixFixWidth<D> flowir(nip, lh);
	
	cfflow -> Evaluate (mir, flowir);
	
	for (int j = 0; j < nip; j++)
	  {
	    Vec<D> normal = Trans (mir[j].GetJacobianInverse()) * normal_ref;       
	    
	    fai.flown(j) = InnerProduct (normal, flowir.Row(j));
	    fai.flown(j) *= ir[j].Weight() * mir[j].GetJacobiDet();
	  }
      }
    







    FlatVector<> vecu = gfu->GetVector().FVDouble();
    Vector<> conv(vecu.Size());
    Vector<> w(vecu.Size());
    Vector<> hu(vecu.Size());
    
#pragma omp parallel
    {
      LocalHeap lh2 = lh.Split();
      

      for (double t = 0; t < tend; t += dt)
        {

#pragma omp single
          cout << "\rt = " << setw(6) << t << flush;
          
          CalcConvection (vecu, conv, lh2);
          SolveM (conv, w, lh2);
          
#pragma omp single
          {
            hu = vecu + (0.5*dt) * w;
          }

          CalcConvection (hu, conv, lh2);
          SolveM (conv, w, lh2);
          
#pragma omp single
          {
            vecu += dt * w;
            
            /*
              cout << " time T/F/M [us] = "
              << 1e6 * timer_element.GetTime()/timer_element.GetCounts()/vecu.Size() << " / "
              << 1e6 * timer_facet.GetTime()/timer_facet.GetCounts()/vecu.Size() << " / "
              << 1e6 * timer_mass.GetTime()/timer_mass.GetCounts()/vecu.Size() 
              << "\r";
            */
            Ng_Redraw();
          }
        }
    }
  }




  void SolveM (FlatVector<double> res, FlatVector<double> vecu,
	       LocalHeap & lh)
  {
    int ne = ma->GetNE();
    const L2HighOrderFESpace & fes = 
      dynamic_cast<const L2HighOrderFESpace&> (*gfu->GetFESpace());

#pragma omp single
    timer_mass.Start();

#pragma omp for
    for (int i = 0; i < ne; i++)
      {
	IntRange dn = fes.GetElementDofs (i);
	vecu.Range (dn) = elementdata[i]->invmass * res.Range (dn);
      }

#pragma omp single
    timer_mass.Stop();
  }






  void CalcConvection (FlatVector<double> vecu, FlatVector<double> conv,
		       LocalHeap & lh)
  {
    
    const L2HighOrderFESpace & fes = 
      dynamic_cast<const L2HighOrderFESpace&> (*gfu->GetFESpace());
    
    
#pragma omp single
    timer_element.Start();
    
    
    int ne = ma->GetNE();
      
#pragma omp for
      for (int i = 0; i < ne; i++)
	{
	  HeapReset hr(lh);
	  
	  const ScalarFiniteElement<D> & fel = 
            static_cast<const ScalarFiniteElement<D>&> (fes.GetFE (i, lh));
	  const IntegrationRule ir(fel.ElementType(), 2*fel.Order());

	  FlatMatrixFixWidth<D> flowip = elementdata[i]->flowip;

	  /*
	  // use this for time-dependent flow
	  MappedIntegrationRule<D,D> mir(ir, ma->GetTrafo (i, 0, lh), lh);
	  FlatMatrixFixWidth<D> flowip(mir.Size(), lh);
	  cfflow -> Evaluate (mir, flowip);
	  for (int j = 0; j < ir.Size(); j++)
	    {
	      Vec<D> flow = mir[j].GetJacobianInverse() * flowip.Row(j);
	      flow *= mir[j].GetWeight();		
	      flowip.Row(j) = flow;
	    }
	  */

	  IntRange dn = fes.GetElementDofs (i);
	  
	  int nipt = ir.Size();
	  FlatVector<> elui(nipt, lh);
	  FlatMatrixFixWidth<D> flowui (nipt, lh);
	  
	  fel.Evaluate (ir, vecu.Range (dn), elui);
	  
	  flowui = flowip;
	  for (int k = 0; k < nipt; k++)
	    flowui.Row(k) *= elui(k);
	  
	  fel.EvaluateGradTrans (ir, flowui, conv.Range(dn));
	}


#pragma omp single
      {
        timer_element.Stop();
        timer_facet.Start();
      }


      int nf = ma->GetNFacets();
      
#pragma omp for
      for (int i = 0; i < nf; i++)
	{
	  HeapReset hr(lh);
	  
	  const FacetData & fai = *facetdata[i];
	  if (fai.elnr[1] != -1)
	    {
	      const DGFiniteElement<D> & fel1 = 
		static_cast<const DGFiniteElement<D>&> (fes.GetFE (fai.elnr[0], lh));
	      const DGFiniteElement<D> & fel2 = 
		static_cast<const DGFiniteElement<D>&> (fes.GetFE (fai.elnr[1], lh));
	      const DGFiniteElement<D-1> & felfacet = 
		static_cast<const DGFiniteElement<D-1>&> (fes.GetFacetFE (i, lh));

	      IntRange dn1 = fes.GetElementDofs (fai.elnr[0]);
	      IntRange dn2 = fes.GetElementDofs (fai.elnr[1]);

	      int ndoffacet = felfacet.GetNDof();
	      int ndof1 = fel1.GetNDof();
	      int ndof2 = fel2.GetNDof();

	      FlatVector<> aelu1(ndof1, lh), aelu2(ndof2, lh);
	      FlatVector<> trace1(ndoffacet, lh), trace2(ndoffacet, lh);

	      fel1.GetTrace (fai.facetnr[0], vecu.Range (dn1), trace1);
	      fel2.GetTrace (fai.facetnr[1], vecu.Range (dn2), trace2);

	      IntegrationRule ir(felfacet.ElementType(), 2*felfacet.Order());
	      int nip = ir.Size();

	      FlatVector<> flown = fai.flown;
	    
	      FlatVector<> tracei1(nip, lh), tracei2(nip, lh);
	      FlatVector<> tracei(nip, lh);

	      felfacet.Evaluate (ir, trace1, tracei1);
	      felfacet.Evaluate (ir, trace2, tracei2);
		    
	      for (int j = 0; j < nip; j++)
		tracei(j) = flown(j) * ( (flown(j) > 0) ? tracei1(j) : tracei2(j) );

	      felfacet.EvaluateTrans (ir, tracei, trace1);
	      fel1.GetTraceTrans (fai.facetnr[0], trace1, aelu1);
	      fel2.GetTraceTrans (fai.facetnr[1], trace1, aelu2);

#pragma omp critical (addres)
	      {
		conv.Range (dn1) -= aelu1;
		conv.Range (dn2) += aelu2;
	      }
	    }
	  else
	    {
	      const DGFiniteElement<D> & fel1 = 
		dynamic_cast<const DGFiniteElement<D>&> (fes.GetFE (fai.elnr[0], lh));
	      const DGFiniteElement<D-1> & felfacet = 
		dynamic_cast<const DGFiniteElement<D-1>&> (fes.GetFacetFE (i, lh));

	      IntRange dn1 = fes.GetElementDofs (fai.elnr[0]);

	      int ndoffacet = felfacet.GetNDof();
	      int ndof1 = fel1.GetNDof();

	      FlatVector<> elu1(ndof1, lh);
	      FlatVector<> trace1(ndoffacet, lh);

	      fel1.GetTrace (fai.facetnr[0], vecu.Range (dn1), trace1);

	      IntegrationRule ir(felfacet.ElementType(), 2*felfacet.Order());
	      int nip = ir.Size();

	      FlatVector<> flown = fai.flown; 
	      FlatVector<> tracei1(nip, lh), tracei(nip, lh);

	      felfacet.Evaluate (ir, trace1, tracei1);
		    
	      for (int j = 0; j < nip; j++)
		tracei(j) = flown(j) * ( (flown(j) > 0) ? tracei1(j) : 0 );

	      felfacet.EvaluateTrans (ir, tracei, trace1);
	      fel1.GetTraceTrans (fai.facetnr[0], trace1, elu1);
	    
#pragma omp critical (addres)
	      {
		conv.Range (dn1) -= elu1;
	      }
	    }
	}

#pragma omp single    
      timer_facet.Stop(); 
      


  }
};




static RegisterNumProc<NumProcLinearHyperbolic<2> > npinit1("linhyp", 2);
static RegisterNumProc<NumProcLinearHyperbolic<3> > npinit2("linhyp", 3);

  










