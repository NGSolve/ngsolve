/*
  
Solver for the linear hyperbolic equation

du/dt  +  div (b u) = 0

by an explicit time-stepping method

*/


#include <solve.hpp>
using namespace ngsolve;


template <int D>
class NumProcLinearHyperbolic : public NumProc
{
protected:
  CoefficientFunction * cfflow;
  GridFunction * gfu;

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
    Array<int> facets;
    ElementData (int ndof, int nip, int nf)
      : flowip(nip), invmass(ndof), facets(nf) { ; }
  };

  Array<FacetData*> facetdata;
  Array<ElementData*> elementdata;

  int ndf;  // num dofs per facet
  ParallelDofs * pardofs;
    
public:
    
  NumProcLinearHyperbolic (PDE & apde, const Flags & flags)
    : NumProc (apde), 
      timer_element("conservation - time element"), 
      timer_facet("conservation - time facet"),
      timer_mass("conservation - time mass")
  {
    gfu = pde.GetGridFunction (flags.GetStringFlag ("gridfunction", "u"));
    cfflow = pde.GetCoefficientFunction (flags.GetStringFlag ("flow", "flow"));

    dt = flags.GetNumFlag ("dt", 0.001);
    tend = flags.GetNumFlag ("tend", 1);
  }



  virtual void Do(LocalHeap & lh)
  {
    cout << "solve conservation equation" << endl;



    // prepare ...

    const L2HighOrderFESpace & fes = 
      dynamic_cast<const L2HighOrderFESpace&> (gfu->GetFESpace());


    int ne = ma.GetNE();
    int nf = ma.GetNFacets();

    elementdata.SetSize (ne);
    facetdata.SetSize (nf);

    ConstantCoefficientFunction one(1);
    MassIntegrator<D> bfi(&one);

    for (int i = 0; i < ne; i++)
      {
	HeapReset hr(lh);
	
	const DGFiniteElement<D> & fel = dynamic_cast<const DGFiniteElement<D>&> (fes.GetFE (i, lh));
	const IntegrationRuleTP<D> ir(ma.GetTrafo(i,0,lh), 2*fel.Order(), &lh);

	const_cast<DGFiniteElement<D>&> (fel).PrecomputeShapes (ir);
	const_cast<DGFiniteElement<D>&> (fel).PrecomputeTrace ();

	Array<int> facets;
	ma.GetElFacets (i, facets);
	MappedIntegrationRule<D,D> mir(ir, ma.GetTrafo (i, 0, lh), lh);
	
	elementdata[i] = new ElementData (fel.GetNDof(), ir.Size(), nf);
	ElementData & edi = *elementdata[i];
	edi.facets = facets;
	cfflow -> Evaluate (mir, FlatMatrix<> (edi.flowip));
			    
	for (int j = 0; j < ir.Size(); j++)
	  {
	    Vec<D> flow = mir[j].GetJacobianInverse() * edi.flowip.Row(j);
	    edi.flowip.Row(j) = mir[j].GetWeight() * flow;
	  }

	FlatMatrix<> mass(fel.GetNDof(), lh);
	bfi.CalcElementMatrix (fel, ma.GetTrafo(i, 0, lh), mass, lh);
	CalcInverse (mass, edi.invmass);
      }



    Array<int> elnums, fnums, vnums;

    MyMPI_Barrier();
    cout << "id = " << MyMPI_GetId() << ", nf = " << nf << endl;
    MyMPI_Barrier();

    for (int i = 0; i < nf; i++)
      {
	cout << "id = " << MyMPI_GetId() << ", facet nr " << i << endl;
	HeapReset hr(lh);
	
	const DGFiniteElement<D-1> & felfacet = 
	  dynamic_cast<const DGFiniteElement<D-1>&> (fes.GetFacetFE (i, lh));
	IntegrationRule ir(felfacet.ElementType(), 2*felfacet.Order());
	const_cast<DGFiniteElement<D-1>&> (felfacet).PrecomputeShapes (ir);
	

	facetdata[i] = new FacetData (ir.Size());
	FacetData & fai = *facetdata[i];

	ma.GetFacetElements (i, elnums);
	
	fai.elnr[1] = -1;
	for (int j = 0; j < elnums.Size(); j++)
	  {
	    fai.elnr[j] = elnums[j];
	    ma.GetElFacets (elnums[j], fnums);
	    for (int k = 0; k < fnums.Size(); k++)
	      if (fnums[k] == i) fai.facetnr[j] = k;
	  }
	
	const DGFiniteElement<D> & fel = 
	  dynamic_cast<const DGFiniteElement<D>&> (fes.GetFE (elnums[0], lh));
	
	ma.GetElVertices (elnums[0], vnums);
	Facet2ElementTrafo transform(fel.ElementType(), vnums); 
	FlatVec<D> normal_ref = 
	  ElementTopology::GetNormals(fel.ElementType())[fai.facetnr[0]];
	
	int nip = ir.Size();

	// transform facet coordinates to element coordinates
	IntegrationRule & irt = transform(fai.facetnr[0], ir, lh);
	MappedIntegrationRule<D,D> mir(irt, ma.GetTrafo(elnums[0], 0, lh), lh);
	
	FlatVector<> flown = fai.flown;
	FlatMatrix<> flowir(nip, D, lh);
	
	cfflow -> Evaluate (mir, flowir);
	
	for (int j = 0; j < nip; j++)
	  {
	    Vec<D> normal = Trans (mir[j].GetJacobianInverse()) * normal_ref;       
	    
	    flown(j) = InnerProduct (normal, flowir.Row(j));
	    flown(j) *= ir[j].Weight() * mir[j].GetMeasure(); 
	  }
      }
    

    MyMPI_Barrier();
    cout << "alle hier, A" << " id = " << MyMPI_GetId() << endl;
    MyMPI_Barrier();


    if (ne)
      {
	ndf = fes.GetFacetFE (0, lh).GetNDof();
	
	Array<Node> dofnodes(ndf*nf);
	for (int i = 0; i < nf; i++)
	  for (int j = 0; j < ndf; j++)
	    dofnodes[i*ndf+j] = Node(NODE_TYPE(D-1), i);
    
	pardofs = new ParallelMeshDofs (ma, dofnodes);
      }
    else
      {
	Array<Node> dofnodes(0);
	pardofs = new ParallelMeshDofs (ma, dofnodes);
      }

    MyMPI_Barrier();
    cout << "alle hier, B" << " id = " << MyMPI_GetId() << endl;
    MyMPI_Barrier();



    FlatVector<> vecu = gfu->GetVector().FVDouble();
    Vector<> conv(vecu.Size());
    Vector<> w(vecu.Size());
    Vector<> hu(vecu.Size());
    

    for (double t = 0; t < tend; t += dt)
      {
	cout << IM(1) << "t = " << setw(6) << t << flush;

	CalcConvection (vecu, conv, lh);
	SolveM (conv, w, lh);

	hu = vecu + (0.5*dt) * w;
	CalcConvection (hu, conv, lh);
	SolveM (conv, w, lh);

	vecu += dt * w;


	cout << " time T/F/M [us] = "
	     << 1e6 * timer_element.GetTime()/timer_element.GetCounts()/vecu.Size() << " / "
	     << 1e6 * timer_facet.GetTime()/timer_facet.GetCounts()/vecu.Size() << " / "
	     << 1e6 * timer_mass.GetTime()/timer_mass.GetCounts()/vecu.Size() 
	     << "\r";

	MyMPI_Barrier();
	Ng_Redraw();
      }
  }




  void SolveM (FlatVector<double> res, FlatVector<double> vecu,
	       LocalHeap & lh)
  {
    int ne = ma.GetNE();
    const L2HighOrderFESpace & fes = 
      dynamic_cast<const L2HighOrderFESpace&> (gfu->GetFESpace());

    timer_mass.Start();

    for (int i = 0; i < ne; i++)
      {
	IntRange dn = fes.GetElementDofs (i);
	vecu.Range (dn) = elementdata[i]->invmass * res.Range (dn);
      }

    timer_mass.Stop();
  }






  void CalcConvection (FlatVector<double> vecu, FlatVector<double> conv,
		       LocalHeap & lh)
  {
    const L2HighOrderFESpace & fes = 
      dynamic_cast<const L2HighOrderFESpace&> (gfu->GetFESpace());

    int ne = ma.GetNE();
    int nf = ma.GetNFacets();
    

    timer_element.Start();
    
#pragma omp parallel 
    {
      LocalHeap & clh = lh, lh = clh.Split();
#pragma omp for
      for (int i = 0; i < ne; i++)
	{
	  HeapReset hr(lh);
	  
	  const ScalarFiniteElement<D> & fel = dynamic_cast<const ScalarFiniteElement<D>&> (fes.GetFE (i, lh));
	  const IntegrationRuleTP<D> ir( ma.GetTrafo(i,0,lh), 2*fel.Order(), &lh);
	  
	  IntRange dn = fes.GetElementDofs (i);
	  
	  int nipt = ir.Size();
	  FlatVector<> elui(nipt, lh);
	  FlatMatrixFixWidth<D> flowui (nipt, lh);

	  fel.Evaluate (ir, vecu.Range (dn), elui);
	  
	  flowui = elementdata[i]->flowip;
	  for (int k = 0; k < nipt; k++)
	    flowui.Row(k) *= elui(k);
	  
	  fel.EvaluateGradTrans (ir, flowui, conv.Range(dn));
	}
    }

    timer_element.Stop();

    

    ParallelVVector<> traces (pardofs->GetNDofLocal(), pardofs);
    traces.SetStatus (DISTRIBUTED);

    timer_facet.Start();


#pragma omp parallel 
    {
      LocalHeap & clh = lh, lh = clh.Split();
#pragma omp for
      for (int i = 0; i < nf; i++)
	{
	  HeapReset hr(lh);
	  
	  const FacetData & fai = *facetdata[i];
	  if (fai.elnr[1] != -1)
	    {
	      const DGFiniteElement<D> & fel1 = 
		dynamic_cast<const DGFiniteElement<D>&> (fes.GetFE (fai.elnr[0], lh));
	      const DGFiniteElement<D> & fel2 = 
		dynamic_cast<const DGFiniteElement<D>&> (fes.GetFE (fai.elnr[1], lh));
	      const DGFiniteElement<D-1> & felfacet = 
		dynamic_cast<const DGFiniteElement<D-1>&> (fes.GetFacetFE (i, lh));

	      IntRange dn1 = fes.GetElementDofs (fai.elnr[0]);
	      IntRange dn2 = fes.GetElementDofs (fai.elnr[1]);

	      int ndoffacet = felfacet.GetNDof();
	      int ndof1 = fel1.GetNDof();
	      int ndof2 = fel2.GetNDof();

	      FlatVector<> aelu1(ndof1, lh), aelu2(ndof2, lh);
	      FlatVector<> trace1(ndoffacet, lh), trace2(ndoffacet, lh);

	      fel1.GetTrace (fai.facetnr[0], vecu.Range (dn1), trace1);
	      fel2.GetTrace (fai.facetnr[1], vecu.Range (dn2), trace2);

	      traces.FV().Range (i*ndf, (i+1)*ndf) = 0.5 * (trace1 + trace2);
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
	      traces.FV().Range(i*ndf, (i+1)*ndf) = 0.5 * trace1;
	    }
	}
    }


    traces.Cumulate();



#pragma omp parallel 
    {
      LocalHeap & clh = lh, lh = clh.Split();
#pragma omp for
      for (int i = 0; i < ne; i++)
	{
	  HeapReset hr(lh);
	  
	  const ElementData & edi = *elementdata[i];

	  const DGFiniteElement<D> & fel = 
	    dynamic_cast<const DGFiniteElement<D>&> (fes.GetFE (i, lh));
	  
	  
	  IntRange dn = fes.GetElementDofs (i);
	  int ndof = fel.GetNDof();


	  FlatVector<> aelu(ndof, lh);

	  for (int j = 0; j < edi.facets.Size(); j++)
	    {
	      int fnr = edi.facets[j];
	      const FacetData & fai = *facetdata[fnr];	      
	      
	      const DGFiniteElement<D-1> & felfacet = 
		dynamic_cast<const DGFiniteElement<D-1>&> (fes.GetFacetFE (fnr, lh));
	      int ndoffacet = felfacet.GetNDof();
	      
	      FlatVector<> trace(ndoffacet, lh);
	      fel.GetTrace (j, vecu.Range (dn), trace);

	      FlatVector<> otrace(ndoffacet, lh);
	      otrace = 2 * traces.FV().Range (fnr*ndf, (fnr+1)*ndf) - trace;
	      

	      IntegrationRule ir(felfacet.ElementType(), 2*felfacet.Order());
	      int nip = ir.Size();

	      FlatVector<> flown = fai.flown;
	    
	      FlatVector<> tracei1(nip, lh), tracei2(nip, lh);
	      FlatVector<> tracei(nip, lh);

	      felfacet.Evaluate (ir, trace, tracei1);
	      felfacet.Evaluate (ir, otrace, tracei2); 
		    
	      for (int k = 0; k < nip; k++)
		{
		  double fn = flown(k);
		  if (fai.elnr[0] != i) fn *= -1;
		  tracei(k) = fn * ( (fn > 0) ? tracei1(k) : tracei2(k) );
		}

	      felfacet.EvaluateTrans (ir, tracei, trace);
	      fel.GetTraceTrans (j, trace, aelu);

	      conv.Range (dn) -= aelu;
	    }
	}
    }
    
    timer_facet.Stop(); 



  }
};




static RegisterNumProc<NumProcLinearHyperbolic<2> > npinit1("linhyp", 2);
static RegisterNumProc<NumProcLinearHyperbolic<3> > npinit2("linhyp", 3);

  










