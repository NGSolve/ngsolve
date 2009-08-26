#include <solve.hpp>


namespace ngsolve
{
  using namespace ngsolve;




  
  

  ///
  class NumProcApplyMat : public NumProc
  {
    const BilinearForm * bfa;
  public: 
    ///
    NumProcApplyMat (PDE & apde, const Flags & flags);
    ///
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcApplyMat (pde, flags);
    }
    ///
    virtual void Do(LocalHeap & lh);
  };
  

  
  NumProcApplyMat ::   
  NumProcApplyMat (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    bfa = apde.GetBilinearForm ("ar");
  }

  void NumProcApplyMat :: Do(LocalHeap & lh)
  {
    const BaseMatrix & mat = bfa->GetMatrix();
    cout << "h = " << mat.Height();
    VVector<Vec<1> > x(mat.Height()), y(mat.Height());
    x = 1;
    y = 0;
    y = mat * x;
    (*testout) << "y = " << y << endl;
  }
  




#ifdef JOACHIM
///
class NumProcTestJS : public NumProc
{
public:
  ///
  NumProcTestJS (PDE & apde, const Flags & flags);
  ///
  virtual ~NumProcTestJS ();
  ///
  virtual void Do();
  
  ///
  void TestHP();
  ///
  void CalcBoundaryCurrent();
  ///
  void CalcDisplacement();
  ///
  void TestBEM();
  ///
  void TestIntegration();
  void TestIE();
  void LaplaceTiming();
};







NumProcTestJS :: NumProcTestJS (PDE & apde, const Flags & flags)
  : NumProc (apde)
{
  ;
}

NumProcTestJS :: ~NumProcTestJS ()
  {
    ;
  }


void NumProcTestJS :: Do()
  {
    LaplaceTiming();
    return;
  }

  void NumProcTestJS :: TestHP()
  {
    cout << "test hp" << endl;
  }





  void NumProcTestJS :: CalcBoundaryCurrent()
  {
    ;
  }




  void NumProcTestJS :: CalcDisplacement()
  {
    ;
  }


  void NumProcTestJS :: TestBEM()
  {
    ;
  }
  void NumProcTestJS :: TestIntegration()
  {
    ;
  }


  void NumProcTestJS::TestIE()
  {
    ;
  }



  double ScalN (int n, const double * pa, const double * pb)
  {
    double sum = 0;
    int i;
    for (i = 0; i < n; i++)
      sum += pa[i] * pb[i];
    return sum;
  }

  template <int N>
  double Scal (const double * pa, const double * pb)
  {
    return *pa * *pb + Scal<N-1> (pa+1, pb+1);
  }
  template <> double Scal<1> (const double * pa, const double * pb)
  {
    return *pa * *pb;
  }



#define RUNTIMEaaa


#ifdef old
  void AssembleLaplace (ElementTransformation & eltrans,
			DenseMatrix & elmat,
			LocalHeap * locheap)
  {
    int i, j, k;
    typedef FE_Trig1 FE;
    enum { nd = FE::NDOF };
    enum { dim = FE::SDIM };
    //    int nd = elmat.Height();
    //    int dim = eltrans.GetElement().SpatialDim();

    const IntegrationRule & ir = GetIntegrationRules().SelectIntegrationRule (ET_TRIG, 1);
    const IntegrationPoint & ip = ir.GetIP(1);
    const ngbla::Mat<2,3> & dshape = FE::ipdata.Get (ip.IPNr()) . dshape;

    const ngbla::FlatMatrix<> & pmat = eltrans.PointMatrix();
    const ngbla::Mat<2,3> & dtrans = FE::ipdata.Get (ip.IPNr()) . dshape;

    //    double fdsh[dim*nd], fjinv[dim*dim], fbmat[dim*nd];
    //    double fdxdxi[dim*dim];
    double * fdxdxi = static_cast<double*> (locheap->Alloc (dim*dim*sizeof(double)));
    double * fdsh = static_cast<double*> (locheap->Alloc (dim*nd*sizeof(double)));
    double * fjinv = static_cast<double*> (locheap->Alloc (dim*dim*sizeof(double)));
    double * fbmat = static_cast<double*> (locheap->Alloc (dim*nd*sizeof(double)));


    for (i = 0; i < dim; i++)
      for (j = 0; j < dim; j++)
	{
	  fdxdxi[dim*i+j] =
#ifdef RUNTIME
	    Scal<nd> (&pmat.Get(i*nd+1), &dtrans.Get(j*nd+1));
#else
	  ScalN (nd, &pmat.Get(i*nd+1), &dtrans.Get(j*nd+1));
#endif
	}

    /*
#ifdef RUNTIME
    fdxdxi[0] = Scal<nd> (&pmat.Get(1), &dtrans.Get(1));
    fdxdxi[1] = Scal<nd> (&pmat.Get(1), &dtrans.Get(nd+1));
    fdxdxi[2] = Scal<nd> (&pmat.Get(nd+1), &dtrans.Get(1));
    fdxdxi[3] = Scal<nd> (&pmat.Get(nd+1), &dtrans.Get(nd+1));
#else
    fdxdxi[0] = ScalN (nd, &pmat.Get(1), &dtrans.Get(1));
    fdxdxi[1] = ScalN (nd, &pmat.Get(1), &dtrans.Get(nd+1));
    fdxdxi[2] = ScalN (nd, &pmat.Get(nd+1), &dtrans.Get(1));
    fdxdxi[3] = ScalN (nd, &pmat.Get(nd+1), &dtrans.Get(nd+1));
#endif
    */

    double det = fdxdxi[0] * fdxdxi[3] - fdxdxi[1] * fdxdxi[2];
    double idet = 1/det;
    fjinv[0] = fdxdxi[3] * idet;
    fjinv[1] = -fdxdxi[1] * idet;
    fjinv[2] = -fdxdxi[2] * idet;
    fjinv[3] = fdxdxi[0] * idet;

    for (i = 0; i < dim; i++)
      for (j = 0; j < nd; j++)
	{
	  fbmat[i+dim*j] = 
#ifdef RUNTIME
	    Scal<dim> (&fjinv[dim*i], &dshape.Get(dim*j+1));  // not ok (transpose)
#else
	    ScalN (dim, &fjinv[dim*i], &dshape.Get(dim*j+1));  // not ok (transpose)
#endif
	}

    for (i = 0; i < nd; i++)
      for (j = 0; j < nd; j++)
	elmat.Elem(i*nd+j+1) =
#ifdef RUNTIME
	  Scal<dim> (fbmat+dim*i, fbmat+dim*j);
#else
	  ScalN (dim, fbmat+dim*i, fbmat+dim*j);
#endif
  }

  using namespace ngbla;

  void AssembleLaplace2 (ElementTransformation & eltrans,
			 DenseMatrix & elmat,
			 LocalHeap * locheap)
  {
    int i, j, k;
    typedef FE_Trig1 FE;
    //    enum { nd = FE::NDOF };
    enum { DIM = FE::SDIM };
    int nd = 3; // elmat.Height();
    int dim = DIM; // eltrans.GetElement().SpatialDim();

    const IntegrationRule & ir = GetIntegrationRules().SelectIntegrationRule (ET_TRIG, 1);
    const IntegrationPoint & ip = ir.GetIP(1);
    const ngbla::Mat<2,3> & dshape = FE::ipdata.Get (ip.IPNr()) . dshape;

    const ngbla::FlatMatrix<> & pmat = eltrans.PointMatrix();
    const ngbla::Mat<2,3> & dtrans = FE::ipdata.Get (ip.IPNr()) . dshape;


    FlatMatrix<double> fpmat (dim, nd, (double*)&pmat.Get(1,1));
    FlatMatrix<double> fdtrans (dim, nd, (double*)&dtrans.Get(1,1));
    FlatMatrix<double> fdsh (dim, nd, (double*)&dshape.Get(1,1));

    Mat<DIM> fdxdxi = fpmat * Trans (fdtrans); 
    Mat<DIM> fjinv = Inv (fdxdxi);

    FlatMatrix<> fbmat (dim, nd, locheap->Alloc<double>(dim*nd));
    FlatMatrix<> felmat (nd, nd, locheap->Alloc<double>(nd*nd));

    double det = 1; // fdxdxi[0] * fdxdxi[3] - fdxdxi[1] * fdxdxi[2];

    fbmat = Trans (fjinv) * fdsh;
    felmat = Trans (fbmat) * fbmat;
  }
#endif









  void NumProcTestJS::LaplaceTiming()
  {
    cout << "laplace timing" << endl;
    const MeshAccess & ma = pde.GetMeshAccess();

    clock_t starttime, endtime;
    starttime = clock();

    LaplaceIntegrator<2,FE_Trig1> lpi (new ConstantCoefficientFunction(1));
    ElasticityIntegrator<2> eli (new ConstantCoefficientFunction(1),
				 new ConstantCoefficientFunction(0.2));

    SourceIntegrator<2> source (new ConstantCoefficientFunction(1));


    //    BlockBilinearFormIntegrator<LaplaceIntegrator<2> > blpi(lpi, 3);
    //    BlockBilinearFormIntegrator<LaplaceIntegrator<2,FE_Trig1> > blpi(lpi, 3);

    BilinearFormIntegrator * bfi = &lpi;
    LinearFormIntegrator * lfi = &source;

    /*
    ElasticityIntegrator eli (2, new ConstantCoefficientFunction(1),
			      new ConstantCoefficientFunction(0.2));
    */
    FE_Trig1 trig;
    ElementTransformation eltrans;
    //    DenseMatrix elmat(trig.GetNDof());
    LocalHeap locheap (100000);

    int n = 1000000;

    ngbla::FlatMatrix<double> elmat;
    ngbla::FlatVector<double> elvec(3, locheap);


    void * heapp = locheap.GetPointer();

    for (int i = 1; i <= n; i++)
      {
	locheap.CleanUp(heapp);
	ma.GetElementTransformation (1, eltrans, locheap);
	bfi->AssembleElementMatrix(trig, eltrans, elmat, locheap);
	//	elvec.Assign (lfi->AssembleElementVector (trig, eltrans, locheap));
	// AssembleLaplace2 (eltrans, elmat, &locheap);
      }

    cout << "elmat = " << elmat << endl;
    cout << "elvec = " << elvec << endl;

    endtime = clock();
    cout << "sec per 1E6 els = " << 1e6 / n * double(endtime - starttime)/CLOCKS_PER_SEC << endl;

    /*
      Laplace, trig1: 169000
      Laplace, tet1:  141000

      Elast, trig1:   103000
      Elast, trig2:    18000
      Elast, tet1:     31000
      Elast, tet2:     2061
    */
   
  }

#endif

  namespace {
    class Init
    { 
    public: 
      Init ();
    };
    
    Init::Init()
    {
      GetNumProcs().AddNumProc ("applymat", NumProcApplyMat::Create);
    }

    Init init;
  }

}
