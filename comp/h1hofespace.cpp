/*********************************************************************/
/* File:   h1hofespace.cpp                                           */
/* Author: Start                                                     */
/* Date:   10. Feb. 2003                                             */
/*********************************************************************/

/**
   High Order Finite Element Space
*/

// #include <comp.hpp>
#include <h1hofespace.hpp>
#include <multigrid.hpp> 
#include "../fem/h1hofe.hpp"
#include "../fem/h1hofefo.hpp"
#include <../fem/hdivhofe.hpp>
#include <../fem/facethofe.hpp>
#include <../fem/nodalhofe.hpp>
#include <l2hofe.hpp>
#include <bdbequations.hpp>
#include <diffop_impl.hpp>
#include <fesconvert.hpp>

using namespace ngmg; 


#ifdef PARALLEL

#include "../parallel/dump.hpp"



template <NODE_TYPE NT, typename TELEM>
class NodalArray
{
  const MeshAccess & ma;
  Array<TELEM> & a;
public:
  NodalArray (const MeshAccess & ama, Array<TELEM> & aa) : ma(ama), a(aa) { ; }
  const MeshAccess & GetMeshAccess() const { return ma; }
  Array<TELEM> & A() { return a; }
};

template <NODE_TYPE NT, typename TELEM>
auto NodalData (const MeshAccess & ama, Array<TELEM> & a) -> NodalArray<NT,TELEM> 
{ return NodalArray<NT,TELEM> (ama, a); }



template <NODE_TYPE NT, typename T> 
Archive & operator & (Archive & archive, NodalArray<NT,T> && a)
{
  auto comm = a.GetMeshAccess().GetCommunicator();
  if (comm.Size() == 1) return archive & a.A();
  
  auto g = [&] (int size) { archive & size; };    

  typedef typename key_trait<NT>::TKEY TKEY;
  auto f = [&] (TKEY key, T val) { archive & val; };
      
  GatherNodalData<NT> (a.GetMeshAccess(), a.A(), g, f);

  return archive;
}



#else
template <NODE_TYPE NT, typename TELEM>
auto NodalData (MeshAccess & ma, Array<TELEM> & a) -> Array<TELEM> & { return a; }
#endif


namespace ngcomp
{



  /// dual operator for H1
  template <int D>
  class DiffOpDual : public DiffOp<DiffOpDual<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 0 };

    static bool SupportsVB (VorB checkvb) { return true; }


    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,BareSliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT && mat, LocalHeap & lh)
    {
      mat.AddSize(DIM_DMAT, fel.GetNDof()) = 0.0;
      static_cast<const ScalarFiniteElement<D>&>(fel).CalcDualShape (mip, mat.Row(0));
    }
    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<!std::is_convertible<MAT,BareSliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT && mat, LocalHeap & lh)
    {
      // fel.CalcDualShape (mip, mat);
      throw Exception(string("DiffOpDual not available for mat ")+typeid(mat).name());
    }

    using DiffOp<DiffOpDual<D>>::AddTransSIMDIR;
    template <typename FEL, class MIR> // , class TVY>
    static void AddTransSIMDIR (const FEL & fel, const MIR & mir,
                                BareSliceMatrix<SIMD<double>> x, BareSliceVector<double> y)
    {
      STACK_ARRAY(SIMD<double>, memx, mir.Size());
      FlatVector<SIMD<double>> hx(mir.Size(), &memx[0]);
      for (size_t i = 0; i < mir.Size(); i++)
        hx(i) = x(0,i) / mir[i].GetMeasure();
      static_cast<const ScalarFiniteElement<D>&>(fel).AddDualTrans (mir.IR(), hx, y);
    }    
  };


  
  template <int _DIM_SPACE, int _DIM_ELEMENT>
  class DiffOpDualH1 : public DiffOpDual<_DIM_SPACE>
  {
  public:
    enum { DIM_SPACE = _DIM_SPACE };
    enum { DIM_ELEMENT = _DIM_ELEMENT };

    typedef DiffOpDualH1<_DIM_SPACE, _DIM_ELEMENT-1> DIFFOP_TRACE;
  };
  
  template <int _DIM_SPACE>
  class DiffOpDualH1<_DIM_SPACE,0> : public DiffOpId<_DIM_SPACE>
  {
  public:    
    enum { DIM_SPACE = _DIM_SPACE };
    enum { DIM_ELEMENT = 0 };

    typedef void DIFFOP_TRACE;
  };


  /// dual operator for VectorH1
  template <int D>
  class DiffOpDualVector : public DiffOp<DiffOpDualVector<D> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = D };
    enum { DIM_ELEMENT = D };
    enum { DIM_DMAT = D };
    enum { DIFFORDER = 0 };

    static bool SupportsVB (VorB checkvb) { return true; }

    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<std::is_convertible<MAT,BareSliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & bfel, const MIP & mip,
                                MAT && mat, LocalHeap & lh)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      mat.AddSize(DIM_DMAT, bfel.GetNDof()) = 0.0;
      for (int i = 0; i < DIM_SPACE; i++)
        {
          auto & feli = static_cast<const ScalarFiniteElement<DIM_ELEMENT>&> (fel[i]);
          feli.CalcDualShape (mip, mat.Row(i).Range(fel.GetRange(i)));
        }
    }
    
    template <typename AFEL, typename MIP, typename MAT,
              typename std::enable_if<!std::is_convertible<MAT,BareSliceMatrix<double,ColMajor>>::value, int>::type = 0>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
                                MAT & mat, LocalHeap & lh)
    {
      throw Exception("DiffOpDual not available for mat ", typeid(mat).name());
    }
  };
  
  template <int _DIM_SPACE, int _DIM_ELEMENT>
  class DiffOpDualVectorH1 : public DiffOpDualVector<_DIM_SPACE>
  {
  public:
    enum { DIM_SPACE = _DIM_SPACE };
    enum { DIM_ELEMENT = _DIM_ELEMENT };

    typedef DiffOpDualVectorH1<_DIM_SPACE, _DIM_ELEMENT-1> DIFFOP_TRACE;
  };
  
  template <int _DIM_SPACE>
  class DiffOpDualVectorH1<_DIM_SPACE,0> : public DiffOpId<_DIM_SPACE>
  {
  public:    
    enum { DIM_SPACE = _DIM_SPACE };
    enum { DIM_ELEMENT = 0 };

    typedef void DIFFOP_TRACE;
  };

  class H1HOProlongation : public Prolongation
  {
    weak_ptr<H1HighOrderFESpace> fes;
    shared_ptr<FESpace> fesL2;

    mutable Array<shared_ptr<BaseMatrix>> convL2toH1;
    mutable Array<shared_ptr<BaseMatrix>> convH1toL2;
  public:
    H1HOProlongation (weak_ptr<H1HighOrderFESpace> afes)
      : fes(afes)
    {
      Flags flagsL2;
      flagsL2.SetFlag ("order", fes.lock()->GetOrder());
      fesL2 = CreateFESpace ("L2", fes.lock()->GetMeshAccess(), flagsL2);
      // fesL2->Update();
      // fesL2->FinalizeUpdate();
      // int levels = fes.lock()->GetMeshAccess()->GetNLevels();
      // convL2toH1.SetSize(levels);
      // convH1toL2.SetSize(levels);

      // LocalHeap lh(10*1000*1000);
      // convL2toH1[levels-1] = ConvertOperator(fesL2, fes.lock(), VOL, lh);
      // convH1toL2[levels-1] = ConvertOperator(fes.lock(), fesL2, VOL, lh);
    }

    virtual void Update (const FESpace & /* fes*/) override
    {
      fesL2->Update();
      fesL2->FinalizeUpdate();
      
      int levels = fes.lock()->GetMeshAccess()->GetNLevels();
      if (convL2toH1.Size() < levels)
        {
          convL2toH1.SetSize(levels);
          convH1toL2.SetSize(levels);
          
          LocalHeap lh(10*1000*1000);
          convL2toH1[levels-1] = ConvertOperator(fesL2, fes.lock(), VOL, lh, nullptr, nullptr, NULL, nullptr, false, true, true, 0, 0, true);
          convH1toL2[levels-1] = ConvertOperator(fes.lock(), fesL2, VOL, lh, nullptr, nullptr, NULL, nullptr, false, true, true, 0, 0, true);
        }

    }

    virtual size_t GetNDofLevel (int level) override
    {
      return fes.lock()->GetNDofLevel(level);
    }
    

    ///
    virtual shared_ptr<SparseMatrix< double >> CreateProlongationMatrix( int finelevel ) const override
    { return NULL; }

    ///
    virtual void ProlongateInline (int finelevel, BaseVector & v) const override
    {
        /*
      if (convL2toH1.Size() < finelevel+1)
        {
          convL2toH1.SetSize(finelevel+1);
          convH1toL2.SetSize(finelevel+1);
          
          LocalHeap lh(10*1000*1000);
          convL2toH1[finelevel] = ConvertOperator(fesL2, fes.lock(), VOL, lh, nullptr, nullptr, NULL, nullptr, false, true, true, 0, 0, true);
          convH1toL2[finelevel] = ConvertOperator(fes.lock(), fesL2, VOL, lh, nullptr, nullptr, NULL, nullptr, false, true, true, 0, 0, true);
        }
        */

      auto vl2 = convL2toH1[finelevel]->CreateRowVector();

      auto shapec = convH1toL2[finelevel-1]->Shape();
      auto shapef = convL2toH1[finelevel]->Shape();

      vl2.Range(get<0>(shapec)) = *convH1toL2[finelevel-1] * v.Range(get<1>(shapec));      
      fesL2->GetProlongation()->ProlongateInline(finelevel, vl2);
      v.Range(get<0>(shapef)) = *convL2toH1[finelevel] * vl2.Range(get<1>(shapef));
    }    

    
    ///
    virtual void RestrictInline (int finelevel, BaseVector & v) const override
    {
      /*
      if (convL2toH1.Size() < finelevel+1)
        {
          convL2toH1.SetSize(finelevel+1);
          convH1toL2.SetSize(finelevel+1);
          
          LocalHeap lh(10*1000*1000);
          convL2toH1[finelevel] = ConvertOperator(fesL2, fes.lock(), VOL, lh, nullptr, nullptr, NULL, nullptr, false, true, true, 0, 0, true);
          convH1toL2[finelevel] = ConvertOperator(fes.lock(), fesL2, VOL, lh, nullptr, nullptr, NULL, nullptr, false, true, true, 0, 0, true);
        }
        */
      auto vl2 = convL2toH1[finelevel]->CreateRowVector();

      auto shapec = convH1toL2[finelevel-1]->Shape();
      auto shapef = convL2toH1[finelevel]->Shape();

      vl2.Range(get<1>(shapef)) = Transpose(*convL2toH1[finelevel]) * v.Range(get<0>(shapef));      
      fesL2->GetProlongation()->RestrictInline(finelevel, vl2);
      v.Range(get<1>(shapec)) = Transpose(*convH1toL2[finelevel-1]) * vl2.Range(get<0>(shapec));
    }    
  };



  
  H1HighOrderFESpace ::  
  H1HighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags)
    : FESpace (ama, flags)
  {
    name = "H1HighOrderFESpace(h1ho)";
    type = "h1ho";
    // define h1ho flags
    DefineDefineFlag("h1ho");
    DefineNumFlag("relorder");
    DefineNumFlag("orderinner");
    DefineNumFlag("orderface");
    DefineNumFlag("orderedge");
    DefineNumFlag("orderquad");
    DefineNumFlag("ordertrig");
    DefineNumFlag("variableorder"); 
    //  DefineNumListFlag("dom_order_min_x");
    //  DefineNumListFlag("dom_order_max_x");
    //  DefineNumListFlag("dom_order_min_y");
    //  DefineNumListFlag("dom_order_max_y");
    //  DefineNumListFlag("dom_order_min_z");
    //  DefineNumListFlag("dom_order_max_z");
    DefineNumFlag("smoothing");
    DefineDefineFlag("wb_withedges");
    if (parseflags) CheckFlags(flags);

    wb_loedge = ma->GetDimension() == 3;
    if (flags.GetDefineFlag("wb_withedges")) wb_loedge = true;
    if (flags.GetDefineFlag("wb_withoutedges") ||
        flags.GetDefineFlagX("wb_withedges").IsFalse()) wb_loedge = false;
    wb_edge = flags.GetDefineFlag ("wb_fulledges");
    
    // Variable order space: 
    //      in case of (var_order && order) or (relorder) 
    var_order = flags.GetDefineFlag("variableorder");  
    fixed_order = flags.GetDefineFlag("fixedorder");  
    order = int (flags.GetNumFlag ("order",1)); 
    if (order < 1) order = 1;

    if(flags.NumFlagDefined("relorder") && !flags.NumFlagDefined("order")) 
      var_order = 1; 
    
    rel_order=int(flags.GetNumFlag("relorder",order-1)); 

    print = flags.GetDefineFlag("print");

    if(flags.NumFlagDefined("order") && flags.NumFlagDefined("relorder")) 
      {
	if(var_order)
	  cerr << " WARNING: H1HoFeSpace: inconsistent flags: variableorder, order and relorder "
	       << "-> variable order space with rel_order " << rel_order << "is used, but order is ignored " << endl; 
	else 
	  cerr << " WARNING: H1HoFeSpace: inconsistent flags: order and rel_order "
	       << "-> uniform order space with order " << order << " is used " << endl; 
      }
    
    uniform_order_inner = int (flags.GetNumFlag ("orderinner", -1));
    uniform_order_face = int (flags.GetNumFlag ("orderface", -1));
    uniform_order_edge = int (flags.GetNumFlag ("orderedge", -1));
    uniform_order_quad = int (flags.GetNumFlag ("orderquad", -1));
    uniform_order_trig = int (flags.GetNumFlag ("ordertrig", -1));
    
    if (flags.NumFlagDefined("smoothing")) 
      throw Exception ("Flag 'smoothing' for fespace is obsolete \n Please use flag 'blocktype' in preconditioner instead");
    nodalp2 = flags.GetDefineFlag ("nodalp2");
    nodal = flags.GetDefineFlag ("nodal");    
    
    highest_order_dc = flags.GetDefineFlag ("highest_order_dc");
    if (highest_order_dc && order < 2)
      throw Exception ("highest_order_dc needs order >= 2");


    test_ho_prolongation = flags.GetDefineFlag("hoprolongation");
    if (test_ho_prolongation)
      no_low_order_space=true;
    
    Flags loflags = flags;
    loflags.SetFlag ("order", 1);

    if (!no_low_order_space)
      low_order_space = make_shared<NodalFESpace> (ma, loflags);
    switch (ma->GetDimension())
      {
      case 1:
        {
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<1>>>();
          flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<1>>>();
          evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<1>>>();
          break;
        }
      case 2:
        {
          // evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdH1<2,2>>>();
          flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
          // evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();
          evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdH1<2,1>>>();
          flux_evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpGradientBoundary<2>>>();
          evaluator[BBND] = make_shared<T_DifferentialOperator<DiffOpIdH1<2,0>>>();
          break;
        }
      case 3:
        {
          // evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdH1<3,3>>>();
          flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<3>>>();
          // evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<3>>>();
          evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdH1<3,2>>>();
          flux_evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpGradientBoundary<3>>>();
          // evaluator[BBND] = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
          evaluator[BBND] = make_shared<T_DifferentialOperator<DiffOpIdH1<3,1>>>();
	  flux_evaluator[BBND] = make_shared<T_DifferentialOperator<DiffOpGradientBBoundary<3>>>();
          // evaluator[BBBND] = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
          evaluator[BBBND] = make_shared<T_DifferentialOperator<DiffOpIdH1<3,0>>>();
          break;
        }
      }
    if (dimension > 1)
      {
        additional_evaluators.Set ("Grad", make_shared<BlockDifferentialOperatorTrans>(flux_evaluator[VOL], dimension));
        if (ma->GetDimension() >= 2)        
          additional_evaluators.Set ("Gradboundary", make_shared<BlockDifferentialOperatorTrans>(flux_evaluator[BND], dimension));
        for (auto vb : { VOL,BND, BBND, BBBND })
          {
            if (evaluator[vb])
              evaluator[vb] = make_shared<BlockDifferentialOperator> (evaluator[vb], dimension);
            if (flux_evaluator[vb])
              flux_evaluator[vb] = make_shared<BlockDifferentialOperator> (flux_evaluator[vb], dimension);
          }
      }
    else
      {
        switch (ma->GetDimension())
          {
          case 1:
            additional_evaluators.Set ("Grad", make_shared<T_DifferentialOperator<DiffOpGradient<1>>>());
            break;
          case 2:
            additional_evaluators.Set ("Grad", make_shared<T_DifferentialOperator<DiffOpGradient<2>>>());
            additional_evaluators.Set ("Gradboundary", make_shared<T_DifferentialOperator<DiffOpGradientBoundary<2>>>());
            break;
          case 3:
            additional_evaluators.Set ("Grad", make_shared<T_DifferentialOperator<DiffOpGradient<3>>>());
            additional_evaluators.Set ("Gradboundary", make_shared<T_DifferentialOperator<DiffOpGradientBoundary<3>>>());
            break;
          default:
            ;
          }
      }

    switch (ma->GetDimension())
      {
      case 1:
        additional_evaluators.Set ("hesse", make_shared<T_DifferentialOperator<DiffOpHesse<1>>> ());
        break;
      case 2:
        additional_evaluators.Set ("hesse", make_shared<T_DifferentialOperator<DiffOpHesse<2>>> ());
        additional_evaluators.Set ("hesseboundary", make_shared<T_DifferentialOperator<DiffOpHesseBoundary<2>>> ());
	if (dimension > 1)
	  { additional_evaluators.Set ("dual", make_shared<BlockDifferentialOperator> (make_shared<T_DifferentialOperator<DiffOpDualH1<2,2>>>(), dimension)); }
	else
	  { additional_evaluators.Set ("dual", make_shared<T_DifferentialOperator<DiffOpDualH1<2,2>>> ()); }
        break;
      case 3:
        additional_evaluators.Set ("hesse", make_shared<T_DifferentialOperator<DiffOpHesse<3>>> ());
	additional_evaluators.Set ("hesseboundary", make_shared<T_DifferentialOperator<DiffOpHesseBoundary<3>>> ());
	if (dimension > 1)
	  { additional_evaluators.Set ("dual", make_shared<BlockDifferentialOperator> (make_shared<T_DifferentialOperator<DiffOpDualH1<3,3>>> (), dimension)); }
	else
	  { additional_evaluators.Set ("dual", make_shared<T_DifferentialOperator<DiffOpDualH1<3,3>>> ()); }
	break;
      default:
        ;
      }

    if (!test_ho_prolongation)
      prol = make_shared<LinearProlongation> (GetMeshAccess());

    needs_transform_vec = false;
  }


  H1HighOrderFESpace :: ~H1HighOrderFESpace ()
  {
    ;
  }

  DocInfo H1HighOrderFESpace :: GetDocu ()
  {
    DocInfo docu = FESpace::GetDocu();
    docu.short_docu = "An H1-conforming finite element space.";
    docu.long_docu =
      R"raw_string(The H1 finite element space consists of continuous and
element-wise polynomial functions. It uses a hierarchical (=modal)
basis built from integrated Legendre polynomials on tensor-product elements,
and Jaboci polynomials on simplicial elements. 

Boundary values are well defined. The function can be used directly on the
boundary, using the trace operator is optional.

The H1 space supports variable order, which can be set individually for edges, 
faces and cells. 

Internal degrees of freedom are declared as local dofs and are eliminated 
if static condensation is on.

The wirebasket consists of all vertex dofs. Optionally, one can include the 
first (the quadratic bubble) edge basis function, or all edge basis functions
into the wirebasket.
)raw_string";      

    
    docu.Arg("wb_withedges") = "bool = true(3D) / false(2D)\n"
      "  use lowest-order edge dofs for BDDC wirebasket";
    docu.Arg("wb_fulledges") = "bool = false\n"
      "  use all edge dofs for BDDC wirebasket";
    docu.Arg("hoprolongation") = "bool = false\n"
      "  (experimental, only trigs) creates high order prolongation,\n"
      "  and switches off low-order space";
    return docu;
  }

  
  void H1HighOrderFESpace :: Update()
  {
    static Timer timer ("H1HighOrderFESpace::Update");
    // static Timer timer1 ("H1HighOrderFESpace::Update 1");
    // static Timer timer2 ("H1HighOrderFESpace::Update 2");
    // static Timer timer3 ("H1HighOrderFESpace::Update 3");
    RegionTimer reg(timer);

    // timer1.Start();
    FESpace::Update();

    if (order_policy == CONSTANT_ORDER)
      fixed_order = true;
    
    TORDER maxorder = 0;
    TORDER minorder = 99; 

    if (low_order_space) low_order_space -> Update();
    
    bool first_update = GetTimeStamp() < ma->GetTimeStamp();
    if (first_update) timestamp = NGS_Object::GetNextTimeStamp();
    
    int dim = ma->GetDimension();
    size_t nv = ma->GetNV();
    size_t ned = (dim <= 1) ? 0 : ma->GetNEdges();
    size_t nfa = (dim <= 2) ? 0 : ma->GetNFaces();
    size_t ne = ma->GetNE();

    if (first_update)
      {
	used_edge.SetSize(ned); 
	used_face.SetSize(nfa); 
	used_vertex.SetSize(nv); 

	used_edge = false; 
	used_face = false; 
	used_vertex = false; 

	// for (FESpace::Element el : Elements (VOL))

	// for (auto vb : { VOL, BND, BBND })
        for (auto vb : Range(VOL, VorB(dim+1)))
          ParallelFor
            (ma->GetNE(vb), [&] (size_t nr)
             {
               ElementId ei(vb, nr);
               Ngs_Element el = (*ma)[ei];
           
	       // if (!DefinedOn (el)) return;
               if (DefinedOnX (el).IsTrue())
                 {
                   used_vertex[el.Vertices()] = true;
                   if (dim >= 2) used_edge[el.Edges()] = true;
                   if (dim == 3) used_face[el.Faces()] = true;
                 }
	     });
    
	/*
	  for (FESpace::Element el : Elements (BND))
	  {
	  used_vertex[el.Vertices()] = true;
	  if (dim >= 2) used_edge[el.Edges()] = true;
	  if (dim == 3) used_face[el.Faces()] = true;
	  }
	*/
    
	ma->AllReduceNodalData (NT_VERTEX, used_vertex, NG_MPI_LOR);
	ma->AllReduceNodalData (NT_EDGE, used_edge, NG_MPI_LOR);
	ma->AllReduceNodalData (NT_FACE, used_face, NG_MPI_LOR);

	// timer1.Stop();
	// timer2.Start();
    
	order_edge.SetSize (ned);
	order_face.SetSize (nfa);
	order_inner.SetSize (ne);

	int p = var_order ?  1 : order; 
    
	order_edge = p; 
	order_face = p; 
	order_inner = p;
	
	if(var_order) 
	  for (Ngs_Element el : ma->Elements<VOL>())
	    {	
	      if (!DefinedOn (el)) continue;
	      int i = el.Nr();
          
	      ELEMENT_TYPE eltype = el.GetType(); 
	      const FACE * faces = ElementTopology::GetFaces (eltype);
	      const EDGE * edges = ElementTopology::GetEdges (eltype);
	      const POINT3D * points = ElementTopology :: GetVertices (eltype);

	      auto vnums = el.Vertices();
	      auto eledges = el.Edges();
	
	      IVec<3,TORDER> el_orders = ma->GetElOrders(i) + IVec<3> (rel_order);

	      maxorder = max2 (maxorder, Max(el_orders));
	      minorder = min2 (minorder, Min(el_orders));
	      // for(int l=0;l<3;l++) maxorder = max2(el_orders[l],maxorder); 
	      // for(int l=0;l<3;l++) minorder = min2(el_orders[l],minorder); 
          
	      order_inner[i] = Max (order_inner[i], el_orders + IVec<3,TORDER>(et_bonus_order[eltype]));
	      // for(int j=0;j<dim;j++) order_inner[i][j] = max2(order_inner[i][j],el_orders[j]);

	      for(int j=0;j<eledges.Size();j++)
		{
		  for(int k=0;k<dim;k++)
		    if(points[edges[j][0]][k] != points[edges[j][1]][k])
		      { 
			order_edge[eledges[j]] = max2(order_edge[eledges[j]],
						      TORDER(el_orders[k]+et_bonus_order[ET_SEGM]));
			k=dim; 
		      }
		}
          
	      if(dim==3)
		{
		  auto elfaces = el.Faces();              
		  for(int j=0;j<elfaces.Size();j++)
		    {
		      // trig_face
		      if(faces[j][3]==-1) 
			{
			  order_face[elfaces[j]][0] = 
			    int(max2(order_face[elfaces[j]][0], 
				     TORDER(el_orders[0]+et_bonus_order[ET_TRIG])));
			  order_face[elfaces[j]][1] = order_face[elfaces[j]][0]; 
			}
		      else //quad_face
			{
			  int fmax = 0;
			  for(int k = 1; k < 4; k++) 
			    if(vnums[faces[j][k]] > vnums[faces[j][fmax]]) fmax = k;   
                      
			  IVec<2> f((fmax+3)%4,(fmax+1)%4); 
			  if(vnums[faces[j][f[1]]] > vnums[faces[j][f[0]]]) swap(f[0],f[1]);
                      
			  for(int l=0;l<2;l++)
			    for(int k=0;k<3;k++)
			      if(points[faces[j][fmax]][k] != points[faces[j][f[l] ]][k])
				{
				  order_face[elfaces[j]][l] = 
				    int(max2(order_face[elfaces[j]][l], 
					     TORDER(el_orders[k]+et_bonus_order[ET_TRIG])));
				  break;
				}
			}
		    }
		}
	    }
    
	else  // not var_order

	  {
	    // for (Ngs_Element el : ma->Elements<VOL>())
	    ParallelFor (ma->GetNE(VOL), [&] (size_t nr)
			 {
			   ElementId ei(VOL, nr);
			   Ngs_Element el = (*ma)[ei];
                     
			   if (!DefinedOn (el)) return; 
                     
			   if (dim >= 2)
			     for (auto e : el.Edges())
			       order_edge[e] = p + et_bonus_order[ET_SEGM];
                     
			   if (dim == 3)
			     for (auto f : el.Faces())
			       order_face[f] = p + et_bonus_order[ma->GetFaceType(f)];

			   order_inner[el.Nr()] = p + et_bonus_order[el.GetType()];
			 });
	  }
	/* 
	   if (ma->GetDimension() == 2 && uniform_order_trig != -1 && uniform_order_quad != -1)
	   {
	   for (int i = 0; i < nel; i++)
	   {
	   if (ma->GetElType(i) == ET_TRIG)
	   order_inner = IVec<3> (uniform_order_trig, uniform_order_trig, uniform_order_trig);
	   else
	   order_inner = IVec<3> (uniform_order_quad, uniform_order_quad, uniform_order_quad);
	   }
	   }
	*/

	// timer2.Stop();
	// timer3.Start();
    
	if(uniform_order_inner > -1)  
	  order_inner = uniform_order_inner;
	if(uniform_order_face > -1 && dim == 3) 
	  order_face = uniform_order_face;
	if(uniform_order_edge > -1)   
	  order_edge = uniform_order_edge; 

	for (auto i : Range(used_edge))
	  if (!used_edge[i]) order_edge[i] = 1; 

	for (auto i : used_face.Range())
	  if (!used_face[i]) order_face[i] = 1; 

	for (ElementId ei : ma->Elements<VOL>())
	  if (!DefinedOn(ei)) order_inner[ei.Nr()] = 1;

	if(print) 
	  {
	    *testout << " H1HoFESpace order " << order << " , var_order " << var_order << " , relorder " << rel_order << endl;  
	    (*testout) << "used_vertex (h1): " << used_vertex << endl;
	    (*testout) << "used_edge (h1): " << used_edge << endl;
	    (*testout) << "used_face (h1): " << used_face << endl;
	
	
	    (*testout) << "order_edge (h1): " << order_edge << endl;
	
	    (*testout) << "order_face (h1): " << order_face <<  endl;
	    (*testout) << "order_inner (h1): " << order_inner << endl;
	  }
    
#ifdef SABINE 
	// fuer EE damit ich mit rel_order auch p-ref machen kann :) 
	order=order_inner[0][0]; 

	if (!var_order) { maxorder = order; minorder = order; }; 
	order = maxorder;

	cout << " H1FESPACE : " << minorder << " <= order <= " << maxorder << endl;  
#endif 
      }
    UpdateDofTables ();
    UpdateCouplingDofArray ();

    if (low_order_space)
      low_order_embedding =
        make_shared<Embedding> (GetNDof(),
                                IntRange(low_order_space->GetNDof()),
                                IsComplex());
    // timer3.Stop();

    if (test_ho_prolongation && !prol)
      {
        prol = make_shared<H1HOProlongation> (dynamic_pointer_cast<H1HighOrderFESpace>(this->shared_from_this()));
        // prol->Update(*this);
      }
  }


  void H1HighOrderFESpace :: UpdateDofTables ()
  {
    static Timer t("H1HighOrderFESpace::UpdateDofTables"); RegionTimer reg(t);

    int dim = ma->GetDimension();
    size_t nv = ma->GetNV();
    size_t ned = (dim <= 1) ? 0 : ma->GetNEdges();
    size_t nfa = (dim <= 2) ? 0 : ma->GetNFaces();
    size_t ne = ma->GetNE();

    int hndof = nv;

    first_edge_dof.SetSize (ned+1);
    for (auto i : Range (ned))
      {
	first_edge_dof[i] = hndof;
        int oe = order_edge[i];
        if (highest_order_dc) oe--;
	if (oe > 1) hndof += oe - 1;
      }
    first_edge_dof[ned] = hndof;

    first_face_dof.SetSize (nfa+1);

    if (nfa)
      ParallelFor
        (nfa, [&] (size_t i)
         {
           int neldof = 0;             
           IVec<2> p = order_face[i];
           switch(ma->GetFaceType(i))
             {
             case ET_TRIG:
               if (p[0] > 2)
                 neldof = (p[0]-1)*(p[0]-2)/2;
               break;
             case ET_QUAD:
               if (p[0] > 1 && p[1] > 1)
                 neldof = (p[0]-1)*(p[1]-1);
               break; 
             default:
               ;
             }
           first_face_dof[i] = neldof;
         });

    // accumulate
    for (auto i : Range(nfa))
      {
        auto neldof = first_face_dof[i];
        first_face_dof[i] = hndof;
        hndof += neldof;
      }
    first_face_dof[nfa] = hndof;
    

    // compute number of element dofs ...
    first_element_dof.SetSize(ne+1);
    ParallelFor
      (ma->GetNE(VOL), [&] (size_t i)
       {
        ElementId ei(VOL,i);
        int neldof = 0;
	IVec<3> p = order_inner[i];
	switch (ma->GetElType(ei))
	  {
	  case ET_TRIG:
	    if(p[0] > 2)
	      neldof = (p[0]-1)*(p[0]-2)/2;
	    break;
	  case ET_QUAD:
	    if(p[0] > 1 && p[1] > 1)
	      neldof = (p[0]-1)*(p[1]-1);
	    break;
	  case ET_TET:
	    if(p[0] > 3)
	      neldof = (p[0]-1)*(p[0]-2)*(p[0]-3)/6;
	    break;
	  case ET_PRISM:
	    if(p[0] > 2 && p[2] > 1)
	      neldof = (p[0]-1)*(p[0]-2)*(p[2]-1)/2;
	    break;
	  case ET_PYRAMID:
	    if(p[0] > 2)
	      neldof = (p[0]-1)*(p[0]-2)*(2*p[0]-3)/6;
	    break;
	  case ET_HEX:
	    if(p[0] > 1 && p[1] > 1 && p[2] > 1) 
	      neldof = (p[0]-1)*(p[1]-1)*(p[2]-1);
	    break;
	  case ET_HEXAMID:
	    if(p[0] > 1 && p[1] > 1 && p[2] > 1) 
	      neldof = (p[0]-1)*(p[1]-1)*(p[2]-1);
	    break;
          case ET_SEGM:
            if (p[0] > 1)
	      neldof = p[0]-1;
            break;
          case ET_POINT:
            neldof = 0;
	    break;
	  }
        if (highest_order_dc)
          neldof += ElementTopology::GetNEdges(ma->GetElType(ei));
        first_element_dof[i] = neldof;        
       });

    // accumulate
    for (auto i : Range(ne))
      {
        auto neldof = first_element_dof[i];
        first_element_dof[i] = hndof;
        hndof += neldof;
      }
    first_element_dof[ne] = hndof;
    // ndof = hndof;
    SetNDof(hndof);
    
    if (print)
      {
        (*testout) << "h1 first edge = " << first_edge_dof << endl;
        (*testout) << "h1 first face = " << first_face_dof << endl;
        (*testout) << "h1 first inner = " << first_element_dof << endl;
      }

    /*
    while (ma->GetNLevels() > ndlevel.Size())
      ndlevel.Append (ndof);
    ndlevel.Last() = ndof;
    */
    if (prol)
      prol->Update(*this);
  }


  void H1HighOrderFESpace :: UpdateCouplingDofArray()
  {
    static Timer t("H1HighOrderFESpace::UpdateCouplingDofArray"); RegionTimer reg(t);    
    ctofdof.SetSize(GetNDof());


    ParallelFor
      (ma->GetNV(), [&] (size_t i)
       {
         ctofdof[i] = used_vertex[i] ? WIREBASKET_DOF : UNUSED_DOF;
       });
    /*
    ctofdof.Range(0,ma->GetNV()) = [&] (size_t i)
      { return used_vertex[i] ? WIREBASKET_DOF : UNUSED_DOF; }  | 1_tasks_per_thread;
    */

      
    int dim = ma->GetDimension();
    size_t ned = (dim <= 1) ? 0 : ma->GetNEdges();
    // for (auto edge : Range (ned))
    ParallelFor
      (ned, [&] (size_t edge)
       {
         IntRange range = GetEdgeDofs (edge);
         if (wb_edge)
           ctofdof[range] = WIREBASKET_DOF;
         else
           {
             ctofdof[range] = INTERFACE_DOF;
             if ( (wb_loedge||nodalp2) && (range.Size() > 0))
               ctofdof[range.First()] = WIREBASKET_DOF;
           }
       });

    if (ma->GetDimension() == 3)
    {
      COUPLING_TYPE face_dof_type = INTERFACE_DOF;
      if (ma->GetNE(VOL) == 0)
        face_dof_type = LOCAL_DOF;
      // for (auto face : Range (ma->GetNFaces()))
      ParallelFor
        (ma->GetNFaces(),
         [&] (size_t face)
         {
           ctofdof[GetFaceDofs(face)] = face_dof_type;
         });
    }
    
    // for (auto el : Range(ma->GetNE()))
    ParallelFor
      (ma->GetNE(),
       [&] (size_t el)
       {
         ctofdof[GetElementDofs(el)] = LOCAL_DOF;
       });
    
    if (print)
      *testout << "ctofdof: " << endl << ctofdof << endl;
  }


  void H1HighOrderFESpace :: DoArchive (Archive & archive)
  {
    low_order_space -> DoArchive (archive);
    FESpace::DoArchive(archive);
    archive & level;

    archive & NodalData<NT_EDGE> (*ma, order_edge);
    archive & NodalData<NT_FACE> (*ma, order_face);
    archive & NodalData<NT_CELL> (*ma, order_inner);

    if (archive.Input())
      UpdateDofTables();
    archive & rel_order & var_order & fixed_order & wb_loedge;
    archive & used_vertex & used_edge & used_face;
    archive & uniform_order_inner & uniform_order_face 
      & uniform_order_edge & uniform_order_quad & uniform_order_trig;
    archive & dom_order_min & dom_order_max;
    // archive & smoother;
    // archive & ndlevel;
    archive & level_adapted_order & nodalp2;
  }

  Array<MemoryUsage> H1HighOrderFESpace :: GetMemoryUsage () const
  {
    auto mu = FESpace::GetMemoryUsage();
    mu += { "H1HighOrder::order_inner", order_inner.Size()*sizeof(IVec<3,TORDER>), 1 };
    mu += { "H1HighOrder::order_face", order_face.Size()*sizeof(IVec<2,TORDER>), 1 };
    mu += { "H1HighOrder::order_edge", order_edge.Size()*sizeof(TORDER), 1 };
    return mu;
  }

  FlatArray<VorB> H1HighOrderFESpace :: GetDualShapeNodes (VorB vb) const
  {
    static VorB nodes[] = { VOL, BND, BBND, BBBND };

    if (first_edge_dof[0] == GetNDof())
      return FlatArray<VorB> (1, &nodes[ma->GetDimension() - int(vb)]);
    if (first_face_dof[0] == GetNDof())
      return FlatArray<VorB> (2, &nodes[ma->GetDimension()-1 - int(vb)]);

    return FlatArray<VorB> (ma->GetDimension()-int(vb)+1, &nodes[0]); 
  }

  
  Timer tgetfe("H1FESpace::GetFE");
  FiniteElement & H1HighOrderFESpace :: GetFE (ElementId ei, Allocator & alloc) const
  {
    // size_t tid = TaskManager::GetThreadId();
    // RegionTimer reg(tgetfe;
    
    Ngs_Element ngel = ma->GetElement(ei);
    ELEMENT_TYPE eltype = ngel.GetType();

    if (!DefinedOn (ei))
      {
        return
          SwitchET (eltype, [&] (auto et) -> FiniteElement&
                      {
                        return *new (alloc) ScalarDummyFE<et.ElementType()> ();
                      });
      }
    
    try
      {
        if (fixed_order && eltype == ET_TRIG)
          {
            switch (order)
              {
              case 1: return *(new (alloc) H1HighOrderFEFO<ET_TRIG,1> ()) -> SetVertexNumbers(ngel.vertices);
              case 2: return *(new (alloc) H1HighOrderFEFO<ET_TRIG,2> ()) -> SetVertexNumbers(ngel.vertices);
              case 3: return *(new (alloc) H1HighOrderFEFO<ET_TRIG,3> ()) -> SetVertexNumbers(ngel.vertices);
              default:
                ; 
              }
          }
        
        if (fixed_order && eltype == ET_TET)
          {
            switch (order)
              {
              case 1: return *(new (alloc) H1HighOrderFEFO<ET_TET,1> ()) -> SetVertexNumbers(ngel.vertices);
              case 2: return *(new (alloc) H1HighOrderFEFO<ET_TET,2> ()) -> SetVertexNumbers(ngel.vertices);
              case 3: return *(new (alloc) H1HighOrderFEFO<ET_TET,3> ()) -> SetVertexNumbers(ngel.vertices);
              default:
                ; 
              }
          }

        if (nodal)
          {
            return SwitchET<ET_SEGM,ET_TRIG,ET_TET> 
              (eltype,
               [&] (auto et) -> FiniteElement&
               {
                 constexpr ELEMENT_TYPE ET = et.ElementType();
                 Ngs_Element ngel = ma->GetElement<et.DIM,BND> (ei.Nr());
                 
                 auto hofe =  new (alloc) NodalHOFE<ET> (order);
                 hofe -> SetVertexNumbers (ngel.vertices);
                 return *hofe;
               });
          }
        
        auto elnr = ei.Nr();
        if (ei.IsVolume())
          {
            return SwitchET
              (eltype, [&] (auto et) -> FiniteElement&
               {
                 constexpr ELEMENT_TYPE ET = et.ElementType();
                 
                 Ngs_Element ngel = ma->GetElement<et.DIM,VOL> (elnr);
                 auto * hofe =  new (alloc) H1HighOrderFE<ET> ();
                 
                 hofe -> SetVertexNumbers (ngel.Vertices());
                 
                 switch (int(ET_trait<ET>::DIM))
                   {
                   case 0:
                     break;
                     
                   case 1:
                     {
                       hofe -> SetOrderEdge (0, order_inner[elnr][0]);
                       break;
                     }
                     
                   case 2:
                     {
                       if(nodalp2)
                         hofe -> SetNodalp2();
                       hofe -> SetOrderEdge (order_edge[ngel.Edges()] );
                       hofe -> SetOrderFace (0, order_inner[elnr]);
                       break;
                     }
                     
                   case 3: default:  
                     {
                       if(nodalp2)
                         hofe -> SetNodalp2();
                       hofe -> SetOrderEdge (order_edge[ngel.Edges()]);
                       hofe -> SetOrderFace (order_face[ngel.Faces()]);
                       hofe -> SetOrderCell (order_inner[elnr]);
                       break;
                     }
                   }
                 
                 hofe -> ComputeNDof();
                 return *hofe;
               });
          }
        else if (ei.IsBoundary())
          {
            return SwitchET<ET_POINT,ET_SEGM,ET_TRIG,ET_QUAD> 
              (eltype,
               [&] (auto et) -> FiniteElement&
               {
                 Ngs_Element ngel = ma->GetElement<et.DIM,BND> (ei.Nr());

                 // MSVC15 gets confused:
                 // auto hofe =  new (alloc) H1HighOrderFE<et.ElementType()> ();
                 
                 constexpr ELEMENT_TYPE ET = et.ElementType();
                 auto hofe =  new (alloc) H1HighOrderFE<ET> ();

                 hofe -> SetVertexNumbers (ngel.vertices);

                 if (et.DIM >= 1)
                   {
                     // hofe -> SetOrderEdge (order_edge[ngel.Edges()]);
                     auto edges = ngel.Edges();
                     for (auto i : Range(edges))
                       hofe->SetOrderEdge (i, order_edge[edges[i]] - (highest_order_dc ? 1 : 0));
                   }
                 
                 if (et.DIM >= 2)
                   hofe -> SetOrderFace (0, order_face[ma->GetSElFace(ei.Nr())]);
                 
                 hofe -> ComputeNDof();
                 return *hofe;
               });
          }
	  else if (ei.VB() == BBND)
	    {
	      switch (eltype)
		{
		case ET_POINT: return T_GetCD2FE<ET_POINT> (elnr, alloc);
		case ET_SEGM: return T_GetCD2FE<ET_SEGM> (elnr, alloc);
		default:
		  throw Exception ("illegal element in H1HoFeSpace::GetCD2FE");
		}
	    }
          else
            {
	      switch (eltype)
		{
		case ET_POINT:
                  {
                    Ngs_Element ngel = ma->GetElement<0,BBBND> (elnr);
                    
                    H1HighOrderFE<ET_POINT> * hofe = new (alloc) H1HighOrderFE<ET_POINT> ();
                    hofe -> SetVertexNumbers (ngel.vertices);
                    hofe -> ComputeNDof();
                    return *hofe;
                  }
		default:
		  throw Exception ("illegal element in H1HoFeSpace::GetCD3FE");
		}
              
            }
      }

    catch (Exception & e)
      {
        e.Append (string("in H1HoFESpace::GetElement, ei = ")+ToString(ei));
        throw;
      }
  }
  
  void H1HighOrderFESpace :: SetOrder (NodeId ni, int order) 
  {
    if (order_policy == CONSTANT_ORDER || order_policy == NODE_TYPE_ORDER)
      throw Exception("In H1HighOrderFESpace::SetOrder. Order policy is constant or node-type!");
    else if (order_policy == OLDSTYLE_ORDER)
      order_policy = VARIABLE_ORDER;
      
    if (order < 1)
      order = 1;
    
    switch (ni.GetType())
      {
      case NT_VERTEX:
        break;
      case NT_EDGE:
        if (ni.GetNr() < order_edge.Size())
          order_edge[ni.GetNr()] = order;
        break;
      case NT_FACE:
        if (ni.GetNr() < order_face.Size())
          order_face[ni.GetNr()] = order;
        break;
      case NT_CELL: case NT_ELEMENT:
        if (ni.GetNr() < order_inner.Size())
          order_inner[ni.GetNr()] = order;
        break;
      case NT_FACET:
        break;
      case NT_GLOBAL:
        break;
      }
  }
  
  int H1HighOrderFESpace :: GetOrder (NodeId ni) const
  {
    switch (ni.GetType())
      {
      case NT_VERTEX:
        return 0;
      case NT_EDGE:
        if (ni.GetNr() < order_edge.Size())
          return order_edge[ni.GetNr()];
        break;
      case NT_FACE:
        if (ni.GetNr() < order_face.Size())
          return order_face[ni.GetNr()][0];
        break;
      case NT_CELL: case NT_ELEMENT:
        if (ni.GetNr() < order_inner.Size())
          return order_inner[ni.GetNr()][0];
        break;
      case NT_FACET:
        break;
      case NT_GLOBAL:
        break;
      }
    return 0;
  }


  template <ELEMENT_TYPE ET>
  FiniteElement & H1HighOrderFESpace :: T_GetFE (int elnr, Allocator & lh) const
  {
    Ngs_Element ngel = ma->GetElement<ET_trait<ET>::DIM,VOL> (elnr);

    H1HighOrderFE<ET> * hofe =  new (lh) H1HighOrderFE<ET> ();
    
    hofe -> SetVertexNumbers (ngel.Vertices());

    switch (int(ET_trait<ET>::DIM))
      {
      case 0:
        break;
        
      case 1:
        {
          hofe -> SetOrderEdge (0, order_inner[elnr][0]);
          break;
        }

      case 2:
        {
	  if(nodalp2)
	    hofe -> SetNodalp2();
          hofe -> SetOrderEdge (order_edge[ngel.Edges()] );
          hofe -> SetOrderFace (0, order_inner[elnr]);
          break;
        }

      case 3: default:  
        {
	  if(nodalp2)
	    hofe -> SetNodalp2();
          hofe -> SetOrderEdge (order_edge[ngel.Edges()]);
          hofe -> SetOrderFace (order_face[ngel.Faces()]);
          hofe -> SetOrderCell (order_inner[elnr]);
          break;
        }
      }

    hofe -> ComputeNDof();
    return *hofe;
  }



  template <ELEMENT_TYPE ET>
  FiniteElement & H1HighOrderFESpace :: T_GetSFE (int elnr, Allocator & lh) const
  {
    Ngs_Element ngel = ma->GetElement<ET_trait<ET>::DIM,BND> (elnr);

    H1HighOrderFE<ET> * hofe =  new (lh) H1HighOrderFE<ET> ();
    
    hofe -> SetVertexNumbers (ngel.vertices);
    
    switch (int (ET_trait<ET>::DIM))
      {
      case 0:
        {
          break;
        }

      case 1:
        {
          hofe -> SetOrderEdge ( order_edge[ngel.Edges()] );
          break;
        }

      case 2: default:  
        {
          hofe -> SetOrderEdge (order_edge[ngel.Edges()]);
	  hofe -> SetOrderFace (0, order_face[ma->GetSElFace(elnr)]);
          break;
        }
      }

    hofe -> ComputeNDof();
    return *hofe;
  }

  template <ELEMENT_TYPE ET>
  FiniteElement & H1HighOrderFESpace :: T_GetCD2FE(int elnr, Allocator & lh) const
  {
    Ngs_Element ngel = ma->GetElement<ET_trait<ET>::DIM,BBND> (elnr);

    H1HighOrderFE<ET> * hofe = new (lh) H1HighOrderFE<ET> ();
    hofe -> SetVertexNumbers (ngel.vertices);
    if(int(ET_trait<ET>::DIM)==1)
      {
	hofe->SetOrderEdge(order_edge[ngel.Edges()]);
      }
    hofe->ComputeNDof();
    return *hofe;
  }
 
  /*
  size_t H1HighOrderFESpace :: GetNDofLevel (int alevel) const
  {
    return ndlevel[alevel];
  }
  */

  void H1HighOrderFESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    Ngs_Element ngel = ma->GetElement(ei);
    dnums.SetSize0();

    if (!DefinedOn (ei))
      return;
    
    dnums = ngel.Vertices();
    if (fixed_order && order==1) return;

    IntRange eldofs;
    if (ei.IsVolume())
      eldofs = GetElementDofs (ei.Nr());
    
    if (ma->GetDimension() >= 2)
      for (auto edge : ngel.Edges())
        {
          dnums += GetEdgeDofs (edge);
          if (ei.IsVolume() && highest_order_dc)
            {
              dnums += eldofs.First();
              eldofs.First()++;
            }
        }

    if (ma->GetDimension() == 3)
      for (auto face : ngel.Faces())
        dnums += GetFaceDofs (face);
    if (ei.IsVolume())
      dnums += eldofs;
  }


  void H1HighOrderFESpace :: 
  GetDofRanges (ElementId ei, Array<IntRange> & dranges) const
  {
    dranges.SetSize(0);

    if (!DefinedOn (ei)) return;
    Ngs_Element ngel = ma->GetElement(ei);

    for (int i = 0; i < ngel.vertices.Size(); i++)
      dranges.Append (IntRange (ngel.vertices[i], ngel.vertices[i]+1));
         
    if (ma->GetDimension() >= 2)
      for (int i = 0; i < ngel.edges.Size(); i++)
        if (GetEdgeDofs (ngel.edges[i]).Size())
          dranges += IntRange (GetEdgeDofs (ngel.edges[i]));

    if (ma->GetDimension() == 3)
      for (int i = 0; i < ngel.faces.Size(); i++)
        if (GetFaceDofs (ngel.faces[i]).Size())
          dranges += IntRange (GetFaceDofs (ngel.faces[i]));

    if (ei.IsVolume())
      if (GetElementDofs (ei.Nr()).Size())
        dranges += IntRange (GetElementDofs (ei.Nr()));
  }



  
  void H1HighOrderFESpace :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  {
    dnums.SetSize(1);
    dnums[0] = vnr;
  }
  
  
  void H1HighOrderFESpace :: GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  {
    dnums = GetEdgeDofs (ednr);
  }

  void H1HighOrderFESpace :: GetFaceDofNrs (int fanr, Array<int> & dnums) const
  {
    dnums.SetSize0();
    if (ma->GetDimension() < 3) return;
    dnums = GetFaceDofs (fanr);
  }

  void H1HighOrderFESpace :: GetInnerDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums = GetElementDofs (elnr);
  }

  
  shared_ptr<Table<int>> H1HighOrderFESpace :: 
  CreateSmoothingBlocks (const Flags & precflags) const
  {
    static Timer t("H1HighOrderFESpace :: CreateSmoothingBlocks");
    RegionTimer reg(t);

    bool eliminate_internal = precflags.GetDefineFlag("eliminate_internal");
    bool subassembled = precflags.GetDefineFlag("subassembled");
    // smoothing_types: 
    // 1: 2d V + E + I 
    // 2: 2d VE + I 
    // 3: V + E + F + I 
    // 4: VE + F + I 
    // 5: VE + FI 
    // 6: VEF + I 
    // 7: VEFI 
    // 8: V + E + FI 
    // 9: V + EF + I 
    // 10:V + EF + FI 

    int smoothing_type = int(precflags.GetNumFlag("blocktype",4)); 
    if (subassembled) smoothing_type = 51;

    size_t nv = ma->GetNV();
    size_t ned = ma->GetNEdges();
    size_t nfa = (ma->GetDimension() == 2) ? 0 : ma->GetNFaces();
    size_t ni = (eliminate_internal) ? 0 : ma->GetNE(); 
   
    cout << IM(4) << " blocktype " << smoothing_type << endl;
    // cout << " Use H1-Block Smoother:  "; 

    FilteredTableCreator creator(GetFreeDofs().get());
    for ( ; !creator.Done(); creator++)
      {
	switch (smoothing_type)
	  {

	  case 1:  // 2d: V + E + I
		
	    if (creator.GetMode() == 1)
              {
                cout << " V + E + I " << endl;
                creator.SetSize(nv+ned+ni);
                continue;
              }
            /*
	    for (size_t i = 0; i < nv; i++)
	      creator.Add (i, i);
	    for (size_t i = 0; i < ned; i++)
	      creator.Add (nv+i, GetEdgeDofs(i));
	    for (size_t i = 0; i < ni; i++)
	      creator.Add (nv+ned+i, GetElementDofs(i));
            */

            ParallelFor (nv, [&creator] (size_t i)
                         { creator.Add (i,i); });
            ParallelFor (ned, [&creator,nv,this] (size_t i)
                         { creator.Add (nv+i, GetEdgeDofs(i)); });
            ParallelFor (ni, [&creator,nv,ned,this] (size_t i)
                         { creator.Add (nv+ned+i, GetElementDofs(i)); });
	    break; 
		
	  case 2: // 2d VE + I

	    if (creator.GetMode() == 1)
              {
                cout << " 2d VE + I " << endl;
                creator.SetSize(nv+ned+ni);
                continue;
              }

	    for (int i = 0; i < nv; i++)
	      creator.Add(i, i);

	    for (auto i : Range(ned))
              for (auto v : ma->GetEdgePNums(i))
                creator.Add (v, GetEdgeDofs(i));
		
	    for (int i = 0; i < ni; i++)
	      creator.Add (nv+ned+i, GetElementDofs(i));
		
	    break;


	  case 3: // V + E + F + I
		
	    if (creator.GetMode() == 1)
	      cout << " V + E + F + I " << endl; 

	    for (int i = 0; i < nv; i++)
	      creator.Add(i, i);

	    for (int i = 0; i < ned; i++)
	      creator.Add (nv+i, GetEdgeDofs(i));

	    for (int i = 0; i < nfa; i++)
	      creator.Add(nv+ned+i, GetFaceDofs(i));

	    for (int i = 0; i < ni; i++)
	      creator.Add (nv+ned+nfa+i, GetElementDofs(i));
	    
	    break; 

	  case 4: // VE + F + I
		
            if (creator.GetMode() == 1)
              {
                cout << IM(4) << " VE + F + I " << endl;
                creator.SetSize(nv+nfa+ni);
                break;
              }
            
            ParallelFor (nv, [&creator] (size_t i)
              { creator.Add(i, i); });
		
            ParallelFor (ned, [&creator,this,&ma=*ma] (size_t i)
              {
                for (auto v : ma.GetEdgePNums(i))
                  creator.Add (v, GetEdgeDofs(i));
              });
	    
            for (int i = 0; i < nfa; i++)
              creator.Add(nv+i, GetFaceDofs(i));
	    
            ParallelFor (ni, [&] (size_t i)
              { creator.Add (nv+nfa+i, GetElementDofs(i)); });
		
	    break; 

	  case 5: // VE + FI

	    if (creator.GetMode() == 1)
	      cout << " VE + FI " << endl; 

	    for (int i = 0; i < nv; i++)
	      creator.Add(i, i);
		
	    for (int i = 0; i < ned; i++)
              for (auto v : ma->GetEdgePNums(i))
                creator.Add (v, GetEdgeDofs(i));
	    
	    for (int i = 0; i < nfa; i++)
	      creator.Add(nv+i, GetFaceDofs(i));

	    for (size_t i = 0; i < ni; i++)
              for (auto f : ma->GetElement( { VOL, i } ).Faces())
                creator.Add (nv+f, GetElementDofs(i));                  

	    break; 


	  case 6: // VEF + I

	    if (creator.GetMode() == 1)
	      cout << " VEF + I " << endl; 

	    for (int i = 0; i < nv; i++)
	      creator.Add (i, i);
		
	    for (int i = 0; i < ned; i++)
	      {
		Ng_Node<1> edge = ma->GetNode<1> (i);
		for (int k = 0; k < 2; k++)
		  creator.Add (edge.vertices[k], GetEdgeDofs(i));
	      }
	    
		
	    for (int i = 0; i < nfa; i++)
	      {
		Ng_Node<2> face = ma->GetNode<2> (i);
		for (int k = 0; k < face.vertices.Size(); k++)
		  creator.Add (face.vertices[k], GetFaceDofs(i));
	      }
	    
	    for (int i = 0; i < ni; i++)
	      creator.Add (nv+i, GetElementDofs(i));
		
	    break; 

	  case 7: // VEFI

	    if (creator.GetMode() == 1)
	      cout << " VEFI " << endl; 

	    for (int i = 0; i < nv; i++)
	      creator.Add (i, i);
		
	    for (int i = 0; i < ned; i++)
	      {
		Ng_Node<1> edge = ma->GetNode<1> (i);
		for (int k = 0; k < 2; k++)
		  creator.Add (edge.vertices[k], GetEdgeDofs(i));
	      }

	    for (int i = 0; i < nfa; i++)
	      {
		Ng_Node<2> face = ma->GetNode<2> (i);
		for (int k = 0; k < face.vertices.Size(); k++)
		  creator.Add (face.vertices[k], GetFaceDofs(i));
	      }
	    
	    for (int i = 0; i < ni; i++)
              for (auto v : ma->GetElement(ElementId(VOL,i)).Vertices())
                creator.Add (v, GetElementDofs(i));

	    break; 
	    
	  case 8: // V + E + FI
	    if (creator.GetMode() == 1)
	      cout << " V + E + FI " << endl; 
		
	    for (int i = 0; i < nv; i++)
	      creator.Add (i, i);
		
	    for (int i = 0; i < ned; i++)
	      creator.Add (nv+i, GetEdgeDofs(i));
		
	    for (int i = 0; i < nfa; i++)
	      creator.Add(nv+ned+i, GetFaceDofs(i));
		
	    for (int i = 0; i < ni; i++)
	      for (int f : ma->GetElement(ElementId(VOL,i)).Faces())
		creator.Add (nv+ned+f, GetElementDofs(i));                  

	    break;


	  case 9: // V + EF + I
	    {
	      if (creator.GetMode() == 1)
                {
                  cout << " V + EF + I " << endl; 
		  creator.SetSize(nv+ned+ni);
                  continue;
                }

              ParallelFor (nv, [&] (size_t i)
                           { creator.Add(i,i); });
              
              ParallelFor (ned, [&] (size_t i) 
                           { creator.Add(nv+i,GetEdgeDofs(i)); });
              
              ParallelFor (nfa, [&] (size_t i) 
                           {
                             for (auto edge : ma->GetFaceEdges(i))
                               creator.Add (nv+edge, GetFaceDofs(i));                                   
                           });
              
              ParallelFor (ni, [&] (size_t i) 
                           { creator.Add(nv+ned+i,GetElementDofs(i)); });
              
              
                  /*
                  for (int i = 0; i < nv; i++)
                    creator.Add (i, i);
                  
                  for (int i = 0; i < ned; i++)
                    creator.Add (nv+i, GetEdgeDofs(i));
                  
                  for (int i = 0; i < nfa; i++)
                    {
                      // Ng_Node<2> face = ma->GetNode<2> (i);
                      // for (int k = 0; k < face.edges.Size(); k++)
                      // creator.Add (face.edges[k], GetFaceDofs(i));
                      
                      ArrayMem<int,4> f2ed;
                      ma->GetFaceEdges (i, f2ed);
                      for (int k = 0; k < f2ed.Size(); k++)
                        creator.Add (nv+f2ed[k], GetFaceDofs(i));
                    }
                  
                  for (int i = 0; i < ni; i++)
                    creator.Add (nv+ned+i, GetElementDofs(i));
                  */
            break; 
          }
        

	  case 10: 
	    if (creator.GetMode() == 1)
	      cout << " V + E + F +I (Cluster)  " << endl; 

	    for (int i = 0; i < nv; i++)
	      creator.Add(ma->GetClusterRepVertex(i), i);

	    for (int i = 0; i < ned; i++)
	      creator.Add (ma->GetClusterRepEdge(i), GetEdgeDofs(i));

	    for (int i = 0; i < nfa; i++)
	      creator.Add(ma->GetClusterRepFace(i), GetFaceDofs(i));

	    for (int i = 0; i < ni; i++)
	      creator.Add (ma->GetClusterRepElement(i), GetElementDofs(i));

	    break; 

	  case 11: 
	    if (creator.GetMode() == 1)
	      cout << " 2d VEI " << endl; 

	    for (int i = 0; i < nv; i++)
		creator.Add (i, i);
		
	    for (int i = 0; i < ned; i++)
	      {
		Ng_Node<1> edge = ma->GetNode<1> (i);
		for (int k = 0; k < 2; k++)
		  creator.Add (edge.vertices[k], GetEdgeDofs(i));
	      }

	    for (int i = 0; i < ni; i++)
	      {
		const Ngs_Element & ngel = ma->GetElement(ElementId(VOL,i));
		for (int j = 0; j < ngel.vertices.Size(); j++)
		  creator.Add (ngel.vertices[j], GetElementDofs(i));
	      }

	    break;


	  case 12: 
	    if (creator.GetMode() == 1)
	      cout << " VEFI Cluster " << endl; 

	    for (int i = 0; i < nv; i++)
	      creator.Add(ma->GetClusterRepVertex(i), i);

	    for (int i = 0; i < ned; i++)
	      {
		Ng_Node<1> edge = ma->GetNode<1> (i);
		int rep[2];
		for (int k = 0; k < 2; k++)
		  rep[k] = ma->GetClusterRepVertex(edge.vertices[k]);
		
		creator.Add (rep[0], GetEdgeDofs(i));
		if (rep[0] != rep[1])
		  creator.Add (rep[1], GetEdgeDofs(i));
	      }

	    for (int i = 0; i < nfa; i++)
	      {
		Ng_Node<2> face = ma->GetNode<2> (i);
		int rep[4];
		
		for (int k = 0; k < face.vertices.Size(); k++)
		  {
		    rep[k] = ma->GetClusterRepVertex(face.vertices[k]);
		    
		    bool ok = true;
		    for (int j = 0; j < k; j++)
		      if (rep[j] == rep[k]) ok = false;
		    if (ok) creator.Add (rep[k], GetFaceDofs(i));
		  }
	      }

	    for (int i = 0; i < ni; i++)
	      {
		Ngs_Element ngel = ma->GetElement (ElementId(VOL,i));
		int rep[8];
		      
		for (int k = 0; k < ngel.vertices.Size(); k++)
		  {
		    rep[k] = ma->GetClusterRepVertex(ngel.vertices[k]);
			
		    bool ok = true;
		    for (int j = 0; j < k; j++)
		      if (rep[j] == rep[k]) ok = false;
		    if (ok) creator.Add (rep[k], GetElementDofs(i));
		  }
	      }
	    break; 

	  case 20: // VE + EF + I
	    {
	      if (creator.GetMode() == 1)
		cout << "VE + EF + I" << endl;
		  
	      for (int i = 0; i < nv; i++)
		creator.Add (i, i);
		  
	      for (int i = 0; i < ned; i++)
		{
		  creator.Add (nv+i, GetEdgeDofs(i));
		  Ng_Node<1> edge = ma->GetNode<1> (i);
		  for (int k = 0; k < 2; k++)
		    creator.Add (edge.vertices[k], GetEdgeDofs(i));
		}
	      
	      for (auto face : Range(nfa))
                for (auto edge : ma->GetFaceEdges(face))
                  creator.Add (nv+edge, GetFaceDofs(face));
	      
	      for (int i = 0; i < ni; i++)
		creator.Add (nv+ned+i, GetElementDofs(i));
	      break;
	    }
	    
	  case 21: // E + F 
	  {
	    if (creator.GetMode() == 1)
	      cout << "Helmholtz" << endl;
		
	    int ds_order = int(precflags.GetNumFlag ("ds_order", 1));
	    cout << "ds_order = " << ds_order << endl;
		
	    for (int i = 0; i < ned; i++)
	      {
		int first = first_edge_dof[i] + ds_order - 1;
		int ndof = first_edge_dof[i+1]-first;
		for (int j = 0; j < ndof; j++)
		  creator.Add (i, first+j);
	      }
	    for (int i = 0; i < nfa; i++)
	      creator.Add(ned+i, GetFaceDofs(i));
		
	    break; 
	  }

	  case 51: 
	  //for BDDC: we have only the condensated (after subassembling) dofs, 
	  //and build patches around each vertex Vertices + Edges
	    
	    if (creator.GetMode() == 1)
	      cout << "BDDC-Edges-around-Vertex-Block" << endl;

	    for (int i = 0; i < nv; i++)
	      creator.Add (i, i);
	    
	    for (int i = 0; i < ned; i++)
	      {
		Ng_Node<1> edge = ma->GetNode<1> (i);
		IntRange edgedofs = GetEdgeDofs(i);

		for (int k = 0; k < 2; k++)
		  for (int l = edgedofs.First(); l < edgedofs.Next(); l++)
		    if (ctofdof[l] == WIREBASKET_DOF)
		      creator.Add (edge.vertices[k], l);
	      }
	    
	    break; 	    
	  }
      }
    // return shared_ptr<Table<int>> (creator.GetTable());
    return make_shared<Table<int>> (creator.MoveTable());
  }
    


  shared_ptr<Array<int>>
  H1HighOrderFESpace :: CreateDirectSolverClusters (const Flags & flags) const
  {
    if (flags.GetDefineFlag("subassembled"))
    {
	cout << "creating bddc-coarse grid(vertices)" << endl;
        auto spclusters = make_shared<Array<int>> (GetNDof());
	Array<int> & clusters = *spclusters;
	clusters = 0;
	int nv = ma->GetNV();
	for (int i = 0; i < nv; i++)
	  if (!IsDirichletVertex(i))
	    clusters[i] = 1;		
	return spclusters;
    }
    
    if (flags.NumFlagDefined ("ds_order"))
      {
	int ds_order = int (flags.GetNumFlag ("ds_order", 1));

        auto spclusters = make_shared<Array<int>> (GetNDof());
	Array<int> & clusters = *spclusters;
	clusters = 0;
	
	int ned = ma->GetNEdges();
	int nv = ma->GetNV();

	// for (int i = 0; i < nv; i++)
        //   clusters[i] = 1;
        clusters.Range(0,nv) = 1;

	for (int i = 0; i < ned; i++)
	  {
	    int first = first_edge_dof[i];
	    int next = first_edge_dof[i+1];
	    for (int j = 0; (j+2 <= ds_order) && (first+j < next) ; j++)
	      clusters[first+j] = 1;
	  }

	/*
	  int nfa = ma->GetNFaces();
	  for (int i = 0; i < nfa; i++)
	  {
	  int first = first_face_dof[i];
	  int next = first_face_dof[i+1];
	  int p = order_face[i][0];
	    
	  // if (usegrad_face[i])
	  int ii = 0;
	  for (int j = 0; j <= p-2; j++)
	  for (int k = 0; k <= p-2-j; k++, ii++)
	  if (j+k+2 <= ds_order)
	  clusters[first+ii] = 1;
	    
	  // other combination
	  for (int j = 0; j <= p-2; j++)
	  for (int k = 0; k <= p-2-j; k++, ii++)
	  if (j+k+2 <= ds_order)
	  clusters[first+ii] = 1;
	    
	  // type 3
	  for (int j = 0; j <= p-2; j++, ii++)
	  if (j+2 <= ds_order)
	  clusters[first+ii] = 1;
	  }
	*/

	return spclusters;
      }




    // return 0;

    // int nv = ma->GetNV();
    // int nd = GetNDof();
    int ne = ma->GetNE();
    auto spclusters = make_shared<Array<int>> (GetNDof());
    Array<int> & clusters = *spclusters;
    clusters = 0;

    // all vertices in global space
    /*
      for (int i = 0; i < nv; i++)
      clusters[i] = 1;

      // all edges in direct solver cluster

      for (i = first_edge_dof[0]; i < first_edge_dof.Last(); i++)
      clusters[i] = 1;
      return &clusters;

    */
   
    // All Vertical Edges in one Cluster for Hex and Prism (-> 2d Problems !) 

    //Array<int> & clusters = *new Array<int> (nd);
    //clusters = 0;


    
    for (int i = 0; i < ne; i++)
      {
        ElementId ei(VOL, i);
	if (ma->GetElType(ei) == ET_PRISM)
	  {
	    auto ednums = ma->GetElEdges (ei);
	    for (int j = 6; j < 9; j++)  //vertical Edges 
              clusters[GetEdgeDofs(ednums[j])] = 2;
	    
	    auto fnums = ma->GetElFaces(ei); // vertical faces 
	    for (int j =2;j<5;j++) 
	      {
                clusters[GetFaceDofs(fnums[j])] = 0;
                /*
		int first = first_face_dof[fnums[j]]; 
		int next = first_face_dof[fnums[j]+1]; 
		
		for (k=first; k < next; k++) 
		  clusters[k]=0; 
                */

		//IVec<2> p = order_face[fnums[j]];
		//for(k=first + 2*(p[0]+1)*(p[1]+1);k<next;k++)
		//  clusters[k]=3;  
	      }
	  }

	else if (ma->GetElType(ei) == ET_HEX)  
	  {
	    auto ednums = ma->GetElEdges (ei);
	    for (int j = 8; j < 12; j++) //vertical edges
              clusters[GetEdgeDofs(ednums[j])] = 2;

	    auto fnums = ma->GetElFaces(ei); // vertical faces 
	    for (int j = 2; j < 6; j++) 
              clusters[GetFaceDofs(fnums[j])] = 3;
	  } 
      }

    Array<int> ednums; 
    for (int i =0; directsolverclustered.Size() > 0 && i<ne; i++)
      {
        ElementId ei(VOL,i);
	if(directsolverclustered[ma->GetElIndex(ei)])
	  {
	    GetDofNrs(i,ednums);
	    for (int k = 0; k<ednums.Size(); k++)
	      {
		clusters[ednums[k]] = 4;
	      }
	  }
      }

   
    clusters[adddirectsolverdofs] = 5;
    /*
    for (int i =0; i< adddirectsolverdofs.Size(); i++)
      {
	clusters[adddirectsolverdofs[i]] = 5;
      }
    */
    const int stdoffset = 6;

    /*
      (*testout) << "directvertexclusters" << directvertexclusters << endl;
      (*testout) << "directedgeclusters" << directedgeclusters << endl;
      (*testout) << "directfaceclusters" << directfaceclusters << endl;
      (*testout) << "directelementclusters" << directelementclusters << endl;
    */

    for (int i = 0; i<directvertexclusters.Size(); i++)
      if(directvertexclusters[i] >= 0)
	clusters[i] = directvertexclusters[i] + stdoffset;

    for (int i = 0; i<directedgeclusters.Size(); i++)
      if(directedgeclusters[i] >= 0)
	for(int j = first_edge_dof[i]; j<first_edge_dof[i+1]; j++)
	  clusters[j] = directedgeclusters[i] + stdoffset;

    for (int i = 0; i<directfaceclusters.Size(); i++)
      if(directfaceclusters[i] >= 0)
	for(int j = first_face_dof[i]; j<first_face_dof[i+1]; j++)
	  clusters[j] = directfaceclusters[i] + stdoffset;
	  
    for (int i = 0; i<directelementclusters.Size(); i++)
      if(directelementclusters[i] >= 0)
	for(int j = first_element_dof[i]; j<first_element_dof[i+1]; j++)
	  clusters[j] = directelementclusters[i] + stdoffset;


    //    (*testout) << "clusters " << clusters << endl;
    
    bool nonzero = false;
    for (int i = 0; !nonzero && i < clusters.Size(); i++)
      if (clusters[i]) nonzero = true;
    if (!nonzero)
      {
	return nullptr;
      }


    // filter with freedofs
    BitArray & free = *GetFreeDofs(flags.GetDefineFlag("eliminate_internal"));
    for (size_t i = 0; i < clusters.Size(); i++)
      if (!free.Test(i))
        clusters[i] = 0;
    
    return spclusters;
  }

  template<>
  shared_ptr<FESpace> RegisterFESpace<H1HighOrderFESpace> :: 
  Create (shared_ptr<MeshAccess> ma, const Flags & flags)
  {
    // we will check the -periodic flag here
    return make_shared<H1HighOrderFESpace> (ma, flags);
  }



  template <class Tx, class T>
  inline void LowEnergyVertexPolynomials2D  (int n, Tx x, T & values)
  {
    JacobiPolynomial (n, x, 0, -1, values);
    Tx sum1 = 0.0, sum2 = 0.0;
    for (int i = 1; i <= n; i++)
      {
	sum1 += 1.0/i;
	sum2 += values[i] / i;
	values[i] = sum2/sum1;
      }
    values[0] = 1;
  }

  template <class Tx, class T>
  inline void LowEnergyVertexPolynomials3D  (int n, Tx x, T & values)
  {
    JacobiPolynomial (n, x, 1, -1, values);
    Tx sum = 0.0;
    for (int i = 1; i <= n; i++)
      {
	sum += (2.0*i+1)/(i+1) * values[i];
	values[i] = 1.0/(i*(i+2)) * sum;
      }
    values[0] = 1;
  }


  class LowEnergyTrig : public ScalarFiniteElement<2>
  {
  public:
    LowEnergyTrig (int order)
      : ScalarFiniteElement<2> (3, order) { ; }

    virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }

    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceVector<> shape) const
    { T_CalcShape (ip(0), ip(1), shape); }
  
    virtual void CalcDShape (const IntegrationPoint & ip, 
                             BareSliceMatrix<> dshape) const
    {
      AutoDiff<2> adx (ip(0), 0);
      AutoDiff<2> ady (ip(1), 1);
      Vector<AutoDiff<2> > shapearray(ndof);
      T_CalcShape<AutoDiff<2>> (adx, ady, shapearray);
      for (int i = 0; i < ndof; i++)
        {
          dshape(i, 0) = shapearray[i].DValue(0);
          dshape(i, 1) = shapearray[i].DValue(1);
        }
    }

  private:
    template <class T>
    void T_CalcShape (const T & x, const T & y, BareSliceVector<T> shape) const
    {
      T lami[3] = { x, y, 1-x-y };
      ArrayMem<T,100> allvalues(order+1);
      for (int i = 0; i < 3; i++)
        {
          LowEnergyVertexPolynomials3D  (order, lami[i], allvalues);
          shape[i] = allvalues[order];
        }
    }
  };


  class LowEnergyTet : public ScalarFiniteElement<3>
  {
  public:
    LowEnergyTet (int order)
      : ScalarFiniteElement<3> (4, order) { ; }

    virtual ELEMENT_TYPE ElementType() const { return ET_TET; }

    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceVector<> shape) const
    { T_CalcShape (ip(0), ip(1), ip(2), shape); }
  
    virtual void CalcDShape (const IntegrationPoint & ip, 
                             BareSliceMatrix<> dshape) const
    {
      AutoDiff<3> adx (ip(0), 0);
      AutoDiff<3> ady (ip(1), 1);
      AutoDiff<3> adz (ip(2), 2);      
      Vector<AutoDiff<3> > shapearray(ndof);
      T_CalcShape<AutoDiff<3>> (adx, ady, adz, shapearray);
      for (int i = 0; i < ndof; i++)
        {
          dshape(i, 0) = shapearray[i].DValue(0);
          dshape(i, 1) = shapearray[i].DValue(1);
          dshape(i, 2) = shapearray[i].DValue(2);
        }
    }

  private:
    template <class T>
    void T_CalcShape (const T & x, const T & y, const T & z, BareSliceVector<T> shape) const
    {
      T lami[4] = { x, y, z, 1-x-y-z };
      ArrayMem<T,100> allvalues(order+1);
      for (int i = 0; i < 4; i++)
        {
          LowEnergyVertexPolynomials3D  (order, lami[i], allvalues);
          shape[i] = allvalues[order];
        }
    }
  };

  

  class LowEnergyVertexFESpace : public FESpace
  {
    int order;
  public:
    LowEnergyVertexFESpace (shared_ptr<MeshAccess> ama, const Flags & flags)
      : FESpace (ama, flags)
    {
      order = int(flags.GetNumFlag ("order", 2));
      if (ma->GetDimension() == 2)
        {
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
          flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
          evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();
        }
      else
        {
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
          flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<3>>>();
          evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<3>>>();
        }
        
    }

    virtual string GetClassName () const
    {
      return "LowEnergyVertexFESpace";
    }

    virtual void Update()
    {
      SetNDof(ma->GetNV());
    }

    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const
    {
      dnums.SetSize0();
      dnums += ma->GetElement(ei).Vertices();
    }
    
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const
    {
      switch (ma->GetElement(ei).GetType())
        {
        case ET_TRIG: return *new(alloc) LowEnergyTrig(order);
        case ET_TET: return *new(alloc) LowEnergyTet(order);
        default: throw Exception("not supported");
        }
    }
  };

  static RegisterFESpace<LowEnergyVertexFESpace> initle ("lowenergyvertex");






  


  template <int DIM_SPC>
  class DiffOpDivFreeReconstructVectorH1 : public DiffOp<DiffOpDivFreeReconstructVectorH1<DIM_SPC> >
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = DIM_SPC };
    enum { DIM_ELEMENT = DIM_SPC };
    enum { DIM_DMAT = DIM_SPC };
    enum { DIFFORDER = 0 };

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const MIP & mip,
                                MAT && bmat, LocalHeap & lh)
    {
      auto & fel = static_cast<const VectorFiniteElement&> (bfel);
      auto & fel_u = static_cast<const ScalarFiniteElement<DIM_SPACE>&> (fel[0]);
      
      int order = fel_u.Order();
      auto & trafo = mip.GetTransformation();
      
      HDivHighOrderFE<ET_TRIG> fel_hdiv(order-1);
      L2HighOrderFE<ET_TRIG> fel_l2(order-2);
      L2HighOrderFE<ET_TRIG> fel_koszul( max(order-3, -1) );
      // FacetFE<ET_TRIG> fel_facet;
      FacetFE<ET_TRIG> & fel_facet = *new (lh) FacetFE<ET_TRIG>;
      fel_facet.SetOrder(order-1);
      Array<int> vnums = { 1, 2, 3 } ;
      fel_facet.SetVertexNumbers(vnums);
      fel_facet.ComputeNDof();

      Matrix<> mat(fel_hdiv.GetNDof());
      Matrix<> rhs(fel_hdiv.GetNDof(), DIM_SPACE*fel_u.GetNDof());

      auto eltype = trafo.GetElementType();
      int nfacet = ElementTopology::GetNFacets(eltype);
      Facet2ElementTrafo transform(eltype); 

      IntRange r_facet(0, fel_facet.GetNDof());
      IntRange r_div(r_facet.Next(), r_facet.Next() + fel_l2.GetNDof() - 1);
      IntRange r_koszul(r_div.Next(), fel_hdiv.GetNDof());

      /*
      *testout << "r_facet = " << r_facet << endl;
      *testout << "r_div = " << r_div << endl;
      *testout << "r_koszul = " << r_koszul << endl;
      */
      
      mat = 0;
      rhs = 0;

      FlatMatrix<> shape_hdiv(fel_hdiv.GetNDof(), 2, lh); 
      FlatVector<> shape_hdivn(fel_hdiv.GetNDof(), lh);
      FlatVector<> shape_u(fel_u.GetNDof(),lh);

      // edge moments
      for (int k = 0; k < nfacet; k++)
        {
          IntRange r_facet_k = fel_facet.GetFacetDofs(k);
          FlatVector<> shape_facet(r_facet_k.Size(), lh);

          HeapReset hr(lh);
          ngfem::ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);
          
          IntegrationRule ir_facet(etfacet, 2*fel.Order());
          IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
          
          MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> mir(ir_facet_vol, trafo, lh);
          mir.ComputeNormalsAndMeasure (eltype, k);
       
          for (int j = 0; j < ir_facet.Size(); j++)
            {
              fel_hdiv.CalcMappedShape (mir[j], shape_hdiv);
              shape_hdivn = shape_hdiv * mir[j].GetNV();
              shape_facet = 0;
              fel_facet.CalcFacetShapeVolIP (k, ir_facet_vol[j], shape_facet);
              mat.Rows(r_facet_k) += mir[j].GetWeight() * shape_facet * Trans(shape_hdivn);

              fel_u.CalcShape(ir_facet_vol[j], shape_u);
              for (int l = 0; l < DIM_SPACE; l++)
                rhs.Rows(r_facet_k).Cols(l*fel_u.GetNDof(), (l+1)*fel_u.GetNDof()) +=
                  (mir[j].GetWeight()*mir[j].GetNV()(l)) * shape_facet * Trans(shape_u);
            }
        }


      IntegrationRule ir(eltype, 2*fel.Order());
      MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> mir(ir, trafo, lh);      

      FlatVector<> div_shape(fel_hdiv.GetNDof(), lh);
      FlatVector<> shape_l2(fel_l2.GetNDof(), lh);
      FlatVector<> shape_koszul(fel_koszul.GetNDof(), lh);
      FlatMatrix<> grad_u(fel_u.GetNDof(), DIM_SPACE, lh);
      
      for (int j = 0; j < ir.Size(); j++)
        {
          fel_hdiv.CalcMappedShape (mir[j], shape_hdiv);
          fel_hdiv.CalcMappedDivShape (mir[j], div_shape);
          fel_u.CalcMappedDShape (mir[j], grad_u);
          fel_u.CalcShape (ir[j], shape_u);
	  
	  //shape_l2 = 0;
	  if(fel_l2.GetNDof() > 0)
	    fel_l2.CalcShape(ir[j], shape_l2);

	  // shape_koszul = 0;
	  if(fel_koszul.GetNDof() > 0)
	    fel_koszul.CalcShape(ir[j], shape_koszul);

          mat.Rows(r_div) += mir[j].GetWeight() * shape_l2.Range(1, shape_l2.Size()) * Trans(div_shape);
          for (int l = 0; l < DIM_SPACE; l++)
            rhs.Rows(r_div).Cols(l*fel_u.GetNDof(), (l+1)*fel_u.GetNDof()) +=
              mir[j].GetWeight() * shape_l2.Range(1, shape_l2.Size()) * Trans(grad_u.Col(l));

          Vec<2> rel_pos = mir[j].Point()-mir[0].Point();
          Vec<2> ymx (rel_pos(1), -rel_pos(0));
          
          mat.Rows(r_koszul) += mir[j].GetWeight() * shape_koszul * Trans(shape_hdiv * ymx);
          for (int l = 0; l < DIM_SPACE; l++)
            rhs.Rows(r_koszul).Cols(l*fel_u.GetNDof(), (l+1)*fel_u.GetNDof()) +=
              (mir[j].GetWeight()*ymx(l)) * shape_koszul * Trans(shape_u);
        }

      // *testout << "mat = " << endl << mat << endl;
      // *testout << "rhs = " << endl << rhs << endl;

      CalcInverse (mat);
      Matrix<> prod = mat * rhs;
      fel_hdiv.CalcMappedShape (mip, shape_hdiv);

      bmat = Trans(shape_hdiv) * prod;
    }
  };



  
  
  template <int D>
  class NGS_DLL_HEADER VectorH1MassIntegrator 
    : public T_BDBIntegrator<DiffOpIdVectorH1<D>, DiagDMat<D> >
  {
  public:
    using T_BDBIntegrator<DiffOpIdVectorH1<D>, DiagDMat<D>>::T_BDBIntegrator;
  };
  


  VectorH1FESpace::VectorH1FESpace (shared_ptr<MeshAccess> ama, const Flags & flags, 
                     bool checkflags)
      : CompoundFESpace(ama, flags)
    {
      type = "VectorH1";
      Array<string> dirichlet_comp;
      string dirnames[] = { "dirichletx", "dirichlety", "dirichletz" };
      interleaved = flags.GetDefineFlag("interleaved");
      for (int i = 0; i <  ma->GetDimension(); i++)
        {
          Flags tmpflags = flags;
          if (flags.StringFlagDefined(dirnames[i]))
            {
              auto dir = flags.GetStringFlag(dirnames[i]);
              if(flags.StringFlagDefined("dirichlet"))
                dir += "|" + flags.GetStringFlag("dirichlet");
              // cout << "dirichlet = " << dir << endl;
              // cout << "dirichlet flag = " << flags.GetStringFlag("dirichlet") << endl;
              tmpflags.SetFlag ("dirichlet", dir);
            }
          if (flags.StringFlagDefined(dirnames[i]+"_bbnd"))
            {
              auto dir = flags.GetStringFlag(dirnames[i]+"_bbnd");
              if(flags.StringFlagDefined("dirichlet_bbnd"))
                dir += "|" + flags.GetStringFlag("dirichlet_bbnd");
              tmpflags.SetFlag ("dirichlet_bbnd", dir);
            }
          if (flags.StringFlagDefined(dirnames[i]+"_bbbnd"))
            {
              auto dir = flags.GetStringFlag(dirnames[i]+"_bbbnd");
              if(flags.StringFlagDefined("dirichlet_bbbnd"))
                dir += "|" + flags.GetStringFlag("dirichlet_bbbnd");
              tmpflags.SetFlag ("dirichlet_bbbnd", dir);
            }
          AddSpace (make_shared<H1HighOrderFESpace> (ama, tmpflags));
        }

      switch (ma->GetDimension())
        {
        case 2:
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdVectorH1<2>>>();
          flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradVectorH1<2>>>();
          evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdVectorH1<2,BND>>>();
          flux_evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpGradBoundaryVectorH1<2>>>();
          additional_evaluators.Set ("div", make_shared<T_DifferentialOperator<DiffOpDivVectorH1<2>>> ()); 
          additional_evaluators.Set ("divfree_reconstruction", make_shared<T_DifferentialOperator<DiffOpDivFreeReconstructVectorH1<2>>> ());
          additional_evaluators.Set ("Grad", make_shared<T_DifferentialOperator<DiffOpGradVectorH1<2>>> ());
          additional_evaluators.Set ("Gradboundary", make_shared<T_DifferentialOperator<DiffOpGradBoundaryVectorH1<2>>>());
          additional_evaluators.Set ("dual", make_shared<T_DifferentialOperator<DiffOpDualVectorH1<2,2>>> ());
          additional_evaluators.Set ("hesse", make_shared<VectorDifferentialOperator>(make_shared<T_DifferentialOperator<DiffOpHesse<2>>>(), 2));
          additional_evaluators.Set ("hesseboundary", make_shared<VectorDifferentialOperator>(make_shared<T_DifferentialOperator<DiffOpHesseBoundary<2>>>(), 2));
          break;
          
        case 3:
          evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdVectorH1<3>>>();
          flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradVectorH1<3>>>();
          evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdVectorH1<3,BND>>>();
          flux_evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpGradBoundaryVectorH1<3>>>();
          //flux_evaluator[BND] = make_shared<VectorDifferentialOperator>(make_shared<T_DifferentialOperator<DiffOpGradientBoundary<3>>>(), 3);
          evaluator[BBND] = make_shared<T_DifferentialOperator<DiffOpIdVectorH1<3,BBND>>>();
          
          additional_evaluators.Set ("div", make_shared<T_DifferentialOperator<DiffOpDivVectorH1<3>>> ());
          additional_evaluators.Set ("Grad", make_shared<T_DifferentialOperator<DiffOpGradVectorH1<3>>> ());
          additional_evaluators.Set ("Gradboundary", make_shared<T_DifferentialOperator<DiffOpGradBoundaryVectorH1<3>>>());
          additional_evaluators.Set ("dual", make_shared<T_DifferentialOperator<DiffOpDualVectorH1<3,3>>> ());
          additional_evaluators.Set ("hesse", make_shared<VectorDifferentialOperator>(make_shared<T_DifferentialOperator<DiffOpHesse<3>>>(), 3));
          additional_evaluators.Set ("hesseboundary", make_shared<VectorDifferentialOperator>(make_shared<T_DifferentialOperator<DiffOpHesseBoundary<3>>>(), 3));          
          break;
        }
    }

  FiniteElement & VectorH1FESpace :: GetFE (ElementId ei, Allocator & alloc) const
  {
    return *new (alloc) VectorFiniteElement (spaces[0]->GetFE(ei, alloc), spaces.Size());
  }

  DocInfo VectorH1FESpace :: GetDocu ()
  {
    DocInfo docu = CompoundFESpace::GetDocu();
    docu.Arg("interleaved") = "bool = False\n"
      "  ordering of dofs changed to x0, y0, z0, x1 ....";
    docu.Arg("dirichletx") = "regexpr\n"
      "  Regular expression string defining the dirichlet boundary\n"
      "  on the first component of VectorH1.\n"
      "  More than one boundary can be combined by the | operator,\n"
      "  i.e.: dirichletx = 'top|right'";
    docu.Arg("dirichlety") = "regexpr\n"
      "  Dirichlet boundary for the second component";
    docu.Arg("dirichletz") = "regexpr\n"
      "  Dirichlet boundary for the third component";
    docu.Arg("dirichletx_bbnd") = "regexpr\n"
      "  Regular expression string defining the dirichlet bboundary,\n"
      "  i.e. points in 2D and edges in 3D, on the first component.\n"
      "  More than one bboundary can be combined by the | operator,\n"
      "  i.e.: dirichletx_bbnd = 'top|right'";
    docu.Arg("dirichlety_bbnd") = "regexpr\n"
      "  Dirichlet bboundary for the second component";
    docu.Arg("dirichletz_bbnd") = "regexpr\n"
      "  Dirichlet bboundary for the third component";
    docu.Arg("dirichletx_bbbnd") = "regexpr\n"
      "  Regular expression string defining the dirichlet bbboundary,\n"
      "  i.e. points in 3D, on the first component.\n"
      "  More than one bbboundary can be combined by the | operator,\n"
      "  i.e.: dirichletx_bbbnd = 'top|right'";
    docu.Arg("dirichlety_bbbnd") = "regexpr\n"
      "  Dirichlet bbboundary for the second component";
    docu.Arg("dirichletz_bbbnd") = "regexpr\n"
      "  Dirichlet bbboundary for the third component";

    return docu;
  }

  void VectorH1FESpace :: FinalizeUpdate()
  {
    CompoundFESpace::FinalizeUpdate();

    if (interleaved)
      {
	free_dofs = make_shared<BitArray> (GetNDof());
	free_dofs->Set();

        size_t dim = spaces.Size();
	for (int i = 0; i < spaces.Size(); i++)
	  {
	    shared_ptr<BitArray> free_dofs_sub = spaces[i]->GetFreeDofs();
            size_t nd = free_dofs_sub->Size();
            for (size_t j = 0, jj = i; j < nd; j++, jj+=dim)
              if (!free_dofs_sub->Test(j))
                free_dofs->Clear (jj);
	  }

        for (size_t i = 0; i < ctofdof.Size(); i++)
          if (ctofdof[i] == UNUSED_DOF)
            free_dofs->Clear(i);

	dirichlet_dofs = *free_dofs;
	dirichlet_dofs.Invert();

        external_free_dofs = make_shared<BitArray>(GetNDof());
        *external_free_dofs = *free_dofs;
        for (int i = 0; i < ctofdof.Size(); i++)
          if (ctofdof[i] & CONDENSABLE_DOF)
            external_free_dofs->Clear(i);
      }
  }
    
  void VectorH1FESpace::SetOrder (ELEMENT_TYPE et, TORDER order)
  {
    FESpace::SetOrder(et, order);
    for (auto & spc : spaces)
      spc->SetOrder (et, order);
  }

  void VectorH1FESpace :: SetOrder (NodeId ni, int order) 
  {
    for (auto & spc : spaces)
      spc->SetOrder (ni, order);
  }
  
  int VectorH1FESpace :: GetOrder (NodeId ni) const
  {
    if (GetNSpaces())
      return spaces[0]->GetOrder (ni);
    return 0;
  }

  void VectorH1FESpace :: GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
    if (interleaved)
      {
        spaces[0]->GetDofNrs(ei, dnums);
        int size1 = dnums.Size();
        int dim = spaces.Size();
        dnums.SetSize(dim*size1);
        for (int i = size1-1; i >= 0; i--)
          {
            DofId di = dnums[i];
            for (int j = 0; j < dim; j++)
              dnums[dim*i+j] = dim*di+j;
          }
      }
    else
      CompoundFESpace::GetDofNrs(ei, dnums);
  }
  
  void VectorH1FESpace :: GetDofNrs (NodeId ni, Array<DofId> & dnums) const
  {
    if (interleaved)
      {
        spaces[0]->GetDofNrs(ni, dnums);
        int size1 = dnums.Size();
        int dim = spaces.Size();
        dnums.SetSize(dim*size1);
        for (int i = size1-1; i >= 0; i--)
          {
            DofId di = dnums[i];
            for (int j = 0; j < dim; j++)
              dnums[dim*i+j] = dim*di+j;
          }
      }
    else
      CompoundFESpace::GetDofNrs(ni, dnums);
  }

  
  
  static RegisterFESpace<H1HighOrderFESpace> init ("h1ho");
  static RegisterFESpace<VectorH1FESpace> initvec ("VectorH1");
}
 
