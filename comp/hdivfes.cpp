/*********************************************************************/
/* File:   hdivfes.cpp                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   27. Jan. 2003                                             */
/*********************************************************************/

/* 
   Finite Element Space for H(div)
*/

#include <comp.hpp>
#include <multigrid.hpp>

#include <../fem/hdivhofe.hpp>  
#include <../fem/hdivlofe.hpp>

namespace ngcomp
{
  using namespace ngcomp;

  /*
  RaviartThomasFESpace :: 
  RaviartThomasFESpace (shared_ptr<MeshAccess> ama,
			int adim, bool acomplex)
    
    : FESpace (ama, 1, adim, acomplex)
  {
    name="RaviartThomasFESpace(hdiv)";
    
    trig    = new FE_RTTrig0;
    segm    = new HDivNormalSegm0;

    if (ma.GetDimension() == 2)
      {
	Array<CoefficientFunction*> coeffs(1);
	coeffs[0] = new ConstantCoefficientFunction(1);
	evaluator = GetIntegrators().CreateBFI("masshdiv", 2, coeffs);
      }
  }
  */

  RaviartThomasFESpace :: RaviartThomasFESpace (shared_ptr<MeshAccess> ama, const Flags& flags, bool parseflags)
  : FESpace (ama, flags)
  {
    name="RaviartThomasFESpace(hdiv)";
    // defined flags
    DefineDefineFlag("hdiv");
    
    // parse flags
    if(parseflags) CheckFlags(flags);
    
    order = 1; // he: see above constructor!
        
    // trig    = new FE_RTTrig0;
    // segm    = new HDivNormalSegm0;

    // SetDummyFE<HDivDummyFE> ();
    
    if (ma->GetDimension() == 2)
    {
      Array<shared_ptr<CoefficientFunction>> coeffs(1);
      coeffs[0] = shared_ptr<CoefficientFunction> (new ConstantCoefficientFunction(1));
      integrator[VOL] = GetIntegrators().CreateBFI("masshdiv", 2, coeffs);
    }
    if (ma->GetDimension() == 3)
    {
      Array<shared_ptr<CoefficientFunction>> coeffs(1);
      coeffs[0] = shared_ptr<CoefficientFunction> (new ConstantCoefficientFunction(1));
      integrator[VOL] = GetIntegrators().CreateBFI("masshdiv", 3, coeffs);
    }
    if (ma->GetDimension() == 2)
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDiv<2>>>();
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdVecHDivBoundary<2>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDiv<2>>>();
      }
    else
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDiv<3>>>();
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdVecHDivBoundary<3>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDiv<3>>>();
      }
  }
      
    RaviartThomasFESpace :: ~RaviartThomasFESpace ()
  {
    ;
  }


  shared_ptr<FESpace> RaviartThomasFESpace :: 
  Create (shared_ptr<MeshAccess> ma, const Flags & flags)
  {
    int order = int(flags.GetNumFlag ("order", 0));

    if (order <= 0)
      return make_shared<RaviartThomasFESpace> (ma, flags, true);      
    else
      return make_shared<HDivHighOrderFESpace> (ma, flags, true);
  }

  
  void RaviartThomasFESpace :: Update()
  {
    shared_ptr<MeshAccess> ma = GetMeshAccess();
    int level = ma->GetNLevels();
    
    if (level == ndlevel.Size())
      return;
    
    if (ma->GetDimension() == 2)
      ndlevel.Append (ma->GetNEdges());
    else
      ndlevel.Append (ma->GetNFaces());

    // FinalizeUpdate (lh);
  }

  FiniteElement & RaviartThomasFESpace :: GetFE (ElementId ei, Allocator & lh) const
  {
    switch(ma->GetElType(ei))
      {
      case ET_TRIG: return *(new (lh) FE_RTTrig0);
      case ET_SEGM: return *(new (lh) HDivNormalSegm0);
      default: throw Exception ("Element type not available for RaviartThomasFESpace::GetFE");
      }
  }
  
  size_t RaviartThomasFESpace :: GetNDof () const throw()
  {
    return ndlevel.Last();
  }
  
  size_t RaviartThomasFESpace :: GetNDofLevel (int level) const
  {
    return ndlevel[level];
  }
  
  
  
  void RaviartThomasFESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    if(ei.VB()==VOL)
      {
	Array<int> forient(6);
	
	if (ma->GetDimension() == 2)
	  ma->GetElEdges (ei.Nr(), dnums, forient);
	else
	  ma->GetElFaces (ei.Nr(), dnums, forient);
	
	if (!DefinedOn (ei))
	  dnums = -1;
        return;
      }

    if(ei.VB()==BND)
      {
	if (ma->GetDimension() == 2)
	  {
	    int eoa[12];
	    Array<int> eorient(12, eoa);
	    ma->GetSElEdges (ei.Nr(), dnums, eorient);
	    
	    if (!DefinedOn(ei))
	      dnums = -1;
            
          }
    // (*testout) << "el = " << elnr << ", dofs = " << dnums << endl;
      }
  
    if(ei.VB()==BBND || ei.VB()==BBBND)
      {
        dnums.SetSize0();
        return;
      }

    /*
      int eoa[12];
      Array<int> eorient(12, eoa);
      GetMeshAccess().GetSElEdges (selnr, dnums, eorient);
      
      if (!DefinedOnBoundary (ma->GetSElIndex (selnr)))
      dnums = -1;
    */
    dnums.SetSize (1);
    dnums = -1;
    
    // (*testout) << "sel = " << selnr << ", dofs = " << dnums << endl;
  }

  /*
  Table<int> * RaviartThomasFESpace :: CreateSmoothingBlocks (int type) const
  {
    return 0;
  }
  */
  
  void RaviartThomasFESpace :: 
  VTransformMR (ElementId ei, 
		SliceMatrix<double> mat, TRANSFORM_TYPE tt) const
  {
    int nd = 3;
    bool boundary = (ei.VB() == BND);
    size_t elnr = ei.Nr();
    if (boundary) return;

    Vector<> fac(nd);
    
    GetTransformationFactors (elnr, fac);
    
    int i, j, k, l;
    
    if (tt & TRANSFORM_MAT_LEFT)
      for (k = 0; k < dimension; k++)
	for (i = 0; i < nd; i++)
	  for (j = 0; j < mat.Width(); j++)
	    mat(k+i*dimension, j) *= fac(i);
	
    if (tt & TRANSFORM_MAT_RIGHT)
      for (l = 0; l < dimension; l++)
	for (i = 0; i < mat.Height(); i++)
	  for (j = 0; j < nd; j++)
	    mat(i, l+j*dimension) *= fac(j);
  }
  
    
  void  RaviartThomasFESpace ::
  VTransformVR (ElementId ei, 
		SliceVector<double> vec, TRANSFORM_TYPE tt) const
  {
    int nd = 3;
    bool boundary = (ei.VB() == BND);
    size_t elnr = ei.Nr();
    
    if (boundary) 
      {
	ArrayMem<int,4> edge_nums, edge_orient;
	ma->GetSElEdges (elnr, edge_nums, edge_orient);
	vec *= edge_orient[0];
	return;
      }

    Vector<> fac(nd);
    
    GetTransformationFactors (elnr, fac);
    
    if ((tt & TRANSFORM_RHS) || (tt & TRANSFORM_SOL) || (tt & TRANSFORM_SOL_INVERSE))
      {
	for (int k = 0; k < dimension; k++)
	  for (int i = 0; i < nd; i++)
	    vec(k+i*dimension) *= fac(i);
      }
  }
  
  void RaviartThomasFESpace ::
  GetTransformationFactors (int elnr, FlatVector<> & fac) const
  {
    Array<int> edge_nums, edge_orient;
    
    fac = 1;

    ma->GetElEdges (elnr, edge_nums, edge_orient);
    for (int i = 0; i < 3; i++)
      fac(i) = edge_orient[i];
  }



  
  class BDM1Prolongation : public Prolongation
  {
    shared_ptr<MeshAccess> ma;
    const FESpace & space;
    
    array<Mat<3,3>, 20> boundaryprol; // 16 + 4
    // array<Mat<3,12>, 6+3+1> innerprol;
    array<Mat<3,12>, 120> innerprol;
    
  public:
    BDM1Prolongation(const FESpace & aspace)
      : ma(aspace.GetMeshAccess()), space(aspace)
    {
      // v0,v1 .. .split edge, v0 in new triangle
      // v2 ... third vertex of coasrse trig
      // v3 ... subdivision vertex
      
      ma->EnableTable("parentfaces");
      // boundary prol (only 8 cases are actually used)
      // add 4 more cases for tet red refinement
      for (int classnr = 0; classnr < 16; classnr++)
        {
          int verts[4] = { 1, 2, 3, 4 };
          if (classnr & 8) Swap(verts[0], verts[1]);
          if (classnr & 4) Swap(verts[1], verts[2]);
          if (classnr & 2) Swap(verts[0], verts[1]);
          if (classnr & 1) Swap(verts[2], verts[3]);

          int vertsc[3] = { verts[1], verts[2], verts[0] };
          int vertsf[3] = { verts[3], verts[2], verts[0] };
          // cout << "coarse: " << vertsc[0] << " " << vertsc[1] << " " << vertsc[2] << endl;
          // cout << "fine:   " << vertsf[0] << " " << vertsf[1] << " " << vertsf[2] << endl;
          
          HDivHighOrderNormalTrig<TrigExtensionMonomial> felc(1);
          felc.SetVertexNumbers (vertsc);

          HDivHighOrderNormalTrig<TrigExtensionMonomial> felf(1);
          felf.SetVertexNumbers (vertsf);
          
          IntegrationRule ir(ET_TRIG, 2);
          Matrix<> massf(3,3), massfc(3,3);
          Vector<> shapef(3), shapec(3);
          massf = 0; massfc = 0;
          
          for (IntegrationPoint ip : ir)
            {
              IntegrationPoint ipc(0.5*ip(0), ip(1));
              // MappedIntegrationPoint<2,2> mipc(ipc, trafoc);
              // MappedIntegrationPoint<2,2> mipf(ip, trafof);

              felc.CalcShape (ipc, shapec);
              felf.CalcShape (ip, shapef);

              massf += ip.Weight() * shapef * Trans(shapef);
              massfc += ip.Weight() * shapef * Trans(shapec);
            }
          CalcInverse (massf);
          boundaryprol[classnr] = 0.5 * massf * massfc;
          //cout << "boundarypol[" << classnr << "] = " << endl << FlatMatrix(boundaryprol[classnr]) << endl;
        }
      
      // 4 red refinement boundary faces
      for (int classnr =16; classnr < 20; classnr++)
        {
          Matrix<> bprol(3,3);
          if (classnr==16){
            bprol(0,0) = 0.25; bprol(0,1) = 0.0625; bprol(0,2) = -0.125;
            bprol(1,0) = 0; bprol(1,1) = 0.125; bprol(1,2) = 0;
            bprol(2,0) = 0; bprol(2,1) = 0.; bprol(2,2) = 0.125;
          }else if (classnr==17){
            bprol(0,0) = -0.25; bprol(0,1) = 0.0625; bprol(0,2) = 0.125;
            bprol(1,0) = 0; bprol(1,1) = 0.125; bprol(1,2) = 0;
            bprol(2,0) = 0; bprol(2,1) = 0.; bprol(2,2) = -0.125;
          }else if (classnr==18){
            bprol(0,0) = 0.25; bprol(0,1) = 0; bprol(0,2) = 0.25;
            bprol(1,0) = 0; bprol(1,1) = -0.0625; bprol(1,2) = 0.375;
            bprol(2,0) = 0; bprol(2,1) = -0.03125; bprol(2,2) = -0.0625;
          }else{
            bprol(0,0) = -0.25; bprol(0,1) = 0; bprol(0,2) = 0;
            bprol(1,0) = 0; bprol(1,1) = 0.0625; bprol(1,2) = 0.375;
            bprol(2,0) = 0; bprol(2,1) = 0.03125; bprol(2,2) = -0.0625;
          }
          boundaryprol[classnr] = bprol;
          //cout << "boundarypol[" << classnr << "] = " << endl << FlatMatrix(boundaryprol[classnr]) << endl;
        }


      /*
      // inner prol 
      for (int classnr = 0; classnr < 10; classnr++)
        {
          vector<int> verts={1,2,3,4,5};// v0, v1, v2, v3, vb
          if (classnr==0){
            verts={1,2,3,4,5};
          }else if (classnr==1){
            verts={1,3,2,4,5};
          }else if (classnr==2){
            verts={1,4,2,3,5};
          }else if (classnr==3){
            verts={2,3,1,4,5};
          }else if (classnr==4){
            verts={2,4,1,3,5};
          }else if (classnr==5){
            verts={3,4,1,2,5};
          }else if (classnr==6){
            verts={1,2,5,3,4};
          }else if (classnr==7){
            verts={1,3,5,2,4};
          }else if (classnr==8){
            verts={2,3,5,1,4};
          }else if (classnr==9){
            verts={1,2,4,5,3};
          }
          // ORIENTATION??
          // Coarse TET
          int vertsc[4] = { verts[0], verts[1], verts[2], verts[3] };
          // Fine TRIG
          int vertsf[3] = { verts[4], verts[2], verts[3] };
          cout << "coarse TET: " << vertsc[0] << " " << vertsc[1] << " " << vertsc[2] << 
            " " << vertsc[3] << endl;
          cout << "bisect TRIG:   " << vertsf[0] << " " << vertsf[1] << " " << vertsf[2] << endl;
          HDivHighOrderFE<ET_TET> felc(1) ;
          felc.SetVertexNumbers (vertsc);
          HDivHighOrderNormalTrig<TrigExtensionMonomial> felf(1);
          felf.SetVertexNumbers (vertsf);
          //cout << felc.GetNDof() << " " << felf.GetNDof() << endl;
          IntegrationRule ir(ET_TRIG, 2);
          Matrix<> massf(3,3), massfc(3,12);
          Matrix<> shapef(3,3), shapec(12,3);
          massf = 0; massfc = 0;
          
          // vertex coordinates of inner face within tet
          // points are images of ref-pnts  (1,0), (0,1), (0,0)
          Matrix<> points = { { 0.5, 0.5, 0 }, { 0, 0, 1 }, { 0, 0, 0 } };
          FE_ElementTransformation<2,3> trig2tet(ET_TRIG, points);
          
          for (IntegrationPoint ip : ir)
            {
              // map ip on trig back to tet
              MappedIntegrationPoint<2,3> mip(ip, trig2tet);
              IntegrationPoint iptet = mip.GetPoint();
              cout << "iptrig = " << ip << ", iptet = " << iptet << endl;
              felc.CalcShape (iptet, shapec);
              felf.CalcMappedShape (mip, shapef);

              massf += ip.Weight() * shapef * Trans(shapef);
              massfc += ip.Weight() * shapef * Trans(shapec);
            }
          CalcInverse (massf);
          Mat<3,12> prolmat = massf * massfc;
          // reorder, since tet-shapes are enumerated all lowest order first
          for (int i = 0; i < 4; i++)
            {
              innerprol[classnr].Col(3*i  ) = prolmat.Col(i);
              innerprol[classnr].Col(3*i+1) = prolmat.Col(4+2*i);
              innerprol[classnr].Col(3*i+2) = prolmat.Col(4+2*i+1);
            }
          cout << "innerpol[" << classnr << "] = " << endl << FlatMatrix(innerprol[classnr]) << endl;
        }
      */


      // inner prol 
      for (int classnr = 0; classnr < 120; classnr++)
        {
          // find a realiztion of permuation classnr

          int hv = classnr;
          int verts[5] = { -1, -1, -1, -1, -1 };
          for (int j = 0; j < 5; j++)
            {
              int nth_el = hv % (j+1);
              hv = hv/(j+1);

              for (int k = 5; k > nth_el; k--)
                verts[k] = verts[k-1];
              verts[nth_el] = j;
            }

          /*
          cout << "realization of class " << classnr << " is ";
          for (int k = 0; k < 5; k++)
            cout << verts[k] << " ";
          cout << endl;
          */
          
          // Coarse TET
          int vertsc[4] = { verts[0], verts[1], verts[2], verts[3] };
          // Fine TRIG
          int vertsf[3] = { verts[4], verts[2], verts[3] };

          /*
          cout << "coarse TET: " << vertsc[0] << " " << vertsc[1] << " " << vertsc[2] << 
            " " << vertsc[3] << endl;
          cout << "bisect TRIG:   " << vertsf[0] << " " << vertsf[1] << " " << vertsf[2] << endl;
          */
          
          HDivHighOrderFE<ET_TET> felc(1) ;
          felc.SetVertexNumbers (vertsc);
          HDivHighOrderNormalTrig<TrigExtensionMonomial> felf(1);
          felf.SetVertexNumbers (vertsf);
          //cout << felc.GetNDof() << " " << felf.GetNDof() << endl;
          IntegrationRule ir(ET_TRIG, 2);
          Matrix<> massf(3,3), massfc(3,12);
          Matrix<> shapef(3,3), shapec(12,3);
          massf = 0; massfc = 0;
          
          // vertex coordinates of inner face within tet
          // points are images of ref-pnts  (1,0), (0,1), (0,0)
          Matrix<> points = { { 0.5, 0.5, 0 }, { 0, 0, 1 }, { 0, 0, 0 } };
          FE_ElementTransformation<2,3> trig2tet(ET_TRIG, points);
          
          for (IntegrationPoint ip : ir)
            {
              // map ip on trig back to tet
              MappedIntegrationPoint<2,3> mip(ip, trig2tet);
              IntegrationPoint iptet = mip.GetPoint();
              // cout << "iptrig = " << ip << ", iptet = " << iptet << endl;
              felc.CalcShape (iptet, shapec);
              felf.CalcMappedShape (mip, shapef);

              massf += ip.Weight() * shapef * Trans(shapef);
              massfc += ip.Weight() * shapef * Trans(shapec);
            }
          CalcInverse (massf);
          Mat<3,12> prolmat = massf * massfc;
          // reorder, since tet-shapes are enumerated all lowest order first
          for (int i = 0; i < 4; i++)
            {
              innerprol[classnr].Col(3*i  ) = prolmat.Col(i);
              innerprol[classnr].Col(3*i+1) = prolmat.Col(4+2*i);
              innerprol[classnr].Col(3*i+2) = prolmat.Col(4+2*i+1);
            }
          // innerprol[classnr] = 0.0;
          // cout << "innerpol[" << classnr << "] = " << endl << FlatMatrix(innerprol[classnr]) << endl;
        }
    }

    
    virtual ~BDM1Prolongation() { }
  
    virtual void Update (const FESpace & fes) { ; }
    virtual shared_ptr<SparseMatrix< double >> CreateProlongationMatrix( int finelevel ) const
    { return nullptr; }

    virtual shared_ptr<BitArray> GetInnerDofs (int finelevel) const
    {
      size_t nc = space.GetNDofLevel (finelevel-1) / 3;
      size_t nf = space.GetNDofLevel (finelevel) / 3;

      BitArray inner(3*nf);
      inner.Clear();
      auto freedofs = space.GetFreeDofs(true);
      
      for (size_t i = nc; i < nf; i++)
        {
          auto [info, nrs] = ma->GetParentFaces(i);
          if (nrs[1] != -1 || info==20)
            for (int j = 0; j < 3; j++)
              if (freedofs->Test(3*i+j))
                inner.SetBit(3*i+j);
        }
      
      cout << IM(5) << "prolongation level " << finelevel
           << " #innerdofs: " << endl << inner.NumSet() << "/" << inner.Size() << endl;
        
      return make_shared<BitArray> (inner);
    }

    virtual void ProlongateInline (int finelevel, BaseVector & v) const
    {
      // cout << "prolongate Hdiv" << endl;
      size_t nc = space.GetNDofLevel (finelevel-1) / 3;
      size_t nf = space.GetNDofLevel (finelevel) / 3;
      
      auto fv = v.FV<double>();
      fv.Range(3*nc, fv.Size()) = 0;

      for (int loop = 0; loop < 5; loop++)
      for (size_t i = nc; i < nf; i++)
        {
          auto [info, nrs] = ma->GetParentFaces(i);
	  int pa1 = nrs[0];
	  int pa2 = nrs[1];
	  int pa3 = nrs[2];
          int pa4 = nrs[3];

          if (pa2 == -1)
            {
              if (info==20) ;// interior face:: do nothing
              else{ // bisect or red face
                Vec<3> fvecc = fv.Range(3*pa1, 3*pa1+3);
                Vec<3> fvecf = boundaryprol[info] * fvecc;
                fv.Range(3*i, 3*i+3) = fvecf;
              }
            }
          else
            {
              Vec<12> fvecc;
              fvecc.Range(0,3) = fv.Range(3*pa1, 3*pa1+3);
              fvecc.Range(3,6) = fv.Range(3*pa2, 3*pa2+3);
              fvecc.Range(6,9) = fv.Range(3*pa3, 3*pa3+3);
              fvecc.Range(9,12) = fv.Range(3*pa4, 3*pa4+3);
              Vec<3> fvecf = innerprol[info] * fvecc;
              fv.Range(3*i, 3*i+3) = fvecf;
            }
        }

      // every face from coarse level got split
      for (size_t i = 0; i < nf; i++)
        {
          auto [info, nrs] = ma->GetParentFaces(i);
          if (nrs[0] != -1 && nrs[1] == -1)
            {
              fv(3*nrs[0]) = 0;
              fv(3*nrs[0]+1) = 0;
              fv(3*nrs[0]+2) = 0;
            }
        }
    }
    
    virtual void RestrictInline (int finelevel, BaseVector & v) const
    {
      size_t nc = space.GetNDofLevel (finelevel-1) / 3;
      size_t nf = space.GetNDofLevel (finelevel) / 3;
      
      auto fv = v.FV<double>();
      fv.Range(3*nf, fv.Size()) = 0;

      // every face from coarse level got split
      for (size_t i = 0; i < nf; i++)
        {
          auto [info, nrs] = ma->GetParentFaces(i);
          if (nrs[0] != -1 && nrs[1] == -1)
            {
              fv(3*nrs[0]) = 0;
              fv(3*nrs[0]+1) = 0;
              fv(3*nrs[0]+2) = 0;
            }
        }

      for (int loop = 0; loop < 5; loop++)      
      for (size_t i = nf; i-- > nc; )
        {
          auto [info, nrs] = ma->GetParentFaces(i);
	  int pa1 = nrs[0];
	  int pa2 = nrs[1];
	  int pa3 = nrs[2];
	  int pa4 = nrs[3];

          if (pa2 == -1)
            {
              if (info==20) ;// interior face:: do nothing
              else{ // bisect or red face
                Vec<3> fvecf = fv.Range(3*i, 3*i+3);
                Vec<3> fvecc = Trans(boundaryprol[info]) * fvecf;
                fv.Range(3*pa1, 3*pa1+3) += fvecc;
                fv.Range(3*i, 3*i+3) = 0.0;           
              }
            }
          else
            {
              Vec<3> fvecf = fv.Range(3*i, 3*i+3);
              Vec<12> fvecc = Trans(innerprol[info]) * fvecf;
              
              fv.Range(3*pa1, 3*pa1+3) += fvecc.Range(0,3);
              fv.Range(3*pa2, 3*pa2+3) += fvecc.Range(3,6);
              fv.Range(3*pa3, 3*pa3+3) += fvecc.Range(6,9);
              fv.Range(3*pa4, 3*pa4+3) += fvecc.Range(9,12);
              
              fv.Range(3*i, 3*i+3) = 0.0;
            }
        }
    }
  };
  
// BDM1 prol in TRIG
  class BDM1ProlongationTRIG : public Prolongation
  {
    shared_ptr<MeshAccess> ma;
    const FESpace & space;
  public:
    BDM1ProlongationTRIG(const FESpace & aspace)
      : ma(aspace.GetMeshAccess()), space(aspace)
    {
      ma->EnableTable("parentedges");
    }
    
    virtual ~BDM1ProlongationTRIG() { }
    
    virtual void Update (const FESpace & fes) { ; }
    virtual shared_ptr<SparseMatrix< double >> CreateProlongationMatrix( int finelevel ) const
    { return nullptr; }
    
    virtual void ProlongateInline (int finelevel, BaseVector & v) const
    {
      size_t nc = space.GetNDofLevel (finelevel-1) / 2;
      size_t nf = space.GetNDofLevel (finelevel) / 2;
      
      auto fv = v.FV<double>();
      fv.Range(2*nf, fv.Size()) = 0;
      
      for (size_t i = nc; i < nf; i++)
        {
          auto [info, nrs] = ma->GetParentEdges(i);
          int pa1 = nrs[0];
          int pa2 = nrs[1];
          int pa3 = nrs[2];
          
          if (pa2 == -1)
            {
              double fac0 = (info & 1) ? 0.5 : -0.5;
              fv(2*i)   = fac0 * fv(2*pa1) + 0.125 * fv(2*pa1+1);
              fv(2*i+1) = 0.25 * fv(2*pa1+1);
            }
          else if (info<8)//bisecting edge
            {
              double fac1 = (info&1) ? 0.5 : -0.5;
              double fac2 = (info&2) ? 0.5 : -0.5;
              double fac3 = (info&4) ? 0.125 : -0.125;
              fv(2*i) = fac1 * fv(2*pa1) + fac2 * fv(2*pa2) + fac3 * fv(2*pa3+1);
              fv(2*i+1) = 0.5 * (fv(2*pa1+1)+fv(2*pa2+1)) - 0.25*fv(2*pa3+1);
            }
          else // info>=8: red edge
            {
              double fac1 = (info&1) ? 0.25 : -0.25;
              double fac2 = (info&2) ? 0.25 : -0.25;
              double fac3 = (info&4) ? 0.25 : -0.25;
              fv(2*i) = fac1 * fv(2*pa1) + fac2 * fv(2*pa2) + fac3 * fv(2*pa3)
                - 0.125 * fv(2*pa1+1) + 0.125 * fv(2*pa2+1);
              fv(2*i+1) = 0.25*fv(2*pa3+1);
            }
        }
      
      // every edge from coarse level got split
      for (size_t i = 0; i < nf; i++)
        {
          auto [info, nrs] = ma->GetParentEdges(i);
          if (nrs[0] != -1 && nrs[1] == -1)
            {
              fv(2*nrs[0]) = 0;
              fv(2*nrs[0]+1) = 0;
            }
        }
    }
    
    virtual void RestrictInline (int finelevel, BaseVector & v) const
    {
      size_t nc = space.GetNDofLevel (finelevel-1) / 2;
      size_t nf = space.GetNDofLevel (finelevel) / 2;
      
      auto fv = v.FV<double>();
      fv.Range(2*nf, fv.Size()) = 0;
      
      // every edge from coarse level got split
      for (size_t i = 0; i < nf; i++)
        {
          auto [info, nrs] = ma->GetParentEdges(i);
          if (nrs[0] != -1 && nrs[1] == -1)
            {
              fv(2*nrs[0]) = 0;
              fv(2*nrs[0]+1) = 0;
            }
        }
      
      
      for (size_t i = nf; i-- > nc; )
        {
          auto [info, nrs] = ma->GetParentEdges(i);
          int pa1 = nrs[0];
          int pa2 = nrs[1];
          int pa3 = nrs[2];
          
          if (pa2 == -1)
            {
              double fac0 = (info & 1) ? 0.5 : -0.5;
              fv(2*pa1) += fac0 * fv(2*i);
              fv(2*pa1+1) += 0.125 * fv(2*i) + 0.25 * fv(2*i+1);
            }
          else if (info<8)//bisecting edge
            {
              double fac1 = (info&1) ? 0.5 : -0.5;
              double fac2 = (info&2) ? 0.5 : -0.5;
              double fac3 = (info&4) ? 0.125 : -0.125;
              fv(2*pa1)   += fac1 * fv(2*i);
              fv(2*pa1+1) += 0.5 * fv(2*i+1);
              fv(2*pa2)   += fac2 * fv(2*i);
              fv(2*pa2+1) += 0.5 * fv(2*i+1);
              fv(2*pa3+1) += fac3 * fv(2*i) - 0.25 * fv(2*i+1);
            }
          else // info>=8: red edge
            {
              double fac1 = (info&1) ? 0.25 : -0.25;
              double fac2 = (info&2) ? 0.25 : -0.25;
              double fac3 = (info&4) ? 0.25 : -0.25;
              fv(2*pa1)   += fac1 * fv(2*i);
              fv(2*pa1+1) -= 0.125 * fv(2*i);
              fv(2*pa2)   += fac2 * fv(2*i);
              fv(2*pa2+1) += 0.125 * fv(2*i);
              fv(2*pa3) += fac3 * fv(2*i);
              fv(2*pa3+1) += 0.25*fv(2*i+1);
            }
        }
      
    }
  };
  
  
  
  class NGS_DLL_HEADER BDM1FESpace : public FESpace
  {
    BitArray active_facets;
    
  public:
    BDM1FESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags=false)
      : FESpace(ama, flags)
      {
        name="BDM1FESpace";
        
        if (ma->GetDimension() == 2)
          {
            evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDiv<2>>>();
            evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdVecHDivBoundary<2>>>();
            flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDiv<2>>>();
            additional_evaluators.Set ("grad", make_shared<T_DifferentialOperator<DiffOpGradientHDiv<2>>> ());
	    additional_evaluators.Set ("dual", make_shared<T_DifferentialOperator<DiffOpHDivDual<2>>> ());
            prol = make_shared<BDM1ProlongationTRIG> (*this);
          }
        else if(ma->GetDimension() == 3) 
          {
            evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDiv<3>>>();
            evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdVecHDivBoundary<3>>>();
            flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDiv<3>>>();
            additional_evaluators.Set ("grad", make_shared<T_DifferentialOperator<DiffOpGradientHDiv<3>>> ());
	    additional_evaluators.Set ("dual", make_shared<T_DifferentialOperator<DiffOpHDivDual<3>>> ());
	    prol = make_shared<BDM1Prolongation> (*this);
          }
      }
    
    virtual ~BDM1FESpace () { }
    virtual const char * GetType()  { return "BDM1"; }
    
    static shared_ptr<FESpace> Create (shared_ptr<MeshAccess> ma, const Flags & flags)
    {
      return make_shared<BDM1FESpace> (ma, flags, true);
    }
    
    void Update() override
    {
      if (ma->GetDimension() == 3)
        {
          size_t nfa = ma->GetNFaces();
          SetNDof (3*nfa);
          active_facets = BitArray(nfa);
          active_facets.Clear();
          for (auto el : ma->Elements(VOL))
            for (auto fa : el.Faces())
              active_facets.SetBit(fa);
      
          ctofdof.SetSize(GetNDof());
          ctofdof = WIREBASKET_DOF;
          for (size_t i = 0; i < nfa; i++)
            if (!active_facets.Test(i))
              ctofdof[3*i] = ctofdof[3*i+1] = ctofdof[3*i+2] = UNUSED_DOF;
          
          // cout << "active faces = " << endl << active_facets << endl;
        }
      else if (ma->GetDimension()==2)
      {
        size_t ned = ma->GetNEdges();
        SetNDof (2*ned);
        active_facets = BitArray(ned);
        active_facets.Clear();
        for (auto el : ma->Elements(VOL))
          for (auto fa : el.Edges())
            active_facets.SetBit(fa);

        ctofdof.SetSize(GetNDof());
        ctofdof = WIREBASKET_DOF;
        for (size_t i = 0; i < ned; i++)
          if (!active_facets.Test(i))
            ctofdof[2*i] = ctofdof[2*i+1] = UNUSED_DOF;
      }
    }
    
    // virtual void DoArchive (Archive & archive) override;
    // virtual void UpdateCouplingDofArray() override;
    
    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const override
    {
      if (ei.VB() == VOL)
        switch (ma->GetElType(ei))
          {
          case ET_TET:   
            {
              Ngs_Element ngel = ma->GetElement<3,VOL> (ei.Nr());
              if (!DefinedOn(ngel)) return * new (lh) HDivDummyFE<ET_TET>();
    
              auto * fe =  new (lh) HDivHighOrderFE<ET_TET> (1);
              fe -> SetVertexNumbers (ngel.Vertices());
              return *fe;
            }

          case ET_TRIG: 
            {
              Ngs_Element ngel = ma->GetElement<2,VOL> (ei.Nr());
              if (!DefinedOn(ngel)) return * new (lh) HDivDummyFE<ET_TRIG>();
              
              auto * fe =  new (lh) HDivHighOrderFE<ET_TRIG> (1);
              fe -> SetVertexNumbers (ngel.Vertices());
              return *fe;
            }
          default:
            ;
          }
      
      else if (ei.VB() == BND)

        switch (ma->GetElType(ei))
          {
          case ET_TRIG:   
            {
              Ngs_Element ngel = ma->GetElement<2,BND> (ei.Nr());
              if (!DefinedOn(ngel)) return * new (lh) HDivNormalDummyFE<ET_TRIG>();
    
              auto fe = new (lh) HDivHighOrderNormalTrig<TrigExtensionMonomial> (1);
              fe -> SetVertexNumbers (ngel.Vertices());
              return *fe;
            }
          case ET_SEGM:
            {
              Ngs_Element ngel = ma->GetElement<1,BND> (ei.Nr());
              if (!DefinedOn(ngel)) return * new (lh) HDivNormalDummyFE<ET_SEGM>();
    
              auto fe = new (lh) HDivHighOrderNormalSegm<TrigExtensionMonomial> (1);
              fe -> SetVertexNumbers (ngel.Vertices());
              return *fe;
            }
          default:
            ;
          }
      
      throw Exception ("Element not available in BDM1 space");
    }
    
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override
    {
      if (ma->GetDimension() == 3)
        {
          auto faces = ma->GetElFaces (ei);
          dnums.SetSize(3*faces.Size());
          for (int i : Range(faces))
            {
              dnums[i] = 3*faces[i];
              dnums[2*i+faces.Size()] = 3*faces[i]+1;
              dnums[2*i+1+faces.Size()] = 3*faces[i]+2;
            }
        }
      else if (ma->GetDimension()==2)
      {
          auto edges = ma->GetElEdges (ei);
          dnums.SetSize(2*edges.Size());
          for (int i : Range(edges))
            {
              dnums[i] = 2*edges[i];
              dnums[i+edges.Size()] = 2*edges[i]+1;
            }

      }
      // cout << "Ei = " << ei << "dnums = " << dnums << endl;
    }

    
    virtual string GetClassName () const override
    {
      return "BDM11FESpace";
    }
    
    virtual void GetVertexDofNrs (int vnr, Array<DofId> & dnums) const override
    { dnums.SetSize0(); }
    virtual void GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const override
    {
      dnums.SetSize0();
      // for 2D ...
      if (ma->GetDimension()==2)
        if (active_facets.Test(ednr))
          {
            dnums.SetSize(2);
            dnums[0] = 2*ednr;
            dnums[1] = 2*ednr+1;
          }
    }
    
    virtual void GetFaceDofNrs (int fanr, Array<DofId> & dnums) const override
    {
      dnums.SetSize0();
      if (ma->GetDimension()==3)      
        if (active_facets.Test(fanr))
          {
            dnums.SetSize(3);
            dnums[0] = 3*fanr;
            dnums[1] = 3*fanr+1;
            dnums[2] = 3*fanr+2;
          }
    }    
    virtual void GetInnerDofNrs (int elnr, Array<DofId> & dnums) const override
    { dnums.SetSize0(); }    
  };

  static RegisterFESpace<BDM1FESpace> initbdm1 ("BDM1");




  

  // register FESpaces
  namespace hdivfes_cpp
  {
    class Init
    { 
    public: 
      Init ();
    };
    
    Init::Init()
    {
      GetFESpaceClasses().AddFESpace ("hdiv", RaviartThomasFESpace::Create,
                                      RaviartThomasFESpace::GetDocu);
    }
    
    Init init;
  }
}
