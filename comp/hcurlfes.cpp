/*********************************************************************/
/* File:   hcurlhdivfes.cpp                                          */
/* Author: Joachim Schoeberl                                         */
/* Date:   12. Jan. 2002                                             */
/*********************************************************************/

/* 
   Finite Element Space
*/

#include "hcurlfes.hpp"
#include <prolongation.hpp>

#include <../fem/hcurllofe.hpp>
#include <../fem/thcurlfe_impl.hpp>
#include <../fem/hcurlfe_utils.hpp>
#include <../fem/hcurl_equations.hpp> 
#include <diffop_impl.hpp>

using namespace ngmg; 


namespace ngcomp
{

  // Nedelec FE Space
  NedelecFESpace :: NedelecFESpace (shared_ptr<MeshAccess> ama, const Flags& flags, bool parseflags)
    : FESpace (ama, flags)
  {
    name="NedelecFESpace(hcurl)";
    DefineDefineFlag("hcurl");
    // parse standard flags
    if(parseflags) CheckFlags(flags);
    
    if( flags.GetDefineFlag("hcurl"))
      cerr << "WARNING: -hcurl flag is deprecated: use -type=hcurl instead" << endl;
    
    // SetDummyFE<HCurlDummyFE> ();

    prol = make_shared<EdgeProlongation> (*this);
    order = 1;

    // auto one = make_shared<ConstantCoefficientFunction>(1);
    // integrator[VOL] = GetIntegrators().CreateBFI("massedge", ma->GetDimension(), one);
    // integrator[BND] = GetIntegrators().CreateBFI("robinedge", ma->GetDimension(), one);

    if (ma->GetDimension() == 2)
      {
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryEdge<2>>>();
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdEdge<2>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpCurlEdge<2>>>();        
      }
    else if(ma->GetDimension() == 3) 
      {
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryEdge<3>>>();
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdEdge<3>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpCurlEdge<3>>>();
        flux_evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpCurlBoundaryEdgeVec<>>>();
	evaluator[BBND] = make_shared<T_DifferentialOperator<DiffOpIdBBoundaryEdge<3>>>();	
      }
    
    discontinuous = flags.GetDefineFlag("discontinuous");
  }

                                    
  NedelecFESpace :: ~NedelecFESpace ()
  {
    ;
  }


  shared_ptr<FESpace> NedelecFESpace :: Create (shared_ptr<MeshAccess> ma, const Flags & flags)
  {
    return make_shared<NedelecFESpace> (ma, flags, true);
  }



  
  void NedelecFESpace :: Update()
  {
    // size_t ne = ma->GetNE();
    // size_t nse = ma->GetNSE();
    size_t ned = ma->GetNEdges();
    
    int level = ma->GetNLevels();
    
    if (level == nelevel.Size())
      return;
    
    nelevel.Append (ned);

    /*
    for (int i=0; i<specialelements.Size(); i++)
      delete specialelements[i];
    specialelements.DeleteAll();
    */
    
    // new implementation of finelevelofedge - array:
    
    size_t oldned = finelevelofedge.Size();
    finelevelofedge.SetSize(ned);
    finelevelofedge.Range (oldned, ned) = -1;

    /*
    for (Ngs_Element el : ma->Elements(VOL))
      if (DefinedOn (el)) 
        finelevelofedge[el.Edges()] = level-1;
    
    for (Ngs_Element el : ma->Elements(BND))
      if (DefinedOn (el)) 
	finelevelofedge[el.Edges()] = level-1;
    */
    for (VorB vb : { VOL, BND })
      ma -> IterateElements
        (vb, [&] (auto el)
         {
           if (this->DefinedOn (el)) 
             finelevelofedge[el.Edges()] = level-1;         
         });


    

    if (ma->HasParentEdges())
      {
        parentedges.SetSize(ned);
        for (int i = 0; i < ned; i++)
          {
            auto [info,nrs] = ma->GetParentEdges(i);
            // cout << "parent of " << i << " : info = " << info
            // << " nrs = " << nrs[0] << "," << nrs[1] << "," << nrs[2] << endl;
            if (nrs[0] > i || nrs[1] > i) cout << "parent is larger" << endl;
            parentedges[i][0] = (nrs[0]==-1) ? -1 : 2*nrs[0]+(info&1);
            parentedges[i][1] = (nrs[1]==-1) ? -1 : 2*nrs[1]+ ( (info&2) / 2);
          }
      }

    // cout << "parentedges, ng: " << parentedges << endl;



    
    // generate edge points, and temporary hash table
    ClosedHashTable<IVec<2>, int> node2edge(5*ned+10);

    edgepoints.SetSize0();
    
    for (size_t i = 0; i < ned; i++)
      {
	IVec<2> edge = ma->GetEdgePNums (i);
	int edgedir = (edge[0] > edge[1]);
	if (edgedir) Swap (edge[0], edge[1]);
	node2edge.Set (edge, i);
	edgepoints.Append (edge);
      }


    // if (!ma->HasParentEdges())
      {
		    static Timer t("build_hierarchy"); RegionTimer reg(t);
    // build edge hierarchy:
    parentedges.SetSize (ned);
    parentedges = IVec<2> (-1,-1);

    for (size_t i = 0; i < ned; i++)
      {
	// cout << "edge " << i << "/" << ned << endl;
	IVec<2> i2 (edgepoints[i][0], edgepoints[i][1]);
	int pa1[2], pa2[2];
	ma->GetParentNodes (i2[0], pa1);
	ma->GetParentNodes (i2[1], pa2);
	
	if (pa1[0] == -1 && pa2[0] == -1)
	  continue;
	
	int issplitedge = 0;
	if (pa1[0] == i2[1] || pa1[1] == i2[1])
	  issplitedge = 1;
	if (pa2[0] == i2[0] || pa2[1] == i2[0])
	  issplitedge = 2;
	
	if (issplitedge)
	  {
	    // edge is obtained by splitting one edge into two parts:
	    IVec<2> paedge;
	    if (issplitedge == 1)
	      paedge = IVec<2> (pa1[0], pa1[1]);
	    else
	      paedge = IVec<2> (pa2[0], pa2[1]);
	    
	    if (paedge[0] > paedge[1]) 
	      Swap (paedge[0], paedge[1]);
	    
	    int paedgenr = node2edge.Get (paedge);
	    int orient = (paedge[0] == i2[0] || paedge[1] == i2[1]) ? 1 : 0;
	    
	    parentedges[i][0] = 2 * paedgenr + orient;
	  }
	else
	  {
	    // edge is splitting edge in middle of triangle:
	    for (int j = 1; j <= 2; j++)
	      {
		IVec<2> paedge1, paedge2;
		if (j == 1)
		  {
		    paedge1 = IVec<2> (pa1[0], i2[1]);
		    paedge2 = IVec<2> (pa1[1], i2[1]);
		  }
		else
		  {
		    paedge1 = IVec<2> (pa2[0], i2[0]);
		    paedge2 = IVec<2> (pa2[1], i2[0]);
		  }
		if (paedge1[0] > paedge1[1]) 
		  Swap (paedge1[0], paedge1[1]);
		if (paedge2[0] > paedge2[1]) 
		  Swap (paedge2[0], paedge2[1]);
		
		int paedgenr1 = 0, paedgenr2 = 0;
		int orient1, orient2;
		
		// if first vertex number is -1, then don't try to find entry in node2edge hash table
		if ( paedge1[0] == -1 || paedge2[0] == -1 )
		  continue;

		if (node2edge.Used (paedge1) && node2edge.Used (paedge2))
		  {
		    // cout << "paedge1 = " << paedge1 << ", i2 = " << i2 << endl;		    
		    // cout << "paedge2 = " << paedge2 << ", i2 = " << i2 << endl;		    
		    paedgenr1 = node2edge.Get (paedge1);
		    orient1 = (paedge1[0] == i2[0] || paedge1[1] == i2[1]) ? 1 : 0;
		    paedgenr2 = node2edge.Get (paedge2);
		    orient2 = (paedge2[0] == i2[0] || paedge2[1] == i2[1]) ? 1 : 0;
		    // cout << "orient1 = " << orient1 << endl;
		    // cout << "orient2 = " << orient2 << endl;		    		    
		    parentedges[i][0] = 2 * paedgenr1 + orient1;	      
		    parentedges[i][1] = 2 * paedgenr2 + orient2;	      
		  }
	      }
	    
	    if (parentedges[i][0] == -1)
	      {
		// quad split
		if (pa1[0] != pa2[0] && 
		    pa1[0] != pa2[1] && 
		    pa1[1] != pa2[0] && 
		    pa1[1] != pa2[1])
		  for (int j = 1; j <= 2; j++)
		    {
		      IVec<2> paedge1, paedge2;
		      if (j == 1)
			{
			  paedge1 = IVec<2> (pa1[0], pa2[0]);
			  paedge2 = IVec<2> (pa1[1], pa2[1]);
			}
		      else
			{
			  paedge1 = IVec<2> (pa1[0], pa2[1]);
			  paedge2 = IVec<2> (pa1[1], pa2[0]);
			}
		      
		      int paedgenr1 = 0, paedgenr2 = 0;
		      int orient1 = 1, orient2 = 1;
		      
		      if (paedge1[0] > paedge1[1]) 
			{
			  Swap (paedge1[0], paedge1[1]);
			  orient1 = 0;
			}
		      if (paedge2[0] > paedge2[1]) 
			{
			  Swap (paedge2[0], paedge2[1]);
			  orient2 = 0;
			}

		      if ( paedge1[0] == -1 || paedge2[0] == -1 )
			continue;
		      
		      if (node2edge.Used (paedge1) && node2edge.Used (paedge2))
			{
			  paedgenr1 = node2edge.Get (paedge1);
			  paedgenr2 = node2edge.Get (paedge2);
			  parentedges[i][0] = 2 * paedgenr1 + orient1;	      
			  parentedges[i][1] = 2 * paedgenr2 + orient2;	      
			}
		    }
	      }
	    
	    if (parentedges[i][0] == -1)
	      {
		// triangle split into quad+trig (from anisotropic pyramids)
		for (int j = 0; j < 2; j++)
		  for (int k = 0; k < 2; k++)
		    {
		      IVec<2> paedge (pa1[1-j], pa2[1-k]);
		      int orientpa = 1;
		      if (paedge[0] > paedge[1]) 
			{
			  Swap (paedge[0], paedge[1]);
			  orientpa = 0;
			}	
		      if (pa1[j] == pa2[k] && node2edge.Used(paedge))
			{
			  int paedgenr = node2edge.Get (paedge);
			  parentedges[i][0] = 2 * paedgenr + orientpa;
			}
		    }
	      }
	  }
      
	if (i > nelevel[0] && parentedges[i][0] == -1)
	  {
	    cerr << "no parent edge found, edge = " 
		 << i2[0] << ", " << i2[1] 
		 << ", pa1 = " << pa1[0] << ", " << pa1[1] 
		 << ", pa2 = " << pa2[0] << ", " << pa2[1]
		 << endl;
	  }
      }
      }

    // cout << "parentedges = " << endl << parentedges << endl;

    prol->Update(*this);
    UpdateCouplingDofArray();
  }


  void NedelecFESpace :: DoArchive(Archive & archive)
  {
    FESpace::DoArchive(archive);
    archive & edgepoints & parentedges;
    archive & finelevelofedge & nelevel;
    archive & discontinuous;
  }

  void  NedelecFESpace :: UpdateCouplingDofArray ()
  {
    int level = ma->GetNLevels()-1;

    ctofdof.SetSize(GetNDof());
    
    for (int edge = 0; edge < ma->GetNEdges(); edge++) 
      ctofdof[edge] = 
	(FineLevelOfEdge(edge) == level) ? WIREBASKET_DOF : UNUSED_DOF; 
  }


  FiniteElement & NedelecFESpace :: GetFE (ElementId ei, Allocator & lh) const
  {
    if(!DefinedOn(ei))
      {
        switch(ma->GetElType(ei))
          {
          case ET_TET:     return * new (lh) HCurlDummyFE<ET_TET>();
          case ET_PRISM:   return * new (lh) HCurlDummyFE<ET_PRISM>();
          case ET_PYRAMID: return * new (lh) HCurlDummyFE<ET_PYRAMID>();
          case ET_TRIG:    return * new (lh) HCurlDummyFE<ET_TRIG>();
          case ET_QUAD:    return * new (lh) HCurlDummyFE<ET_QUAD>();
          case ET_SEGM:    return * new (lh) HCurlDummyFE<ET_SEGM>();
          case ET_HEX:     return * new (lh) HCurlDummyFE<ET_HEX>();
          default:
            throw Exception ("Inconsistent element type in NedelecFESpace::GetFE");
          }
      }
    switch (ma->GetElType(ei))
      {
      case ET_TET:     return * new (lh) FE_NedelecTet1;
      case ET_PRISM:   return * new (lh) FE_NedelecPrism1;
      case ET_PYRAMID: return * new (lh) FE_NedelecPyramid1;
      case ET_TRIG:    return * new (lh) FE_NedelecTrig1;
      case ET_QUAD:    return * new (lh) FE_NedelecQuad1;
      case ET_SEGM:    return * new (lh) FE_NedelecSegm1;
      case ET_HEX:     return * new (lh) FE_NedelecHex1;
      default:
        throw Exception ("Inconsistent element type in NedelecFESpace::GetFE");
      }
  }

  size_t NedelecFESpace :: GetNDof () const throw()
  {
    return nelevel.Last();
  }

  size_t NedelecFESpace :: GetNDofLevel (int level) const
  {
    return nelevel[level];
  }


  void NedelecFESpace :: GetDofRanges (ElementId ei, Array<IntRange> & dranges) const
  {
    dranges.SetSize (0);
    if (!DefinedOn (ei)) return;

    Ngs_Element ngel = ma->GetElement(ei);
    for (int i = 0; i < ngel.edges.Size(); i++)
      dranges.Append (IntRange (ngel.edges[i], ngel.edges[i]+1));
  }

  
  void NedelecFESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    if (DefinedOn (ei))
      dnums = ma->GetElEdges (ei);
    else
      dnums.SetSize0();
  }



  template <class T>
  void NedelecFESpace::T_TransformMat (ElementId ei,
                                       SliceMatrix<T> mat, TRANSFORM_TYPE tt) const
  {
    Ngs_Element ngel = ma->GetElement(ei);
    ELEMENT_TYPE eltype = ngel.GetType();
    
    int ned = ElementTopology::GetNEdges (eltype);
    auto edges = ElementTopology::GetEdges (eltype);
    ArrayMem<int,12> eorient(ned);
    for (int i = 0; i < ned; i++)
      eorient[i] = 
        ngel.vertices[edges[i][0]] < ngel.vertices[edges[i][1]]  
        ? 1 : -1;

    if (tt & TRANSFORM_MAT_LEFT)
      for (int i = 0; i < ned; i++)
        for (int k = 0; k < dimension; k++)
          mat.Row(k+i*dimension) *= eorient[i];

    if (tt & TRANSFORM_MAT_RIGHT)
      for (int j = 0; j < ned; j++)
        for (int l = 0; l < dimension; l++)
          mat.Col(l+j*dimension) *= eorient[j];
  }


  template <class T>
  void NedelecFESpace::T_TransformVec (ElementId ei, 
                                       SliceVector<T> vec, TRANSFORM_TYPE tt) const
  {
    /*
    int nd;
    ArrayMem<int,12> enums, eorient;
    LocalHeapMem<1000> lh("NedelecFESpace - transformvec");

    if (boundary)
      {
	GetMeshAccess().GetSElEdges (elnr, enums, eorient);
	nd = GetSFE (elnr, lh).GetNDof();
      }
    else
      {
	GetMeshAccess().GetElEdges (elnr, enums, eorient);
	nd = GetFE (elnr, lh).GetNDof();
      }
    */


    Ngs_Element ngel = ma->GetElement(ei);
    ELEMENT_TYPE eltype = ngel.GetType();
    
    int ned = ElementTopology::GetNEdges (eltype);
    auto edges = ElementTopology::GetEdges (eltype);
    ArrayMem<int,12> eorient(ned);
    for (int i = 0; i < ned; i++)
      eorient[i] = 
        ngel.vertices[edges[i][0]] < ngel.vertices[edges[i][1]]  
        ? 1 : -1;


    if ((tt & TRANSFORM_RHS) || (tt & TRANSFORM_SOL) || (tt & TRANSFORM_SOL_INVERSE))
      {
	for (int k = 0; k < dimension; k++)
	  for (int i = 0; i < ned; i++)
	    vec(k+i*dimension) *= eorient[i];
      }
  }

  
  shared_ptr<Table<int>> NedelecFESpace :: CreateSmoothingBlocks (const Flags & precflags) const
  {
    return CreateSmoothingBlocks (int (precflags.GetNumFlag ("loblocktype", SB_AFW)));
  }



  shared_ptr<Table<int>> NedelecFESpace :: CreateSmoothingBlocks (int type) const
  {
    cout << IM(3) << "NedelecFESpace::CreateSmoothingBlocks" << endl;

    int nd = GetNDof();
    int nv = ma->GetNV();
    int level = ma->GetNLevels()-1;
    auto freedofs = GetFreeDofs();

    Table<int> *node2edge = 0;
    //type = SB_AFW;  
    switch (type)
      {
      case SB_AFW:
	{
	  cout << IM(3) << " ******** Low-order H(Curl) Smoother: AFW" << endl;
	  Array<int> cnts(nv);
	  for (int k = 1; k <= 2; k++)
	    {
	      if (k == 2)
		node2edge = new Table<int>(cnts);
	    
	      cnts = 0;

	    
	      for (int j = 0; j < nd; j++)
		{
		  if (FineLevelOfEdge(j) < level) continue;
		  if (freedofs && !freedofs->Test(j)) continue;

		  int ep1 = EdgePoint1(j);
		  int ep2 = EdgePoint2(j);

		
		  // for anisotropic connections:
		  int cep1 = ma->GetClusterRepVertex(ep1);
		  int cep2 = ma->GetClusterRepVertex(ep2);

		  if (k == 2)
		    {
		      (*node2edge)[cep1][cnts[cep1]] = j;
		      cnts[cep1]++;
		    
		      if (cep1 != cep2)
			{
			  (*node2edge)[cep2][cnts[cep2]] = j;
			  cnts[cep2]++;
			}
		    }
		  else
		    {
		      cnts[cep1]++;
		      if (cep1 != cep2)
			cnts[cep2]++;
		    }
		
		}
	    }
	  //(*testout) << "node2egde: " << *node2edge << endl;
	  break;
	}
      case SB_JAC: // only for getting bad condition numbers  ... 
	{
	  cout << " Jacobi Smoother for Low-order H(Curl) --> bad conditoning" << endl;  
	  Array<int> cnts(nd);
	  for (int k = 1; k <= 2; k++)
	    {
	      if (k == 2)
		node2edge = new Table<int>(cnts);
	    
	      cnts = 0;

	    
	      for (int j = 0; j < nd; j++)
		{
		  if (FineLevelOfEdge(j) < level) continue;
		  if (freedofs && !freedofs->Test(j)) continue;

		  if (k == 2)
		    {
		      (*node2edge)[j][0] = j;
		   
		    }
		  else
		    {
		      cnts[j]=1; 
		    }
		
		}
	    }
	  (*testout) << "node2egde: " << *node2edge << endl;
	  break;
	}
	//     case SB_AFW:
	//       {
	// 	Array<int> cnts(nv);
	// 	for (int k = 1; k <= 2; k++)
	// 	  {
	// 	    if (k == 2)
	// 	      node2edge = new Table<int>(cnts);
	    
	// 	    cnts = 0;
	    
	// 	    for (int j = 0; j < nd; j++)
	// 	      {
	// 		if (FineLevelOfEdge(j) < level) continue;
		
	// 		int ep1 = EdgePoint1(j);
	// 		int ep2 = EdgePoint2(j);

		
	// 		// for anisotropic connections:
	// 		int cep1 = ma->GetClusterRepVertex(ep1);
	// 		int cep2 = ma->GetClusterRepVertex(ep2);

	// 		if (k == 2)
	// 		  {
	// 		    (*node2edge)[cep1][cnts[cep1]] = j;
	// 		    cnts[cep1]++;
		    
	// 		    if (cep1 != cep2)
	// 		      {
	// 			(*node2edge)[cep2][cnts[cep2]] = j;
	// 			cnts[cep2]++;
	// 		      }
	// 		  }
	// 		else
	// 		  {
	// 		    cnts[cep1]++;
	// 		    if (cep1 != cep2)
	// 		      cnts[cep2]++;
	// 		  }
		
	// 	      }
	// 	  }
	// 	//	(*testout) << "node2egde: " << *node2edge << endl;
	// 	break;
	//       }



      case SB_HIPTMAIR:
	{
	  Array<int> cnts(nd);
	  for (int k = 1; k <= 2; k++)
	    {
	      if (k == 2)
		node2edge = new Table<int>(cnts);
	    
	      cnts = 0;
	    
	      for (int j = 0; j < nd; j++)
		{
		  if (FineLevelOfEdge(j) < level) continue;
		  if (freedofs && !freedofs->Test(j)) continue;

		  int ecl = ma->GetClusterRepEdge (j);
		  if (ecl < nv)
		    ecl = j;
		  else
		    ecl -= nv;

		  if (k == 2)
		    {
		      (*node2edge)[ecl][cnts[ecl]] = j;
		      cnts[ecl]++;
		    }
		  else
		    {
		      cnts[ecl]++;
		    }
		
		}
	    }
	  break;
	}
      case SB_POTENTIAL:
	{
	  Array<int> cnts(nv);
	  for (int k = 1; k <= 2; k++)
	    {
	      if (k == 2)
		node2edge = new Table<int>(cnts);
	    
	      cnts = 0;
	    
	      for (int j = 0; j < nv; j++)
		{
		  int vcl = ma->GetClusterRepVertex (j);
		  if (k == 2)
		    {
		      (*node2edge)[vcl][cnts[vcl]] = j;
		      cnts[vcl]++;
		    }
		  else
		    {
		      cnts[vcl]++;
		    }
		}
	    }
	  break;
	}
      }
  
    return shared_ptr<Table<int>> (node2edge);
  }

  SparseMatrix<double> * 
  NedelecFESpace :: CreateGradient() const
  {
    int i;
    int ned = GetNDof();
    int level = ma->GetNLevels()-1;

    Array<int> cnts(ned);
    for (i = 0; i < ned; i++)
      cnts[i] = (FineLevelOfEdge(i) == level) ? 2 : 0;

    SparseMatrix<double> & grad = *new SparseMatrix<double>(cnts, ma->GetNV());

    for (i = 0; i < ned; i++)
      {
	if (FineLevelOfEdge(i) < level) continue;
	grad.CreatePosition (i, edgepoints[i][0]);
	grad.CreatePosition (i, edgepoints[i][1]);
      }
    for (i = 0; i < ned; i++)
      {
	if (FineLevelOfEdge(i) < level) continue;
	grad(i, edgepoints[i][0]) = 1;
	grad(i, edgepoints[i][1]) = -1;
      }

    return &grad;
  }



  void NedelecFESpace :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
  }

  void NedelecFESpace :: GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  {
    dnums.SetSize(1);
    dnums[0] = ednr;
  }

  void NedelecFESpace :: GetFaceDofNrs (int fanr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
  }

  void NedelecFESpace :: GetInnerDofNrs (int elnr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
  }



  
  
  class EdgeP1Prolongation : public Prolongation
  {
    shared_ptr<MeshAccess> ma;
    const FESpace & space;
  public:
    EdgeP1Prolongation(const FESpace & aspace)
      : ma(aspace.GetMeshAccess()), space(aspace)
    {
      ma->EnableTable("parentedges");
    }
    
    virtual ~EdgeP1Prolongation() { }
  
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
              fv(2*i)   = fac0 * fv(2*pa1) - 0.125 * fv(2*pa1+1);
              fv(2*i+1) = 0.25 * fv(2*pa1+1);
            }
          else if (info<8)//bisecting edge
            {
              double fac1 = (info&1) ? 0.5 : -0.5;
              double fac2 = (info&2) ? 0.5 : -0.5;
              double fac3 = (info&4) ? -0.125 : 0.125;
              fv(2*i) = fac1 * fv(2*pa1) + fac2 * fv(2*pa2) + fac3 * fv(2*pa3+1);
              fv(2*i+1) = 0.5 * (fv(2*pa1+1)+fv(2*pa2+1)) - 0.25*fv(2*pa3+1);
            }
          else // info>=8: red edge
            {
              double fac1 = (info&1) ? 0.25 : -0.25;
              double fac2 = (info&2) ? 0.25 : -0.25;
              double fac3 = (info&4) ? 0.25 : -0.25;
              fv(2*i) = fac1 * fv(2*pa1) + fac2 * fv(2*pa2) + fac3 * fv(2*pa3)
                + 0.125 * fv(2*pa1+1) - 0.125 * fv(2*pa2+1);
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
              fv(2*pa1+1) += -0.125 * fv(2*i) + 0.25 * fv(2*i+1);
            }
          else if (info<8)//bisecting edge
            {
              double fac1 = (info&1) ? 0.5 : -0.5;
              double fac2 = (info&2) ? 0.5 : -0.5;
              double fac3 = (info&4) ? -0.125 : 0.125;
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
              fv(2*pa1+1) += 0.125 * fv(2*i);
              fv(2*pa2)   += fac2 * fv(2*i);
              fv(2*pa2+1) -= 0.125 * fv(2*i);
              fv(2*pa3) += fac3 * fv(2*i);
              fv(2*pa3+1) += 0.25*fv(2*i+1);
            }
        }
      
    }
  };

}

namespace ngfem {
  class NedelecP1Trig : public T_HCurlFiniteElementFO<NedelecP1Trig,ET_TRIG,6,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (TIP<2,Tx> ip, TFA & shape) 
    {
      // Tx x = hx[0], y = hx[1];
      Tx x = ip.x, y = ip.y;
      Tx lami[3] = { x, y, 1-x-y };
      
      auto edges = ElementTopology::GetEdges (ET_TRIG);
      for (int i = 0; i < 3; i++)
        {
          shape[i] = uDv_minus_vDu (lami[edges[i][0]], lami[edges[i][1]]);
          shape[i+3] = Du (-0.5*lami[edges[i][0]]*lami[edges[i][1]]);
        }
    }
  };

  class NedelecP1Tet : public T_HCurlFiniteElementFO<NedelecP1Tet,ET_TET,12,1>
  {
  public:
    template<typename Tx, typename TFA>  
    static void T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
    {
      Tx lami[4] = { ip.x, ip.y, ip.z, 1-ip.x-ip.y-ip.z };
      
      auto edges = ElementTopology::GetEdges (ET_TET);
      for (int i = 0; i < 6; i++)
        {
          shape[i] = uDv_minus_vDu (lami[edges[i][0]], lami[edges[i][1]]);
          shape[i+6] = Du (-0.5*lami[edges[i][0]]*lami[edges[i][1]]);
        }
    }
  };
  
  // template class T_HCurlHighOrderFiniteElement<ET_TRIG,NedelecP1Trig>;
  // template class T_HCurlHighOrderFiniteElement<ET_TET,NedelecP1Tet>;
}

namespace ngcomp {
  
  class NGS_DLL_HEADER NedelecP1FESpace : public FESpace
  {
    BitArray active_edges;
  public:
    NedelecP1FESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags=false)
      : FESpace(ama, flags)
      {
        name="NedelecP1FESpace";
        
        if (ma->GetDimension() == 2)
          {
            evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryEdge<2>>>();
            evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdEdge<2>>>();
            flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpCurlEdge<2>>>();        

            additional_evaluators.Set ("grad", make_shared<T_DifferentialOperator<DiffOpGradientHCurl<2>>> ());
          }
        else if(ma->GetDimension() == 3) 
          {
            evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryEdge<3>>>();
            evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdEdge<3>>>();
            flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpCurlEdge<3>>>();
            flux_evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpCurlBoundaryEdgeVec<>>>();
            evaluator[BBND] = make_shared<T_DifferentialOperator<DiffOpIdBBoundaryEdge<3>>>();

            additional_evaluators.Set ("grad", make_shared<T_DifferentialOperator<DiffOpGradientHCurl<3>>> ());            
          }
        prol = make_shared<EdgeP1Prolongation> (*this);
      }
    
    virtual ~NedelecP1FESpace () { }
    virtual const char * GetType()  { return "NedelecP1"; }

    /*
    static shared_ptr<FESpace> Create (shared_ptr<MeshAccess> ma, const Flags & flags)
    {
      return make_shared<NedelecFESpace2> (ma, flags, true);
    }
    */
    void Update() override
    {
      size_t ned = ma->GetNEdges();
      SetNDof (2*ned);
      active_edges = BitArray(ned);
      active_edges.Clear();
      for (auto el : ma->Elements(VOL))
        for (auto ed : el.Edges())
          active_edges.SetBit(ed);
      
      ctofdof.SetSize(GetNDof());
      ctofdof = WIREBASKET_DOF;
      for (size_t i = 0; i < ned; i++)
        if (!active_edges.Test(i))
          ctofdof[2*i] = ctofdof[2*i+1] = UNUSED_DOF;
      // cout << "active edges = " << endl << active_edges << endl;
    }
    
    // virtual void DoArchive (Archive & archive) override;
    // virtual void UpdateCouplingDofArray() override;
    
    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const override
    {
      switch (ma->GetElType(ei))
        {
        case ET_TET:     return * new (lh) NedelecP1Tet;
        case ET_TRIG:    return * new (lh) NedelecP1Trig;
        default:
          throw Exception ("Inconsistent element type in NedelecFESpace::GetFE");
        }
    }
    
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override
    {
      auto edges = ma->GetElEdges (ei);
      dnums.SetSize(2*edges.Size());
      for (int i : Range(edges))
        {
          dnums[i] = 2*edges[i];
          dnums[i+edges.Size()] = 2*edges[i]+1;
        }
    }

    // virtual shared_ptr<Table<int>> CreateSmoothingBlocks (const Flags & precflags) const override;

    template <class T>
      NGS_DLL_HEADER void T_TransformMat (ElementId ei, 
                                          SliceMatrix<T> mat, TRANSFORM_TYPE tt) const
    {
      Ngs_Element ngel = ma->GetElement(ei);
      ELEMENT_TYPE eltype = ngel.GetType();
      
      int ned = ElementTopology::GetNEdges (eltype);
      auto edges = ElementTopology::GetEdges (eltype);
      ArrayMem<int,12> eorient(ned);
      for (int i = 0; i < ned; i++)
        eorient[i] = 
          ngel.vertices[edges[i][0]] < ngel.vertices[edges[i][1]]  
                                       ? 1 : -1;
      
      if (tt & TRANSFORM_MAT_LEFT)
        for (int i = 0; i < ned; i++)
          for (int k = 0; k < dimension; k++)
            mat.Row(k+i*dimension) *= eorient[i];
      
      if (tt & TRANSFORM_MAT_RIGHT)
        for (int j = 0; j < ned; j++)
          for (int l = 0; l < dimension; l++)
            mat.Col(l+j*dimension) *= eorient[j];
    }
    
    template <class T>
      NGS_DLL_HEADER void T_TransformVec (ElementId ei, 
                                          SliceVector<T> vec, TRANSFORM_TYPE tt) const
    {
      Ngs_Element ngel = ma->GetElement(ei);
      ELEMENT_TYPE eltype = ngel.GetType();
      
      int ned = ElementTopology::GetNEdges (eltype);
      auto edges = ElementTopology::GetEdges (eltype);
      ArrayMem<int,12> eorient(ned);
      for (int i = 0; i < ned; i++)
        eorient[i] = 
          ngel.vertices[edges[i][0]] < ngel.vertices[edges[i][1]]  
        ? 1 : -1;
      
      
      if ((tt & TRANSFORM_RHS) || (tt & TRANSFORM_SOL) || (tt & TRANSFORM_SOL_INVERSE))
        {
          for (int k = 0; k < dimension; k++)
            for (int i = 0; i < ned; i++)
              vec(k+i*dimension) *= eorient[i];
        }
    }
    
    
    virtual void VTransformMR (ElementId ei, 
                               SliceMatrix<double> mat, TRANSFORM_TYPE tt) const override
    {
      T_TransformMat (ei, mat, tt);
    }
    
    virtual void VTransformMC (ElementId ei, 
                               SliceMatrix<Complex> mat, TRANSFORM_TYPE tt) const override
    {
      T_TransformMat (ei, mat, tt);
    }
    
    virtual void VTransformVR (ElementId ei, 
                               SliceVector<double> vec, TRANSFORM_TYPE tt) const override
    {
      T_TransformVec (ei, vec, tt);
    }
    
    virtual void VTransformVC (ElementId ei, 
                               SliceVector<Complex> vec, TRANSFORM_TYPE tt) const override
    {
      T_TransformVec (ei, vec, tt);
    }
    
    virtual string GetClassName () const override
    {
      return "NedelecP1FESpace";
    }
    
    virtual void GetVertexDofNrs (int vnr, Array<DofId> & dnums) const override
    { dnums.SetSize0(); }
    virtual void GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const override
    {
      if (active_edges.Test(ednr))
        {
          dnums.SetSize(2);
          dnums[0] = 2*ednr;
          dnums[1] = 2*ednr+1;
        }
      else
        dnums.SetSize0();
    }
    virtual void GetFaceDofNrs (int fanr, Array<DofId> & dnums) const override
    { dnums.SetSize0(); }    
    virtual void GetInnerDofNrs (int elnr, Array<DofId> & dnums) const override
    { dnums.SetSize0(); }    
  };

  static RegisterFESpace<NedelecP1FESpace> initnedelec ("HCurlP1");






  // register FESpaces
  namespace hcurlhdives_cpp
  {
    class Init
    { 
    public: 
      Init ();
    };
    
    Init::Init()
    {
      GetFESpaceClasses().AddFESpace ("hcurl", NedelecFESpace::Create,
                                      NedelecFESpace::GetDocu);
    }

    
    Init init;
  }



}
