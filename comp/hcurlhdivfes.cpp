/*********************************************************************/
/* File:   hcurlhdivfes.cpp                                          */
/* Author: Joachim Schoeberl                                         */
/* Date:   12. Jan. 2002                                             */
/*********************************************************************/

/* 
   Finite Element Space
*/

#include <comp.hpp>
#include <multigrid.hpp>

#include <../fem/hcurllofe.hpp>
#include <../fem/thcurlfe_impl.hpp>


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

    auto one = make_shared<ConstantCoefficientFunction>(1);
    integrator[VOL] = GetIntegrators().CreateBFI("massedge", ma->GetDimension(), one);
    integrator[BND] = GetIntegrators().CreateBFI("robinedge", ma->GetDimension(), one);

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
    int order = int(flags.GetNumFlag ("order", 1)); 
    if (order >= 2) 
      return make_shared<NedelecFESpace2> (ma, flags, true);
    else
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
    ClosedHashTable<INT<2>, int> node2edge(5*ned+10);

    edgepoints.SetSize0();
    
    for (size_t i = 0; i < ned; i++)
      {
	INT<2> edge = ma->GetEdgePNums (i);
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
    parentedges = INT<2> (-1,-1);

    for (size_t i = 0; i < ned; i++)
      {
	// cout << "edge " << i << "/" << ned << endl;
	INT<2> i2 (edgepoints[i][0], edgepoints[i][1]);
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
	    INT<2> paedge;
	    if (issplitedge == 1)
	      paedge = INT<2> (pa1[0], pa1[1]);
	    else
	      paedge = INT<2> (pa2[0], pa2[1]);
	    
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
		INT<2> paedge1, paedge2;
		if (j == 1)
		  {
		    paedge1 = INT<2> (pa1[0], i2[1]);
		    paedge2 = INT<2> (pa1[1], i2[1]);
		  }
		else
		  {
		    paedge1 = INT<2> (pa2[0], i2[0]);
		    paedge2 = INT<2> (pa2[1], i2[0]);
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
		      INT<2> paedge1, paedge2;
		      if (j == 1)
			{
			  paedge1 = INT<2> (pa1[0], pa2[0]);
			  paedge2 = INT<2> (pa1[1], pa2[1]);
			}
		      else
			{
			  paedge1 = INT<2> (pa1[0], pa2[1]);
			  paedge2 = INT<2> (pa1[1], pa2[0]);
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
		      INT<2> paedge (pa1[1-j], pa2[1-k]);
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
    const EDGE * edges = ElementTopology::GetEdges (eltype);
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
    const EDGE * edges = ElementTopology::GetEdges (eltype);
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
      : ma(aspace.GetMeshAccess()), space(aspace) { ; }
    virtual ~EdgeP1Prolongation() { }
  
    virtual void Update (const FESpace & fes) { ; }
    virtual SparseMatrix< double >* CreateProlongationMatrix( int finelevel ) const
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
          else
            {
              double fac1 = (info&1) ? 0.5 : -0.5;
              double fac2 = (info&2) ? 0.5 : -0.5;
              double fac3 = (info&4) ? -0.125 : 0.125;
              fv(2*i) = fac1 * fv(2*pa1) + fac2 * fv(2*pa2) + fac3 * fv(2*pa3+1);
              fv(2*i+1) = 0.5 * (fv(2*pa1+1)+fv(2*pa2+1)) - 0.25*fv(2*pa3+1);
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
          else
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
      
      const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
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
      
      const EDGE * edges = ElementTopology::GetEdges (ET_TET);
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
    
    static shared_ptr<FESpace> Create (shared_ptr<MeshAccess> ma, const Flags & flags)
    {
      return make_shared<NedelecFESpace2> (ma, flags, true);
    }
    
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
      const EDGE * edges = ElementTopology::GetEdges (eltype);
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
      const EDGE * edges = ElementTopology::GetEdges (eltype);
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




  

  NedelecFESpace2 :: NedelecFESpace2 (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags)
    : FESpace (ama, flags)
  {
    name="NedelecFESpace2(hcurl)";
    // defined flags
    DefineDefineFlag("hcurl");
    DefineNumFlag("zorder");
    DefineNumListFlag("gradientdomains");
    DefineNumListFlag("gradientboundaries");
    DefineNumListFlag("direcsolverdomains");
    DefineNumListFlag("directsolvermaterials");
  
    if(parseflags) CheckFlags(flags);
  
    order = (int) flags.GetNumFlag("order",2);
    zorder = int(flags.GetNumFlag ("zorder", order));

    ned = 0;

    n_plane_edge_dofs = order;
    n_z_edge_dofs = zorder;
    n_edge_dofs = max2 (n_z_edge_dofs, n_plane_edge_dofs);
    n_pyramid_el_dofs = 0;


    curltet = 0;
    curlprism = 0;
    curlpyramid = 0;
  
    switch (order)
      {
      case 2:
	{
	  segm    = new FE_NedelecSegm2;
	  trig    = new FE_NedelecTrig2;
	  tet     = new FE_NedelecTet2;

	  n_trig_face_dofs = 0;
	  n_tet_el_dofs = 0;

	  n_quad_face_dofs = 3*zorder-2;
	  n_prism_el_dofs = 0;
	  n_prism_nograd_el_dofs = 0;
	  n_pyramid_el_dofs = 0;

	  switch (zorder)
	    {
	    case 1:
	      {
		quad    = new FE_TNedelecQuad<2,1>;
		prism   = new FE_TNedelecPrism2<1>;
		//	      pyramid = new FE_NedelecPyramid2;
		break;
	      }
	    case 2:
	      {
		quad    = new FE_TNedelecQuad<2,2>;
		prism   = new FE_TNedelecPrism2<2>;
		pyramid = new FE_NedelecPyramid2;
		break;
	      }
	    case 3:
	      {
		quad    = new FE_TNedelecQuad<2,3>;
		prism   = new FE_TNedelecPrism2<3>;
		break;
	      }
	    case 4:
	      {
		quad    = new FE_TNedelecQuad<2,4>;
		prism   = new FE_TNedelecPrism2<4>;
		break;
	      }
	    }
	  break;
	}
      case 3:
	{
	  segm    = new FE_NedelecSegm3;
	  trig    = new FE_NedelecTrig3;
	  tet     = new FE_NedelecTet3;
	  curltet     = new FE_NedelecTet3NoGrad;

	  n_trig_face_dofs = 3;
	  n_tet_el_dofs = 0;
	  n_quad_face_dofs = 5*zorder-3;
	  n_prism_el_dofs = 4*zorder-3;
	  n_prism_nograd_el_dofs = 3*zorder-2;
	  n_pyramid_el_dofs = 9;

	  switch (zorder)
	    {
	    case 1:
	      {
		quad    = new FE_TNedelecQuad<3,1>;
		prism   = new FE_TNedelecPrism3<1>;	      
		curlprism   = new FE_TNedelecPrism3NoGrad<1>;
		break;
	      }
	    case 2:
	      {
		quad    = new FE_TNedelecQuad<3,2>;
		prism   = new FE_TNedelecPrism3<2>;
		curlprism   = new FE_TNedelecPrism3NoGrad<2>;
		break;
	      }
	    case 3:
	      {
		quad    = new FE_TNedelecQuad<3,3>;
		prism   = new FE_TNedelecPrism3<3>;
		// pyramid = new FE_NedelecPyramid3;
		curlprism   = new FE_TNedelecPrism3NoGrad<3>;
		break;
	      }
	    }
	  break;
	}
      }

    if (!curltet) curltet = tet;
    if (!curlprism) curlprism = prism;


    gradientdomains.SetSize (ma->GetNDomains());
    gradientdomains.Set();
    gradientboundaries.SetSize (ma->GetNBoundaries());
    gradientboundaries.Set();

    if (flags.NumListFlagDefined("gradientdomains"))
      {
        cout << "has gradientdomains" << endl;
	const Array<double> & graddomains = flags.GetNumListFlag ("gradientdomains");
	for (int i = 0; i < gradientdomains.Size(); i++)
	  if (!graddomains[i])
	    gradientdomains.Clear(i);
      }

    if (flags.NumListFlagDefined("gradientboundaries"))
      {
	const Array<double> & gradbounds = flags.GetNumListFlag ("gradientboundaries");
	for (int i = 0; i < gradientboundaries.Size(); i++)
	  if (!gradbounds[i])
	    gradientboundaries.Clear(i);
      }


    cout << "gradientdomains = " << gradientdomains << endl;

    Flags loflags;
    loflags.SetFlag ("order", 1);
    loflags.SetFlag ("dim", dimension);
    if (iscomplex) loflags.SetFlag ("complex");
    low_order_space = make_shared<NedelecFESpace> (ma, loflags);
    // low_order_space = new NedelecFESpace (ama, 1, dimension, iscomplex);

    prol = make_shared<EdgeProlongation> (*static_cast<NedelecFESpace*> (low_order_space.get()));
    /*
      CreateVecObject1(prol, EdgeProlongation,
      dimension, iscomplex, 
      *static_cast<NedelecFESpace*> (low_order_space));
      */

    // Evaluator for shape tester 
    if (ma->GetDimension() == 2)
      {
	Array<shared_ptr<CoefficientFunction>> coeffs(1);
	coeffs[0] = shared_ptr<CoefficientFunction> (new ConstantCoefficientFunction(1));
	integrator[VOL] = GetIntegrators().CreateBFI("massedge", 2, coeffs);
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryEdge<2>>>();
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdEdge<2>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpCurlEdge<2>>>();        
      }
    else if(ma->GetDimension() == 3) 
      {
	Array<shared_ptr<CoefficientFunction>> coeffs(1); 
	coeffs[0] = shared_ptr<CoefficientFunction> (new ConstantCoefficientFunction(1)); 
	integrator[VOL] = GetIntegrators().CreateBFI("massedge",3,coeffs); 
	integrator[BND] = GetIntegrators().CreateBFI("robinedge",3,coeffs);
        
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundaryEdge<3>>>();
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdEdge<3>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpCurlEdge<3>>>();
        flux_evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpCurlBoundaryEdgeVec<>>>();
	evaluator[BBND] = make_shared<T_DifferentialOperator<DiffOpIdBBoundaryEdge<3>>>();	
      }


    if(flags.NumListFlagDefined("directsolverdomains"))
      {
	directsolverclustered.SetSize(ama->GetNDomains());
	directsolverclustered = false;
	Array<double> clusters(flags.GetNumListFlag("directsolverdomains"));
	for(int i=0; i<clusters.Size(); i++) 
	  directsolverclustered[static_cast<int>(clusters[i])-1] = true; // 1-based!!
      }
    
    if(flags.StringListFlagDefined("directsolvermaterials"))
      {
	directsolvermaterials.SetSize(flags.GetStringListFlag("directsolvermaterials").Size());
	for(int i=0; i<directsolvermaterials.Size(); i++)
	  directsolvermaterials[i] = flags.GetStringListFlag("directsolvermaterials")[i];
      }

  }


  NedelecFESpace2 :: ~NedelecFESpace2 ()
  {
    //delete low_order_space;
  }

  void NedelecFESpace2 :: Update()
  {
    if (low_order_space)
      low_order_space -> Update();

    int i, j;

    int level = ma->GetNLevels();
    int ne = ma->GetNE();
    int nse = ma->GetNSE();
    // int np = ma->GetNP();

    ned = ma->GetNEdges();
    nfa = ma->GetNFaces(); 
    nel = ma->GetNE(); 

    if (ma->GetDimension() == 2)
      nfa = nel;

    if (gradientedge.Size() == ned)
      return;

    // new definition of gradient edges

    Array<int> enums, fnums, forient;

    gradientedge.SetSize(ned);
    gradientedge = 1;
    gradientface.SetSize(nfa);
    gradientface = 1;

    /*
      for (i = 0; i < nse; i++)
      {
      ma->GetSElEdges (i, enums, eorient);
      for (j = 0; j < enums.Size(); j++)
      gradientedge[enums[j]] = 1;
      }
    */


    first_face_dof.SetSize (nfa+1);
    first_el_dof.SetSize (nel+1);

    first_face_dof[0] = n_edge_dofs * ned;
    for (i = 0; i < nfa; i++)
      {
	auto pnums = ma->GetFacePNums (i);
	if (pnums.Size() == 3)
	  first_face_dof[i+1] = first_face_dof[i] + n_trig_face_dofs;
	else
	  first_face_dof[i+1] = first_face_dof[i] + n_quad_face_dofs;
      }

    first_el_dof[0] = first_face_dof[nfa];
    for (i = 0; i < ne; i++)
      {
        ElementId ei(VOL,i);
	bool gradel = gradientdomains[ma->GetElIndex(ei)];
	switch (ma->GetElType (ei))
	  {
	  case ET_TET:
	    first_el_dof[i+1] = first_el_dof[i] + n_tet_el_dofs; 
	    break;
	  case ET_PRISM:
	    {
	      if (gradel) first_el_dof[i+1] = first_el_dof[i] + n_prism_el_dofs; 
	      else first_el_dof[i+1] = first_el_dof[i] + n_prism_nograd_el_dofs; 
	      break;
	    }
	  case ET_PYRAMID:
	    first_el_dof[i+1] = first_el_dof[i] + n_pyramid_el_dofs; 
	    break;
	  default:
	    throw Exception ("unhandled case in NedelecFESpace2::Update");
	  }
      }

    if (level != ndlevel.Size())
      ndlevel.Append (first_el_dof[nel]);
  

    if (gradientdomains.Size())
      {
	gradientedge = 0;
	gradientface = 0;

	for (i = 0; i < ne; i++)
	  {
            ElementId ei(VOL,i);
	    if (gradientdomains[ma->GetElIndex(ei)])
	      {
		auto enums = ma->GetElEdges (ei);
		for (j = 0; j < enums.Size(); j++)
		  gradientedge[enums[j]] = 1;
		auto fnums = ma->GetElFaces (ei);
		for (j = 0; j < fnums.Size(); j++)
		  gradientface[fnums[j]] = 1;
	      }
	  }
      }

    fnums.SetSize(1);
    forient.SetSize(1);
    (*testout) << "gradientboundaries = " << endl << gradientboundaries << endl;
    if (gradientboundaries.Size())
      for (i = 0; i < nse; i++)
	{
          ElementId sei(BND,i);
	  if (gradientboundaries[ma->GetElIndex(sei)])
	    {
	      auto enums = ma->GetElEdges (sei);
	      for (j = 0; j < enums.Size(); j++)
		gradientedge[enums[j]] = 1;
	      ma->GetSElFace (i, fnums[0], forient[0]);
	      gradientface[fnums[0]] = 1;
	    }
	}

    //  gradientface = 0;
    //  (*testout) << "gradientedges = " << endl << gradientedge << endl;


    // FinalizeUpdate (lh);
  }





  size_t NedelecFESpace2 :: GetNDof () const throw()
  {
    return ndlevel.Last();
  }

  size_t NedelecFESpace2 :: GetNDofLevel (int level) const
  {
    return ndlevel[level];
  }



  FiniteElement & NedelecFESpace2 :: GetFE (ElementId ei, Allocator & lh) const
  {
    FiniteElement * fe = 0;
    ELEMENT_TYPE typ = ma->GetElType(ei);

    switch (typ)
      {
      case ET_TET:
	fe = tet; break;
      case ET_PYRAMID:
	fe = pyramid; break;
      case ET_PRISM:
	fe = prism; break;
      case ET_HEX:
	fe = hex; break;
      case ET_TRIG:
	fe = trig; break;
      case ET_QUAD:
	fe = quad; break;
      default:
	fe = 0;
      }
    
    if (!gradientdomains[ma->GetElIndex(ei)])
      {
	switch (typ)
	  {
	  case ET_TET:
	    fe = curltet; break;
	  case ET_PRISM:
	    fe = curlprism; break;
	    //	case ET_PYRAMID:
	    //	  fe = curlpyramid; break;
	  default:
	    ;
	  }
      }
  
    if (!fe)
      {
	stringstream str;
	str << "FESpace " << GetClassName() 
	    << ", undefined eltype " 
	    << ElementTopology::GetElementName(ma->GetElType(ei))
	    << ", order = " << order << endl;
	throw Exception (str.str());
      }
  
    return *fe;
  }

  
  void NedelecFESpace2 :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    if(ei.VB()==VOL)
      {
	// int eled, elfa;
	int j;

	ArrayMem<int,6> fnums, forient;
	// ArrayMem<int,12> enums;
	
	auto enums = ma->GetElEdges (ei);
	ma->GetElFaces (ei.Nr(), fnums, forient);
	
	LocalHeapMem<1000> lh("NedelecFESpace2, GetDofNrs");
	int nd = GetFE (ei, lh).GetNDof();
	dnums.SetSize(nd);
	dnums = -1;
	
	int index = ma->GetElIndex (ei);
	
	if (!DefinedOn (VOL, index)) return;
	
	bool graddom = gradientdomains[index];
	
	switch (ma->GetElType(ei))
	  {
	  case ET_TRIG:
	    {
	      switch (order)
		{
		case 2:
		  {
		    for (j = 0; j < 3; j++)
		      {
			dnums[j] = enums[j];
			if (gradientedge[enums[j]])
			  dnums[j+3] = enums[j] + ned;
		      }
		    break;
		  }	
		case 3:
		  {
		    for (j = 0; j < 3; j++)
		      {
			int edgenr = enums[j];
			dnums[j] = edgenr;
			if (gradientedge[edgenr])
			  {
			    dnums[j+3] = edgenr + ned;
			    dnums[j+6] = edgenr + 2*ned;
			  }
		      }
		    int nfd = gradientface[ei.Nr()] ? 3 : 2;
		    for (j = 0; j < nfd; j++)
		      dnums[9+j] = first_el_dof[ei.Nr()];
		    break;
		  }
		}
	      break;
	    }
	  case ET_TET:
	    {
	      switch (order)
		{
		case 2:
		  {
		    for (j = 0; j < 6; j++)
		      {
			dnums[j] = enums[j];
			if (gradientedge[enums[j]])
			  dnums[j+6] = enums[j] + ned;
		      }
		    break;
		  }
		case 3:
		  {
		    for (j = 0; j < 6; j++)
		      dnums[j] = enums[j];
		    
		    int base = 6;
		    
		    if (nd == 30)
		      {
			for (j = 0; j < 6; j++)
			  {
			    int edgenr = enums[j];
			    if (gradientedge[edgenr] && nd == 30)
			      {
				dnums[base+j]  = edgenr + ned;
			    dnums[base+6+j] = edgenr + 2*ned;
			      }
			  }
			base += 12;
		      }
		    
		    for (j = 0; j < 4; j++)
		      {
			int facedir = forient[j];
			
			static const int reorder[8][3] =
			  { { 1, 2, 3 }, { 2, 1, 3 }, { 1, 3, 2 }, { 2, 3, 1 },
			    { 2, 1, 3 }, { 1, 2, 3 }, { 3, 1, 2 }, { 3, 2, 1 } };
			
			int facebase = first_face_dof[fnums[j]];
			
			dnums[base + 3*j+reorder[facedir][0]-1] = facebase;
			dnums[base + 3*j+reorder[facedir][1]-1] = facebase + 1;
			if (gradientface[fnums[j]])
			  dnums[base + 3*j+reorder[facedir][2]-1] = facebase + 2;
		      }	  
		    break;
		  }
		}
	      break;
	    }
	  case ET_PRISM:
	    {
	      switch (order)
		{
		case 2:
		  {
		    int j, k, ii = 0;
		    // all edges
		    for (j = 0; j < 9; j++)
		      dnums[ii++] = enums[j];
		    
		    // horizontal edges
		    for (j = 0; j < 6; j++)
		      dnums[ii++] = 
			(gradientedge[enums[j]]) ? enums[j] + ned : -1;
		    
		    // vertical edges
		    for (j = 6; j < 9; j++)
		      for (k = 0; k < zorder-1; k++)
			dnums[ii++] = 
			  (gradientedge[enums[j]]) ? enums[j] + ned*(k+1) : -1;
		    
		    // quad faces:
		    for (j = 0; j < 3; j++)
		      {
			int facebase = first_face_dof[fnums[j+2]];
			for (k = 0; k < n_quad_face_dofs; k++)
			  dnums[ii++] = facebase + k;
		      }		    
		    break;
		  }
		  
		case 3:
		  {
		    int j, k, ii = 0;
		    
		    // all edges
		    for (j = 0; j < 9; j++)
		      dnums[ii++] = enums[j];
		    
		    if (graddom)
		      {
			// horizontal edges
			for (j = 0; j < 6; j++)
			  if (gradientedge[enums[j]])
			    {
			      dnums[ii++] = enums[j] + ned;
			      dnums[ii++] = enums[j] + 2*ned;
			    }
			  else ii+=2;
			
			// vertical edges
			for (j = 6; j < 9; j++)
			  if (gradientedge[enums[j]])
			    for (k = 0; k < zorder-1; k++)
			      dnums[ii++] = enums[j] + ned*(k+1);
			  else
			    ii += zorder-1;
		      }
		    
		    // trig faces:
		    for (j = 0; j < 2; j++)
		      {
			int facedir = forient[j];
			
			static const int reorder[8][3] =
			  { { 1, 2, 3 }, { 2, 1, 3 }, { 1, 3, 2 }, { 2, 3, 1 },
			    { 2, 1, 3 }, { 1, 2, 3 }, { 3, 1, 2 }, { 3, 2, 1 } };
			
			int facebase = first_face_dof[fnums[j]];
			
			dnums[ii+reorder[facedir][0]-1] = facebase;
			dnums[ii+reorder[facedir][1]-1] = facebase + 1;
			if (gradientface[fnums[j]])
			  dnums[ii+reorder[facedir][2]-1] = facebase + 2;
			ii += 3;
		      }
		    
		    // quad faces:
		    for (j = 0; j < 3; j++)
		      {
			int facebase = first_face_dof[fnums[j+2]];
			for (k = 0; k < n_quad_face_dofs; k++)
			  dnums[ii++] = facebase + k;
			if (!gradientface[fnums[j+2]]) dnums[ii-1] = -1;
		      }		    
		    
		    // vol dofs:
		    int elbase = first_el_dof[ei.Nr()];
		    int next = first_el_dof[ei.Nr()+1];
		    for (k = elbase; k < next; k++)
		      dnums[ii++] = k;
		    
		    break;
		  }
		}
	      break;
	    }
	    
	    
	  case ET_PYRAMID:
	    {
	      switch (order)
		{
		case 2:
		  {
		    int j, k; 
		    // all edges
		    for (j = 0; j < 8; j++)
		      dnums[j] = enums[j];
		    for (j = 0; j < 8; j++)
		      if (gradientedge[enums[j]])
			dnums[8+j] = enums[j] + ned;
		    
		    // quad face:
		    int facebase = first_face_dof[fnums[4]];
		    for (k = 0; k < 4; k++)
		      dnums[16+k] = facebase + k;
		    
		    break;
		  }
		case 3:
		  {
		    int j, k;
		    for (j = 0; j < 8; j++)
		      {
			int edgenr = enums[j];
			dnums[j] = edgenr;
			if (gradientedge[edgenr])
			  {
			    dnums[j+8] = edgenr + ned;
			    dnums[j+16] = edgenr + 2*ned;
			  }
		      }
		    int ii = 24;
		    
		    for (j = 0; j < 4; j++)
		      {
			int facedir = forient[j];
			
			static const int reorder[8][3] =
			  { { 1, 2, 3 }, { 2, 1, 3 }, { 1, 3, 2 }, { 2, 3, 1 },
			    { 2, 1, 3 }, { 1, 2, 3 }, { 3, 1, 2 }, { 3, 2, 1 } };
			
			int facebase = first_face_dof[fnums[j]];
			
			dnums[ii + reorder[facedir][0]-1] = facebase;
			dnums[ii + reorder[facedir][1]-1] = facebase + 1;
			if (gradientface[fnums[j]])
			  dnums[ii + reorder[facedir][2]-1] = facebase + 2;
			ii += 3;
		      }	  
		    
		    // quad face:
		    int facebase = first_face_dof[fnums[4]];
		    for (k = 0; k < n_quad_face_dofs; k++)
		      dnums[ii++] = facebase + k;
		    if (!gradientface[fnums[4]]) dnums[ii-1] = -1;
		    
		    for (k = 0; k < n_pyramid_el_dofs; k++)
		      dnums[ii++] = first_el_dof[ei.Nr()]+k;
		    
		    break;
		  }
		}
	      break;
	    }
	  default:
	    {
	      cerr << "NedelecFE2, GetDofNrs, unknown element" << endl;
	    }
	  }
      }
    if(ei.VB()==BND)
      {

    int fnum, forient;

    ma->GetSElFace (ei.Nr(), fnum, forient);
    auto enums = ma->GetElEdges (ei);

    LocalHeapMem<1000> lh("NedelecFESpace2, GetSDofNrs");
    int nd = GetFE (ElementId(BND,ei.Nr()), lh).GetNDof();
    dnums.SetSize(nd);
    dnums = -1;

    if (!DefinedOn (BND, ma->GetElIndex (ei)))
      return;

    switch (ma->GetElType(ei))
      {
      case ET_TRIG:
	{
	  switch (order)
	    {
	    case 2:
	      {
		for (int j = 0; j < 3; j++)
		  {
		    dnums[j] = enums[j];
		    if (gradientedge[enums[j]])
		      dnums[j+3] = enums[j] + ned;
		  }
		break;
	      }
	    case 3:
	      {
		int j;
		for (j = 0; j < 3; j++)
		  {
		    int edgenr = enums[j];
		    dnums[j] = edgenr;
		    if (gradientedge[edgenr])
		      {
			dnums[j+3] = edgenr + ned;
			dnums[j+6] = edgenr + 2*ned;
		      }
		  }

		int facedir = forient;
	      
		static const int reorder[8][3] =
		  { { 1, 2, 3 }, { 2, 1, 3 }, { 1, 3, 2 }, { 2, 3, 1 },
		    { 2, 1, 3 }, { 1, 2, 3 }, { 3, 1, 2 }, { 3, 2, 1 } };
	      
		int facebase = first_face_dof[fnum];
		dnums[9+reorder[facedir][0]-1] = facebase;
		dnums[9+reorder[facedir][1]-1] = facebase + 1;
		if (gradientface[fnum])
		  dnums[9+reorder[facedir][2]-1] = facebase + 2;
		break;
	      }
	    }
	  break;
	}
      case ET_QUAD:
	{
	  switch (order)
	    {
	    case 2:
	      {
		int j, k, ii = 0;
		// all edges
		for (j = 0; j < 4; j++)
		  dnums[j] = enums[j];
		// horizontal edges
		for (j = 0; j < 2; j++)
		  if (gradientedge[enums[j]])
		    dnums[4+j] = enums[j] + ned;
		ii = 6;
		// vertical edges
		for (j = 2; j < 4; j++)
		  if (gradientedge[enums[j]])
		    for (k = 0; k < zorder-1; k++)
		      dnums[ii++] = enums[j] + ned*(k+1);
		  else
		    ii += zorder-1;

		int facebase = first_face_dof[fnum];
		for (k = 0; k < n_quad_face_dofs; k++)
		  dnums[ii++] = facebase + k;
		break;
	      }
	    case 3:
	      {
		int j, k, ii = 0;
		// all edges
		for (j = 0; j < 4; j++)
		  dnums[ii++] = enums[j];

		// horizontal edges
		for (j = 0; j < 2; j++)
		  if (gradientedge[enums[j]])
		    {
		      dnums[ii++] = enums[j] + ned;
		      dnums[ii++] = enums[j] + 2*ned;
		    }
		  else 
		    ii+=2;

		// vertical edges
		for (j = 2; j < 4; j++)
		  if (gradientedge[enums[j]])
		    for (k = 0; k < zorder-1; k++)
		      dnums[ii++] = enums[j] + ned*(k+1);
		  else
		    ii+=zorder-1;

		int facebase = first_face_dof[fnum];
		for (k = 0; k < n_quad_face_dofs; k++)
		  dnums[ii++] = facebase + k;
		if (!gradientface[fnum]) dnums[ii-1] = -1;
		break;
	      }
	    }
	  break;
	}
      case ET_SEGM:
	{
	  for (int k = 0; k < order; k++)
	    dnums[k] = enums[0] + k * ned;
	  break;
	}
      default:
	{
	  throw Exception ("Unhandled Element in GetSDofNrs");
	}
      }

      }
    if(ei.VB()==BBND)
      throw Exception("Nedelec FESpace2 not implemented for BBND elements");
    //  (*testout) << "el = " << elnr << ", dnums = " << dnums << endl;
  }





  void NedelecFESpace2 ::GetTransformation (ELEMENT_TYPE eltype, 
					    int elnr,
					    const Array<int> & eorient,
					    const Array<int> & forient,
					    FlatVector<double> & fac) const
  {
    bool graddom = gradientdomains[ma->GetElIndex(ElementId(VOL,elnr))];

    fac = 1.0;
    switch (eltype)
      {
      case ET_SEGM:
	{
	  fac(0) = eorient[0];
	  if (order >= 3) fac(2) = eorient[0];
	  break;
	}

      case ET_TRIG:
	{
	  for (int i = 0; i < 3; i++)
	    {
	      fac(i) = eorient[i];
	      if (order >= 3) fac(i+6) = eorient[i];
	    }
	  break;
	}

      case ET_QUAD:
	{
	  int i, j, k, ii;

	  for (i = 0; i < 4; i++)
	    fac(i) = eorient[i];
	  ii = 4;

	  for (i = 0; i < 2; i++)
	    {
	      if (order >= 3)
		fac(ii+1) = eorient[i];
	      ii += order-1;
	    }

	  for (i = 0; i < 2; i++)
	    {
	      if (zorder >= 3)
		fac(ii+1) = eorient[i+2];
	      ii += zorder-1;
	    }

	
	  int nmx = order * (zorder-1);
	  // int nmy = (order-1) * zorder;

	  // vertical flip
	  if (forient[0] & 1)
	    { 
	      // horizontal moments
	      for (j = 0; j < order; j++)
		for (k = 1; k < zorder-1; k+=2)
		  fac(ii+j*(zorder-1)+k) *= -1;

	      // vertical moments
	      for (j = 0; j < order-1; j++)
		for (k = 0; k < zorder; k+=2)
		  fac(ii+nmx+j*zorder+k) *= -1;
	    }

	  // horizontal flip
	  if (forient[0] & 2)
	    {
	      // horizontal moments
	      for (j = 0; j < order; j+=2)
		for (k = 0; k < zorder-1; k++)
		  fac(ii+j*(zorder-1)+k) *= -1;

	      // vertical moments
	      for (j = 1; j < order-1; j+=2)
		for (k = 0; k < zorder; k++)
		  fac(ii+nmx+j*zorder+k) *= -1;
	    }
	  break;
	}
      
      case ET_TET:
	{
	  if (graddom)
	    for (int i = 0; i < 6; i++)
	      {
		fac(i) = eorient[i];
		if (order >= 3)
		  fac(i+12) = eorient[i];
	      }
	  else
	    for (int i = 0; i < 6; i++)
	      fac(i) = eorient[i];
	  break;
	}
      

      case ET_PRISM:
	{
	  int i, j, k, ii;
	  for (i = 0; i < 9; i++)
	    fac(i) = eorient[i];
	  ii = 9;

	  if (graddom)
	    {
	      for (i = 0; i < 6; i++)
		{
		  if (order >= 3)
		    fac(ii+1) = eorient[i];
		  ii += order-1;
		}
	    
	      for (i = 0; i < 3; i++)
		{
		  if (zorder >= 3)
		    fac(ii+1) = eorient[i+6];
		  ii += zorder-1;
		}
	    }

	  // trig faces:
	  if (order == 3)
	    ii += 6;


	
	  int nmx = order * (zorder-1);
	  // int nmy = (order-1) * zorder;

	  for (i = 2; i < 5; i++)
	    {
	      if (forient[i] & 1)
		{
		  for (j = 0; j < order; j++)
		    for (k = 1; k < zorder-1; k+=2)
		      fac(ii+j*(zorder-1)+k) *= -1;
		
		  for (j = 0; j < order-1; j++)
		    for (k = 0; k < zorder; k+=2)
		      fac(ii+nmx+j*zorder+k) *= -1;
		}
	    
	      if (forient[i] & 2)
		{
		  for (j = 0; j < order; j+=2)
		    for (k = 0; k < zorder-1; k++)
		      fac(ii+j*(zorder-1)+k) *= -1;
		
		  for (j = 1; j < order-1; j+=2)
		    for (k = 0; k < zorder; k++)
		      fac(ii+nmx+j*zorder+k) *= -1;
		}
	      ii += n_quad_face_dofs;
	    }
	  break;
	}



      case ET_PYRAMID:
	{
	  int i, ii;
	  for (i = 0; i < 8; i++)
	    {
	      fac(i) = eorient[i];
	      if (order >= 3)
		fac(i+16) = eorient[i];
	    }

	  ii = 8*order + 4 * n_trig_face_dofs;

	  // quad face:
	  if (order == 2 && zorder == 1)
	    {
	      if (forient[4] & 1)
		fac(ii) *= -1;
	      ii += n_quad_face_dofs;
	    }
	
	  if (order == 2 && zorder == 2)
	    {
	      if (forient[4] & 1)
		fac(ii+2) = -1;
	      if (forient[4] & 2)
		fac(ii) = -1;
	      ii += n_quad_face_dofs;
	    }

	  if (order == 2 && zorder == 3)
	    {
	      // x-components: 0,1,2,3
	      // y-components: 4,5,6
	      if (forient[4] & 1)
		{
		  fac(ii+1) *= -1;
		  fac(ii+3) *= -1;
		  fac(ii+4) *= -1;
		  fac(ii+6) *= -1;
		}
	      if (forient[4] & 2)
		{
		  fac(ii) *= -1;
		  fac(ii+1) *= -1;
		}
	      ii += n_quad_face_dofs;
	    }	

	
	  if (order == 3 && zorder == 1)
	    {
	      if (forient[4] & 1)
		{ // vertical flip
		  fac(ii) *= -1;
		  fac(ii+1) *= -1;
		}
	    
	      if (forient[4] & 2)
		{ // horizontal flip
		  fac(ii+1) *= -1;
		}
	      ii += n_quad_face_dofs;
	    }


	  if (order == 3 && zorder == 2)
	    {
	      if (forient[4] & 1)
		{ // vertical flip
		  fac(ii+3) *= -1;
		  fac(ii+5) *= -1;
		}
		
	      if (forient[4] & 2)
		{ // horizontal flip
		  fac(ii+0) *= -1;
		  fac(ii+2) *= -1;
		  fac(ii+4) *= -1;
		  fac(ii+5) *= -1;
		}
	      ii += n_quad_face_dofs;
	    }
	
	  if (order == 3 && zorder == 3)
	    {
	      if (forient[4] & 1)
		{ // vertical flip
		  fac(ii+1) *= -1;
		  fac(ii+3) *= -1;
		  fac(ii+5) *= -1;
		  fac(ii+6) *= -1;
		  fac(ii+8) *= -1;
		  fac(ii+9) *= -1;
		  fac(ii+11) *= -1;
		}
	    
	      if (forient[4] & 2)
		{ // horizontal flip
		  fac(ii+0) *= -1;
		  fac(ii+1) *= -1;
		  fac(ii+4) *= -1;
		  fac(ii+5) *= -1;
		  fac(ii+9) *= -1;
		  fac(ii+10) *= -1;
		  fac(ii+11) *= -1;
		}
	      ii += n_quad_face_dofs;
	    }
	  break;
	}

      default:
        {
          cerr << "unhandled case 152345" << endl;
          break;
        }

      }
  }




  template <class T>
  void NedelecFESpace2::TransformMat (ElementId ei, 
				      SliceMatrix<T> mat, TRANSFORM_TYPE tt) const
  {
    int nd;
    ELEMENT_TYPE et;
    ArrayMem<int,12> enums, eorient;
    ArrayMem<int,6> fnums, forient;
    LocalHeapMem<1000> lh("NedelecFESpace2 - TransformMat");

    bool boundary = ei.VB() == BND;
    size_t elnr = ei.Nr();
    if (boundary)
      {
	nd = GetFE (ei, lh).GetNDof();
	et = ma->GetElType (ei);
	ma->GetSElEdges (elnr, enums, eorient);
	ma->GetSElFace (elnr, fnums[0], forient[0]);
      }
    else
      {
	nd = GetFE (ei, lh).GetNDof();
	et = ma->GetElType (ei);
	ma->GetElEdges (elnr, enums, eorient);
	ma->GetElFaces (elnr, fnums, forient);
      }

  
    ArrayMem<double, 100> mem(nd);
    FlatVector<double> fac(nd, &mem[0]);

    GetTransformation (et, elnr, eorient, forient, fac);

    if (tt & TRANSFORM_MAT_LEFT)
      for (int k = 0; k < dimension; k++)
	for (int i = 0; i < nd; i++)
	  for (int j = 0; j < mat.Width(); j++)
	    mat(k+i*dimension, j) *= fac(i);
  
    if (tt & TRANSFORM_MAT_RIGHT)
      for (int l = 0; l < dimension; l++)
	for (int i = 0; i < mat.Height(); i++)
	  for (int j = 0; j < nd; j++)
	    mat(i, l+j*dimension) *= fac(j);
  }







  template <class T>
  void NedelecFESpace2::TransformVec (ElementId ei,
				      SliceVector<T> vec, TRANSFORM_TYPE tt) const
  {
    int nd;
    ELEMENT_TYPE et;
    /*
      int ena[12], eoa[12];
      int fna[12], foa[12];
      Array<int> enums(12, ena), eorient(12, eoa);
      Array<int> fnums(6, fna), forient(6, foa);
    */
    ArrayMem<int,12> enums, eorient;
    ArrayMem<int,6> fnums, forient;
    LocalHeapMem<1000> lh ("Nedelecfespace2, transformvec");
    size_t elnr = ei.Nr();
    if (ei.VB()==BND)
      {
	nd = GetFE (ei, lh).GetNDof();
	et = ma->GetElType (ei);
	ma->GetSElEdges (elnr, enums, eorient);
	ma->GetSElFace (elnr, fnums[0], forient[0]);
      }
    else
      {
	nd = GetFE (ei, lh).GetNDof();
	et = ma->GetElType (ei);
	ma->GetElEdges (elnr, enums, eorient);
	ma->GetElFaces (elnr, fnums, forient);
      }

    ArrayMem<double, 100> mem(nd);
    FlatVector<double> fac(nd, &mem[0]);

    GetTransformation (et, elnr, eorient, forient, fac);

    for (int k = 0; k < dimension; k++)
      for (int i = 0; i < nd; i++)
	vec(k+i*dimension) *= fac(i);
  }



  template
  void NedelecFESpace2::TransformVec<double> 
  (ElementId ei, SliceVector<double> vec, TRANSFORM_TYPE tt) const;
  template
  void NedelecFESpace2::TransformVec<Complex> 
  (ElementId ei, SliceVector<Complex> vec, TRANSFORM_TYPE tt) const;

  /*
  template
  void NedelecFESpace2::TransformMat<const FlatMatrix<double> > 
  (int elnr, bool boundary, const FlatMatrix<double> & mat, TRANSFORM_TYPE tt) const;
  template
  void NedelecFESpace2::TransformMat<const FlatMatrix<Complex> > 
  (int elnr, bool boundary, const FlatMatrix<Complex> & mat, TRANSFORM_TYPE tt) const;
  */
  
  template
  void NedelecFESpace2::TransformMat<double> 
  (ElementId ei, SliceMatrix<double> mat, TRANSFORM_TYPE tt) const;
  template
  void NedelecFESpace2::TransformMat<Complex> 
  (ElementId ei, SliceMatrix<Complex> mat, TRANSFORM_TYPE tt) const;








  void NedelecFESpace2 ::
  SetGradientDomains (const BitArray & adoms)
  {
    gradientdomains = adoms;
  }

  void NedelecFESpace2 ::
  SetGradientBoundaries (const BitArray & abnds)
  {
    gradientboundaries = abnds;
  }



  shared_ptr<Table<int>> NedelecFESpace2 ::
  CreateSmoothingBlocks (const Flags & precflags) const
  // CreateSmoothingBlocks (int type) const
  {
    int type = int(precflags.GetNumFlag("blocktype", NedelecFESpace::SB_AFW));
    cout << "Ned2, CreateSmoothingblocks, type = " << type << endl;
    int ne = ma->GetNE();
    // int nse = ma->GetNSE();
    int nv = ma->GetNV();
    // int nd = nv+ned + nfa + ne;
    // int level = ma->GetNLevels()-1;

    int i, j, k, l;

    ArrayMem<int,12> pnums;
    ArrayMem<int,37> dnums, dcluster;
    ArrayMem<int,12> eorient, fnums, forient, fpnum1, fpnum2;
    ArrayMem<int,12> ecluster, fcluster;

    // int elcluster;
    Table<int> * it = NULL;
    // const NedelecFESpace & fe1 = 
    // dynamic_cast<const NedelecFESpace&> (*low_order_space);

  
  
    switch (type)
      {
     
      case NedelecFESpace::SB_AFW:
	{
	  cout << " Take Smoothing Block Type SB_AFW " << endl; 
	  // all edge-dofs in common with vertex-representant
	  // all vertically aligned edges and faces
	  // all vertically aligned faces and elements
	
	  // non-overlapping small blocks
	  Array<int> cnts(nv+ned+nfa+ne);
	  for (k = 1; k <= 2; k++)
	    {
	      if (k == 2)
		it = new Table<int>(cnts);
	      cnts = 0;
	      for (i = 0; i < ned; i++)
		{
		  int nd, cl;
		  int ecl = ma->GetClusterRepEdge(i);
		  if (ecl < 0) continue;
		  if (ecl < nv)
		    { // vertical edge -> assign to edge-block (not ecluster = vertex)
		      cl = i+nv;
		      nd = n_z_edge_dofs;
		    }
		  else
		    { // plane edge -> assign to edge-cluster
		      cl = ecl;
		      nd = n_plane_edge_dofs;
		    }
	
		  if (!gradientedge[i] && nd > 1)
		    nd = 1;
		
		  if (k == 1)
		    cnts[cl] += nd;
		  else
		    for (l = 0; l < nd; l++)
		      (*it)[cl][cnts[cl]++] = i + l * ned;
		}
	    
	      for (i = 0; i < nfa; i++)
		{
		  int cl = ma->GetClusterRepFace (i);
		  if (cl < 0) continue;

		  int nd = first_face_dof[i+1] - first_face_dof[i];
		  if (!gradientface[i])
		    {
		      if (nd == 3) nd = 2;
		      if (nd == 4) nd = 3;
		    }
		  // if (nd == 3 && !gradientface[i]) nd = 2;

		  if (k == 1)
		    cnts[cl] += nd;
		  else
		    for (l = 0; l < nd; l++)
		      (*it)[cl][cnts[cl]++] = first_face_dof[i] + l;

		
		}
	      for (i = 0; i < nel; i++)
		{
		  int cl = ma->GetClusterRepElement (i);
		  int nd = first_el_dof[i+1] - first_el_dof[i];
		
		  if (k == 1)
		    cnts[cl] += nd;
		  else
		    for (l = 0; l < nd; l++)
		      (*it)[cl][cnts[cl]++] = first_el_dof[i] + l;

		  // test ... 
		  // if (ma->GetElType(i) == ET_PYRAMID)
		  {
		    ma->GetElFaces (i, fnums, forient);
		    for (j = 0; j < forient.Size(); j++)
		      {		
			int fcl = ma->GetClusterRepFace (fnums[j]);
			// int fcl = nv+ned+fnums[j];
			if (fcl != cl)
			  {
			    if (k == 1)
			      cnts[fcl] += nd;
			    else
			      for (l = 0; l < nd; l++)
				(*it)[fcl][cnts[fcl]++] = first_el_dof[i] + l;
			  }
		      }
		  }
		}

	      for (i = 0; i < ned; i++)
		{
		  int ecl = ma->GetClusterRepEdge (i);
		  if (ecl < 0) continue;

		  // int pi1, pi2;
		  auto pnts = ma->GetEdgePNums (i);
		  int pi1 = ma->GetClusterRepVertex (pnts[0]);
		  int pi2 = ma->GetClusterRepVertex (pnts[1]);

		  int nd = (pi1 == pi2) 
		    ? n_z_edge_dofs
		    : n_plane_edge_dofs;

		  if (!gradientedge[i] && nd > 1)
		    nd = 1;

		  if (k == 1)
		    {
		      cnts[pi1] += nd;
		      if (pi1 != pi2)
			cnts[pi2] += nd;
		    }
		  else
		    {
		      for (l = 0; l < nd; l++)
			{
			  (*it)[pi1][cnts[pi1]++] = i + l * ned;
			  if (pi1 != pi2)
			    (*it)[pi2][cnts[pi2]++] = i + l * ned;
			}
		    }
		}


	      for (i = 0; i < nfa; i++)
		{
		  int cl = ma->GetClusterRepFace (i);
		  if (cl < 0) continue;
		
		  // ma->GetFaceEdges (i, enums);
                  auto enums = ma->GetFaceEdges (i);

		  int nd = first_face_dof[i+1] - first_face_dof[i];
		  if (nd == n_quad_face_dofs) continue;
		  // if (nd == 3 && !gradientface[i]) nd = 2;
		  if (!gradientface[i])
		    {
		      if (nd == 3) nd = 2;
		      if (nd == 4) nd = 3;
		    }

		  for (j = 0; j < enums.Size(); j++)
		    {
		      //		    int pi = nv + enums[j];
		      int ecl = ma->GetClusterRepEdge (enums[j]);

		      if (ecl == cl) continue;


		      /*
		      // edges on face cluster
		      int pi1, pi2;
		      ma->GetEdgePNums (enums[j], pi1, pi2);
		      pi1 = ma->GetClusterRepVertex (pi1);
		      pi2 = ma->GetClusterRepVertex (pi2);

		      int nedof;
		      if (pi1 == pi2)
		      nedof = n_z_edge_dofs;
		      else
		      nedof = n_plane_edge_dofs;

		      if (k == 1)
		      cnts[cl] += nedof;
		      else
		      for (l = 0; l < nedof; l++)
		      (*it)[cl][cnts[cl]++] = enums[j] + l * ned;
		      */

		      // face on eclusters:
		      if (k == 1)
			cnts[ecl] += nd;
		      else 
			for (l = 0; l < nd; l++)
			  (*it)[ecl][cnts[ecl]++] = first_face_dof[i] + l;
		    }
		}


	      /*	    
		BitArray prism_faces(nfa);
		prism_faces.Clear();
		for (i = 0; i < -nel; i++)
		{
		ma->GetElFaces (i, fnums, forient);
		if (fnums.Size() == 4)
		for (j = 0; j < 4; j++)
		prism_faces.Set(fnums[j]);
		}
		for (i = 0; i < -nse; i++)
		{
		int fnr, fori;
		ma->GetSElFace (i, fnr, fori);
		prism_faces.Set(fnr);
		}

		for (i = 0; i < nfa; i++)
		{
		if (!prism_faces[i]) continue;

		int cl = ma->GetClusterRepFace (i);
		if (cl < 0) continue;
		
		ma->GetFacePNums (i, pnums);
		int nd = first_face_dof[i+1] - first_face_dof[i];
		
		if (pnums.Size() == 4)
		{ // quad... use diagonal
		pnums[1] = pnums[2];
		pnums.SetSize(2);
		}

		for (j = 0; j < pnums.Size(); j++)
		{
		int pi = ma->GetClusterRepVertex (pnums[j]);
		if (k == 1)
		cnts[pi] += nd;
		else
		for (l = 0; l < nd; l++)
		(*it)[pi][cnts[pi]++] = first_face_dof[i] + l;
		}
		}
	      */
	    
	      /*
		Array<int> fpnums;
		for (i = 0; i < nfa; i++)
		{
		ma->GetFacePNums (i, fpnums);
		for (j = 0; j < fpnums.Size(); j++)
		fpnums[j] = ma->GetClusterRepVertex (fpnums[j]);
	      
		for (j = 0; j < fpnums.Size(); j++)
		{
		bool dup = 0;
		for (l = 0; l < j; l++)
		if (fpnums[j] == fpnums[l])
		dup = 1;

		if (!dup)
		{
		int nd = first_face_dof[i+1] - first_face_dof[i];
		int pi = fpnums[j];
		if (k == 1)
		{
		cnts[pi] += nd;
		}
		else
		{
		for (l = 0; l < nd; l++)
		{
		(*it)[pi][cnts[pi]] = first_face_dof[i] + l;
		cnts[pi]++;
		}
		}
		}
		}
		}
	      */
	    }
	
	  //(*testout) << "AFW Blocks: " << (*it) << endl;
	  break;
	}


      case  NedelecFESpace::SB_HIPTMAIR:
	{
	  // all edge-dofs in common with vertex-representant
	  // all vertically aligned edges and faces
	  // all vertically aligned faces and elements
	
	  // non-overlapping small blocks
	  Array<int> cnts(nv+ned+nfa+ne);
	  cout << " Take Smoothing Block Type SB_HIPTMAIR " << endl; 
	  for (k = 1; k <= 2; k++)
	    {
	      if (k == 2)
		it = new Table<int>(cnts);
	      cnts = 0;
	    
	      for (i = 0; i < ned; i++)
		{
		  int nd, cl;
		  int ecl = ma->GetClusterRepEdge(i);
		  if (ecl < 0) continue;
		  if (ecl < nv)
		    {
		      cl = i+nv;
		      nd = n_z_edge_dofs;
		    }
		  else
		    {
		      cl = ecl;
		      nd = n_plane_edge_dofs;
		    }
		
		  if (k == 1)
		    {
		      cnts[cl] += nd;
		    }
		  else
		    {
		      for (l = 0; l < nd; l++)
			(*it)[cl][cnts[cl]++] = i + l * ned;
		    }
		}
	    
	      for (i = 0; i < nfa; i++)
		{
		  int cl = ma->GetClusterRepFace (i);
		  if (cl < 0) continue;

		  int nd = first_face_dof[i+1] - first_face_dof[i];

		  if (k == 1)
		    {
		      cnts[cl] += nd;
		    }
		  else
		    {
		      for (l = 0; l < nd; l++)
			(*it)[cl][cnts[cl]++] = first_face_dof[i] + l;
		    }
		}
	    
	      for (i = 0; i < nel; i++)
		{
		  int cl = ma->GetClusterRepElement (i);
		  int nd = first_el_dof[i+1] - first_el_dof[i];
		
		  if (k == 1)
		    {
		      cnts[cl] += nd;
		    }
		  else
		    {
		      for (l = 0; l < nd; l++)
			(*it)[cl][cnts[cl]++] = first_el_dof[i] + l;
		    }
		}
	    }



	  break;
	}


      case  NedelecFESpace::SB_POTENTIAL:
	{
	  Array<int> cnts(nv+ned);
	  for (int k = 1; k <= 2; k++)
	    {
	      if (k == 2)
		it = new Table<int>(cnts);
	    
	      cnts = 0;
	    
	      for (int j = 0; j < nv; j++)
		{
		  int vcl = ma->GetClusterRepVertex (j);
		  if (k == 2)
		    {
		      (*it)[vcl][cnts[vcl]++] = j;
		    }
		  else
		    {
		      cnts[vcl]++;
		    }
		}


	      for (i = 0; i < ned; i++)
		{
		  int nd, cl;
		  int ecl = ma->GetClusterRepEdge(i);
		  if (ecl < 0) continue;

		  if (ecl < nv)
		    { // vertical edge
		      cl = i+nv;
		      nd = n_z_edge_dofs-1;
		    }
		  else
		    {
		      cl = ecl;
		      nd = n_plane_edge_dofs-1;
		    }

		  if (k == 1)
		    {
		      cnts[cl] += nd;
		    }
		  else
		    {
		      for (l = 0; l < nd; l++)
			{
			  (*it)[cl][cnts[cl]] = i + l * ned + nv;
			  cnts[cl]++;
			}
		    }
		}
	    }
	  break;
	}
      }
  
    // (*testout) << "Nedelec2, Smoothingblocks type = " << type << endl << (*it) << endl;
    return shared_ptr<Table<int>> (it);
  }

  BitArray * NedelecFESpace2 :: 
  CreateIntermediatePlanes (int type) const
  {
    BitArray & ba = *new BitArray (GetNDof());
    ba.Clear();

    int i;
    for (i = 0; i < ned; i++)
      {
	auto pnts = ma->GetEdgePNums (i);
	  
	if (ma->GetClusterRepVertex (pnts[0]) ==
	    ma->GetClusterRepVertex (pnts[1]))
	  {
	    for (int l = 1; l < n_z_edge_dofs; l++)
	      ba.SetBit (i + l * ned);
	  }
      }


    for (i = 0; i < nfa; i++)
      {
	int first = first_face_dof[i];
	int nd = first_face_dof[i+1] - first;
	if (nd == n_quad_face_dofs)
	  {

	    if (order == 2 && zorder == 1)
	      {
		ba.SetBit (first);
	      }
	  
	    if (order == 2 && zorder == 2)
	      {
		ba.SetBit (first + 0);
		ba.SetBit (first + 2);
		ba.SetBit (first + 3);
	      }

	    if (order == 2 && zorder == 3)
	      {
		ba.SetBit (first + 0);
		ba.SetBit (first + 1);
		ba.SetBit (first + 5);
		ba.SetBit (first + 6);
	      }

	    if (order == 3 && zorder == 1)
	      {
		ba.SetBit (first + 0);
		ba.SetBit (first + 1);
	      }

	    if (order == 3 && zorder == 2)
	      {
		ba.SetBit (first + 0);
		ba.SetBit (first + 3);
		ba.SetBit (first + 4);

		/*
		  for (int l = 0; l < 7; l++)
		  ba.Set (first+l);
		*/
	      }

	    if (order == 3 && zorder == 3)
	      {
		ba.SetBit (first + 0);
		ba.SetBit (first + 1);
		ba.SetBit (first + 6);
		ba.SetBit (first + 7);
		ba.SetBit (first + 8);
	      }
	  }
      }

    /*
      for (i = first_el_dof[0]; i < first_el_dof[nel]; i++)
      ba.Set(i);
    */
    return &ba;
  }


  shared_ptr<Array<int>>
  NedelecFESpace2 :: CreateDirectSolverClusters (const Flags & flags) const
  {
    (*testout) << "CreateDirectSolverClusters" << endl;

    // int nv = ma->GetNV();
    int nd = GetNDof();
    int ne = ma->GetNE();


    auto spclusters = make_shared<Array<int>> (nd);
    Array<int> & clusters = *spclusters;
    clusters = 0;


    //lo 
    //for(i=0;i<ma->GetNEdges();i++)
    //  clusters[i]=1; 

    //
    for (size_t i = 0; i < ned; i++)
      {
	auto pts = ma->GetEdgePNums (i);
	  
	if (ma->GetClusterRepVertex (pts[0]) ==
	    ma->GetClusterRepVertex (pts[1]))
	  {
	    for (int l = 1; l < n_z_edge_dofs; l++)
	      clusters[i + l * ned] = 1;
	  }
      }

  
    for (size_t i = 0; i < nfa; i++)
      {
	int first = first_face_dof[i];
	int nd = first_face_dof[i+1] - first;
	if (nd == n_quad_face_dofs)
	  {

	    if (order == 2 && zorder == 1)
	      {
		clusters[first] = 1;
	      }
	  
	    if (order == 2 && zorder == 2)
	      {
		clusters[first + 0] = 1;
		clusters[first + 2] = 1;
		clusters[first + 3] = 1;
	      }

	    if (order == 2 && zorder == 3)
	      {
		clusters[first + 0] = 1;
		clusters[first + 1] = 1;
		clusters[first + 5] = 1;
		clusters[first + 6] = 1;
	      }

	    if (order == 3 && zorder == 1)
	      {
		clusters[first + 0] = 1;
		clusters[first + 1] = 1;
	      }

	    if (order == 3 && zorder == 2)
	      {
		clusters[first + 0] = 1;
		clusters[first + 3] = 1;
		clusters[first + 4] = 1;

		/*
		  for (int l = 0; l < 7; l++)
		  ba.Set (first+l);
		*/
	      }

	    if (order == 3 && zorder == 3)
	      {
		clusters[first + 0] = 1;
		clusters[first + 1] = 1;
		clusters[first + 6] = 1;
		clusters[first + 7] = 1;
		clusters[first + 8] = 1;
	      }
	  }
      }


    //


    Array<int> dnums;
  

    for(size_t i=0; i<ne && (directsolverclustered.Size() > 0 || directsolvermaterials.Size() > 0); i++)
      {
        ElementId ei(VOL,i);
	if((directsolverclustered.Size() > 0 && directsolverclustered[ma->GetElIndex(ei)]) || 
	   directsolvermaterials.Contains(ma->GetMaterial(ei)))
	  {     
	    ELEMENT_TYPE eltype = ma->GetElType(ei);
	    if(eltype != ET_PRISM) continue;

	    GetDofNrs(i,dnums);
	    for(size_t k=0; k<dnums.Size(); k++)
	      if(dnums[k] >= 0) clusters[dnums[k]] = 2;

	  }
      }
    
    for(size_t i=0; i< adddirectsolverdofs.Size(); i++)
      {
	clusters[adddirectsolverdofs[i]] = 2;
      }


    //   int clusternum = 2;
    //   for(i=0; i<ne; i++)
    //     {
    //       ELEMENT_TYPE eltype = ma->GetElType(i);
    //       if(eltype == ET_PRISM)
    // 	{
    // 	  GetDofNrs(i,dnums);
    // 	  for(k=0; k<dnums.Size(); k++)
    // 	    {
    // 	      if(dnums[k] >= 0) clusters[dnums[k]] = clusternum;
    // 	    }
    // 	  //clusternum++;
    // 	}
    //     }
	  
	 
    return spclusters;
    

  }



  SparseMatrix<double> * 
  NedelecFESpace2 :: CreateGradient() const
  {
    cout << "update gradient, N2" << endl;
    int j;
    size_t nv = ma->GetNV();
    int level = ma->GetNLevels()-1;
    const NedelecFESpace & fe1 = 
      dynamic_cast<const NedelecFESpace&> (*low_order_space);

    Array<int> cnts(GetNDof());
    cnts = 0;
    for (size_t i = 0; i < ned; i++)
      {
	if (fe1.FineLevelOfEdge(i) == level)
	  {
	    cnts[i] = 2;
	    for (j = 1; j < n_edge_dofs; j++)
	      cnts[i+j*ned] = 1;
	  }
      }
  

    SparseMatrix<double> & grad = *new SparseMatrix<double>(cnts);

    for (size_t i = 0; i < ned; i++)
      {
	if (fe1.FineLevelOfEdge(i) < level) continue;
	auto pts = ma->GetEdgePNums (i);
	grad.CreatePosition (i, pts[0]);
	grad.CreatePosition (i, pts[1]);
      }

    for (size_t i = 0; i < ned; i++)
      {
	if (fe1.FineLevelOfEdge(i) < level) continue;
	auto pts = ma->GetEdgePNums (i);
	grad(i, pts[0]) = 1;
	grad(i, pts[1]) = -1;
      }


    for (size_t i = 0; i < ned; i++)
      {
	if (fe1.FineLevelOfEdge(i) == level)
	  {
	    for (j = 1; j < n_edge_dofs; j++)
	      grad.CreatePosition(i+j*ned, i+(j-1)*ned+nv);
	  }
      }
    for (size_t i = 0; i < ned; i++)
      {
	if (fe1.FineLevelOfEdge(i) == level)
	  {
	    for (j = 1; j < n_edge_dofs; j++)
	      grad(i+j*ned, i+(j-1)*ned+nv) = 1;
	  }
      }

    (*testout) << "grad, p2 = " << grad << endl;
    return &grad;
  }





  void NedelecFESpace2 :: 
  LockSomeDofs (BaseMatrix & mat) const
  {
  
    cout << "Lock hanging dofs" << endl;

    // int eled, elfa;
    int i, j, k;
    int ne = ma->GetNE();
    // eled = 8;
    // elfa = 5;


    Matrix<double> elmat(1);
    Matrix<Complex> elmatc(1);
    elmat(0,0) = 1e15;
    elmatc(0,0) = 1e15;
    Array<int> dnums(1);

    Array<int> fnums, forient;
    Array<int> lock;

    cout << "type is " << typeid(mat).name() << endl;
    SparseMatrixSymmetric<Mat<1,1,double> > & smat =
      dynamic_cast<SparseMatrixSymmetric<Mat<1,1,double> > &> (mat);
  
    for (i = 0; i < ne; i++)
      {
        ElementId ei(VOL, i);
	lock.SetSize (0);
	switch (ma->GetElType(ei))
	  {
	  case ET_PRISM:
	    {
	      ma->GetElFaces (i, fnums, forient);
	      auto enums = ma->GetElEdges (ei);
	    
	      if (order == 3)
		{
		  // lock 3rd vert. edge dofs and trig face dofs
		  for (j = 6; j < 9; j++)
		    {
		      lock.Append (3 * enums[j]);
		    }
		  for (j = 0; j < 2; j++)
		    {
		      int base = first_face_dof[fnums[j]];
		      for (k = 0; k < n_trig_face_dofs; k++)
			lock.Append (base+k);
		    }
		}
	      break;
	    }
	  default:
	    { 
	      ;
	    }
	  }
      
	for (k = 0; k < lock.Size(); k++)
	  {
	    smat(lock[k], lock[k]) (0,0) += 1e15;
	  }
      }



    /*
      if (ma->GetElType(elnr) == ET_PYRAMID)
      {

      switch (type)
      {
      case N2:
      {
      for (j = 1; j <= eled; j++)
      {
      //		dnums.Elem(1) = abs (elementedges.Get(elnr)[j-1]) + ned;
      dnums.Elem(1) = enums.Elem(j) + ned;
      mat.AddElementMatrix (dnums, elmat);
      }
      for (j = 1; j <= elfa; j++)
      {
      dnums.Elem(1) =  2*ned + 2*fnums.Get(j) - 1;
      //		dnums.Elem(1) =  2*ned + 2*elementfaces.Get(elnr)[j-1] - 1;
      mat.AddElementMatrix (dnums, elmat);
      dnums.Elem(1)++;
      mat.AddElementMatrix (dnums, elmat);
      }  
      break;
      }
      case BDM1:
      {
      for (j = 1; j <= elfa; j++)
      {
      dnums.Elem(1) =  2*ned + fnums.Get(j);
      mat.AddElementMatrix (dnums, elmat);
      }  
      break;
      }
      case BDM2:
      {
      for (j = 1; j <= eled; j++)
      {
      //		int edge = abs (elementedges.Get(elnr)[j-1]);
      int edge = enums.Get(j);
      dnums.Elem(1) = edge + 2 * ned;
      mat.AddElementMatrix (dnums, elmat);
      }
      for (j = 1; j <= elfa; j++)
      {
      dnums.Elem(1) =  3*ned + 3*fnums.Get(j) - 2;
      //		dnums.Elem(1) =  3*ned + 3*elementfaces.Get(elnr)[j-1] - 2;
      mat.AddElementMatrix (dnums, elmat);
      dnums.Elem(1)++;
      mat.AddElementMatrix (dnums, elmat);
      dnums.Elem(1)++;
      mat.AddElementMatrix (dnums, elmat);
      }  
      break;
      }
      }
      }
    */
  }


  void NedelecFESpace2 :: GetVertexDofNrs (int vnr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
  }

  void NedelecFESpace2 :: GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  {
    cout << "EdgeDofs vom Nedelec2 space: SABINE FRAGEN.." << endl;
    dnums.SetSize(0);
  }

  void NedelecFESpace2 :: GetFaceDofNrs (int fanr, Array<int> & dnums) const
  {
    cout << "FaceDofs vom Nedelec2 space: SABINE FRAGEN.." << endl;
    dnums.SetSize(0);
  }

  void NedelecFESpace2 :: GetInnerDofNrs (int elnr, Array<int> & dnums) const
  {
    cout << "InnerDofs vom Nedelec2 space: SABINE FRAGEN.." << endl;
    dnums.SetSize(0);
  }

  // void NedelecFESpace2 :: 
  // AddGradient (double fac, const BaseVector & pot, BaseVector & grad) const
  // {
  //   shared_ptr<MeshAccess> ma = GetMeshAccess();

  //   int ned = ma->GetNEdges();
  
  //   const BaseSystemVector & svpot = 
  //     dynamic_cast<const BaseSystemVector&> (pot);
  //   BaseSystemVector & svgrad = 
  //     dynamic_cast<BaseSystemVector&> (grad);

  //   int sdim = svpot.SystemDim();
  //   if (sdim != 1)
  //     {
  //       cerr << "NedelecFESpace2::AddGradient not implemented for sdim != 1" << endl;
  //       exit(1);
  //     }

  //   int i, j;

  //   const SystemVector<SysVector1d> & svpot1 = 
  //     dynamic_cast<const SystemVector<SysVector1d>&> (pot);
  //   SystemVector<SysVector1d> & svgrad1 = 
  //     dynamic_cast<SystemVector<SysVector1d>&> (grad);

  //   BitArray usededges(ned);
  
  
  //   for (i = 1; i <= ned; i++)
  //     if (space.FineLevelOfEdge (i) >= level)
  //       {
  // 	int ep1 = space.EdgePoint1(i);
  // 	int ep2 = space.EdgePoint2(i);
	
  // 	for (j = 1; j <= sdim; j++)
  // 	  svgrad1.Elem(i, j) += fac * (svpot1.Get (ep1, j) - svpot1.Get(ep2, j));
  //       }

  // }
  // void NedelecFESpace2 ::  
  // ApplyGradientT (int level, const BaseVector & gradt, BaseVector & pott) const
  // {
  // }




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
