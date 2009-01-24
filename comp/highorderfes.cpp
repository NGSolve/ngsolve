#ifdef OLD_HIGHORDERFES

// the first high order fes in Ngsolve
// uses a nodal basis

/*********************************************************************/
/* File:   highorderfes.cc                                           */
/* Author: Joachim Schoeberl                                         */
/* Date:   28. Oct. 2000                                             */
/*********************************************************************/

/* 
   High Order Finite Element space
*/

#include <comp.hpp>

namespace ngcomp
{
  using namespace ngcomp;

  NodalFESpaceP ::   
  NodalFESpaceP (const MeshAccess & ama,
		 int ap, int adim, bool acomplex)
    : FESpace (ama, ap, adim, acomplex)
  {
    augmented = 0;
    p = ap;

    Flags loflags;
    loflags.SetFlag ("order", 1);
    loflags.SetFlag ("dim", dimension);
    if (iscomplex) loflags.SetFlag ("complex");
    low_order_space = new NodalFESpace (ma, loflags);
    // low_order_space = new NodalFESpace (ama, 1, adim, acomplex);

    /*
    if (augmented)
      {
	segm = new FE_Augmented_SegmP(p);
	trig = new FE_Augmented_TrigP(p);
	tet  = new FE_Augmented_TetP(p);
      }
    else
    */
      {
	segm = new FE_SegmP(p);
	trig = new FE_TrigP(p);
	tet  = new FE_TetP(p);
      }

#ifdef NONE
    quad = new FE_QuadP(p);
    prism = new FE_PrismP(p);
    hex = new FE_HexP(p);
#endif

    n_edge_dofs = p-1;
    // with quads
    //  n_face_dofs = (p-1)*(p-1);
    // trig only
    n_face_dofs = ((p-1)*(p-2)) / 2;

    if (ma.GetDimension() == 2)
      n_el_dofs = sqr (p-1);
    else
      n_el_dofs = ((p-4)*(p-4)*(p-4) + 6*(p-4)*(p-4)+11*(p-4)+6) / 6;


    if (ma.GetDimension() == 2)
      {
	evaluator = 
	  new MassIntegrator<2> (new ConstantCoefficientFunction(1));
	boundary_evaluator = 0;
      }
    else
      {
	evaluator = 
	  new MassIntegrator<3> (new ConstantCoefficientFunction(1));
	boundary_evaluator = 
	  new RobinIntegrator<3> (new ConstantCoefficientFunction(1));
      }

    if (dimension > 1)
      evaluator = new BlockBilinearFormIntegrator (*evaluator, dimension);      

  }


  NodalFESpaceP :: ~NodalFESpaceP ()
  {
    ;
  }

  FESpace * NodalFESpaceP :: Create (const MeshAccess & ama, const Flags & flags)
  {
    int aorder = int(flags.GetNumFlag ("order", 2));
    int adim = int(flags.GetNumFlag ("dim", 1));
    bool aiscomplex = flags.GetDefineFlag ("complex");

    return new NodalFESpaceP (ama, aorder, adim, aiscomplex);
  }




  void NodalFESpaceP :: Update()
  {
    if (low_order_space)
      low_order_space -> Update();

    nv = ma.GetNV();
    nedge = ma.GetNEdges();
    nface = ma.GetNFaces();
    ne = ma.GetNE();

    if (ma.GetDimension() == 2)
      ndof = nv +
	n_edge_dofs * nedge +
	n_el_dofs * ne;
    else
      ndof = nv + 
	n_edge_dofs * nedge + 
	n_face_dofs * nface +
	n_el_dofs * ne;

    if (augmented)
      ndof += (p-1) * nv;

#ifdef NONE
    delete edgedofs;
    delete facedofs;

    //  static int nfdoftab[] = { 0, 0, 0, 1, 3, 6, 10 };
    //  static int nidoftab[] = { 0, 0, 0, 0, 1, 4, 10 };

    const MeshAccess & ma = GetMeshAccess();


    int nedof = p - 1;

    int pp = p-3;
    int nfdof = (pp*pp+3*pp+2)/2;
    pp = p-4;
    int nidof =  (pp*pp*pp + 6 * pp * pp + 11 * pp + 6) / 6;

    int nquadfdof = (p-1)*(p-1);
    int nprismidof = nfdof * (p-1);
    int nhexidof = (p-1)*(p-1)*(p-1);


    int i, j, k, l;
    Array<int> pnums;

    ndof = ma.GetNP();

    edgedofs = new INT_2_HASHTABLE<int> (ndof+1);
    facedofs = new INT_3_HASHTABLE<int> (ndof+1);

    int ne = ma.GetNE();
    for (i = 1; i <= ne; i++)
      {
	ELEMENT_TYPE type = ma.GetElType (i);
	ma.GetElPNums (i, pnums);

	// define edges, trig-faces and quad-faces:
	int tetedges[][2] = 
	  { { 1, 2 }, { 1, 3 }, { 1, 4 }, { 2, 3 }, { 2, 4 }, { 3, 4 }, 
	    { 0, 0 } };
	int tettrigs[][3] = 
	  { { 1, 2, 3 }, { 1, 2, 4 }, { 1, 3, 4 }, { 2, 3, 4 },
	    { 0, 0, 0 } };
	int tetquads[][4] =
	  { { 0, 0, 0, 0 } };

	int prismedges[][2] = 
	  { { 1, 2 }, { 1, 3 }, { 2, 3 }, 
	    { 4, 5 }, { 4, 6 }, { 5, 6 },
	    { 1, 4 }, { 2, 5 }, { 3, 6 }, 
	    { 0, 0 } };
	int prismtrigs[][3] = 
	  { { 1, 2, 3 }, { 4, 5, 6 }, { 0, 0, 0 } };
	int prismquads[][4] =
	  { { 1, 2, 5, 4 }, { 2, 3, 6, 5 }, { 3, 1, 4, 6 }, { 0, 0, 0, 0 } };

	int hexedges[][2] =
	  { { 1, 2 }, { 2, 3 }, { 3, 4 }, { 4, 1 }, 
	    { 5, 6 }, { 6, 7 }, { 7, 8 }, { 8, 5 }, 
	    { 1, 5 }, { 2, 6 }, { 3, 7 }, { 4, 8 },
	    { 0, 0 } };
	int textrigs[][3] = 
	  { { 0, 0, 0 } };
	int texquads[][4] = 
	  { { 1, 2, 3, 4 }, { 5, 6, 7, 8 },
	    { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 4, 8, 7 }, { 4, 1, 5, 8 },
	    { 0, 0, 0, 0} };


	int (*ep)[2];
	int (*tfp)[3];
	int (*qfp)[4];
	switch (type)
	  {
	  case ET_TET:
	    ep = tetedges; 
	    tfp = tettrigs;
	    qfp = tetquads;
	    break;
	  case ET_PRISM:
	    ep = prismedges; 
	    tfp = prismtrigs;
	    qfp = prismquads;
	    break;
	  case ET_HEX:
	    ep = hexedges; 
	    tfp = textrigs;
	    qfp = texquads;
	    break;
	  }
	j = 0;
	while (ep[j][0])
	  {
	    INT_2 edge(pnums.Get(ep[j][0]), pnums.Get(ep[j][1]));
	    edge.Sort();
	    if (!edgedofs->Used (edge))
	      {
		edgedofs->Set (edge, ndof+1);
		ndof += nedof;
	      }
	    j++;
	  }

	/*
	  for (j = 1; j <= 3; j++)
	  for (k = j+1; k <= 4; k++)
	  {
	  INT_2 edge(pnums.Get(j), pnums.Get(k));
	  edge.Sort();
	  if (!edgedofs->Used (edge))
	  {
	  edgedofs->Set (edge, ndof+1);
	  ndof += nedof;
	  }
	  }
	*/



	j = 0;
	while (tfp[j][0])
	  {
	    INT_3 tface(pnums.Get(tfp[j][0]), 
			pnums.Get(tfp[j][1]),
			pnums.Get(tfp[j][2]));
	    tface.Sort();
	    if (!facedofs->Used (tface))
	      {
		facedofs->Set (tface, ndof+1);
		ndof += nfdof;
	      }
	    j++;
	  }


	j = 0;
	while (qfp[j][0])
	  {
	    INT_4Q qface(pnums.Get(qfp[j][0]), 
			 pnums.Get(qfp[j][1]),
			 pnums.Get(qfp[j][2]),
			 pnums.Get(qfp[j][3]));
	    qface.Sort();
	    INT_3 qf3 (qface.I1(), qface.I2(), qface.I3());
	    if (!facedofs->Used (qf3))
	      {
		facedofs->Set (qf3, ndof+1);
		ndof += nquadfdof;
	      }
	    j++;
	  }


	/*
	  for (j = 1; j <= 2; j++)
	  for (k = j+1; k <= 3; k++)
	  for (l = k+1; l <= 4; l++)
	  {
	  INT_3 face(pnums.Get(j), pnums.Get(k), pnums.Get(l));
	  face.Sort();
	  if (!facedofs->Used (face))
	  {
	  facedofs->Set (face, ndof+1);
	  ndof += nfdof;
	  }
	  }
	*/
      }

    eldofs.SetSize(ne);
    for (i = 1; i <= ne; i++)
      {
	eldofs.Elem(i) = ndof+1;

	ELEMENT_TYPE type = ma.GetElType (i);
	switch (type)
	  {
	  case ET_TET: ndof += nidof; break;
	  case ET_PRISM: ndof += nprismidof; break;
	  case ET_HEX: ndof += nhexidof; break;
	  }
      }
#endif
  }


  void NodalFESpaceP :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    Array<int> pnums;
    Array<int> ednums, edorient, fanums, faorient;

    ma.GetElPNums (elnr, pnums);
    ma.GetElEdges (elnr, ednums, edorient);
    if (ma.GetDimension() == 3)
      ma.GetElFaces (elnr, fanums, faorient);

    /*
      cout << endl << endl <<  "elnr = " << elnr << endl;
      cout << "pnums = " << pnums << endl;
      cout << "ednums = " << ednums << endl;
    */

    ELEMENT_TYPE type = ma.GetElType (elnr);

    int ndof;
    switch (type)
      {
      case ET_TRIG:  ndof = trig->GetNDof(); break;
      case ET_TET:   ndof = tet->GetNDof(); break;
      case ET_PRISM: ndof = prism->GetNDof(); break;
	//    case ET_HEX:   ndof = hex->GetNDof(); break;
      }
    dnums.SetSize (ndof);
  
    int i = 0, di;
    int j, k, l;
    int lami[4];

    switch (type)
      {
      case ET_TRIG:
	{
	  for (lami[0] = 0; lami[0] <= p; lami[0]++)
	    for (lami[1] = 0; lami[1] <= p-lami[0]; lami[1]++)
	      {
		lami[2] = p - lami[0] - lami[1];
		di = -1;
	      
		// vertex dofs
		for (j = 0; j < 3; j++)
		  if (lami[j] == p)
		    {
		      di = pnums[j];
		    }
	      
		// edge dofs
		if (di == -1)
		  {
		    for (j = 0; j < 2; j++)
		      for (k = j+1; k < 3; k++)
			if (lami[j] + lami[k] == p)
			  {
			    int ednum = 
			      ednums[ElementTopology::GetEdgeNr(ET_TRIG, j, k)];
			    di = GetEdgeDof (ednum, pnums[j], pnums[k],
					     lami[j], lami[k]);
			  }
		  }
	      
		// element dofs
		if (di == -1)
		  di = GetElementDof (elnr, lami[0], lami[1], lami[2], 0);
	      
		dnums[i] = di;
		i++;
	      }

	  if (augmented)
	    {
	      int base = nv +
		n_edge_dofs * nedge +
		n_el_dofs * ne;

	      if (augmented == 1)
		for (j = 0; j < 3; j++)
		  for (k = 0; k < p-1; k++)
		    {
		      dnums[i] = base + pnums[j] * (p-1) + k;
		      i++;
		    }
	      else
		for (j = 0; j < 3; j++)
		  {
		    dnums[i] = base + pnums[j] * (p-1) + p-2;
		    i++;
		    for (k = 1; k < p-1; k++)
		      {
			dnums[i] = -1;
			i++;
		      }
		  }
	    }

	  break;
	}


      case ET_TET:
	{
	  for (lami[0] = 0; lami[0] <= p; lami[0]++)
	    for (lami[1] = 0; lami[1] <= p-lami[0]; lami[1]++)
	      for (lami[2] = 0; lami[2] <= p-lami[0]-lami[1]; lami[2]++)
		{
		  lami[3] = p - lami[0] - lami[1] - lami[2];
		  di = -1;

		  // vertex dofs
		  for (j = 0; j < 4; j++)
		    if (lami[j] == p)
		      {
			di = pnums[j];
		      }

		  // edge dofs
		  if (di == -1)
		    {
		      for (j = 0; j < 3; j++)
			for (k = j+1; k < 4; k++)
			  if (lami[j] + lami[k] == p)
			    {
			      int ednum = 
				ednums[ElementTopology::GetEdgeNr(ET_TET, j, k)];
			      di = GetEdgeDof (ednum, pnums[j], pnums[k],
					       lami[j], lami[k]);
			    }
		    }

		  // face dofs
		  if (di == -1)
		    {
		      for (j = 0; j < 2; j++)
			for (k = j+1; k < 3; k++)
			  for (l = k+1; l < 4; l++)
			    if (lami[j] + lami[k] + lami[l] == p)
			      {
				int fanum = 
				  fanums[ElementTopology::GetFaceNr(ET_TET, j, k, l)];
				di = GetFaceDof (fanum, pnums[j], pnums[k], pnums[l],
						 lami[j], lami[k], lami[l]);
			      }
		    }

		  // element dofs
		  if (di == -1)
		    di = GetElementDof (elnr, lami[0], lami[1], lami[2], lami[3]);
		
		  dnums[i] = di;
		  i++;
		}

	  if (augmented)
	    {
	      int base = nv +
		n_edge_dofs * nedge +
		n_face_dofs * nface +
		n_el_dofs * ne;

	      if (augmented == 1)
		for (j = 0; j < 4; j++)
		  for (k = 0; k < p-1; k++)
		    {
		      dnums[i] = base + pnums[j] * (p-1) + k;
		      i++;
		    }
	      else
		for (j = 0; j < 4; j++)
		  {
		    dnums[i] = base + pnums[j] * (p-1) + p-2;
		    i++;
		    for (k = 1; k < p-1; k++)
		      {
			dnums[i] = -1;
			i++;
		      }
		  }
	      /*
	      for (j = 0; j < 4; j++)
		for (k = 0; k < p-1; k++)
		  {
		    dnums[i] = base + pnums[j] * (p-1) + k;
		    i++;
		  }
	      */
	    }

	  break;
	}
#ifdef NONE
      case ET_PRISM:
	{
	  Array<int> hpnums(6);
	  hpnums.Elem(1) = pnums.Get(3);
	  hpnums.Elem(2) = pnums.Get(2);
	  hpnums.Elem(3) = pnums.Get(1);
	  hpnums.Elem(4) = pnums.Get(6);
	  hpnums.Elem(5) = pnums.Get(5);
	  hpnums.Elem(6) = pnums.Get(4);

	  int l1, l2, l3, lz;
	  for (lz = 0; lz <= p; lz++)
	    for (l2 = 0; l2 <= p; l2++)
	      for (l1 = 0; l1 <= p-l2; l1++)
		{
		  l3 = p - l1 - l2;
		  //		(*testout) << "l1,2,3 = " << l1 << ", " << l2 << ", " << l3 << ", lz = " << lz << endl;
		  di = 0;
		  // corner dofs
		  if (lz == 0 && l3 == p) { di = hpnums.Get(1); }
		  if (lz == 0 && l1 == p) { di = hpnums.Get(2); }	
		  if (lz == 0 && l2 == p) { di = hpnums.Get(3); }
		  if (lz == p && l3 == p) { di = hpnums.Get(4); }
		  if (lz == p && l1 == p) { di = hpnums.Get(5); }	
		  if (lz == p && l2 == p) { di = hpnums.Get(6); }

		  if (!di)
		    { // edge dofs
		      if (lz == 0 && l2 == 0) di = GetEdgeDof (hpnums.Get(2), hpnums.Get(1), l1, p-l1);
		      if (lz == 0 && l1 == 0) di = GetEdgeDof (hpnums.Get(3), hpnums.Get(1), l2, p-l2);
		      if (lz == 0 && l3 == 0) di = GetEdgeDof (hpnums.Get(2), hpnums.Get(3), l1, p-l1);

		      if (lz == p && l2 == 0) di = GetEdgeDof (hpnums.Get(5), hpnums.Get(4), l1, p-l1);
		      if (lz == p && l1 == 0) di = GetEdgeDof (hpnums.Get(6), hpnums.Get(4), l2, p-l2);
		      if (lz == p && l3 == 0) di = GetEdgeDof (hpnums.Get(5), hpnums.Get(6), l1, p-l1);

		      if (l3 == p) di = GetEdgeDof (hpnums.Get(4), hpnums.Get(1), lz, p-lz);
		      if (l1 == p) di = GetEdgeDof (hpnums.Get(5), hpnums.Get(2), lz, p-lz);
		      if (l2 == p) di = GetEdgeDof (hpnums.Get(6), hpnums.Get(3), lz, p-lz);
		    }
		

		  if (!di)
		    { // face dofs

		      if (lz == 0)
			di = GetFaceDof (hpnums.Get(1), hpnums.Get(2), hpnums.Get(3), l3, l1, l2);
		      if (lz == p)
			di = GetFaceDof (hpnums.Get(4), hpnums.Get(5), hpnums.Get(6), l3, l1, l2);


		      if (l1 == 0)
			di = GetQuadFaceDof (hpnums.Get(1), hpnums.Get(3), hpnums.Get(6), hpnums.Get(4),
					     l2, lz);
		      if (l2 == 0)
			di = GetQuadFaceDof (hpnums.Get(1), hpnums.Get(2), hpnums.Get(5), hpnums.Get(4),
					     l1, lz);
		      if (l3 == 0)
			di = GetQuadFaceDof (hpnums.Get(2), hpnums.Get(3), hpnums.Get(6), hpnums.Get(5),
					     l2, lz);
		    }

		  if (!di)
		    { // internal dof
		      di = GetElementDof (elnr, l1, l2, l3, lz);
		    }
		    
		  i++;
		  dnums.Elem(i) = di;
		}
	  break;
	}




      case ET_HEX:
	{
	  int l1, l2, l3;
	  for (l3 = 0; l3 <= p; l3++)
	    for (l2 = 0; l2 <= p; l2++)
	      for (l1 = 0; l1 <= p; l1++)
		{
		  di = 0;
		  // corner dofs
		  if (l1 == 0 && l2 == 0 && l3 == 0) { di = pnums.Get(1); }
		  if (l1 == p && l2 == 0 && l3 == 0) { di = pnums.Get(2); }
		  if (l1 == p && l2 == p && l3 == 0) { di = pnums.Get(3); }
		  if (l1 == 0 && l2 == p && l3 == 0) { di = pnums.Get(4); }
		  if (l1 == 0 && l2 == 0 && l3 == p) { di = pnums.Get(5); }
		  if (l1 == p && l2 == 0 && l3 == p) { di = pnums.Get(6); }
		  if (l1 == p && l2 == p && l3 == p) { di = pnums.Get(7); }
		  if (l1 == 0 && l2 == p && l3 == p) { di = pnums.Get(8); }

		  if (!di)
		    { // edge dofs
		      if (l2 == 0 && l3 == 0) di = GetEdgeDof (pnums.Get(2), pnums.Get(1), l1, p-l1);
		      if (l2 == p && l3 == 0) di = GetEdgeDof (pnums.Get(3), pnums.Get(4), l1, p-l1);
		      if (l2 == 0 && l3 == p) di = GetEdgeDof (pnums.Get(6), pnums.Get(5), l1, p-l1);
		      if (l2 == p && l3 == p) di = GetEdgeDof (pnums.Get(7), pnums.Get(8), l1, p-l1);
		    
		      if (l1 == 0 && l3 == 0) di = GetEdgeDof (pnums.Get(4), pnums.Get(1), l2, p-l2);
		      if (l1 == p && l3 == 0) di = GetEdgeDof (pnums.Get(3), pnums.Get(2), l2, p-l2);
		      if (l1 == 0 && l3 == p) di = GetEdgeDof (pnums.Get(8), pnums.Get(5), l2, p-l2);
		      if (l1 == p && l3 == p) di = GetEdgeDof (pnums.Get(7), pnums.Get(6), l2, p-l2);
		    
		      if (l1 == 0 && l2 == 0) di = GetEdgeDof (pnums.Get(5), pnums.Get(1), l3, p-l3);
		      if (l1 == p && l2 == 0) di = GetEdgeDof (pnums.Get(6), pnums.Get(2), l3, p-l3);
		      if (l1 == 0 && l2 == p) di = GetEdgeDof (pnums.Get(8), pnums.Get(4), l3, p-l3);
		      if (l1 == p && l2 == p) di = GetEdgeDof (pnums.Get(7), pnums.Get(3), l3, p-l3);
		    }
		
		  if (!di)
		    { // face dofs
		      if (l1 == 0)
			di = GetQuadFaceDof (pnums.Get(1), pnums.Get(4), pnums.Get(8), pnums.Get(5),
					     l2, l3);
		      if (l1 == p)
			di = GetQuadFaceDof (pnums.Get(2), pnums.Get(3), pnums.Get(7), pnums.Get(6),
					     l2, l3);

		      if (l2 == 0)
			di = GetQuadFaceDof (pnums.Get(1), pnums.Get(2), pnums.Get(6), pnums.Get(5),
					     l1, l3);
		      if (l2 == p)
			di = GetQuadFaceDof (pnums.Get(4), pnums.Get(3), pnums.Get(7), pnums.Get(8),
					     l1, l3);

		      if (l3 == 0)
			di = GetQuadFaceDof (pnums.Get(1), pnums.Get(2), pnums.Get(3), pnums.Get(4),
					     l1, l2);
		      if (l3 == p)
			di = GetQuadFaceDof (pnums.Get(5), pnums.Get(6), pnums.Get(7), pnums.Get(8),
					     l1, l2);
		    }

		  if (!di)
		    { // internal dof
		      di = GetElementDof (elnr, l1, l2, l3, -1);
		    }
		    
		  i++;
		  dnums.Elem(i) = di;
		}
	  break;
	}
#endif
      }
    (*testout) << "el " << elnr << ", dnums = " << dnums << endl;
  }


  void NodalFESpaceP :: GetSDofNrs (int selnr, Array<int> & dnums) const
  {
    //  cout << "getsdofnrs" << endl;

    Array<int> pnums, ednums, edorient;
    int fanum, faorient;
    ma.GetSElPNums (selnr, pnums);
    ma.GetSElEdges (selnr, ednums, edorient);
    if (ma.GetDimension() == 3)
      ma.GetSElFace (selnr, fanum, faorient);


    /*
      cout << "pnums = " << pnums << endl;
      cout << "ednums = " << ednums << endl;
    */

    ELEMENT_TYPE type = ma.GetSElType (selnr);

    int ndof;
    switch (type)
      {
      case ET_SEGM: ndof = segm->GetNDof(); break;
      case ET_TRIG: ndof = trig->GetNDof(); break;
      case ET_QUAD: ndof = quad->GetNDof(); break;
      }


    dnums.SetSize (ndof);


    int i = 0, di;
    int j, k, l;
    int lami[3];

    switch (type)
      {
      case ET_SEGM:
	{
	  dnums[0] = pnums[0];
	  dnums[p] = pnums[1];
	
	  i = 1;
	  for (lami[0] = 1; lami[0] < p; lami[0]++)
	    {
	      lami[1] = p - lami[0];
	      dnums[i] = GetEdgeDof (ednums[0], pnums[0], pnums[1], lami[0], lami[1]);
	      i++;
	    }
	
	  i = p+1;
	  if (augmented)
	    {
	      int base = nv +
		n_edge_dofs * nedge +
		n_el_dofs * ne;
	      for (j = 0; j < 2; j++)
		for (k = 0; k < p-1; k++)
		  {
		    dnums[i] = base + pnums[j] * (p-1) + k;
		    i++;
		  }
	    }

	  break;
	}

      case ET_TRIG:
	{
	  for (lami[0] = 0; lami[0] <= p; lami[0]++)
	    for (lami[1] = 0; lami[1] <= p-lami[0]; lami[1]++)
	      {
		lami[2] = p - lami[0] - lami[1];
		di = -1;

		// vertex dof
		for (j = 0; j < 3; j++)
		  if (lami[j] == p)
		    {
		      di = pnums[j];
		    }
	      
		// edge dof
		if (di == -1)
		  {
		    for (j = 0; j < 2; j++)
		      for (k = j+1; k < 3; k++)
			if (lami[j] + lami[k] == p)
			  {
			    int ednum = 
			      ednums[ElementTopology::GetEdgeNr(ET_TRIG, j, k)];
			    di = GetEdgeDof (ednum, pnums[j], pnums[k],
					     lami[j], lami[k]);
			  }
		  }
	      
		if (di == -1)
		  {
		    di = GetFaceDof (fanum, pnums[0], pnums[1], pnums[2],
				     lami[0], lami[1], lami[2]);
		  }
		dnums[i] = di;
		i++;
	      }

	  if (augmented)
	    {
	      int base = nv +
		n_edge_dofs * nedge +
		n_face_dofs * nface +
		n_el_dofs * ne;
	    
	      if (augmented == 1)
		for (j = 0; j < 3; j++)
		  for (k = 0; k < p-1; k++)
		    {
		      dnums[i] = base + pnums[j] * (p-1) + k;
		      i++;
		    }
	      else
		for (j = 0; j < 3; j++)
		  {
		    dnums[i] = base + pnums[j] * (p-1) + p-2;
		    i++;
		    for (k = 1; k < p-1; k++)
		      {
			dnums[i] = -1;
			i++;
		      }
		  }
	      /*
	      for (j = 0; j < 3; j++)
		for (k = 0; k < p-1; k++)
		  {
		    dnums[i] = base + pnums[j] * (p-1) + k;
		    i++;
		  }
	      */
	    }

	  break;
	}


#ifdef NONE

      case ET_QUAD:
	{
	  for (lami[1] = 0; lami[1] <= p; lami[1]++)
	    for (lami[0] = 0; lami[0] <= p; lami[0]++)
	      {
		di = 0;

		// nodal dofs
		if (lami[0] == 0 && lami[1] == 0) di = pnums.Get(1);
		if (lami[0] == p && lami[1] == 0) di = pnums.Get(2);
		if (lami[0] == p && lami[1] == p) di = pnums.Get(3);
		if (lami[0] == 0 && lami[1] == p) di = pnums.Get(4);
	    
		if (!di)
		  { // edge dofs
		    if (lami[0] == 0) di = GetEdgeDof (pnums.Get(4), pnums.Get(1),
						       lami[1], p-lami[1]);
		    if (lami[0] == p) di = GetEdgeDof (pnums.Get(3), pnums.Get(2),
						       lami[1], p-lami[1]);
		    if (lami[1] == 0) di = GetEdgeDof (pnums.Get(2), pnums.Get(1),
						       lami[0], p-lami[0]);
		    if (lami[1] == p) di = GetEdgeDof (pnums.Get(3), pnums.Get(4),
						       lami[0], p-lami[0]);
		  }
	    
		if (!di)
		  di = GetQuadFaceDof (pnums.Get(1), pnums.Get(2), 
				       pnums.Get(3), pnums.Get(4),
				       lami[0], lami[1]);
	    
		i++;
		dnums.Elem(i) = di;
	      }
	}
#endif
      }

    //  cout << "surf el " << selnr << ", ednum = " << ednums << " dnums = " << dnums << endl;
    /*
      (*testout) << "surf, dnums = ";
      for (i = 1; i <= dnums.Size(); i++)
      (*testout) << dnums.Get(i) << " ";
      (*testout) << endl;
    */
  }

  int NodalFESpaceP :: GetEdgeDof (int ednr, int v1, int v2, int lam1, int lam2) const
  {
    int vi1, vi2;

    int base = nv + n_edge_dofs * ednr;
    ma.GetEdgePNums (ednr, vi1, vi2);

    if (vi1 == v1)
      base += lam1 - 1;
    else
      base += lam2 - 1;
    
    return base;
  }

  int NodalFESpaceP :: GetFaceDof (int fnr, int i1, int i2, int i3, int lam1, int lam2, int lam3) const
  {
    int base = nv + n_edge_dofs * nedge + n_face_dofs * fnr;
    Array<int> pnums;
    ma.GetFacePNums (fnr, pnums);

    int lr1, lr2, lr3;

    if (pnums[0] == i1) lr1 = lam1;
    else if (pnums[0] == i2) lr1 = lam2;
    else if (pnums[0] == i3) lr1 = lam3;
    else
      cerr << "GetFaceDof, errror" << endl;

    if (pnums[1] == i1) lr2 = lam1;
    else if (pnums[1] == i2) lr2 = lam2;
    else if (pnums[1] == i3) lr2 = lam3;
    else
      cerr << "GetFaceDof, errror 2" << endl;

    base += lr1-1 + (lr1+lr2-1) * (lr1+lr2-2) / 2;

    return base;

#ifdef NONE
    INT_3 face(i1, i2, i3);
    face.Sort();
    int base = facedofs->Get(face);
  
    int lr1, lr2;
    if (face.I1() == i1)
      lr1 = lam1;
    else if (face.I1() == i2)
      lr1 = lam2;
    else
      lr1 = lam3;

    if (face.I2() == i1)
      lr2 = lam1;
    else if (face.I2() == i2)
      lr2 = lam2;
    else
      lr2 = lam3;

    // lr2 = p - lr2 - 1;
    //  int loc = (lr1-1) + lr2 * (lr2-1)/2;
    int loc = lr1-1 + (lr2-1) * (p-1) - (lr2-1) * lr2 / 2;
    //  (*testout) << "lr1,2 = " << lr1 << ", " << lr2 << ", Loc = " << loc << endl;
    return base+loc;
#endif
  }
 


  int NodalFESpaceP :: GetQuadFaceDof (int fnr, int i1, int i2, int i3, int i4,
				       int lam1, int lam2) const
  {
    return 0;
#ifdef NONE
    INT_4Q face(i1, i2, i3, i4);
    //  (*testout) << "get quad face, " << face;

    if (min2 (face.I2(), face.I3()) < min2 (face.I1(), face.I4()))
      { 
	swap (face.I1(), face.I2()); swap (face.I3(), face.I4());
	lam1 = p - lam1;
      }
    if (face.I4() < face.I1())
      { 
	swap (face.I1(), face.I4()); swap (face.I2(), face.I3());
	lam2 = p - lam2;
      }
    if (face.I4() < face.I2())
      { 
	swap (face.I2(), face.I4()); 
	swap (lam1, lam2);
      }
    //  (*testout) << ", sorted, " << face;

    INT_3 f3 (face.I1(), face.I2(), face.I3());

    face.Sort();
    //  (*testout) << " == " << face << endl;

    int base = facedofs->Get(f3);
    int loc = lam1-1 + (lam2-1) * (p-1);
    return base+loc;
#endif
  }



  int NodalFESpaceP :: GetElementDof (int elnr, int lam1, int lam2, 
				      int lam3, int lam4) const
  {
    ELEMENT_TYPE type = ma.GetElType (elnr);
    int base;
    if (ma.GetDimension() == 2)
      base = nv + n_edge_dofs * nedge + 
	n_el_dofs * elnr;
    else
      base = nv + n_edge_dofs * nedge + 
	n_face_dofs * nface + n_el_dofs * elnr;
  
    int loc;

    loc = 0;
    switch (type)
      {
      case ET_TRIG:
	{
	  loc = lam1-1 + (lam1+lam2-1) * (lam1+lam2-2) / 2;
	  break;
	} 
      case ET_TET:
	{
	  lam1--;
	  lam2--;
	  lam3--;
	
	  loc = lam1 + lam2 * (p-3-lam3) - lam2*(lam2-1)/2;
	  int j;
	  for (j = 0; j < lam3; j++)
	    loc += (p-3-j) * (p-2-j) / 2;

	  //	loc = lam1-1 + (p-1) * (lam2-1) + (p-1)*(p-1) * (lam3-1);
	  break;
	}
      case ET_HEX:
	{
	  loc = (lam1-1) + (lam2-1)*(p-1) + (lam3-1) * (p-1) * (p-1);
	  break;
	}
      case ET_PRISM:
	{
	  int pp = p-3;
	  int nfdof = (pp*pp+3*pp+2)/2;

	  loc = (lam4-1) * nfdof +  lam1-1 + (lam2-1) * (p-1) - (lam2-1) * lam2 / 2; 
	  break;
	}
      default:
	{
	  cerr << "NodalFESpaceP::GetElementDof not implemented for el " << type << endl;
	}
      }
    return base+loc;
  }




  template <class MAT>
  void NodalFESpaceP :: TransformMat (int elnr, bool boundary,
				      MAT & mat, TRANSFORM_TYPE tt) const
  {
    ;
  }




  Table<int> * NodalFESpaceP :: CreateSmoothingBlocks (int type) const
  {
    int i, j, k, l;
    Array<int> eledges, orient, pnums;
    cout << "NodalFESpaceP::CreateSmoothingBlocks" << endl;

    if (ma.GetDimension() == 2)
      {
	// V,E,C, 2D
	Array<int> cnts(nv+nedge+ne+1);
	cnts = 0;

	if (!augmented)
	  for (i = 0; i < nv; i++)
	    cnts[i] = 1;
	else
	  for (i = 0; i < nv; i++)
	    cnts[i] = p;
      
	for (i = 0; i < nedge; i++)
	  cnts[nv+i] = n_edge_dofs;

	for (i = 0; i < ne; i++)
	  cnts[nv+nedge+i] = n_el_dofs;
      
	cnts[nv+nedge+ne] = 0; // nv;
      
	Table<int> * it = new Table<int> (cnts);
      
	cnts = 0;
      
	for (i = 0; i < nv; i++)
	  (*it)[i][0] = i;
      
	if (augmented)
	  {
	    int base = nv +
	      n_edge_dofs * nedge +
	      n_el_dofs * ne;
	  
	    for (i = 0; i < nv; i++)
	      for (j = 0; j < p-1; j++)
		(*it)[i][j+1] = base+(p-1)*i+j;
	  }
      
	k = nv;
	for (i = 0; i < nedge; i++)
	  {
	    for (j = 0; j < n_edge_dofs; j++)
	      {
		(*it)[nv+i][j] = k; 
		k++; 
	      }
	    cnts[nv+i] = n_edge_dofs;
	  }
	for (i = 0; i < ne; i++)
	  {
	    for (l = 0; l < n_el_dofs; l++)
	      {
		(*it)[nv+nedge+nface+i][l] = k;
		k++;
	      }
	  }

	for (i = 0; i < cnts[nv+nedge+ne]; i++)
	  (*it)[nv+nedge+ne][i] = i;
      
	(*testout) << "table: " << (*it) << endl;
	return it;
	/*
	// V,E+C, 2D
	Array<int> cnts(nv+nedge+1);
	cnts = 0;
      
	if (!augmented)
	for (i = 0; i < nv; i++)
	cnts[i] = 1;
	else
	for (i = 0; i < nv; i++)
	cnts[i] = p;
      
	for (i = 0; i < nedge; i++)
	cnts[nv+i] = n_edge_dofs;
	for (i = 0; i < ne; i++)
	{
	ma.GetElEdges (i, eledges, orient);
	for (j = 0; j < eledges.Size(); j++)
	cnts[nv+eledges[j]] += n_el_dofs;
	}
      
	cnts[nv+nedge] = nv;
      
	Table<int> * it = new Table<int> (cnts);
      
	cnts = 0;

	for (i = 0; i < nv; i++)
	(*it)[i][0] = i;

	if (augmented)
	{
	int base = nv +
	n_edge_dofs * nedge +
	n_el_dofs * ne;
	  
	for (i = 0; i < nv; i++)
	for (j = 0; j < p-1; j++)
	(*it)[i][j+1] = base+(p-1)*i+j;
	}
      
      
	k = nv;
	for (i = 0; i < nedge; i++)
	{
	for (j = 0; j < n_edge_dofs; j++)
	{
	(*it)[nv+i][j] = k; 
	k++; 
	}
	cnts[nv+i] = n_edge_dofs;
	}
      
	for (i = 0; i < ne; i++)
	{
	ma.GetElEdges (i, eledges, orient);
	for (int enr = 0; enr < eledges.Size(); enr++)
	for (j = 0; j < n_el_dofs; j++)
	{
	int cl = nv+eledges[enr];
	(*it)[cl][cnts[cl]] = k + j;
	cnts[cl]++;
	}
	k += n_el_dofs;
	}
      
	for (i = 0; i < nv; i++)
	(*it)[nv+nedge][i] = i;

	(*testout) << "table: " << (*it) << endl;
	return it;
	*/
      }
    else
      {
	/*
	// V,E,F,C, 3D
	Array<int> cnts(nv+nedge+nface+ne+1);
	cnts = 0;

	if (!augmented)
	for (i = 0; i < nv; i++)
	cnts[i] = 1;
	else
	for (i = 0; i < nv; i++)
	cnts[i] = p;
      
	for (i = 0; i < nedge; i++)
	cnts[nv+i] = n_edge_dofs;

	for (i = 0; i < nface; i++)
	cnts[nv+nedge+i] = n_face_dofs;

	for (i = 0; i < ne; i++)
	cnts[nv+nedge+nface+i] = n_el_dofs;
      
	cnts[nv+nedge+nface+ne] = nv;
      
	Table<int> * it = new Table<int> (cnts);
      
	cnts = 0;
      
	for (i = 0; i < nv; i++)
	(*it)[i][0] = i;
      
	if (augmented)
	{
	int base = nv +
	n_edge_dofs * nedge +
	n_face_dofs * nface +
	n_el_dofs * ne;
	  
	for (i = 0; i < nv; i++)
	for (j = 0; j < p-1; j++)
	(*it)[i][j+1] = base+(p-1)*i+j;
	}
      
      
	k = nv;
	for (i = 0; i < nedge; i++)
	{
	for (j = 0; j < n_edge_dofs; j++)
	{
	(*it)[nv+i][j] = k; 
	k++; 
	}
	cnts[nv+i] = n_edge_dofs;
	}
      
	for (i = 0; i < nface; i++)
	{
	for (l = 0; l < n_face_dofs; l++)
	{
	(*it)[nv+nedge+i][l] = k;
	k++;
	}
	}

	for (i = 0; i < ne; i++)
	{
	for (l = 0; l < n_el_dofs; l++)
	{
	(*it)[nv+nedge+nface+i][l] = k;
	k++;
	}
	}

	for (i = 0; i < nv; i++)
	(*it)[nv+nedge+nface+ne][i] = i;
      
	(*testout) << "table: " << (*it) << endl;
	return it;
	*/


	// V,E,F+C, 3D
	Array<int> cnts(nv+nedge+nface+1);
	cnts = 0;

	if (!augmented)
	for (i = 0; i < nv; i++)
	cnts[i] = 1;
	else
	for (i = 0; i < nv; i++)
	cnts[i] = p;
      
	for (i = 0; i < nedge; i++)
	cnts[nv+i] = n_edge_dofs;

	for (i = 0; i < nface; i++)
	cnts[nv+nedge+i] = n_face_dofs;

	for (i = 0; i < ne; i++)
	{
	ma.GetElFaces (i, eledges, orient);
	for (j = 0; j < eledges.Size(); j++)
	cnts[nv+nedge+eledges[j]] += n_el_dofs;
	}
      
	cnts[nv+nedge+nface] = nv;
      
	Table<int> * it = new Table<int> (cnts);
      
	cnts = 0;
      
	for (i = 0; i < nv; i++)
	(*it)[i][0] = i;
      
	if (augmented)
	{
	int base = nv +
	n_edge_dofs * nedge +
	n_face_dofs * nface +
	n_el_dofs * ne;
	  
	for (i = 0; i < nv; i++)
	for (j = 0; j < p-1; j++)
	(*it)[i][j+1] = base+(p-1)*i+j;
	}
      
      
	k = nv;
	for (i = 0; i < nedge; i++)
	{
	for (j = 0; j < n_edge_dofs; j++)
	{
	(*it)[nv+i][j] = k; 
	k++; 
	}
	cnts[nv+i] = n_edge_dofs;
	}
      
	for (i = 0; i < nface; i++)
	{
	for (l = 0; l < n_face_dofs; l++)
	{
	(*it)[nv+nedge+i][l] = k;
	k++;
	}
	}

	cnts = n_face_dofs;
	for (i = 0; i < ne; i++)
	{
	ma.GetElFaces (i, eledges, orient);
	for (int enr = 0; enr < eledges.Size(); enr++)
	for (j = 0; j < n_el_dofs; j++)
	{
	int cl = nv+nedge+eledges[enr];
	(*it)[cl][cnts[cl]] = k + j;
	cnts[cl]++;
	}
	k += n_el_dofs;
	}
      

	for (i = 0; i < nv; i++)
	(*it)[nv+nedge+nface][i] = i;
      
	(*testout) << "table: " << (*it) << endl;
	return it;


	/*
	// V,E+F+C, 3D
	HashTable<INT<2>, int> vert2edge(nedge+1);
	for (i = 0; i < nedge; i++)
	  {
	    INT<2> edge;
	    ma.GetEdgePNums (i, edge[0], edge[1]);
	    edge.Sort();
	    vert2edge.Set (edge, i);
	  }
      
	// V,E,V+C, 3D
	Array<int> cnts(nv+nedge+1);
	cnts = 0;

	if (augmented == 0)
	  for (i = 0; i < nv; i++)
	    cnts[i] = 1;
	else if (augmented == 1)
	  for (i = 0; i < nv; i++)
	    cnts[i] = p;
	else
	  for (i = 0; i < nv; i++)
	    cnts[i] = 2;
      
	for (i = 0; i < nedge; i++)
	  cnts[nv+i] = n_edge_dofs;
	for (i = 0; i < nface; i++)
	  {
	    ma.GetFacePNums (i, pnums);
	    for (j = 0; j < 3; j++)
	      {
		INT<2> edge;
		edge[0] = pnums[j];
		edge[1] = pnums[(j+1)%3];
		edge.Sort();
		int edgenr = vert2edge.Get(edge);
		cnts[nv+edgenr] += n_face_dofs;
	      }
	  }
	for (i = 0; i < ne; i++)
	  {
	    ma.GetElEdges (i, eledges, orient);
	    for (j = 0; j < eledges.Size(); j++)
	      cnts[nv+eledges[j]] += n_el_dofs;
	  }
      
	cnts[nv+nedge] = nv;
      
	Table<int> * it = new Table<int> (cnts);
      

	cnts = 0;
      
	for (i = 0; i < nv; i++)
	  (*it)[i][0] = i;

	if (augmented)
	  {
	    int base = nv +
	      n_edge_dofs * nedge +
	      n_face_dofs * nface +
	      n_el_dofs * ne;
	  
	    if (augmented == 1)
	      for (i = 0; i < nv; i++)
		for (j = 0; j < p-1; j++)
		  (*it)[i][j+1] = base+(p-1)*i+j;
	    else
	      for (i = 0; i < nv; i++)
		(*it)[i][1] = base+(p-1)*i + p-2;
	  }
      
      
	k = nv;
	for (i = 0; i < nedge; i++)
	  {
	    for (j = 0; j < n_edge_dofs; j++)
	      {
		(*it)[nv+i][j] = k; 
		k++; 
	      }
	    cnts[nv+i] = n_edge_dofs;
	  }
      
	for (i = 0; i < nface; i++)
	  {
	    ma.GetFacePNums (i, pnums);
	    for (j = 0; j < 3; j++)
	      {
		INT<2> edge;
		edge[0] = pnums[j];
		edge[1] = pnums[(j+1)%3];
		edge.Sort();
		int edgenr = vert2edge.Get(edge);

		for (l = 0; l < n_face_dofs; l++)
		  {
		    int cl = nv+edgenr;
		    (*it)[cl][cnts[cl]] = k+l;
		    cnts[cl]++;
		  }
	      }
	    k += n_face_dofs;
	  }

	for (i = 0; i < ne; i++)
	  {
	    ma.GetElEdges (i, eledges, orient);
	    for (int enr = 0; enr < eledges.Size(); enr++)
	      for (j = 0; j < n_el_dofs; j++)
		{
		  int cl = nv+eledges[enr];
		  (*it)[cl][cnts[cl]] = k + j;
		  cnts[cl]++;
		}
	    k += n_el_dofs;
	  }
      

	for (i = 0; i < nv; i++)
	  (*it)[nv+nedge][i] = i;
      
	//  (*testout) << "table: " << (*it) << endl;
	return it;
	*/

      }    
    



#ifdef VEFC_precond
    Array<int> cnts(nv+nedge+nface+ne);
    cnts = 0;
    for (i = 0; i < nv; i++)
      cnts[i] = 1;
    for (i = 0; i < nedge; i++)
      cnts[nv+i] = n_edge_dofs;
    for (i = 0; i < nface; i++)
      cnts[nv+nedge+i] = n_face_dofs;
    for (i = 0; i < ne; i++)
      cnts[nv+nedge+nface+i] = n_el_dofs;

    Table<int> * it = new Table<int> (cnts);

    for (i = 0; i < nv; i++)
      (*it)[i][0] = i;
  
    k = nv;
    for (i = 0; i < nedge; i++)
      for (j = 0; j < n_edge_dofs; j++)
	{
	  (*it)[nv+i][j] = k; 
	  k++; 
	}

    for (i = 0; i < nface; i++)
      for (j = 0; j < n_face_dofs; j++)
	{
	  (*it)[nv+nedge+i][j] = k; 
	  k++; 
	}

    for (i = 0; i < ne; i++)
      for (j = 0; j < n_el_dofs; j++)
	{
	  (*it)[nv+nedge+nface+i][j] = k; 
	  k++; 
	}

    (*testout) << "table: " << (*it) << endl;
    return it;

#endif


#ifdef NONE
    int i, j, k;
    INT_2 edata;
    int enr;
    INT_3 fdata;
    int fnr;

    int nedof = p - 1;
    int pp = p-3;
    int nfdof = (pp*pp+3*pp+2)/2;
    pp = p-4;
    int nidof =  (pp*pp*pp + 6 * pp * pp + 11 * pp + 6) / 6;

  
    blocks.SetSize (ndof);
    for (i = 1; i <= ndof; i++)
      blocks.AddUnique (i, i);
    /*
      int np = ma.GetNP();
      for (i = 1; i <= np; i++)
      blocks.AddUnique (1, i);
    */
    for (i = 1; i <= edgedofs->GetNBags(); i++)
      for (j = 1; j <= edgedofs->GetBagSize(i); j++)
	{
	  edgedofs->GetData (i, j, edata, enr);
	  (*testout) << "i = " << i << ", j = " << j << ",  enr = " << enr << endl;
	  for (k = 0; k < nedof; k++)
	    blocks.AddUnique (enr, enr+k);
	}

    for (i = 1; i <= facedofs->GetNBags(); i++)
      for (j = 1; j <= facedofs->GetBagSize(i); j++)
	{
	  facedofs->GetData (i, j, fdata, fnr);
	  (*testout) << "i = " << i << ", j = " << j << ",  enr = " << fnr << endl;
	  for (k = 0; k < nfdof; k++)
	    blocks.AddUnique (fnr, fnr+k);
	}

    for (i = 1; i <= eldofs.Size(); i++)
      {
	int base = eldofs.Get(i);
	for (k = 0; k < nidof; k++)
	  blocks.AddUnique (base, base+k);
      }
#endif
  }




// register FESpaces
namespace {
  class Init
  { 
  public: 
    Init ();
    };
    
    Init::Init()
    {
      GetFESpaceClasses().AddFESpace ("h1nodal", NodalFESpaceP::Create);
    }

    
    Init init;
}



}


#endif
