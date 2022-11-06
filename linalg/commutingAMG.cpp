/*********************************************************************/
/* File:   commutingamg.cpp                                          */
/* Author: Joachim Schoeberl                                         */
/* Date:   15. Aug. 2002                                             */
/*********************************************************************/

// #define DEBUG
#include <la.hpp>

namespace ngla
{
  
  AMG_H1 :: AMG_H1 (const BaseMatrix & sysmat,
		    Array<INT<2> > & e2v,
		    Array<double> & weighte,
		    int levels)
  {
    // find number of vertices
    int ne = e2v.Size();
    int nv = 0;
    for (int i = 0; i < ne; i++)
      for (int j = 0; j < 2; j++)
	nv = max2 (nv, e2v[i][j]);
    nv++;

    cout << "ne = " << ne << ", nv = " << nv << endl;

    jacobi = 0;
    coarsemat = nullptr;
    inv = 0;

    if (nv < 20 || levels == 0)
      {
	recAMG = 0;
	prol = 0;
	return;
      }

    Array<int> vcoarse(nv), connected(nv);
    Array<bool> edge_collapse(ne);
    Array<double> edge_collapse_weight(ne);
    edge_collapse = 0;
    Array<int> v2edge(nv);
    
    /*
    for (int i = 0; i < weighte.Size(); i++)
      for (int j = 0; j < 2; j++)
	if (e2v[i][j] == -1)
	  cout << "unused edge" << endl;
    */


    /*
    // compute weight to collapse edge
    Array<double> vstrength(nv);
    vstrength = 0.0;
    for (i = 0; i < weighte.Size(); i++)
      for (j = 0; j < 2; j++)
	vstrength[e2v[i][j]] += sqr (weighte[i]);

    for (i = 0; i < ne; i++)
      edge_collapse_weight[i] = 
	sqr (weighte[i]) / min2 (vstrength[e2v[i][0]], vstrength[e2v[i][1]]);
	//	sqr (weighte[i]) / sqrt (vstrength[e2v[i][0]] * vstrength[e2v[i][1]]);

	*/
    // compute weight to collapse edge
    Array<double> vstrength(nv);
    vstrength = 0.0;
    for (int i = 0; i < weighte.Size(); i++)
      if (e2v[i][0] != -1)
	for (int j = 0; j < 2; j++)
	  vstrength[e2v[i][j]] += weighte[i];
    
    for (int i = 0; i < ne; i++)
      if (e2v[i][0] != -1)
	{
	  double vstr1 = vstrength[e2v[i][0]];
	  double vstr2 = vstrength[e2v[i][1]];
	  edge_collapse_weight[i] = 
	    weighte[i] * (vstr1+vstr2) / (vstr1 * vstr2);
	}

    
    // figure out best edges to collapse (iterative improvement)
    v2edge = -1;
    edge_collapse = false;
    bool changed;
    do
      {
	changed = 0;
	for (int i = 0; i < ne; i++)
	  if (e2v[i][0] != -1 && edge_collapse_weight[i] > 0.1)
	    {
	      if (v2edge[e2v[i][0]] == -1 && v2edge[e2v[i][1]] == -1)
		
		{
		  edge_collapse[i] = true;
		  v2edge[e2v[i][0]] = i;
		  v2edge[e2v[i][1]] = i;
		  changed = 1;
		}	    
	      
	      else
		
		{
		  for (int j = 0; j < 2; j++)
		    {
		      int pi1 = e2v[i][j];
		      int pi2 = e2v[i][1-j];
		      
		      if (v2edge[pi1] != -1 && v2edge[pi2] == -1 &&
			  edge_collapse_weight[i] > 
			  edge_collapse_weight[v2edge[pi1]])
			{
			  int remove = v2edge[pi1];
			  edge_collapse[i] = true;
			  edge_collapse[remove] = false;
			  v2edge[e2v[remove][0]] = -1;
			  v2edge[e2v[remove][1]] = -1;
			  v2edge[pi1] = i;
			  v2edge[pi2] = i;
			  changed = 1;
			  break;
			}
		    }
		}
	    }
      }
    while (changed);

    // compute fine vertex to coarse vertex map (vcoarse)
    for (int i = 0; i < nv; i++)
      connected[i] = i;

    for (int i = 0; i < ne; i++)
      if (e2v[i][0] != -1)
      if (edge_collapse[i])
	{
	  int pi1 = e2v[i][0], pi2 = e2v[i][1];
	  if (vstrength[pi1] > vstrength[pi2])
	    connected[pi2] = pi1;
	  else
	    connected[pi1] = pi2;
	}

    int ncv = 0;
    for (int i = 0; i < nv; i++)
      if (connected[i] == i)
	vcoarse[i] = ncv++;

    for (int i = 0; i < nv; i++)
      if (connected[i] != i)
	vcoarse[i] = vcoarse[connected[i]];




    // compute fine edge to coarse edge map (ecoarse)
    HashTable<INT<2>, int> ht_ecoarse(e2v.Size());
    for (int i = 0; i < e2v.Size(); i++)
      if (e2v[i][0] != -1)
      {
	INT<2> ce;
	for (int j = 0; j < 2; j++)
	  ce[j] = vcoarse[e2v[i][j]];
	ce.Sort();
	if (ce[0] != ce[1])
	  ht_ecoarse.Set (ce, -1);
      }
    
    Array<INT<2> > ce2v;
    for (int i = 0; i < ht_ecoarse.Size(); i++)
      for (int j = 0; j < ht_ecoarse.EntrySize(i); j++)
	{
	  INT<2> ce;
	  int efi;
	  ht_ecoarse.GetData (i, j, ce, efi);
	  efi = ce2v.Size();
	  ce2v.Append (ce);
	  ht_ecoarse.SetData (i, j, ce, efi);
	}

    Array<int> ecoarse(ne);
    for (int i = 0; i < e2v.Size(); i++)
      if (e2v[i][0] != -1)
      {
	INT<2> ce;
	for (int j = 0; j < 2; j++)
	  ce[j] = vcoarse[e2v[i][j]];
	ce.Sort();
	if (ce[0] != ce[1])
	  ecoarse[i] = ht_ecoarse.Get(ce);
	else
	  ecoarse[i] = -1;
      }

    // coarse edge weights:
    Array<double> cweighte(ce2v.Size());
    cweighte = 0;
    for (int i = 0; i < e2v.Size(); i++)
      if (e2v[i][0] != -1)
	if (ecoarse[i] != -1)
	  cweighte[ecoarse[i]] += weighte[i];
    

    // compute prolongation matrix 
    // prol  ... for matrix projection

    // piecewise constant:
    Array<int> nne(nv);
    nne = 1;

    prol = new SparseMatrix<double> (nne, ncv);
    for (int i = 0; i < nv; i++)
      (*prol)(i, vcoarse[i]) = 1;

    recAMG = new AMG_H1 (sysmat, ce2v, cweighte, levels-1);
  }
  

  AMG_H1 :: ~AMG_H1 ()
  {
    delete prol;
    delete recAMG;

    // delete jacobi;
    // delete coarsemat;
    // delete inv;
  }
  

  void AMG_H1 :: ComputeMatrices (const BaseSparseMatrix & mat)
  {
    pmat = &mat;

    if (0) // if (prol)
      {
	// define smoothing blocks
	int nv = prol->Width();
	Array<int> cnt(nv);
	cnt = 0;
	for (int i = 0; i < prol->Height(); i++)
	  for (int j = 0; j < prol->GetRowIndices(i).Size(); j++)
	    cnt[prol->GetRowIndices(i)[j]]++;
	
	Table<int> smblocks(cnt);
	cnt = 0;
	for (int i = 0; i < prol->Height(); i++)
	  for (int j = 0; j < prol->GetRowIndices(i).Size(); j++)
	    {
	      int jj = prol->GetRowIndices(i)[j];
	      smblocks[jj][cnt[jj]] = i;
	      cnt[jj]++;
	    }
	// *testout << "prol = " << endl << *prol << endl;
	// *testout << "smoothing blocks = " << endl << *smblocks << endl;
	
	bjacobi = mat.CreateBlockJacobiPrecond (make_shared<Table<int>> (smblocks));
      }
    jacobi = mat.CreateJacobiPrecond ();

    if (recAMG)
      {
	coarsemat = shared_ptr<BaseSparseMatrix>(mat.Restrict (*prol));
	recAMG -> ComputeMatrices (*coarsemat);
	inv = 0;
      }
    else 
      {
        // const_cast<BaseSparseMatrix&> (mat).SetInverseType ( SPARSECHOLESKY );
	mat.SetInverseType (SPARSECHOLESKY);
        inv = mat.InverseMatrix();
      }
  }

  size_t AMG_H1 :: NZE() const
  {
    size_t nze = pmat->NZE();
    if (recAMG)
      nze += recAMG->NZE();
    return nze;
  }

  
  void AMG_H1 :: Mult (const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("H1-AMG::Mult");
    RegionTimer reg (timer);

    if (inv)
      {
	y = (*inv) * x;
	return;
      }

    auto hv = pmat->CreateColVector();
    auto wc = coarsemat->CreateColVector();
    auto dc = coarsemat->CreateColVector();

    y = 0;
    jacobi->GSSmooth (y, x);

    if (recAMG)
      {
	hv = x - (*pmat) * y;
	dc = Transpose (*prol) * hv;

	if (recAMG)
	  recAMG -> Mult (dc, wc);

	y += (*prol) * wc;
      }

    jacobi->GSSmoothBack (y, x);
  }
      



  /* ************************* HCurl AMG ********************************* */




  void SortFace (INT<4> & face)
  {
    if (face[3] == -1)
      {
	if (face[0] > face[1]) swap (face[0], face[1]);
	if (face[1] > face[2]) swap (face[1], face[2]);
	if (face[0] > face[1]) swap (face[0], face[1]);
      }
    else
      {
	while (face[1] < face[0] || face[2] < face[0] || face[3] < face[0])
	  {
	    int hi = face[0];
	    face[0] = face[1];
	    face[1] = face[2];
	    face[2] = face[3];
	    face[3] = hi;
	  }

	if(face[1] > face[3]) // sensible?
	  swap(face[1],face[3]);
      }
  }


  
  AMG_HCurl :: AMG_HCurl (const BaseMatrix & sysmat,
			  const Array<Vec<3> > & vertices,
                          Array<INT<2> > & e2v,
			  Array<INT<4> > & f2v,
			  Array<double> & weighte,
			  Array<double> & weightf,
			  int levels)
  {
    static int timer1 = NgProfiler::CreateTimer ("AMG, reg1");
    static int timer2 = NgProfiler::CreateTimer ("AMG, reg2");
    static int timer2b = NgProfiler::CreateTimer ("AMG, reg2b");
    static int timer3 = NgProfiler::CreateTimer ("AMG, reg3");
    static int timer4 = NgProfiler::CreateTimer ("AMG, reg4");
    
    NgProfiler::StartTimer (timer1);



    if (vertices.Size())   // visualization of faces
      {
        char buf[20];
        snprintf (buf, 20, "amg.level%d.surf", levels);
        cout << "filename = " << buf << endl;
        ofstream out (buf);
        out << "surfacemesh" << endl;
        out << vertices.Size() << endl;
        for (int i = 0; i < vertices.Size(); i++)
          out << vertices[i](0) << " " 
              << vertices[i](1) << " " 
              << vertices[i](2) << "\n";

        double minweight=1e10;
        int cntf = 0;
        for (int i = 0; i < f2v.Size(); i++)
          if (weightf[i] < minweight) cntf++;

        out << f2v.Size() << endl;
        for (int i = 0; i < f2v.Size(); i++)
	  if (weightf[i] < minweight) 
	    {
	      out << f2v[i][0]+1 << " "
		  << f2v[i][1]+1 << " "
		  << f2v[i][2]+1;
	      if(f2v[i][3] != -1)
		out << " " << f2v[i][3]+1;
	      out << endl;
	    }
      }



    for (int i = 0; i < e2v.Size(); i++)
      if (e2v[i][0] == -1)
        cout << "es gibt sie doch !!!" << endl;




    // find number of vertices
    int ne = e2v.Size();
    int nf = f2v.Size();
    int nv = 0;
    for (int i = 0; i < ne; i++)
      for (int j = 0; j < 2; j++)
	if (e2v[i][j] > nv)
	  nv = e2v[i][j];
    nv++;

    cout << "nfa = " << nf << ", ned = " << ne << ", nv = " << nv << endl;

    Array<bool> edge_collapse(ne);
    Array<double> edge_collapse_weight(ne);
    edge_collapse = 0;


    // find loops of 3 edges without face
    HashTable<INT<3>, int> faceht (f2v.Size());
    HashTable<INT<4>, int> qfaceht (f2v.Size());

    // bool quadfacesexist = false;
    
    for (int i = 0; i < f2v.Size(); i++)
      {
       	if(f2v[i][3] == -1)
	  {
	    INT<3> hv  (f2v[i][0], f2v[i][1], f2v[i][2]);
	    hv.Sort();
	    faceht.Set (hv, i);
	  }
	else
	  {
	    INT<4> hv  (f2v[i][0], f2v[i][1], f2v[i][2], f2v[i][3]);
	    hv.Sort();
	    qfaceht.Set (hv, i);
	    // quadfacesexist = true;
	  }
      }

    Array<int> cnt_v2e(nv);
    cnt_v2e = 0;
    for (int i = 0; i < ne; i++)
      {
        cnt_v2e[e2v[i][0]]++;
        cnt_v2e[e2v[i][1]]++;
      }
    Table<int> v2e(cnt_v2e);
    cnt_v2e = 0;
    for (int i = 0; i < ne; i++)
      {
        v2e[e2v[i][0]][cnt_v2e[e2v[i][0]]++] = i;
        v2e[e2v[i][1]][cnt_v2e[e2v[i][1]]++] = i;
      }

    NgProfiler::StopTimer (timer1);
    
    NgProfiler::StartTimer (timer2);


    
    // for (int xx = 0; xx < 5; xx++)
    bool appended = true;
    while(appended)
      {
	appended = false;
	for (int i = 0; i < nv; i++)
	  {
	    for (int j = 0; j < v2e[i].Size(); j++)
	      for (int k = 0; k < j; k++)
		{
		  int ej = v2e[i][j];
		  int ek = v2e[i][k];
		  int vi2 = e2v[ej][0]+e2v[ej][1]-i;
		  int vi3 = e2v[ek][0]+e2v[ek][1]-i;
		  INT<3> f123(i, vi2, vi3);
		  f123.Sort();
		  if (faceht.Used(f123))
		    for (int l = 0; l < k; l++)
		      {
			int el = v2e[i][l];
			int vi4 = e2v[el][0]+e2v[el][1]-i;
			INT<3> f124(i, vi2, vi4);
			INT<3> f134(i, vi3, vi4);
			f124.Sort();
			f134.Sort();
			if (faceht.Used(f124) && faceht.Used(f134))
			  {
			    
			    // 			double w = 1.0 / (1.0 / weightf[faceht.Get(f123)] + 
			    // 					  1.0 / weightf[faceht.Get(f124)] + 
			    // 					  1.0 / weightf[faceht.Get(f134)]);
			    
			    
			    double w =  min3 (weightf[faceht.Get(f123)],
					      weightf[faceht.Get(f124)],
					      weightf[faceht.Get(f134)]);
			    INT<3> f234(vi2, vi3, vi4);
			    f234.Sort();
			    if (faceht.Used(f234))
			      {
				int fi = faceht.Get(f234);
				weightf[fi] = max2 (weightf[fi], w);
			      }
			    else
			      {
				f2v.Append (INT<4> (vi2, vi3, vi4, -1));
				weightf.Append (w);
				faceht.Set (f234, weightf.Size()-1);
				nf++;
				appended = true;
			      }
			  }
		      }
		}
	  }
      }


    NgProfiler::StopTimer (timer2);
    
    NgProfiler::StartTimer (timer2b);


    /*
      for (int i = 0; i < ne; i++)
      {
        int vi1 = e2v[i][0];
        int vi2 = e2v[i][1];
        for (int j = 0; j < v2e[vi1].Size(); j++)
          for (int k = 0; k < v2e[vi2].Size(); k++)
            {
              int e2 = v2e[vi1][j];
              int e3 = v2e[vi2][k];
              int v3 = e2v[e2][0]+e2v[e2][1]-vi1;
              int v4 = e2v[e3][0]+e2v[e3][1]-vi2;

              if (i < e2 && i < e3 && v3 == v4)
                {
                  INT<3> face(vi1, vi2, v3);
                  face.Sort();

                  if (!faceht.Used(face))
                    {
                      *testout << "level = " << levels 
                               << " has no face " << face << endl;
                      cout << "loop found !!" << endl;
                      f2v.Append (INT<4> (vi1, vi2, v3, -1));
                      weightf.Append (0);
                      nf++;
                    }
                  else
                    *testout << "has face " << face << endl;
                }
                
            }
      }
    */


    // compute face 2 edge table
    HashTable<INT<2>, int> ht_edge(e2v.Size());
    Array<INT<4> > f2e(nf);

    for (int i = 0; i < ne; i++)
      if (e2v[i][0] != -1)
	{
	  INT<2> ce = e2v[i];
	  ce.Sort();
	  ht_edge.Set (ce, i);
	}

    for (int i = 0; i < f2v.Size(); i++)
      if (f2v[i][0] != -1)
	{
	  int nfv = (f2v[i][3] == -1) ? 3 : 4;
	  f2e[i][3] = -1;
	  for (int j = 0; j < nfv; j++)
	    {
	      INT<2> ce;
	      ce[0] = f2v[i][j];
	      ce[1] = f2v[i][(j+1)%nfv];
	      ce.Sort();
	      if (!ht_edge.Used (ce))
		{
		  /*
		    cout << "Err: unused edge, " 
		    << "face = " << f2v[i] << endl;
		*/
		  (*testout) << "Err: unused edge, " 
			     << "face = " << f2v[i] << endl;
		  cout << "Err: unused edge, " 
                       << "face = " << f2v[i] << endl;
		  f2e[i][j] = -1;
		}
	      else
		f2e[i][j] = ht_edge.Get(ce);
	    }
	}
      else
	for (int j = 0; j < 4; j++)
	  f2e[i][j] = -1;
	  
    
    //    (*testout) << "weightf = " << weightf << endl;
    


    Array<double> sume(ne);
    sume = 0;
    for (int i = 0; i < nf; i++)
      for (int j = 0; j < 4; j++)
        if (f2e[i][j] != -1)
          sume[f2e[i][j]] += sqr (weightf[i]);




    Array<double> face_collapse_weight(nf);
    for (int i = 0; i < nf; i++)
      {
        double mine = 1e99;
        for (int j = 0; j < 4; j++)
          if (f2e[i][j] != -1)
            mine = min2 (mine, sume[f2e[i][j]]);
        face_collapse_weight[i] = sqr (weightf[i]) / mine;
  
        //       double maxe = 0;
        //       for (j = 0; j < 4; j++)
        //       if (f2e[i][j] != -1)
        //       maxe = max2 (maxe, sume[f2e[i][j]]);
        //       face_collapse_weight[i] = sqr (weightf[i]) / maxe;

      }

    /*
    edge_collapse_weight = 1e99;
    for (i = 0; i < nf; i++)
      for (j = 0; j < 4; j++)
        if (f2e[i][j] != -1)
          edge_collapse_weight[f2e[i][j]] = 
            min2 (edge_collapse_weight[f2e[i][j]],
                  face_collapse_weight[i]);
    */

    Array<int> cnt_e2f(ne);
    cnt_e2f = 0;
    for (int i = 0; i < nf; i++)
      for (int j = 0; j < 4; j++)
        if (f2e[i][j] != -1)
          cnt_e2f[f2e[i][j]]++;

    Table<int> e2f(cnt_e2f);
    cnt_e2f = 0;
    for (int i = 0; i < nf; i++)
      for (int j = 0; j < 4; j++)
        if (f2e[i][j] != -1)
          e2f[f2e[i][j]][cnt_e2f[f2e[i][j]]++] = i;
  

    for (int i = 0; i < ne; i++)
      {
        double mine = 1e99;
        for (int j = 0; j < e2f[i].Size(); j++)
          {
            int fnr = e2f[i][j];
            double maxf = 0;
            for (int k = 0; k < 4; k++)
              {
                int enr = f2e[fnr][k];
                if (enr != -1  && enr != i)
                  maxf = max2 (maxf, sqr(weightf[fnr])/sume[enr]);
              }
            mine = min2(mine, maxf);
          }
        edge_collapse_weight[i] = mine;
      }


    /*
    for (double fac = 1; fac > 1e-10; fac *= 0.1)
      {
        int cnt = 0;
        for (int i = 0; i < ne; i++)
          if (edge_collapse_weight[i] < fac) cnt++;
        cout << cnt << " edges have weight < " << fac << endl;
      }
    */


    NgProfiler::StopTimer (timer2b);
    NgProfiler::StartTimer (timer3);



    // figure out best edges to collapse (iterative improvement)
    Array<int> v2edge(nv);
    v2edge = -1;
    edge_collapse = 0;

    for (int cnt = 0, changed = 1; cnt < 3 && changed; cnt++)
      {
        changed = 0;

	for (int i = 0; i < ne; i++)
          {
            if (e2v[i][0] == -1) continue;
            if (edge_collapse_weight[i] < 0.05) continue;
       
       
            bool ok = true;
            int vi1 = e2v[i][0];
            int vi2 = e2v[i][1];

            // check, whether weak holes will collapse:
            for (int j = 0; j < v2e[vi1].Size(); j++)
              {
                int e2 = v2e[vi1][j];
                if (e2 == i) continue;
                int vi3 = e2v[e2][0]+e2v[e2][1]-vi1;


                // check whether weak faces will be closed
                INT<2> ep2(vi2, vi3);
                ep2.Sort();
                if (ht_edge.Used (ep2))
                  {
                    INT<3> face(vi1, vi2, vi3);
                    face.Sort();
                    
                    double val = (faceht.Used(face)) ? face_collapse_weight[faceht.Get(face)] : 0;
                    if (val < 0.01) { ok = false; break; }
		  }
		else
		  {
// 		    bool foundquad = false;
// 		    for(int k=0; ok && k < v2e[vi2].Size(); k++)
// 		      {
// 			int e3 = v2e[vi2][k];
// 			if(e3 == i) continue;
// 			int vi4 = e2v[e3][0]+e2v[e3][1]-vi2;
			
// 			INT<2> ep2(vi3,vi4);
// 			ep2.Sort();
// 			if (ht_edge.Used (ep2))
// 			  {
// 			    INT<4> face(vi1,vi2,vi3,vi4);
// 			    face.Sort();

// 			    if(qfaceht.Used(face))
// 			      {
// 				foundquad = true;
// 				if(face_collapse_weight[qfaceht.Get(face)] < 0.01)
// 				  ok = false;
// 			      }
// 			  }
// 		      }
// 		    if(!foundquad && v2edge[vi3] != -1)
		    if(v2edge[vi3] != -1)
		      {
			int eop = v2edge[vi3];
			int vi4 = e2v[eop][0]+e2v[eop][1]-vi3;
			if(vi4 == vi1 || vi4 == vi2) continue; // sensible?
			
			INT<2> epair(vi2, vi4);
			epair.Sort();
			if (ht_edge.Used (epair))
			  {
			    INT<4> f1234(vi1,vi2,vi3,vi4);  // sensible?
			    f1234.Sort();
			    double v1234 = (qfaceht.Used(f1234)) ? face_collapse_weight[qfaceht.Get(f1234)] : 0;
			    if (v1234 > 0.01) continue;	
			    

			    // not possible
			    //INT<3> f123(vi1, vi2, vi3);
			    //f123.Sort();
			    //if(faceht.Used(f123))
			    //  cout << endl << "POSSIBLE!!" << endl;
			    //double v123 = (faceht.Used(f123)) ? face_collapse_weight[faceht.Get(f123)] : 0;
			    //if (v123 > 0.01) continue;
			    
			    INT<3> f124(vi1, vi2, vi4);
			    f124.Sort();
			    double v124 = (faceht.Used(f124)) ? face_collapse_weight[faceht.Get(f124)] : 0;
			    if (v124 > 0.01) continue;
			    
			    ok = false; 
			    break;
			  }
		      }
		  }
              }


            if (!ok) continue;

            if (v2edge[e2v[i][0]] == -1 && v2edge[e2v[i][1]] == -1)
              
              {
                edge_collapse[i] = 1;
                v2edge[e2v[i][0]] = i;
                v2edge[e2v[i][1]] = i;
                changed = 1;
              }	    
            else
              {  // iterative improvement (really improving ?)
                for (int j = 0; j < 2; j++)
                  {
                    int pi1 = e2v[i][j];
                    int pi2 = e2v[i][1-j];
		    
                    if (v2edge[pi1] != -1 && v2edge[pi2] == -1 &&
                        edge_collapse_weight[i] > 
                        edge_collapse_weight[v2edge[pi1]])
                      {
                        int remove = v2edge[pi1];
                        edge_collapse[i] = 1;
                        edge_collapse[remove] = 0;
                        v2edge[e2v[remove][0]] = -1;
                        v2edge[e2v[remove][1]] = -1;
                        v2edge[pi1] = i;
                        v2edge[pi2] = i;
                        changed = 1;
                        break;
                      }
                  }
              }
          }
      }

    NgProfiler::StopTimer (timer3);
    NgProfiler::StartTimer (timer4);




    // compute fine vertex to coarse vertex map (vcoarse)
    Array<int> vcoarse(nv), connected(nv);
    for (int i = 0; i < nv; i++)
      connected[i] = i;

    for (int i = 0; i < ne; i++)
      if (edge_collapse[i])
	{
	  int pi1 = e2v[i][0];
	  int pi2 = e2v[i][1];
	  int mini = min2(pi1, pi2);
	  connected[pi1] = mini;
	  connected[pi2] = mini;
	}

    int ncv = 0;
    for (int i = 0; i < nv; i++)
      if (connected[i] == i)
	{
	  vcoarse[i] = ncv;
	  ncv++;
	}
      else
	vcoarse[i] = vcoarse[connected[i]];


    Array<Vec<3> > cvertices;
    if (vertices.Size())
      {
        cvertices.SetSize(ncv);
        cvertices = Vec<3> (0,0,0);
        for (int i = 0; i < nv; i++)
          cvertices[vcoarse[i]] = vertices[i];
      }




    // compute fine edge to coarse edge map (ecoarse)
    HashTable<INT<2>, int> ht_ecoarse(ne);
    for (int i = 0; i < ne; i++)
      if (e2v[i][0] != -1)
	{
	  INT<2> ce;
	  for (int j = 0; j < 2; j++)
	    ce[j] = vcoarse[e2v[i][j]];
	  if (ce[0] != ce[1])
	    {
	      if (ce[0] > ce[1]) Swap (ce[0], ce[1]);
	      ht_ecoarse.Set (ce, -1);
	    }
	}
    
    Array<INT<2> > ce2v;
    for (int i = 0; i < ht_ecoarse.Size(); i++)
      for (int j = 0; j < ht_ecoarse.EntrySize(i); j++)
	{
	  INT<2> ce;
	  int efi;
	  ht_ecoarse.GetData (i, j, ce, efi);
	  efi = ce2v.Size();
	  ce2v.Append (ce);
	  ht_ecoarse.SetData (i, j, ce, efi);
	}

    Array<int> ecoarse(ne);
    Array<short> ecoarsesign(ne);
    ecoarsesign = 1;
    for (int i = 0; i < ne; i++)
      {
	if (e2v[i][0] == -1)
	  {
	    ecoarsesign[i] = 0;
	    ecoarse[i] = -1;
	    continue;
	  }

	INT<2> ce;
	for (int j = 0; j < 2; j++)
	  ce[j] = vcoarse[e2v[i][j]];
	
	if (ce[0] == ce[1])
	  {
	    ecoarsesign[i] = 0;
	    ecoarse[i] = -1;
	  }
	else 
	  {
	    if (ce[0] > ce[1]) 
	      {
		Swap (ce[0], ce[1]);
		ecoarsesign[i] = -1;
	      }
	    ecoarse[i] = ht_ecoarse.Get(ce);
	  }
      }

    // coarse edge weights:
    Array<double> cweighte(ce2v.Size());
    cweighte = 0;
    for (int i = 0; i < e2v.Size(); i++)
      if (ecoarse[i] != -1)
	cweighte[ecoarse[i]] += weighte[i];






    // compute fine face to coarse face map (fcoarse)
    HashTable<INT<4>, int> ht_fcoarse(nf);
    for (int i = 0; i < nf; i++)
      {
	bool valid = true;
	for (int j = 0; j < 3; j++)
	  if (f2e[i][j] == -1)
	    valid = false;
	if(!valid)
	  {
	    cout << "not valid" << endl;
	    continue;
	  }

	INT<4> cf;
	for (int j = 0; j < 4; j++)
	  if (f2v[i][j] != -1)
	    cf[j] = vcoarse[f2v[i][j]];
	  else
	    cf[j] = -1;
	SortFace(cf);
	bool degenerated = 0;
	for (int j = 0; j < 3; j++)
	  for (int k = j+1; k < 4; k++)
	    if (cf[j] == cf[k])
	      degenerated = 1;
	//	if (cf[0] != cf[1] && cf[1] != cf[2] && cf[2] != cf[3] && cf[3] != cf[0])
	if (!degenerated)
	  ht_fcoarse.Set (cf, -1);
      }
    
    Array<INT<4> > cf2v;
    for (int i = 0; i < ht_fcoarse.Size(); i++)
      for (int j = 0; j < ht_fcoarse.EntrySize(i); j++)
	{
	  INT<4> cf;
	  int cfi;
	  ht_fcoarse.GetData (i, j, cf, cfi);
	  cfi = cf2v.Size();
	  cf2v.Append (cf);
	  ht_fcoarse.SetData (i, j, cf, cfi);
	}

    Array<int> fcoarse(nf);
    fcoarse = -1;
    for (int i = 0; i < nf; i++)
      {
	bool valid = 1;
	for (int j = 0; j < 3; j++)
	  if (f2e[i][j] == -1)
	    valid = 0;
	if(!valid)
	  {
	    cout << "not valid" << endl;
	    continue;
	  }

	INT<4> cf;
	for (int j = 0; j < 4; j++)
	  if (f2v[i][j] != -1)
	    cf[j] = vcoarse[f2v[i][j]];
	  else
	    cf[j] = -1;
	SortFace(cf);
	bool degenerated = 0;
	for (int j = 0; j < 3; j++)
	  for (int k = j+1; k < 4; k++)
	    if (cf[j] == cf[k])
	      degenerated = 1;
	// if (cf[0] != cf[1] && cf[1] != cf[2] && cf[2] != cf[3] && cf[3] != cf[0])
	if (!degenerated)
	  fcoarse[i] = ht_fcoarse.Get(cf);
      }

    // coarse face weights:
    Array<double> cweightf(cf2v.Size());
    cweightf = 0;
    for (int i = 0; i < nf; i++)
      if (fcoarse[i] != -1)
	cweightf[fcoarse[i]] += weightf[i];


    // compute prolongation matrix
    Array<int> nne(ne);
    nne = 0;
    for (int i = 0; i < ne; i++)
      if (ecoarse[i] != -1)
	nne[i] = 1;

    prol = new SparseMatrix<double> (nne);
    // prol = dynamic_cast< SparseMatrixTM<double>* >(sysmat.CreateMatrix(nne));
    for (int i = 0; i < ne; i++)
      if (ecoarse[i] != -1)
	{
	  prol->CreatePosition (i, ecoarse[i]);
	  (*prol)(i, ecoarse[i]) = ecoarsesign[i];
	}

    // compute gradient matrix:
    for (int i = 0; i < ne; i++)
      if (e2v[i][0] == -1)
	nne[i] = 0;
      else
	nne[i] = 2;

    grad = new SparseMatrix<double> (nne);
    // grad = dynamic_cast< SparseMatrixTM<double>* >(sysmat.CreateMatrix(nne));
    for (int i = 0; i < ne; i++)
      if (e2v[i][0] != -1)
      {
	grad->CreatePosition (i, e2v[i][0]);
	grad->CreatePosition (i, e2v[i][1]);
	(*grad)(i, e2v[i][0]) = 1;
	(*grad)(i, e2v[i][1]) = -1;
      }

    NgProfiler::StopTimer (timer4);

    cout << "nv = " << nv << ", ncv = " << ncv << endl;
    if (nv > 50 && ncv < 0.9*nv && levels != 0)
      {
	recAMG = new AMG_HCurl (sysmat, cvertices, ce2v, cf2v, cweighte, cweightf, levels-1);
        h1AMG = new AMG_H1 (sysmat, e2v, weighte, levels);
        // h1AMG = new AMG_H1 (e2v, weighte, 0);
      }
    else
      {
	recAMG = 0;
	h1AMG = 0;
      }
  }


  AMG_HCurl :: ~AMG_HCurl ()
  {
    delete prol;
    delete recAMG;
  }
  



  void AMG_HCurl :: ComputeMatrices (const BaseSparseMatrix & mat)
  {
    cout << "compute HCurl matrices" << endl;

    pmat = &mat;
    coarsemat = shared_ptr<BaseSparseMatrix>(mat.Restrict (*prol));
    jacobi = mat.CreateJacobiPrecond ();

    h1mat = mat.Restrict (*grad);
    //    dynamic_cast<SparseMatrixSymmetric<Mat<1,1> >&> (*h1mat) (0,0) += 1;
    dynamic_cast<SparseMatrixTM<double>&> (*h1mat)(0,0) += 1;

    if (recAMG)
      {
	recAMG -> ComputeMatrices (*coarsemat);
	h1AMG -> ComputeMatrices (*h1mat);
	inv = 0;
	/*
	EigenSystem eigen (mat, *this);
	eigen.Calc();
	eigen.PrintEigenValues (cout);
	*/
      }
    else 
      {
	cout << "cal inverse, size = " << mat.Height() << endl;
        const_cast<BaseSparseMatrix&> (mat).SetInverseType ( SPARSECHOLESKY );
	inv = mat.InverseMatrix();
      }
  }



  size_t AMG_HCurl :: NZE() const
  {
    size_t nze = pmat->NZE() + h1mat->NZE();
    if (recAMG) 
      nze += recAMG->NZE() + h1AMG->NZE();

    return nze;
  }


  
  void AMG_HCurl :: Mult (const BaseVector & x, BaseVector & y) const
  {
    if (inv)
      {
	y = (*inv) * x;
	return;
      }

    auto hv = pmat->CreateColVector();
    auto res = pmat->CreateColVector();
    auto wc = coarsemat->CreateColVector();
    auto dc = coarsemat->CreateColVector();

    auto wh1 = h1mat->CreateColVector();
    auto dh1 = h1mat->CreateColVector();
    
    y = 0;
    res = x;
    hv = 0;

    // for (int k = 0; k < 5; k++)
    jacobi->GSSmooth (y, x, res /* , hv */);
    pmat -> MultAdd1 (-1, y, res);   // update residual

    // do potential smoothing:
    hv = res;
    // hv = x - (*pmat) * y;
    dh1 = Transpose (*grad) * hv;
    wh1 = (*h1AMG) * dh1;
    y += (*grad) * wh1;

    if (recAMG)
      {
	hv = x - (*pmat) * y;
	dc = Transpose (*prol) * hv;

	if (recAMG)
	  recAMG -> Mult (dc, wc);

	y += (*prol) * wc;
      }

    hv = x - (*pmat) * y;
    dh1 = Transpose (*grad) * hv;
    wh1 = (*h1AMG) * dh1;
    y += (*grad) * wh1;

    // for (int k = 0; k < 5; k++)
      jacobi->GSSmoothBack (y, x);
  }
}


