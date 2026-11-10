/* *************************************************************************/
/* File:   order.cc                                                        */
/* Author: Joachim Schoeberl                                               */
/* Date:   25. Mar. 2000                                                   */
/* *************************************************************************/



// #include <la.hpp>
#include <bit>
#include <ngstd.hpp>
#include "order.hpp"

namespace ngla
{
  

  /* 

  Finds ordering for sparse Cholesky Factorization
  Algorithm: Minimum Degree
              
  See:       The evolution of the 
  minimum degree ordering algorithm

  A. George and J.W.H. Liu
  SIAM Review, Vol31, 1989, pp1-19

  */


  // ngstd::BlockAllocator CliqueEl :: ball(sizeof (CliqueEl));

  void MDOVertex::DoArchive(Archive& ar)
  {
    ar & master & nextminion & numminions & numcliques & eliminated
      & used & flag & nconnected;
    ar.Do(connected, nconnected);
  }
  
  MinimumDegreeOrdering :: MinimumDegreeOrdering (int an)
    :  n(an), 
       cliques(an), order(an), blocknr(an), vertices(an), 
       priqueue(an, an+1),
       ball(sizeof (CliqueEl), 1000)
  {
    static Timer t("MinimumDegreeOrdering::ctor"); RegionTimer r(t);
    /*
    cliques = NULL;
    blocknr = 0;
    order = 0;
    */
    /*
    for (int i = 0; i < n; i++)
      vertices[i].Init(i);
    */
    ParallelForRange (n, [&] (IntRange r)
                      {
                        cliques.Range(r) = NULL;
                        blocknr.Range(r) = 0;
                        order.Range(r) = 0;
                        
                        for (auto i : r)
                          vertices[i].Init(i);
                      });
  }

  void MinimumDegreeOrdering::DoArchive(Archive & ar)
  {
    ar & n & nused & order & blocknr & vertices & priqueue;
    ar & edges & adjstart & adjdata & adjlen & adjw & vmark & markid & approx & curstamp;
    if(ar.Output())
      {
        ar << cliques.Size();
        for(auto clique : cliques)
          ar << clique->vnr << clique->eliminate << clique->Flag();
        for(auto clique : cliques)
          if(clique)
            {
              ar << cliques.Pos(clique->next);
              ar << cliques.Pos(clique->nextcl);
              ar << cliques.Pos(clique->clmaster);
            }
      }
    else
      {
        size_t clique_size;
        ar & clique_size;
        cliques.SetSize(clique_size);
        for(auto i : Range(clique_size))
          {
            int vnr;
            bool flag, eliminate;
            ar & vnr & eliminate & flag;
            cliques[i] = new(ball) CliqueEl(vnr);
            cliques[i]->SetFlag(flag);
            cliques[i]->eliminate = eliminate;
          }
        for(auto clique : cliques)
          if(clique)
            {
              size_t next, nextcl, clmaster;
              ar & next & nextcl & clmaster;
              clique->next = cliques[next];
              clique->nextcl = cliques[nextcl];
              clique->clmaster = cliques[clmaster];
            }
      }
  }

  void MinimumDegreeOrdering :: AddEdge (int v1, int v2)
  {
    if (v1 == v2) return;
    edges.Append (v1);
    edges.Append (v2);
  }


  void MinimumDegreeOrdering :: PrintCliques () 
  {
    for (int i = 0; i < n; i++)
      if (!vertices[i].Eliminated())
	{
	  (*testout) << "Vertex " << i << ", degree = " 
		     << CalcDegree (i) 
		     << ", adjacency = " << Adj(i) << endl;
	  
	  for (CliqueEl * p1 = cliques[i]; p1; p1 = p1->nextcl)
	    {
	      CliqueEl * p2 = p1;
	      (*testout) << "( ";
	      do
		{
		  if (!vertices[p2->GetVertexNr()].Eliminated())
		    (*testout) << p2->GetVertexNr() << " ";
		  p2 = p2->next;
		}
	      while (p2 != p1);
	      (*testout) << ")";
	    }
	  (*testout) << endl;
	}
  }



  int MinimumDegreeOrdering :: CalcDegree (int v1)
  {
    // distinct vertices of the elements and the adjacency of v1, weighted
    // with their minions; v1 itself is counted too
    for (CliqueEl * p1 = cliques[v1]; p1; p1 = p1->nextcl)
      {
	CliqueEl * p2 = p1;
	do { vertices[p2->GetVertexNr()].SetUsed(0); p2 = p2->next; } while (p2 != p1);
      }
    for (int w : Adj(v1)) vertices[w].SetUsed(0);
    vertices[v1].SetUsed(0);

    int deg = 0;
    auto count = [&] (int v2)
    {
      if (!vertices[v2].Used() && !vertices[v2].Eliminated() && IsMaster(v2))
        {
          deg += 1+NumMinions(v2);
          vertices[v2].SetUsed(1);
        }
    };
    count (v1);
    for (CliqueEl * p1 = cliques[v1]; p1; p1 = p1->nextcl)
      {
	CliqueEl * p2 = p1;
	do { count (p2->GetVertexNr()); p2 = p2->next; } while (p2 != p1);
      }
    for (int w : Adj(v1)) count (w);
    return deg;
  }
 


  void MinimumDegreeOrdering :: EliminateMasterVertex (int v)
  {
    static Timer t("MDO::EliminateMaster", NoTracing);
    RegionTimer reg(t);
    curstamp++;
    ReleaseChain (v);

    // the elements of v go, their members and the adjacency of v form the
    // new element (masters only, each once; v itself excluded via Used)
    for (CliqueEl * p1 = cliques[v]; p1; p1 = p1->nextcl)
      {
	CliqueEl * p2 = p1;
	do { vertices[p2->GetVertexNr()].SetUsed (0); p2->eliminate = 1; p2 = p2->next; } while (p2 != p1);
      }
    for (int w : Adj(v)) vertices[w].SetUsed (0);
    vertices[v].SetUsed (1);

    CliqueEl * newp = NULL;
    auto addmember = [&] (int w)
    {
      if (vertices[w].Used() || vertices[w].Eliminated() || !IsMaster(w)) return;
      CliqueEl * p3 = new (ball) CliqueEl (w);
      p3->next = newp;
      p3->clmaster = newp ? newp->clmaster : p3;
      p3->eliminate = 0;
      p3->nextcl = NULL;
      p3->SetFlag (false);
      newp = p3;
      vertices[w].SetUsed (1);
    };
    // in the order the clique lists had them (elements, then the edges
    // last added first): the ring order decides ties in the priority queue
    for (CliqueEl * p1 = cliques[v]; p1; p1 = p1->nextcl)
      {
	CliqueEl * p2 = p1;
	do { addmember (p2->GetVertexNr()); p2 = p2->next; } while (p2 != p1);
      }
    {
      auto a = Adj(v);
      for (int k = a.Size()-1; k >= 0; k--) addmember (a[k]);
    }
    adjlen[v] = 0;

    if (!newp)
      {
        // no neighbours left: the block is v and its minions
	CliqueEl * p1 = cliques[v];
	while (p1)
	  {
	    CliqueEl * p2 = p1->GetNextClique();
	    ball.Free(p1);
	    p1 = p2;
	  }
	cliques[v] = NULL;	
        delete [] vertices[v].connected;
        int cnt = NumMinions(v);
	vertices[v].nconnected = cnt;
	vertices[v].connected = cnt ? new int[cnt] : NULL;
        cnt = 0;
        for (int hv = NextMinion(v); hv != -1; hv = NextMinion(hv))
          vertices[v].connected[cnt++] = hv;
	return;
      }

    // close the new element, its weight, flag its members
    newp -> clmaster -> next = newp;
    {
      int wp = 0;
      CliqueEl * p3 = newp;
      do { wp += 1 + NumMinions (p3->Nr()); vertices[p3->Nr()].SetFlag (1); p3 = p3->next; } while (p3 != newp);
      newp->clmaster->weight = wp;
    }

    // dominated elements: all members in the new element
    if (approx)
      {
        // one pass over the members: w(e) = weight of e outside the new
        // element, zero means dominated
        touched.SetSize0();
        CliqueEl * p3 = newp;
        do
          {
            int wi = 1 + NumMinions (p3->Nr());
            for (CliqueEl * p1 = cliques[p3->Nr()]; p1; p1 = p1->nextcl)
              {
                if (p1->eliminate) continue;
                CliqueEl * cm = p1->clmaster;
                if (cm->stamp != curstamp) { cm->stamp = curstamp; cm->w = cm->weight; touched.Append (cm); }
                cm->w -= wi;
              }
            p3 = p3->next;
          }
        while (p3 != newp);

        for (CliqueEl * cm : touched)
          if (cm->w == 0)
            {
              CliqueEl * p2 = cm;
              do { p2->eliminate = 1; p2 = p2->next; } while (p2 != cm);
            }
      }
    else
      {
        CliqueEl * p3 = newp;
        do
          {
            for (CliqueEl * p1 = cliques[p3->Nr()]; p1; p1 = p1->nextcl)
              {
                if (p1->clmaster != p1 || p1->eliminate) continue;
                bool dominated = true;
                CliqueEl * p2 = p1;
                do
                  {
                    if (!vertices[p2->Nr()].Flag()) { dominated = false; break; }
                    p2 = p2->next;
                  }
                while (p1 != p2);
                if (dominated)
                  {
                    p2 = p1;
                    do { p2->eliminate = 1; p2 = p2->next; } while (p1 != p2);
                  }
              }
            p3 = p3->next;
          }
        while (p3 != newp);
      }

    {
      // delete the dominated elements from the members' lists, the new
      // element becomes the head of every member's list
      CliqueEl * p3 = newp;
      do
        {
          p3->nextcl = cliques[p3->Nr()];
          while (p3->nextcl && p3->nextcl->eliminate)
            p3->nextcl = p3->nextcl->nextcl;

          CliqueEl hcel(-1);
          hcel.nextcl = cliques[p3->Nr()];
          CliqueEl * p1 = &hcel;
          while (p1)
            {
              while (p1->nextcl && p1->nextcl->eliminate)
                {
                  CliqueEl * hp = p1->nextcl;
                  p1->nextcl = p1->nextcl->nextcl;
                  ball.Free (hp);
                  vertices[p3->Nr()].numcliques--;
                }
              p1 = p1->nextcl;
            }
          
          cliques[p3->Nr()] = p3;
          p3 = p3->next;
          vertices[p3->Nr()].numcliques++;
        }
      while (p3 != newp);
    }

    {
      // delete the elements of v
      CliqueEl * p1 = cliques[v];
      while (p1)
        {
          CliqueEl * p2 = p1->GetNextClique();
          ball.Free (p1);
          p1 = p2;
        }
      cliques[v] = NULL;
    }

    {
      // prune the adjacency of the members: entries in the new element are
      // absorbed by it, eliminated and minion entries are stale
      CliqueEl * p3 = newp;
      do
        {
          int i = p3->Nr();
          auto a = Adj(i);
          int cnt = 0, wsum = 0;
          for (int w : a)
            if (!vertices[w].Eliminated() && IsMaster(w) && !vertices[w].Flag())
              { a[cnt++] = w; wsum += 1 + NumMinions(w); }
          adjlen[i] = cnt;
          adjw[i] = wsum;
          p3 = p3->next;
        }
      while (p3 != newp);
    }

    // indistinguishable vertices: only the new element and no adjacency
    // left -> merged into v; equal elements and adjacency -> merged
    CliqueEl * p3 = newp;
    bool merged;
    do
      {
        // a merge may release waiting vertices: next pass
        merged = false;
        do
          {
            if (IsMaster (p3->GetVertexNr()) && NumCliques (*p3) == 1 && adjlen[*p3] == 0 && nwait[*p3] == 0)
              {
                ReleaseChain (*p3);
                SetMaster (v, *p3);
                merged = true;
              }
            p3 = p3->next;
          }
        while (p3 != newp);
      }
    while (merged && level.Size());

    p3 = newp;
    do
      {
	if (IsMaster (p3->GetVertexNr()) && nwait[*p3] == 0)
	  {
	    int nclp3 = NumCliques (*p3), nadj3 = adjlen[*p3];
            SetFlagCliques (p3->GetVertexNr());
            markid++;
            for (int w : Adj(*p3)) vmark[w] = markid;

            for (CliqueEl * p4 = p3->next; p4 != newp; p4 = p4->next)
              if (IsMaster (*p4) && NumCliques(*p4) == nclp3 && adjlen[*p4] == nadj3 && nwait[*p4] == 0)
                { 
                  bool same = true;
                  for (CliqueEl * p1 = cliques[p4->Nr()]; p1; p1 = p1->nextcl)
                    if (!p1->Flag()) { same = false; break; }
                  if (same)
                    for (int w : Adj(*p4))
                      if (vmark[w] != markid) { same = false; break; }
                  if (same)
                    SetMaster (*p3, *p4);
                }

            ClearFlagCliques (p3->GetVertexNr());
          }
	p3 = p3->next;
      }
    while (p3 != newp);

    {
      // the elimination structure of v: its minions, then the members of
      // the new element with their minions
      int cnt = NumMinions(v);
      CliqueEl * p3 = newp;
      do 
        {
          if (IsMaster(*p3))
            cnt += 1+NumMinions(*p3);
          p3 = p3->next;
        }
      while (p3 != newp);

      delete [] vertices[v].connected;
      vertices[v].nconnected = cnt;
      vertices[v].connected = new int[cnt];
      cnt = 0;
      
      int hv = NextMinion(v);
      while (hv != -1)
        {
          vertices[v].connected[cnt++] = hv;
          hv = NextMinion (hv);
        }
      
      do 
        {
          if (IsMaster(*p3))
            {
              int hv = *p3;
              do
                {
                  vertices[v].connected[cnt++] = hv;
                  hv = NextMinion (hv);
                }
              while (hv != -1);
            }
          p3 = p3->next;
        }
      while (p3 != newp);
    }

    CliqueEl * anymaster = nullptr;
    {
      // disconnect the minions from all elements (their adjacency entries
      // elsewhere are skipped as stale)
      int cnt = 0;
      CliqueEl * p3 = newp;
      do 
        {
          if (IsMaster(*p3)) anymaster = p3;
          cnt++;
          p3 = p3->next;
        }
      while (p3 != newp);

      p3 = newp;
      for (int i = 0; i < cnt; i++)
        {
          int v3o = p3->GetVertexNr();
          p3 = p3->next;
          
          if (!IsMaster (v3o))
            {
              vertices[v3o].SetFlag (0);
              adjlen[v3o] = 0;
              for (CliqueEl * p1 = cliques[v3o]; p1; )
                {
                  if (p1->clmaster == p1)
                    {
                      CliqueEl * newmaster = p1->next;
                      newmaster->weight = p1->weight;
                      newmaster->w = p1->w;
                      newmaster->stamp = p1->stamp;
                      for (CliqueEl * p2 = p1->next; p2 != p1; p2=p2->next)
                        p2->clmaster = newmaster;
                      p1->clmaster = newmaster;
                    }

                  for (CliqueEl * p2 = p1; true; p2=p2->next)
                    if (p2->next == p1)  
                      {
                        p2->next = p2->next->next;
                        break;
                      }

                  CliqueEl * hp = p1;
                  p1 = p1->nextcl;
                  ball.Free (hp);
                }
              cliques[v3o] = nullptr;
            }
        }
    }
    
    // degrees of the remaining (master) members
    if (anymaster)
      {
        int wp = 0;
        CliqueEl * p3 = anymaster;
        do { wp += 1 + NumMinions (*p3); p3 = p3->next; } while (p3 != anymaster);
        anymaster->clmaster->weight = wp;

        if (approx)
          {
            // deg(i) <= 1 + |L_p \ i| + |A_i| + sum over the other elements of w(e),
            // and the old degree plus the growth is a bound too (AMD)
            p3 = anymaster;
            do
              {
                int v3 = *p3;
                int growth = wp - (1 + NumMinions (v3));
                int deg = 1 + growth + adjw[v3];
                for (CliqueEl * p1 = cliques[v3]; p1; p1 = p1->nextcl)
                  if (p1 != p3 && p1->clmaster->stamp == curstamp)
                    deg += p1->clmaster->w;
                int olddeg = GetDeg (v3);
                if (olddeg > 0 && olddeg < n) deg = min (deg, olddeg + growth);
                SetDeg (v3, min (deg, n-1));
                p3 = p3->next;
              }
            while (p3 != anymaster);
          }
        else
          {
            // exact: the new element counts once, the other elements and the
            // adjacency only outside it
            p3 = anymaster;
            do
              {
                int v3 = *p3;
                int deg = wp;
                for (CliqueEl * p1 = cliques[v3]; p1; p1 = p1->nextcl)
                  if (p1 != p3)
                    {
                      CliqueEl * p2 = p1;
                      do { vertices[*p2].SetUsed (0); p2 = p2->next; } while (p2 != p1);
                    }
                for (int w : Adj(v3)) vertices[w].SetUsed (0);
                for (CliqueEl * p1 = cliques[v3]; p1; p1 = p1->nextcl)
                  if (p1 != p3)
                    {
                      CliqueEl * p2 = p1;
                      do
                        {
                          int w = *p2;
                          if (!vertices[w].Flag() && !vertices[w].Used())
                            { deg += 1 + NumMinions (w); vertices[w].SetUsed (1); }
                          p2 = p2->next;
                        }
                      while (p2 != p1);
                    }
                // adjacency: pruned of the new element, but may overlap the other elements
                for (int w : Adj(v3))
                  if (!vertices[w].Used())
                    { deg += 1 + NumMinions (w); vertices[w].SetUsed (1); }
                SetDeg (v3, deg - NumMinions (v3));
                p3 = p3->next;
              }
            while (p3 != anymaster);
          }

        p3 = anymaster;
        do { vertices[*p3].SetFlag (0); p3 = p3->next; } while (p3 != anymaster);
      }
  }







  void MinimumDegreeOrdering :: EliminateMinionVertex (int v)
  {
    if (cliques[v]) // NumCliques(v) != 0)
      throw Exception ("Eliminate Minion should have exactly no clique");

    vertices[v].nconnected = 0;
    delete [] vertices[v].connected;    
    vertices[v].connected = NULL;
    adjlen[v] = 0;
  }







  void MinimumDegreeOrdering :: Order()
  {
    static Timer reorder_timer("MinimumDegreeOrdering::Order");
    RegionTimer reg(reorder_timer);

    cout << IM(4) << "start order" << endl;

    if (auto *tm = GetTaskManager(); tm) tm -> StopWorkers();

    // adjacency arrays from the collected edges
    adjstart.SetSize (n+1);
    adjlen.SetSize (n);
    adjw.SetSize (n);
    vmark.SetSize (n);
    adjstart = 0;
    for (size_t k = 0; k < edges.Size(); k += 2)
      { adjstart[edges[k]+1]++; adjstart[edges[k+1]+1]++; }

    int locked_dofs = 0;
    for (int i = 0; i < n; i++)
      if (vertices[i].Eliminated()) locked_dofs++;
    nused = n-locked_dofs;

    // dense vertices leave the graph, each costs O(n) per elimination
    densev.SetSize0();
    denseedges.SetSize0();
    if (densefac > 0)
      {
        size_t limit = max (16.0, densefac * sqrt (double(nused)));
        for (int i = 0; i < n; i++)
          if (adjstart[i+1] > limit) densev.Append (i);
        if (densev.Size() > 256) densev.SetSize0();   // too many
      }
    for (int d : lastdofs)
      if (!vertices[d].Eliminated() && !densev.Contains(d)) densev.Append (d);
    if (densev.Size() >= nused) densev.SetSize0();
    if (densev.Size())
      {
        vmark = 0;
        for (int d : densev) vmark[d] = 1;
        adjstart = 0;
        size_t keep = 0;
        for (size_t k = 0; k < edges.Size(); k += 2)
          {
            int v1 = edges[k], v2 = edges[k+1];
            if (vmark[v1] || vmark[v2])
              { denseedges.Append (v1); denseedges.Append (v2); }
            else
              {
                edges[keep++] = v1; edges[keep++] = v2;
                adjstart[v1+1]++; adjstart[v2+1]++;
              }
          }
        edges.SetSize (keep);
      }
    for (int i = 0; i < n; i++) adjstart[i+1] += adjstart[i];
    adjdata.SetSize (adjstart[n]);
    adjlen = 0;
    for (size_t k = 0; k < edges.Size(); k += 2)
      {
        int v1 = edges[k], v2 = edges[k+1];
        adjdata[adjstart[v1] + adjlen[v1]++] = v2;
        adjdata[adjstart[v2] + adjlen[v2]++] = v1;
      }
    edges.DeleteAll();
    vmark = 0;
    adjw = 0;

    nwait.SetSize (n);
    nwait = 0;
    relstart.SetSize (n+1);
    relstart = 0;
    if (level.Size())
      {
        // levels of the -1 vertices, breadth first
        int maxfixed = 0;
        for (int l : level) maxfixed = max (maxfixed, l);
        Array<int> front, next;
        for (int k = 0; ; k++)
          {
            if (k <= maxfixed)
              {
                front.SetSize0();
                for (int i = 0; i < n; i++)
                  if (level[i] == k) front.Append (i);
              }
            else if (!front.Size()) break;
            next.SetSize0();
            for (int v : front)
              for (int w : Adj(v))
                if (level[w] == -1) { level[w] = k+1; next.Append (w); }
            Swap (front, next);
          }
        for (int & l : level)
          if (l == -1) l = 0;

        for (int v = 0; v < n; v++)
          for (int w : Adj(v))
            if (level[w] > level[v]) { relstart[v+1]++; nwait[w]++; }
        if (waitfrac < 1)
          for (int & nw : nwait) nw = int (ceil (waitfrac * nw));
        for (int i = 0; i < n; i++) relstart[i+1] += relstart[i];
        reldata.SetSize (relstart[n]);
        for (int v = 0; v < n; v++)
          {
            size_t pos = relstart[v];
            for (int w : Adj(v))
              if (level[w] > level[v]) reldata[pos++] = w;
          }
        waitdeg.SetSize (n);
      }

    for (int j = 0; j < n; j++)
      SetDeg (j, 1+adjlen[j]);

    int minj = -1;
    int lastel = -1;

    if (n > 5000)
      cout << IM(4) << "order " << flush;

    for (int i = 0; i < n; i++)
      if (vertices[i].Eliminated())
        priqueue.SetDegree (i, n);
    for (int d : densev)
      priqueue.SetDegree (d, n);
    int nreg = nused - densev.Size();

    for (int i = 0; i < nreg; i++)
      {
	if (n > 5000 && i % 1000 == 999)
	  {
	    if (i % 10000 == 9999)
	      cout << IM(4) << "+" << flush;
	    else
	      cout << IM(4) << "." << flush;
	  }
      
	if (lastel != -1 && vertices[lastel].NextMinion() != -1)
	  {
	    minj = vertices[lastel].NextMinion();

	    if (vertices[minj].Eliminated())
	      cerr << "alread eliminated !!!" << endl;
	    priqueue.Invalidate(minj);

	    blocknr[i] = blocknr[i-1];
	    EliminateMinionVertex (minj);
	  }

	else
	  {
	    // find new master vertex
	    do
	      {
		minj = priqueue.MinDegree();
		priqueue.Invalidate(minj); 
		if (vertices[minj].Master() != minj)
		  priqueue.SetDegree (minj, n);
	      }
	    while (vertices[minj].Master() != minj);

	    blocknr[i] = i;
	    EliminateMasterVertex (minj);
	  }

	order[i] = minj;
	vertices[minj].SetEliminated (1);
	lastel = minj;
      }

    if (level.Size())
      {
        // inside a supernode: lower levels first. Root block: waiting vertices
        // with a neighbour outside go first (multiplier fixing a kernel)
        Array<int> pos(n);
        pos = -1;
        for (int i = 0; i < nreg; i++) pos[order[i]] = i;
        Array<unsigned char> outside(n);
        outside = 0;
        for (int i = 0; i < nreg; i++)
          for (int w : reldata.Range (relstart[order[i]], relstart[order[i]+1]))
            if (pos[w] >= 0 && blocknr[pos[w]] != blocknr[i]) outside[w] = 1;

        for (int i0 = 0; i0 < nreg; )
          {
            int i1 = i0+1;
            while (i1 < nreg && blocknr[i1] == i0) i1++;
            int oldfirst = order[i0];
            bool root = vertices[oldfirst].nconnected == i1-i0-1;
            auto key = [&] (int v) -> std::pair<int,int>
            {
              if (level[v] == 0) return { 1, 0 };
              return { root && outside[v] ? 0 : 2, level[v] };
            };
            auto block = order.Range (i0, i1);
            std::stable_sort (block.begin(), block.end(),
                              [&] (int a, int b) { return key(a) < key(b); });
            if (int newfirst = order[i0]; newfirst != oldfirst)
              {
                std::swap (vertices[newfirst].connected, vertices[oldfirst].connected);
                std::swap (vertices[newfirst].nconnected, vertices[oldfirst].nconnected);
                for (int j = 0; j < vertices[newfirst].nconnected; j++)
                  if (vertices[newfirst].connected[j] == newfirst)
                    vertices[newfirst].connected[j] = oldfirst;
              }
            i0 = i1;
          }
      }

    if (densev.Size())
      AppendDense (nreg);

    if (auto *tm = GetTaskManager(); tm) tm -> StartWorkers();
  }


  void MinimumDegreeOrdering :: AppendDense (int nreg)
  {
    static Timer t("MinimumDegreeOrdering::AppendDense");
    RegionTimer reg(t);
    int m = densev.Size();

    // dense vertices get position nreg+j
    Array<int> pos(n);
    pos = -1;
    for (int i = 0; i < nreg; i++) pos[order[i]] = i;
    for (int j = 0; j < m; j++) pos[densev[j]] = nreg+j;

    // parent in the elimination tree: block of the first external vertex
    Array<int> bstart, parent(nreg);
    for (int i0 = 0; i0 < nreg; )
      {
        int i1 = i0+1;
        while (i1 < nreg && blocknr[i1] == i0) i1++;
        bstart.Append (i0);
        MDOVertex & v = vertices[order[i0]];
        int minpos = nreg;
        for (int j = 0; j < v.nconnected; j++)
          if (int p = pos[v.connected[j]]; p >= i1) minpos = min (minpos, p);
        parent[i0] = minpos < nreg ? blocknr[minpos] : -1;
        i0 = i1;
      }

    // structure of the dense vertices: adjacent blocks and their ancestors
    int nchunks = (m+63)/64;
    Array<uint64_t> mask(size_t(nchunks)*nreg);
    mask = 0;
    Array<int> firstnb(m);
    firstnb = nreg;
    for (size_t k = 0; k < denseedges.Size(); k += 2)
      {
        int p1 = pos[denseedges[k]], p2 = pos[denseedges[k+1]];
        if (p1 > p2) Swap (p1, p2);
        if (p1 < 0 || p1 >= nreg || p2 < nreg) continue;
        int j = p2-nreg;
        mask[size_t(j/64)*nreg + blocknr[p1]] |= uint64_t(1) << (j%64);
        firstnb[j] = min (firstnb[j], p1);
      }
    for (int c = 0; c < nchunks; c++)
      {
        auto cm = mask.Range (size_t(c)*nreg, size_t(c+1)*nreg);
        for (int b : bstart)
          if (cm[b] && parent[b] >= 0) cm[parent[b]] |= cm[b];
      }
    for (int b : bstart)
      {
        int cnt = 0;
        for (int c = 0; c < nchunks; c++)
          cnt += std::popcount (mask[size_t(c)*nreg + b]);
        if (!cnt) continue;
        MDOVertex & v = vertices[order[b]];
        int * newc = new int[v.nconnected + cnt];
        for (int j = 0; j < v.nconnected; j++) newc[j] = v.connected[j];
        delete [] v.connected;
        v.connected = newc;
        for (int j = 0; j < m; j++)
          if (mask[size_t(j/64)*nreg + b] >> (j%64) & 1)
            v.connected[v.nconnected++] = densev[j];
      }

    // final block: dense vertices, then the coupled root block. A constraint
    // fixing a kernel (mean value, rigid body motions) must not come last.
    // Several roots: only their last vertices of each level
    Array<int> ntail(nreg);
    ntail = 0;
    size_t sumtail = 0;
    int nroots = 0;
    for (size_t ib = 0; ib < bstart.Size(); ib++)
      if (int b = bstart[ib]; parent[b] == -1)
        {
          int cnt = 0;
          for (int c = 0; c < nchunks; c++)
            cnt += std::popcount (mask[size_t(c)*nreg + b]);
          if (!cnt) continue;
          int bend = ib+1 < bstart.Size() ? bstart[ib+1] : nreg;
          ntail[b] = min (2*cnt, bend-b);
          sumtail += ntail[b];
          nroots++;
        }
    if (nroots == 1)
      {
        for (int & nt : ntail)
          if (nt) nt = nreg;
      }
    else if (sumtail > size_t (max (64, 4*m)))
      ntail = 0;               // many components: plain append

    Array<int> neworder(nused), newblocknr(nused), tail;
    Array<unsigned char> intail(nreg);
    intail = 0;
    int cnt = 0;
    for (size_t ib = 0; ib < bstart.Size(); ib++)
      {
        int b = bstart[ib];
        int bend = ib+1 < bstart.Size() ? bstart[ib+1] : nreg;
        int nb = cnt;
        if (ntail[b])
          for (int i = bend-1, curlevel = -1, taken = 0; i >= b; i--)
            {
              int l = level.Size() ? level[order[i]] : 0;
              if (l != curlevel) { curlevel = l; taken = 0; }
              if (taken < ntail[b]) { intail[i] = 1; taken++; }
            }
        int oldfirst = order[b], newfirst = -1;
        for (int i = b; i < bend; i++)
          if (intail[i])
            tail.Append (order[i]);
          else
            {
              if (newfirst == -1) newfirst = order[i];
              neworder[cnt] = order[i]; newblocknr[cnt++] = nb;
            }
        if (newfirst == -1)
          {
            MDOVertex & v = vertices[oldfirst];
            delete [] v.connected;
            v.connected = nullptr;
            v.nconnected = 0;
          }
        else if (newfirst != oldfirst)
          {
            std::swap (vertices[newfirst].connected, vertices[oldfirst].connected);
            std::swap (vertices[newfirst].nconnected, vertices[oldfirst].nconnected);
            for (int j = 0; j < vertices[newfirst].nconnected; j++)
              if (vertices[newfirst].connected[j] == newfirst)
                vertices[newfirst].connected[j] = oldfirst;
          }
      }
    int f0 = cnt;
    // dense vertices coupling only to the tails go to the very end
    Array<unsigned char> before(m);
    for (int j = 0; j < m; j++) before[j] = firstnb[j] == nreg;
    for (size_t k = 0; k < denseedges.Size(); k += 2)
      {
        int p1 = pos[denseedges[k]], p2 = pos[denseedges[k+1]];
        if (p1 > p2) Swap (p1, p2);
        if (p1 >= 0 && p1 < nreg && p2 >= nreg && !intail[p1]) before[p2-nreg] = 1;
      }
    Array<int> atend;
    for (int j = 0; j < m; j++)
      if (before[j]) neworder[cnt++] = densev[j]; else atend.Append (densev[j]);
    for (int v : tail) neworder[cnt++] = v;
    for (int v : atend) neworder[cnt++] = v;
    for (int i = f0; i < nused; i++) newblocknr[i] = f0;

    order.Range(0, nused) = neworder;
    blocknr.Range(0, nused) = newblocknr;

    for (int d : densev) vertices[d].nconnected = 0;
    MDOVertex & first = vertices[order[f0]];
    delete [] first.connected;
    first.nconnected = nused-f0-1;
    first.connected = new int[first.nconnected];
    for (int i = f0+1; i < nused; i++)
      first.connected[i-f0-1] = order[i];
  }



  MinimumDegreeOrdering:: ~MinimumDegreeOrdering ()
  {
    // cout << "~MDO: all data should be deleted, please double-check" << endl;
    for (int i = 0; i < vertices.Size(); i++)
      delete [] vertices[i].connected;
  }


  void MinimumDegreeOrdering :: SetFlagNodes (int v)
  {
    for (CliqueEl * p1 = cliques[v]; p1; p1 = p1->nextcl)
      {
	CliqueEl * p2 = p1;
	do
	  {
	    vertices[p2->GetVertexNr()].SetFlag (1);
	    p2 = p2->next;
	  }
	while (p2 != p1);
      }
  }

  void MinimumDegreeOrdering :: ClearFlagNodes (int v)
  {
    for (CliqueEl * p1 = cliques[v]; p1; p1 = p1->nextcl)
      {
	CliqueEl * p2 = p1;
	do
	  {
	    vertices[p2->GetVertexNr()].SetFlag (0);
	    p2 = p2->next;
	  }
	while (p2 != p1);
      }
  }

  void MinimumDegreeOrdering :: SetFlagCliques (int v)
  {
    for (CliqueEl * p1 = cliques[v]; p1; p1 = p1->nextcl)
      p1->SetFlag (true);
  }


  void MinimumDegreeOrdering :: ClearFlagCliques (int v)
  {
    for (CliqueEl * p1 = cliques[v]; p1; p1 = p1->nextcl)
      p1->SetFlag (false);
  }


  void MinimumDegreeOrdering :: SetMaster (int master, int minion)
  {
    int hv = master;
    while (vertices[hv].NextMinion() != -1)
      hv = vertices[hv].NextMinion();
			    
    vertices[hv].SetNextMinion (minion);
    while (hv != -1)
      {
	vertices[hv].SetMaster (master);
	hv = vertices[hv].NextMinion();
      }    

    vertices[master].numminions += 1+vertices[minion].numminions;
    priqueue.SetDegree(minion, n);
  }




  MDOPriorityQueue :: MDOPriorityQueue (int size, int maxdeg)
    : list(size), first_in_class(maxdeg)
  {
    /*
    for (int i = 0; i < size; i++)
      list[i].degree = 0;
    for (int i = 0; i < maxdeg; i++)
      first_in_class[i] = -1;
    */
    ParallelFor(size, [&](size_t i)
                { this->list[i].degree = 0; });
    ParallelFor(maxdeg, [&](size_t i)
                { first_in_class[i] = -1; });
  }

  MDOPriorityQueue :: ~MDOPriorityQueue ()
  {
    ;
  }

  void MDOPriorityQueue::DoArchive(Archive &ar)
  {
    ar & list & first_in_class;
  }

  int MDOPriorityQueue :: MinDegree () const
  {
    for (int i = 0; i < first_in_class.Size(); i++)
      {
	if (first_in_class[i] != -1)
	  return first_in_class[i];
      }
    return 0;
  }

  void MDOPriorityQueue :: SetDegree (int nr, int deg)
  {
    if (deg == 0)
      deg++;

    if (list[nr].degree > 0)
      Invalidate (nr);

    if (first_in_class[deg] != -1)
      {
	int next = first_in_class[deg];
	int prev = list[next].prev;

	list[nr].next = next;
	list[nr].prev = prev;
	list[next].prev = nr;
	list[prev].next = nr;
      }
    else
      {
	list[nr].next = nr;
	list[nr].prev = nr;
	first_in_class[deg] = nr;
      }
    list[nr].degree = deg;
  }

  void MDOPriorityQueue :: Invalidate (int nr)
  {
    if (!list[nr].degree)
      cerr << "already eliminated" << endl;

    if (list[nr].next == nr)
      { // just one element in class
	first_in_class[list[nr].degree] = -1;
      }
    else
      {
	int next = list[nr].next;
	int prev = list[nr].prev;
      
	list[prev].next = next;
	list[next].prev = prev;
	first_in_class[list[nr].degree] = next;
      }
    list[nr].degree = 0;
  }

}
