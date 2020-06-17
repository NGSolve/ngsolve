/* *************************************************************************/
/* File:   order.cc                                                        */
/* Author: Joachim Schoeberl                                               */
/* Date:   25. Mar. 2000                                                   */
/* *************************************************************************/



#include <la.hpp>
namespace ngla
{
  using namespace ngla;


  /* 

  Finds ordering for sparse Cholesky Factorization
  Algorithm: Minimum Degree
              
  See:       The evolution of the 
  minimum degree ordering algorithm

  A. George and J.W.H. Liu
  SIAM Review, Vol31, 1989, pp1-19

  */


  // ngstd::BlockAllocator CliqueEl :: ball(sizeof (CliqueEl));


  
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



  void MinimumDegreeOrdering :: AddEdge (int v1, int v2)
  {
    if (v1 == v2) return;

    CliqueEl * p1 = new (ball) CliqueEl(v1);
    CliqueEl * p2 = new (ball) CliqueEl(v2);

    p1->eliminate = 0;
    p2->eliminate = 0;

    p1->next = p2;
    p2->next = p1;

    p1->clmaster = p1;
    p2->clmaster = p1;

    p1->SetFlag (false);
    p2->SetFlag (false);

    p1->nextcl = cliques[v1];
    cliques[v1] = p1;
    p2->nextcl = cliques[v2];
    cliques[v2] = p2;

    vertices[v1].numcliques++;
    vertices[v2].numcliques++;
  }


  void MinimumDegreeOrdering :: PrintCliques () 
  {
    for (int i = 0; i < n; i++)
      if (!vertices[i].Eliminated())
	{
	  (*testout) << "Vertex " << i << ", degree = " 
		     << CalcDegree (i) 
		     << endl;
	  
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
    // static Timer t("MDO::CalcDegree");
    // RegionTimer reg(t);

    int deg = 0;

    // clear flags to count vertices just once:
    for (CliqueEl * p1 = cliques[v1]; p1; p1 = p1->nextcl)
      {
	CliqueEl * p2 = p1;
	do
	  {
	    vertices[p2->GetVertexNr()].SetUsed(0);   
	    p2 = p2->next;
	  }
	while (p2 != p1);
      }

    
    // count members of all cliques of vertex
    
    for (CliqueEl * p1 = cliques[v1]; p1; p1 = p1->nextcl)
      {
	CliqueEl * p2 = p1;
	do
	  {
	    int v2 = p2->GetVertexNr();
	    if (!vertices[v2].Used())
	      {
		if (IsMaster(v2)) 
                  deg += 1+NumMinions(v2);
                else
                  cerr << "we still have minions" << endl;
		vertices[p2->GetVertexNr()].SetUsed(1);
	      }
	    p2 = p2->next;
	  }
	while (p2 != p1);
      }

    return deg;
  }
 


  void MinimumDegreeOrdering :: EliminateMasterVertex (int v)
  {
    static Timer t("MDO::EliminateMaster", 2);
    /*
    static Timer t1("MDO::EliminateMaster 1", 2);
    static Timer t2("MDO::EliminateMaster 2", 2);
    static Timer t2a("MDO::EliminateMaster 2a", 2);
    static Timer t2b("MDO::EliminateMaster 2b", 2);
    static Timer t2c("MDO::EliminateMaster 2c", 2);
    static Timer t3("MDO::EliminateMaster 3", 2);
    static Timer t4("MDO::EliminateMaster 4 (calcdeg)", 2);
    */
    RegionTimer reg(t);

    // t1.Start();
    // (*testout) << "Eliminate Master Vertex " << v  << endl;

    // int numminions = NumMinions (v);

    for (CliqueEl * p1 = cliques[v]; p1; p1 = p1->nextcl)
      {
	CliqueEl * p2 = p1;
	do
	  {
	    vertices[p2->GetVertexNr()].SetUsed (0);
	    p2->eliminate = 1;
	    p2 = p2->next;
	  }
	while (p2 != p1);
      }


    //  new clique is union of cliques of vertex v
    CliqueEl * newp = NULL;

    for (CliqueEl * p1 = cliques[v]; p1; p1 = p1->nextcl)
      {
	CliqueEl * p2 = p1;
	do
	  {
	    if (!vertices[p2->GetVertexNr()].Used() && p2->GetVertexNr() != v)
	      {
                CliqueEl * p3 = new (ball) CliqueEl (p2->GetVertexNr());

		p3 -> next = newp;
                p3 -> clmaster = newp ? newp->clmaster : p3;

		p3 -> eliminate = 0;
		p3 -> nextcl = NULL;
		p3 -> SetFlag (false);

		newp = p3;

		vertices[p2->GetVertexNr()].SetUsed(1);
	      }
	    p2 = p2->next;
	  }
	while (p2 != p1);
      }

    // t1.Stop();


    if (!newp)
      {
        // was not in any clique with more than 1 member
	CliqueEl * p1 = cliques[v];
	while (p1)
	  {
	    CliqueEl * p2 = p1->GetNextClique();
	    ball.Free(p1);
	    p1 = p2;
	  }
	cliques[v] = NULL;	
	vertices[v].nconnected = 0;
        delete [] vertices[v].connected;
	vertices[v].connected = NULL;
	
	return;
      }


	
    // close new clique
    newp -> clmaster -> next = newp;
      
    // t2.Start();
    // t2a.Start();


    // find dominated cliques

    // flag all vertices belonging to new clique
    {
      CliqueEl * p3 = newp;
      do
        {
          vertices[p3->Nr()].SetFlag (1);
          p3 = p3->next;
        }
      while (p3 != newp);
    }

    

    {
      // check all cliques of all members of new clique
      CliqueEl * p3 = newp;
      do
        {
          for (CliqueEl * p1 = cliques[p3->Nr()]; p1; p1 = p1->nextcl)
            {
              if (p1->clmaster != p1) continue;

              if (!p1->eliminate)
                {
                  bool dominated = 1;
		  
                  CliqueEl * p2 = p1;
                  do
                    {
                      if (!vertices[p2->Nr()].Flag())
                        {
                          dominated = 0;
                          break;
                        }
                      p2 = p2->next;
                    }
                  while (p1 != p2);
                  
                  // if clique is dominated, set eliminate flag
                  if (dominated)
                    {
                      CliqueEl * p2 = p1;
                      do
                        {
                          p2->eliminate = 1;
                          p2 = p2->next;
                        }
                      while (p1 != p2);
                    }
                } 
            }
          p3 = p3->next;
        }
      while (p3 != newp);
    }
    
    // t2a.Stop();
    // t2b.Start();
    
    {
      // clear flag
      CliqueEl * p3 = newp;
      do
        {
          vertices[p3->GetVertexNr()].SetFlag (0);
          p3 = p3->next;
        }
      while (p3 != newp);
    }
    
  
    {
      // delete cliques with eliminate flag set
      CliqueEl * p3 = newp;
      do
        {
          // p3->next is first remaining clique:
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
      // delete also clique-elements of v
      CliqueEl * p1 = cliques[v];
      while (p1)
        {
          CliqueEl * p2 = p1->GetNextClique();
          ball.Free (p1);
          p1 = p2;
        }
      cliques[v] = NULL;
    }
    
    // t2b.Stop();
    // t2c.Start();


    // find equivalent cliques
    CliqueEl * p3 = newp;
    do
      {
	if (IsMaster (p3->GetVertexNr()))
	  {
	    int nclp3 = NumCliques (*p3);
	    if ( nclp3 == 1)
	      {
		// only in new clique ==> connect to v
		SetMaster (v, *p3);
	      }
	  }
        
	p3 = p3->next;
      }
    while (p3 != newp);


    p3 = newp;
    do
      {
	if (IsMaster (p3->GetVertexNr()))
	  {
	    int nclp3 = NumCliques (*p3);

            SetFlagCliques (p3->GetVertexNr());

            for (CliqueEl * p4 = p3->next; p4 != newp; p4 = p4->next)
              {
                // have p3 and p4 equivalent cliques ?
                if (IsMaster (*p4) && NumCliques(*p4) == nclp3)
                  { 
                    bool samecl = true;

                    for (CliqueEl * p1 = cliques[p4->Nr()]; p1; p1 = p1->nextcl)
                      if (!p1->Flag()) 
                        {
                          samecl = false;
                          break;
                        }
                    
                    if (samecl)
                      SetMaster (*p3, *p4);
                  }
              }

            ClearFlagCliques (p3->GetVertexNr());
          }

	p3 = p3->next;
      }
    while (p3 != newp);

    // t2c.Stop();
    // t2.Stop();


    // t3.Start();

    {
      // setup elimination data structures for vertex v

      int cnt = NumMinions(v);
      CliqueEl * p3 = newp;
      do 
        {
          if (IsMaster(*p3))
            cnt += 1+NumMinions(*p3);
          p3 = p3->next;
        }
      while (p3 != newp);

      delete [] vertices[v].connected;  // necessary ?       
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
      // disconnect minions 
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
              for (CliqueEl * p1 = cliques[v3o]; p1; )
                {
                  // disconnect p1
                  if (p1->clmaster == p1)
                    {
                      CliqueEl * newmaster = p1->next;
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
      // t3.Stop();
    }
    
    // t4.Start();
    // calc master degrees in new clique
    if (anymaster)
      {
        CliqueEl * p3 = anymaster;
        do
          {
            priqueue.SetDegree (*p3, CalcDegree (*p3) - NumMinions (*p3));
            p3 = p3->next;
          }
        while (p3 != anymaster);      
      }
    // t4.Stop();
  }







  void MinimumDegreeOrdering :: EliminateMinionVertex (int v)
  {
    if (cliques[v]) // NumCliques(v) != 0)
      throw Exception ("Eliminate Minion should have exactly no clique");

    vertices[v].nconnected = 0;
    delete [] vertices[v].connected;    
    vertices[v].connected = NULL;
  }







  void MinimumDegreeOrdering :: Order()
  {
    static Timer reorder_timer("MinimumDegreeOrdering::Order");
    RegionTimer reg(reorder_timer);

    cout << IM(4) << "start order" << endl;

    if (task_manager) task_manager -> StopWorkers();

    for (int j = 0; j < n; j++)
      {
	// priqueue.SetDegree(j, CalcDegree(j));
	priqueue.SetDegree(j, 1+NumCliques(j));
      }

    int minj = -1;
    int lastel = -1;

    if (n > 5000)
      cout << IM(4) << "order " << flush;

    int locked_dofs = 0;
    for (int i = 0; i < n; i++)
      if (vertices[i].Eliminated())                
        {
          locked_dofs++;
          priqueue.SetDegree (i, n);                  
        }
    nused = n-locked_dofs;

    for (int i = 0; i < nused; i++)
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
    // PrintCliques();
    // NgProfiler::Print (cout);
    if (task_manager) task_manager -> StartWorkers();
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
