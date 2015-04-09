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
    cliques = NULL;
    blocknr = 0;
    order = 0;

    for (int i = 0; i < n; i++)
      vertices[i].Init(i);
  }



  void MinimumDegreeOrdering :: AddEdge (int v1, int v2)
  {
    // *testout << "Add edge " << v1 << " - " << v2 << endl;
    if (v1 == v2) return;

    CliqueEl *p1, *p2;
    // p1 = new CliqueEl;
    // p2 = new CliqueEl;
    
    p1 = (CliqueEl*)ball.Alloc();
    p2 = (CliqueEl*)ball.Alloc();
  
    p1->vnr = v1;
    p2->vnr = v2;

    p1->eliminate = 0;
    p2->eliminate = 0;

    p1->next = p2;
    p2->next = p1;

    p1->flag = 0;
    p2->flag = 0;

    p1->nextcl = cliques[v1];
    cliques[v1] = p1;
    p2->nextcl = cliques[v2];
    cliques[v2] = p2;
  }


  void MinimumDegreeOrdering :: PrintCliques () 
  {
    CliqueEl *p1, *p2;

    for (int i = 0; i < n; i++)
      if (!vertices[i].Eliminated())
	{
	  (*testout) << "Vertex " << i << ", degree = " 
		     << CalcDegree (i) 
		     << endl;
	  
	  p1 = cliques[i];
	  while (p1)
	    {
	      p2 = p1;
	      (*testout) << "( ";
	      do
		{
		  if (!vertices[p2->GetVertexNr()].Eliminated())
		    (*testout) << p2->GetVertexNr() << " ";
		  p2 = p2->next;
		}
	      while (p2 != p1);
	      (*testout) << ")";
	      p1 = p1->nextcl;
	    }
	  (*testout) << endl;
	}
  }





  int MinimumDegreeOrdering :: CalcDegree (int v1)
  {
    static Timer t("MDO::CalcDegree");
    RegionTimer reg(t);

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
                  deg += 1+NumSlaves(v2);
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
    static Timer t("MDO::EliminateMaster (incl calcdeg)");
    static Timer t1("MDO::EliminateMaster 1");
    static Timer t2("MDO::EliminateMaster 2");
    static Timer t2a("MDO::EliminateMaster 2a");
    static Timer t2b("MDO::EliminateMaster 2b");
    static Timer t2c("MDO::EliminateMaster 2c");
    static Timer t2ca("MDO::EliminateMaster 2ca");
    static Timer t2cb("MDO::EliminateMaster 2cb");
    static Timer t2cc("MDO::EliminateMaster 2cc");
    static Timer t3("MDO::EliminateMaster 3");
    RegionTimer reg(t);

    t1.Start();
    // (*testout) << "Eliminate Master Vertex " << v  << endl;

    // int numslaves = NumSlaves (v);

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
		CliqueEl * p3 = (CliqueEl*)ball.Alloc();

		p3 -> next = newp;
		p3 -> vnr = p2->GetVertexNr();
		p3 -> eliminate = 0;
	      
		p3 -> nextcl = NULL;
		p3 -> flag = 0;

		newp = p3;

		vertices[p2->GetVertexNr()].SetUsed(1);
	      }
	    p2 = p2->next;
	  }
	while (p2 != p1);
      }



    if (!newp)
      {
	CliqueEl * p1 = cliques[v];
	while (p1)
	  {
	    CliqueEl * p2 = p1->GetNextClique();
	    // delete p1;
	    ball.Free(p1);
	    p1 = p2;
	  }
	cliques[v] = NULL;	
	vertices[v].nconnected = 0;
	vertices[v].connected = NULL;
	
	return;
      }


	
    // close new clique
    {
      CliqueEl * p3 = newp;
      while (p3->next) p3 = p3->next;
      p3->next = newp;
    }

    t1.Stop();
    t2.Start();
    t2a.Start();



    // find dominated cliques

    // set flag for new clique
    {
      CliqueEl * p3 = newp;
      do
        {
          vertices[p3->GetVertexNr()].SetFlag (1);
          p3 = p3->next;
        }
      while (p3 != newp);
    }

    

    {
      // check all cliques of all members of new clique
      CliqueEl * p3 = newp;
      do
        {
          for (CliqueEl * p1 = cliques[p3->GetVertexNr()]; p1; p1 = p1->nextcl)
            {
              if (!p1->eliminate)
                {
                  bool dominated = 1;
		  
                  CliqueEl * p2 = p1;
                  do
                    {
                      if (!vertices[p2->GetVertexNr()].Flag())
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
    
    t2a.Stop();
    t2b.Start();
    
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
      // delete old cliques
      CliqueEl * p3 = newp;
      do
        {
          // p3->next is first remaining clique:
          
          p3->nextcl = cliques[p3->GetVertexNr()];
          while (p3->nextcl && p3->nextcl->eliminate)
            p3->nextcl = p3->nextcl->nextcl;
          
          
          CliqueEl hcel;
          hcel.nextcl = cliques[p3->GetVertexNr()];
          CliqueEl * p1 = &hcel;
	  
          while (p1)
            {
              while (p1->nextcl && p1->nextcl->eliminate)
                {
                  CliqueEl * hp = p1->nextcl;
                  p1->nextcl = p1->nextcl->nextcl;
                  ball.Free (hp);
                }
              p1 = p1->nextcl;
            }
          
          cliques[p3->GetVertexNr()] = p3;
          
          p3 = p3->next;
        }
      while (p3 != newp);
    }
      

    {
      
      CliqueEl * p1 = cliques[v];
      while (p1)
        {
          CliqueEl * p2 = p1->GetNextClique();
          ball.Free (p1);
          p1 = p2;
        }
      
      cliques[v] = NULL;
    }
    
    t2b.Stop();
    t2c.Start();

    // find connections
    CliqueEl * p3 = newp;
    do
      {
	// if (vertices[p3->GetVertexNr()].Master() == p3->GetVertexNr())
	if (IsMaster (p3->GetVertexNr()))
	  {
	    int nclp3 = NumCliques (p3->GetVertexNr());
	    if ( nclp3 == 1)
	      {
		// only in new clique, connect to v
		SetMaster (v, p3->GetVertexNr());
	      }

	    else

	      {
                t2ca.Start();
		SetFlagCliques (p3->GetVertexNr());
                t2ca.Stop();
                t2cb.Start();
		for (CliqueEl * p4 = p3->next; p4 != newp; p4 = p4->next)
		  {
		    // are p3 and p4 connected ?
		    
		    if (IsMaster (p4->GetVertexNr()))
		      {
			bool samecl = 1;
                        int nclp4 = 0;
                        for (CliqueEl * p1 = cliques[p4->GetVertexNr()]; p1; p1 = p1->nextcl)
                          {
                            if (!p1->flag) 
                              {
                                samecl = 0;
                                break;
                              }
                            nclp4++;
                          }

                        if (nclp4 != nclp3) samecl = 0;

			if (samecl)
			  SetMaster (p3->GetVertexNr(), p4->GetVertexNr());
		      }
		  }
                t2cb.Stop();
                t2cc.Start();
		ClearFlagCliques (p3->GetVertexNr());
                t2cc.Stop();
	      }
	  }

	p3 = p3->next;
      }
    while (p3 != newp);

    t2c.Stop();
    t2.Stop();

    // calc master degrees in new clique
    {
      CliqueEl * p3 = newp;
      do
        {
          int v3 = p3->GetVertexNr();
          if (IsMaster(v3))
            priqueue.SetDegree (v3, CalcDegree (v3) - NumSlaves (v3));
          p3 = p3->next;
        }
      while (p3 != newp);      
    }
    

    t3.Start();

    {
      int cnt = NumSlaves(v);
      CliqueEl * p3 = newp;
      do 
        {
          int v3 = p3->GetVertexNr();
          if (IsMaster(v3))
            cnt += 1+NumSlaves(v3);
          p3 = p3->next;
        }
      while (p3 != newp);
      
      vertices[v].nconnected = cnt;
      vertices[v].connected = new int[cnt];
      cnt = 0;
      
      int hv = NextSlave(v);
      while (hv != -1)
        {
          vertices[v].connected[cnt++] = hv;
          hv = NextSlave (hv);
        }
      
      do 
        {
          int v3 = p3->GetVertexNr();
          if (IsMaster(v3))
            {
              int hv = v3;
              do
                {
                  vertices[v].connected[cnt++] = hv;
                  hv = NextSlave (hv);
                }
              while (hv != -1);
            }
          p3 = p3->next;
        }
      while (p3 != newp);
    }
      


    // disconnect slaves 
    {
      int cnt = 0;
      CliqueEl * p3 = newp;
      do 
        {
          cnt++;
          p3 = p3->next;
        }
      while (p3 != newp);


      p3 = newp;
      for (int i = 0; i < cnt; i++)
        {
          int v3 = p3->GetVertexNr();
          
          p3 = p3->next;
          
          if (!IsMaster (v3))
            {
              for (CliqueEl * p1 = cliques[v3]; p1; )
                {
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
              cliques[v3] = 0;
            }
        }
      t3.Stop();
    }
  }







  void MinimumDegreeOrdering :: EliminateSlaveVertex (int v)
  {
    if (cliques[v]) // NumCliques(v) != 0)
      throw Exception ("Eliminate Slave should have exactly no clique");

    vertices[v].nconnected = 0; 
    vertices[v].connected = NULL;
  }







  void MinimumDegreeOrdering :: Order()
  {
    static Timer reorder_timer("MinimumDegreeOrdering::Order");
    RegionTimer reg(reorder_timer);

    cout << IM(4) << "start order" << endl;

    for (int j = 0; j < n; j++)
      {
	// priqueue.SetDegree(j, CalcDegree(j));
	priqueue.SetDegree(j, 1+NumCliques(j));
      }

    int minj = -1;
    int lastel = -1;

    if (n > 5000)
      cout << IM(4) << "order " << flush;

    for (int i = 0; i < n; i++)
      {
	if (n > 5000 && i % 1000 == 999)
	  {
	    if (i % 10000 == 9999)
	      cout << IM(4) << "+" << flush;
	    else
	      cout << IM(4) << "." << flush;
	    //	    cout << "allocated els = " << ball.NumElements() << endl;

	    // NgProfiler::Print (cout);
	  }
      
	if (lastel != -1 && vertices[lastel].NextSlave() != -1)
	  {
	    minj = vertices[lastel].NextSlave();

	    if (vertices[minj].Eliminated())
	      cerr << "alread eliminated !!!" << endl;
	    priqueue.Invalidate(minj);

	    blocknr[i] = blocknr[i-1];
	    EliminateSlaveVertex (minj);
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
  }



  MinimumDegreeOrdering:: ~MinimumDegreeOrdering ()
  {
    for (int i = 0; i < vertices.Size(); i++)
      delete [] vertices[i].connected;
  }




  int MinimumDegreeOrdering :: NumCliques (int v) const
  {
    static Timer t("MinimumDegreeOrdering::NumCliques");
    RegionTimer reg(t);

    int cnt = 0;
    for (CliqueEl * p1 = cliques[v]; p1; p1 = p1->nextcl)
      cnt++;
    return cnt;
  }


  void MinimumDegreeOrdering :: SetFlagNodes (int v)
  {
    CliqueEl * p1 = cliques[v];
    while (p1)
      {
	CliqueEl * p2 = p1;
	do
	  {
	    vertices[p2->GetVertexNr()].SetFlag (1);
	    p2 = p2->next;
	  }
	while (p2 != p1);
	p1 = p1->nextcl;
      }
  }

  void MinimumDegreeOrdering :: ClearFlagNodes (int v)
  {
    CliqueEl * p1 = cliques[v];
    while (p1)
      {
	CliqueEl * p2 = p1;
	do
	  {
	    vertices[p2->GetVertexNr()].SetFlag (0);
	    p2 = p2->next;
	  }
	while (p2 != p1);
	p1 = p1->nextcl;
      }
  }

  void MinimumDegreeOrdering :: SetFlagCliques (int v)
  {
    for (CliqueEl * p1 = cliques[v]; p1; p1 = p1->nextcl)
      {
	CliqueEl * p2 = p1;
	do
	  {
	    p2->flag = 1;
	    p2 = p2->next;
	  }
	while (p2 != p1);
      }
  }

  void MinimumDegreeOrdering :: ClearFlagCliques (int v)
  {
    for (CliqueEl * p1 = cliques[v]; p1; p1 = p1->nextcl)
      {
	CliqueEl * p2 = p1;
	do
	  {
	    p2->flag = 0;
	    p2 = p2->next;
	  }
	while (p2 != p1);
      }
  }



  void MinimumDegreeOrdering :: SetMaster (int master, int slave)
  {
    int hv = master;
    while (vertices[hv].NextSlave() != -1)
      hv = vertices[hv].NextSlave();
			    
    vertices[hv].SetNextSlave (slave);
    while (hv != -1)
      {
	vertices[hv].SetMaster (master);
	hv = vertices[hv].NextSlave();
      }    

    vertices[master].numslaves += 1+vertices[slave].numslaves;
    priqueue.SetDegree(slave, n);
  }




  MDOPriorityQueue :: MDOPriorityQueue (int size, int maxdeg)
    : list(size), first_in_class(maxdeg)
  {
    for (int i = 0; i < size; i++)
      list[i].degree = 0;
    for (int i = 0; i < maxdeg; i++)
      first_in_class[i] = -1;
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
