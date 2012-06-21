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
    CliqueEl *p1, *p2;
    int deg = 0;

    // clear flags to count vertices just once:
    p1 = cliques[v1];
    while (p1)
      {
	p2 = p1;
	do
	  {
	    vertices[p2->GetVertexNr()].SetUsed(0);   
	    p2 = p2->next;
	  }
	while (p2 != p1);
	p1 = p1->nextcl;
      }


    // count members of all cliques of vertex
    p1 = cliques[v1];
    while (p1)
      {
	p2 = p1;
	do
	  {
	    int v2 = p2->GetVertexNr();
	    if (!vertices[v2].Used())
	      {
		if (vertices[v2].Master() == v2)
		  {
		    deg += 1+NumSlaves(v2);
		  }
		vertices[p2->GetVertexNr()].SetUsed(1);
	      }
	    p2 = p2->next;
	  }
	while (p2 != p1);
	p1 = p1->nextcl;
      }

    return deg;
  }
 


  void MinimumDegreeOrdering :: EliminateMasterVertex (int v)
  {
    // (*testout) << "Eliminate Master Vertex " << v  << endl;

    /*
    if (v == 361611)
      PrintCliques();
    if (v == 397582)
      PrintCliques();
    */

    // int numslaves = NumSlaves (v);


    // clear used-field
    /*
      CliqueEl * p1 = cliques[v];
      while (p1)
      {
      CliqueEl * p2 = p1;
      do
      {
      vertices[p2->GetVertexNr()].SetUsed (0);
      p2->eliminate = 1;
      p2 = p2->next;
      }
      while (p2 != p1);
      p1 = p1->nextcl;
      }
    */
    
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

    /*
    p1 = cliques[v];
    while (p1)
      {
	p2 = p1;
	do
	  {
	    if (!vertices[p2->GetVertexNr()].Used() && p2->GetVertexNr() != v)
	      {
		// p3 = new CliqueEl;
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
	p1 = p1->nextcl;
      }
    */




    for (CliqueEl * p1 = cliques[v]; p1; p1 = p1->nextcl)
      {
	CliqueEl * p2 = p1;
	do
	  {
	    if (!vertices[p2->GetVertexNr()].Used() && p2->GetVertexNr() != v)
	      {
		// CliqueEl * p3 = new CliqueEl;
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


    CliqueEl *p1, *p2, *p3, *p4;



    {
    // find dominated cliques

    // set flag for new clique
    p3 = newp;
    do
      {
	vertices[p3->GetVertexNr()].SetFlag (1);
	p3 = p3->next;
      }
    while (p3 != newp);
    

    // check all cliques of all members of new clique
    p3 = newp;
    do
      {
	p1 = cliques[p3->GetVertexNr()];
	while (p1)
	  {
	    if (!p1->eliminate)
	      {
		bool dominated = 1;
		    
		p2 = p1;
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
		    p2 = p1;
		    do
		      {
			p2->eliminate = 1;
			p2 = p2->next;
		      }
		    while (p1 != p2);
		  }
	      } 
	    p1 = p1->nextcl;
	  }
	p3 = p3->next;
      }
    while (p3 != newp);
     
    // clear flag
    p3 = newp;
    do
      {
	vertices[p3->GetVertexNr()].SetFlag (0);
	p3 = p3->next;
      }
    while (p3 != newp);
  
  
  
    // delete old cliques
    p3 = newp;
    do
      {
	// p3->next is first remaining clique:

	p3->nextcl = cliques[p3->GetVertexNr()];
	while (p3->nextcl && p3->nextcl->eliminate)
	  p3->nextcl = p3->nextcl->nextcl;


	CliqueEl hcel;
	hcel.nextcl = cliques[p3->GetVertexNr()];
	p1 = &hcel;
	  
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

  
    p1 = cliques[v];
    while (p1)
      {
	p2 = p1->GetNextClique();
	ball.Free (p1);
	p1 = p2;
      }

    cliques[v] = NULL;
    }


    // find connections
    p3 = newp;
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
		SetFlagCliques (p3->GetVertexNr());
		p4 = p3->next;
		while (p4 != newp)
		  {
		    // are p3 and p4 connected ?
		    
		    if (IsMaster (p4->GetVertexNr()))
		      {
			bool samecl = 1;
			
                        /*
			if (NumCliques (p4->GetVertexNr()) == nclp3)
			  {
                            p1 = p4;
			    while (p1)
			      {
				if (!p1->flag) 
                                  {
                                    samecl = 0;
                                    break;
                                  }
				p1 = p1->nextcl;
			      }
			  }
			else
			  samecl = 0;
                        */

                        p1 = cliques[p4->GetVertexNr()];
                        int nclp4 = 0;
                        while (p1)
                          {
                            if (!p1->flag) 
                              {
                                samecl = 0;
                                break;
                              }
                            p1 = p1->nextcl;
                            nclp4++;
                          }
                        if (nclp4 != nclp3) samecl = 0;

			if (samecl)
			  SetMaster (p3->GetVertexNr(), p4->GetVertexNr());
		      }
		    p4 = p4->next;
		  }
		ClearFlagCliques (p3->GetVertexNr());
	      }
	  }

	p3 = p3->next;
      }
    while (p3 != newp);



    // calc master degrees in new clique
    {
      p3 = newp;
      do
        {
          int v3 = p3->GetVertexNr();
          if (IsMaster(v3))
            priqueue.SetDegree (v3, CalcDegree (v3) - NumSlaves (v3));
          p3 = p3->next;
        }
      while (p3 != newp);      
    }
    

    int cnt = NumSlaves(v);
    p3 = newp;
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




    // disconnect slaves 
    cnt = 0;
    p3 = newp;
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
  }








  void MinimumDegreeOrdering :: EliminateSlaveVertex (int v)
  {
    if (cliques[v]) // NumCliques(v) != 0)
      throw Exception ("Eliminate Slave should have exactly no clique");


    // if (NumCliques(v) != 1)
    //      throw Exception ("Eliminate Slave should have exactly one clique");

    // (*testout) << "Eliminate Slave Vertex " << v << endl;
    // PrintCliques();

    CliqueEl *p1, *newp;
    
    if (cliques[v])
      {
	newp = cliques[v]->next;
	
	if (newp != cliques[v])
	  {
	    p1 = cliques[v];
	    while (1)
	      {
		if (p1->next == cliques[v])
		  {
		    p1->next = p1->next->next;
		    break;
		  }
		p1 = p1->next;
	      }
	  }
	else
	  {
	    newp = NULL;
	  }
	
	// delete cliques[v];
	ball.Free (cliques[v]);
	cliques[v] = NULL;
      }




    // PrintCliques();

    /*
    if (newp)
      {
	// eliminated node was slave
	
	p3 = newp;
	do
	  {
	    int deg = priqueue.GetDegree (p3->GetVertexNr())-1;
	    priqueue.SetDegree (p3->GetVertexNr(), deg);
	    p3 = p3->next;
	  }
	while (p3 != newp); 
      }
    */

    
    /*
    int cnt = 0;

    if (newp)
      {
	p3 = newp;
	do 
	  {
	    cnt++;
	    p3 = p3->next;
	  }
	while (p3 != newp);
      
	vertices[v].nconnected = cnt;
	vertices[v].connected = new int[cnt];
	cnt = 0;
      
	do 
	  {
	    vertices[v].connected[cnt] = p3->vnr;
	    cnt++;
	    p3 = p3->next;
	  }
	while (p3 != newp);
      }
    else
      {
	vertices[v].nconnected = 0;
	vertices[v].connected = NULL;
      }

    (*testout) << "slave " << v << ", connected: ";
    for (int i = 0; i < vertices[v].nconnected; i++)
      (*testout) << vertices[v].connected[i] << " ";
    (*testout) << endl;
*/

    /*
    int master = vertices[v].Master();
    int cnt = vertices[master].nconnected;
    int * master_connected = vertices[master].connected;
    int cnt2 = 0;
    for (int i = 0; i < cnt; i++)
      if (!vertices[master_connected[i]].Eliminated() && master_connected[i] != v)
	cnt2++;
    */

    /*
    if (cnt2 != 0)
      {
	vertices[v].nconnected = cnt2;
	vertices[v].connected = new int[cnt2];
	cnt2 = 0;
	for (int i = 0; i < cnt; i++)
	  if (!vertices[master_connected[i]].Eliminated() && master_connected[i] != v)
	    {
	      vertices[v].connected[cnt2] = master_connected[i];
	      cnt2++;
	    }
      }
    else
      {
	vertices[v].nconnected = 0;
	vertices[v].connected = NULL;
      }
    */

    vertices[v].nconnected = 0; // cnt2;
    vertices[v].connected = NULL;


    /*
    (*testout) << "slave " << v << ", connected new: ";
    for (int i = 0; i < vertices[v].nconnected; i++)
      (*testout) << vertices[v].connected[i] << " ";
    (*testout) << endl;
    */

    //  cliques.Elem(v) = savep;
    // cliques[v] = NULL;
  }



















  void MinimumDegreeOrdering :: Order()
  {
    static Timer reorder_timer("MinimumDegreeOrdering::Order");
    RegionTimer reg(reorder_timer);

    cout << "start order" << endl;

    for (int j = 0; j < n; j++)
      {
	// priqueue.SetDegree(j, CalcDegree(j));
	priqueue.SetDegree(j, 1+NumCliques(j));
	// cout << "deg = " << CalcDegree(j) << ", numcl = " << NumCliques(j) << endl;
      }

    int minj = -1;
    int lastel = -1;

    if (n > 5000)
      cout << "order " << flush;

    for (int i = 0; i < n; i++)
      {
	if (n > 5000 && i % 1000 == 999)
	  {
	    if (i % 10000 == 9999)
	      cout << "+" << flush;
	    else
	      cout << "." << flush;
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
    int cnt = 0;

    CliqueEl * p1 = cliques[v];
    while (p1)
      {
	p1 = p1->nextcl;
	cnt++;
      }
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
    CliqueEl * p1 = cliques[v];
    while (p1)
      {
	CliqueEl * p2 = p1;
	do
	  {
	    p2->flag = 1;
	    p2 = p2->next;
	  }
	while (p2 != p1);
	p1 = p1->nextcl;
      }
  }

  void MinimumDegreeOrdering :: ClearFlagCliques (int v)
  {
    CliqueEl * p1 = cliques[v];
    while (p1)
      {
	CliqueEl * p2 = p1;
	do
	  {
	    p2->flag = 0;
	    p2 = p2->next;
	  }
	while (p2 != p1);
	p1 = p1->nextcl;
      }
  }


  /*
  int MinimumDegreeOrdering :: GetNZE() const
  {
    int cnt = 0;
    for (int i = 0; i < n; i++)
      cnt += vertices[i].nconnected;
    return cnt;
  }
  */


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
