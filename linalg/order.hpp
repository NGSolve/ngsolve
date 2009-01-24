#ifndef FILE_ORDER
#define FILE_ORDER

/* *************************************************************************/
/* File:   order.hh                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   18. Jun. 97                                                    */
/* *************************************************************************/


/*
  reordering for sparse cholesky factoriztion
*/

///
class MDOPriorityQueue
{
  struct entry { 
    int degree, prev, next;
  };
  Array<entry> list;
  Array<int> first_in_class;
public:
  MDOPriorityQueue (int size, int maxdeg);
  ~MDOPriorityQueue ();
  int MinDegree () const;
  int GetDegree (int nr) const { return list[nr].degree; }
  void SetDegree (int nr, int deg);
  void Invalidate (int nr);
};


///
class MDOVertex
{
protected:
  int master;         /// master of node
  int nextslave;      /// linked list of slaves
  int numslaves;      /// number of slaves
  bool eliminated;    /// node is eliminated
  bool used;          /// temporary field (used in calcorder)
  bool flag;          


public:
  ///
  int * connected;
  ///
  int nconnected;

  ///
  MDOVertex(int ma=0)  { Init (ma); }
  ///
  ~MDOVertex()   { ; }

  /// (it is a POD !!!)
  void Init (int ma)
  {
    master = ma;
    nextslave = -1;
    numslaves = 0;
    flag = 0;
    eliminated = 0;
    used = 0;
  }

  ///
  int Master() const { return master; };
  ///
  void SetMaster(int ma) { master = ma; };
  ///
  int NextSlave () const {return nextslave; };
  ///
  void SetNextSlave( int ns ) { nextslave = ns; };
  ///
  bool Eliminated() const {return eliminated; };
  ///
  void SetEliminated(bool el) {eliminated = el; };
  ///
  bool Used() const {return used; };
  ///
  void SetUsed(bool us) {used = us; } ;
  ///
  bool Flag() const {return flag; };
  ///
  void SetFlag(bool fl) {flag = fl;};

  friend class MinimumDegreeOrdering;
};

/// 
class CliqueEl
{
public:
  /// 
  CliqueEl *next;
  CliqueEl *nextcl;
  ///
  int vnr:30;
  ///
  bool eliminate;
  ///
  bool flag;
  
  ///
  CliqueEl () {
    next = NULL;
    nextcl = NULL;
    eliminate = 0;
    flag = 0;
  }

  ///
  CliqueEl * GetNext()
  { return next; }

  ///
  CliqueEl * GetNextClique()
  { return nextcl; }

  ///
  int GetVertexNr() const
  {
    return vnr; 
  }
  
  /*  
  ///
  static ngstd::BlockAllocator ball;
  ///

  void * operator new(size_t)
  {
    cout << "call own new" << endl;
    return ball.Alloc();
  }

  ///
  void operator delete (void * p, size_t)
  {
    cout << "call own delete" << endl;
    ball.Free (p);
  }
  */
};
  


///
class MinimumDegreeOrdering
{
public:
  ///
  int n;
  ///
  Array<CliqueEl*> cliques;
  ///
  Array<int> order;
  ///
  Array<int> blocknr;
  ///
  Array<MDOVertex> vertices;
  ///
  MDOPriorityQueue priqueue;
  ///
  ngstd::BlockAllocator ball;
public:
  ///
  MinimumDegreeOrdering (int an);
  ///
  void AddEdge (int v1, int v2);
  ///
  void PrintCliques ();

  ///
  int CalcDegree (int v1);
  ///
  void EliminateMasterVertex (int v);
  void EliminateSlaveVertex (int v);
  ///
  void Order();
  /// 
  ~MinimumDegreeOrdering();

  ///
  int NumCliques (int v) const;

  /// set/clear flag for all nodes in clique
  void SetFlagNodes (int v);
  ///
  void ClearFlagNodes (int v);
  
  /// set/clear flag in all cliques of node
  void SetFlagCliques (int v);
  ///
  void ClearFlagCliques (int v);
  ///
  int Size () const { return n; }
  /// number of non-zero elements
  //  int GetNZE() const;
  //  friend class SparseCholesky;

  int NextSlave (int vnr) const
  {
    return vertices[vnr].NextSlave();
  }
  ///
  int NumSlaves (int vnr) const
  {
    return vertices[vnr].numslaves;
    /*
    int next = vertices[vnr].NextSlave();
    int cnt = 0;
    while (next != -1)
      {
	next = vertices[next].NextSlave();
	cnt++;
      }
    return cnt;
    */
  }
  ///
  bool IsMaster (int vnr) const
  {
    return vertices[vnr].Master() == vnr;
  }

  void SetMaster (int master, int slave);
};






#endif
