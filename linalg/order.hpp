#ifndef FILE_ORDER
#define FILE_ORDER

/* *************************************************************************/
/* File:   order.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   18. Jun. 97                                                    */
/* *************************************************************************/


namespace ngla
{
  using namespace ngcore;
  using namespace ngstd;  

  /*
    reordering for sparse cholesky factoriztion
  */

  ///
  class MDOPriorityQueue
  {
    struct entry { 
      int degree, prev, next;
      void DoArchive(Archive& ar) { ar & degree & prev & next; }
    };
    Array<entry> list;
    Array<int> first_in_class;
  public:
    MDOPriorityQueue (int size, int maxdeg);
    MDOPriorityQueue () {}
    ~MDOPriorityQueue ();

    void DoArchive(Archive& ar);

    int MinDegree () const;
    int GetDegree (int nr) const { return list[nr].degree; }
    void SetDegree (int nr, int deg);
    void Invalidate (int nr);

    const MemoryTracer& GetMemoryTracer() const
    {
      return mt;
    }

    void StartMemoryTracing() const
    {
      mt.Track(list, "list", first_in_class, "first_in_class");
    }
  private:
    MemoryTracer mt;
  };


  ///
  class MDOVertex
  {
  protected:
    int master;         /// master of node
    int nextminion;      /// linked list of minions
    int numminions;      /// number of minions
    int numcliques;     /// number of cliques
    bool eliminated;    /// node is eliminated
    bool used;          /// temporary field (used in calcorder)
    bool flag;          
    

  public:
    ///
    int * connected = nullptr;
    ///
    int nconnected;

    ///
    MDOVertex() = default;
    MDOVertex(int ma)  { Init (ma); }
    ///
    ~MDOVertex()   { ; }

    void DoArchive(Archive& ar);

    /// 
    void Init (int ma)
    {
      master = ma;
      nextminion = -1;
      numminions = 0;
      numcliques = 0;
      flag = 0;
      eliminated = 0;
      used = 0;
    }

    ///
    int Master() const { return master; };
    ///
    void SetMaster(int ma) { master = ma; };
    ///
    int NextMinion () const {return nextminion; };
    ///
    void SetNextMinion( int ns ) { nextminion = ns; };
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
    bool flag = false;
  public:
    CliqueEl * next = nullptr;
    CliqueEl * nextcl = nullptr;
    CliqueEl * clmaster;
    int vnr;
    bool eliminate = false;
    // in the clique master: weighted size (vertices incl. minions), and the
    // per-elimination counter of the approximate degree update
    int weight = 0, w = 0, stamp = -1;
  
    CliqueEl (int avnr) : vnr(avnr) { ; }

    CliqueEl * GetNext() { return next; }
    CliqueEl * GetNextClique() { return nextcl; }
    int GetVertexNr() const { return vnr; }
    
    void SetFlag (bool aflag) { clmaster->flag = aflag; }
    bool Flag () const { return clmaster->flag; }
    
    operator int() const { return vnr; }
    int Nr() const { return vnr; }
  };
  

  ///
  class MinimumDegreeOrdering
  {
  public:
    ///
    int n, nused;
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
    /// approximate degrees (the AMD bound) instead of exact ones: about twice as fast, similar fill
    bool approx = false;
    int curstamp = 0;
    Array<CliqueEl*> touched;   // elements met in the current elimination

    /*
      original edges: adjacency arrays (csr), built in Order() from the
      edges collected by AddEdge; entries only get removed (absorbed into
      an element, or stale: eliminated / minion vertices, skipped lazily)
    */
    Array<int> edges;           // pairs, until Order()
    Array<size_t> adjstart;
    Array<int> adjdata, adjlen, adjw;   // current length, weighted size outside the new clique
    Array<int> vmark;           // per-vertex marks of the equivalence test
    int markid = 0;
    FlatArray<int> Adj (int v) { return adjdata.Range (adjstart[v], adjstart[v]+adjlen[v]); }
  public:
    ///
    MinimumDegreeOrdering (int an);
    /// for archive
    MinimumDegreeOrdering() : ball(sizeof(CliqueEl), 1000) {}

    void DoArchive(Archive& archive);

    void AddEdge (int v1, int v2);
    ///
    void PrintCliques ();

    ///
    int CalcDegree (int v1);
    ///
    void EliminateMasterVertex (int v);
    void EliminateMinionVertex (int v);
    ///
    void Order();
    /// 
    ~MinimumDegreeOrdering();

    ///
    int NumCliques (int v) const
    {
      return vertices[v].numcliques;
    }

    ///
    void SetUnusedVertex (int v) { vertices[v].SetEliminated(true); order[v] = -1; }
      
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

    int NextMinion (int vnr) const
    {
      return vertices[vnr].NextMinion();
    }
    ///
    int NumMinions (int vnr) const 
    {
      return vertices[vnr].numminions;
    }

    ///
    bool IsMaster (int vnr) const
    {
      return vertices[vnr].Master() == vnr;
    }

    void SetMaster (int master, int minion);


    void StartMemoryTracing() const
    {
      mt.Track(cliques, "cliques",
               order, "order",
               blocknr, "blocknr",
               vertices, "vertices",
               priqueue, "priqueue",
               ball, "ball");
    }
    const MemoryTracer& GetMemoryTracer() const { return mt; }

  private:
    MemoryTracer mt;
  };


}








#endif
