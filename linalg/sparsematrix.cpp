/*********************************************************************/
/* File:   sparsematrix.cpp                                          */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


#define FILE_SPARSEMATRIX_CPP

#include <la.hpp>

#include "pardisoinverse.hpp"
#include "umfpackinverse.hpp"
#include "superluinverse.hpp"
#include "mumpsinverse.hpp"


namespace ngla
{

  int UsedBits (size_t nr)
  {
    int cnt = 0;
    size_t bit = 1;
    while (bit < nr)
      {
        cnt++;
        bit *= 2;
      }
    return cnt;
  }
  
  
  auto Split (size_t i, int bits)
  {
    auto hi = i >> bits;
    auto lo = i & ((1 << bits) -1);
    return make_pair(hi,lo);
  }  
  

  MatrixGraph :: MatrixGraph (const Array<int> & elsperrow, int awidth)
  {
    size = elsperrow.Size();
    width = awidth;
    owner = true;

    firsti.SetSize (size+1);
    nze = 0;
    for (int i = 0; i < size; i++)
      {
	firsti[i] = nze;
	nze += elsperrow[i];
      }
    firsti[size] = nze;
    
    // colnr.SetSize (nze+1);
    colnr = NumaDistributedArray<int> (nze);   // nze+1

    /*
    for (size_t i = 0; i < nze; i++)
      colnr[i] = -1;
    */
    
    FlatArray<int> hcolnr = colnr;
    ParallelForRange (nze,
		      [hcolnr] (auto myrange)
		      {
			hcolnr.Range(myrange) = -1;
		      });

    CalcBalancing ();
  }
                                                                                                                                                                                                                  
  MatrixGraph :: MatrixGraph (int as, int max_elsperrow) 
  {
    size = as;
    width = as;
    nze = as * max_elsperrow;

    // colnr.SetSize (as*max_elsperrow+1);
    colnr = NumaDistributedArray<int> (as*max_elsperrow+1);
    firsti.SetSize (as+1);
    owner = true;
    
    for (int i = 0; i < as*max_elsperrow; i++)
      colnr[i] = -1;
    colnr[as*max_elsperrow] = 0;
    
    for (int i = 0; i < as+1; i++)
      firsti[i] = i*max_elsperrow;

    CalcBalancing ();
  }
  

  
  MatrixGraph :: MatrixGraph (const MatrixGraph & agraph, bool stealgraph)
  {
    MatrixGraph & graph = const_cast<MatrixGraph&> (agraph);
    size = graph.size;
    width = graph.width;
    nze = graph.nze;
    owner = false;

    if (stealgraph)
      {
	firsti.Swap (graph.firsti);
	colnr.Swap (graph.colnr);
      }
    else
      {
	firsti.SetSize (size+1);
	// colnr.SetSize (nze);
        colnr = NumaDistributedArray<int> (nze);
        
	for (int i = 0; i < size+1; i++)
	  firsti[i] = graph.firsti[i];
	for (size_t i = 0; i < nze; i++)
	  colnr[i] = graph.colnr[i];
      }
    // inversetype = agraph.GetInverseType();
    CalcBalancing ();
  }

  MatrixGraph :: MatrixGraph (MatrixGraph && graph)
  {
    if (!graph.owner) {
      throw Exception("Matrix-Graph Move-constructor with graph that is now owner ... is that valid?");
    }
    size = move(graph.size);
    width = move(graph.width);
    nze = move(graph.nze);
    owner = true;
    firsti.Swap (graph.firsti);
    colnr.Swap (graph.colnr);
    CalcBalancing ();
  }


  /*
  template <typename FUNC>
  INLINE void MergeArrays (FlatArray<int*> ptrs,
                           FlatArray<int> sizes,
                           FUNC f)
  {
    int nactive = 0;
    for (auto s : sizes)
      if (s) nactive++;
    
    while (nactive)
      {
        int minval = numeric_limits<int>::max();
        for (int * ptr : ptrs)
          if (ptr)
            minval = min2(minval, *ptr);

        f(minval);
        for (int i : sizes.Range())
          if (ptrs[i] && (*ptrs[i] == minval))
            {
              ptrs[i]++;
              sizes[i]--;
              if (sizes[i] == 0)
                {
                    nactive--;
                    ptrs[i] = nullptr;
                  }
            }
      }
  }  
  */


  template <typename FUNC>
  INLINE void MergeArrays1 (FlatArray<int*> ptrs,
                           FlatArray<int> sizes,
                           // FlatArray<int> minvals,
                           FUNC f)
  {
    STACK_ARRAY(int, minvals, sizes.Size());
    int nactive = 0;
    for (auto i : sizes.Range())
      if (sizes[i]) 
        {
          nactive++;
          minvals[i] = *ptrs[i];
        }
      else
        minvals[i] = numeric_limits<int>::max();
    
    while (nactive)
      {
        int minval = minvals[0];
        for (int i = 1; i < sizes.Size(); i++)
          minval = min2(minval, minvals[i]);

        
        f(minval);

        for (int i : sizes.Range())
          if (minvals[i] == minval)
            {
              ptrs[i]++;
              sizes[i]--;
              if (sizes[i] == 0)
                {
                  nactive--;
                  minvals[i] = numeric_limits<int>::max();
                }
              else
                minvals[i] = *ptrs[i];
            }
      }
  }  


  template <typename FUNC>
  INLINE void MergeArrays (FlatArray<int*> ptrs,
                           FlatArray<int> sizes,
                           FUNC f)
  {
    if (ptrs.Size() <= 16)
      {
        MergeArrays1 (ptrs, sizes, f);
        return;
      }
    // static Timer tall("merge arrays"); RegionTimer reg(tall);
    
    struct valsrc { int val, src; };
    struct trange { int idx, val;  };
    
    int nactive = 0;
    int nrange = 0;
    /*
    STACK_ARRAY(valsrc, minvals, sizes.Size());
    STACK_ARRAY(trange, ranges, sizes.Size()+1);
    */
    ArrayMem<valsrc,1024> minvals(sizes.Size());
    ArrayMem<trange,1024> ranges(sizes.Size()+1);
    
    constexpr int nhash = 1024; // power of 2
    int hashes[nhash];

    for (int i = 0; i < nhash; i++)
      hashes[i] = -1;
  
    // take first value from every array ...
    for (auto i : sizes.Range())
      while (sizes[i]) 
        {
          auto val = *ptrs[i];
          sizes[i]--;
          ptrs[i]++;
          if (hashes[val&(nhash-1)] == val) continue;
          minvals[nactive].val = val;
          hashes[val&(nhash-1)] = val;
          minvals[nactive].src = i;
          nactive++;
          break;
        }
    // cout << "starting values " << endl << minvals << endl;

    // presort minvals:
    // values from ranges[i] up are all smaller or equal rangevals[i]
    // minvals[j] <= rangevals[i]  <==> j >= ranges[i]   \forall i >= 0

    while (nactive > 0)
      {
        // partial quicksort
        int lower = 0;
        if (nrange > 0) lower = ranges[nrange-1].idx;

        // cout << "start quicksort in range [" << lower << ", " << nactive-1 << "]" << endl;
        while (true)
          {          
            int firstval = minvals[lower].val;
            int otherval = firstval;

            for (int i = lower+1; i < nactive; i++)
              {
                if (minvals[i].val != firstval)
                  {
                    otherval = minvals[i].val;
                    break;
                  }
              }

            if (firstval == otherval)
              { // all values in last range are the same -> presorting commplete
                if (nrange == 0)
                  {
                    ranges[0].idx = 0;
                    ranges[0].val = firstval;
                    nrange = 1;
                  }
                break;
              }

            int midval = (firstval+otherval)/2;
            // midval is not the largest value, so we find a new separation ...
            int l = lower, r = nactive-1;
            while (l <= r)
              {
                while (minvals[l].val > midval) l++;
                while (minvals[r].val <= midval) r--;
                if (l < r)
                  Swap (minvals[l++], minvals[r--]);
              }
          
            // elements from l up are <= midval
            ranges[nrange].idx = l;
            ranges[nrange].val = midval;
            nrange++;
          
            lower = l;
          }

        nrange--;
        int last = ranges[nrange].idx; 
        f(minvals[last].val);
      
        // insert new values
        FlatArray<valsrc> tmp(nactive-last, &minvals[last]);
        nactive = last;

        for (valsrc vs : tmp)
          while (sizes[vs.src])
            {
              vs.val = *ptrs[vs.src];
              sizes[vs.src]--;
              ptrs[vs.src]++;

              // take next value if already in queue
              if (hashes[vs.val&(nhash-1)] == vs.val) continue;

              int prevpos = nactive;
              for (int i = nrange-1; i >= 0; i--)
                {
                  if (vs.val <= ranges[i].val)
                    break;
                  int pos = ranges[i].idx;
                  minvals[prevpos] = minvals[pos];
                  prevpos = pos;
                  ranges[i].idx++;
                }

              minvals[prevpos] = vs;
              hashes[vs.val&(nhash-1)] = vs.val;            
              nactive++;
              break;
            }
      }
  }




  inline void MergeSortedArrays (FlatArray<int> in1, FlatArray<int> in2,
                                 Array<int> & out)
  {
    out.SetSize(in1.Size()+in2.Size());

    int i1 = 0, i2 = 0, io = 0;
    while (i1 < in1.Size() && i2 < in2.Size())
      {
        int newel;
        if (in1[i1] == in2[i2])
          {
            newel = in1[i1++]; i2++;
          }
        else if (in1[i1] < in2[i2])
          newel = in1[i1++];
        else
          newel = in2[i2++];
        out[io++] = newel; 
      }
    
    while (i1 < in1.Size())
      out[io++] = in1[i1++];
    while (i2 < in2.Size())
      out[io++] = in2[i2++];
                      
    out.SetSize(io);
  }


  MatrixGraph :: MatrixGraph (int asize, int awidth, const Table<int> & rowelements, 
                              const Table<int> & colelements, 
                              bool symmetric)
  {
    // make sure that taskmanager is up ...
    /*
    RunWithTaskManager 
      ([&]() 
       {
    */

    static Timer timer("MatrixGraph");
    static Timer timer_dof2el("MatrixGraph - build dof2el table");
    static Timer timer_prefix("MatrixGraph - prefix");    
    RegionTimer reg (timer);

    bool includediag = (&rowelements == &colelements);
     
    int ndof = asize;
    TableCreator<int> creator(ndof);


    ParallelFor (Range(colelements.Size()), 
                 [&] (int i) { QuickSort (colelements[i]); });
    
    timer_dof2el.Start();
    for ( ; !creator.Done(); creator++)
      {    
        ParallelFor (Range(rowelements.Size()),
                     [&] (int i)
                     {
                       for (auto e : rowelements[i])
                         creator.Add(e, i);
                     },
                     TasksPerThread(10));
      }
    timer_dof2el.Stop();

    Table<int> dof2element = creator.MoveTable();

    // #define NEWDOF2EL
#ifdef NEWDOF2EL
    {
    static Timer timer_newdof2el("MatrixGraph - new build dof2el table");
    static Timer timer_newdof2elb("MatrixGraph - new build dof2el table b");
    static Timer timer_newdof2el_1("MatrixGraph - new build dof2el table 1");
    static Timer timer_newdof2el_2("MatrixGraph - new build dof2el table 2");

            
    // atomic-free dof2element via COO format
    timer_newdof2elb.Start();
    int rowbits = UsedBits (rowelements.Size());
    int rowbits_hi = rowbits / 2;
    int rowbits_lo = rowbits - rowbits_hi;
    int rows_hi = 1 << rowbits_hi;
    int rows_lo = 1 << rowbits_lo;
    int nrows_hi = (rowelements.Size()+rows_lo-1) / rows_lo;
    
    int dofbits = UsedBits (ndof);
    int dofbits_hi = dofbits/2;
    int dofbits_lo = dofbits - dofbits_hi;
    int dofs_hi = 1 << dofbits_hi;
    int dofs_lo = 1 << dofbits_lo;    
    int ndofs_hi = (ndof+dofs_lo-1)/dofs_lo;

    { // build a global table
      Array<int> cnt_coo(nrows_hi*ndofs_hi);
      // ParallelFor (nrows_hi, [&] (int row_hi)
      ParallelJob ([&] (TaskInfo &ti)
                   {
		     int row_hi = ti.task_nr;
                     auto mycnt = cnt_coo.Range(row_hi*ndofs_hi, (row_hi+1)*ndofs_hi);
                     mycnt = 0;
                     for (int row_lo = 0; row_lo < rows_lo; row_lo++)
                       {
                         int row = (row_hi << rowbits_lo)+row_lo;
                         if (row < rowelements.Size())
                           {
                             for (auto d : rowelements[row])
                               {
                                 int dof_hi, dof_lo;
                                 tie(dof_hi, dof_lo) = Split (d, dofbits_lo);
                                 mycnt[dof_hi]++;
                               }
                           }
                       }
		     // }, TasksPerThread(4));
		   }, nrows_hi);
      Table<int> dof2element_coo(cnt_coo);

      /*
      ParallelJob ([&] (TaskInfo &ti)
                   {
		     int row_hi = ti.task_nr;
		     auto mytable = dof2element_coo.Range(row_hi*ndofs_hi, (row_hi+1)*ndofs_hi);
		     mytable.AsArray() = 0;
		     // for (auto i = row_hi*ndofs_hi; i < (row_hi+1)*ndofs_hi; i++)
		     // dof2element_coo[i] = 0;
		   }, nrows_hi);
      */
      /*
      ParallelForRange (dof2element_coo.Size(), [&] (IntRange r)
			{
			  dof2element_coo.Range(r).AsArray() = 0;
			}, 1+0.3*TaskManager::GetNumThreads());
      */
      // ParallelFor (nrows_hi, [&] (int row_hi)
      ParallelJob ([&] (TaskInfo &ti)
                   {
		     int row_hi = ti.task_nr;
                     auto mycnt = cnt_coo.Range(row_hi*ndofs_hi, (row_hi+1)*ndofs_hi);
                     auto mytable = dof2element_coo.Range(row_hi*ndofs_hi, (row_hi+1)*ndofs_hi);
                     mycnt = 0;
                     for (int row_lo = 0; row_lo < rows_lo; row_lo++)
                       {
                         int row = (row_hi << rowbits_lo)+row_lo;
                         if (row < rowelements.Size())
                           {
                             for (auto d : rowelements[row])
                               {
                                 int dof_hi, dof_lo;
                                 tie(dof_hi, dof_lo) = Split (d, dofbits_lo);
                                 mytable[dof_hi][mycnt[dof_hi]++] = row_lo*dofs_lo+dof_lo;
                               }
                           }
                       }
		   }, nrows_hi);
      // }, TasksPerThread(4));


    Array<int> newcnt(ndof);
    // ParallelFor(ndofs_hi, [&] (int dof_hi)
    ParallelJob ([&] (TaskInfo &ti)
		 {
		   int dof_hi = ti.task_nr;
                  auto first = dof_hi << dofbits_lo;
                  if (first >= ndof) return;
                  auto next = min ( (dof_hi+1) << dofbits_lo, ndof);
                  newcnt.Range(first, next) = 0;
                  
                  // for (int row_hi = 0; row_hi < rows_hi; row_hi++)
                  for (int row_hi : Range(nrows_hi))
                    for (auto p : dof2element_coo[dof_hi+row_hi*ndofs_hi])
                      {
                        int dof_lo, row_lo;
                        tie(row_lo, dof_lo) = Split (p, dofbits_lo);
                        int dof = (dof_hi << dofbits_lo) + dof_lo;
                        newcnt[dof]++;
                      }
		 }, ndofs_hi);
		  // }, TasksPerThread(4));
    
    Table<int> newdof2element(newcnt);
    // ParallelFor (ndofs_hi, [&] (int dof_hi)
    ParallelJob ([&] (TaskInfo &ti)
                 {
		   int dof_hi = ti.task_nr;
                   auto first = dof_hi << dofbits_lo;
                   if (first >= ndof) return;
                   auto next = min ( (dof_hi+1) << dofbits_lo, ndof);
                   newcnt.Range(first, next) = 0;
                   
                   // for (int row_hi = 0; row_hi < rows_hi; row_hi++)
                   for (int row_hi : Range(nrows_hi))
                     for (auto p : dof2element_coo[dof_hi+row_hi*ndofs_hi])                     
                       {
                         int dof_lo, row_lo;
                         tie(row_lo, dof_lo) = Split (p, dofbits_lo);
                         int dof = (dof_hi << dofbits_lo) + dof_lo;
                         newdof2element[dof][newcnt[dof]++] = row_hi*rows_lo+row_lo;
                       }
		 }, ndofs_hi);
    // }, TasksPerThread(4));
    timer_newdof2elb.Stop();

    if (dof2element.Size() < 100)
      {
        cout << "dof2el = " << dof2element << endl;
        cout << "newdof2el = " << newdof2element << endl;
      }
    }
    
    timer_newdof2el.Start();
    
    Array<Table<int>> entries(nrows_hi);
    ParallelFor (nrows_hi, [&] (int row_hi)
                 {
                   Array<int> cnt_entries(ndofs_hi);
                   {
                   RegionTracer reg(TaskManager::GetThreadId(), timer_newdof2el_1);
                   cnt_entries = 0;
                   
                   for (int row_lo = 0; row_lo < rows_lo; row_lo++)
                     {
                       int row = (row_hi << rowbits_lo)+row_lo;
                       if (row < rowelements.Size())
                         {
                           for (auto d : rowelements[row])
                             {
                               int dof_hi, dof_lo;
                               tie(dof_hi, dof_lo) = Split (d, dofbits_lo);
                               cnt_entries[dof_hi]++;
                             }
                         }
                     }
                   }
                   Table<int> loctable(cnt_entries);
                   cnt_entries = 0;

                   {
                   RegionTracer reg(TaskManager::GetThreadId(), timer_newdof2el_2);
                   for (int row_lo = 0; row_lo < rows_lo; row_lo++)
                     {
                       int row = (row_hi << rowbits_lo)+row_lo;
                       if (row < rowelements.Size())
                         {
                           for (auto d : rowelements[row])
                             {
                               int dof_hi, dof_lo;
                               tie(dof_hi, dof_lo) = Split (d, dofbits_lo);
                               loctable[dof_hi][cnt_entries[dof_hi]++] = row_lo*dofs_lo+dof_lo;
                             }
                         }
                     }
                   }
                   entries[row_hi] = move(loctable);
                 },
                 TasksPerThread(4));
    
    
    Array<int> newcnt(ndof);
    ParallelFor(ndofs_hi, [&] (int dof_hi)
                {
                  auto first = dof_hi << dofbits_lo;
                  if (first >= ndof) return;
                  auto next = min ( (dof_hi+1) << dofbits_lo, ndof);
                  newcnt.Range(first, next) = 0;
                  
                  // for (int row_hi = 0; row_hi < rows_hi; row_hi++)
                  for (int row_hi : Range(entries))
                    for (auto p : entries[row_hi][dof_hi])
                      {
                        int dof_lo, row_lo;
                        tie(row_lo, dof_lo) = Split (p, dofbits_lo);
                        int dof = (dof_hi << dofbits_lo) + dof_lo;
                        newcnt[dof]++;
                      }
                }, TasksPerThread(4));
    
    Table<int> newdof2element(newcnt);
    ParallelFor (ndofs_hi, [&] (int dof_hi)
                 {
                   auto first = dof_hi << dofbits_lo;
                   if (first >= ndof) return;
                   auto next = min ( (dof_hi+1) << dofbits_lo, ndof);
                   newcnt.Range(first, next) = 0;
                   
                   // for (int row_hi = 0; row_hi < rows_hi; row_hi++)
                   for (int row_hi : Range(entries))                     
                     for (auto p : entries[row_hi][dof_hi])                         
                       {
                         int dof_lo, row_lo;
                         tie(row_lo, dof_lo) = Split (p, dofbits_lo);
                         int dof = (dof_hi << dofbits_lo) + dof_lo;
                         newdof2element[dof][newcnt[dof]++] = row_hi*rows_lo+row_lo;
                       }
                 }, TasksPerThread(4));
    timer_newdof2el.Stop();    

    if (dof2element.Size() < 100)
      {
        cout << "dof2el = " << dof2element << endl;
        cout << "newdof2el = " << newdof2element << endl;
      }

    // delete table-memory on same thread as it was created
    ParallelFor (nrows_hi, [&] (int row_hi)
                 { 
                   entries[row_hi] = Table<int>();
                 });
    

    // just for testing: with numa local array:
    ParallelFor (nrows_hi, [&] (int row_hi)
                 {
                   Array<int> cnt_entries(ndofs_hi);
                   cnt_entries = 0;
                   
                   for (int row_lo = 0; row_lo < rows_lo; row_lo++)
                     {
                       int row = (row_hi << rowbits_lo)+row_lo;
                       if (row < rowelements.Size())
                         {
                           for (auto d : rowelements[row])
                             {
                               int dof_hi, dof_lo;
                               tie(dof_hi, dof_lo) = Split (d, dofbits_lo);
                               cnt_entries[dof_hi]++;
                             }
                         }
                     }
                 });

    ParallelFor (nrows_hi, [&] (int row_hi)
                 {
                   NumaLocalArray<int> cnt_entries(ndofs_hi);
                   cnt_entries = 0;
                   
                   for (int row_lo = 0; row_lo < rows_lo; row_lo++)
                     {
                       int row = (row_hi << rowbits_lo)+row_lo;
                       if (row < rowelements.Size())
                         {
                           for (auto d : rowelements[row])
                             {
                               int dof_hi, dof_lo;
                               tie(dof_hi, dof_lo) = Split (d, dofbits_lo);
                               cnt_entries[dof_hi]++;
                             }
                         }
                     }
                 });
                 
    }

    
#endif
    
    Array<int> cnt(ndof);
    // cnt = 0;
    ParallelJob ([&] (TaskInfo ti)
                 {
                   auto r = Range(ndof).Split(ti.task_nr, ti.ntasks);
                   cnt[r] = 0;
                 });

    /*
    class ProfileData
    {
    public:
      double tstart, tend;
      int size;
    };
    Array<ProfileData> prof(ndof);
    */
    


    for (int loop = 1; loop <= 2; loop++)
      {
        if (!symmetric)
          {
            ParallelForRange 
              (Range(ndof), [&](IntRange myr) 
               {
                 ArrayMem<int, 50> sizes;
                 ArrayMem<int*, 50> ptrs;

                 for (int i : myr)
                   {
                     sizes.SetSize(dof2element[i].Size());
                     ptrs.SetSize(dof2element[i].Size());

                     for (int j : dof2element[i].Range())
                       {
                         sizes[j] = colelements[dof2element[i][j]].Size();
                         ptrs[j] = colelements[dof2element[i][j]].Addr(0);
                       }
                     
                     if (loop == 1)
                       {
                         int cnti = 0;
                         MergeArrays(ptrs, sizes, [&cnti] (int col) { cnti++; } );
                         cnt[i] = cnti;
                       }
                     else
                       {
                         auto ptr = colnr.Data()+firsti[i];
                         MergeArrays(ptrs, sizes, [&ptr] (int col) 
                                     {
                                       *ptr = col;
                                       ptr++;
                                     } );
                       }
                   }
               },
               TasksPerThread(20));
            // 20 * task_manager->GetNumThreads());
          }
        else
          {
            ParallelForRange 
              (Range(ndof),[&](IntRange myr)
               {
                 Array<int> rowdofs;
                 Array<int> rowdofs1;
                 
                 // for (int i : sl)
                 for (int i : myr)
                   {
                     rowdofs.SetSize0();
                     if (includediag) rowdofs += i;
                     
                     for (auto elnr : dof2element[i])
                       {
                         rowdofs.Swap (rowdofs1);
                         auto row = colelements[elnr];
                         
                         rowdofs.SetSize(rowdofs1.Size()+row.Size());
                         
                         int i1 = 0, i2 = 0, i3 = 0;
                         while (i1 < rowdofs1.Size() && i2 < row.Size() && row[i2] <= i)
                           {
                             int newel;
                             if (rowdofs1[i1] == row[i2])
                               {
                                 newel = rowdofs1[i1++]; i2++;
                               }
                             else if (rowdofs1[i1] < row[i2])
                               newel = rowdofs1[i1++];
                             else
                               newel = row[i2++];
                             rowdofs[i3++] = newel; 
                           }
                         
                         while (i1 < rowdofs1.Size())
                           rowdofs[i3++] = rowdofs1[i1++];
                         while (i2 < row.Size() && row[i2] <= i)
                           rowdofs[i3++] = row[i2++];
                         
                         rowdofs.SetSize(i3);
                       }
                  
                  
                     if (loop == 1)
                       cnt[i] = rowdofs.Size();
                     else
                       colnr.Range(firsti[i], firsti[i+1]) = rowdofs;
                   }
               }, TasksPerThread(5));
            
          }
        
        
        if (loop == 1)
          {
            size = ndof;
            width = awidth;
            owner = true;
            
            firsti.SetSize (size+1);
            /*
            nze = 0;
            for (int i = 0; i < size; i++)
              {
                firsti[i] = nze;
                nze += cnt[i];
              }
            firsti[size] = nze;
            */
            
            timer_prefix.Start();
            Array<size_t> partial_sums(TaskManager::GetNumThreads()+1);
            partial_sums[0] = 0;
            ParallelJob
              ([&] (TaskInfo ti)
               {
                 IntRange r = IntRange(size).Split(ti.task_nr, ti.ntasks);
                 size_t mysum = 0;
                 for (size_t i : r)
                   mysum += cnt[i];
                 partial_sums[ti.task_nr+1] = mysum;
               });

            for (size_t i = 1; i < partial_sums.Size(); i++)
              partial_sums[i] += partial_sums[i-1];

            ParallelJob
              ([&] (TaskInfo ti)
               {
                 IntRange r = IntRange(size).Split(ti.task_nr, ti.ntasks);
                 size_t mysum = partial_sums[ti.task_nr];
                 for (size_t i : r)
                   {
                     firsti[i] = mysum;
                     mysum += cnt[i];
                   }
               });
            nze = partial_sums[partial_sums.Size()-1];
            firsti[size] = nze;
            timer_prefix.Stop();
            
            colnr = NumaDistributedArray<int> (nze);  // nze+1

	    CalcBalancing ();

            // first touch memory (numa!)
            ParallelFor (balance, [&](int row) 
                         {
                           colnr.Range(firsti[row], firsti[row+1]) = 0;
                         });
          }
        else
          {
            ;
          }
      }
    /*
    ofstream out("creategraph.out");
    double sumtime = 0;
    for (int i = 0; i < prof.Size(); i++)
      {
        out << i << " " << prof[i].size << " ";
        double t = (prof[i].tend-prof[i].tstart) * 1e6;
        if (t > 100) out << "           ";
        out << t << endl;
        sumtime += t;
      }
    cout << "sum time = " << sumtime << endl;
    cout << "sum time per thread = " << sumtime/48 << endl;
    */
    // });
  }

  
  /*
  MatrixGraph :: MatrixGraph (const Table<int> & dof2dof, 
			      bool symmetric)
  {
    static Timer timer ("MatrixGraph");
    RegionTimer reg (timer);

    int ndof = dof2dof.Size();

    Array<int> cnt(ndof);
    cnt = 0;
    for (int i = 0; i < dof2dof.Size(); i++)
      {
        FlatArray<int> dof = dof2dof[i];
        for (int j = 0; j < dof.Size(); j++)
	  if (!symmetric || dof[j]>=i)
	    cnt[dof[j]]++;
      }

    size = ndof;
    owner = true;
    
    firsti.Alloc (size+1);
    firsti.SetName ("matrix graph, table 1");

    nze = 0;
    for (int i = 0; i < size; i++)
      {
	firsti[i] = nze;
	nze += cnt[i];
      }
    firsti[size] = nze;
    
    colnr.Alloc (nze+1);
    colnr.SetName ("matrix graph");
    
    Array<int> mark(ndof);

//     cnt = 0;

    mark = -1;

    if (!symmetric)

      for (int i = 0; i < ndof; i++)
        {
          int cnti = firsti[i];
          mark[i] = i;
          colnr[cnti++] = i;
          
              for (int k = 0; k < dof2dof[i].Size(); k++)
                {
                  int d2 = dof2dof[i][k];
		  if (mark[d2] != i)
		    {
		      mark[d2] = i;
		      colnr[cnti++] = d2;
		    }
                }
        }
    
    else
      for (int i = 0; i < ndof; i++)
        {
          int cnti = firsti[i];
          mark[i] = i;
          colnr[cnti++] = i;
          
	  for (int k = 0; k < dof2dof[i].Size(); k++)
	    {
	      int d2 = dof2dof[i][k];
	      if(d2<i)
		if (mark[d2] != i)
		  {
		    mark[d2] = i;
		    colnr[cnti++] = d2;
		  }
	    }
        }
    
    for (int i = 0; i < ndof; i++)
      QuickSort (GetRowIndices(i));
    
    colnr[nze] = 0;
  }
  */

  
  MatrixGraph :: ~MatrixGraph ()
  {
    ;
  }
  
  void MatrixGraph :: Compress()
  {
    cout << "compress not implemented" << endl; 
  }
  

  /// returns position of Element (i, j), exception for unused
  size_t MatrixGraph :: GetPosition (int i, int j) const
  {
    /*
      for (int k = firsti[i]; k < firsti[i+1]; k++)
      if (colnr[k] == j) return k;
    */
    
    size_t first = firsti[i];
    size_t last = firsti[i+1];
    while (last > first + 5)
      {
	size_t mid = (first+last) / 2;
	if (colnr[mid] > j)
	  last = mid;
	else
	  {
	    if (colnr[mid] == j) return mid;
	    first = mid+1;
	  }
      }
    for (size_t k = first; k < last; k++)
      if (colnr[k] == j) return k;
    
    stringstream err;
    err << "illegal position: " << i << ", " << j << endl;
    throw Exception (err.str());
  }
  
  
  /// returns position of Element (i, j), -1 for unused
  size_t MatrixGraph :: GetPositionTest (int i, int j) const
  {
    /*
      for (int k = firsti[i]; k < firsti[i+1]; k++)
      if (colnr[k] == j) return k;
    */
    
    size_t first = firsti[i];
    size_t last = firsti[i+1];
    while (last > first + 5)
      {
	size_t mid = (first+last) / 2;
	if (colnr[mid] > j)
	  last = mid;
	else
	  {
	    if (colnr[mid] == j) return mid;
	    first = mid+1;
	  }
      }
    for (size_t k = first; k < last; k++)
      if (colnr[k] == j) return k;


    return numeric_limits<size_t>::max();
  }
  
  size_t MatrixGraph :: CreatePosition (int i, int j)
  {
    size_t first = firsti[i]; 
    size_t last = firsti[i+1];
    /*
      (*testout) << "row = " << i << ", col = " << j << endl;
      (*testout) << "first = " << first << ", last = " << last << endl;
      for (int k = first; k < last; k++)
      (*testout) << colnr[k] << " ";
      (*testout) << endl;
    */
    // while (last > first + 5)
    while (last > first + 2)
      {
	size_t mid = (first+last) / 2;
	// (*testout) << "first = " << first << ", last = " << last << ", mid = " << mid << ", colnr[mid] = " << colnr[mid] << endl;

        if (colnr[mid] == j) return mid;

	if (colnr[mid] > j || colnr[mid] == -1)
	  last = mid+1;
	else
          first = mid+1;
      }


    for (size_t k = first; k < last; k++)
      {
	if (colnr[k] == -1)
	  {
	    colnr[k] = j;
	    return k;
	  }
	
	if (colnr[k] == j) return k;
	
	if (colnr[k] > j)
	  {
	    if (colnr[firsti[i+1]-1] != -1)
	      throw Exception ("sparse matrix row full 1 !");
	    
	    for (size_t l = firsti[i+1]-1; l > k; l--)
	      colnr[l] = colnr[l-1];

	    colnr[k] = j;
	    return k;
	  }
      }


    /*
      for (int k = firsti[i]; k < firsti[i+1]; k++)
      {
      if (colnr[k] == -1)
      {
      colnr[k] = j;
      return k;
      }
	
      if (colnr[k] == j) return k;
	
      if (colnr[k] > j)
      {
      if (colnr[firsti[i+1]-1] != -1)
      throw Exception ("sparse matrix row full 1 !");
	    
      for (int l = firsti[i+1]-1; l > k; l--)
      colnr[l] = colnr[l-1];

      colnr[k] = j;
      return k;
      }
      }
    */
    throw Exception ("sparse matrix row full 2 !");
  }
  


  void MatrixGraph :: 
  GetPositionsSorted (int row, int n, int * pos) const
  {
    if (n == 1)
      {
	pos[0] = GetPosition (row, pos[0]);
	return;
      }
    
    int i = 0;
    int posi = pos[i];
    size_t endk = firsti[row+1];
    for (size_t k = firsti[row]; k < endk; k++)
      {
	if (colnr[k] == posi)
	  {
	    pos[i] = k;
	    i++;
	    if (i == n) return;
	    posi = pos[i];
	  }
      }

    throw Exception ("GetPositionSorted: not matching");
  }

  
  
  template <typename Tarray>
  int BinSearch(const Tarray & v, size_t i) {
    int n = v.Size();
    if (n == 0) return 0;
    
    int first = 0;
    int last = n-1;
    if(v[0]>i) return 0;
    if(v[n-1] <= i) return n;
    while(last-first>1) {
      int m = (first+last)/2;
      if(v[m]<i)
            first = m;
      else
	last = m;
    }
    return first;
  }
  

  void MatrixGraph :: CalcBalancing ()
  {
    static Timer timer ("MatrixGraph - CalcBalancing");
    RegionTimer reg (timer);

    balance.Calc (size, [&] (int row) { return 1 + GetRowIndices(row).Size(); });
  }
  
  void MatrixGraph :: FindSameNZE()
  {
    return;

    same_nze.SetSize (size);

    same_nze[0] = 0;
    for (int i = 1; i < size; i++)
      if (GetRowIndices(i) == GetRowIndices(i-1))
	same_nze[i] = same_nze[i-1];
      else
	same_nze[i] = i;

    (*testout) << "same_nze = " << endl << same_nze << endl;

    int sum = 0;
    for (int i = 0; i < size; i++)
      if (same_nze[i] != i) sum++;
    cout << "same_nze = " << sum << "out of " << size << endl;
  }
  
  ostream & MatrixGraph :: Print (ostream & ost) const
  {
    for (int i = 0; i < size; i++)
      {
	ost << "Row " << i << ":";
	
	for (size_t j = firsti[i]; j < firsti[i+1]; j++)
	  ost << " " << colnr[j];
	ost << "\n";
      }
    return ost;
  }


  Array<MemoryUsage> MatrixGraph :: GetMemoryUsage () const
  {
    return { { "MatrixGraph", (nze+size)*sizeof(int), 1 } };
  }









  BaseSparseMatrix :: ~BaseSparseMatrix ()
  { 
    ;
  }

  INVERSETYPE BaseSparseMatrix ::
  SetInverseType (string ainversetype) const
  {
    INVERSETYPE old_invtype = inversetype;

    if      (ainversetype == "pardiso" )      SetInverseType ( PARDISO );
    else if (ainversetype == "pardisospd")    SetInverseType ( PARDISOSPD );
    else if (ainversetype == "superlu")       SetInverseType ( SUPERLU );
    else if (ainversetype == "superlu_dist")  SetInverseType ( SUPERLU_DIST );
    else if (ainversetype == "mumps")         SetInverseType ( MUMPS );
    else if (ainversetype == "masterinverse") SetInverseType ( MASTERINVERSE );
    else if (ainversetype == "sparsecholesky") SetInverseType ( SPARSECHOLESKY );
    else if (ainversetype == "umfpack")       SetInverseType ( UMFPACK );
    else
      {
        throw Exception (ToString("undefined inverse ")+ainversetype+
                         "\nallowed is: 'sparsecholesky', 'pardiso', 'pardisospd', 'mumps', 'masterinverse', 'umfpack'");
      }
    return old_invtype;
  }
}


#include "sparsematrix_impl.hpp"


namespace ngla
{


  // #ifdef OLD  leave it for a while
  shared_ptr<SparseMatrixTM<double>>
  TransposeMatrix (const SparseMatrixTM<double> & mat)
  {
    return dynamic_pointer_cast<SparseMatrixTM<double>> (mat.CreateTranspose());
    /*
    static Timer t1("TransposeMatrix 1");
    static Timer t2("TransposeMatrix 2");
    t1.Start();
    Array<int> cnt(mat.Width());
    cnt = 0;
    ParallelFor (mat.Height(), [&] (int i)
                 {
                   for (int c : mat.GetRowIndices(i))
                     AsAtomic (cnt[c]) ++;
                 });
    t1.Stop();
    t2.Start();

    shared_ptr<SparseMatrixTM<double>> trans;
    if (dynamic_cast<const SparseMatrix<double,double,double>*>(&mat))
      trans = make_shared<SparseMatrix<double>>(cnt, mat.Height());
    else if (dynamic_cast<const SparseMatrix<double,Complex,Complex>*>(&mat))
      trans = make_shared<SparseMatrix<double,Complex,Complex>>(cnt, mat.Height());
    else
      throw Exception(string("SparseMatrix-Transpose, not available for type")+typeid(mat).name()); 

    cnt = 0;
    ParallelFor (mat.Height(), [&] (int i)
                 {
                   for (int ci : Range(mat.GetRowIndices(i)))
                     {
                       int c = mat.GetRowIndices(i)[ci];
                       int pos = AsAtomic(cnt[c])++;
                       trans -> GetRowIndices(c)[pos] = i;
                       trans -> GetRowValues(c)[pos] = mat.GetRowValues(i)[ci];
                     }
                 });

    ParallelFor (trans->Height(), [&] (int r)
                 {
                   auto rowvals = trans->GetRowValues(r);
                   BubbleSort (trans->GetRowIndices(r),
                               FlatArray<double> (rowvals.Size(), rowvals.Data()));
                 });

    t2.Stop();
    return trans;
    */
  }
  // #endif
  
  shared_ptr<SparseMatrix<double,double>>
  MakeFullMatrix (const SparseMatrix<double, double> & mat)
  {
    Array<int> cnt(mat.Width());
    cnt = 0;
    for (int i = 0; i < mat.Height(); i++)
      {
        cnt[i] += mat.GetRowIndices(i).Size();
        for (int c : mat.GetRowIndices(i))
          if (c < i) cnt[c]++;
      }
    
    auto full = make_shared<SparseMatrix<double>>(cnt);
    cnt = 0;

    ParallelFor (mat.Height(), [&] (int i)
                 {
                   for (int ci : Range(mat.GetRowIndices(i)))
                     {
                       full -> GetRowIndices(i)[cnt[i]] = mat.GetRowIndices(i)[ci];
                       full -> GetRowValues(i)[cnt[i]] = mat.GetRowValues(i)[ci];
                       cnt[i] ++;
                     }
                 });
    
    for (int i = 0; i < mat.Height(); i++)
      for (int ci : Range(mat.GetRowIndices(i)))
        {
          int c = mat.GetRowIndices(i)[ci];
          if (c == i) continue;
          full -> GetRowIndices(c)[cnt[c]] = i;
          full -> GetRowValues(c)[cnt[c]] = mat.GetRowValues(i)[ci];
          cnt[c] ++;
        }
    return full;    
  }

  shared_ptr<SparseMatrixSymmetric<double,double>>
  GetSymmetricMatrix (SparseMatrixTM<double> & mat)
  {
    Array<int> cnt(mat.Width());
    cnt = 0;
    for (int i = 0; i < mat.Height(); i++)
      for (int c : mat.GetRowIndices(i))
        if (c <= i)
          cnt[i]++;

    auto symm = make_shared<SparseMatrixSymmetric<double>>(cnt);

    for (int i = 0; i < mat.Height(); i++)
      for (int ci : Range(symm->GetRowIndices(i)))
        {
            symm -> GetRowIndices(i)[ci] = mat.GetRowIndices(i)[ci];
            symm -> GetRowValues(i)[ci] = mat.GetRowValues(i)[ci];
        }

    return symm;    
  }


  template <typename TM_Res, typename TM1, typename TM2>
  shared_ptr<SparseMatrixTM<TM_Res>>
  MatMult (const SparseMatrixTM<TM1> & mata, const SparseMatrixTM<TM2> & matb)
  {
    static Timer t ("sparse matrix multiplication");
    static Timer t1a ("sparse matrix multiplication - setup a");
    static Timer t1b ("sparse matrix multiplication - setup b");
    static Timer t1b1 ("sparse matrix multiplication - setup b1");
    static Timer t2 ("sparse matrix multiplication - mult"); 
    RegionTimer reg(t);

    t1a.Start();

    // find graph of product
    Array<int> cnt(mata.Height());
    cnt = 0;

    ParallelForRange
      (mata.Height(), [&] (IntRange r)
       {
         Array<int*> ptrs;
         Array<int> sizes;
         for (int i : r)
           {
             auto mata_ci = mata.GetRowIndices(i);
             ptrs.SetSize(mata_ci.Size());
             sizes.SetSize(mata_ci.Size());
             for (int j : Range(mata_ci))
               {
                 ptrs[j] = matb.GetRowIndices(mata_ci[j]).Addr(0);
                 sizes[j] = matb.GetRowIndices(mata_ci[j]).Size();
               }
             int cnti = 0;
             MergeArrays(ptrs, sizes, [&cnti] (int col) { cnti++; } );
             cnt[i] = cnti;
           }
       },
       TasksPerThread(10));

    t1a.Stop();
    t1b.Start();
    t1b1.Start();
    // auto prod = make_shared<SparseMatrix<TM_Res>>(cnt, matb.Width());
    shared_ptr<SparseMatrixTM<TM_Res>> prod;
    
    if constexpr (is_same<TM_Res,double>()) {
        if (dynamic_cast<const SparseMatrix<double,double,double>*>(&mata) &&
            dynamic_cast<const SparseMatrix<double,double,double>*>(&matb))
          prod = make_shared<SparseMatrix<TM_Res,double,double>>(cnt, matb.Width());
        else if (dynamic_cast<const SparseMatrix<double,Complex,Complex>*>(&mata) ||
                 dynamic_cast<const SparseMatrix<double,Complex,Complex>*>(&matb))
          prod = make_shared<SparseMatrix<TM_Res,Complex,Complex>>(cnt, matb.Width());
        else
          prod = make_shared<SparseMatrix<TM_Res>>(cnt, matb.Width());  // as it was, no complex supported
      }
    else
      prod = make_shared<SparseMatrix<TM_Res>>(cnt, matb.Width());  // as it was, no complex supported      
                      
        
        
    prod->AsVector() = 0.0;
    t1b1.Stop();
    // fill col-indices
    ParallelForRange
      (mata.Height(), [&] (IntRange r)
       {
         Array<int*> ptrs;
         Array<int> sizes;
         for (int i : r)
           {
             auto mata_ci = mata.GetRowIndices(i);
             ptrs.SetSize(mata_ci.Size());
             sizes.SetSize(mata_ci.Size());
             for (int j : Range(mata_ci))
               {
                 ptrs[j] = matb.GetRowIndices(mata_ci[j]).Addr(0);
                 sizes[j] = matb.GetRowIndices(mata_ci[j]).Size();
               }
             int * ptr = prod->GetRowIndices(i).Addr(0);
             MergeArrays(ptrs, sizes, [&ptr] (int col)
                         {
                           *ptr = col;
                           ptr++;
                         } );
           }
       },
       TasksPerThread(10));


    t1b.Stop();
    t2.Start();
    
    ParallelForRange
      (mata.Height(), [&] (IntRange r)
       {
         struct thash { int idx; int pos; };

         size_t maxci = 0;
         for (auto i : r)
           maxci = max2(maxci, size_t (prod->GetRowIndices(i).Size()));

         size_t nhash = 2048;
         while (nhash < 2*maxci) nhash *= 2;
         ArrayMem<thash,2048> hash(nhash);
         size_t nhashm1 = nhash-1;

         /*
         // constexpr unsigned int nhash = 1024; // power of 2, fits well onto stack and cash ..
         // constexpr unsigned int nhash = 16384;
         size_t nhash = 16384;
         thash hash[nhash];
         */
         
         for (auto i : r)
           {
             auto mata_ci = mata.GetRowIndices(i);
             auto matc_ci = prod->GetRowIndices(i);
             auto matc_vals = prod->GetRowValues(i);
             
             for (int k = 0; k < matc_ci.Size(); k++)
               {
                 size_t hashval = size_t(matc_ci[k]) & nhashm1; // % nhash;
                 hash[hashval].pos = k;
                 hash[hashval].idx = matc_ci[k];
               }
             
             for (int j : Range(mata_ci))
               {
                 auto vala = mata.GetRowValues(i)[j];
                 int rowb = mata.GetRowIndices(i)[j];
                 
                 auto matb_ci = matb.GetRowIndices(rowb);
                 auto matb_vals = matb.GetRowValues(rowb);
                 for (int k = 0; k < matb_ci.Size(); k++)
                   {
                     auto colb = matb_ci[k];
                     unsigned hashval = unsigned(colb) & nhashm1; // % nhash;
                     if (hash[hashval].idx == colb)
                       { // lucky fast branch
                        matc_vals[hash[hashval].pos] += vala * matb_vals[k]; 
                       }
                     else
                      { // do the binary search
                        (*prod)(i,colb) += vala * matb_vals[k];
                      }
                   }
               }
           }
       },
       TasksPerThread(10));

    t2.Stop();
    return prod;
  }

  shared_ptr<SparseMatrixTM<double>> MatMult (const SparseMatrixTM<double> & mata,
                                              const SparseMatrixTM<double> & matb)
  {
    return MatMult<double, double, double>(mata, matb);
  }

  template <class TM, class TV>
  shared_ptr<BaseSparseMatrix>
  SparseMatrixSymmetric<TM,TV> :: Restrict (const SparseMatrixTM<double> & prol,
					    shared_ptr<BaseSparseMatrix> acmat ) const
  {
    static Timer t ("sparsematrix - restrict");
    static Timer tbuild ("sparsematrix - restrict, build matrix");
    static Timer tcomp ("sparsematrix - restrict, compute matrix");
    RegionTimer reg(t);
    
    int n = this->Height();

    auto cmat = dynamic_pointer_cast<SparseMatrixSymmetric<TM,TV>>(acmat);
 
    // if no coarse matrix, build up matrix-graph!
    if ( !cmat )
      {
        RegionTimer reg(tbuild);

	Array<int> marks(n);
	Array<INT<2> > e2v;
	for (int i = 0; i < n; i++)
	  for (int j = 0; j < this->GetRowIndices(i).Size(); j++)
	    {
	      int col = this->GetRowIndices(i)[j];
	      FlatArray<int> prol_rowind = prol.GetRowIndices(i);
	      FlatArray<int> prol_colind = prol.GetRowIndices(col);

	      for (int k = 0; k < prol_rowind.Size(); k++)
		for (int l = 0; l < prol_colind.Size(); l++)
		  {
		    int kk = prol_rowind[k];
		    int ll = prol_colind[l];
		    
		    if (kk >= ll) swap (kk,ll);
		    e2v.Append (INT<2> (kk,ll));
		  }
	    }

	int nc = 0;
	for (int i = 0; i < e2v.Size(); i++)
	  nc = max2 (nc, e2v[i][1]);
	nc++;

	// *testout << "e2v = " << endl << e2v << endl;
        
        // count all entries in row with multiplicity
	Array<int> cnt(nc);
	cnt = 0;
	for (int i = 0; i < e2v.Size(); i++)
	  cnt[e2v[i][1]]++;

	Table<int> v2e(cnt);
	cnt = 0;
	for (int i = 0; i < e2v.Size(); i++)
	  {
	    int v1 = e2v[i][1];
	    v2e[v1][cnt[v1]++] = i;
	  }
	
	cnt = 0;
	marks = -1;

        // count all entries in row withOUT multiplicity
	for (int i = 0; i < nc; i++)
	  for (int j = 0; j < v2e[i].Size(); j++)
	    {
	      int jj = v2e[i][j];
	      int v0 = e2v[jj][0];
	      if (marks[v0] != i) 
		{
		  cnt[i]++;
		  marks[v0] = i;
		}
	    }

	cmat = make_shared<SparseMatrixSymmetric<TM,TV>> (cnt);

	marks = -1;
	for (int i = 0; i < nc; i++)
	  for (int j = 0; j < v2e[i].Size(); j++)
	    {
	      int jj = v2e[i][j];
	      int v0 = e2v[jj][0];
	      if (marks[v0] != i) 
		{
		  marks[v0] = i;
		  cmat -> CreatePosition (i, v0);
		}
	    }
      }

    *cmat = 0.;
    RegionTimer reg2(tcomp);
	  
    for (int i = 0; i < n; i++)
      {
        FlatArray<int> mat_ri = this->GetRowIndices(i);
        FlatVector<TM> mat_rval = this->GetRowValues(i);

        for (int j = 0; j < mat_ri.Size(); j++)
          {
            int col = mat_ri[j];
            TM mat_val = mat_rval[j]; 

            FlatArray<int> prol_ri_i = prol.GetRowIndices(i);
            FlatArray<int> prol_ri_col = prol.GetRowIndices(col);
            FlatVector<double> prol_rval_i = prol.GetRowValues(i);
            FlatVector<double> prol_rval_col = prol.GetRowValues(col);

            for (int k = 0; k < prol_ri_i.Size(); k++)
              for (int l = 0; l < prol_ri_col.Size(); l++)
                {
                  int kk = prol_ri_i[k];
                  int ll = prol_ri_col[l];
                  
                  if ( kk>=ll && kk < cmat->Height() )
                    {
                      (*cmat)(kk,ll) += 
                        prol_rval_i[k] * prol_rval_col[l] * mat_val; 
                    }
                  
                  if (ll >= kk && i != col && ll < cmat->Height() )
                    {
                      (*cmat)(ll,kk) += 
                        prol_rval_col[l] * prol_rval_i[k] * Trans(mat_val); 
                    }
                  
                }
          }
      }
    return cmat;
  }
  



  template <> shared_ptr<BaseSparseMatrix>
  SparseMatrix<double> :: Restrict (const SparseMatrixTM<double> & prol,
                                    shared_ptr<BaseSparseMatrix> acmat ) const
  {
    static Timer t ("sparsematrix - restrict");
    RegionTimer reg(t);

    // auto prolT = TransposeMatrix(prol);
    auto prolT = dynamic_pointer_cast<SparseMatrixTM<double>> (prol.CreateTranspose());

    auto prod1 = MatMult<double, double, double>(*this, prol);
    auto prod = MatMult<double, double, double>(*prolT, *prod1);
    return prod;
  }

  template <> shared_ptr<BaseSparseMatrix>
  SparseMatrix<std::complex<double>> :: Restrict (const SparseMatrixTM<double> & prol,
                                                  shared_ptr<BaseSparseMatrix> acmat ) const
  {
    static Timer t ("sparsematrix - restrict");
    RegionTimer reg(t);
    // new version
    // auto prolT = TransposeMatrix(prol);
    auto prolT = dynamic_pointer_cast<SparseMatrixTM<double>> (prol.CreateTranspose());
    
    auto prod1 = MatMult<std::complex<double>, std::complex<double>, double>(*this, prol);
    auto prod = MatMult<std::complex<double>, double, std::complex<double>>(*prolT, *prod1);
    return prod;
  }





  template <> shared_ptr<BaseSparseMatrix>
  SparseMatrixSymmetric<double,double> :: Restrict (const SparseMatrixTM<double> & prol,
                                                    shared_ptr<BaseSparseMatrix> acmat ) const
  {
    static Timer t ("sparsematrixsymmetric - restrict");
    RegionTimer reg(t);
    // new version
    // auto prolT = TransposeMatrix(prol);
    auto prolT = dynamic_pointer_cast<SparseMatrixTM<double>> (prol.CreateTranspose());    
    auto full = MakeFullMatrix(*this);

    auto prod1 = MatMult<double, double, double>(*full, prol);
    auto prod = MatMult<double, double, double>(*prolT, *prod1);

    auto prodhalf = GetSymmetricMatrix (*prod);
    return prodhalf;


#ifdef OLD

    static Timer t ("sparsematrix - restrict");
    static Timer tbuild ("sparsematrix - restrict, build matrix");
    static Timer tbuild1 ("sparsematrix - restrict, build matrix1");
    static Timer tbuild2 ("sparsematrix - restrict, build matrix2");
    static Timer tbuild3 ("sparsematrix - restrict, build matrix3");
    static Timer tbuild4 ("sparsematrix - restrict, build matrix4");
    static Timer tcomp ("sparsematrix - restrict, compute matrix");
    RegionTimer reg(t);
    
    int n = this->Height();

    SparseMatrixSymmetric<double,double>* cmat = 
      dynamic_cast< SparseMatrixSymmetric<double,double>* > ( acmat );
 
    // if no coarse matrix, build up matrix-graph!
    if ( !cmat )
      {
        RegionTimer reg(tbuild);
        tbuild1.Start();
	Array<INT<2> > e2v;
	for (int i = 0; i < n; i++)
	  for (int j = 0; j < this->GetRowIndices(i).Size(); j++)
	    {
	      int col = this->GetRowIndices(i)[j];

	      for (int kk : prol.GetRowIndices(i))
		for (int ll : prol.GetRowIndices(col))
		  {
                    /*
		    if (kk >= ll) swap (kk,ll);
		    e2v.Append (INT<2> (kk,ll));
                    */
                    INT<2> i2 = (kk<ll) ? INT<2>(kk,ll) : INT<2>(ll,kk);
                    e2v.Append (i2);
		  }
	    }

	int nc = 0;
        for (auto vpair : e2v)
          nc = max2 (nc, vpair[1]);
	nc++;

        tbuild1.Stop();
        tbuild2.Start();

        // count all entries in row with multiplicity
	Array<int> cnt(nc);
	cnt = 0;
	for (int i = 0; i < e2v.Size(); i++)
	  cnt[e2v[i][1]]++;

	Table<int> v2e(cnt);
	cnt = 0;
	for (int i = 0; i < e2v.Size(); i++)
	  {
	    int v1 = e2v[i][1];
	    v2e[v1][cnt[v1]++] = i;
	  }

        tbuild2.Stop();
        tbuild3.Start();
	
	cnt = 0;
	Array<int> marks(n);
	marks = -1;

        // count all entries in row withOUT multiplicity
	for (int i = 0; i < nc; i++)
	  for (int j = 0; j < v2e[i].Size(); j++)
	    {
	      int jj = v2e[i][j];
	      int v0 = e2v[jj][0];
	      if (marks[v0] != i) 
		{
		  cnt[i]++;
		  marks[v0] = i;
		}
	    }

        tbuild3.Stop();
        tbuild4.Start();

	cmat = new SparseMatrixSymmetric<double,double> (cnt);

	marks = -1;
	for (int i = 0; i < nc; i++)
	  for (int j = 0; j < v2e[i].Size(); j++)
	    {
	      int jj = v2e[i][j];
	      int v0 = e2v[jj][0];
	      if (marks[v0] != i) 
		{
		  marks[v0] = i;
		  cmat -> CreatePosition (i, v0);
		}
	    }

        tbuild4.Stop();
      }

    *cmat = 0.;
    RegionTimer reg2(tcomp);

    // #pragma omp parallel for	  
    for (int i = 0; i < n; i++)
      {
        FlatArray<int> mat_ri = this->GetRowIndices(i);
        FlatVector<double> mat_rval = this->GetRowValues(i);

        for (int j = 0; j < mat_ri.Size(); j++)
          {
            int col = mat_ri[j];
            double mat_val = mat_rval[j]; 

            FlatArray<int> prol_rowind = prol.GetRowIndices(i);
            FlatArray<int> prol_colind = prol.GetRowIndices(col);
            FlatVector<double> prol_rowval = prol.GetRowValues(i);
            FlatVector<double> prol_colval = prol.GetRowValues(col);

            for (int k = 0; k < prol_rowind.Size(); k++)
              for (int l = 0; l < prol_colind.Size(); l++)
                {
                  int kk = prol_rowind[k];
                  int ll = prol_colind[l];
                  
                  if ( kk>=ll && kk < cmat->Height() )
                    {
                      auto val = prol_rowval[k] * prol_colval[l] * mat_val;
                      // #pragma omp atomic                      
                      (*cmat)(kk,ll) += val;
                    }
                  
                  if (ll >= kk && i != col && ll < cmat->Height() )
                    {
                      auto val = prol_colval[l] * prol_rowval[k] * Trans(mat_val);
                      // #pragma omp atomic                      
                      (*cmat)(ll,kk) += val;
                    }
                  
                }
          }
      }
    return cmat;
#endif

  }







  

  template <class TM, class TV>
  shared_ptr<BaseMatrix> SparseMatrixSymmetric<TM,TV> :: InverseMatrix (shared_ptr<BitArray> subset) const
  {
    if ( this->GetInverseType() == SUPERLU_DIST )
      throw Exception ("SparseMatrix::InverseMatrix:  SuperLU_DIST_Inverse not available");

    if (  BaseSparseMatrix :: GetInverseType()  == SUPERLU )
      {
#ifdef USE_SUPERLU
	return new SuperLUInverse<TM,TV_ROW,TV_COL> (*this, subset, nullptr, 1);
#else
	throw Exception ("SparseMatrix::InverseMatrix:  SuperLUInverse not available");
#endif
      }
    else if ( BaseSparseMatrix :: GetInverseType()  == PARDISO ||  BaseSparseMatrix :: GetInverseType()  == PARDISOSPD)
      {
        if(is_pardiso_available)
          return make_shared<PardisoInverse<TM,TV_ROW,TV_COL>> (*this, subset, nullptr, 1);
        else
          throw Exception ("SparseMatrix::InverseMatrix:  PardisoInverse not available");
      }
    else if (  BaseSparseMatrix :: GetInverseType()  == UMFPACK)
      {
#ifdef USE_UMFPACK
	return make_shared<UmfpackInverse<TM,TV_ROW,TV_COL>> (*this, subset, nullptr, 1);
#else
	throw Exception ("SparseMatrix::InverseMatrix:  UmfpackInverse not available");
#endif
      }
    else if ( BaseSparseMatrix :: GetInverseType()  == MUMPS )
      {
#ifdef USE_MUMPS
	return make_shared<MumpsInverse<TM,TV_ROW,TV_COL>> (*this, subset, nullptr, 1);
#else
	throw Exception ("SparseMatrix::InverseMatrix:  MumpsInverse not available");
#endif
      }
    else
      return make_shared<SparseCholesky<TM,TV_ROW,TV_COL>> (*this, subset);
  }

  template <class TM, class TV>
  shared_ptr<BaseMatrix> SparseMatrixSymmetric<TM,TV> :: InverseMatrix (shared_ptr<const Array<int>> clusters) const
  {
    if ( this->GetInverseType() == SUPERLU_DIST )
      throw Exception ("SparseMatrix::InverseMatrix:  SuperLU_DIST_Inverse not available");

    if (  BaseSparseMatrix :: GetInverseType()  == SUPERLU )
      {
#ifdef USE_SUPERLU
	return make_shared<SuperLUInverse<TM,TV_ROW,TV_COL>> (*this, nullptr, clusters, 1);
#else
	throw Exception ("SparseMatrix::InverseMatrix:  SuperLUInverse not available");
#endif
      }
    else if (  BaseSparseMatrix :: GetInverseType()  == PARDISO ||  BaseSparseMatrix :: GetInverseType()  == PARDISOSPD)
      {
        if(is_pardiso_available)
          return make_shared<PardisoInverse<TM,TV_ROW,TV_COL>> (*this, nullptr, clusters, 1);
        else
          throw Exception ("SparseMatrix::InverseMatrix:  PardisoInverse not available");
      }
    else if (  BaseSparseMatrix :: GetInverseType()  == UMFPACK)
      {
#ifdef USE_UMFPACK
	return make_shared<UmfpackInverse<TM,TV_ROW,TV_COL>> (*this, nullptr, clusters, 1);
#else
	throw Exception ("SparseMatrix::InverseMatrix:  UmfpackInverse not available");
#endif
      }
    else if (  BaseSparseMatrix :: GetInverseType()  == MUMPS )
      {
#ifdef USE_MUMPS
	return make_shared<MumpsInverse<TM,TV_ROW,TV_COL>> (*this, nullptr, clusters, 1);
#else
	throw Exception ("SparseMatrix::InverseMatrix:  MumpsInverse not available");
#endif
      }
    else
      return make_shared<SparseCholesky<TM,TV_ROW,TV_COL>> (*this, nullptr, clusters);
  }


  template <>
  void SparseMatrix<double,double,double> ::
  MultAdd (FlatVector<double> alpha, const MultiVector & x, MultiVector & y) const
  {
    // BaseMatrix::MultAdd (alpha, x, y);

    static Timer t("SparseMatrix::MultAdd Multivec"); RegionTimer reg(t);
    t.AddFlops (this->NZE()*x.Size());

    /*
    task_manager -> CreateJob
      ([&] (TaskInfo & ti)
       {
         int tasks_per_part = ti.ntasks / balance.Size();
         int mypart = ti.task_nr / tasks_per_part;
         int num_in_part = ti.task_nr % tasks_per_part;

         auto myrange = balance[mypart].Split (num_in_part, tasks_per_part);
         cout << "myrange = " << myrange << endl;
         cout << "tasknr = " << ti.task_nr << "  ntasks = " << ti.ntasks << endl;
    */

    ParallelForRange
      (balance, [&] (IntRange myrange)
       {
         int i = 0;
         for ( ; i + 4 <= x.Size(); i += 4)
           {
             auto fx0 = x[i+0]->FVDouble();
             auto fx1 = x[i+1]->FVDouble();
             auto fx2 = x[i+2]->FVDouble();
             auto fx3 = x[i+3]->FVDouble();
             auto fy0 = y[i+0]->FVDouble();
             auto fy1 = y[i+1]->FVDouble();
             auto fy2 = y[i+2]->FVDouble();
             auto fy3 = y[i+3]->FVDouble();
             double a0 = alpha[i+0];
             double a1 = alpha[i+1];
             double a2 = alpha[i+2];
             double a3 = alpha[i+3];
             for (auto row : myrange)
               {
                 double sum0 = 0;
                 double sum1 = 0;
                 double sum2 = 0;
                 double sum3 = 0;
                 for (size_t j = firsti[row]; j < firsti[row+1]; j++)
                   {
                     sum0 += data[j] * fx0(colnr[j]);
                     sum1 += data[j] * fx1(colnr[j]);
                     sum2 += data[j] * fx2(colnr[j]);
                     sum3 += data[j] * fx3(colnr[j]);
                   }
                 fy0(row) += a0 * sum0;
                 fy1(row) += a1 * sum1;
                 fy2(row) += a2 * sum2;
                 fy3(row) += a3 * sum3;
               }
           }

         for ( ; i+1 <= x.Size(); i++)
           {
             auto fx0 = x[i+0]->FVDouble();
             auto fy0 = y[i+0]->FVDouble();
             double a0 = alpha[i+0];
             for (auto row : myrange)
               {
                 double sum0 = 0;
                 for (size_t j = firsti[row]; j < firsti[row+1]; j++)
                   {
                     sum0 += data[j] * fx0(colnr[j]);
                   }
                 fy0(row) += a0 * sum0;
               }
           }

       });
  }


  //  template class SparseMatrix<double>;
#define NONExx
#ifdef NONExx


  template class SparseMatrixTM<double>;
  template class SparseMatrixTM<Complex>;
  template class SparseMatrix<double>;
  template class SparseMatrix<Complex>;
  template class SparseMatrix<double, Complex, Complex>;
  template class SparseMatrixSymmetric<double>;
  template class SparseMatrixSymmetric<Complex>;
  template class SparseMatrixSymmetric<double, Complex>;


#define INST_A_SPM(N, M, SCAL)			  \
  template class SparseMatrixTM<Mat<N, M, SCAL>>; \
  template class SparseMatrix<Mat<N, M, SCAL>>;	  \

#define INST_SPMS(N)					      \
  INST_A_SPM(N, N, double);				      \
  INST_A_SPM(1, N, double);				      \
  INST_A_SPM(N, 1, double);				      \
  INST_A_SPM(N, N, Complex);				      \
  INST_A_SPM(1, N, Complex);				      \
  INST_A_SPM(N, 1, Complex);				      \
  template class SparseMatrixSymmetric<Mat<N, N, double>>;    \
  template class SparseMatrixSymmetric<Mat<N, N, Complex>>;

  /** Do this case by hand - otherwise duplicate instantiation N,N and 1,N etc**/
  INST_A_SPM(1, 1, double);
  INST_A_SPM(1, 1, Complex);
  template class SparseMatrixSymmetric<Mat<1, 1, double>>;	\
  template class SparseMatrixSymmetric<Mat<1, 1, Complex>>;

  
#if MAX_SYS_DIM >= 2
  INST_SPMS(2);
#endif
#if MAX_SYS_DIM >= 3
  INST_SPMS(3);
#endif
#if MAX_SYS_DIM >= 4
  INST_SPMS(4);
#endif
#if MAX_SYS_DIM >= 5
  INST_SPMS(5);
#endif
#if MAX_SYS_DIM >= 6
  INST_SPMS(6);
#endif
#if MAX_SYS_DIM >= 7
  INST_SPMS(7);
#endif
#if MAX_SYS_DIM >= 8
  INST_SPMS(8);
#endif

  /**
#if MAX_SYS_DIM >= 1
  template class SparseMatrixTM<Mat<1,1,double> >;
  template class SparseMatrixTM<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class SparseMatrixTM<Mat<2,2,double> >;
  template class SparseMatrixTM<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class SparseMatrixTM<Mat<3,3,double> >;
  template class SparseMatrixTM<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class SparseMatrixTM<Mat<4,4,double> >;
  template class SparseMatrixTM<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class SparseMatrixTM<Mat<5,5,double> >;
  template class SparseMatrixTM<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class SparseMatrixTM<Mat<6,6,double> >;
  template class SparseMatrixTM<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class SparseMatrixTM<Mat<7,7,double> >;
  template class SparseMatrixTM<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class SparseMatrixTM<Mat<8,8,double> >;
  template class SparseMatrixTM<Mat<8,8,Complex> >;
#endif
  **/

  /*
  template class SparseMatrixSymmetricTM<double>;
  template class SparseMatrixSymmetricTM<Complex>;


#if MAX_SYS_DIM >= 1
  template class SparseMatrixSymmetricTM<Mat<1,1,double> >;
  template class SparseMatrixSymmetricTM<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class SparseMatrixSymmetricTM<Mat<2,2,double> >;
  template class SparseMatrixSymmetricTM<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class SparseMatrixSymmetricTM<Mat<3,3,double> >;
  template class SparseMatrixSymmetricTM<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class SparseMatrixSymmetricTM<Mat<4,4,double> >;
  template class SparseMatrixSymmetricTM<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class SparseMatrixSymmetricTM<Mat<5,5,double> >;
  template class SparseMatrixSymmetricTM<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class SparseMatrixSymmetricTM<Mat<6,6,double> >;
  template class SparseMatrixSymmetricTM<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class SparseMatrixSymmetricTM<Mat<7,7,double> >;
  template class SparseMatrixSymmetricTM<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class SparseMatrixSymmetricTM<Mat<8,8,double> >;
  template class SparseMatrixSymmetricTM<Mat<8,8,Complex> >;
#endif
  */








  /**
#if MAX_SYS_DIM >= 1
  template class SparseMatrix<Mat<1,1,double> >;
  template class SparseMatrix<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class SparseMatrix<Mat<2,2,double> >;
  template class SparseMatrix<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class SparseMatrix<Mat<3,3,double> >;
  template class SparseMatrix<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class SparseMatrix<Mat<4,4,double> >;
  template class SparseMatrix<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class SparseMatrix<Mat<5,5,double> >;
  template class SparseMatrix<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class SparseMatrix<Mat<6,6,double> >;
  template class SparseMatrix<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class SparseMatrix<Mat<7,7,double> >;
  template class SparseMatrix<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class SparseMatrix<Mat<8,8,double> >;
  template class SparseMatrix<Mat<8,8,Complex> >;
#endif


  template class SparseMatrixSymmetric<double>;
  template class SparseMatrixSymmetric<Complex>;
  template class SparseMatrixSymmetric<double, Complex>;


#if MAX_SYS_DIM >= 1
  template class SparseMatrixSymmetric<Mat<1,1,double> >;
  template class SparseMatrixSymmetric<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class SparseMatrixSymmetric<Mat<2,2,double> >;
  template class SparseMatrixSymmetric<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class SparseMatrixSymmetric<Mat<3,3,double> >;
  template class SparseMatrixSymmetric<Mat<3,3,Complex> >;
#endif

#if MAX_SYS_DIM >= 4
  template class SparseMatrixSymmetric<Mat<4,4,double> >;
  template class SparseMatrixSymmetric<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class SparseMatrixSymmetric<Mat<5,5,double> >;
  template class SparseMatrixSymmetric<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class SparseMatrixSymmetric<Mat<6,6,double> >;
  template class SparseMatrixSymmetric<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class SparseMatrixSymmetric<Mat<7,7,double> >;
  template class SparseMatrixSymmetric<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class SparseMatrixSymmetric<Mat<8,8,double> >;
  template class SparseMatrixSymmetric<Mat<8,8,Complex> >;
#endif
  **/

#ifdef CACHEBLOCKSIZE
  template class SparseMatrixSymmetric<double, Vec<CACHEBLOCKSIZE> >;
#endif
  
#if MAX_CACHEBLOCKS >= 2
  template class SparseMatrixSymmetric<double, Vec<2> >;
#ifdef GOLD
  // template class SparseMatrixSymmetric<double, MD<2> >;
#endif
#endif
#if MAX_CACHEBLOCKS >= 3
  template class SparseMatrixSymmetric<double, Vec<3> >;
  template class SparseMatrixSymmetric<double, Vec<4> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class SparseMatrixSymmetric<double, Vec<5> >;
  template class SparseMatrixSymmetric<double, Vec<6> >;
  template class SparseMatrixSymmetric<double, Vec<7> >;
  template class SparseMatrixSymmetric<double, Vec<8> >;
  template class SparseMatrixSymmetric<double, Vec<9> >;
  template class SparseMatrixSymmetric<double, Vec<10> >;
  template class SparseMatrixSymmetric<double, Vec<11> >;
  template class SparseMatrixSymmetric<double, Vec<12> >;
  template class SparseMatrixSymmetric<double, Vec<13> >;
  template class SparseMatrixSymmetric<double, Vec<14> >;
  template class SparseMatrixSymmetric<double, Vec<15> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class SparseMatrixSymmetric<Complex, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class SparseMatrixSymmetric<Complex, Vec<3,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class SparseMatrixSymmetric<Complex, Vec<5,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<6,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<7,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<8,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<9,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<10,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<11,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<12,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<13,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<14,Complex> >;
  template class SparseMatrixSymmetric<Complex, Vec<15,Complex> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class SparseMatrixSymmetric<double, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class SparseMatrixSymmetric<double, Vec<3,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class SparseMatrixSymmetric<double, Vec<5,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<6,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<7,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<8,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<9,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<10,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<11,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<12,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<13,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<14,Complex> >;
  template class SparseMatrixSymmetric<double, Vec<15,Complex> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class SparseMatrix<double, Vec<2>, Vec<2> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class SparseMatrix<double, Vec<3>, Vec<3> >;
  template class SparseMatrix<double, Vec<4>, Vec<4> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class SparseMatrix<double, Vec<5>, Vec<5> >;
  template class SparseMatrix<double, Vec<6>, Vec<6> >;
  template class SparseMatrix<double, Vec<7>, Vec<7> >;
  template class SparseMatrix<double, Vec<8>, Vec<8> >;
  template class SparseMatrix<double, Vec<9>, Vec<9> >;
  template class SparseMatrix<double, Vec<10>, Vec<10> >;
  template class SparseMatrix<double, Vec<11>, Vec<11> >;
  template class SparseMatrix<double, Vec<12>, Vec<12> >;
  template class SparseMatrix<double, Vec<13>, Vec<13> >;
  template class SparseMatrix<double, Vec<14>, Vec<14> >;
  template class SparseMatrix<double, Vec<15>, Vec<15> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class SparseMatrix<Complex, Vec<2,Complex>, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class SparseMatrix<Complex, Vec<3,Complex>, Vec<3,Complex> >;
  template class SparseMatrix<Complex, Vec<4,Complex>, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class SparseMatrix<Complex, Vec<5,Complex>, Vec<5,Complex> >;
  template class SparseMatrix<Complex, Vec<6,Complex>, Vec<6,Complex> >;
  template class SparseMatrix<Complex, Vec<7,Complex>, Vec<7,Complex> >;
  template class SparseMatrix<Complex, Vec<8,Complex>, Vec<8,Complex> >;
  template class SparseMatrix<Complex, Vec<9,Complex>, Vec<9,Complex> >;
  template class SparseMatrix<Complex, Vec<10,Complex>, Vec<10,Complex> >;
  template class SparseMatrix<Complex, Vec<11,Complex>, Vec<11,Complex> >;
  template class SparseMatrix<Complex, Vec<12,Complex>, Vec<12,Complex> >;
  template class SparseMatrix<Complex, Vec<13,Complex>, Vec<13,Complex> >;
  template class SparseMatrix<Complex, Vec<14,Complex>, Vec<14,Complex> >;
  template class SparseMatrix<Complex, Vec<15,Complex>, Vec<15,Complex> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class SparseMatrix<double, Vec<2,Complex>, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class SparseMatrix<double, Vec<3,Complex>, Vec<3,Complex> >;
  template class SparseMatrix<double, Vec<4,Complex>, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class SparseMatrix<double, Vec<5,Complex>, Vec<5,Complex> >;
  template class SparseMatrix<double, Vec<6,Complex>, Vec<6,Complex> >;
  template class SparseMatrix<double, Vec<7,Complex>, Vec<7,Complex> >;
  template class SparseMatrix<double, Vec<8,Complex>, Vec<8,Complex> >;
  template class SparseMatrix<double, Vec<9,Complex>, Vec<9,Complex> >;
  template class SparseMatrix<double, Vec<10,Complex>, Vec<10,Complex> >;
  template class SparseMatrix<double, Vec<11,Complex>, Vec<11,Complex> >;
  template class SparseMatrix<double, Vec<12,Complex>, Vec<12,Complex> >;
  template class SparseMatrix<double, Vec<13,Complex>, Vec<13,Complex> >;
  template class SparseMatrix<double, Vec<14,Complex>, Vec<14,Complex> >;
  template class SparseMatrix<double, Vec<15,Complex>, Vec<15,Complex> >;
#endif
#endif

}
