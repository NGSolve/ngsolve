/**************************************************************************/
/* File:   table.cpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   25. Mar. 2000                                                  */
/**************************************************************************/

/* 
   Abstract data type Table
*/

#include <ngstd.hpp>

namespace ngstd
{

  /*
  BaseTable :: BaseTable (int asize, int entrysize)
  {
    size = asize;
    index = new int[size+1];
    for (int i = 0; i <= size; i++)
      index[i] = i*entrysize;
  }


  BaseTable :: BaseTable (const FlatArray<int> & entrysize)
  {
    int i, cnt = 0;
    size  = entrysize.Size();

    index = new int[size+1];
    for (i = 0; i < size; i++)
      {
	index[i] = cnt;
	cnt += entrysize[i];
      }
    index[i] = cnt;
  }
  */

  /*
  BaseTable :: ~BaseTable ()
  {
    delete [] index;
  }
  */

#ifdef PARALLEL_TABLE
  template <typename TI> 
  size_t * TablePrefixSum2 (FlatArray<TI> entrysize)
  {
    /*
    size_t size  = entrysize.Size();
    size_t * index = new size_t[size+1];

    size_t cnt = 0;
    for (size_t i = 0; i < size; i++)
      {
	index[i] = cnt;
	cnt += entrysize[i];
      }
    index[size] = cnt;
    return index;
    */

    size_t size  = entrysize.Size();
    size_t * index = new size_t[size+1];

    Array<size_t> partial_sums(TaskManager::GetNumThreads()+1);
    partial_sums[0] = 0;
    ParallelJob
      ([&] (TaskInfo ti)
       {
         IntRange r = IntRange(size).Split(ti.task_nr, ti.ntasks);
         size_t mysum = 0;
         for (size_t i : r)
           mysum += entrysize[i];
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
             index[i] = mysum;
             mysum += entrysize[i];
           }
       });
    index[size] = partial_sums.Last();

    return index;
  }

  DLL_HEADER size_t * TablePrefixSum32 (FlatArray<unsigned int> entrysize)
  { return TablePrefixSum2 (entrysize); }
  DLL_HEADER size_t * TablePrefixSum64 (FlatArray<size_t> entrysize)
  { return TablePrefixSum2 (entrysize); }
  /*
  DLL_HEADER template size_t * TablePrefixSum<int> (FlatArray<int> entrysize);
  DLL_HEADER template size_t * TablePrefixSum<unsigned int> (FlatArray<unsigned int> entrysize);
  DLL_HEADER template size_t * TablePrefixSum<size_t> (FlatArray<size_t> entrysize);
  DLL_HEADER template size_t * TablePrefixSum<atomic<int>> (FlatArray<atomic<int>> entrysize);
  */
#endif
  
  BaseDynamicTable :: BaseDynamicTable (int size)
    : data(size)
  {
    for (int i = 0; i < size; i++)
      {
	data[i].maxsize = 0;
	data[i].size = 0;
	data[i].col = NULL;
      }
    oneblock = NULL;
  }

  BaseDynamicTable :: BaseDynamicTable (const Array<int> & entrysizes, int elemsize)
    : data(entrysizes.Size())
  {
    int cnt = 0;
    int n = entrysizes.Size();
    
    for (int i = 0; i < n; i++)
      cnt += entrysizes[i];
    oneblock = new char[elemsize * cnt];

    cnt = 0;
    for (int i = 0; i < n; i++)
      {
	data[i].maxsize = entrysizes[i];
	data[i].size = 0;

	data[i].col = &oneblock[elemsize * cnt];
	cnt += entrysizes[i];
      }
  }


  BaseDynamicTable :: ~BaseDynamicTable ()
  {
    if (oneblock)
      delete [] oneblock;
    else
      for (int i = 0; i < data.Size(); i++)
	delete [] static_cast<char*> (data[i].col);
  }

  void BaseDynamicTable :: SetSize (int size)
  {
    for (int i = 0; i < data.Size(); i++)
      delete [] static_cast<char*> (data[i].col);

    data.SetSize(size);
    for (int i = 0; i < size; i++)
      {
	data[i].maxsize = 0;
	data[i].size = 0;
	data[i].col = NULL;
      }    
  }

  void BaseDynamicTable :: IncSize (int i, int elsize)
  {
    if (i < 0 || i >= data.Size())
      {
	cerr << "BaseDynamicTable::Inc: Out of range, i = " << i << ", size = " << data.Size() << endl;
	return;
      }
    
    linestruct & line = data[i];
    
    if (line.size == line.maxsize)
      {
	void * p = new char [(2*line.maxsize+5) * elsize];
      
	memcpy (p, line.col, line.maxsize * elsize);
	delete [] static_cast<char*> (line.col);
	line.col = p;
	line.maxsize = 2*line.maxsize+5;
      }
  
    line.size++;
  }

  void BaseDynamicTable :: DecSize (int i)
  {
    if (i < 0 || i >= data.Size())
      {
	cerr << "BaseDynamicTable::Dec: Out of range" << endl;
	return;
      }
  
    linestruct & line = data[i];
  
    if (line.size == 0)
      {
	cerr << "BaseDynamicTable::Dec: EntrySize < 0" << endl;
	return;      
      }
  
    line.size--;
  }
  
  void FilteredTableCreator::Add (size_t blocknr, int data)
  {
    if (!takedofs||takedofs->Test(data))
      TableCreator<int>::Add(blocknr,data);
  }

  void FilteredTableCreator::Add (size_t blocknr, IntRange range)
  {
    for (size_t i=range.First(); i<range.Next();i++)
      if (!takedofs||takedofs->Test(i))
	TableCreator<int>::Add(blocknr,i);
  }  
  
  void FilteredTableCreator::Add (size_t blocknr, FlatArray<int> dofs)
  {
    for (size_t i = 0; i < dofs.Size(); i++)
      if (!takedofs||takedofs->Test(dofs[i]))
	TableCreator<int>::Add(blocknr,dofs[i]);
  }  
}
