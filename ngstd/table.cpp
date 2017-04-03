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
