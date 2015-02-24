/*
  Standard data types demos
 */


// ng-soft header files
#include <ngstd.hpp>
using namespace ngstd;


int main ()
{
  // ****************** Using Arrays *********************

  // an array with memory allocation/deallocation
  Array<int> ia(10);

  ia = -1;
  for (int i = 0; i < 5; i++)
    ia[i] = i;
  cout << "array = " << ia << endl;

  // an array without memory handling, initialize with size and pointer
  FlatArray<int> sub(3, &ia[2]);
  cout << "sub-array = " << sub << endl;

  // access a semi-open sub-array [first,next)
  cout << "sub-array = " << ia.Range (2,5) << endl;



  // ******************** Tables ****************************
  
  DynamicTable<int> tab(10);
  tab.Add (1, 5);
  tab.Add (1, 8);
  tab.Add (1, 5);
  tab.Add (7, 3);

  cout << "table = " << endl << tab << endl;

  // a compact table requires a priori knowledge of entrysizes
  // copy dynamic table into a compact table
  Array<int> entry_sizes(tab.Size());
  for (int i = 0; i < tab.Size(); i++)
    entry_sizes[i] = tab[i].Size();

  Table<int> tab2(entry_sizes);
  for (int i = 0; i < tab.Size(); i++)
    for (int j = 0; j < tab[i].Size(); j++)
      tab2[i][j] = tab[i][j];
  
  cout << "table2 = " << endl << tab2 << endl;

  // ******************** Local Heaps ************************

  // initialize heap memory handler with 1000 bytes
  LocalHeap lh(1000, "demo - localheap");

  cout << "available: " << lh.Available() << endl;

  int * ip = new (lh) int [5];
  cout << "available: " << lh.Available() << endl;
  
  // stores pointer to current position
  void * heap_pointer = lh.GetPointer();

  int * ip2 = lh.Alloc<int> (5);
  cout << "available: " << lh.Available() << endl;

  // resets heap pointer to stored position -> ip2 array is invalid
  lh.CleanUp (heap_pointer);
  cout << "available: " << lh.Available() << endl;

  try 
    {
      int * ip = lh.Alloc<int>(500);
    }
  catch (Exception & e)
    {
      cout << "Exception " << e.What() << endl;
    }

  cout << "End of Tests" << endl;
  return 0;
}
