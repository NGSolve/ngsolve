// ngscxx -shared sumup.cpp -lngcomp -o mymip.so

#include <solve.hpp>
using namespace ngsolve;

/*
  
we sum-up the number of elements per vertex

if the mesh is distributed, we have to accumulate the sum across processes

*/


class NumProcSumUp : public NumProc
{
protected:
  shared_ptr<GridFunction> gfu;

public:

  NumProcSumUp (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    gfu = pde.GetGridFunction (flags.GetStringFlag ("gridfunction", "u"));
  }

  // solve at one level
  virtual void Do(LocalHeap & lh)
  {
    static Timer t("Numroc SumUp :: Do");
    RegionTimer reg(t);


    cout << "solve SumUp" << endl;

    FlatVector<> vec = gfu->GetVector().FV<double>();
    vec = 0.0;

    for (Ngs_Element el : ma->Elements(VOL))
      for (int v : el.Vertices())
	vec[v]++;
 



    Vector<> orig_vec = vec;

    // and now, we sum up across processes:


    // convenience function wrappers in ngstd/mpiwrapper.hpp
    int id    = MyMPI_GetId();  
    int ntask = MyMPI_GetNTasks();  
    
 
    // Version 1:
    /*
    // everyone sends to everyone
    // blocking communication
    // if me < you, first send, otherwise first receive
    for (int other = 0; other < ntask; other++)
      {
	if (other == id) continue;   // don't send to me

	Array<int> procs;
	Array<double> send_data;
	Array<double> recv_data;

	for (int i = 0; i < ma->GetNV(); i++)
	  {
	    ma->GetDistantProcs (Node(NT_VERTEX, i), procs);
	    if (procs.Contains(other))
	      send_data.Append (orig_vec[i]);
	  }
	
	if (other > id)
	  {
	    MyMPI_Send (send_data, other);
	    MyMPI_Recv (recv_data, other);
	  }
	else
	  {
	    MyMPI_Recv (recv_data, other);
	    MyMPI_Send (send_data, other);
	  }


	for (int i = 0, cnt = 0; i < ma->GetNV(); i++)
	  {
	    ma->GetDistantProcs (Node(NT_VERTEX, i), procs);
	    if (procs.Contains(other))
	      vec[i] += recv_data[cnt++];
	  }
      }
    */
  

    // Version 2:
    // post all communication jobs simultaneously:

    DynamicTable<int> send_table (ntask);
    DynamicTable<int> recv_table (ntask);

    for (int i = 0; i < ma->GetNV(); i++)
      for (int p : ma->GetDistantProcs (Node(NT_VERTEX, i)))
	{
	  send_table.Add (p, vec(i));
	  recv_table.Add (p, -1);  // dummy value 
	}

    Array<MPI_Request> requests;
    for (int other = 0; other < ntask; other++)
      if (other != id)
	{
	  requests += MyMPI_ISend (send_table[other], other);
	  requests += MyMPI_IRecv (recv_table[other], other);
	}

    // wait till all requests are fulfilled
    MyMPI_WaitAll (requests);

    Array<int> cnt(ntask); cnt = 0;

    for (int i = 0; i < ma->GetNV(); i++)
      for (int p : ma->GetDistantProcs (Node(NT_VERTEX, i)))
	vec(i) += recv_table[p][cnt[p]++];

    
    /*
    // Version 3:
    // a ngsolve-library function for this purpose:
    Array<double> data(vec.Size(), &vec(0));
    AllReduceNodalData (NT_VERTEX, data, MPI_SUM, *ma);
   */
  }

};



// register the numproc 'SumUp' 
static RegisterNumProc<NumProcSumUp> npinit1("SumUp");
