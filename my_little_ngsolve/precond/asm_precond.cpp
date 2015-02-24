/*

Additive / Multiplicative two-level Schwarz preconditioner

To define such a preconditioner, the finite element space 
has to provide
* the blocks used for an additive Schwartz smoother
* the dofs which build the coarse system

*/


#include <comp.hpp>
using namespace ngcomp;



class MyH1HighOrderFESpace : public H1HighOrderFESpace
{
public:
  MyH1HighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & flags)
    : H1HighOrderFESpace (ama, flags)
  {
    cout << "H1HighOrderSpace with smoothing-blocks tutorial" << endl;
  }

  virtual Table<int> * CreateSmoothingBlocks (const Flags & flags) const
  {
    cout << "My CreateSmoothingBlocks" << endl;

    /*
      We have to return a table representing the block-smoother.
      Each table-entry contains all the dofs of one block.
      A dof may belong to several blocks. This leads to an 
      block-smoother with overlaps
    */

    
    // are we eliminating internal dofs (static condensation) ?
    bool elim = flags.GetDefineFlag("eliminate_internal");
    cout << "elim = " << elim << endl;


    /*
      A table-creator creates a table in compressed form.
      It first determines the size of the table, then the entry-sizes,
      and finally it stores the entries.
      To obtain that, we have to provide the desired 
      entries multiple times.

      A filtered table-creator also filters out Dirichlet-dofs, and
      eliminated internal dofs
     */

    FilteredTableCreator creator(GetFreeDofs(elim));

    for ( ; !creator.Done(); creator++)
      {

        if (flags.GetDefineFlag ("ebe"))  
          { 
            // element-by-element preconditioner
            // put all dofs of one element into one block

            Array<int> dnums;
            for (int i = 0; i < ma->GetNE(); i++)
              {
                GetDofNrs(i, dnums);
                for (int j = 0; j < dnums.Size(); j++)
                  creator.Add (i, dnums[j]);  
              }
          }
        
        // here you can put alternatives ...
        // else if ( ...)

        else
          {
            throw Exception ("\n\nMyH1HighOrderFESpace::CreateSmoothingBlocks\n"
                             "undefined smoothing blocks\n"
                             "Possibilities are:  -ebe\n");
          }
      }

    // *testout << "smoothing block: " << endl << *creator.GetTable() << endl;
    return creator.GetTable();
  }


  virtual Array<int> * CreateDirectSolverClusters (const Flags & flags) const
  {
    /*
      Dofs marked in direct-solve-clusters are solved for directly.
      We have to return an int-array, one entry per dof.
      A value 0 means that we do not use a direct solver for that dof.
      All dofs with the same positive value are taken into the same cluster.
    */

    Array<int> & clusters = * new Array<int>(GetNDof());

    // all lowest-order dofs (i.e. vertex dofs) are solved directly
    clusters = 0;
    clusters.Range (0, ma->GetNV()) = 1;

    // but note: take only non-Dirichlet dofs
    const BitArray & freedofs = *GetFreeDofs();
    for (int i = 0; i < clusters.Size(); i++)
      if (!freedofs[i]) clusters[i] = 0;

    return &clusters;
  }
};



static RegisterFESpace<MyH1HighOrderFESpace> initifes ("myh1ho");
