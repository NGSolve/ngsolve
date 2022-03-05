
/*
  
  Assemble the system matrix

  The input is
  - a finite element space, which provides the basis functions
  - an integrator, which computes the element matrices
  
  The result is a sparse matrix
*/

#include <comp.hpp>

namespace ngcomp
{
  shared_ptr<BaseSparseMatrix> MyAssembleMatrix(shared_ptr<FESpace> fes,
                                                shared_ptr<BilinearFormIntegrator> bfi,
                                                bool parallel)
  {
    cout << "We assemble matrix" << endl;
    Timer tgraph("MyAssemble - buildmatrixgraph");
    Timer tassemble("MyAssemble - assemble matrix");
    Timer telement("MyAssemble - calc element matrix");        
    auto ma = fes->GetMeshAccess();
    
    int ndof = fes->GetNDof();
    int ne = ma->GetNE();

    // we build a sparse matrix
    // the non-zero pattern is given by the connectivity pattern provided by the FESpace
      
    // setup element-to-dof table:
    // first we get the number of dofs per element ...
    tgraph.Start();
    Array<int> dnums;
    Array<int> cnt(ne);

    for (int i = 0; i < ma->GetNE(VOL); i++) 
      {
        fes->GetDofNrs (ElementId(VOL,i), dnums);
        cnt[i] = dnums.Size();
      }
      
    // allocate the table in compressed form ...
    Table<int> el2dof(cnt);

    // and fill it
    for (auto ei : ma->Elements(VOL))
      {
        fes->GetDofNrs (ei, dnums);
        el2dof[ei.Nr()] = dnums;
      }
    // cout << "el2dof - table: " << el2dof << endl;

    // generate sparse matrix of size ndof x ndof
    // from element-to-dof table for rows and columns
    auto mat = make_shared<SparseMatrix<double>> (ndof, ndof, el2dof, el2dof, false);
    mat -> SetZero();
    tgraph.Stop();
    tassemble.Start();
    
    if (!parallel)
      {
        cout << "sequential assembling" << endl;

        LocalHeap lh(1000*1000); // reserve 1MB 
        
        // loop over all volume elements
        for (int i = 0; i < ma->GetNE(VOL); i++)
          {
            HeapReset hr(lh);  // cleanup heap at end of scope
            ElementId ei(VOL, i);
            
            // let FESpace generate the finite element
            FiniteElement & fel = fes->GetFE (ei, lh);
            
            // the global dof numbers of the element
            fes->GetDofNrs (ei, dnums);
            
            // the mesh knows the geometry of the element
            const ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
            
            // compute the element matrix
            FlatMatrix<> elmat (fel.GetNDof(), lh);
            bfi->CalcElementMatrix (fel, eltrans, elmat, lh);
            
            mat->AddElementMatrix (dnums, dnums, elmat);
          }
      }
    else
      {
        cout << "parallel assembling" << endl;
        
        /*
          Parallel iteration over elements.
          The lambda-function will be called for all elements, 
          parallelized over available threads.
          
          To avoid race conditions, the FESpace provides a coloring, such that
          two elements of the same color have independent dof-numbers. These elements
          can be processed savely in parallel.
          
          The memory of the Localheap is split into peaces, one for each thread 
        */
        
        LocalHeap lh(1000000, "mylocalheap", true); // 1MB per thread
        IterateElements(*fes, VOL, lh, [&] (ElementId ei, LocalHeap &lh)
                        {
                          const ElementTransformation & eltrans = ma->GetTrafo(ei,lh);
                          const FiniteElement & fel = fes->GetFE(ei,lh);
                          
                          Array<int> dnums(fel.GetNDof(), lh);
                          fes->GetDofNrs (ei, dnums);
                          
                          telement.Start();
                          FlatMatrix<> elmat (fel.GetNDof(), lh);
                          bfi->CalcElementMatrix (fel, eltrans, elmat, lh);
                          telement.Stop();
                          
                          mat->AddElementMatrix (dnums, dnums, elmat);
                        });
      }
    
    tassemble.Stop();
    return mat;
  }
}
