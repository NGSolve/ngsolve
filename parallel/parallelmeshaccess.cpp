#ifdef PARALLEL
 
#include <comp.hpp>
#include <parallelngs.hpp>

#include <parallelinterface.hpp>


namespace ngparallel
{
  using namespace netgen;

  ParallelMeshAccessxxx :: ParallelMeshAccessxxx ( const MeshAccess & ama )
    : ma( ama )
  {
   ;
  }

  /*
  int ParallelMeshAccessxxx ::GetDistantNodeNum (int proc, NODE_TYPE nt, int locnum) const
  {
    switch (nt)
      {
      case NT_VERTEX: return NgPar_GetDistantPNum ( proc, locnum ); 
      case NT_EDGE: return NgPar_GetDistantEdgeNum ( proc, locnum ); 
      case NT_FACE: return NgPar_GetDistantFaceNum ( proc, locnum ); 
      case NT_CELL: return NgPar_GetDistantElNum ( proc, locnum ); 
      }
    return -1;
  }
  */

  int ParallelMeshAccessxxx ::GetGlobalNodeNum (NODE_TYPE nt, int locnum) const
  {
    switch (nt)
      {
      case NT_VERTEX: return NgPar_GetDistantPNum ( 0, locnum ); 
      case NT_EDGE: return NgPar_GetDistantEdgeNum ( 0, locnum ); 
      case NT_FACE: return NgPar_GetDistantFaceNum ( 0, locnum ); 
      case NT_CELL: return NgPar_GetDistantElNum ( 0, locnum ); 
      }
    return -1;
  }


  int ParallelMeshAccessxxx :: GetDistantNodeNums ( NODE_TYPE nt, int locnum, 
                                                 ngstd::Array<int[2]> & distnums ) const
  {
    distnums.SetSize( NgPar_GetNDistantNodeNums(nt, locnum) );
    return NgPar_GetDistantNodeNums ( nt, locnum, &distnums[0][0] );
  }
}


#endif
