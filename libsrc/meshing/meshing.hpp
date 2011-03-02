#ifndef FILE_MESHING
#define FILE_MESHING



#include "../include/myadt.hpp"
#include "../include/gprim.hpp"
#include "../include/linalg.hpp"
#include "../include/opti.hpp"



namespace netgen
{
  // extern int printmessage_importance;

  class CSGeometry;
  class NetgenGeometry;
}
  
  
#include "msghandler.hpp"
#include "meshtype.hpp"
#include "localh.hpp"
#include "meshclass.hpp"
#include "global.hpp"


namespace netgen
{
#include "meshtool.hpp"
#include "ruler2.hpp"
#include "adfront2.hpp"
#include "meshing2.hpp"
#include "improve2.hpp"


#include "geomsearch.hpp"
#include "adfront3.hpp"
#include "ruler3.hpp"

#define _INCLUDE_MORE


#include "meshing3.hpp"
#include "improve3.hpp"

#include "findip.hpp"
#include "findip2.hpp"

#include "topology.hpp"
#include "curvedelems.hpp"
#include "clusters.hpp"

#include "meshfunc.hpp"

#include "bisect.hpp"
#include "hprefinement.hpp"
#include "boundarylayer.hpp"
#include "specials.hpp"
}

#include "validate.hpp"
#include "basegeom.hpp"

#ifdef PARALLEL
#include "paralleltop.hpp"
#endif


#endif
