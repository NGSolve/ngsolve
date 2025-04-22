#ifndef FILE_COEFFICIENT_IMPL
#define FILE_COEFFICIENT_IMPL


namespace ngfem
{


  struct GenericIdentity
  {
    template <typename T> T operator() (T x) const { return x; }
    static string Name() { return  " "; }
    template <typename T> T Diff (T x) const { throw Exception("no GenericIdentity::Diff"); }        
    void DoArchive(Archive& ar) {}
  };

}

#endif
