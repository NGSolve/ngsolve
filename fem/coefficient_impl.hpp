#ifndef FILE_COEFFICIENT_IMPL
#define FILE_COEFFICIENT_IMPL


namespace ngfem
{


  struct GenericIdentity
  {
    template <typename T> T operator() (T x) const { return x; }
    static string Name() { return  " "; }
    void DoArchive(Archive& ar) {}
  };

}

#endif
