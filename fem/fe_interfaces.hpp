
#ifndef FILE_FEINTERFACES_
#define FILE_FEINTERFACES_

namespace ngfem
{
  template<ELEMENT_TYPE ET>
  class VertexOrientedFE
  {
  protected:
    enum { N_VERTEX = ET_trait<ET>::N_VERTEX };
    int vnums[N_VERTEX];

  public:
    template <typename TA>
    INLINE void SetVertexNumbers (const TA & avnums)
    { for (int i = 0; i < N_VERTEX; i++) vnums[i] = avnums[i]; }
    /// assign vertex number
    INLINE void SetVertexNumber (int nr, int vnum) { vnums[nr] = vnum; }
    INLINE FlatArray<const int> GetVertexNumbers () const { return { N_VERTEX, vnums }; }
    auto GetVertexOrientedEdge (int nr) const
    {
      return ET_trait<ET>::GetEdgeSort (nr, vnums);
    }
    auto GetVertexOrientedFace (int nr) const
    {
      return ET_trait<ET>::GetFaceSort (nr, vnums);
    }
    /// vertex-order class: bit k is set if vnums[i] > vnums[j] for the k-th pair i < j.
    /// Shape functions depend on the vertex numbers only through this class.
    /// -1 if two vertex numbers coincide, since the class then does not fix the orientation.
    int VertexOrderClass () const
    {
      int cls = 0, bit = 0;
      for (int i = 0; i < N_VERTEX; i++)
        for (int j = i+1; j < N_VERTEX; j++, bit++)
          {
            if (vnums[i] == vnums[j]) return -1;
            if (vnums[i] > vnums[j]) cls |= 1 << bit;
          }
      return cls;
    }
  };
}

#endif // FILE_FEINTERFACES_
