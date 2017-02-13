
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

    auto GetVertexOrientedEdge (int nr) const
    {
      return ET_trait<ET>::GetEdgeSort (nr, vnums);
    }
    auto GetVertexOrientedFace (int nr) const
    {
      return ET_trait<ET>::GetFaceSort (nr, vnums);
    }
  };
}

#endif // FILE_FEINTERFACES_
