NGX_INLINE DLL_HEADER Ng_Point Ngx_Mesh :: GetPoint (int nr) const
{
  return Ng_Point (&mesh->Point(nr + PointIndex::BASE)(0));
}


template <>
NGX_INLINE DLL_HEADER int Ngx_Mesh :: GetElementIndex<1> (int nr) const
{
  return (*mesh)[SegmentIndex(nr)].si;
}
  
template <>
NGX_INLINE DLL_HEADER int Ngx_Mesh :: GetElementIndex<2> (int nr) const
{
  int ind = (*mesh)[SurfaceElementIndex(nr)].GetIndex(); 
  return mesh->GetFaceDescriptor(ind).BCProperty();
}

template <>
NGX_INLINE DLL_HEADER int Ngx_Mesh :: GetElementIndex<3> (int nr) const
{
  return (*mesh)[ElementIndex(nr)].GetIndex();
}

template <> 
NGX_INLINE DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<1> (int nr) const
{
  const Segment & el = mesh->LineSegment (SegmentIndex(nr));

  Ng_Element ret;
  ret.type = NG_ELEMENT_TYPE(el.GetType());

  ret.points.num = el.GetNP();
  ret.points.ptr = (int*)&(el[0]);

  ret.vertices.num = 2;
  ret.vertices.ptr = (int*)&(el[0]);

  ret.edges.num = 1;
  ret.edges.ptr = (T_EDGE2*)mesh->GetTopology().GetSegmentElementEdgesPtr (nr);

  ret.faces.num = 0;
  ret.faces.ptr = NULL;

  return ret;
}

template <> 
NGX_INLINE DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<2> (int nr) const
{
  const Element2d & el = mesh->SurfaceElement (SurfaceElementIndex (nr));
  
  Ng_Element ret;
  ret.type = NG_ELEMENT_TYPE(el.GetType());
  ret.points.num = el.GetNP();
  ret.points.ptr  = (int*)&el[0];

  ret.vertices.num = el.GetNV();
  ret.vertices.ptr = (int*)&(el[0]);

  ret.edges.num = MeshTopology::GetNEdges (el.GetType());
  ret.edges.ptr = (T_EDGE2*)mesh->GetTopology().GetSurfaceElementEdgesPtr (nr);

  ret.faces.num = MeshTopology::GetNFaces (el.GetType());
  ret.faces.ptr = (T_FACE2*)mesh->GetTopology().GetSurfaceElementFacesPtr (nr);

  return ret;
}

template <> 
NGX_INLINE DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<3> (int nr) const
{
  const Element & el = mesh->VolumeElement (ElementIndex (nr));
  
  Ng_Element ret;
  ret.type = NG_ELEMENT_TYPE(el.GetType());
  ret.points.num = el.GetNP();
  ret.points.ptr = (int*)&el[0];

  ret.vertices.num = el.GetNV();
  ret.vertices.ptr = (int*)&(el[0]);

  ret.edges.num = MeshTopology::GetNEdges (el.GetType());
  ret.edges.ptr = (T_EDGE2*)mesh->GetTopology().GetElementEdgesPtr (nr);

  ret.faces.num = MeshTopology::GetNFaces (el.GetType());
  ret.faces.ptr = (T_FACE2*)mesh->GetTopology().GetElementFacesPtr (nr);

  return ret;
}
