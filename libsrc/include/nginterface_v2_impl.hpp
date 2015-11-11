NGX_INLINE DLL_HEADER Ng_Point Ngx_Mesh :: GetPoint (int nr) const
{
  return Ng_Point (&mesh->Point(nr + PointIndex::BASE)(0));
}


template <>
NGX_INLINE DLL_HEADER int Ngx_Mesh :: GetElementIndex<0> (int nr) const
{
  return (*mesh).pointelements[nr].index;
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
NGX_INLINE DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<0> (int nr) const
{
  const Element0d & el = mesh->pointelements[nr];
  
  Ng_Element ret;
  ret.type = NG_PNT;
  ret.index = el.index;
  
  ret.points.num = 1;
  ret.points.ptr = (int*)&el.pnum;
  
  ret.vertices.num = 1;
  ret.vertices.ptr = (int*)&el.pnum;
  
  ret.edges.num = 0;
  ret.edges.ptr = NULL;
  
  ret.faces.num = 0;
  ret.faces.ptr = NULL;
  
  return ret;
}



template <> 
NGX_INLINE DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<1> (int nr) const
{
  const Segment & el = mesh->LineSegment (SegmentIndex(nr));

  Ng_Element ret;
  ret.type = NG_ELEMENT_TYPE(el.GetType());
  ret.index = el.si;
  ret.points.num = el.GetNP();
  ret.points.ptr = (int*)&(el[0]);

  ret.vertices.num = 2;
  ret.vertices.ptr = (int*)&(el[0]);

  ret.edges.num = 1;
  ret.edges.ptr = (T_EDGE2*)mesh->GetTopology().GetSegmentElementEdgesPtr (nr);

  ret.faces.num = 0;
  ret.faces.ptr = NULL;

  // ret.is_curved = mesh->GetCurvedElements().IsSegmentCurved(nr);
  ret.is_curved = el.IsCurved();

  return ret;
}

template <> 
NGX_INLINE DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<2> (int nr) const
{
  const Element2d & el = mesh->SurfaceElement (SurfaceElementIndex (nr));
  
  Ng_Element ret;
  ret.type = NG_ELEMENT_TYPE(el.GetType());
  ret.index = mesh->GetFaceDescriptor(el.GetIndex()).BCProperty();
  ret.points.num = el.GetNP();
  ret.points.ptr  = (int*)&el[0];

  ret.vertices.num = el.GetNV();
  ret.vertices.ptr = (int*)&(el[0]);

  ret.edges.num = MeshTopology::GetNEdges (el.GetType());
  ret.edges.ptr = (T_EDGE2*)mesh->GetTopology().GetSurfaceElementEdgesPtr (nr);

  ret.faces.num = MeshTopology::GetNFaces (el.GetType());
  ret.faces.ptr = (T_FACE2*)mesh->GetTopology().GetSurfaceElementFacesPtr (nr);

  ret.is_curved = el.IsCurved();
  return ret;
}

template <> 
NGX_INLINE DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<3> (int nr) const
{
  const Element & el = mesh->VolumeElement (ElementIndex (nr));
  
  Ng_Element ret;
  ret.type = NG_ELEMENT_TYPE(el.GetType());
  ret.index = el.GetIndex();
  ret.points.num = el.GetNP();
  ret.points.ptr = (int*)&el[0];

  ret.vertices.num = el.GetNV();
  ret.vertices.ptr = (int*)&(el[0]);

  ret.edges.num = MeshTopology::GetNEdges (el.GetType());
  ret.edges.ptr = (T_EDGE2*)mesh->GetTopology().GetElementEdgesPtr (nr);

  ret.faces.num = MeshTopology::GetNFaces (el.GetType());
  ret.faces.ptr = (T_FACE2*)mesh->GetTopology().GetElementFacesPtr (nr);

  // ret.is_curved = mesh->GetCurvedElements().IsElementCurved(nr);
  ret.is_curved = el.IsCurved();

  return ret;
}
