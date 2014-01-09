namespace ngfem
{


class DefaultHash
{
  int classnr;
  int order;
  int intorder;
public:
  DefaultHash (int aclassnr, int aorder, int aintorder)
    : classnr(aclassnr), order(aorder), intorder(aintorder) { ; }

  int ClassNr() const { return classnr; }
  int Order() const { return order; }
  int IntOrder() const { return intorder; }
};

inline int HashValue (const DefaultHash & hash, int size)
{
  return (hash.ClassNr() + 32 * (hash.Order()+hash.IntOrder()) ) % size;
}

inline bool operator== (const DefaultHash & hash1, const DefaultHash & hash2) 
{
  return (hash1.ClassNr() == hash2.ClassNr() && 
	  hash1.Order() == hash2.Order() && 
	  hash1.IntOrder() == hash2.IntOrder());
}





template <typename T, typename T_HASH = DefaultHash>
class PrecomputedShapesContainer
{
  HashTable<T_HASH, T*> shapes;
  
public:
  PrecomputedShapesContainer () 
    : shapes(640)        // 32 classes, roughly up to order+intorder 20
  { ; }

  void Add (int classnr, int order, int intorder, T * pre)
  {
    shapes.Set (T_HASH (classnr, order, intorder), pre);
  }
  
  T * Get (int classnr, int order, int intorder)
  {
    if (shapes.Used (T_HASH (classnr, order, intorder)))
      return shapes.Get(T_HASH (classnr, order, intorder));
    else
      return NULL;
  }
};


}
