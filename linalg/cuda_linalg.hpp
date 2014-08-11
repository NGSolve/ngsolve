#ifdef CUDA

namespace ngla
{

  class UnifiedVector : public BaseVector
  {
    int size;
    double * host_data;
    double * dev_data;
    bool host_uptodate;
    bool dev_uptodate;
    
  public:
    UnifiedVector (int asize);
    
    BaseVector & operator= (double d);
    BaseVector & operator= (BaseVector & v2);

    void UpdateHost ();    
    void UpdateDevice ();
  };

}

#endif
