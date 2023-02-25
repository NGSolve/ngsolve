/*********************************************************************/
/* File:   basematrix.cpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

/* 
   base class in matrix hierarchy
*/

#include <la.hpp>
#include <../parallel/parallel_matrices.hpp>

namespace ngla
{
  BaseMatrix :: BaseMatrix()
    : paralleldofs (NULL)
  {
    ;
  }

  BaseMatrix :: BaseMatrix (shared_ptr<ParallelDofs> aparalleldofs)
    : paralleldofs ( aparalleldofs )
  {     
    ;
  }
  
  BaseMatrix :: ~BaseMatrix ()
  {
    ;
  }

  int BaseMatrix :: VHeight() const
  {
    throw Exception (string("BaseMatrix::VHeight not overloaded, type = ")+typeid(*this).name());
  }
  
  int BaseMatrix :: VWidth() const
  {
    throw Exception (string("BaseMatrix::VWidth not overloaded, type = ")+typeid(*this).name());
  }

  BaseVector & BaseMatrix :: AsVector()
  {
    throw Exception (string("BaseMatrix::AsVector not overloaded, type = ")+typeid(*this).name());
  }

  const BaseVector & BaseMatrix :: AsVector() const
  {
    throw Exception (string("BaseMatrix::AsVector const not overloaded, type = ")+typeid(*this).name());    
  }
  
  void BaseMatrix :: SetZero()
  {
    AsVector() = 0;
  }

  ostream & BaseMatrix :: Print (ostream & ost) const
  {
    return (ost << "Print base-matrix" << endl);
  }

  Array<MemoryUsage> BaseMatrix :: GetMemoryUsage () const { return Array<MemoryUsage>(); } 

  size_t BaseMatrix :: NZE () const
  {
    throw Exception(string("NZE not overloaded for matrix-type")
                    +typeid(*this).name());
  }
  
  shared_ptr<BaseMatrix> BaseMatrix :: CreateMatrix () const
  {
    throw Exception (string("BaseMatrix::CreateMatrix not overloaded, type = ")+typeid(*this).name());        
  }

  AutoVector BaseMatrix :: CreateVector () const
  {
    throw Exception (string("BaseMatrix::CreateVector not overloaded, type = ")+typeid(*this).name());            
  }

  // AutoVector BaseMatrix :: CreateRowVector () const
  // {
  //   return CreateVector();
  // }

  // AutoVector BaseMatrix :: CreateColVector () const
  // {
  //   return CreateVector();
  // }




  
  void BaseMatrix :: Mult (const BaseVector & x, BaseVector & y) const
  {
    // y = 0;
    if(safety_check & 1)
      throw Exception("Mult or MultAdd must be implemented for BaseMatrix!");
    y.SetZero();
    MultAdd (1, x, y);
  }

  void BaseMatrix :: MultTrans (const BaseVector & x, BaseVector & y) const
  {
    if(IsSymmetric().IsTrue())
      { Mult(x, y); return; }
    if(safety_check & 2)
      throw Exception("MultTransAdd or MultTrans must be implemented for (maybe) not symmetric BaseMatrix!");

    // y = 0;
    y.SetZero();
    MultTransAdd (1, x, y);
  }

  void BaseMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    //    cout << "Warning: BaseMatrix::MultAdd(double), this = " << typeid(*this).name() << endl;
    auto temp = y.CreateVector();
    safety_check |= 1;
    Mult (x, *temp);
    y += s * *temp;
  }

  void BaseMatrix :: MultAdd (Complex s, const BaseVector & x, BaseVector & y) const 
  {
    /*
    stringstream err;
    err << "BaseMatrix::MultAdd (Complex) called, type = " 
	<< typeid(*this).name();
    throw Exception (err.str());
    */
    auto temp = y.CreateVector();
    safety_check |= 1;
    Mult (x, *temp);
    y += s * *temp;
  }
  
  void BaseMatrix :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    if(IsSymmetric().IsTrue())
      { MultAdd(s, x, y); return; }

    auto temp = y.CreateVector();
    safety_check |= 2;
    MultTrans (x, *temp);
    y += s * *temp;
    /*
    cout << "warning: BaseMatrix::MultTransAdd(double) calls MultAdd, ";
    cout << "type = " << typeid(*this).name() << endl;
    MultAdd (s, x, y);
    return;
    */
    /*
    stringstream err;
    err << "BaseMatrix::MultTransAdd (double) called, type = " 
	<< typeid(*this).name();
    throw Exception (err.str());
    */
  }

  void BaseMatrix :: MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const
  {
    if(IsSymmetric().IsTrue())
      { MultAdd(s, x, y); return; }

    safety_check |= 2;
    auto temp = y.CreateVector();
    MultTrans (x, *temp);
    y += s * *temp;

    /*
    //    cout << "warning: BaseMatrix::MultTransAdd(complex) calls MultAdd" << endl;
    MultAdd (s, x, y);
    return;

    stringstream err;
    err << "BaseMatrix::MultTransAdd (Complex) called, type = " 
	<< typeid(*this).name();
    throw Exception (err.str());
    */
  }

  void BaseMatrix :: MultConjTransAdd (Complex s, const BaseVector & x, BaseVector & y) const
  {
    auto tmpx = x.CreateVector();
    auto tmpy = y.CreateVector();
    tmpx.FV<Complex>() = Conj(x.FV<Complex>());
    tmpy.FV<Complex>() = Conj(y.FV<Complex>());
    MultTransAdd (Conj(s), tmpx, tmpy);
    y.FV<Complex>() = Conj(tmpy.FV<Complex>());
    // throw Exception(string("MultHermitianAdd not overloaded for type ")+typeid(*this).name());
  }

  void BaseMatrix :: MultAdd (FlatVector<double> alpha, const MultiVector & x, MultiVector & y) const
  {
    for (int i = 0; i < alpha.Size(); i++)
      MultAdd (alpha[i], *x[i], *y[i]);
  }


  
   // to split mat x vec for symmetric matrices
  void BaseMatrix :: MultAdd1 (double s, const BaseVector & x, BaseVector & y,
			       const BitArray * ainner,
			       const Array<int> * acluster) const
  {
    MultAdd (s, x, y);
  }

  void BaseMatrix :: MultAdd2 (double s, const BaseVector & x, BaseVector & y,
			       const BitArray * ainner,
			       const Array<int> * acluster) const
  {
    ;
  }
  



  shared_ptr<BaseMatrix> BaseMatrix :: InverseMatrix (shared_ptr<BitArray> subset) const 
  {
    cerr << "BaseMatrix::InverseMatrix not available" << endl;
    return NULL;
  }
  
  shared_ptr<BaseMatrix> BaseMatrix :: InverseMatrix (shared_ptr<const Array<int>> clusters) const
  {
    cerr << "BaseMatrix::InverseMatrix not available" << endl;
    return NULL;
  }
  
  INVERSETYPE BaseMatrix :: SetInverseType ( INVERSETYPE ainversetype ) const
  {
    cerr << "BaseMatrix::SetInverseType not available" << endl;
    return SPARSECHOLESKY;
  }
  
  INVERSETYPE BaseMatrix :: SetInverseType ( string ainversetype ) const
  {
    cerr << "BaseMatrix::SetInverseType not available" << endl;
    return SPARSECHOLESKY;
  }
  
  INVERSETYPE BaseMatrix :: GetInverseType () const
  {
    cerr << "BaseMatrix::GetInverseType not available" << endl;
    return SPARSECHOLESKY;
  }

  void BaseMatrix :: DoArchive (Archive & ar)
  {
    ;
  }

  template <typename TSCAL>
  Matrix<TSCAL> BaseMatrix :: ToDense() const
  {
    auto vecx = CreateRowVector();
    auto vecy = CreateColVector();
    
    Matrix<TSCAL> dmat(Height(), Width());
    auto fx = vecx.FV<TSCAL>();
    auto fy = vecy.FV<TSCAL>();
    for (int i = 0; i < fx.Size(); i++)
      {
        fx = 0;
        fx(i) = 1;
        Mult (vecx, vecy);
        dmat.Col(i) = fy;
      }
    return std::move(dmat);
  }

  template Matrix<double> BaseMatrix :: ToDense<double>() const;
  template Matrix<Complex> BaseMatrix :: ToDense<Complex>() const;


  double BaseMatrix::Timing (int runs) const
  {
    Timer t("timing");
    auto vx = CreateRowVector();
    auto vy = CreateColVector();

    vx = 0;
    t.Start();
    for (int i = 0; i < runs; i++)
      Mult (vx, vy);
    t.Stop();
    return t.GetTime()/runs;
  }
  

  template<>
  S_BaseMatrix<double> :: S_BaseMatrix () 
  { ; }


  template<>
  S_BaseMatrix<double> :: ~S_BaseMatrix () 
  { ; }

  S_BaseMatrix<Complex> :: S_BaseMatrix () 
  { ; }

  S_BaseMatrix<Complex> :: ~S_BaseMatrix () 
  { ; }


  void S_BaseMatrix<Complex> :: 
  MultAdd (double s, const BaseVector & x, BaseVector & y) const 
  {
    MultAdd (Complex(s), x, y);
  }

  void S_BaseMatrix<Complex> :: 
  MultAdd (Complex s, const BaseVector & x, BaseVector & y) const 
  {
    stringstream err;
    err << "S_BaseMatrix<Complex>::MultAdd (Complex) called, type = " 
	<< typeid(*this).name();
    throw Exception (err.str());
  }
  
  void S_BaseMatrix<Complex> :: 
  MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    MultTransAdd (Complex(s), x, y);
  }

  void S_BaseMatrix<Complex> :: 
  MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const
  {
    stringstream err;
    err << "S_BaseMatrix<Complex>::MultTransAdd (Complex) called, type = " 
	<< typeid(*this).name();
    throw Exception (err.str());
  }


  void VMatVecExpr :: CheckSize (BaseVector & dest_vec) const
  {
    if (m.Height() != dest_vec.Size() || m.Width() != x.Size())
      throw Exception (ToString ("matrix-vector: size does not fit\n") +
                       "matrix-type = " + typeid(m).name() + "\n" +
                       "matrix:     " + ToString(m.Height()) + " x " + ToString(m.Width()) + "\n"
                       "vector in : " + ToString(x.Size()) + "\n"
                       "vector res: " + ToString(dest_vec.Size()));

  }



  string GetInverseName (INVERSETYPE type)
  {
    switch (type)
      {
      case PARDISO:         return "pardiso";
      case PARDISOSPD:      return "pardisospd";
      case SPARSECHOLESKY:  return "sparsecholesky";
      case SUPERLU:         return "superlu";
      case SUPERLU_DIST:    return "superlu_dist";
      case MUMPS:           return "mumps";
      case MASTERINVERSE:   return "masterinverse";
      case UMFPACK:         return "umfpack";
      }
    return "";
  }

  std::map<type_index, function<shared_ptr<BaseMatrix>(const BaseMatrix&)>> BaseMatrix::devmatcreator;

  shared_ptr<BaseMatrix> BaseMatrix :: CreateDeviceMatrix() const
  {
    auto it = devmatcreator.find(typeid(*this));
    if (it == devmatcreator.end())
      {
        cout << IM(7) << "No deviceMatrix creator function for type " << typeid(*this).name() << endl;
        return const_cast<BaseMatrix*>(this)->shared_from_this();
        // return nullptr;
      }
    cout << IM(7) << "DeviceMatrix creator function for type " << typeid(*this).name() << endl;
    return (*it).second(*this);
  }


  
  BaseMatrix::OperatorInfo BaseMatrix :: GetOperatorInfo () const
  {
    OperatorInfo info;
    info.name = typeid(*this).name();
    info.height = Height();
    info.width = Width();
    return info;
  }

  BaseMatrix::OperatorInfo IdentityMatrix :: GetOperatorInfo () const
  {
    OperatorInfo info;
    if (has_format)
      {
        info.name = "Identity";
        info.height = Height();
        info.width = Width();
      }
    else
      {
        info.name = "Identity (any format)";
        info.height = 0;
        info.width = 0;
      }
    return info;
  }

  BaseMatrix::OperatorInfo Transpose :: GetOperatorInfo () const
  {
    OperatorInfo info;
    info.name = "Transpose";
    try
      {
        info.height = Height();
        info.width = Width();
      }
    catch (Exception &)
      {
        cerr << "Transpose::GetOperatorInfo, got exception for H/W" << endl;        
      };
    info.childs += &bm;
    return info;
  }

  
  BaseMatrix::OperatorInfo SumMatrix :: GetOperatorInfo () const
  {
    OperatorInfo info;
    info.name = "SumMatrix";
    try
      {
        info.height = Height();
        info.width = Width();
      }
    catch (Exception &)
      {
        cerr << "SumMatrix::GetOperatorInfo, got exception for H/W" << endl;        
      };
    info.childs += &bma;
    info.childs += &bmb;
    return info;
  }

  
  
  BaseMatrix::OperatorInfo ProductMatrix :: GetOperatorInfo () const
  {
    OperatorInfo info;
    info.name = "ProductMatrix";
    try
      {
        info.height = Height();
        info.width = Width();
      }
    catch (Exception &)
      {
        cerr << "ProductMatrix::GetOperatorInfo, got exception for H/W" << endl;
      }
    info.childs += &bma;
    info.childs += &bmb;
    return info;
  }


  
  void BaseMatrix :: PrintOperatorInfo (ostream & ost, int level) const
  {
    auto info = GetOperatorInfo();
    
    ost << string(2*level, ' ');
    ost << info.name << ", h = " << info.height << ", w = " << info.width;
    if (IsComplex()) ost << " complex";
    ost << endl;
    for (auto c : info.childs)
      {
        try
          {
            c->PrintOperatorInfo (ost, level+1);
          }
        catch (Exception & e)
          {
            ost << "got exception, child type is " << typeid(*c).name() << endl;
          } 
      }
  }





  shared_ptr<BaseMatrix> ComposeOperators (shared_ptr<BaseMatrix> a,
                                           shared_ptr<BaseMatrix> b)
  {
    if (auto embb = dynamic_pointer_cast<EmbeddingTranspose> (b))
      {
        // cout << "embeddingT optimization" << endl;
        return make_shared<EmbeddedTransposeMatrix> (embb->Width(), embb->GetRange(), a);
      }
    
    if (auto emba = dynamic_pointer_cast<Embedding> (a))
      {
        // cout << "embedding optimization" << endl;        
        return make_shared<EmbeddedMatrix> (emba->Height(), emba->GetRange(), b);
      }

    auto para = dynamic_pointer_cast<ParallelMatrix> (a);
    auto parb = dynamic_pointer_cast<ParallelMatrix> (b);

    if (para && parb)
      {
        if (RowType(para->GetOpType()) == ColType(parb->GetOpType()))
          {
            // cout << "combining parallel matrices" << endl;
            auto localprod = ComposeOperators (para->GetMatrix(), parb->GetMatrix());
            return make_shared<ParallelMatrix> (localprod,
                                                parb->GetRowParallelDofs(),
                                                para->GetColParallelDofs(),
                                                ParallelOp(RowType(parb->GetOpType()), ColType(para->GetOpType())));
          }
        else
          {
            cerr << "illegal operator composition" << endl;
            cerr << "optyp A = " << int(para->GetOpType()) << endl;
            cerr << "optyp B = " << int(parb->GetOpType()) << endl;
            auto & locmata = *para->GetMatrix();
            auto & locmatb = *parb->GetMatrix();
            cerr << "type a parallelmat of = " << typeid(locmata).name()
                 << ", type b = " <<typeid(locmatb).name() << endl;
          }
      }
    
    return make_shared<ProductMatrix> (a, b);
  }

  shared_ptr<BaseMatrix> AddOperators (shared_ptr<BaseMatrix> a,
                                       shared_ptr<BaseMatrix> b,
                                       double faca, double facb)
  {
    auto para = dynamic_pointer_cast<ParallelMatrix> (a);
    auto parb = dynamic_pointer_cast<ParallelMatrix> (b);
    if (para && parb)
      {
        if (para->GetOpType() == parb->GetOpType())
          return make_shared<ParallelMatrix> (AddOperators (para->GetMatrix(), parb->GetMatrix(), faca, facb), 
                                              para->GetRowParallelDofs(),
                                              para->GetColParallelDofs(),
                                              para->GetOpType());

        cerr << "Adding parallel matrices of different types, type a = "
             << int(para->GetOpType()) << ", type b = " << int(parb->GetOpType()) << endl;
      }
    return make_shared<SumMatrix> (a, b, faca, facb);
  }

  shared_ptr<BaseMatrix> TransposeOperator (shared_ptr<BaseMatrix> mat)
  {
    if (auto emb = dynamic_pointer_cast<Embedding> (mat))
      return make_shared<EmbeddingTranspose> (emb->Height(), emb->GetRange(), emb->IsComplex());

    if (auto emb = dynamic_pointer_cast<EmbeddingTranspose> (mat))
      return make_shared<Embedding> (emb->Width(), emb->GetRange(), emb->IsComplex());

    if (auto parmat = dynamic_pointer_cast<ParallelMatrix> (mat))
      return make_shared<ParallelMatrix> (TransposeOperator(parmat->GetMatrix()),
                                          parmat->GetColParallelDofs(),
                                          parmat->GetRowParallelDofs(),
                                          ParallelOp(InvertType(ColType(parmat->GetOpType())),
                                                     InvertType(RowType(parmat->GetOpType()))));

    
    return make_shared<Transpose> (mat);
  }
  
    

}
