/* *************************************************************************/
/* File:   umfpackinverse.hpp                                              */
/* Author: Matthias Hochsteger                                             */
/* Date:   Nov. 15                                                         */
/* *************************************************************************/

#include <la.hpp>
#include "umfpackinverse.hpp"

#ifdef USE_UMFPACK

#ifdef USE_MKL
  extern "C" void mkl_free_buffers (void);
#endif // USE_MKL

using ngbla::integer;

#include <umfpack.h>

#include "umfpackinverse.hpp"

template<typename T> double GetReal( const T& x ) { return x.real(); }
template<> double GetReal( const double& x ) { return x; }

template<typename T> double GetImag( const T& x ) { return x.imag(); }
template<> double GetImag( const double& x ) { throw("cannot access imaginary part of double"); }

template<typename T> void SetReal( T& x, double val ) { x.real(val); }
template<> void SetReal( double& x, double val ) { x=val; }

template<typename T> void SetImag( T& x, double val ) { x.imag(val); }
template<> void SetImag( double& x, double val ) { }

namespace ngla
{

  class SubsetAll
  {
  public:
    bool Used(int i) { return true; }
    bool Used(int i, int j) { return true; }
  };

  class SubsetFree
  {
    const BitArray & free;
  public:
    SubsetFree (const BitArray & afree) : free(afree) { ; }
    bool Used(int i) { return free.Test(i); }
    bool Used(int i, int j) { return free.Test(i) && free.Test(j); }
  };

  class SubsetCluster
  {
    const Array<int> & cluster;
  public:
    SubsetCluster (const Array<int> & acluster) : cluster(acluster) { ; }
    bool Used(int i) { return cluster[i]; }
    bool Used(int i, int j) { return cluster[i] && (cluster[i] == cluster[j]); }
  };


  template<class TM>
  UmfpackInverseTM<TM> ::
  UmfpackInverseTM (shared_ptr<const SparseMatrixTM<TM>> a,
		    shared_ptr<BitArray> ainner,
		    shared_ptr<const Array<int>> acluster,
		    int asymmetric)
    : SparseFactorization (a, ainner, acluster)
  {
    static Timer timer("Umfpack Inverse");
    static Timer timer_sym("Umfpack Inverse symbolic");
    static Timer timer_num("Umfpack Inverse numeric");
    RegionTimer reg (timer);


    symmetric = asymmetric;
    //is_complex = mat_traits<TM>::IS_COMPLEX;
    is_complex = ngbla::IsComplex<TM>();

    (*testout) << "Umfpack, symmetric = " << symmetric << endl;

    if (inner && cluster)
      throw Exception("UmfpackInverse: Cannot use inner and cluster");

    if ( (inner && inner->Size() < a->Height()) ||
	 (cluster && cluster->Size() < a->Height() ) )
      {
	cout << "UmfpackInverse: Size of inner/cluster does not match matrix size!" << endl;
	throw Exception("Invalid parameters inner/cluster. Thrown by UmfpackInverse.");
      }

    if ( int( mat_traits<TM>::WIDTH) != int(mat_traits<TM>::HEIGHT) )
      {
	cout << "UmfpackInverse: Each entry in the square matrix has to be a square matrix!" << endl;
	throw Exception("No Square Matrix. Thrown by UmfpackInverse.");
      }


    entrysize = mat_traits<TM>::HEIGHT;
    height = a->Height() * entrysize;

    *testout << "matrix.InverseTpye = " <<  a->GetInverseType() << endl;

    if (inner)
      GetUmfpackMatrix (*a, SubsetFree (*inner));
    else if (cluster)
      GetUmfpackMatrix (*a, SubsetCluster (*cluster));
    else
      GetUmfpackMatrix (*a, SubsetAll());

    nze = rowstart[compressed_height];


    cout << IM(3) << "call umfpack ..." << flush;

    if (task_manager) task_manager -> StopWorkers();

    int status;
    double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
    umfpack_dl_defaults(Control);

    double *data = reinterpret_cast<double *>(values.Data());
    if (values.Size() > 0)
      {
        try
          {
            if(is_complex)
              {
                status = umfpack_zl_symbolic ( compressed_height, compressed_height, rowstart.Data(), indices.Data(), data, nullptr, &Symbolic, Control, Info );
                umfpack_zl_report_status( nullptr, status );
                if( status!= UMFPACK_OK ) throw Exception("UmfpackInverse: Symbolic factorization failed.");
                
                status = umfpack_zl_numeric (rowstart.Data(), indices.Data(), data, nullptr, Symbolic, &Numeric, nullptr, nullptr );
                umfpack_zl_report_status( nullptr, status );
                if( status!= UMFPACK_OK ) throw Exception("UmfpackInverse: Numeric factorization failed.");
              }
            else
              {
                timer_sym.Start();
                status = umfpack_dl_symbolic ( compressed_height, compressed_height, rowstart.Data(), indices.Data(), data, &Symbolic, Control, Info );
                timer_sym.Stop();
                umfpack_dl_report_status( nullptr, status );
                if( status!= UMFPACK_OK ) throw Exception("UmfpackInverse: Symbolic factorization failed.");

                timer_num.Start();                
                status = umfpack_dl_numeric (rowstart.Data(), indices.Data(), data, Symbolic, &Numeric, nullptr, nullptr );
                timer_num.Stop();                
                
                umfpack_dl_report_status( nullptr, status );
                if( status!= UMFPACK_OK ) throw Exception("UmfpackInverse: Numeric factorization failed.");
              }
          }
        catch(Exception & e)
          {
            if (task_manager) task_manager -> StartWorkers();
            throw;
          }
      }

    if (task_manager) task_manager -> StartWorkers();

    cout << IM(3) << " done" << endl;


  }


  template<class TM> template <typename TSUBSET>
  void UmfpackInverseTM<TM> ::
  GetUmfpackMatrix (const SparseMatrixTM<TM> & a, TSUBSET subset)
  {
    Array<int> icompress (a.Height());
    icompress = -1;

    compress.SetSize(0);

    int cnt = 0;
    for (int i = 0; i < a.Height(); i++)
      if (subset.Used(i))
	{
	  icompress[i] = cnt++;
	  compress.Append(i);
	}
    compressed = true;
    compressed_height = cnt * entrysize;

    rowstart.SetSize (compressed_height+1);
    rowstart = 0;
    if (symmetric)
      {
	// --- transform lower left to upper right triangular matrix ---
	// 1.) build array 'rowstart':
	// (a) get nr. of entries for each row

	for (int i = 0; i < a.Height(); i++)
	  {
	    if (!subset.Used(i)) continue;

	    FlatArray<int> ind = a.GetRowIndices(i);
	    for (int j = 0; j < ind.Size(); j++ )
	      {
		if (!subset.Used(i,ind[j])) continue;
		int ccol = entrysize*icompress[ind[j]];
		int rrow = entrysize*icompress[i];

		if (i != ind[j])
		  for (int k = 0; k < entrysize; k++ )
                    {
                      rowstart[ccol+k+1] += entrysize;
                      rowstart[rrow+k+1] += entrysize;
                    }
		else
		  for (int k = 0; k < entrysize; k++ )
		    rowstart[ccol+k+1] += entrysize; // -k;
	      }
	  }

	// (b) accumulate
	rowstart[0] = 0;
	for (int i = 1; i <= compressed_height; i++)
	  rowstart[i] += rowstart[i-1];


        // 2.) build whole (unsymmetrically stored) matrix:
	indices.SetSize(rowstart[compressed_height]);
	values.SetSize (rowstart[compressed_height]);

	Array<int> counter(compressed_height);
	counter = 0;
	for (int i = 0; i < a.Height(); i++ )
	  {
	    if (!subset.Used(i)) continue;
	    int ri = entrysize*icompress[i];

	    FlatArray<int> ind = a.GetRowIndices(i);
	    FlatVector<TM> values = a.GetRowValues(i);

	    for (int j = 0; j < ind.Size(); j++)
	      {
		if (!subset.Used (i,ind[j])) continue;
		int ci = entrysize*icompress[ind[j]];

		for (int k = 0; k < entrysize; k++)
		  for (int l = 0; l < entrysize; l++ )
                    if ( (i != ind[j]) || (k >= l) )
                      {
                        int row = ri+k;
                        int col = ci+l;
                        // set entry
                        size_t index = rowstart[row]+counter[row];
                        indices[index] = col;
                        this->values[index] = Access(values[j],k,l);
                        counter[row]++;

                        // set transposed entry
                        if (row != col)
                          {
                            size_t index = rowstart[col]+counter[col];
                            indices[index] = row;
                            this->values[index] = Access(values[j],k,l);
                            counter[col]++;
                          }
                      }

              }
	  }
      }
    else
      {
	int counter = 0;
	for (int i = 0; i < a.Height(); i++ )
	  {
	    if (!subset.Used(i)) continue;
	    FlatArray<int> ind = a.GetRowIndices(i);

	    int rowelems = 0;
	    for (int j = 0; j < ind.Size(); j++)
	      if (subset.Used(i,ind[j]))
		rowelems += entrysize;

	    int ci = entrysize * icompress[i];
	    for (int k = 0; k < entrysize; k++, counter += rowelems )
	      rowstart[ci+k] = counter;
	  }
	rowstart[compressed_height] = counter;

	indices.SetSize(counter);
	values.SetSize (counter);

	for (int i = 0; i < a.Height(); i++ )
	  {
	    if (!subset.Used(i)) continue;
	    FlatArray<int> ind = a.GetRowIndices(i);
	    FlatVector<TM> values = a.GetRowValues(i);

	    int ci = entrysize*icompress[i];

	    int counter = 0;
	    for (int j = 0; j < ind.Size(); j++ )
	      if (subset.Used(i,ind[j]))
		{
		  int ccol = entrysize*icompress[ind[j]];
		  for (int k=0; k<entrysize; k++ )
		    for (int l=0; l<entrysize; l++ )
		      {
			indices[rowstart[ci+k]+counter+l] = ccol+l;
			this->values[rowstart[ci+k]+counter+l] = Access(values[j],k,l);
		      }
		  counter+=entrysize;
		}
	  }
      }
  }




  template<class TM, class TV_ROW, class TV_COL>
  void UmfpackInverse<TM,TV_ROW,TV_COL> ::
  Mult (const BaseVector & x, BaseVector & y) const
  {
    static Timer timer(string("Umfpack Solve, mat = ") + typeid(TM).name() + ", vec = " + typeid(TV_ROW).name());
    RegionTimer reg (timer);

    FlatVector<TVX> fx = x.FV<TVX> ();
    FlatVector<TVX> fy = y.FV<TVX> ();

    int nrhs = fx.Size() / (height/entrysize);

    if(nrhs>1) throw Exception("UmfpackInverse: Multiple right-hand sides not supported.");

    // bool is_vector_complex = mat_traits<TVX>::IS_COMPLEX;
    bool is_vector_complex = ngbla::IsComplex<TVX>();
    if(is_complex && !is_vector_complex) throw Exception("UmfpackInverse: Cannot solve with complex matrix and real vector.");

    if (fx.Size() != fy.Size())
      {
	cout << "UmfpackInverse::Mult .. sizes don't match" << endl;
	cout << "type<TVX> = " << typeid(TVX).name() << endl;
	cout << "type<TM> = " << typeid(TM).name() << endl;
	cout << "fx.size = " << fx.Size() << endl;
	cout << "fy.size = " << fy.Size() << endl;
	cout << "size(x) = " << x.Size() << endl;
	cout << "size(y) = " << y.Size() << endl;
	cout << "height = " << height/entrysize << endl;
      }

    if (this->values.Size() == 0)
      {
        y = 0.0;
        return;
      }

    double *data = reinterpret_cast<double *>(this->values.Data());

    if (is_complex)
      {
        // complex matrix and vectors
        FlatVector<TSCAL> mx(height, (TSCAL*)fx.Data());
        FlatVector<TSCAL> my(height, (TSCAL*)fy.Data());

        Vector<TSCAL> hx(compressed_height);
        Vector<TSCAL> hy(compressed_height);

        double *data_x = reinterpret_cast<double *>(hx.Data());
        double *data_y = reinterpret_cast<double *>(hy.Data());

        for (int i : Range(compress.Size()) )
           hx(i) = mx(compress[i]);

        int status = umfpack_zl_solve ( UMFPACK_Aat, rowstart.Data(), indices.Data(), data, nullptr, data_y, nullptr, data_x, nullptr, this->Numeric, nullptr, nullptr );
        umfpack_zl_report_status( nullptr, status );
        if( status!= UMFPACK_OK ) throw Exception("UmfpackInverse: Solve failed.");

        my = 0;
        for (int i : Range(compress.Size()) )
           my(compress[i]) = hy(i);

      }
    else
      {
        // real matrix
        FlatVector<double> mx(height, (double*)fx.Data());
        FlatVector<double> my(height, (double*)fy.Data());
        Vector<double> hx(compressed_height);
        Vector<double> hy(compressed_height);

        for (int i : Range(compress.Size()) )
          for (int j = 0; j < entrysize; j++)
            hx(i*entrysize+j) = GetReal(mx(compress[i]*entrysize+j));

        int status = umfpack_dl_solve ( UMFPACK_Aat, rowstart.Data(), indices.Data(), data, hy.Data(), hx.Data(), this->Numeric, nullptr, nullptr );
        umfpack_dl_report_status( nullptr, status );
        if( status!= UMFPACK_OK ) throw Exception("UmfpackInverse: Solve failed.");

        my = 0;
        for (int i : Range(compress.Size()) )
          for (int j = 0; j < entrysize; j++)
            my(compress[i]*entrysize+j) = hy(i*entrysize+j);

        if(is_vector_complex)
          {
            // complex vectors
            for (int i : Range(compress.Size()) )
              for (int j = 0; j < entrysize; j++)
                hx(i*entrysize+j) = GetImag(mx(compress[i]*entrysize+j));

            int status = umfpack_dl_solve ( UMFPACK_Aat, rowstart.Data(), indices.Data(), data, hy.Data(), hx.Data(), this->Numeric, nullptr, nullptr );
            umfpack_dl_report_status( nullptr, status );
            if( status!= UMFPACK_OK ) throw Exception("UmfpackInverse: Solve failed.");

            for (int i : Range(compress.Size()) )
              for (int j = 0; j < entrysize; j++)
                SetImag(my(i*entrysize+j), hy(compress[i]*entrysize+j));
          }
      }
  }




  template<class TM, class TV_ROW, class TV_COL>
  void UmfpackInverse<TM,TV_ROW,TV_COL> ::
  MultTrans (const BaseVector & x, BaseVector & y) const
  {
    static Timer timer(string("Umfpack Solve, mat = ") + typeid(TM).name() + ", vec = " + typeid(TV_ROW).name());
    RegionTimer reg (timer);

    FlatVector<TVX> fx = x.FV<TVX> ();
    FlatVector<TVX> fy = y.FV<TVX> ();

    int nrhs = fx.Size() / (height/entrysize);

    if(nrhs>1) throw Exception("UmfpackInverse: Multiple right-hand sides not supported.");

    // bool is_vector_complex = mat_traits<TVX>::IS_COMPLEX;
    bool is_vector_complex = ngbla::IsComplex<TVX>();
    if(is_complex && !is_vector_complex) throw Exception("UmfpackInverse: Cannot solve with complex matrix and real vector.");

    if (fx.Size() != fy.Size())
      {
	cout << "UmfpackInverse::Mult .. sizes don't match" << endl;
	cout << "type<TVX> = " << typeid(TVX).name() << endl;
	cout << "type<TM> = " << typeid(TM).name() << endl;
	cout << "fx.size = " << fx.Size() << endl;
	cout << "fy.size = " << fy.Size() << endl;
	cout << "size(x) = " << x.Size() << endl;
	cout << "size(y) = " << y.Size() << endl;
	cout << "height = " << height/entrysize << endl;
      }


    double *data = reinterpret_cast<double *>(this->values.Data());

    if (is_complex)
      {
        // complex matrix and vectors
        FlatVector<TSCAL> mx(height, (TSCAL*)fx.Data());
        FlatVector<TSCAL> my(height, (TSCAL*)fy.Data());

        Vector<TSCAL> hx(compressed_height);
        Vector<TSCAL> hy(compressed_height);

        double *data_x = reinterpret_cast<double *>(&hx[0]);
        double *data_y = reinterpret_cast<double *>(&hy[0]);

        for (int i : Range(compress.Size()) )
           hx(i) = mx(compress[i]);

        int status = umfpack_zl_solve ( UMFPACK_A, &rowstart[0], &indices[0], data, nullptr, data_y, nullptr, data_x, nullptr, this->Numeric, nullptr, nullptr );
        umfpack_zl_report_status( nullptr, status );
        if( status!= UMFPACK_OK ) throw Exception("UmfpackInverse: Solve failed.");

        my = 0;
        for (int i : Range(compress.Size()) )
           my(compress[i]) = hy(i);

      }
    else
      {
        // real matrix
        FlatVector<double> mx(height, (double*)fx.Data());
        FlatVector<double> my(height, (double*)fy.Data());
        Vector<double> hx(compressed_height);
        Vector<double> hy(compressed_height);

        for (int i : Range(compress.Size()) )
          for (int j = 0; j < entrysize; j++)
            hx(i*entrysize+j) = GetReal(mx(compress[i]*entrysize+j));

        int status = umfpack_dl_solve ( UMFPACK_Aat, &rowstart[0], &indices[0], data, &hy(0), &hx(0), this->Numeric, nullptr, nullptr );
        umfpack_dl_report_status( nullptr, status );
        if( status!= UMFPACK_OK ) throw Exception("UmfpackInverse: Solve failed.");

        my = 0;
        for (int i : Range(compress.Size()) )
          for (int j = 0; j < entrysize; j++)
            my(compress[i]*entrysize+j) = hy(i*entrysize+j);

        if(is_vector_complex)
          {
            // complex vectors
            for (int i : Range(compress.Size()) )
              for (int j = 0; j < entrysize; j++)
                hx(i*entrysize+j) = GetImag(mx(compress[i]*entrysize+j));

            int status = umfpack_dl_solve ( UMFPACK_A, &rowstart[0], &indices[0], data, &hy(0), &hx(0), this->Numeric, nullptr, nullptr );
            umfpack_dl_report_status( nullptr, status );
            if( status!= UMFPACK_OK ) throw Exception("UmfpackInverse: Solve failed.");

            for (int i : Range(compress.Size()) )
              for (int j = 0; j < entrysize; j++)
                SetImag(my(i*entrysize+j), hy(compress[i]*entrysize+j));
          }
      }
  }








  
  template<class TM>
  ostream & UmfpackInverseTM<TM> :: Print (ostream & ost) const
  {
    cout << "UmfpackInverse::Print not implemented!" << endl;
    return ost;
  }


  template<class TM>
  UmfpackInverseTM<TM> :: ~UmfpackInverseTM()
  {
    if(is_complex)
      {
        umfpack_zl_free_symbolic ( &Symbolic );
        umfpack_zl_free_numeric ( &Numeric );
      }
    else
      {
        umfpack_dl_free_symbolic ( &Symbolic );
        umfpack_dl_free_numeric ( &Numeric );
      }

#ifdef USE_MKL
    mkl_free_buffers();
#endif // USE_MKL
  }



  template<class TM>
  void UmfpackInverseTM<TM> :: SetMatrixType() // TM entry)
  {
    // if (mat_traits<TM>::IS_COMPLEX)
    if (ngbla::IsComplex<TM>())
      is_complex = true;
    else
      is_complex = false;

    if (print)
      cout << ", sym = " << int(symmetric)
	   << ", complex = " << is_complex << endl;
  }

  template<class TM>
  void UmfpackInverseTM<TM> :: Update()
  {
    cout << IM(3) << "call umfpack update..." << flush;

    auto castmatrix = dynamic_pointer_cast<const SparseMatrix<TM>>(matrix.lock());

    if (inner)
      GetUmfpackMatrix (*castmatrix, SubsetFree (*inner));
    else if (cluster)
      GetUmfpackMatrix (*castmatrix, SubsetCluster (*cluster));
    else
      GetUmfpackMatrix (*castmatrix, SubsetAll());
    
    if (task_manager) task_manager -> StopWorkers();

    int status;
    double *data = reinterpret_cast<double *>(values.Data());
    try
      {
        if(is_complex)
          {
            umfpack_zl_free_numeric ( &Numeric );

            status = umfpack_zl_numeric (&rowstart[0], &indices[0], data, nullptr, Symbolic, &Numeric, nullptr, nullptr );
            umfpack_zl_report_status( nullptr, status );
            if( status!= UMFPACK_OK ) throw Exception("UmfpackInverse: Numeric factorization failed.");
          }
        else
          {
            umfpack_dl_free_numeric ( &Numeric );

            status = umfpack_dl_numeric (&rowstart[0], &indices[0], data, Symbolic, &Numeric, nullptr, nullptr );
            umfpack_dl_report_status( nullptr, status );
            if( status!= UMFPACK_OK ) throw Exception("UmfpackInverse: Numeric factorization failed.");
          }
      }
    catch(Exception & e)
      {
        if (task_manager) task_manager -> StartWorkers();
        throw;
      }

    if (task_manager) task_manager -> StartWorkers();

    cout << IM(3) << " done" << endl;
  }




  template class UmfpackInverseTM<double>;
  template class UmfpackInverseTM<Complex>;
  template class UmfpackInverse<double>;
  template class UmfpackInverse<Complex>;
  template class UmfpackInverse<double,Complex,Complex>;
#if MAX_SYS_DIM >= 1
  template class UmfpackInverseTM<Mat<1,1,double> >;
  template class UmfpackInverseTM<Mat<1,1,Complex> >;
  template class UmfpackInverse<Mat<1,1,double> >;
  template class UmfpackInverse<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class UmfpackInverseTM<Mat<2,2,double> >;
  template class UmfpackInverseTM<Mat<2,2,Complex> >;
  template class UmfpackInverse<Mat<2,2,double> >;
  template class UmfpackInverse<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class UmfpackInverseTM<Mat<3,3,double> >;
  template class UmfpackInverseTM<Mat<3,3,Complex> >;
  template class UmfpackInverse<Mat<3,3,double> >;
  template class UmfpackInverse<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class UmfpackInverseTM<Mat<4,4,double> >;
  template class UmfpackInverseTM<Mat<4,4,Complex> >;
  template class UmfpackInverse<Mat<4,4,double> >;
  template class UmfpackInverse<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class UmfpackInverseTM<Mat<5,5,double> >;
  template class UmfpackInverseTM<Mat<5,5,Complex> >;
  template class UmfpackInverse<Mat<5,5,double> >;
  template class UmfpackInverse<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class UmfpackInverseTM<Mat<6,6,double> >;
  template class UmfpackInverseTM<Mat<6,6,Complex> >;
  template class UmfpackInverse<Mat<6,6,double> >;
  template class UmfpackInverse<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class UmfpackInverseTM<Mat<7,7,double> >;
  template class UmfpackInverseTM<Mat<7,7,Complex> >;
  template class UmfpackInverse<Mat<7,7,double> >;
  template class UmfpackInverse<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class UmfpackInverseTM<Mat<8,8,double> >;
  template class UmfpackInverseTM<Mat<8,8,Complex> >;
  template class UmfpackInverse<Mat<8,8,double> >;
  template class UmfpackInverse<Mat<8,8,Complex> >;
#endif


}



#endif

