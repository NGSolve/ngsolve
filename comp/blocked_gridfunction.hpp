#ifndef FILE_BLOCKEDGRIDFUNCTION
#define FILE_BLOCKEDGRIDFUNCTION

namespace ngcomp
{

  // class for cache optimized handling of gridfunctions (blocksize around 10)

#define HARDBLOCKSIZE 10


  template <class SCAL>
  class SmallBlockGridFunction : public S_GridFunction<SCAL>
  {
  protected:
    int softblocksize;
    int restvecblocksize;

    Array< S_BaseVector<SCAL> * > blockvec;


  public:
    SmallBlockGridFunction (const FESpace & afespace, const string & aname, const Flags & flags)
      : S_GridFunction<SCAL> (afespace,aname,flags)
    {
      softblocksize = flags.GetNumFlag("size",1);

      restvecblocksize = softblocksize % HARDBLOCKSIZE;

      if(restvecblocksize == 0)
	blockvec.SetSize(softblocksize/HARDBLOCKSIZE);
      else
	blockvec.SetSize(softblocksize/HARDBLOCKSIZE + 1);

      blockvec = NULL;
    

      Visualize (this->name);
    }


    ~SmallBlockGridFunction()
    {
      for(int i=0; i<blockvec.Size(); i++)
	delete blockvec[i];
    }

    
    virtual void Update ()
    {
      try
	{
	  NgLock lock(this -> mutex, 1);
	
	  int ndof = this->GetFESpace().GetNDof();
	
	  for(int i=0; i<blockvec.Size(); i++)
	    {
	      if(blockvec[i] && ndof = blockvec[i]->Size())
		break;

	      S_BaseVector<SCAL> * ovec = blockvec[i];

	      if(i == blockvec.Size()-1 && restvecblocksize != 0)
		{
		  switch(restvecblocksize)
		    {
		    case 1:
		      blockvec[i] = new VVector<SCAL>(ndof);
		      break;
		    case 2:
		      blockvec[i] = new VVector < Vec<2,SCAL> >(ndof);
		      break;
		    case 3:
		      blockvec[i] = new VVector < Vec<3,SCAL> >(ndof);
		      break;
		    case 4:
		      blockvec[i] = new VVector < Vec<4,SCAL> >(ndof);
		      break;
		    case 5:
		      blockvec[i] = new VVector < Vec<5,SCAL> >(ndof);
		      break;
		    case 6:
		      blockvec[i] = new VVector < Vec<6,SCAL> >(ndof);
		      break;
		    case 7:
		      blockvec[i] = new VVector < Vec<7,SCAL> >(ndof);
		      break;
		    case 8:
		      blockvec[i] = new VVector < Vec<8,SCAL> >(ndof);
		      break;
		    case 9:
		      blockvec[i] = new VVector < Vec<9,SCAL> >(ndof);
		      break;
		    case 10:
		      blockvec[i] = new VVector < Vec<10,SCAL> >(ndof);
		      break;
		    case 11:
		      blockvec[i] = new VVector < Vec<11,SCAL> >(ndof);
		      break;
		    case 12:
		      blockvec[i] = new VVector < Vec<12,SCAL> >(ndof);
		      break;
		    default:
		      throw Exception(string("SmallBlockGridfunction size ") + softblocksize + string(" not possible. Contact developers."));
		      break;				    
		    }
		}
	      else
		blockvec[i] = new VVector<HARDBLOCKSIZE,SCAL>(ndof);

	      if (this->nested && ovec && this->GetFESpace().GetProlongation())
		{
		  *testout << "do prolongation" << endl;

		  *blockvec[i] = SCAL(0);
		
		  *(blockvec[i]->Range (0, ovec->Size())) += (*ovec);
		
		
		  const_cast<ngmg::Prolongation&> (*this->GetFESpace().GetProlongation()).Update();
		
		  cout << "prolongate gridfunction" << endl;
		  this->GetFESpace().GetProlongation()->ProlongateInline
		    (this->GetMeshAccess().GetNLevels()-1, *blockvec[i]);
		}
	      else
		{
		  *vec[i] = SCAL(0);
		}

	      Visualize (this->name);
	    
	      delete ovec;
	    }
	

	  this -> level_updated = this -> ma.GetNLevels();
	}
      catch (exception & e)
	{
	  Exception e2 (e.what());
	  e2.Append ("\nIn SmallBlockGridFunction::Update()\n");
	  throw e2;
	}
      catch (Exception & e)
	{
	  e.Append ("In SmallBlockGridFunction::Update()\n");
	  throw e;
	}    
    }
	
	


    virtual BaseVector & GetVector (int comp = 0);
    virtual const BaseVector & GetVector (int comp = 0) const;

    ///
    virtual void GetElementVector (const Array<int> & dnums,
				   FlatVector<TSCAL> & elvec) const;

    ///
    virtual void SetElementVector (const Array<int> & dnums,
				   const FlatVector<TSCAL> & elvec);


    ///
    virtual void GetElementVector (int comp,
				   const Array<int> & dnums,
				   FlatVector<TSCAL> & elvec) const;

    ///
    virtual void SetElementVector (int comp,
				   const Array<int> & dnums,
				   const FlatVector<TSCAL> & elvec);



  
    


  };




}



#endif // FILE_BLOCKEDGRIDFUNCTION
