#ifndef FILE_MSGHANDLER
#define FILE_MSGHANDLER

/**************************************************************************/
/* File:   msghandler.hh                                                  */
/* Author: Johannes Gerstmayr                                             */
/* Date:   20. Nov. 99                                                    */
/**************************************************************************/


namespace netgen
{

  extern void PrintDot(char ch = '.');


  //Message Pipeline:

  //importance: importance of message: 1=very important, 3=middle, 5=low, 7=unimportant
  extern DLL_HEADER void PrintMessage(int importance, 
			   const MyStr& s1, const MyStr& s2=MyStr());
  extern DLL_HEADER void PrintMessage(int importance, 
			   const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4=MyStr());
  extern DLL_HEADER void PrintMessage(int importance, 
			   const MyStr& s1, const MyStr& s2, const MyStr& s3, const MyStr& s4, 
			   const MyStr& s5, const MyStr& s6=MyStr(), const MyStr& s7=MyStr(), const MyStr& s8=MyStr());

  // CR without line-feed
  extern DLL_HEADER void PrintMessageCR(int importance, 
			     const MyStr& s1, const MyStr& s2="", const MyStr& s3="", const MyStr& s4="", 
			     const MyStr& s5="", const MyStr& s6="", const MyStr& s7="", const MyStr& s8="");
  extern DLL_HEADER void PrintFnStart(const MyStr& s1, const MyStr& s2="", const MyStr& s3="", const MyStr& s4="", 
			   const MyStr& s5="", const MyStr& s6="", const MyStr& s7="", const MyStr& s8="");
  extern DLL_HEADER void PrintWarning(const MyStr& s1, const MyStr& s2="", const MyStr& s3="", const MyStr& s4="", 
			   const MyStr& s5="", const MyStr& s6="", const MyStr& s7="", const MyStr& s8="");
  extern DLL_HEADER void PrintError(const MyStr& s1, const MyStr& s2="", const MyStr& s3="", const MyStr& s4="", 
			 const MyStr& s5="", const MyStr& s6="", const MyStr& s7="", const MyStr& s8="");
  extern DLL_HEADER void PrintFileError(const MyStr& s1, const MyStr& s2="", const MyStr& s3="", const MyStr& s4="", 
			     const MyStr& s5="", const MyStr& s6="", const MyStr& s7="", const MyStr& s8="");
  extern DLL_HEADER void PrintSysError(const MyStr& s1, const MyStr& s2="", const MyStr& s3="", const MyStr& s4="", 
			    const MyStr& s5="", const MyStr& s6="", const MyStr& s7="", const MyStr& s8="");
  extern DLL_HEADER void PrintUserError(const MyStr& s1, const MyStr& s2="", const MyStr& s3="", const MyStr& s4="", 
			     const MyStr& s5="", const MyStr& s6="", const MyStr& s7="", const MyStr& s8="");
  extern DLL_HEADER void PrintTime(const MyStr& s1="", const MyStr& s2="", const MyStr& s3="", const MyStr& s4="", 
			const MyStr& s5="", const MyStr& s6="", const MyStr& s7="", const MyStr& s8="");
  extern DLL_HEADER void SetStatMsg(const MyStr& s);

  extern DLL_HEADER void PushStatus(const MyStr& s);
  extern DLL_HEADER void PushStatusF(const MyStr& s);
  extern DLL_HEADER void PopStatus();
  extern DLL_HEADER void SetThreadPercent(double percent);
  extern DLL_HEADER void GetStatus(MyStr & s, double & percentage);
}


#endif

