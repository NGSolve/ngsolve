#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

#define maxlen 1000

int main (int argc, char ** argv)
{
  if (argc != 3)
    {
      cout << "use:  makerlsfile  infile outfile" << endl;
      exit(1);
    }
  

  char line[maxlen], infile[maxlen], outfile[maxlen];\
  char ch;
  int i, j;

  /*
  cout << "infile: ";
  cin >> infile;
  cout << "outfile: ";
  cin >> outfile;
  */

  ifstream inf (argv[1]);
  ofstream outf (argv[2]);

  outf << "const char * ngscript[] = {" << endl;
  while (inf.good())
    {
      i = 0;

      inf.get(ch);
      while (ch != '\n' && inf.good() && i < maxlen)
	{
	  if (ch == '\"')
	    {
	      line[i] = '\\';
	      line[i+1] = '\"';
	      i+= 2;
	    }
	  else if (ch == '\\')
	    {
	      line[i] = '\\';
	      line[i+1] = '\\';
	      i+= 2;
	    }
	  else
	    {
	      line[i] = ch;
	      i++;
	    }
	  inf.get(ch);
	}
      line[i] = 0;
      cout << line << endl;
      outf << "\"" << line << "\\n\",\\" << endl;
    }
  outf << "0};" << endl;
  return 0;
}
