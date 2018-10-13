/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
   global.c

   global variables and error handling routines for the bwtext package

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.

   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include "bwtext.h"


// ========== global variables ====================================
filter *_bwtext_Text_filter, *_bwtext_Bwt_filter; // filters 
char *_bwtext_Tname, *_bwtext_Bwtname;            // file names
int _bwtext_Verbose=0;                            // verbosity level
int _bwtext_Use_divsufsort=1;      // use divsufsort for suffix sorting
/* ================================================================
   parameters that can be passed to the filter (usually from the 
   command line). these parameters apply only to the filter used for 
   writing. it is assumed that such parameters are hard-coded in the 
   output file so no parameters need to be specified for reading
   ================================================================ */  
double _bwtext_FilterPar1;          // parameter #1
int _bwtext_FilterPar1_flag;        // !=0 if parameter #1 was defined
double _bwtext_FilterPar2;          // parameter #2
int _bwtext_FilterPar2_flag;        // !=0 if parameter #2 was defined



void _bwtext_fatal_error(char *s)
{
  fprintf(stderr,"Fatal error in function: %s\n",s);
  exit(1);
}

void _bwtext_ioerror(char *file, char *fun)
{
  fprintf(stderr,"Error working on file %s in function: %s\n",file,fun);
  exit(_bwtext_IOERROR);
}


void _bwtext_out_of_mem(char *f)
{
  fprintf(stderr, "Out of memory in function %s!\n", f);
  exit(2);
}



// ========= only useful for debigging 

#if 0
// ----- this function prints any char in a readable form
void _bwtext_pretty_putchar(int c)
{
  
  if(c>=32 && c<127)      // printable char
    fprintf(stderr,"  %c", c);
  else if(c=='\n')
    fprintf(stderr," \\n");        // \n
  else if(c=='\t')
    fprintf(stderr," \\t");        // \t
  else     
    fprintf(stderr," %02x", c);      // print hex code
}


void _bwtext_print_block(uint8 *t, uint64 n, char *s)
{
  uint64 i;

  fprintf(stderr,"------- begin: %s ---------\n",s);
  for(i=0;i<n;i++)
    _bwtext_pretty_putchar(t[i]);
  fprintf(stderr,"------- end: %s ---------\n",s);

}


void _bwtext_print_string(uint8 *t, int64 n) 
{
  int k;
  for(k=0;k<n;k++)
    _bwtext_pretty_putchar(t[k]);
  fprintf(stderr,"\n");
}

#endif
