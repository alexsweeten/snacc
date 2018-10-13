/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
   misc.c

   miscellaneous routines for the bwtext package

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.

   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include "bwtext.h"




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



// ========= useful for debigging 

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
