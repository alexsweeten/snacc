#include "filters.h"
#include <unistd.h>


// open file fname for reading or writing  according to mode
// return NULL on failure
void *_bwtext_plain_open(char *fname, int32 mode)
{ 
  if(mode=='r')
    return fopen(fname,"rb");       
  if(mode=='w')
    return fopen(fname,"wb");       
  _bwtext_fatal_error("illegal mode (_bwtext_plain_open)");
  return NULL;
}

// close file. If possible return the number of bytes read/written 
// so far (it is used only for providing info in Verbose mode)
// otherwise return 0. Return -1 if an error occurs
int64 _bwtext_plain_close(void *f)
{
  int64 pos = ftello((FILE *)f); 
 
  if(fclose((FILE *)f)!=0) return -1;
  return pos;
}

// read up to n bytes from f and store them in buffer b[]
// return number of bytes actually read 
int32 _bwtext_plain_read(void *f, uint8 *b, int32 n)
{
  return fread(b,1,n,f);  // alternative: return read(fileno((FILE *)f), b, n);
}

// write to f n bytes from buffer b[]
// return number of bytes actually written  
int32 _bwtext_plain_write(void *f, uint8 *b, int32 n)
{
  return write(fileno((FILE *)f),b,n);
}

// rewind file f (set pointer at the beginning of the file)
// this operation has to be supported only for files opened in read mode
void _bwtext_plain_rewind(void *f)
{
  rewind((FILE *)f);
}

// stores the functions in the filter data structure 
void _bwtext_register_plain_filter(filter *t)
{
  t->name = "no compression";
  t->extension = "";
  t->open = _bwtext_plain_open;
  t->read = _bwtext_plain_read;
  t->write = _bwtext_plain_write;
  t->rewind = _bwtext_plain_rewind;
  t->close = _bwtext_plain_close;
}

