#include "filters.h"
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <zlib.h>


void *_bwtext_gz_open(char *fname, int32 mode)
{
  int mask=0;
  char *smode="";

  if(sizeof(off_t)<8) 
    _bwtext_fatal_error("_bwtext_gz_open (>4GB files not supported)");

  // does a very rudimental conversion
  if (mode=='r')
    {mask = O_RDONLY; smode="r";}
  else if(mode=='w')
    {mask = O_WRONLY | O_CREAT; smode="w";}
  else
    _bwtext_fatal_error("_bwtext_gz_open (illegal mode)");

  int fd = open(fname,mask, S_IRUSR|S_IWUSR);   // was: open64
  if(fd==-1) return NULL;
  return gzdopen(fd,smode); 
}


// close the file f return the current position in file
// if available
int64 _bwtext_gz_close(void *f)
{
    if(gzclose((gzFile *)f)!=0) return -1;
    return 0; // position not available sorry!
}

int32 _bwtext_gz_read(void *f, uint8 *b, int32 n)
{
  return gzread((gzFile *)f, b, n);
}


int32 _bwtext_gz_write(void *f, uint8 *b, int32 n)
{
  //assert(n>0);
  if(n==0) return 0;
  int ris=gzwrite((gzFile *)f,b,n);
  return ris==0 ? -1 : ris;   // returns -1 on error
}

void _bwtext_gz_rewind(void *f)
{
  gzrewind((gzFile *)f);
}


void _bwtext_register_gz_filter(filter *t)
{
  t->name = "gzip (via libz)";
  t->extension = ".gz";
  t->open = _bwtext_gz_open;
  t->read = _bwtext_gz_read;
  t->write = _bwtext_gz_write;
  t->rewind = _bwtext_gz_rewind;
  t->close = _bwtext_gz_close;
}

