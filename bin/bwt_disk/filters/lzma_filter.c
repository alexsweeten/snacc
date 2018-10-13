#include "filters.h"
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include "lzmadec.h"


typedef struct {
  int32 type;           // 0 for reading, 1 for writing
  // ---- field for reading operations  
  lzmadec_FILE *input;
  // ----- fields for writing operations
  FILE *stream;
  char *command;   // command used for opening the pipe (not necessary here)
} lzma_file;



/* **************************************************************
   open a file for reading or writing according to mode
   Reading is done usin the libz-like functions of liblzmadec,
   writing uses pipes (ii is very slow so it should not be used)
   ************************************************************** */
void *_bwtext_lzma_open(char *fname, int32 mode)
{
  lzma_file *f;

  if(sizeof(off_t)<8) 
    _bwtext_fatal_error("_bwtext_lzma_open (>4GB files not supported)");

  f = (lzma_file *) malloc(sizeof(lzma_file));
  if(f==NULL) _bwtext_out_of_mem("_bwtext_lzma_open");


  // distinguish between reading and writing
  if (mode=='r') {
    f->type = 0;
    int fd = open(fname,O_RDONLY, S_IRUSR|S_IWUSR); 
    if(fd==-1) return NULL;
    f->input = lzmadec_dopen(fd);

  } else if(mode=='w') {
    f->type = 1;
    asprintf(&(f->command), "lzma -1 -c > %s",fname);
    if(f->command==NULL) _bwtext_out_of_mem("_bwtext_lzma_open");
    // fprintf(stderr,"***** Opening (W): %s\n",f->command); 
    f->stream = popen(f->command,"w");  
    if (f->stream==NULL)
      _bwtext_fatal_error("_bwtext_lzma_open: unable to open pipe");

  } else
    _bwtext_fatal_error("_bwtext_lzma_open (illegal mode)");

  return f; 
}


// close the file f return the current position in file
// if available
int64 _bwtext_lzma_close(void *w)
{
  lzma_file *f = (lzma_file *) w;

  if(f->type==0) {
    if(lzmadec_close(f->input)!=0) return -1;
  }
  else if(f->type==1) {
    if(pclose(f->stream)!=0) return -1;
    free(f->command);
  }
  else
    _bwtext_fatal_error("_bwtext_lzma_close (illegal mode)");

  free(f);
  return 0; // position not available sorry!
}

int32 _bwtext_lzma_read(void *w, uint8 *b, int32 n)
{
  lzma_file *f = (lzma_file *) w;

  return lzmadec_read(f->input, b, n);
}


int32 _bwtext_lzma_write(void *w, uint8 *b, int32 n)
{
  if(n==0) return 0;

  lzma_file *f = (lzma_file *) w;
  return write(fileno((FILE *)f->stream),b,n);
}

void _bwtext_lzma_rewind(void *w)
{
  lzma_file *f = (lzma_file *) w;
  if(f->type!=0)
    _bwtext_fatal_error("_bwtext_lzma_rewind (llegal request)");
  lzmadec_rewind(f->input);
}


void _bwtext_register_lzma_filter(filter *t)
{
  t->name = "lzma (via liblzamdec)";
  t->extension = ".lzma";
  t->open = _bwtext_lzma_open;
  t->read = _bwtext_lzma_read;
  t->write = _bwtext_lzma_write;
  t->rewind = _bwtext_lzma_rewind;
  t->close = _bwtext_lzma_close;
}

