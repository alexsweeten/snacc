#include "filters.h"
#include <unistd.h>
#include <string.h>

extern int _bwtext_Verbose;                // verbosity level
extern double _bwtext_FilterPar1;          // parameter #1
extern int _bwtext_FilterPar1_flag;        // !=0 if parameter #1 was defined
extern double _bwtext_FilterPar2;          // parameter #2
extern int _bwtext_FilterPar2_flag;        // !=0 if parameter #2 was defined


typedef struct {
  FILE *stream;
  int32 type;           // 0 for reading, 1 for writing
  char *command;        // command used for opening the pipe
  char *filename;       // useful for rewinding 
} piped_file;


/* **************************************************************
   open a file for reading or writing according to mode
   initialize the data structure containg pipe information 
   ************************************************************** */
void *_bwtext_pipedlzma_open(char *fname, int32 mode)
{ 
  piped_file *f;

  if(sizeof(off_t)<8) 
    _bwtext_fatal_error("_bwtext_pipedlzma_open: " 
			"your system does't support >4GB files");

  if(mode=='r') {
    f = (piped_file *) malloc(sizeof(piped_file));
    if(f==NULL) _bwtext_out_of_mem("_bwtext_pipedlzma_open");
    f->filename = strdup(fname);
    if(f->filename==NULL) _bwtext_out_of_mem("_bwtext_pipedlzma_open");
    asprintf(&(f->command),"lzcat -S .bwt %s",fname);
    if(f->command==NULL) _bwtext_out_of_mem("_bwtext_pipedlzma_open");
    // fprintf(stderr,"***** Opening (R): %s\n",f->command); 
    f->stream = popen(f->command,"r");  
    if (f->stream==NULL)
    _bwtext_fatal_error("_bwtext_pipedlzma_open: unable to open pipe");
    f->type = 0;
    return f;
  }
   
  if(mode=='w') {
    f = (piped_file *) malloc(sizeof(piped_file));
    if(f==NULL) _bwtext_out_of_mem("_bwtext_pipedlzma_open");
    f->filename = strdup(fname);
    if(f->filename==NULL) _bwtext_out_of_mem("_bwtext_pipedlzma_open");
    asprintf(&(f->command), "lzma -c > %s",fname);
    if(f->command==NULL) _bwtext_out_of_mem("_bwtext_pipedlzma_open");
    // fprintf(stderr,"***** Opening (W): %s\n",f->command); 
    f->stream = popen(f->command,"w");  
    if (f->stream==NULL)
    _bwtext_fatal_error("_bwtext_pipedlzma_open: unable to open pipe");
    f->type = 1;
    return f;
  }
  _bwtext_fatal_error("_bwtext_pipedlzma_open (illegal mode)");
  return NULL;
}


int64 _bwtext_pipedlzma_close(void *w)
{

  piped_file *f = (piped_file *) w;

  if(pclose(f->stream)!=0) return -1;
  free(f->command);
  free(f->filename);
  free(f);

  return 0;//pos;
}


int32 _bwtext_pipedlzma_read(void *w, uint8 *b, int32 n)
{
  piped_file *f = (piped_file *) w;
  return fread(b,1,n,f->stream);
}


int32 _bwtext_pipedlzma_write(void *w, uint8 *b, int32 n)
{
  piped_file *f = (piped_file *) w;
  return write(fileno((FILE *)f->stream),b,n);
}


void _bwtext_pipedlzma_rewind(void *w)
{
  piped_file *f = (piped_file *) w;

  // fprintf(stderr,"***** Closing: %s\n",f->command);// !!!!!!!!!
  int ris =  pclose(f->stream);
  if(ris!=0) {
    ;
    //error(0,ris,"");
    // fprintf(stderr,"[pclose %d] ",ris);
    // _bwtext_fatal_error("_bwtext_pipedlzma_rewind: error closing pipe");
  }

  f->stream = popen(f->command,"r");  
  if (f->stream==NULL)
    _bwtext_fatal_error("_bwtext_pipedlzma_rewind: unable to open pipe");
}


void _bwtext_register_pipedlzma_filter(filter *t)
{
  t->name = "piped lzma";
  t->extension = ".lzma";
  t->open = _bwtext_pipedlzma_open;
  t->read = _bwtext_pipedlzma_read;
  t->write = _bwtext_pipedlzma_write;
  t->rewind = _bwtext_pipedlzma_rewind;
  t->close = _bwtext_pipedlzma_close;
}

