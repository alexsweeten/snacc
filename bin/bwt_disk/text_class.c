/* *************************************************************
   functions handling reading/writing data from a (huge) text file
   in a sequential manner. For text file we mean one for which we 
   want to compute the bwt.
   All functions in the bwtext library use these function to access 
   the text file. The functions here access to the actual file
   using filters (eg plain, gzip). Thus, the bwtext library can use
   differen formats in a transparent way. The functions here provide a 
   buffering mechanism so no further buffering should be 
   necessary (eg the filters will only be requested to read/write
   in large chunks). 
   ************************************************************* */
#include "bwtext.h"


// external variables giving access to the text file
extern char *_bwtext_Tname;
extern filter *_bwtext_Text_filter;


#define TEXT_BUF_SIZE _bwtext_TEXT_BUFFER_SIZE


// local global variable representing the current open text file
static textfile *t=NULL;

// open text for reading
void _bwtext_text_open() 
{

  assert(t==NULL);
  t = (textfile *) malloc(sizeof(*t));
  if(!t) _bwtext_out_of_mem("_bwtext_text_open");

  t->filename = _bwtext_Tname;                // set filename to default
  t->fil = _bwtext_Text_filter;               // set filter to default filter
  t->filep = t->fil->open(t->filename,'r');   // open file for reading 
  if(t->filep==NULL) {
    perror(t->filename);
    _bwtext_ioerror(t->filename," _bwtext_text_open");
  }

  t->textbuffer = (uint8 *) malloc(TEXT_BUF_SIZE*sizeof(*t->textbuffer));
  if(!t->textbuffer) 
    _bwtext_out_of_mem("_bwtext_text_open");
  t->text_pos=t->buf_cur=t->buf_size=0;
}

void _bwtext_text_close()
{

  t->fil->close(t->filep);
  free(t->textbuffer);
  free(t);
  t=NULL;        // this is a static variable 
}

void _bwtext_text_rewind()
{
  t->fil->rewind(t->filep);
  t->text_pos=t->buf_cur=t->buf_size=0;
}

// return next char o EOF is no more available
int32 _bwtext_text_getc()
{

  if(t->buf_size>0 && t->buf_cur<t->buf_size) {
    // buffer not empty 
    t->text_pos++;
    return t->textbuffer[t->buf_cur++];
  }

  // refill buffer
  t->buf_size = t->fil->read(t->filep, t->textbuffer,TEXT_BUF_SIZE);
  t->buf_cur=0;
 
  if(t->buf_size==0)
    return EOF;    // no more char

  t->text_pos++;
  return t->textbuffer[t->buf_cur++];
}

boolean _bwtext_text_moretocome()
{
 
  if(t->buf_size>0 && t->buf_cur<t->buf_size) {
    // buffer not empty 
    return 1;
  }

  // refill buffer
  t->buf_size = t->fil->read(t->filep, t->textbuffer,TEXT_BUF_SIZE);
  t->buf_cur=0;
 
  if(t->buf_size==0)
    return 0;    // no more char
  else 
    return 1;
}

// copy up to size uint8's in buf starting from buf[size-1] and going
// backwards. return the number of uint8's copied 
uint32 _bwtext_text_copy_reverse(uint8 *buf, uint32 size)
{
  uint32 to_be_written = size;
  uint32 written=0;

  while(to_be_written>0) {
    if(t->buf_cur<t->buf_size) {
      buf[--to_be_written] = t->textbuffer[t->buf_cur++];
      written++;
    }
    else {
      // buffer empty: try to get new data 
      t->buf_size = t->fil->read(t->filep, t->textbuffer,TEXT_BUF_SIZE);
      t->buf_cur=0;
      if(t->buf_size==0) {
	t->text_pos += written;
	return written;
      }
    }
  }
  assert(to_be_written==0);
  assert(written==size);
  t->text_pos += written;
  return written;
}

boolean _bwtext_text_checkpos(uint64 p)
{
  return (boolean) (p==t->text_pos);
}


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   functions providing write access to a text file 
   These are used for example by the unbwt program 
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */

// open text for writing
textfile *_bwtext_text_open_write(char *s, filter *tf, int32 n) 
{
  textfile *t;

  t = (textfile *) malloc(sizeof(*t));
  if(!t) _bwtext_out_of_mem("_bwtext_text_open_write");

  t->buf_size = n;
  t->filename = s;
  t->fil = tf;
  t->filep = t->fil->open(s,'w');
  if(t->filep==NULL) {
    perror(s);
    _bwtext_ioerror(s,"_bwtext_text_open_write");
  }

  assert(t->buf_size>0);
  t->textbuffer = (uint8 *) malloc(t->buf_size*sizeof(*t->textbuffer));
  if(!t->textbuffer) 
    _bwtext_out_of_mem("_bwtext_text_open_write");
  t->text_pos=t->buf_cur=0;
  return t;
}

// write a single char
void _bwtext_text_putc(textfile *t, int32 c)
{

  assert(c<256 && c>=0);

  if(t->buf_cur==t->buf_size) {
    if(t->fil->write(t->filep,t->textbuffer,t->buf_size)!=t->buf_size) {
      perror(t->filename);
      _bwtext_ioerror(t->filename,"_bwtext_text_putc");
    }
    t->buf_cur=0;
  }
  t->text_pos++;
  t->textbuffer[t->buf_cur++] = (uint8) c;

}

// flush and close write buffer
void _bwtext_text_close_write(textfile *t)
{

  if(t->buf_cur>0) {
    if(t->fil->write(t->filep,t->textbuffer,t->buf_cur)!=t->buf_cur) {
      perror(t->filename);
      _bwtext_ioerror(t->filename,"_bwtext_text_close_write");
    }
    t->buf_cur=0;
  }
  // close and deallocate
  if(t->fil->close(t->filep)<0)
    _bwtext_ioerror(t->filename,"_bwtext_text_close_write");
  free(t->textbuffer);
  free(t);
}
