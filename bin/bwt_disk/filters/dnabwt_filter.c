#include "filters.h"


// maximum value for a run of N or ACGT: runs longer than this are splitted
// it must be close to 1<<21 in order not to waste space in the 7x8 encoding
// and it has to be a multiple of 4 to avoid overflow problems in 
// dna encoding (see flush_A())    
#define Run_limit ((1<<21)-4)

// array for decoding of ACGT values
static uint8 ACGTdecode[] = {'A','C','G','T'};


typedef struct {
  FILE *in;       // file from where the compressed data is read
  FILE *out;      // file where the compressed data is written
  int32 Arun;     // flag indicating whether we are in the middle of a 
                  // ACGT run or a N run 
  int32 Acount;   // number of pending ACGT chars (for r/w)
  int32 Ncount;   // number of pending N chars (for r/w)
  uint8 *Abuffer; // buffer where pending ACGT chars are stored 
  int32 Apos;     // current position inside Abuffer (for reading only)
  int64 tot_written; // number of bytes written (for writing only)
  int32 raw_bytes;// number of bytes to be sent/received raw 
  int32 raw_init; // initial value for the field raw_bytes 
} dna_file;


// =========== local prototypes ==============
static void flush_A(dna_file *f);
static void flush_N(dna_file *f);
static int32 fwrite7x8(uint32 n, FILE *f);
static int32 fread7x8(FILE *f);



/* **************************************************************
   open a file for reading or writing according to mode
   initialize the data structures required for range (de)coding
   _bwtext_ can be safely removed
   ************************************************************** */
static void *_bwtext_dnabwt_open(char *fname, int32 mode, int32 raw)
{ 
  int i;
  dna_file *f;

  if(sizeof(off_t)<8) 
    _bwtext_fatal_error("_bwtext_dnabwt_open: your system does't support >4GB files");

  f = (dna_file *) malloc(sizeof(dna_file));
  if(f==NULL) _bwtext_out_of_mem("_bwtext_dnabwt_open");
  f->Abuffer = (uint8 *) malloc(Run_limit);
  if(f->Abuffer==NULL) _bwtext_out_of_mem("_bwtext_dnabwt_open");
  for(i=0;i<Run_limit;i++) f->Abuffer[i] = 'A'; // init with legal values
  f->Apos=f->Acount=f->Ncount=0;
  f->raw_bytes = f->raw_init = raw;

  if(mode=='r') {
    f->in = fopen(fname,"rb");   // was: fopen64
    if(f->in==NULL) {
      perror(fname);
      _bwtext_ioerror(fname,"_bwtext_dnabwt_open");
    }
    f->out = NULL;
    f->Arun=0;  // a sequence of N was just completed
  }   
  else if(mode=='w') {
    // fprintf(stderr,"Opening for writing file %s\n",fname);
    f->in = NULL;
    f->out = fopen(fname,"wb"); // was: fopen64
    if(f->out==NULL) {
      perror(fname);
      _bwtext_ioerror(fname,"_bwtext_dnabwt_open");
    }
    f->Arun=1; // we are currently in a sequence of ACGT's
    f->tot_written=0;
  }
  else 
    _bwtext_fatal_error("_bwtext_dnabwt_open (illegal mode)");

  return f;
}


void *_bwtext_dnabwt_open0(char *fname, int32 mode)
{
  return _bwtext_dnabwt_open(fname, mode, 0);
}

void *_bwtext_dnabwt_open16(char *fname, int32 mode)
{
  return _bwtext_dnabwt_open(fname, mode, 16);
}

/* deallocate data structure and close file */ 
int64 _bwtext_dnabwt_close(void *w)
{
  dna_file *f = (dna_file *) w;
  int64 pos = 0; 

  if(f->out!=NULL) {
    if(f->Acount>0) {
      assert(f->Arun==1);
      flush_A(f);
    } else if(f->Ncount>0) {
      assert(f->Arun==0);
      flush_N(f);
    }
    pos = f->tot_written;
    if(fclose(f->out)!=0) return -1;
  }
  else {
    if(fclose(f->in)!=0) return -1;
  }
  free(f->Abuffer);
  free(f);
  return pos;
}


/* ***************************************************************
   read up to n symbols from *w and write them to b[]
   return number of symbols read
   *************************************************************** */   
int32 _bwtext_dnabwt_read(void *w, uint8 *b, int32 n) 
{
  dna_file *f = (dna_file *) w;
  int32 i,j,s,written=0;  // number of symbols written to *b

  // fprintf(stderr,"#### %ld %d\n",f->Apos,n);

  // check if we already reached the end of the file
  if(f->Acount<0 || f->Ncount<0) return 0;

  while(written<n && f->raw_bytes>0) {
    s = getc(f->in);
    if(s==EOF) _bwtext_fatal_error("_bwtext_dnabwt_open (unexpected EOF)");
    b[written++] = s;
    f->raw_bytes--;
  }

  while(written<n) {
    if(f->Ncount>0) { // there are N's pending 
      assert(f->Acount==0);
      assert(f->Arun==0);
      while(written<n && f->Ncount>0) {
	b[written++] = 'N';
	f->Ncount--; 
      }
    }
    else if(f->Acount>0) { // there are ACGT's pending
      assert(f->Ncount==0);
      assert(f->Arun==1);
      while(written<n && f->Acount>0) {
	b[written++] = f->Abuffer[f->Apos++];
	f->Acount--; 
      }
    }
    else if(f->Arun==1) { 
      // no symbols pending and just completed a run of ACGT's
      assert(f->Ncount==0);
      assert(f->Acount==0);
      f->Arun=0;         // a run of N's is starting 
      f->Ncount = fread7x8(f->in);
      // fprintf(stderr,"<-- N %d\n",f->Ncount);
      if(f->Ncount<0) 
	return written; // there are no other runs
    }
    else {  // no symbols pending and we just completed a run of N's
      assert(f->Ncount==0);
      assert(f->Acount==0);
      assert(f->Arun==0);
      f->Arun=1;               // a run of ACGT's is starting 
      f->Acount = fread7x8(f->in);
      // fprintf(stderr,"<-- A %d\n",f->Acount);
      if(f->Acount<0) 
	return written; // there are no other runs
      for(j=i=0;i<(f->Acount+3)/4;i++) {
	s = getc(f->in);
	if(s==EOF) _bwtext_fatal_error("_bwtext_dnabwt_open (unexpected EOF)");
	f->Abuffer[j++]=ACGTdecode[s&0x03];
	f->Abuffer[j++]=ACGTdecode[(s>>2)&0x03];
	f->Abuffer[j++]=ACGTdecode[(s>>4)&0x03];
	f->Abuffer[j++]=ACGTdecode[(s>>6)&0x03];
        // Note: we are writing to Abuffer up to 3 "garbage" symbols: it is
        // not a problem becaues Abuffer-size is a multiple of 4 and
        // Acount keeps track of the correct buffer size 
      }
      assert(j>=f->Acount);
      f->Apos = 0; 	     
    }
  }
  return written;
}


// encode b[0,n) and write the result to *w
int32 _bwtext_dnabwt_write(void *w, uint8 *b, int32 n)
{
  dna_file *f = (dna_file *) w;
  int32 written=0,next;

  // fprintf(stderr,"%%%% %llx %d\n",f->tot_written,n);

  while(written<n && f->raw_bytes>0) {
    next =b[written++];
    putc(next,f->out);
    f->tot_written++;
    f->raw_bytes--;
  }

  // fprintf(stderr,"%%%% %llx %d\n",f->tot_written,n);

  while(written<n) {
    next=b[written++];
    if(next!='N' && next!='A' && next!='C' && next!='G' && next!='T')
      _bwtext_fatal_error("_bwtext_dnabwt_write (invalid symbol)");
    if(f->Arun==0 && next=='N') { 
      if(++f->Ncount >= Run_limit) {
	assert(f->Ncount==Run_limit);
	flush_N(f);
      }
    }
    else if(f->Arun==1 && next=='N') {
      flush_A(f);
      f->Ncount++;
    }
    // next is A|C|G|T
    else if(f->Arun==0) { // end the current run of N's
      flush_N(f);
      f->Abuffer[f->Acount++] = (uint8) next;
    }
    else {
      f->Abuffer[f->Acount++] = (uint8) next;
      if(f->Acount>=Run_limit) {
	assert(f->Acount==Run_limit);
	flush_A(f);
      }
    }
  }
  return written;
}

void _bwtext_dnabwt_rewind(void *w)
{
  dna_file *f = (dna_file *) w;

  if(f->out!=NULL)
    _bwtext_fatal_error("_bwtext_dnabwt_rewind (illegal request)");
  assert(f->in!=NULL);
  rewind(f->in);
  f->Arun=0;
  f->Apos=f->Acount=f->Ncount=0;
  f->raw_bytes = f->raw_init;
}



void _bwtext_register_dnabwt_filter(filter *t)
{
  t->name = "dna five symbols encoding (ACGTN)";
  t->extension = ".atn";
  t->open = _bwtext_dnabwt_open16;
  t->read = _bwtext_dnabwt_read;
  t->write = _bwtext_dnabwt_write;
  t->rewind = _bwtext_dnabwt_rewind;
  t->close = _bwtext_dnabwt_close;
}


void _bwtext_register_dna_filter(filter *t)
{
  t->name = "pure dna five symbols encoding  (A, C, G, T, N)";
  t->extension = ".2bn";
  t->open = _bwtext_dnabwt_open0;
  t->read = _bwtext_dnabwt_read;
  t->write = _bwtext_dnabwt_write;
  t->rewind = _bwtext_dnabwt_rewind;
  t->close = _bwtext_dnabwt_close;
}



// ================= local functions ===============

// write to disk the current sequence of ACGT
static void flush_A(dna_file *f)
{
  int32 i,j,k,b,t=0; 

  assert(f->Arun==1);
  assert(f->Acount <= Run_limit);
  // fprintf(stderr,"(%llx) A --> %d\n",f->tot_written,f->Acount);
  f->tot_written += fwrite7x8(f->Acount, f->out);
  // write ceil(f->Acount/4) bytes: we write a number of symbols which 
  // is ==0 mod 4. This could be more than f->Acount but there is no overflow
  // since Run_limit is ==0 mod 4 
  for(j=i=0;i<(f->Acount+3)/4;i++) {
    b=0;
    for(k=0; k<=6; k+=2) { // do 4 iterations 
      switch(f->Abuffer[j++]) {
      case 'A': t=0; break;
      case 'C': t=1; break;
      case 'G': t=2; break;
      case 'T': t=3; break;
      default:
	_bwtext_fatal_error("flush_A (invalid symbol in buffer)");
      }
      b += t << k;
    }
    assert(b<256);
    fputc(b,f->out);
  }
  f->tot_written += (f->Acount+3)/4;
  f->Arun=0;
  f->Ncount=0;
  return;
}

static void flush_N(dna_file *f)
{

  assert(f->Arun==0);
  assert(f->Ncount <= Run_limit);
  // fprintf(stderr,"(%llx) N --> %d\n",f->tot_written,f->Ncount);
  f->tot_written += fwrite7x8(f->Ncount, f->out);
  f->Arun=1;
  f->Acount=0;
  return;
}




// write n using the 7x8 scheme. The most significant bit of each byte 
// indicates if the representation extends to the successive bytes.
// return the number of bytes written
static int32 fwrite7x8(uint32 n, FILE *f)
{
  uint32 t;
  int32 written=0;

  do {
    t = n & 0x7f;       // takes last seven bits
    n = n >> 7;
    if(n==0) t |= 0x80;  // set 8th bit
    putc(t,f);
    // fprintf(stderr,"%x ",t);
    written++;
  } while(n>0);
  return written;
}


static int32 fread7x8(FILE *f)
{
  int32 i,c,ris=0;

  c = getc(f);
  if(c==EOF) return -1;

  for(i=0;  ;i+=7) {
    assert(i<=28);
    ris |= (c& 0x7f) << i;     // add 7 bits 
    if((c&0x80)!=0) break;       // more bits ?
    c = getc(f);
    if(c==EOF) _bwtext_fatal_error("_bwtext_read7x8 (unexpected EOF)");
  }
  return ris;
}

