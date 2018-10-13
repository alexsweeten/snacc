#include "filters.h"
#include "range.c"

#define RLE_symbol (_BW_ALPHA_SIZE)
#define EOF_symbol (_BW_ALPHA_SIZE+1)
#define Illegal_symbol (-1)
#define default_Max_frequency 16383     // maximum allowed frequency count
#define default_Increment 1             // freq. increment for each occurrence


typedef struct {
  FILE *file;
  int32 type;           // 0 for reading, 1 for writing
  rc_encoder *rce;      // encoder 
  rc_decoder *rcd;      // decode
  rc_model   *rcm;      // model
  Uint *freqT;          // frequency table used by the model 
  uint64 pending;       // number of pending symbols  
  int32  pending_symb;  // current pending symbol
  int32  next_symb;     // symbol after the pending one     
} rc_file;


// =========== local prototypes ==============
static uint32 uint32_read(FILE *f);
static void uint32_write(FILE *f, uint32 u);
static int32 int_log2(uint64 u);
static void rcrle_write_pending_run(rc_file *f);
static void init_rcfilter_parameters(int32 *maxf, int32 *incr);
static void rc_init_decoding(rc_file *f);
static void rc_init_encoding(rc_file *f);
static void rc_stop_decoding(rc_file *f);
static void rc_stop_encoding(rc_file *f);

/* **************************************************************
   open a file for reading or writing according to mode
   initialize the data structures required for range (de)coding
   ************************************************************** */
void *_bwtext_rc_open(char *fname, int32 mode)
{ 
  rc_file *f;

  if(sizeof(off_t)<8) 
    _bwtext_fatal_error("_bwtext_rc_open: your system does't support >4GB files");

  if(mode=='r') {
    f = (rc_file *) malloc(sizeof(rc_file));
    if(f==NULL) _bwtext_out_of_mem("_bwtext_rc_open");
    f->file = fopen(fname,"rb");   // was: fopen64
    f->type = 0;
    rc_init_decoding(f);
    return f;
  }   
  if(mode=='w') {
    f = (rc_file *) malloc(sizeof(rc_file));
    if(f==NULL) _bwtext_out_of_mem("_bwtext_rc_open");
    f->file = fopen(fname,"wb"); // was: fopen64
    f->type = 1;
    rc_init_encoding(f);
    return f;
  }
  _bwtext_fatal_error("_bwtext_rc_open (illegal mode)");
  return NULL;
}

/* deallocate rc-related data structures and close file */ 
int64 _bwtext_rc_close(void *w)
{
  rc_file *f = (rc_file *) w;
  FILE *tmpfile = f->file;


  if(f->type==0)
    rc_stop_decoding(f);
  else if(f->type==1)
    rc_stop_encoding(f);
  else
    _bwtext_fatal_error("_bwtext_rc_close (illegal mode)");
  free(f);

  // current position (??? it may be changed below due to buffering and eof)
  int64 pos = ftello(tmpfile);  // was: (int64) ftello64(tmpfile); 

  if(fclose(tmpfile)!=0) return -1;
  return pos;
}

/* ***************************************************************
   read up to n symbols from w and write them to b[]
   return number of symbols read
   *************************************************************** */   
int32 _bwtext_rcrle_read(void *w, uint8 *b, int32 n)
{
  rc_file *f = (rc_file *) w;
  int32 next,read=0;
  
  if(f->pending_symb==EOF_symbol) return 0;

  while(1) {
    // ---- take care of pending symbols 
    while(f->pending>0) {
      if(read==n) return n;  // done 
      assert(f->pending_symb>=0 && f->pending_symb<_BW_ALPHA_SIZE);
      b[read++] = (uint8) f->pending_symb;
      f->pending--;
    }
    if(read==n) return n;  // done
    // ---- start a new run
    assert(f->pending==0);
    if(f->next_symb==Illegal_symbol)
      f->pending_symb = DecodeSymbol(f->rcd,f->rcm);
    else {
      f->pending_symb = f->next_symb;
      f->next_symb = Illegal_symbol;
    }
    f->pending=1;
    if(f->pending_symb==EOF_symbol) {
      assert(read<n);
      return read;         // sorry we do not have n symbols
    }
    assert(f->pending_symb!=RLE_symbol);
    // --- compute length of new run
    next=DecodeSymbol(f->rcd, f->rcm);
    while(next==f->pending_symb || next==RLE_symbol) {
      if(next==RLE_symbol)
	f->pending = (f->pending<<1);        // RLE stands for bit 0 
      else
	f->pending = (f->pending<<1)+1;      // current symbol stands for bit 1
      next = DecodeSymbol(f->rcd, f->rcm);
    }
    assert(f->pending>=1);
    f->next_symb = next; 
  }
  // it is an error if we exit the main loop
  _bwtext_fatal_error(" _bwtext_rcrle_read");
  return 0;
}


// encode b[0,n) using rle+rc 
int32 _bwtext_rcrle_write(void *w, uint8 *b, int32 n)
{
  rc_file *f = (rc_file *) w;
  int32 written=0,next;

  while(written<n) {
    next=b[written++];
    if(next==f->pending_symb) {
      f->pending++; continue;
    }
    else { 
      if(f->pending_symb!=Illegal_symbol)
	rcrle_write_pending_run(f);
      f->pending_symb = next;
      f->pending = 1;
    }
  }
  return written;
}

void _bwtext_rc_rewind(void *w)
{
  rc_file *f = (rc_file *) w;

  if(f->type!=0)
    _bwtext_fatal_error("_bwtext_rc_rewind (illegal request)");
  rc_stop_decoding(f);
  rewind(f->file);
  rc_init_decoding(f);
}

void _bwtext_register_rcrle_filter(filter *t)
{
  t->name = "rle + range coding";
  t->extension = ".rrc";
  t->open = _bwtext_rc_open;
  t->read = _bwtext_rcrle_read;
  t->write = _bwtext_rcrle_write;
  t->rewind = _bwtext_rc_rewind;
  t->close = _bwtext_rc_close;
}



// ===================== static functions ===================

// create and initialize the model and the encoder
static void rc_init_decoding(rc_file *f)
{
  int32 increm, max_freq, r;
  uint32 t;

  // make sure f has been created for reading
  assert(f->type==0); 
  // -------- create and init model         
  t = uint32_read(f->file);    // read and check model parameters
  increm = t & 0xFFFF;
  max_freq = 1+((t>>16) & 0xFFFF);
  if(max_freq>65536 || increm>65535)
    _bwtext_fatal_error("rc_init_decoding (invalid update parameters)\n");
  // allocate freq table for alphabet symbols and EOF
  f->freqT = (Uint *) malloc(sizeof(Uint)*(_BW_ALPHA_SIZE+2));
  // allocate model
  f->rcm = (rc_model *) malloc(sizeof(rc_model));
  if(f->freqT==NULL || f->rcm==NULL) _bwtext_out_of_mem("rc_init_decoding");
  // init model
  r=ModelInit(f->rcm, _BW_ALPHA_SIZE+2,f->freqT,NULL,increm,max_freq,ADAPT);
  if (r != RC_OK)
    _bwtext_fatal_error("rc_init_decoding (ModelInit)\n");

  // ---- create and init decoder 
  f->rcd = (rc_decoder *) malloc(sizeof(rc_decoder));
  if(f->rcd==NULL) _bwtext_out_of_mem("rc_init_decoding");
  StartDecode(f->rcd, f->file, NULL, 0);

  f->rce=NULL;            // encoder not used   
  f->pending = 0;         // no pending symbol 
  f->pending_symb = Illegal_symbol;   
  f->next_symb    = Illegal_symbol;
}


// create and initialize the model and the decoder
static void rc_init_encoding(rc_file *f)
{
  int32 max_freq, increm, r;

  // make sure f has been created for writing
  assert(f->type==1); 
  // ------- write update parameters to file
  init_rcfilter_parameters(&max_freq, &increm);
  if(max_freq>65536 || increm>65535)     // check parameters
       _bwtext_fatal_error("rc_init_encoding (invalid parameters)\n");
  uint32_write(f->file,((max_freq-1)<<16)+increm);

  // allocate freq table for alphabet symbols and EOF
  f->freqT = (Uint *) malloc(sizeof(Uint)*(_BW_ALPHA_SIZE+2));
  // allocate model
  f->rcm = (rc_model *) malloc(sizeof(rc_model));
  if(f->freqT==NULL || f->rcm==NULL) _bwtext_out_of_mem("rc_init_encoding");
  // init model
  r=ModelInit(f->rcm, _BW_ALPHA_SIZE+2,f->freqT,NULL,increm,max_freq,ADAPT);
  if (r != RC_OK)
    _bwtext_fatal_error("rc_init_encoding (ModelInit)");

  // ---- create and init encoder 
  f->rce = (rc_encoder *) malloc(sizeof(rc_encoder));
  if(f->rce==NULL) _bwtext_out_of_mem("rc_init_encoding");
  StartEncode(f->rce, f->file, NULL, 0);

  f->rcd=NULL;            // decoder not used   
  f->pending = 0;         // no pending symbol 
  f->pending_symb = Illegal_symbol;   
  f->next_symb    = Illegal_symbol;    
}

static void rc_stop_decoding(rc_file *f)
{
  // assert(f->pending==0); // all symbols have been read ???
  assert(f->type==0);
  assert(f->rce==NULL);
  // deallocate structures allocated in rc_init_decoding
  free(f->rcd);
  free(f->rcm);
  free(f->freqT); 
}

static void rc_stop_encoding(rc_file *f)
{
  assert(f->type==1);
  assert(f->rcd==NULL);
  // flush remaining bits and add EOF
  rcrle_write_pending_run(f);
  EncodeSymbol(f->rce, f->rcm, EOF_symbol);  
  FinishEncode(f->rce);
  // deallocate structures allocated in rc_init_encoding
  free(f->rce);
  free(f->rcm);
  free(f->freqT);  
}

static void rcrle_write_pending_run(rc_file *f)
{
  uint64 count = f->pending;
  uint64 mask;
  int32 b, bits, symbol;

  bits=int_log2(count);   	      //  number of bit to write count
  mask = ( ((uint64) 1) << (bits-1) );
  for (b=bits-1; b>=0; b--) {         // write count in binary ("bits" bit)
    if((count & mask) == 0 )          // get b-th bit of count 
      symbol = RLE_symbol;            // RLE_symbol stands for 0
    else
      symbol = f->pending_symb;       // cur symbol stands for 1
    EncodeSymbol(f->rce,f->rcm,symbol); 
    if(f->rce->error!=RC_OK)
	_bwtext_fatal_error("rcrle_write_pending_run (Encoding Error)");
    count <<= 1;                  // avoid dangerous right shift of mask
  }
}

static void init_rcfilter_parameters(int32 *max_frequency, int32 *increment)
{
  // init order0 model parameters
  *max_frequency = _bwtext_FilterPar1_flag ? 
                 (int32) _bwtext_FilterPar1 : default_Max_frequency;
  *increment = _bwtext_FilterPar2_flag ? 
                 (int32) _bwtext_FilterPar2 : default_Increment;
  if(_bwtext_Verbose>2) 
    fprintf(stderr,"Range_coding options: Max_frequency=%d, Incr.=%d\n",
	    *max_frequency, *increment); 

}

/* Read an uint32 from f */
static uint32 uint32_read(FILE *f)
{
  uint32 u;

  u =  ((uint32) getc(f)) <<24;
  u |= ((uint32) getc(f)) <<16;
  u |= ((uint32) getc(f))  <<8;
  u |= ((uint32) getc(f));
  return u;
}

/* write an uint32 to f */
static void uint32_write(FILE *f, uint32 u)
{
  assert(f!=NULL);
  putc((u>>24) & 0xffL,f);
  putc((u>>16) & 0xffL,f);
  putc((u>> 8) & 0xffL,f);
  putc(      u & 0xffL,f);
}

static int32 int_log2(uint64 u)    
{
  int32 i = 1;
  uint64 r = 1;

  assert(u>0);
  while((i<=64) && (r<u)){
    r=2*r+1;
    i = i+1;
  }
  assert(i<=64);
  return i;
}
