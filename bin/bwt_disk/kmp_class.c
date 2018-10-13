#include "bwtext.h"

#define simple_buffer 0
#define long_run_limit 16

// global prototypes
void _bwtext_fatal_error(char *s);
void _bwtext_out_of_mem(char *f);


extern int _bwtext_Verbose;


// local prototypes
static void write_run(boolean, uint64 n, uint32 *pos, uint32 size, uint32 *buf);
static uint64 read_run(boolean *, uint32 *pos, uint32 size, uint32 *buf);


/* **************************************************************
   create a kmp_buffer for s[0,size) allocating and initializing 
   b->string, b->kmp_shift and allocating b->bit_buffer. 
   ************************************************************** */
_bwtext_kbuf *_bwtext_kb_create(uint8 *s, uint32 size, uint32 bbsize)
{
  uint32 k,q;
  _bwtext_kbuf *b;

  assert(size<65535);
  b = (_bwtext_kbuf *) malloc(sizeof(*b));
  if(!b)  _bwtext_out_of_mem("_bwtext_kb_create");
  b->string = (uint8 *) malloc(size*sizeof(uint8));
  b->kmp_shift = (uint16 *) malloc((size+1)*sizeof(uint16));
  b->bit_buffer = (uint32 *) calloc( bbsize/32, sizeof(uint32));  
  if(!b->string || !b->kmp_shift ||  !b->bit_buffer)
    _bwtext_out_of_mem("_bwtext_kb_create");
  b->size=size;                      // size of the string
  b->bit_buffer_size=32*(bbsize/32); // size of the bit_buffer in bits 
  b->cur_bit = b->current=0;
  b->pending0 = b->pending1 = 0;
  b->chars_seen=0;

  // copy string to b->string
  for(k=0;k<size;k++)
    b->string[k] = s[size-1-k];
  
  // compute kmp shift 
  // kmp_shift[i]=j => string[j-1...] pfx of string[i-1...]
  b->kmp_shift[0] = 0;        // this entry is never used
  b->kmp_shift[1] = 0;
  k=0;
  for(q=1;q<size;q++) {
    // invariant: string[k-1...]  pfx of string[q-1...] 
    assert(k<q);
    if(b->string[q] == b->string[k]) {
      // it is string[k...]  pfx of string[q...]
      b->kmp_shift[q+1] = ++k;
      continue;
    }
    while(k>0) {
      assert(k<q);
      assert(b->string[q] != b->string[k]);
      k = b->kmp_shift[k];
      if(b->string[q] == b->string[k]) {
	b->kmp_shift[q+1] = ++k;
	break;
      }
    }
    if(k==0)
      b->kmp_shift[q+1] = 0;	
  }
  // redundant check on kmp_shift
  for(k=1;k<=size;k++) {
    assert(k>b->kmp_shift[k]);
    for(q=0;q<b->kmp_shift[k];q++)
      assert(b->string[q]==b->string[q+k-b->kmp_shift[k]]);
  }
  return b;
}

// add one char to the virtual buffer. if buffer matches string 
// record bit gt in bit_buffer 
void _bwtext_kb_add_char(_bwtext_kbuf *b, uint32 c, boolean gt)
{
  assert(c<256);
  assert(gt==0 || gt==1);
  b->chars_seen++;
  if(c==b->string[b->current]) {
    if(++b->current==b->size) {
       b->current = b->kmp_shift[b->current];
#if simple_buffer
      if(gt)
	b->bit_buffer[b->cur_bit/32] |=  1<<(b->cur_bit%32);
      if(++b->cur_bit==b->bit_buffer_size)
	_bwtext_fatal_error("kmp bit buffer full (_bwtext_kb_add_char)");
#else
      if(gt) {  // store 1 bit
	if(b->pending0>0) {
	  assert(b->pending1==0);
	  write_run(0,b->pending0,&(b->cur_bit),
		    b->bit_buffer_size,b->bit_buffer);
	  b->pending0 = 0;
	  b->pending1 = 1;
	}
	else b->pending1++;
      }
      else { // gt==0
	if(b->pending1>0) {
	  assert(b->pending0==0);
	  write_run(1,b->pending1,&(b->cur_bit),
		    b->bit_buffer_size,b->bit_buffer);
	  b->pending1 = 0;
	  b->pending0 = 1;
	}
	else b->pending0++;
      }
#endif
    }
    return;
  }
  // c does not match
  while(b->current >0) { 
    assert(c!=b->string[b->current]);   // loop invariant 
    b->current = b->kmp_shift[b->current];
    if(c==b->string[b->current]) {
      b->current++;
      return;
    }
  }
  assert(c!=b->string[b->current]);
  assert(b->current==0);
  return;
}




// add one char to the virtual buffer. if buffer matches string 
// return true and  store bit from bit_buffer in gt 
boolean _bwtext_kb_revisit_char(_bwtext_kbuf *b, uint32 c, boolean *gt)
{

  assert(c<256);
  b->chars_seen--; 
  if(c==b->string[b->current]) {
    if(++b->current==b->size) {
      b->current = b->kmp_shift[b->current];
      // get bit from bit_buffer
#if simple_buffer
      if(b->cur_bit>=b->stored_bits)
	_bwtext_fatal_error("kmp bit buffer empty (_bwtext_kb_revisit_char)");
      *gt = (b->bit_buffer[b->cur_bit/32]>>(b->cur_bit%32)) & 1;
      b->cur_bit++;
#else
      if(b->pending0>0) {
	assert(b->pending1==0);
	b->pending0--;
	*gt=0;
      }
      else if(b->pending1>0) {
	assert(b->pending0==0);
	b->pending1--;
	*gt=1;
      }
      else { // get next run from buffer
	assert(b->pending0==0 && b->pending1==0);
	uint64 len=read_run(gt,&(b->cur_bit),b->stored_bits,b->bit_buffer);
	assert(len>0);
	if(*gt)   b->pending1 = len-1;
	else      b->pending0 = len-1;
      }
#endif
      assert(*gt==1 || *gt==0);
      return 1;
    }
    return 0;
  }
  // c does not match
  while(b->current >0) { 
    assert(c!=b->string[b->current]);   // loop invariant 
    b->current = b->kmp_shift[b->current];
    if(c==b->string[b->current]) {
      b->current++;
      return 0;
    }
  }
  assert(c!=b->string[b->current]);
  assert(b->current==0);
  return 0;
}


// prepare the structure for the revisting phase
void _bwtext_kb_rewind(_bwtext_kbuf *b)
{
#if !simple_buffer
  if(b->pending0>0) {
    assert(b->pending1==0);
    write_run(0,b->pending0,&(b->cur_bit),b->bit_buffer_size,b->bit_buffer);
  }
  else if(b->pending1>0) {
    assert(b->pending0==0);
    write_run(1,b->pending1,&(b->cur_bit),b->bit_buffer_size,b->bit_buffer);
  }
#endif
  b->pending0=b->pending1=0;
  b->stored_bits = b->cur_bit;
  b->cur_bit = b->current = 0;
  if (_bwtext_Verbose>1 &&  b->stored_bits>0)
    fprintf(stderr,"!%u bit(s) stored in bit buffer!\n", b->stored_bits);
}

// check that every written bit has been read and deallocate kbuffer
_bwtext_kbuf *_bwtext_kb_close(_bwtext_kbuf *b)
{
  assert(b->cur_bit==b->stored_bits);
  assert(b->chars_seen==0);
  free(b->bit_buffer); 
  free(b->kmp_shift);  
  free(b->string);     
  return NULL;
}


// ################### static functions ###########

// compute the len of the gamma code for n
inline int32 gammacode_len(uint64 n)
{
  int32 b=0;

  assert(n>0);
  while (n) { b++; n >>= 1; }
  assert(b>0 && b <=64);
  return 2*b-1;
}


static 
void write_run(boolean val, uint64 n, uint32 *pos, uint32 size, uint32 *buffer)
{
  int32 i;
  assert(val==0 || val==1);
  assert(n>0);

  // fprintf(stderr,"---> [%d %lld] (bit buffer)\n",val,n);

  uint32 run =  (n<long_run_limit) ? ((uint32) n) : long_run_limit;

  // write run bits
  if(*pos+run>=size) 
    _bwtext_fatal_error("kmp bit buffer full (write_run)");
  if(val==0)
    *pos += run; // only advance pointer
  else  //val==1
    for(i=run;i>0;i--) {
      buffer[(*pos)/32] |=  1<<((*pos)%32);
      *pos += 1; 
    }

  // if n<long_run_limit we are done, otherwise we have to write 
  // 1+n-long_run_limit using gamma code 
  if(n>=long_run_limit) { 
    uint64 v = 1+(n-long_run_limit);
    int32 len = gammacode_len(v);
    if(*pos+len>=size) 
      _bwtext_fatal_error("kmp bit buffer full (write_run)");
    *pos += (len-1)/2; // "write" prefix of zeroes 
    for(i=(len-1)/2;i>=0;i--) { // write binary digits 
      if((v>>i)&1)
	buffer[(*pos)/32] |=  1<<((*pos)%32);
      *pos += 1; 
    }
  }
}


static 
uint64 read_run(boolean *val, uint32 *pos, uint32 size, uint32 *buffer)
{
  uint64 count;

  if( *pos >= size)
    _bwtext_fatal_error("kmp bit buffer empty (read_run)");
  *val = (buffer[(*pos)/32] >> ((*pos)%32)) &1;
  (*pos)++;
  count=1;

  while(count<long_run_limit && *pos < size) {
    if(*val!= ((buffer[(*pos)/32] >>((*pos)%32)) &1)) { 
      // fprintf(stderr,"<--- [%d %lld] (bit buffer)\n",*val,count);
      return count;  // exit if we find a bit different from the first one
    }
    count++;
    (*pos)++;
  }
  if(count==long_run_limit) {  // there is a gammacode to add
    int32 bit,nbits=0;
    // get number of digits
    while(1) { 
      if( *pos >= size)
	_bwtext_fatal_error("kmp bit buffer empty (read_run)");
      bit = ((buffer[(*pos)/32]>>((*pos)%32)) &1);
      (*pos)++;
      if(bit==1) break;
      nbits++;
      assert(nbits<64);
    }
    // read actual digits
    uint64 gamma=1;            // most significant bit 
    for( ;nbits>0;nbits--) {
      if( *pos >= size)
	_bwtext_fatal_error("kmp bit buffer empty (read_run)");
      bit = ((buffer[(*pos)/32]>>((*pos)%32)) &1);
      (*pos)++;
      gamma = 2*gamma+bit;
    }
    count += gamma-1;
  }
  // fprintf(stderr,"<--- [%d %lld] (bit buffer)\n",*val,count);
  return count;
}
















// =========== non usati

#if 0
static void write7x8(uint64 n, uint32 *pos, uint32 size, uint8 *buffer)
{
  uint64 t;

  do {
    t = n & 0x7f;       // takes last seven bits
    n = n >> 7;
    if(n==0) t |= 0x80;  // set 8th bit
    buffer[(*pos)++] = (uint8) t;    
    if(*pos>= size) 
      _bwtext_fatal_error("kmp bit buffer full (kmp_class::write7x8)");
  } while(n>0);
}

static uint64 read7x8(uint32 *pos, uint32 size, uint8 *buffer)
{
  int32 i;
  uint64 c,ris=0;

  c =  buffer[(*pos)++];
  if(*pos>=size)
    _bwtext_fatal_error("kmp bit buffer empty (kmp_class::read7x8)");
  for(i=0;  ;i+=7) {
    assert(i<=63);
    ris |= (c& 0x7f) << i;      // add 7 bits 
    if((c&0x80)!=0) break;      // more bits ?
    c =  buffer[(*pos)++];
    if(*pos>=size)
      _bwtext_fatal_error("kmp bit buffer empty (kmp_class::read7x8)");
  }
  return ris;
}

#endif


