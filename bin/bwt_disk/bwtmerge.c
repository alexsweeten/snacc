/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
   bwtmerge.c

   external memory sa construction algorithm based on the algorithm by
   Baeza-Yates, Gonnet and Snir, modified as suggested by Crauser and 
   Ferragina and with the further modification of using the rank-next array 

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include "bwtext.h"
#include <time.h>


// local prototypes
static void shell_sufsort(uint8 *t, uint32 *sa, uint32 n, uint32 *tn_rank);
static int cmp_with_gt_tn(uint8 *t, uint32 i,uint32 j,uint32 n,uint32 *gt_tn);
static int long_suffix_compare(uint64 i, uint8 *local_pfx, uint8 *tlast);
static uint32 get_char_from_first_column(uint32 index, uint32 *first_col);
static uint32 compute_first_segment_bwt(uint8 *text, uint8 *bwt, 
					uint32 n, uint32 *sa, uint32 *aux);
static void suf_insert_bwt(uint32 *bwt32,uint32 n,int32 last,uint32 *gaps,
			   uint64 n_oldsuf,uint32 *start_sa_range,uint32 rk0,uint32);
static void init_sb_b_counts(uint32 *bwt32, uint32 n, uint32 rk0,
			     uint32 sb_counts[][_BW_ALPHA_SIZE]);
static uint32 compute_prefix_rank(uint32 *,uint32,uint32,uint32 rk0,uint32 c, 
				  uint32 sb_counts[][_BW_ALPHA_SIZE]);
static uint32 isa2bwt(uint32 *isa, uint32 *bwt32, uint32 n, uint32 *first_col);
//static uint32 compute_prefix_rank_stupid(uint32 *bwt32,uint32 i,uint32 rk0,uint32 c);


// prototypes for kmp buffers
_bwtext_kbuf *_bwtext_kb_create(uint8 *s, uint32 size, uint32 bbsize);
void _bwtext_kb_add_char(_bwtext_kbuf *b, uint32 c, boolean gt);
boolean _bwtext_kb_revisit_char(_bwtext_kbuf *b, uint32 c, boolean *gt);
void _bwtext_kb_rewind(_bwtext_kbuf *b);
_bwtext_kbuf *_bwtext_kb_close(_bwtext_kbuf *b);

// prototype for internal suffix sorting
// void _bwtext_isufsort(uint8 *t, uint32 *sa, uint32 n, 
//		      uint32 *gt_tn, uint32 *gt_eof, uint32 *isa);
void _bwtext_isufsort_plain(uint8 *t, uint32 *sa, uint32 n,  uint32 *isa);
void _bwtext_isufsort_rank(uint8 *t, uint32 n, uint32 *gt_tn, uint32 *gt_eof, 
			   uint32 *isa, uint32 *aux); // computes rank



// external global variables
extern int _bwtext_Verbose;

// local global variables
#define PFX_BUFFER_SIZE _bwtext_PFX_BUFFER_SIZE
#define KMP_BIT_BUF_SIZE _bwtext_KMP_BIT_BUF_SIZE

static _bwtext_kbuf *kmp_in=NULL;
static _bwtext_kbuf *kmp_out=NULL;


static time_t Tot_saint_time=0;
static time_t Tot_sufins_time=0;
static time_t Tot_merge_time=0;
static time_t Tot_housekeeping_time=0;
uint64 Long_suffix_comparisons=0;

 
// read bit at position i from packed bit vector vec of uint32's
#define bit_read32(vec,i)  ( (vec[(i)/32] & (1u<<((i)%32))) ? 1u : 0 ) 
// write bit 1 to position i in logical bit vector vec
#define bit_set32(vec,i)  ( vec[(i)/32] |= ( 1u<<((i)%32)) )  



/* 
  we use uint32 integers to represent sa, rank, rank_next and so on.
  So internally we can compute sa of segments of size n with 
  n < 2^{32}-1 (we want that n+1 is still a legal uint32). 
  so every variable indexing t, sa, etc must be an uint32
*/


// name and filters of text and bwt files are passed trough 
// global variables
int64 _bwtext_sufsortmain(int32 imem_mb)
{
  uint32 *sa, *ranknext, *gaps, *isa;  // internal memory arrays with alias
  uint32 *gt_tn, *bit_aux, *bwt32;
  uint8 *t;
  uint64 last=0, first=0;         // portion of the input text considered
  uint64 occ_global[_BW_ALPHA_SIZE]; // # of occ. of each char in whole file  
  uint32 occ[_BW_ALPHA_SIZE];     // emulation of first column bwt matrix 
  uint32 start_sa_range[_BW_ALPHA_SIZE+1];
  uint32 size,n;                  // max block size and current block size
  uint32 i,j,new_rank0,gt_tn_size,tot; 
  uint32 last_symbol, new_rank_last;
  boolean more_to_come;           // true iff the current is not last block
  uint64 imem_size;               // size of usable internal memory
  int64 cbwt_size=0;              // size of the compressed (partial) bwt
  time_t tstart,tstart0;          // tstart0 never changes

  tstart0 = tstart = time(NULL);
  // ------------ open file containing text: if this fails ....
  _bwtext_text_open(); 
  
  // ------------ compute constants for memory allocation
  imem_size =((uint64)imem_mb)*1048576; // size of available internal memory
  assert( (imem_size/8) < 4294967295U); //1+(imem_size/8) must fit in a uint32
  size = (uint32) (imem_size/8);        // this is the block size
  assert(size>=PFX_BUFFER_SIZE);        // very short blocks do not make sense
  assert(size%256==0);                  // required by compute_prefix_rank
  uint64 mallocsize = sizeof(uint32) * (((uint64) size)+1);
  assert(mallocsize<SIZE_MAX);          // otherwise malloc wont work

  if(_bwtext_Verbose>0)
    fprintf(stderr,"Using blocks of size %u\n",size);

  // ------------ allocate internal memory
  sa = (uint32 *) malloc((size_t) mallocsize);
  isa = (uint32 *) malloc((size_t) mallocsize);
  if(!sa || !isa)  _bwtext_out_of_mem("bwtext_main");
  // define alias for the arrays that use the same memory
  bwt32 = ranknext = sa; gaps=isa;
  // allocate t and gt_tn inside isa 
  // make sure that all the array we want to fit into isa do not exceed
  // the available space. This could be a problem for larger alphabets
  gt_tn_size = (size+31)/32;           // = ceil(size/32)
  gt_tn = isa;                         // bit vector gt_tn[0..ceil(size/32)]
  bit_aux = (gt_tn + gt_tn_size);      // bit vector bit_aux[0..ceil(size/32)]
  t = (uint8 *) (gt_tn + 2*gt_tn_size);// t[0,2*size]
  assert(2*gt_tn_size*sizeof(uint32)+(2*size+1)*sizeof(uint8) <= 
         mallocsize*sizeof(uint32));

  // ------- init occ_global[]
  for(i=0;i<_BW_ALPHA_SIZE;i++) occ_global[i]=0;

  // ---------- processing of the first block
  assert(_bwtext_text_checkpos(0));  // check position in text
  // copy up to size uint8's starting with t[size-1], going backwards
  n=_bwtext_text_copy_reverse(t,size); 
  assert(n<=size);
  if (n==0) _bwtext_fatal_error("Empty input file");
  assert(_bwtext_text_checkpos(n));  // check position in text
  more_to_come = _bwtext_text_moretocome();
  for(i=size-n;i<size;i++) occ_global[t[i]]++;  // update occ_global[]
  Tot_housekeeping_time += time(NULL) - tstart;

  if(_bwtext_Verbose>1)
    fprintf(stderr,"Beginning working on block %d-%lld\n",0,(long long) n);

  // compute bwt and sa of the first block 
  tstart = time(NULL);
  new_rank0 = compute_first_segment_bwt(t+(size-n), t+size,n,sa,isa); 
  Tot_saint_time += time(NULL)-tstart;

  tstart = time(NULL);
  cbwt_size = _bwtext_bwt_init_temp(n+1,new_rank0+1,t+size); // init bwt file
  if(_bwtext_Verbose>2)
    fprintf(stderr,"Bwt file succesfully initialized\n");

  // ------------- prepare for next block 
  if(more_to_come) { 
    for(i=0;i<gt_tn_size;i++) gt_tn[i]=0;   // clean gt_tn vector
    for(i=new_rank0+1,j=0;i<n;i++) {        // init gt_tn vector
      assert(sa[i]<n && sa[i]!=0);
      //gt_tn[sa[i]/32] |= ((uint32) 1) <<(sa[i]%32); 
      bit_set32(gt_tn,sa[i]);
      j++;
    }
    // a final check
    assert(j==n-1-new_rank0);
    
    // ------------ copy first block into second half of t
    for(i=0;i<size;i++)
      (t+size)[i] = t[i];

    // ------------ init kmp buffer
    assert(kmp_in==NULL);
    kmp_in = _bwtext_kb_create(t, PFX_BUFFER_SIZE, KMP_BIT_BUF_SIZE);

  }
  Tot_housekeeping_time += time(NULL) - tstart;
  

  // ---------------  main loop in which we add the bwt of a new block 
  first = 0; last=n;      
  while(more_to_come) {

    /* invariants at the beginning of each iteration:
         last is the next position of the text we have to read
         t[size,2*size) contains the previous segment
         gt_tn[0,size) is properly initialized */
    assert(_bwtext_text_checkpos(last));  // check position in text

    tstart = time(NULL);
    // copy up to size uint8's starting with t[size-1], going backwards
    n=_bwtext_text_copy_reverse(t,size); 
    assert(n<=size && n>0);
    more_to_come = _bwtext_text_moretocome();

    first = last; last+=n;
    assert(_bwtext_text_checkpos(last));  // check position in text
    if(_bwtext_Verbose>1) {               // block we are beginnig to process
      fprintf(stderr,"-> %llu ", (unsigned long long) last);
      fprintf(stderr, 
	      "isa: %ld, sufi: %ld, mrg: %ld, hsk %ld, cbwt %lld,",
	      Tot_saint_time,Tot_sufins_time,Tot_merge_time,
	      Tot_housekeeping_time, (long long) cbwt_size);
      // report elapsed time in microsecond per symbol
      fprintf(stderr," muXs: %.2lf\n", (time(NULL)-tstart0)*(1000000.0/first));
    }

    // ---------- adjust kmp buffers
    // complete filling of kmp_in 
    for(i=size-1;i>0;i--)
      _bwtext_kb_add_char(kmp_in, (t+size)[i], (gt_tn[i/32]>>(i%32))&1);
    assert(kmp_in->chars_seen==first-1);
    // rewind and swap kmp_in <-> kmp_out
    _bwtext_kb_rewind(kmp_in);
    assert(kmp_out==NULL);
    kmp_out=kmp_in;
    kmp_in = NULL;
    assert(kmp_in==NULL);
    if(more_to_come)
      kmp_in = _bwtext_kb_create(t, PFX_BUFFER_SIZE, KMP_BIT_BUF_SIZE);

    // init start_sa_range for later emulation of first column bwt matrix;
    for(i=0;i<_BW_ALPHA_SIZE;i++) occ[i]=0;
    for(i=size-n;i<size;i++)       
      occ[t[i]]++;
    for(i=0,tot=0;i<_BW_ALPHA_SIZE;i++) {
      start_sa_range[i] = tot; tot+= occ[i]; occ_global[i] += occ[i];
    }
    start_sa_range[_BW_ALPHA_SIZE] = tot;  
    assert(tot==n);
    last_symbol = t[size-1];
    Tot_housekeeping_time += time(NULL) - tstart;

    // sort the suffixes in t[size-n ... size) 
    // Range: isa[0,n) -> [0,n)
    tstart = time(NULL);
    if(n<16) shell_sufsort(t+(size-n),sa,n,gt_tn);
    else _bwtext_isufsort_rank(t+(size-n),n,gt_tn,bit_aux,isa,sa);      
    Tot_saint_time += time(NULL)-tstart;
    if(_bwtext_Verbose>2)
      fprintf(stderr,"Suffix sorting done\n");

    // --------- compute current bwt 
    tstart = time(NULL);
    new_rank0 = isa[0];   // rank of t[first...]: we do not want to loose this
    new_rank_last = isa[n-1];  // rank of t[n-1...]
    j = isa2bwt(isa, bwt32, n, start_sa_range);
    assert(j==last_symbol);
    Tot_housekeeping_time += time(NULL) - tstart;
    // --------- suffix insertion by "backward-search"
    tstart = time(NULL);
    suf_insert_bwt(bwt32,n,last_symbol,gaps,first,start_sa_range,
		   new_rank0,new_rank_last);
    Tot_sufins_time += time(NULL) - tstart;

    
    // ------------- bwt merging
    // update global bwt given the gap array and the current bwt 
    tstart = time(NULL);
    cbwt_size = _bwtext_bwtmerge_temp(gaps,bwt32,n,new_rank0,last_symbol);
    Tot_merge_time += time(NULL) - tstart;

    // if this was the last segment stop
    if(!more_to_come) break;     
    
    tstart = time(NULL);
    // -------------- prepare for next block
    // compute rankprev in place (replacing ranknext==bwt)
    // Range: rankprev[0,n) -> [0,n]
    // rankprev[j]=n when j is the rank of t[0..]
    uint32 *rankprev=ranknext;
    for(i=0;i<n;i++)                        // compute rank_prev
      if(i!=new_rank0) {
	j = bwt32[i] & 255;                 // this is bwt[i]
	rankprev[i] = start_sa_range[j]++;
	if(rankprev[i]==new_rank_last)
	  rankprev[i]=start_sa_range[j]++;  //  we have to skip this rank
      }
      else
	rankprev[i] = n;

    for(i=0;i<gt_tn_size;i++) gt_tn[i]=0;  // clear bit vector
    i=new_rank_last;
    for(j=n-1;j>0;j--) {
      // invariant: i is the rank of t[j..]
      assert(i!=new_rank0);
      assert(i<n);
      if(i>new_rank0)
	bit_set32(gt_tn,j); //gt_tn[j/32] |= ( ((uint32) 1) << (j%32) );
      i = rankprev[i];
    }
    assert(i==new_rank0);
    
    // ---------- copy current block into second half of t
    n=_bwtext_text_copy_reverse(t+size,size); 
    assert(n==size);
    Tot_housekeeping_time += time(NULL) - tstart;

  }   // ---------- end main loop

  _bwtext_text_close();       // close text file
  _bwtext_bwt_finalize_temp();  // rename temp bwt file 
  _bwtext_bwt_write_auxfile(occ_global);  // write auxiliary info 

  if(_bwtext_Verbose>1) {

    fprintf(stderr, 
	    "Text size: %llu.\n", (long long) last);
    fprintf(stderr, 
	    "Compressed bwt size: (not exact) %lld.\n", (long long) cbwt_size);
    fprintf(stderr, 
	    "Long suffix comparisons: %llu.\n",
	    (unsigned long long) Long_suffix_comparisons);
    fprintf(stderr, 
	    "Total time spent building internal sa's: %ld secs.\n",
	    Tot_saint_time);
    fprintf(stderr, 
	    "Total time spent for suffix insertion: %ld secs.\n",
	    Tot_sufins_time);
    fprintf(stderr, 
	    "Total time spent merging: %ld secs.\n",
	    Tot_merge_time);
    fprintf(stderr, 
	    "Total time spent housekeeping: %ld secs.\n",
	    Tot_housekeeping_time);
  }

  // ----- deallocate
  assert(kmp_out==NULL && kmp_in==NULL);
  free(isa);
  free(sa);
  return last;
}


/* ****************************************************************
   compute the bwt of the first segment of the text
   \param text: text[0..n) is the text for which we need the bwt
                we have the guarantee that text[n...2n) is allocated
                and available for use.
   \param bwt: bwt[0..n] is the location in which we have to write 
               the bwt. Note that bwt[0...n] can overlap with t[n,2n)
   \param sa: array of size at least n+1 where we can write the bwt
   \param aux: array of size at least n+1 that can be used during 
               the computation. Note that aux may overlap t and bwt
   Output:
   the returned value is the rank of t[0..n)
   ************************************************************* */   
static uint32 compute_first_segment_bwt(uint8 *text, uint8 *bwt, 
				      uint32 n, uint32 *sa, uint32 *aux)
{
  uint32 i, new_rank0;

  _bwtext_isufsort_plain(text,sa,n,aux);      	
  if(_bwtext_Verbose>2)
    fprintf(stderr,"Suffix sorting done\n");

  bwt[0] = text[n-1];                  // text[n-1] is the first bwt char
  new_rank0=n;
  for(i=0;i<n;i++) 
    if(sa[i]==0) {
      new_rank0=i;                        // rank of text[0..]
      bwt[i+1] =bwt[i];                   // to help compression
    }
    else
      bwt[i+1]=text[sa[i]-1];
  assert(new_rank0<n);

  return new_rank0;
}

/* *******************************************************************
   fill the gap array computing the number of old suffixes 
   between two consecutive new suffixes

   For each suffix in t[last,bign) we find its lexicografic position
   within the suffix array of t[n-size,size) aka t[first,last)
   this is done using the bwt

   In gaps[i] we store the number of suffixes of t[last,bign) which
   are smaller that the one in sa[i] and larger than sa[i-1]. In 
   gaps[n] we store the number of suffixes which are greater than sa[n-1]
   
   The position of t[bign-1,bign) is simply the one
   immediately before the first suffix starting with t[bign-1]

   Note on kmp_buffers: 
     kmp_out has been just rewinded and is ready for revisitng
     kmp_in has been just initialized with a new string
   ******************************************************************* */
static void suf_insert_bwt(uint32 *bwt32, uint32 n, int32 last, uint32 *gaps, 
		       uint64 num_oldsuf, uint32 *start_sa_range,uint32 rk0,uint32 rklst)
{
  uint32 sb_counts[1+((n-1)>>24)][_BW_ALPHA_SIZE],j,cfirst,cur_rank,old;
  uint8 pfx_buffer[PFX_BUFFER_SIZE];
  uint64 i;
  int32 c;
  boolean long_lcp,gt;

  init_sb_b_counts(bwt32, n, rk0, sb_counts); // init counts for rank queries
  _bwtext_text_rewind();           // rewind the input text 
  for(j=0;j<=n;j++) gaps[j]=0;     // clear gaps

  // -------- eof symbol 
  gaps[0]++;                      
  // -------- last symbol (it is actually the first since we work on T^R)
  c = _bwtext_text_getc();               // last symbol of t
  assert(0<=c && c<_BW_ALPHA_SIZE);
  pfx_buffer[0] = (uint8) c;  
  cur_rank = start_sa_range[c];    // determining the rank of c$ is easy
  gaps[cur_rank] += 1;
  if(kmp_in) _bwtext_kb_add_char(kmp_in,c,cur_rank>rk0);
  // -------- all others symbols
  for(i=1;i<num_oldsuf;i++) {         // i = # symbols read from t
    long_lcp = _bwtext_kb_revisit_char(kmp_out,c,&gt); // revisit previous char
    c = _bwtext_text_getc();          // next symbol
    assert(0<=c && c<_BW_ALPHA_SIZE);
    cfirst = start_sa_range[c];
    old = cur_rank;
    if(cur_rank==0)
      cur_rank=cfirst;
    else {
      cur_rank=cfirst+compute_prefix_rank(bwt32,cur_rank-1,n,rk0,c,sb_counts);
      assert(cur_rank<=n);
    }
    if(c==last) {
      if(cur_rank==rklst) {
	if(long_lcp) {if(gt) cur_rank++;}    // t[i-1,bign) > t[last,bign)
	else if(long_suffix_compare(i-1,pfx_buffer,kmp_out->string)>0)
	  cur_rank++;
      }
      else if(cur_rank>rklst) {
	  cur_rank++;
	  if(long_lcp) assert(gt);    // t[i-1,bign) > t[last,bign)
	  else assert(long_suffix_compare(i-1,pfx_buffer,kmp_out->string)>0);
          // the last assertion should be costly to compute so is should
          // be removed. However apparently it makes execution faster 
      }
      else {
	if(long_lcp) assert(gt==0);    // t[i-1,bign) > t[last,bign)
        else assert(long_suffix_compare(i-1,pfx_buffer,kmp_out->string)<0);
          // the last assertion should be costly to computed so is should
          // be removed. However apparently it makes execution faster
      } 
    }
    pfx_buffer[i%PFX_BUFFER_SIZE] = (uint8) c;
    assert(cur_rank<=n);
    gaps[cur_rank] += 1;
    assert(gaps[cur_rank]>0);
    if(kmp_in) _bwtext_kb_add_char(kmp_in,c,cur_rank>rk0);    
  }
  // check that gaps accounts for all old suffixes
  for(i=0,j=0;j<=n;j++)
    i+= gaps[j];
  assert(i==num_oldsuf+1);             // +1 is for the empty suffix
  kmp_out = _bwtext_kb_close(kmp_out); // deallocate and set kmp_out=NULL
}


/* ===================================================================
   init_sb_b_counts() and compute_prefix_rank() are used together to
   support rank quesries on a prefix of the bwt. Their space usage 
   is 4n + (256*4*n)/2^24 = 4n + n/2^16 bytes.
   Bucket size is 256, superbucket size is 2^24 independently of n
   The bwt is stored in the low 8 bits of a uint32 array, the 24
   highest bits are used for the bucket table. See below for details.

   Important: this scheme requires that the bwt is stored in an array
   whose size is a multiple of 256 (n can be not a multiple of 256 but
   in this case there must be free space in bwt up to the next position
   whose index i = 255 mod 256). Unfortunately we cannot check this 
   here, so it is checked in _bwtext_sufsortmain()
   =================================================================== */
   

/* ****************************************************************
   init the bucket and superbucket counts for the array bwt32[0,n)
   excluding the position rk0
   at exit superbucket[i][c] contains the number of occurrences 
   of character c in bwt[0,i*2^24). 
   in addition if x is a multiple of 256, at exit 
   bwt32[x+c]>>8 is the number of occurrences of the character c
   from the beginning of the superbucket containing x up to position x-1
   **************************************************************** */
#define Two224 (1u<<24)        // superbucket size =2^24
#define Mask24 (Two224-1)      // 24 low order bits set to one
#define Mask8 ((uint32) 255)
static void init_sb_b_counts(uint32 *bwt32, uint32 n, uint32 rk0,
		      uint32 sb_counts[][_BW_ALPHA_SIZE])
{
  uint32 count[_BW_ALPHA_SIZE];
  uint32 i,j,num,sb,tot;

  // -------- the scheme requires that bwt32 has size =0 mod 256
  for(j=n;(j&Mask8)!=0;j++)
    bwt32[j] = 0;             // clear up to the next position multiple of 256
  for(j=0;j<_BW_ALPHA_SIZE;j++) count[j]=0; // init global counts
  for(i=sb=0;i<n;i++) {
    if((i&Mask24)==0) {              // i multiple of 2^24
      for(j=0;j<_BW_ALPHA_SIZE;j++)
	sb_counts[sb][j] = count[j];
      sb++;                          // increment number of superbucket
    }
    if((i&Mask8)==0) {               // i multiple of 2^8
      for(j=0;j<_BW_ALPHA_SIZE;j++) {
	num = count[j]-sb_counts[sb-1][j];
	if(num> Two224-256) 
	  _bwtext_fatal_error("init_sb_b_counts");
	bwt32[i+j] += num<<8; // we write up to 255 positions after i
      }
    }
    if(i!=rk0) count[bwt32[i]&Mask8]++; // increment count  
  }
  // check that the total number of occurrences is n-1
  tot=0;
  for(j=0;j<_BW_ALPHA_SIZE;j++) 
    tot += count[j];
  assert(tot==n-1);
  assert( ((n-1)>>24) == sb-1 ); 
}


/* ****************************************************************
    compute the number of occurrences of charater c in bwt32[0,i]
    the position bwt32[rk0] should not be considered because
    it should contain the eof symbol
   **************************************************************** */
static uint32 compute_prefix_rank(uint32 *bwt32,uint32 i,uint32 n,uint32 rk0, 
				  uint32 c,uint32 sb_counts[][_BW_ALPHA_SIZE])
{
  uint32 cur_bucket_start, next_bucket_start, j;
  uint32 sb_count;                         // count up to the beginning of sb
  uint32 b_count;                            // count up to the beginning of b
  uint32 bend_count;                         // count up to the end of b
  uint32 this;                               // count in this bucket  
  
  sb_count = sb_counts[i>>24][c];            // up to the beginning of sb
  cur_bucket_start = i & (~Mask8);
  b_count = bwt32[cur_bucket_start+c]>>8;    // up to the beginning of bucket
  next_bucket_start = cur_bucket_start+256; 
  if( next_bucket_start< n) {
    // get count at the end of current bucket
    if( (next_bucket_start & Mask24) == 0) { // next b is in different sb
      assert(i>>24 != next_bucket_start>>24);
      bend_count = sb_counts[next_bucket_start>>24][c]-sb_count;
    }
    else
      bend_count = bwt32[next_bucket_start+c]>>8;
    assert(bend_count>=b_count);
    if(bend_count==b_count)
      return sb_count + b_count;
    if((i& 128) != 0) {              // equivalet to i%256>=128
      for(this=0,j=i+1;j<next_bucket_start;j++)// note: next_bucket_start<n
	if( (j!=rk0) && ((bwt32[j] & Mask8)==c) ) 
	  this++;
      return sb_count + bend_count - this;
    } 
  }
  // get contribution of this bucket
  for(this=0, j=cur_bucket_start; j<=i; j++)
    if( (j!=rk0) && ((bwt32[j] & Mask8)==c) ) 
      this++;
  return sb_count+b_count+this;
}

#if 0
static 
uint32 compute_prefix_rank_stupid(uint32 *bwt32,uint32 i,uint32 rk0,uint32 c)
{
  uint32 j,tot=0;
  for(j=0;j<=i;j++)
    if( (j!=rk0) && ((bwt32[j] & Mask8)==c) ) 
      tot++;
  return tot;
}
#endif


/* ********************************************************
   compare suffix with t[i...] with t[last...]
   returns <0 if t[i...] <t[last] and >0 if t[i...]>t[last...]
   by construction the lcp between these two suffixes is at 
   most Pfx_buffer_size symbols
   ******************************************************** */ 
static int long_suffix_compare(uint64 i, uint8 *local_pfx, uint8 *tlast)
{
  uint32 k;

  Long_suffix_comparisons++;
  for(k=PFX_BUFFER_SIZE-1; ; ) {
    if(tlast[k]>local_pfx[i%PFX_BUFFER_SIZE])
      return -1;   // t[last..] is larger
    if(tlast[k]<local_pfx[i%PFX_BUFFER_SIZE])
      return 1;    // t[last..] is smaller
    // the current symbols are equal: step by one
    if(k-- ==0)
      _bwtext_fatal_error("Illegal lcp in long_suffix_compare()");
    if(i-- ==0) 
      return -1;  // t[i..] is smaller
  }
  assert(0);   // the execution should never reach here
  return 0;    // does not matter
}


/* *******************************************************************
   compute the bwt given the isa array 
   return the last text character t[n-1]
   \param isa:       the inverse suffix array
   \param bwt32:     array where the bwt will be written
   \param n:         number of suffixes
   \param first_col: compressed representation of the first column
                     of the bwt matrix. first_col[i] is the number of
                     symbols smaller than i
   ******************************************************************* */
static uint32 isa2bwt(uint32 *isa, uint32 *bwt32, uint32 n, uint32 *first_col)
{
  uint32 i, c;

  c=0;
  for(i=0;i<n;i++) {
    bwt32[isa[i]] = c;
    c = get_char_from_first_column(isa[i],first_col);
  }
  // take care of the position where the eof symbol should go
  // writing a nearby symbol to help (run-length) compression 
  if(isa[0]>0) 
    bwt32[isa[0]]=bwt32[isa[0]-1]; // preceding symbol
  else if(isa[0]<n-1)
    bwt32[isa[0]]=bwt32[isa[0]+1]; // following symbol 
  else 
    bwt32[isa[0]]=c;  // this is certainly in t
  return c;           // this is t[n-1]
}

// given an index returns the corresponding char in the bwt matrix
// to be improved using binary search
static uint32 get_char_from_first_column(uint32 index, uint32 *first_col)
{
  uint32 i;

  for(i=0;i<_BW_ALPHA_SIZE;i++)
    if(first_col[i]<=index && index < first_col[i+1])
      return i;
  _bwtext_fatal_error("get_char_from_first_column");
  return 0;  
}


/* *******************************************************************
   given a string t[0...] (the endpoint is unknown) compute the 
   ordering of the suffixes t[0...] up to t[n-1...], (n in total),
   and write it to sa[0,n). The sorting is possible since 
   in addition to t[0,n) we are provided with t[n,n+size) and with a 
   (packed) bit vector gt_tn[] providing the rank of t[n...] with
   respect to t[j...] for j=n ... size-1 (see below). 
   The comparison between two suffixes t[i...] and t[j...] goes as follows:
     1. compare character by character. if a difference is found OK
     2. as soon as min(i+k,j+k)=n use gt_tn
   input 
     \param t:    input text. length=n+size 
     \param sa:   space for the sa.  length = n
     \param n:    we have to sort the suffixes t[0...] to t[n-1...]
     \param gt_tn: (reads: larger_than_t[n...])
                   logically it is a bit vector of length size,
                   Provides the rank of the suffix t[n...] within
                   the suffixes t[j...] for j=n ... n+size-1 in the sense
                   t[j...] > t[n...] iff gt_tn[j-n]==1
                   In practice it is a uint32 array of length ceil(size/32).
   output
     the sorted order of suffixes  t[0 ...] to t[n-1...] is stored sa
   ******************************************************************* */
static void shell_sufsort(uint8 *t, uint32 *sa, uint32 n, uint32 *gt_tn) 
{
  uint32 i, j, h, k, tmp;
  int incs[16] = { 1391376, 463792, 198768, 86961, 33936,
		   13776, 4592, 1968, 861, 336, 112, 48, 21, 7, 3, 1 };

  assert(n>0);

  // init suffix array
  for(i=0;i<n;i++) sa[i]=i;

  // shellsort
  for ( k = 0; k < 16; k++) {
    h = incs[k]; 
    for (i = h; i < n; i++) { 
      tmp = sa[i]; 
      for(j=i; j>= h; j -= h) {
	if(cmp_with_gt_tn(t, tmp, sa[j-h],n,gt_tn) < 0) 
	  sa[j] = sa[j-h];
	else
	  break;
      }
      sa[j] = tmp;
    } 
  }
}


// compare t[i..] and t[j..] using the bit vector gt_tn[]
// returns segnum(t[i..]-t[j..])
static int cmp_with_gt_tn(uint8 *t, uint32 i,uint32 j,uint32 n,uint32 *gt_tn) 
{
  uint32 r;
  
  assert(i!=j);
  assert(i<n && j<n);
  if(i<j) {
    while(i<n) {
      if(t[i]<t[j]) return -1;
      else if(t[i]>t[j]) return 1;
      i++; j++;
    }
    assert(i==n);
    r = j-n;
    assert(r>0 && r<n);
    if( bit_read32(gt_tn,r) ) // gt_tn[r/32] & (bit << (r%32))
      return -1;
    return 1;
  }
  else {  // here j<i
    while(j<n) {
      if(t[i]<t[j]) return -1;
      else if(t[i]>t[j]) return 1;
      i++; j++;
    }
    assert(j==n);
    r = i-n;
    assert(r>0 && r<n);
    if( bit_read32(gt_tn,r) ) //  gt_tn[r/32] & (bit << (r%32))
      return 1;
    return -1;
  }
}

