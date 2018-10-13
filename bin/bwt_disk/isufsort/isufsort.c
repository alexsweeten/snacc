#include <string.h>
#include <assert.h>
#include "../bwtext.h"
#include "libdivsufsort-lite/divsufsort.h"



// read bit at position i from packed bit vector vec
#define bit_read32(vec,i)  ( (vec[(i)/32] & (1u<<((i)%32))) ? 1u : 0 ) 
// write bit 1 to position i in logical bit vector vec
#define bit_set32(vec,i)  ( vec[(i)/32] |= ( 1u<<((i)%32)) )  



// global prototypes
void _bwtext_fatal_error(char *s);
void _bwtext_out_of_mem(char *f);
void _bwtext_ls_suffixsort(int *x, int *p, int n, int asize);
void _bwtext_ls_suffixsort_rank(int *x, int *p, int n, int asize);
void _bwtext_ls_suffixsort_plain(int *x, int *p, int n, int asize);


// local prototypes
static uint32 remap_alphabet(uint8 *t, uint32 n, uint32 *newt, uint32 *gt_eof);
static void compute_gt_eof(uint8 *, uint32, uint32 *, uint32 *, uint32 *);
static void compute_kmp_shift(uint8 *t, uint32 n, uint32 *shift);
static void ls_sufsort_rank(uint32 *isa, uint32 *sa, uint32 n, uint32 asize);
static void ls_sufsort_plain(uint32 *isa, uint32 *sa, uint32 n, uint32 asize);



extern int _bwtext_Verbose;
extern int _bwtext_Use_divsufsort;


/* *******************************************************************
   given a string t[0..n-1] compute the its suffix array, that is the 
   ordering of the suffixes t[0...n-1] up to t[n-1...n-1], (n in total),
   and write it to sa[0,n). 
   input 
     \param t:      input text. length=n
     \param sa:     space for the sa.  length = n+1
     \param n:      we have to sort the suffixes of t[0,n)
     \param isa:    available array of n+1 uint32's. In LS it temporary stores
                    the inverse suffix (=rank) array  
                    Important note:
                    the memory used by isa overlaps with the memory used by
                    t so they cannot be used simultaneously. 
   output
     the sorted order of suffixes  t[0 ...] to t[n-1...] is stored sa[0..n-1]
     t[0,n) is restored
   ******************************************************************* */
void _bwtext_isufsort_plain(uint8 *t, uint32 *sa, uint32 n,  uint32 *isa) 
{

  if(_bwtext_Use_divsufsort) {
    int ris = divsufsort(t, (int32 *) sa, n);
    if(ris!=0) 
      _bwtext_fatal_error("Divsufsort Error (_bwtext_isufsort_plain)"); 
    return;
  }
  else {
    uint32 occ[_BW_ALPHA_SIZE], map[_BW_ALPHA_SIZE];
    uint32 i,j,c,asize;
    
    // count occurrences of each symbol and remap alphabet
    memset( occ, 0, _BW_ALPHA_SIZE* sizeof(uint32) );
    for(i=0;i<n;i++) occ[t[i]]++;
    
    // compute the mapping between new and old symbols
    asize=0;
    for(i=0;i<_BW_ALPHA_SIZE;i++) {
      if(occ[i]>0) {
	map[i] = asize++;
	// fprintf(stderr,"%c (%d occ) -> %d\n",i,occ[i],map[i]);
      } else
	map[i] = 2+_BW_ALPHA_SIZE;
    }
    
    // do the mapping writing the remapped text to sa[]
    for(i=0;i<n;i++) {
      c = t[i];
      assert(map[c]< 2+_BW_ALPHA_SIZE );
      sa[i] = map[c]; // actual remapping
    }
    
    // now we no longer need t
    // copy the remapped text to isa
    for(i=0;i<n;i++) isa[i] = sa[i];
    
    // do the sorting (for now with larsson-sadakane)
    ls_sufsort_plain(isa,sa,n,asize);
    
    // now sa contains the suffix array and isa its inverse
    // recover the text in t (overwriting isa)
    j=0;
    for(i=0;i<n;i++) {
      while(occ[j]==0) {         // search next nonzero in occ 
	j++; assert(j<_BW_ALPHA_SIZE);
      }
      t[sa[i]] = (uint8) j;
      occ[j]--;
    }

    // check that all symbols have been assigned 
    assert(occ[j]==0);
    for(i=j+1;i<_BW_ALPHA_SIZE;i++)
      assert(occ[i]==0);
    // now t[0,n) contains the text and sa[0,n) the suffix array
    return;
  }
}

/* ***************************************************************
   plain (with usual eof rule) suffix sorting algorithm 
   using larsson-sadakane.
   \param: isa[0,n): the input string over the alphabet [0,asize)
           the algorithm append isa[n]=0 so isa must be of length n+1 
   \param: sa[0,n): array where the sa is written.
   *************************************************************** */
static void ls_sufsort_plain(uint32 *isa, uint32 *sa, uint32 n, uint32 asize)
{
  uint32 i;
  assert(n+1< (1u<<31));
  assert(n>0);

  _bwtext_ls_suffixsort_plain((int32 *)isa,(int32 *)sa,(int) n,(int) asize);
  assert(sa[0]==n);
  for(i=0;i<n;i++) // shift by one position to get rid of the eof
    sa[i] = sa[i+1];

}



/* *******************************************************************
   given a string t[0...] (the endpoint is unknown) compute the 
   ordering of the suffixes t[0...] up to t[n-1...], (n in total),
   and write the relative rank of these suffixes to isa[0] up to isa[n-1]
   Since isa[] and t[] share the same memory space the isa[] overwrites t[].

   The sorting is possible since in addition to t[0,n) we are provided with
   t[n,2n) and with a (packed) bit vector gt_tn[] providing us with the rank
   of t[n...] with respect to t[j...] for j=n ... 2n-1 (see below).
 
   The comparison between two suffixes t[i...] and t[j...] goes as follows:
     1. compare character by character. if a difference is found OK
     2. as soon as min(i+k,j+k)=n use gt_tn
   input 
     \param t:      input text. length=2n
     \param n:      we have to sort the suffixes t[0...] to t[n-1...]
     \param gt_tn:  (reads: larger_than_t[n...])
                    logically it is a bit vector of length n,
                    Provides the rank of the suffix t[n...] within
                    the suffixes t[j...] for j=n ... 2n-1 in the sense
                    t[j...] > t[n...] iff gt_tn[j-n]==1
                    In practice it is a uint32 array of length ceil(n/32).
                    If gt_tn==NULL then t[n..] should be considered
                    the lex. smallest suffix (used for first sa computation)
     \param gt_eof: available (packed) bit array of size n
     \param isa:    array of n+1 uint32's where the ranks are written.
                    Important note:
                    the memory used by isa overlaps with the memory used by
                    t, gt_tn and gt_eof so they cannot be used simultaneously. 
     \param aux:    auxiliary uint32 array of length = n+1

   output
     the relative ranks of suffixes  t[0 ...] to t[n-1...] is stored 
     in isa[0..n-1], t, gt_tn and aux are overwritten
   ******************************************************************* */
void _bwtext_isufsort_rank(uint8 *t, uint32 n, uint32 *gt_tn, uint32 *gt_eof, 
			    uint32 *isa, uint32 *aux) 
{
  uint32 i, newasize;


  // compute the packed bit vector gt_eof[] such that
  // for i=0,..,n-1  gt_eof[i]==1 iff t[i..] > t[n..]
  compute_gt_eof(t,n,aux,gt_tn,gt_eof);

  // remap symbols in t[] according to gt_eof bits: the new string is written
  // to sa[] and the remapping is such that suffix-sorting of sa[0,n-1]
  // is equivalent to sorting the suffixes t[0..], t[n-1...]
  newasize = remap_alphabet(t,n,aux,gt_eof);


  if(newasize<256 && _bwtext_Use_divsufsort) {
    // ------ it is possible to use divsufsort
    if(_bwtext_Verbose>2) fprintf(stderr,"DV ");
    assert(n< (1u<<31));
    for(i=0;i<n;i++) {
      assert(aux[i]<256);
      t[i] = (uint8) aux[i];
    }
    int ris = divsufsort(t,(int *) aux,n); // compute sa in aux
    if(ris!=0) 
      _bwtext_fatal_error("Divsufsort Error (_bwtext_isufsort_rank)"); 
    for(i=0;i<n;i++) {
      assert(aux[i]<n);
      isa[aux[i]] = i; 
    }
  }
  else {
    // ------- we have to use larsson-sadakane
    if(_bwtext_Verbose>1) fprintf(stderr,"LS ");

    // now we no longer need t, gt_tn, and gt_eos
    // copy the remapped text to isa
    for(i=0;i<n;i++) isa[i] = aux[i];

    // do the sorting (for now with larsson-sadakane)
    ls_sufsort_rank(isa,aux,n,newasize);
  }

  return;
}

static void ls_sufsort_rank(uint32 *isa, uint32 *aux, uint32 n, uint32 asize)
{
  assert(n< (1u<<31));
  assert(n>0);

  _bwtext_ls_suffixsort_rank((int32 *) isa,(int32 *)aux,(int) n-1,(int) asize);
}





static uint32 remap_alphabet(uint8 *t, uint32 n, uint32 *newt, uint32 *gt_eof)
{
  uint32 occ[_BW_ALPHA_SIZE+2],i,c;
  uint32 map[_BW_ALPHA_SIZE+2], asize;

  // count # of occurrences of each new symbol
  memset( occ, 0, (2+_BW_ALPHA_SIZE)* sizeof(uint32));
  for(i=0;i<n-1;i++) {
    if (t[i]<t[n-1] || (t[i]==t[n-1] && bit_read32(gt_eof,i+1)==0)) {
      occ[t[i]]++;
      assert(occ[t[i]]>0);
    } else if (t[i]>t[n-1] || (t[i]==t[n-1] && bit_read32(gt_eof,i+1)==1)) {
      occ[t[i]+2]++;
      assert(occ[t[i]+2]>0);
    }
    else 
      assert(0);
  }
  occ[t[n-1]+1]++;
  assert(occ[t[n-1]+1]==1);

  // compute the mapping between new and old symbols
  asize=0;
  for(i=0;i<2+_BW_ALPHA_SIZE;i++) {
    //fprintf(stderr,">> %x: %d\n",i,occ[i]);//!!!!!!!!!!
    if(occ[i]>0) {
      map[i] = asize++;
      // fprintf(stderr,"%c (%d occ) -> %d\n",i,occ[i],map[i]);
    } else
      map[i] = 2+_BW_ALPHA_SIZE;
  }

  // do the mapping
  for(i=0;i<n;i++) {
    if(i==n-1)
      c = t[i]+1;
    else if(t[i]<t[n-1] || (t[i]==t[n-1] && bit_read32(gt_eof,i+1)==0))
      c = t[i];
    else if (t[i]>t[n-1] || (t[i]==t[n-1] && bit_read32(gt_eof,i+1)==1)) 
      c = t[i]+2;
    else 
      assert(0);
    assert(map[c]< 2+_BW_ALPHA_SIZE );
    newt[i] = map[c]; // actual remapping
    // fprintf(stderr,"%d %d\n",t[i], newt[i]);
  }
  return asize;
}

#if 0
  // compute cumulative sums 
  sum=0;
  for(i=0;i<_BW_ALPHA_SIZE+2;i++) {
    old_sum = sum;
    sum += occ[i];
    assert(sum>=old_sum);
    occ[i]=old_sum;
  }
  assert(sum==n);
  // occ[] now contains the remapping
#endif

#if 0
static uint32 remap_alphabet_plain(uint8 *t, uint32 n, uint32 *newt)
{
  uint32 occ[_BW_ALPHA_SIZE],i,c;
  uint32 map[_BW_ALPHA_SIZE], asize;

  // count # of occurrences of each new symbol
  memset( occ, 0, (_BW_ALPHA_SIZE)* sizeof(uint32));
  for(i=0;i<n;i++) 
    occ[t[i]]++;

  // compute the mapping between new and old symbols
  asize=0;
  for(i=0;i<_BW_ALPHA_SIZE;i++) {
    if(occ[i]>0) {
      map[i] = asize++;
      // fprintf(stderr,"%c (%d occ) -> %d\n",i,occ[i],map[i]);
    } else
      map[i] = 2+_BW_ALPHA_SIZE;
  }

  // do the mapping
  for(i=0;i<n;i++) {
    c = t[i];
    assert(map[c]< 2+_BW_ALPHA_SIZE );
    newt[i] = map[c]; // actual remapping
  }

  return asize;
}
#endif


/* *******************************************************************
   given a string t[0...] (the endpoint is unknown) compute the
   (packed) bit array gt_eof[] such that for i=0,...,n-1
   gt_eof[i]==1 iff t[i...] > t[n...].
   
   This is useful becuase we need to compute the ordering of the suffixes
   t[i..] for i=0,n-1 and for doing this we want to use a classical 
   suffix array construction algorithm with the variant that the last symbol
   has a known rank (instead of having the smallest lex. rank). 
   Once we have the gt_eof bit array we are able to compute the
   rank of t[n-1...] within t[0...] ... t[n-1...].

   The computation of gt_eof is possible since in addition to 
   t[0,n) we are provided with t[n,2n) and the  
   (packed) bit vector gt_tn[] providing us with the rank of t[n...] with
   respect to t[j...] for j=n ... 2n-1 (see below). 
   The comparison between two suffixes t[i...] and t[j...] goes as follows:
     1. compare character by character. if a difference is found OK
     2. as soon as min(i+k,j+k)=n use gt_tn
   input 
     \param t:    input text. length=2n
     \param aux:  working space of n+1 uint32's
     \param n:    we have to sort the suffixes t[0...] to t[n-1...]
     \param gt_tn: (reads: larger_than_t[n...])
                   packed bit vector of length n,
                   Provides the rank of the suffix t[n...] within
                   the suffixes t[j...] for j=n ... 2n-1 in the sense
                   t[j...] > t[n...] iff gt_tn[j-n]==1
                   In practice it is a uint32 array of length ceil(n/32).
     \param gt_eof: address packed bit vector of size n where the
                    result is stored ( like it is a uint32 array of 
                    length ceil(n/32).)                     
   output
     the desired bit_vector is stored in gt_eof
   ******************************************************************* */
static void 
compute_gt_eof(uint8 *t, uint32 n, uint32 *aux, uint32 *gt_tn, uint32 *gt_eof)
{
  uint32 i,j,k,h,startj;
  uint32 *kmp_shift = aux;

  // compute kmp_shift for t[n,2n-1] and write it to kmp_shift=sa
  compute_kmp_shift(t+n,n,kmp_shift);
  //for(i=0;i<100;i++) fprintf(stderr,"%d ",kmp_shift[i]); //!!!!!!!

  // clear the gt_eof packed bit array (write all 0's)
  for(i=0;i<(n+31)/32;i++) gt_eof[i]=0;

  // compute the rank of t[n...] among the suffixes t[i...]  i=0..n-1
  // and write it to gt_eof so that t[i..n] > t[n..] iff gt_eof[i]==1
  for(i=startj=0;i<n;i++) {
    for(j=startj; ;j++) {
      if(i+j==n) { // we have reached t[n]
	assert(j>0 && j<n);
	if(bit_read32(gt_tn,j)==0) 
	  bit_set32(gt_eof,i); // gt_eof[i] = !gt_tn[j]
	break;
      }
      else if(t[i+j]!=t[n+j]) {
	if(t[i+j]>t[n+j])   bit_set32(gt_eof,i);
	break;
      }
      assert(t[i+j] == t[n+j]);// if t[i+j] == t[n+j] just continue
    }
    // mismatch at t[n+j]!=t[j]
    if(j==0) startj=0;
    else {                  // there was a match t[0..j-1] == t[n .. n+j-1]
      assert(j>0 && j<=n);  // j is the lenght of the match 
      k = kmp_shift[j-1];   // t[0..k-1] is a suffix of t[0..j-1] 
      assert(k<j );
      startj=k;             // t[n+k] is the next char to consider
      for(h=1;h<j-k;h++)    // copy  bits
	if(bit_read32(gt_tn,h)==1) 
	  bit_set32(gt_eof,i+h); // gt_eof[i+t] = gt_tn[t]
      i += j-k-1;           // skip already written positions 
    }
  }
  // check correctness of gt_eof with an expensive test 
  //for(i=0;i<n;i++)
  //  assert(bit_read(gt_eof,i)==tn_compare(t,i,n,gt_tn));

  return;

}





/* ***************************************************************
   compute the kmp shift for t[0..n-1].
   fill array shift[0...n-1] as follows:
     for i=1...n
       shift[i-1] = the largest j, 0<=j<i, such that t[0,j-1] == t[i-j,i-1]
       that is shift[i-1]=j  iff t[0,j-1] is a suffix of t[0,i-1] and it 
       is the longest prefix with this property.
       (note: j is the lenght of the match)         
   *************************************************************** */ 
static void compute_kmp_shift(uint8 *t, uint32 n, uint32 *shift)
{
  uint32 k,q;

  shift[0] = 0;
  k=0;
  for(q=1;q<n;q++) {
    // at iterarion q compute shift[q] 
    while(1) {
      // invariant: t[0,k-1]  suffix of t[0,q-1] that is
      // there is a lenght-k match with a suffix of t[0,q-1] 
      assert(k<q);
      if(t[q] == t[k]) {
	shift[q] = ++k; // it is t[0,k]  suffix of t[0,q]
	break;
      }
      if(k==0) {shift[q] = 0; break;}
      assert(k>0);
      assert(k>shift[k-1]);
      k = shift[k-1];
    }
  }
}








/* 
  **** alternative version of kmp_shift, correct but longer ****

  shift[1] = k = 0;
  for(q=1;q<n;q++) {
    // here we compute shift[q+1] 
    // invariant: t[0,k-1]  suffix of t[0,q-1] 
    assert(k<q);
    if(t[q] == t[k]) {
      // it is t[0,k]  suffix of t[0,q]
      shift[q+1] = ++k;
      continue;
    }
    while(k>0) { 
      assert(k<q);
      assert(t[q] != t[k]);
      k = shift[k];
      if(t[q] == t[k]) {
	shift[q+1] = ++k;
	break;
      }
    }
    if(k==0) shift[q+1] = 0;	
  }
*/




// ============================ No longer used

#if 0

static int check_suffixes(uint32 *t,uint32 i, uint32 j, uint32 n)
{
  uint32 k;

  for(k=0;k<n;k++) {
    assert(i+k<n && j+k<n);
    if (t[i+k]>t[j+k]) return 1;
    else if (t[i+k]<t[j+k]) return 0;
  }
  assert(0);
  return 2;
}


/* ****************************************************************
   compare t[i..] with t[n..] (useful only for debugging)
   return 1 iff t[i ..] > t[n ..] 
   **************************************************************** */ 
static uint32 tn_compare(uint8 *t, int i, int n, uint32 *gt_tn)
{
  int j;

  assert(i>=0 && i<n);
  for(j=0; ;j++) {
    if(i+j==n) return 1 - bit_read32(gt_tn,j);
    if(t[i+j]>t[n+j]) return 1;
    if(t[i+j]<t[n+j]) return 0;
  }
  assert(0);
  return 0;
}




/* *******************************************************************
   given a string t[0...] (the endpoint is unknown) compute the 
   ordering of the suffixes t[0...] up to t[n-1...], (n in total),
   and write it to sa[0,n). The sorting is possible since 
   in addition to t[0,n) we are provided with t[n,2n) and with a 
   (packed) bit vector gt_tn[] providing us with the rank of t[n...] with
   respect to t[j...] for j=n ... 2n-1 (see below). 
   The comparison between two suffixes t[i...] and t[j...] goes as follows:
     1. compare character by character. if a difference is found OK
     2. as soon as min(i+k,j+k)=n use gt_tn
   input 
     \param t:      input text. length=2n
     \param sa:     space for the sa.  length = n+1
     \param n:      we have to sort the suffixes t[0...] to t[n-1...]
     \param gt_tn:  (reads: larger_than_t[n...])
                    logically it is a bit vector of length n,
                    Provides the rank of the suffix t[n...] within
                    the suffixes t[j...] for j=n ... 2n-1 in the sense
                    t[j...] > t[n...] iff gt_tn[j-n]==1
                    In practice it is a uint32 array of length ceil(n/32).
                    If gt_tn==NULL then t[n..] should be considered
                    the lex. smallest suffix (used for first sa computation)
     \param gt_eof: available (packed) bit array of size n
     \param isa:    available array of n+1 uint32's. Important note:
                    the memoru used by isa overlaps with the memory used by
                    t, gt_tn and gt_eof so they cannot be used simultaneously. 
   output
     the sorted order of suffixes  t[0 ...] to t[n-1...] is stored sa[0..n-1]
     t[0,n) is restored, but t[n,2n) and gt_tn are destroyed (overwritten)
   ******************************************************************* */
void _bwtext_isufsort(uint8 *t, uint32 *sa, uint32 n, 
		      uint32 *gt_tn, uint32 *gt_eof, uint32 *isa) 
{
  uint32 occ[_BW_ALPHA_SIZE];
  uint32 i,j, newasize;

  // count occurrences of each symbol
  memset( occ, 0, _BW_ALPHA_SIZE* sizeof(uint32) );
  for(i=0;i<n;i++) occ[t[i]]++;

  // compute the packed bit vector gt_eof[] such that
  // for i=0,..,n-1  gt_eof[i]==1 iff t[i..] > t[n..]
  if(gt_tn==NULL) 
    for(i=0;i<n;i++)
      bit_set32(gt_eof,i);
  else
    compute_gt_eof(t,n,sa,gt_tn,gt_eof);


  // !!!!!!!!!!!!!!!!!!!!!
#if 0
  for(i=0;i<10;i++)
    fprintf(stderr,"%d",bit_read32(gt_tn,i));
  fprintf(stderr,"xxx\n");
  for(i=n-10;i<n;i++)
    fprintf(stderr,"%d",bit_read32(gt_eof,i));
  fprintf(stderr,"xxx\n");
#endif

  // remap symbols in t[] according to gt_eof bits: the new string is written
  // to sa[] and the remapping is such that suffix-sorting of sa[0,n-1]
  // is equivalent to sorting the suffixes t[0..], t[n-1...]
  newasize = remap_alphabet(t,n,sa,gt_eof);

  // now we no longer need t, gt_tn, and gt_eos
  // copy the remapped text to isa
  for(i=0;i<n;i++) isa[i] = sa[i];

  // do the sorting (for now with larsson-sadakane)
  ls_sufsort(isa,sa,n,newasize);

  // now sa contains the suffix array and isa its inverse
  // recover the text in t (overwriting isa)
  j=0;
  for(i=0;i<n;i++) {
    while(occ[j]==0) {         // search next nonzero in occ 
      j++; assert(j<_BW_ALPHA_SIZE);
    }
    t[sa[i]] = (uint8) j;
    occ[j]--;
  }

  // check that all symbols have been assigned 
  assert(occ[j]==0);
  for(i=j+1;i<_BW_ALPHA_SIZE;i++)
    assert(occ[i]==0);

  // now t[0,n) contains the text and sa[0,n) the suffix array
  return;
}

static void ls_sufsort(uint32 *isa, uint32 *sa, uint32 n, uint32 asize)
{
  //uint32 *salvatesto; //!!!!!!!!!!!
  assert(n< (1u<<31));
  assert(n>0);

  //salvatesto = (uint32 *) malloc(n*4);
  //assert(salvatesto!=NULL);

  #if 0
  {uint32 i;//!!!!!!!!!!!
  for(i=0;i<n;i++) { 
    salvatesto[i]=isa[i];
    if(i<n-1) assert(isa[i]!=isa[n-1]);
  }
  }
  #endif

  _bwtext_ls_suffixsort((int32 *) isa, (int32 *) sa, (int) n-1,(int) asize);

  #if 0
  {uint32 i;
  for(i=0;i<n-1;i++) 
    if(check_suffixes(salvatesto,sa[i],sa[i+1],n)!=0)
      assert(1);
  }
  #endif

  // free(salvatesto); !!!!!!!!!!!!
}

#endif
