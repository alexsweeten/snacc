/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   functions for handling a bwt file via a filter
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include "bwtext.h"
#include <string.h>

// external variables providing access to the bwt file
extern char *_bwtext_Bwtname;
extern filter *_bwtext_Bwt_filter;

// local prototypes
static char *generate_bwt_temp_filename(char *fname);
static char *generate_bwt_final_filename(char *fname);
static void _bwtext_bwt_deallocate(bwtfile *f, boolean close, boolean delete);
static void _bwtext_bwt_create_temp(uint64 size, uint64 eof);
static int64 _bwtext_bwt_close_temp();
static void write64(bwtfile *, uint64);
static uint64 read64(bwtfile *);
static void block_transfer(bwtfile *in, bwtfile *out, uint64 n,uint32 c);
static int write64_le(uint64 u, FILE *f);


#define BSIZE _bwtext_BWT_BUFFER_SIZE

// structure representing the current bwt file
static bwtfile *b=NULL;

/* >>>>>>>>>>> functions used for the I/O of the bwt >>>>>>>>>>>>> */


/* ***************************************************************
   create a temporary bwt files with lenght=n, eofpos=eof and data
   s[0,n-1]. the bwt file is assigned to the static variable b
   return the number of bytes written in the bwt file
   *************************************************************** */
int64 _bwtext_bwt_init_temp(uint32 n, uint32 eof, uint8 *s)
{
  _bwtext_bwt_create_temp(n,eof);        // create temp

  assert(n<(1u<<31));
  if(b->fil->write(b->filep,s,n)!= ((int32) n)) {       // write data
    perror(b->filename);
    _bwtext_ioerror(b->filename,"_bwtext_bwt_init_temp");
  }

  return _bwtext_bwt_close_temp();              // close temp file

} 

/* *********************************************************************
   reopen in read only mode the current bwt file. closing and reopening is
   equivalent to a rewind, but for some filters (eg gzip) a rewind is only
   possible for read only files. 
   read size and eof and set pointer to the beginning of bwt data
   ********************************************************************* */
static void _bwtext_bwt_reopen(uint64 *size, uint64 *eof)
{

  b->filep = b->fil->open(b->filename,'r');
  if(b->filep==NULL) {
    perror(b->filename);
    _bwtext_ioerror(b->filename,"_bwtext_bwt_reopen");
  }
  *size = read64(b);
  *eof  = read64(b);
}


/* ******************************************************************
   rename the current temp bwt file to the final bwt name 
   also deallocate space for the current temp bwt file
   ****************************************************************** */
void _bwtext_bwt_finalize_temp() 
{
  char *s= _bwtext_Bwtname;
  char *final_name;

  final_name = generate_bwt_final_filename(s);
  if(rename(b->filename,final_name)==-1) {
    perror(final_name);
    _bwtext_ioerror(b->filename,"_bwtext_bwt_finalize_temp");
  }
  free(final_name);
  _bwtext_bwt_deallocate(b,0,0);
} 

/* *****************************************************************
   write a file containing the number of occurrences of each symbol
   in the input file written in little-endian format
   ***************************************************************** */ 
void _bwtext_bwt_write_auxfile(uint64 *occ_global)
{
  char *s= _bwtext_Bwtname, *final_name, *aux_name;
  FILE *aux; int i;

  // create aux file name
  final_name = generate_bwt_final_filename(s);
  asprintf(&aux_name,"%s.aux",final_name);
  // open aux file
  aux = fopen(aux_name,"wb");
  if(aux==NULL) {
    perror(aux_name);
    _bwtext_ioerror(aux_name,"_bwtext_write_auxfile");
  }
  // write occ_global[] to file aux in little-endian format
  for(i=0;i<_BW_ALPHA_SIZE;i++) {
    uint64 u = occ_global[i]; 
    if(write64_le(u, aux)) { 
      perror(aux_name);
      _bwtext_ioerror(aux_name,"_bwtext_write_auxfile");
    }
  }
  fclose(aux);
  free(aux_name);
  free(final_name);
}


/* ********************************************************************
   merge the bwt stored in b with the one in cur_bwt[0,n). The 
   merging is possible since gaps[0,n] contains the number of suffixes 
   in b between any pair of suffixes in cur_bwt[].
   cur_rank0 is the lex. of the longest suffix in cur_bwt[] and it is
   used to determine the position of the eof symbol in the merged bwt
   cur_last is the length-1 suffix of cur_bwt[] and it is used to replace
   the eof sumbol of the bwt stored in b  
   ******************************************************************** */
int64 _bwtext_bwtmerge_temp(uint32 *gaps, uint32 *cur_bwt, uint32 n, 
				  uint32 cur_rank0, uint32 cur_last)
{
  uint64 old_size, old_eof, new_size, new_eof, i,tot,gi;
  uint32 next_bwt_char;
  bwtfile *old;

  // rewind current bwt file and save old values
  _bwtext_bwt_reopen(&old_size,&old_eof);
  assert(old_size==b->size && old_eof==b->eofpos);
  // set old to current and current to NULL;
  old=b; b=NULL;

  // compute new size and eof  
  new_size = old_size+n;
  tot=n;                     // do a redundant check on size_new
  for(i=0;i<=n;i++)
    tot += gaps[i];
  assert(tot==new_size); 

  new_eof = 0;
  for(i=0;i<cur_rank0;i++)
    new_eof += gaps[i]+1;    // we skip bwt[i] + gaps[i] preceding suffixes
  new_eof += gaps[i];        // skip suffixes immediately before cur_rank0
  assert(new_eof>0);
  
  // create new bwt file and make it the current one
  _bwtext_bwt_create_temp(new_size,new_eof);

  // copy from old bwt to new bwt
  assert(cur_last<256);
  for(tot=i=0;i<=n;i++) {
    // convert gaps[i] to a uint64 (in the future gaps[] could be more complex)
    gi = gaps[i];
    // compute next bwt char taking into account new_eof and last gap 
    if(i<n) {
      if(i==cur_rank0) assert(tot+i+gi==new_eof);
      next_bwt_char = cur_bwt[i] & 255;  // delete high oreder bits
      assert(next_bwt_char<256);
    }
    else next_bwt_char = 256;

    if(tot>old_eof || (tot + gi <= old_eof)) 
      block_transfer(old,b,gi,next_bwt_char);
    else { // eof was in this block 
      block_transfer(old,b,old_eof-tot,cur_last);
      _bwtext_Bwt_filter->read(old->filep,(uint8 *)&cur_last,1); // skip eof
      block_transfer(old,b, gi-(old_eof-tot)-1,next_bwt_char);     
    }
    tot += gi;
  }   
  assert(tot==old_size);
  int64 written = _bwtext_bwt_close_temp();
  _bwtext_bwt_deallocate(old, 1, 1); // close and delete old bwt file
  return written;
} 

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   functions providing read access to a bwt file
   These are used by the unbwt program
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */

// open a bwt file for reading
bwtfile *_bwtext_bwt_open_read(char *fname, filter *f)
{
  bwtfile *b;

  // allocate new structure for file
  b = (bwtfile *) malloc(sizeof(bwtfile));
  if(!b) _bwtext_out_of_mem("_bwtext_bwt_open_read");

  // set filename
  b->filename = strdup(fname);
  assert(b->filename!=NULL);
  // set filter
  b->fil = f;

  // open file
  b->filep = b->fil->open(b->filename,'r');
  if(b->filep==NULL) {
    perror(b->filename);
    _bwtext_ioerror(b->filename,"_bwtext_bwt_open_read");
  }
  
  // store in b file size and eofpos
  b->size = read64(b);
  b->eofpos = read64(b);
  return b;
}


// read from the current bwt file up to n chars in a user-supplied buffer 
int32 _bwtext_bwt_read(bwtfile *b, uint8 *buf, int32 n)
{
  int32 m;

  if( (m=b->fil->read(b->filep,buf,n))==-1) {
    perror(b->filename);
    _bwtext_ioerror(b->filename,"_bwtext_bwt_read");
  }
  assert(m<=n);
  return m;
}

// close a bwt file 
void _bwtext_bwt_close(bwtfile *b)
{
  if(b->fil->close(b->filep)<0)
    _bwtext_ioerror(b->filename,"_bwtext_bwt_close");
  free(b->filename);
  free(b);
}




/* >>>>>>>>>>>>>>>>>>>>>>>> static functions >>>>>>>>>>>>>>>>>>>> */


/* *****************************************************************
   create a new temporary bwt file and store its data in b
   ***************************************************************** */
static void _bwtext_bwt_create_temp(uint64 size, uint64 eof)
{
  char *s= _bwtext_Bwtname;

  // allocate new structure for current file
  assert(b==NULL);
  b = (bwtfile *) malloc(sizeof(bwtfile));
  if(!b) _bwtext_out_of_mem("_bwtext_bwt_create_temp");

  // set filter to default filter
  b->fil = _bwtext_Bwt_filter;
  // create appropriate filename
  b->filename = generate_bwt_temp_filename(s);
  assert(b->filename!=NULL);

  // open file
  b->filep = b->fil->open(b->filename,'w');
  if(b->filep==NULL) {
    perror(b->filename);
    _bwtext_ioerror(b->filename,"_bwtext_bwt_create_temp");
  }
  
  // store in b and write to file size and eof
  b->size=size; b->eofpos=eof;
  write64(b,size);
  write64(b,eof);
  return;
}

// close the current bwt file return the number of bytes written so far
// (this could be less than the actual file size due to buffering )
int64 _bwtext_bwt_close_temp()
{
  int64 ris = b->fil->close(b->filep); 
  if(ris<0) 
    _bwtext_ioerror(b->filename,"_bwtext_bwt_close");
  return ris;
}


/* ************************************************************************
   transfer n bytes from in to out. Then, if c<256, c is written to out
   ************************************************************************ */
static void block_transfer(bwtfile *in, bwtfile *out, uint64 n, uint32 c)
{
  uint8 buffer[BSIZE];
  int64 r;

  assert(BSIZE<4294967295U);
  //buffer = (uint8 *) malloc(BSIZE*sizeof(*buffer));
  if(buffer==NULL) _bwtext_out_of_mem("block_transfer");
  while(n>=BSIZE) {
    if( in->fil->read(in->filep,buffer,BSIZE)!=BSIZE) {
      perror(in->filename);
      _bwtext_ioerror(in->filename,"block_transfer");
    }
    if( out->fil->write(out->filep,buffer,BSIZE)!=BSIZE) {
      perror(out->filename);
      _bwtext_ioerror(out->filename,"block_transfer");
    } 
    n -= BSIZE;
  }
  // read a block of less than BSIZE bytes
  assert(n<BSIZE);
  r = (int64) n;
  if(r>0)
    if(in->fil->read(in->filep,buffer,r)!=r) {
      perror(in->filename);
      _bwtext_ioerror(in->filename,"block_transfer");
    }
  if(c<256)  buffer[r++] = (uint8) c; // add c if < 256
  if( out->fil->write(out->filep,buffer,r)!=r) {
    perror(out->filename);
    _bwtext_ioerror(out->filename,"block_transfer");
  }
  // free(buffer);
}

// free space used by a bwtfile and optionally delete the associated file
static void _bwtext_bwt_deallocate(bwtfile *f, boolean close, boolean delete)
{
 
  assert(f!=NULL);
  // delete=0;      // uncomment this to keep the bwt tmp files 
  if(close) {
    if(f->fil->close(f->filep)<0) 
      _bwtext_ioerror(f->filename,"_bwtext_bwt_deallocate");
  }
  if(delete) {
    if(remove(f->filename)==-1) {
      perror(f->filename);      
      _bwtext_ioerror(b->filename,"_bwtext_bwt_deallocate");
    }
  }
  free(f->filename);
  free(f);
}

static char *generate_bwt_temp_filename(char *fname) 
{
  static char *base=NULL;
  static int count=0;
  char *name;

  if(base==NULL) base = fname;
  assert(base!=NULL);
  asprintf(&name,"%s.%d.tmp.bwt", base, count++);
  return name;
}

static char *generate_bwt_final_filename(char *fname) 
{
  char *name;

  assert(fname!=NULL);
  asprintf(&name,"%s", fname);
  return name;
}

// read and write 64 bits from f
static uint64 read64(bwtfile *f)
{
  uint8 buffer[8]; int i;
  uint64 u;

  if(f->fil->read(f->filep,buffer,8)!=8) {
    perror(f->filename);
    _bwtext_ioerror(f->filename,"read64");
  }
  for(u=0,i=7;i>=0;i--) {
    u = u<<8;         // shift old content
    u += buffer[i];   // add one byte 
  }
  return u;
}
static void write64(bwtfile *f, uint64 u)
{
  uint8 buffer[8]; int i;

  for(i=0;i<8;i++) {
    buffer[i] = (uint8) (u & 0xFF); // take lsb 
    u = u>>8;                       // throw lsb away
  }
  if(f->fil->write(f->filep,buffer,8)!=8) {
    perror(f->filename);
    _bwtext_ioerror(f->filename,"read64");
  }
}


// write a uint64 in little-endian format to file f
// return 1 on error, 0 if successful  
static int write64_le(uint64 u, FILE *f) 
{
  int j,c;

  for(j=0;j<8;j++) {
    c = fputc((int) (u & 0xFF),f);    // output lsb  
    u = u>>8;                         // throw lsb away
    if(c==EOF) return 1;
  }
  return 0;
} 


