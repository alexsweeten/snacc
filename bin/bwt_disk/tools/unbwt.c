/* *********************************************************************
   In memory inversion of a bwt computed using the bwtext library 
   Giovanni Manzini  (giovanni.manzini@unipmn.it)
   27 Nov 2007
   04 Mar 2009
   27 Jan 2010

    Copyright (C) 2007,2008,2009,2010  Giovanni Manzini

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

   ********************************************************************* */
#include "../bwtext.h"
#include <unistd.h>

// parameters to be obtained from the command line 
extern int _bwtext_Verbose;                // verbosity level
extern double _bwtext_FilterPar1;          // parameter #1
extern int _bwtext_FilterPar1_flag;        // !=0 if parameter #1 was defined
extern double _bwtext_FilterPar2;          // parameter #2
extern int _bwtext_FilterPar2_flag;        // !=0 if parameter #2 was defined

static  uint32 get_char_from_first_column64(uint64 i, uint64 *first_col);


int main(int argc, char *argv[])
{
  void _bwtext_invert_bwt(char *bnam, filter *bf, char *tnam, filter *tf);
  char *bnam=NULL, *tnam=NULL; 
  filter *bf,*tf;
  extern char *optarg;
  extern int optind, opterr, optopt;
  int c, text_filter, bwt_filter, num_opt;

  /* ------------- read options from command line ----------- */
  num_opt = opterr = 0;
  text_filter = 0; bwt_filter = -1;
  _bwtext_Verbose=0;
  _bwtext_FilterPar1_flag = _bwtext_FilterPar2_flag = 0;
  while ((c=getopt(argc, argv, "vb:t:o:y:z:")) != -1) {
    switch (c) 
      {
      case 't':
        text_filter = atoi(optarg); break;
      case 'b':
        bwt_filter = atoi(optarg); break;
      case 'o':
        tnam = optarg; break;
      case 'v':
        _bwtext_Verbose++; break;
      case 'y':
        _bwtext_FilterPar1 = atof(optarg);
	_bwtext_FilterPar1_flag =1; break;
      case 'z':
        _bwtext_FilterPar2 = atof(optarg);
	_bwtext_FilterPar2_flag =1; break;
      case '?':
        fprintf(stderr,"Unknown option: %c -unbwt/main-\n", optopt);
        exit(1);
      }
    num_opt++;
  }
  if(optind<argc && argv[optind]!=NULL) {
    bnam=argv[optind];
  }
  else {
    fprintf(stderr, 
	    "Usage:  %s [-v][-b fil][-t fil][-o textfile]", argv[0]);
    fprintf(stderr,"[-y val][-z val] bwtfile\n");
    fprintf(stderr,"\t-v          produces a verbose output\n");
    fprintf(stderr,"\t-b fil      bwt file filter [def. autodetect]\n");
    fprintf(stderr,"\t-t fil      text file filter [def. none]\n");
    fprintf(stderr,"\t-o textfile write text to textfile\n");
    fprintf(stderr,"\t-y val      text file filter parameter #1\n");
    fprintf(stderr,"\t-z val      text file filter parameter #2\n\n");
    _bwtext_print_available_filters(stderr);
    fprintf(stderr,"\n");
    return 0;
  }
  if(_bwtext_Verbose) {
    fprintf(stderr,"Command line: ");
    for(c=0;c<argc;c++)
      fprintf(stderr,"%s ",argv[c]);
    fprintf(stderr,"\n");
  }

  // initialize filters
  if(!(tf = _bwtext_get_filter(text_filter,tnam)))
    _bwtext_fatal_error("Invalid text filter (unbwt/main)");
  if(!(bf=_bwtext_get_filter(bwt_filter,bnam)))
    _bwtext_fatal_error("Invalid bwt filter (unbwt/main)");
  if(_bwtext_Verbose) {
    fprintf(stderr,"Using bwt filter: %s\n",bf->name);
    fprintf(stderr,"Using text filter: %s\n",tf->name);
  }
  // init text file name 
  if(tnam==NULL) 
    asprintf(&tnam,"%s.twb%s",bnam,tf->extension);


  // do the inversion
  _bwtext_invert_bwt(bnam,bf,tnam,tf);
  return 0;
}



#define BSIZE 1000000
void _bwtext_invert_bwt(char *bnam, filter *bf, char *tnam, filter *tf)
{
  bwtfile *b;
  textfile *t;
  uint64 occ[_BW_ALPHA_SIZE];      // emulation of first column bwt matrix 
  uint64 start_sa_range[_BW_ALPHA_SIZE+1];
  uint64 i,tot,read,textsize;
  uint32 *rankprev, *bwt;
  uint8 buffer[BSIZE];
  int32 n,j;

  b = _bwtext_bwt_open_read(bnam,bf);
  // fprintf(stderr,"bwt size: %lld, eofpos: %lld\n",(long long)b->size,
  //                                                 (long long) b->eofpos);
  if(b->size>4294967294U)
    _bwtext_out_of_mem("BWT too large for internal memory inversion");
   bwt = rankprev = (uint32 *) malloc(b->size*sizeof(*rankprev));
  if(rankprev==NULL)
    _bwtext_out_of_mem("BWT too large for internal memory inversion");
  
  // read bwt from file
  for(read=0;read<b->size;read+=n) {
    n = min(BSIZE,b->size-read);
    j = _bwtext_bwt_read(b,buffer,n);
    assert(j==n);
    for(j=0;j<n;j++)
      bwt[read+j] = buffer[j];
  }

  for(i=0;i<_BW_ALPHA_SIZE;i++)
    occ[i]=0;
  for(i=0;i<b->size;i++)
    if(i!=b->eofpos) occ[bwt[i]]++;
  
  for(i=0,tot=0;i<_BW_ALPHA_SIZE;i++) {
    start_sa_range[i] = tot; tot+= occ[i];
  }
  start_sa_range[_BW_ALPHA_SIZE] = tot;  
  assert(tot==b->size-1);

  // bwt -> rankprev inplace
  for(i=0;i<b->size;i++)
    if(i!=b->eofpos)
      rankprev[i] = 1+start_sa_range[bwt[i]]++; // 1 is for the eof
    else 
      rankprev[i] = 0;
  assert(start_sa_range[_BW_ALPHA_SIZE-1]==start_sa_range[_BW_ALPHA_SIZE]);
 
  // recover the original start_sa_range values
  tot=0;
  for(i=0;i<_BW_ALPHA_SIZE;i++) {
    uint64 temp = start_sa_range[i];
    start_sa_range[i]=tot;
    tot=temp;
  }
  assert(tot==start_sa_range[_BW_ALPHA_SIZE]);
  textsize = b->size-1;  // save text size
  _bwtext_bwt_close(b);  // close bwt file
  
  // output original text 
  t = _bwtext_text_open_write(tnam,tf,1000000);
  for(i=rankprev[0],tot=0;tot<textsize;i=rankprev[i],tot++) {
    assert(i>0);
    uint32 c = get_char_from_first_column64(i-1,start_sa_range);
    _bwtext_text_putc(t,c);
  }
  assert(i==0);
  _bwtext_text_close_write(t);
  free(rankprev);
}


// given an index returns the corresponding char in the bwt matrix
// to be improved using binary search
static uint32 get_char_from_first_column64(uint64 index, uint64 *first_col)
{
  uint64 i;

  for(i=0;i<_BW_ALPHA_SIZE;i++)
    if(first_col[i]<=index && index < first_col[i+1])
      return i;
  _bwtext_fatal_error("get_char_from_first_column64");
  return 0;  
}

