/* *********************************************************************
   contruction of the auxiliary file containing the number of occurrences
   of each character in the input file. 
   Giovanni Manzini  (giovanni.manzini@unipmn.it)
   24 Jan 2010
   ********************************************************************* */
#include "../bwtext.h"
#include <unistd.h>

// parameters to be obtained from the command line 
extern int _bwtext_Verbose;                // verbosity level



static void count(char *fnam, filter *fil, char *auxnam);
static int write64_le(uint64 u, FILE *f);


int main(int argc, char *argv[])
{
  char *inam=NULL, *outnam=NULL; 
  filter *inf;
  extern char *optarg;
  extern int optind, opterr, optopt;
  int c, input_filter, num_opt;

  /* ------------- read options from command line ----------- */
  num_opt = opterr = 0;
  input_filter = -1;
  _bwtext_Verbose=0;
  while ((c=getopt(argc, argv, "vi:")) != -1) {
    switch (c) 
      {
      case 'i':
        input_filter = atoi(optarg); break;
      case 'v':
        _bwtext_Verbose++; break;
      case '?':
        fprintf(stderr,"Unknown option: %c -text_count/main-\n", optopt);
        exit(1);
      }
    num_opt++;
  }
  if(optind+1<argc && argv[optind]!=NULL && argv[optind+1]!=NULL) {
    inam=argv[optind];
    outnam=argv[optind+1];
  }
  else {
    fprintf(stderr, 
	    "Usage:  %s [-i fil][-v] infile outfile\n", argv[0]);
    fprintf(stderr,"\t-i fil      input filter [def. autodetect]\n");
    fprintf(stderr,"\t-v          produces a verbose output\n\n");
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

  if(!(inf=_bwtext_get_filter(input_filter,inam)))
    _bwtext_fatal_error("Invalid input filter (text_count/main)");
  if(_bwtext_Verbose) {
    fprintf(stderr,"Using input filter: %s\n",inf->name);
  }

  // do the conversion
  count(inam,inf,outnam);
  return 0;
}



#define BSIZE 1000000
static void count(char *innam, filter *inf, char *outnam)
{
  void *infile;
  uint8 buffer[BSIZE];
  int32 n,i;
  int64 usize=0;
  uint64 occ[_BW_ALPHA_SIZE];  

  infile = inf->open(innam,'r');     // open input file
  if(infile==NULL) {
    perror(innam);
    _bwtext_fatal_error("Invalid input file (text_count/count)");
  }

  for(i=0;i<_BW_ALPHA_SIZE;i++) occ[i]=0; // clear occ[]

  do {  // compute occ[]
    n = inf->read(infile,buffer,BSIZE);
    for(i=0;i<n;i++)  occ[buffer[i]]++;
    usize += n;
  } while(n==BSIZE);

  inf->close(infile);

  // open output file
  FILE *out = fopen(outnam,"wb");
  if(out==NULL) {
    perror(outnam);
    _bwtext_fatal_error("Invalid output file (text_count/count)");
  }
  // write occ to output file
  for(i=0;i<_BW_ALPHA_SIZE;i++) {
    if(write64_le(occ[i],out)) {
      perror(outnam);
      _bwtext_fatal_error("Write file error (text_count/count)");
    }
  }
  fclose(out);

  if(_bwtext_Verbose) 
    fprintf(stderr,"Original file size: %jd\n", (intmax_t) usize);

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

