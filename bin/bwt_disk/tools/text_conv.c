/* *********************************************************************
   conversion between the various formats supported by the bwtext library 
   Giovanni Manzini  (giovanni.manzini@unipmn.it)
   16 Jun 2009
   ********************************************************************* */
#include "../bwtext.h"
#include <unistd.h>

// parameters to be obtained from the command line 
extern int _bwtext_Verbose;                // verbosity level
extern double _bwtext_FilterPar1;          // parameter #1
extern int _bwtext_FilterPar1_flag;        // !=0 if parameter #1 was defined
extern double _bwtext_FilterPar2;          // parameter #2
extern int _bwtext_FilterPar2_flag;        // !=0 if parameter #2 was defined

// note: these filter paramaters are used by the output filter 
// or, more in general, by all filters used in write mode 


static void convert(char *bnam, filter *inf, char *tnam, filter *outf);

int main(int argc, char *argv[])
{
  char *inam=NULL, *outnam=NULL; 
  filter *inf,*outf;
  extern char *optarg;
  extern int optind, opterr, optopt;
  int c, input_filter, output_filter, num_opt;

  /* ------------- read options from command line ----------- */
  num_opt = opterr = 0;
  output_filter = 0; input_filter = -1;
  _bwtext_Verbose=0;
  _bwtext_FilterPar1_flag = _bwtext_FilterPar2_flag = 0;
  while ((c=getopt(argc, argv, "vi:o:y:z:")) != -1) {
    switch (c) 
      {
      case 'i':
        input_filter = atoi(optarg); break;
      case 'o':
        output_filter = atoi(optarg); break;
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
  if(optind+1<argc && argv[optind]!=NULL && argv[optind+1]!=NULL) {
    inam=argv[optind];
    outnam=argv[optind+1];
  }
  else {
    fprintf(stderr, 
	    "Usage:  %s [-i fil][-o fil]", argv[0]);
    fprintf(stderr,"[-y val][-z val][-v] infile outfile\n");
    fprintf(stderr,"\t-v          produces a verbose output\n");
    fprintf(stderr,"\t-i fil      input filter [def. autodetect]\n");
    fprintf(stderr,"\t-o fil      output filter [def. none]\n");
    fprintf(stderr,"\t-y val      output filter parameter #1\n");
    fprintf(stderr,"\t-z val      output filter parameter #2\n\n");
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
  if(!(outf = _bwtext_get_filter(output_filter,outnam)))
    _bwtext_fatal_error("Invalid output filter (text_conv/main)");
  if(!(inf=_bwtext_get_filter(input_filter,inam)))
    _bwtext_fatal_error("Invalid input filter (text_conv/main)");
  if(_bwtext_Verbose) {
    fprintf(stderr,"Using input filter: %s\n",inf->name);
    fprintf(stderr,"Using output filter: %s\n",outf->name);
  }

  // do the conversion
  convert(inam,inf,outnam,outf);
  return 0;
}



#define BSIZE 1000000
static void convert(char *innam, filter *inf, char *outnam, filter *outf)
{
  void *infile;
  void *outfile;
  uint8 buffer[BSIZE];
  int32 n;
  int64 usize=0, csize;

  infile = inf->open(innam,'r');     // open input file
  if(infile==NULL) {
    perror(innam);
    _bwtext_fatal_error("Invalid input file (text_conv/convert)");
  }
  outfile = outf->open(outnam,'w');  // open output file
  if(outfile==NULL) {
    perror(outnam);
    _bwtext_fatal_error("Invalid output file (text_conv/convert)");
  }

  do {
    n = inf->read(infile,buffer,BSIZE);
    outf->write(outfile,buffer,n);
    usize += n;
  } while(n==BSIZE);

  csize = outf->close(outfile);
  inf->close(infile);

  if(_bwtext_Verbose) {
    fprintf(stderr,"Original file size: %lld\n", (long long) usize);
    fprintf(stderr,"Compressed file size (extimated): %lld\n", 
	    (long long) csize);
  }

}
 


