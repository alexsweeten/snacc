/* *********************************************************************
   tool for reverting the content of a file 
   Giovanni Manzini  (giovanni.manzini@unipmn.it)
   26 Jan 2010
   ********************************************************************* */
#include "../bwtext.h"
#include <unistd.h>

// parameters to be obtained from the command line 
extern int _bwtext_Verbose;                // verbosity level
extern double _bwtext_FilterPar1;          // parameter #1
extern int _bwtext_FilterPar1_flag;        // !=0 if parameter #1 was defined
extern double _bwtext_FilterPar2;          // parameter #2
extern int _bwtext_FilterPar2_flag;        // !=0 if parameter #2 was defined


static void reverse(char *inam, char *onam, filter *outf);


int main(int argc, char *argv[])
{
  char *inam=NULL, *outnam=NULL; 
  filter *outf;
  extern char *optarg;
  extern int optind, opterr, optopt;
  int c, output_filter, num_opt;

  /* ------------- read options from command line ----------- */
  num_opt = opterr = 0;
  output_filter = 0;
  _bwtext_Verbose=0;
  _bwtext_FilterPar1_flag = _bwtext_FilterPar2_flag = 0;
  while ((c=getopt(argc, argv, "vo:y:z:")) != -1) {
    switch (c) 
      {
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
    fprintf(stderr,"Usage:  %s [-o fil]", argv[0]);
    fprintf(stderr,"[-y val][-z val][-v] infile outfile\n");
    fprintf(stderr,"\t-v          produces a verbose output\n");
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
    _bwtext_fatal_error("Invalid output filter (text_rev/main)");
  if(_bwtext_Verbose) {
    fprintf(stderr,"Using output filter: %s\n",outf->name);
  }

  // do the conversion
  reverse(inam,outnam,outf);
  return 0;
}



#define BSIZE 1000000
static void reverse(char *innam, char *outnam, filter *outf)
{
  FILE *infile;
  void *outfile;
  uint8 tmp, buffer[BSIZE];
  int32 n,i;
  int64 read=0, csize, usize, offset;

  infile = fopen(innam,"rb");        // open input file
  if(infile==NULL) {
    perror(innam);
    _bwtext_fatal_error("Invalid input file (text_rev/reverse)");
  }
  outfile = outf->open(outnam,'w');  // open output file
  if(outfile==NULL) {
    perror(outnam);
    _bwtext_fatal_error("Invalid output file (text_rev/reverse)");
  }

  if(fseek(infile,0,SEEK_END)) {     // get input file size
    perror(innam);
    _bwtext_fatal_error("Error in fseek (text_rev/reverse)");
  }
  usize = ftell(infile);
  if(usize == -1) {
    perror(innam);
    _bwtext_fatal_error("Error in ftell (text_rev/reverse)");
  }

  // --------------- reverse content 
  do {                               
    offset = usize - read - BSIZE;
    if(offset>=0) {                  // there are at least BSIZE chars left
      if(fseek(infile,offset,SEEK_SET)) {
	perror(innam);
	_bwtext_fatal_error("Error in fseek (text_rev/reverse)");
      }
      n = fread(buffer,1,BSIZE,infile);
      if(n!=BSIZE) {
	perror(innam);
	_bwtext_fatal_error("Error in fread (text_rev/reverse)");
      }
    }
    else {                         // less than BSIZE chars left
      if(fseek(infile,0,SEEK_SET)) {
	perror(innam);
	_bwtext_fatal_error("Error in fseek (text_rev/reverse)");
      }
      n = fread(buffer,1,BSIZE+offset,infile);
      if(n!=BSIZE+offset) {
	perror(innam);
	_bwtext_fatal_error("Error in fread (text_rev/reverse)");
      }
    }
    // reverse chars inside buffer[]
    for(i=0;i<n-1-i;i++) {
      tmp = buffer[i];
      buffer[i] = buffer[n-1-i];
      buffer[n-1-i] = tmp;
    }
    // write buffer
    outf->write(outfile,buffer,n);
    read += n;
  } while(n==BSIZE);
  assert(read==usize);

  csize = outf->close(outfile);
  fclose(infile);

  if(_bwtext_Verbose) {
    fprintf(stderr,"Original file size: %jd\n", (intmax_t) usize);
    fprintf(stderr,"Compressed file size (extimated): %jd\n", 
	    (intmax_t) csize);
  }

}
 


